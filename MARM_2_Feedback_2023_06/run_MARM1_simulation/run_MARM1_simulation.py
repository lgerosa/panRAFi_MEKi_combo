#this Python script runs MARM1 simulations

#import standard support functions
import sys
import numpy as np
import pandas as pd
import seaborn as sns

#import custom support functions
from equilibrate import equilibrate
from get_species_index import get_species_index

#impor  model 
from pysb.simulator import ScipyOdeSimulator
from pysb.core import as_complex_pattern
from MARM1 import model

def run_MARM1_simulation(params_file, param_set_index):

    #define file output name
    output_file='./sim_output_files/sim_results_param_%d.csv' % param_set_index

    #define conditions
    n_doses=11;
    lb_RAFi=-4
    ub_RAFi=1
    lb_MEKi=-5
    ub_MEKi=0
    RAFi_concentration = np.append(np.array([0.0]), np.logspace(lb_RAFi, ub_RAFi, n_doses))
    MEKi_concentration = np.append(np.array([0.0]), np.logspace(lb_MEKi, ub_MEKi, n_doses))
    t_pretrt = [24]
    EGF_concentration = [0, 100]
    t_trt = [2]
    param_set_index = [param_set_index]
    N_time_points = 97
 
    #generate condition settings for multiple doses
    settings_list=[]
    for param in param_set_index:
        for rafi in RAFi_concentration:
            for meki in MEKi_concentration:
                for pretrt in t_pretrt:
                    for trt in t_trt:
                        for egfc in EGF_concentration:
                            settings_list.append([param, pretrt, rafi, meki, trt, egfc]) 

    #generate condition settings for EGF stimulations after varying time of RAFi inhibition
    for param in param_set_index:
        RAFi_concentration = [1.0]
        MEKi_concentration = [0.0]
        t_pretrt = [0.83, 0.25, 0.5, 1 , 2, 4, 8]
        for rafi in RAFi_concentration:
            for meki in MEKi_concentration:
                for pretrt in t_pretrt:
                    for trt in t_trt:
                        for egfc in EGF_concentration:
                            settings_list.append([param, pretrt, rafi, meki, trt, egfc]) 

    param_prev = -1
    #run a simulation for each selected condition
    for iset in range(len(settings_list)):
   
        #unload the settings and prepare the unperturbed model
        [param, pretrt, rafi, meki, trt, egfc] =  settings_list[iset]
    
        #run a simulation of the unperturbed model to obtain initial steady state (if not run before)
        if not (param == param_prev):
           param_sets = pd.read_csv(params_file, index_col=0)
           param_sets = param_sets.drop('chi2', axis=1)
           params = param_sets.iloc[param].to_dict()
           ScipyOdeSimulator._use_inline = True
           sim = ScipyOdeSimulator(model, param_values=params,  atol=1E-50)
           df_eq = equilibrate(sim, None)
        

        #run a time-course simulation for the pretreatment phase
        RAFi_index = get_species_index(model, model.monomers.RAFi(raf=None))
        MEKi_index = get_species_index(model, model.monomers.MEKi(mek=None))
        EGF_index = get_species_index(model, model.monomers.EGF(rtk=None))
        initials_pre = df_eq.iloc[-1, :len(model.species)].copy()
        initials_pre[RAFi_index] = rafi
        initials_pre[MEKi_index] = meki
        initials_pre[EGF_index] = 0.0
        tspan_pretrt = np.linspace(0, pretrt, N_time_points)
        df_pre= sim.run(tspan=tspan_pretrt, initials=initials_pre.to_list()).dataframe
        df_pre['time'] = df_pre.index
        df_pre['time'] = df_pre['time']-pretrt
        df_pre['time'].iloc[-1] = 0
        df_pre.reset_index(drop=True, inplace=True)
        df_pre.set_index('time', inplace=True)
    
        #run a time-course simulation for the EGF perturbaion phase
        tspan_trt = np.linspace(0, trt, N_time_points)
        initials_trt = df_pre.iloc[-1, :len(model.species)].copy()
        initials_trt[RAFi_index] = rafi
        initials_trt[MEKi_index] = meki
        initials_trt[EGF_index] = egfc / model.expressions['m_Da_EGF'].get_value()
        df_trt = sim.run(tspan=tspan_trt, initials=initials_trt.to_list()).dataframe
    
        #concatenate pretreatment and EGFR perturnations and settings 
        #obs_names = [x.name for x in model.observables]
        #obs = pd.concat([df_pre, res_trt.dataframe])[obs_names]
        # obs = pd.concat([df_pre, df_trt.iloc[1:]])[df_pre.keys()[1007:-1]]
        obs = pd.concat([df_pre, df_trt.iloc[1:]])[df_pre.keys()[len(model.species):]]
        obs.loc[:, (obs < 1e-10).all()] = 0
   
        #append simulations results on file 
        #df_settings=pd.DataFrame(obs.shape[0]*[[param, pretrt, rafi, meki, trt, egfc]], columns=['param_set_index', 'time_pre_treatment' , 'RAFi_concentration' , 'MEKi_concentration' , 'time_treatment', 'EGF_concentration'] )  
        df_settings=pd.DataFrame(obs.shape[0]*[['A375_sim', param, 'Vemurafenib', 'Cobimetinib', rafi, meki, pretrt, pretrt, egfc, trt]], columns=['Cell_line', 'Parameter_set', 'Drug A', 'Drug B', 'Concentration A (uM)', 'Concentration B (uM)', 'Time A (h)', 'Time B (h)', 'EGF (ng/mL)', 'EGF total duration (h)'] )    
        obs = obs.join(df_settings.set_index(obs.index))
        #obs = pd.concat([df_settings,obs.index], axis=1)
        if (iset == 0): 
           obs.to_csv(output_file, mode='w', header=True)
        else:
           obs.to_csv(output_file, mode='a', header=False)
        
        #update param used
        param_prev=param

if __name__ == "__main__":
    run_MARM1_simulation(sys.argv[1],int(sys.argv[2]))


