#code to generate a template for the DDoses file

#set working directory
cwd <- "/gstore/home/gerosal/projects/work/panRAFi_MEKi_combo"
setwd(cwd)

#set directory names
data_dir = 'data'
dda_dir = 'Dose_projections'

#### CREATE FuDrugs datasets for RAF/MEK inhibitors
FuDrugs = data.frame(DrugNamePlot = character(),
                     type = character(),
                     fu_FBS  = numeric(), 
                     min_fu_FBS  = numeric(),
                     max_fu_FBS  = numeric(),
                     FBS_perc  = numeric()
                     )
#Belva
FuDrugs[nrow(FuDrugs)+1,] <- c('panRAFi_Belvarafenib', 'Cell_culture', 0.034, 0.034, 0.034, 10) #measured, Michael Dolton email
FuDrugs[nrow(FuDrugs)+1,] <- c('panRAFi_Belvarafenib', 'Cell_culture', 0.068, 0.068, 0.068, 5) #assumed, Michael Dolton email
#Cobi
FuDrugs[nrow(FuDrugs)+1,] <- c('MEKi_Cobimetinib', 'Cell_culture',  0.1955,  0.1955, 0.1955, 10) #measured, Michael Dolton email
FuDrugs[nrow(FuDrugs)+1,] <- c('MEKi_Cobimetinib', 'Cell_culture', 0.3,  0.3, 0.3, 5) #assumed, Michael Dolton email

#manual conversion to numeric 
FuDrugs[,  c(3:6)] <- sapply(FuDrugs[, c(3:6)], as.numeric)
write.csv(FuDrugs, file.path(cwd, data_dir, dda_dir, 'FuDrugs.csv'), quote=F, row.names=FALSE)


##### COBIMETINIB plasma free drug calculations ####
#From Michal Dolton's email:
#Steady-state trough total concentration (median, 10th, 90th): 16.7 ng/mL (3.53 - 53.9 ng/mL)
#Cobimetinib Steady-state trough free concentration (median, 10th, 90th): 0.87 ng/mL (0.18 - 2.80 ng/mL)
#Cobi MW: 531.3 g/mol
#other doses from https://www.accessdata.fda.gov/drugsatfda_docs/nda/2015/206192orig1s000clinpharmr.pdf
#table 7
drug_ng_ml_cnames <- c('value','min','max')
Cobi_ng_ml <- c( c(16.7, 3.53, 53.9) , #20_mg_TIW
                 c(37, 10^(log10(37)-(log10(58)-log10(37))), 58), #20_mg_QD average = AUC 886/24 h
                 c(160, 10^(log10(160)-(log10(272)-log10(160))), 272), #40_mg_QD average = AUC 3840/24 h
                 c(233, 10^(log10(233)-(log10(364)-log10(233))), 364) #60_mg_QD average = AUC 5600/24 h
                 #c(335, 10^(log10(335)-(log10(525)-log10(335))), 525) #80_mg_QD average = AUC 8060/24 h
                )
Cobi_rnames <- c('20mg_TIW', '20mg_QD', '40mg_QD', '60mg_QD') #, '80mg_QD')
#calculate free conc
fu_human_Cobi <- 0.052
Cobi_MW <- 531.3
Cobi_uM_free_conc <- fu_human_Cobi * (Cobi_ng_ml / Cobi_MW) # (ng/mL / g/mol) = 10^-6 mol/L = uMol
#create matrix
Cobi_uM_free_conc <- t(matrix(Cobi_uM_free_conc, nrow=3, ncol=length(Cobi_rnames)))
rownames(Cobi_uM_free_conc) <- Cobi_rnames
colnames(Cobi_uM_free_conc) <- drug_ng_ml_cnames

### Belvarafenib plasma free drug calculations ####
# From Nature paper: Supp Fig. 5b
# from Michale Dolton's email: 200 mg - 300 mg BID we are considering
Belva_ng_ml <- c(c(789, 10^(log10(789)-(log10(1115)-log10(789))), 1115),  # 50mg_QD  AUC 18931 /24
                 c(1748, 10^(log10(1748)-(log10(1986)-log10(1748))), 1986), # 100mg_QD  AUC 41952 /24
                 c(2009, 10^(log10(2009)-(log10(2213)-log10(2009))), 2213), # 200mg_QD AUC 48215 /24
                 c(2871, 10^(log10(2871)-(log10(3213)-log10(2871))), 3213) # 400mg_BID AUC 68914 / 24
)
Belva_rnames <- c('50mg_QD', '100mg_QD', '200mg_QD', '400mg_BID')
#calculate free conc
fu_human_Belva <- 0.00258
Belva_MW <- 478.93
Belva_uM_free_conc <- fu_human_Belva * (Belva_ng_ml / Belva_MW) # (ng/mL / g/mol) = 10^-6 mol/L = uMol
#create matrix
Belva_uM_free_conc <- t(matrix(Belva_uM_free_conc, nrow=3, ncol=length(Belva_rnames)))
rownames(Belva_uM_free_conc) <- Belva_rnames
colnames(Belva_uM_free_conc) <- drug_ng_ml_cnames


#### CREATE DOSES  ####
DDoses = data.frame(DrugNamePlot = character(), 
                    dd = numeric(), 
                    dd_min = numeric(),
                    dd_max = numeric(),
                    title = character())
for (i in 1:nrow(Belva_uM_free_conc)){
  DDoses[nrow(DDoses)+1,] <- c('panRAFi_Belvarafenib', Belva_uM_free_conc[i,1], Belva_uM_free_conc[i,2], Belva_uM_free_conc[i,3], rownames(Belva_uM_free_conc)[i])
}
for (i in 1:nrow(Cobi_uM_free_conc)){
  DDoses[nrow(DDoses)+1,] <- c('MEKi_Cobimetinib', Cobi_uM_free_conc[i,1], Cobi_uM_free_conc[i,2], Cobi_uM_free_conc[i,3], rownames(Cobi_uM_free_conc)[i])
}

#manual conversion to numeric 
DDoses[, c(2:4)] <- sapply(DDoses[, c(2:4)], as.numeric)
write.csv(DDoses, file.path(cwd, data_dir, dda_dir, 'DDoses.csv'), quote=F, row.names=FALSE)
