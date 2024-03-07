
#Simulate a PySB model from given initial conditions until it reaches steady state
import numpy as np
import pandas as pd

def equilibrate(simulator, initials):
    """Simulate a model from given initial conditions until it reaches steady state"""
    scale = 10
    t_start = 1e-4
    df = None
    tspan = np.geomspace(t_start, t_start * scale)
    while True:
        #print(f"    at t={tspan[-1]:<5.3g} ... ", end='', flush=True)
        res = simulator.run(tspan=tspan, initials=initials)
        df = pd.concat([df, res.dataframe.iloc[1:]])
        initials = res.species[-1]
        close = np.isclose(
            *res.species[[-1,-2]].view(float).reshape(2,-1),
            rtol=1e-3
        )
        cs = np.sum(close)
        n = len(simulator.model.species)
        ##print(f"{cs}/{n} species converged")
        if np.all(close):
            break
        tspan *= scale
    return df
