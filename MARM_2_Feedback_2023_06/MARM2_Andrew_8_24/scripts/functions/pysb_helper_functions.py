import numpy as np
import pandas as pd
import itertools
from pysb.core import as_complex_pattern

def equilibrate(simulator, initials,verbose = True):
    """Simulate a model from given initial conditions until it reaches steady state"""
    scale = 10
    t_start = 10
    df = None
    tspan = np.geomspace(t_start, t_start * scale)
    while True:
        if verbose:
            print(f"    at t={tspan[-1]:<5.3g} ... ", end='', flush=True)
        res = simulator.run(tspan=tspan, initials=initials)
        df = pd.concat([df, res.dataframe.iloc[1:]])
        initials = res.species[-1]
        close = np.isclose(
            *res.species[[-1,-2]].view(float).reshape(2,-1),
            rtol=1e-3
        )
        cs = np.sum(close)
        n = len(simulator.model.species)
        if verbose:
            print(f"{cs}/{n} species converged")
        if np.all(close):
            break
        tspan *= scale
    return df

def get_species_index(model, pattern):
    """Return the integer species number for a given species in the model"""
    pattern = as_complex_pattern(pattern)
    matches = [
        i for i, s in enumerate(model.species)
        if s.is_equivalent_to(pattern)
    ]
    n = len(matches)
    assert n == 1, f"Expected exactly one match, got {n}"
    return matches[0]