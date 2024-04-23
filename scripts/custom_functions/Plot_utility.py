import numpy as np
def generate_heatmap_pixel_map(log_conc_values):
    """gives mapping between pixel location in imag and concentration of drug"""
    values = np.sort(np.unique(log_conc_values))
    count = len(values)
    m = np.array([[values[1],1],[values[-1],1]])
    b = np.array([1.5/count,1-1/(count*2)])
    lin_tf = np.linalg.solve(m,b)
    return lin_tf

def signif(x, p):
    x = np.asarray(x)
    x_positive = np.where(np.isfinite(x) & (x != 0), np.abs(x), 10**(p-1))
    mags = 10 ** (p - 1 - np.floor(np.log10(x_positive)))
    return np.round(x * mags) / mags