import numpy as np

def threshold(I, T_min, T_max):
    I_thr = np.interp(I, (T_min, T_max), (0, np.iinfo(I.dtype).max)).astype(I.dtype)
    I_thr[I < T_min] = 0
    return I_thr