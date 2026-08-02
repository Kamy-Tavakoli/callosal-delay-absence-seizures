import numpy as np

def get_parameters():

    N = 2
    dt = 0.025

    T = (40000000 + 24000) 
    beta = 25
    tau = np.array([0.5 /dt, 0.5 /dt], dtype=np.int64)

    a_e = 0.7
    a_i = 0.7
    a_th = 0.07
    a_rtn = 0.035

    i_e = -0.5
    i_i = -0.61
    i_th = 0.16
    i_rtn = -0.514

    wee = 0.9
    wei = 1.5
    wie = -2.2
    wii = -0.5
    weth = 1.9
    wthi = 1.8
    wthe = 1.5
    wrtnth = -1.6
    wthrtn = 1.44
    wertn = 2

    return (
        T, N, beta, a_e, a_i, a_th, a_rtn,
        i_e, i_i, i_th, i_rtn,
        wee, wei, wie, wii,
        weth, wthi, wthe, wrtnth, wthrtn, wertn,
        tau, dt
    )
