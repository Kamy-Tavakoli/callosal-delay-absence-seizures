import numpy as np
from numba import njit

@njit
def sigmoid(u, beta):
    return 1 / (1 + np.exp(-beta * u))


@njit
def run_simulation(
    num_delays, cv, T, N, wcc, beta, a_e, a_i, a_th, a_rtn,
    i_e, i_i, i_th, i_rtn, wee, wei, wie, wii, weth, wthi, wthe, wrtnth, wthrtn, wertn,
    tau, noise_intensity, dt, seed
):
    noise_e = np.zeros((N, T))
    noise_i = np.zeros((N, T))
    noise_th = np.zeros((N, T))
    noise_rtn = np.zeros((N, T))
    np.random.seed(seed)
    for n in range(N):
        noise_e[n] = np.random.normal(0, 1, T)
        noise_i[n] = np.random.normal(0, 1, T)
        noise_th[n] = np.random.normal(0, 1, T)
        noise_rtn[n] = np.random.normal(0, 1, T)

    D_e = 1e-5 * noise_intensity
    D_i = 1e-5 * noise_intensity
    D_th = 1e-7 * noise_intensity
    D_rtn = 1e-7 * noise_intensity

    e = np.zeros((N, T))
    i_pop = np.zeros((N, T))
    th = np.zeros((N, T))
    rtn = np.zeros((N, T))
    callosal = np.zeros((N, T))


    max_delay = int(round(2.5 /dt))
    tau_callosal = np.zeros(num_delays, dtype=np.int64)
    for k in range(num_delays):
        while True:
            tau_old = np.round(4.5 / np.abs(cv + np.random.normal() * 0.2))
            tau_new = int(round(tau_old*0.1 / dt))
            if tau_new < max_delay:
                tau_callosal[k] = max(1, tau_new)
                break



    #Euler-Maruyama integration
    for t in range(1, T - 1):

        #delayed callosal connections
        for u in range(num_delays):
            idx = t - tau_callosal[u]
            if idx >= 0:
                callosal[0, t] += wcc * sigmoid(e[1, idx], beta) / num_delays
                callosal[1, t] += wcc * sigmoid(e[0, idx], beta) / num_delays

        for n in range(N):
            idx_delay = t - tau[n]
            e_delay = e[n, idx_delay] if idx_delay >= 0 else 0.0
            th_delay = th[n, idx_delay] if idx_delay >= 0 else 0.0

            e[n, t + 1] = (
                e[n, t]
                + dt * a_e * (
                    -e[n, t]
                    + wee * sigmoid(e[n, t], beta)
                    + wie * sigmoid(i_pop[n, t], beta)
                    + wthe * sigmoid(th_delay, beta)
                    + callosal[n, t]
                    + i_e
                )
                + np.sqrt(2 * D_e * dt) * noise_e[n, t]
            )

            i_pop[n, t + 1] = (
                i_pop[n, t]
                + dt * a_i * (
                    -i_pop[n, t]
                    + wei * sigmoid(e[n, t], beta)
                    + wii * sigmoid(i_pop[n, t], beta)
                    + wthi * sigmoid(th_delay, beta)
                    + i_i
                )
                + np.sqrt(2 * D_i * dt) * noise_i[n, t]
            )

            th[n, t + 1] = (
                th[n, t]
                + dt * a_th * (
                    -th[n, t]
                    + weth * sigmoid(e_delay, beta)
                    + wrtnth * sigmoid(rtn[n, t], beta)
                    + i_th
                )
                + np.sqrt(2 * D_th * dt) * noise_th[n, t]
            )

            rtn[n, t + 1] = (
                rtn[n, t]
                + dt * a_rtn * (
                    -rtn[n, t]
                    + wthrtn * sigmoid(th[n, t], beta)
                    + wertn * sigmoid(e_delay, beta)
                    + i_rtn
                )
                + np.sqrt(2 * D_rtn * dt) * noise_rtn[n, t]
            )

    return e