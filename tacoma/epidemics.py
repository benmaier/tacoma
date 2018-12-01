# -*- coding: utf-8 -*-
"""
This module provides functions related to the simulation and
measurement of epidemics.
"""

import numpy as np
import tacoma as tc


def measure_i_inf(temporal_network_or_model,epidemic_object,t_equilibrate,is_static=False,verbose=False):
    """Get the equilibrium ratio of infected. 

    Parameters
    ----------
    N : int
        Number of nodes
    P : float
        Probability to reconnect

    Returns
    -------
    numpy.ndarray
        Group size distribution of this configuration. The m-th entry
        of this array contains the expected number of nodes of groups
        of size m.
    """

    tn = temporal_network_or_model
    eo = epidemic_object

    N = tn.N
    t_eq = t_equilibrate
    t_run = eo.t_simulation

    t_run_total = t_eq + t_run

    tc.gillespie_epidemics(tn,eo,is_static=is_static,verbose=verbose)

    t = np.array(eo.time)
    I = np.array(eo.I,dtype=float) / N
    r0 = np.array(eo.R0)
    t0 = t[0]

    if t[-1]>t0+t_eq:
        ndcs = np.where(t>=t_eq+t0)[0]
        ti, i = t[ndcs], I[ndcs]
        i_inf = tc.time_average(ti,i,tmax=t0+t_run_total)
        i_inf_std = np.sqrt(tc.time_average(ti,(i-i_inf)**2,tmax=t0+t_run_total))

        r0 = r0[t>=t_eq+t0]
        this_t = t[t>=t_eq+t0]
        R0 = tc.time_average(this_t,r0,tmax=t0+t_run_total)
    else:
        i_inf = 0
        i_inf_std = 0.
        R0 = r0[-1]

    result = ( i_inf, i_inf_std, R0 )

    return result


if __name__ == "__main__":

    N = 5
    k = 2
    omega = 1.6
    recovery_rate = 0.1
    R0 = 10
    t_run_total = 1000

    AM = tc.ActivityModel(N,
                           k/(N-1.),
                           omega,
                           t0 = 2000
                          )
    infection_rate = R0 / k * recovery_rate


    SIS = tc.SIS(N,t_run_total,infection_rate,recovery_rate,
            number_of_initially_infected=N,
            sampling_dt=0.0,
            )

    print(measure_i_inf(AM, SIS,t_equilibrate=900))
