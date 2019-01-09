# -*- coding: utf-8 -*-
"""
This module provides functions related to the simulation and
measurement of epidemics.
"""

import numpy as np
import tacoma as tc

from _tacoma import gillespie_QS_SIS_on_EdgeActivityModel, gillespie_QS_SIS_on_edge_lists


def simulate_and_measure_i_inf(temporal_network_or_model,epidemic_object,t_equilibrate,is_static=False,verbose=False):
    """Get the equilibrium ratio of infected. 

    Parameters
    ----------
    temporal_network_or_model : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, :class:`_tacoma.edge_lists_with_histograms`, or :class:`_tacoma.EdgeActivityModel`.
        An instance of a temporal network or network model.
    epidemic_object : :class:`_tacoma.SI`, :class:`_tacoma.SIS`, :class:`_tacoma.SIR`, :class:`_tacoma.SIRS`, :class:`_tacoma.node_based_SIS`
        An initialized epidemic object.
    t_equilibrate: float
        Time passed after t0 after which to start measuring.
    is_static : bool, default : False
        The algorithm works a bit differently if it knows that the network is actually static.
        It works only with instances of :class:`_tacoma.edge_lists`.
    verbose: bool, optional
        Be chatty.


    Returns
    -------
    i_inf: float
        Temporal average over the ratio of infected after equilibration.
    i_inf_std: float
        RMSE of the ratio of infected.
    R0: float
        As measured after equilibration
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

def simulate_quasi_stationary_SIS_on_model(model, qs_sis, verbose=False):

    while True:

        gillespie_QS_SIS_on_EdgeActivityModel(model, qs_sis, verbose=verbose)

        if qs_sis.ended_in_absorbing_state():
            node_status, G = qs_sis.get_random_configuration()
            qs_sis.set_node_configuration(qs_sis.last_active_time, node_status)
            model.set_initial_configuration(qs_sis.last_active_time, G)
        else:
            break

    return qs_sis.get_infection_observables()
        

def simulate_quasi_stationary_SIS_on_static_network(network, qs_sis, verbose=False):

    t = float(qs_sis.last_active_time)

    while True:

        gillespie_QS_SIS_on_edge_lists(network, qs_sis, is_static=True, verbose=verbose)

        if qs_sis.ended_in_absorbing_state():

            delta_t = qs_sis.last_active_time - t

            node_status, G = qs_sis.get_random_configuration()
            qs_sis.set_initial_configuration(qs_sis.last_active_time, node_status)


            qs_sis.t_simulation -= delta_t
            #print(qs_sis.last_active_time, t, qs_sis.t_simulation)
            #sys.exit(1)
            new_t = [qs_sis.last_active_time]
            network.t = new_t
            t = qs_sis.last_active_time
            print("t = ", t,": died")
        else:
            break

    return qs_sis.get_infection_observables()

def get_SIS_critical_infection_rate(tn, recovery_rate):

    pass

def get_SIS_critical_recovery_rate(tn, infection_rate):
    pass

"""
def get_SIS_epidemic_threshold(tn, infection_rate=None, recovery_rate=None):

    if infection_rate is None and recovery_rate is None:
        raise ValueError('Please provide either an infection rate (to find the critical recovery rate) (x)or a recovery rate (to find the critical infection rate).')
    elif infection_rate is not None and recovery_rate is not None:
        raise ValueError('Please provide either an infection rate (to find the critical recovery rate) (x)or a recovery rate (to find the critical infection rate), not both.')
"""

    

if __name__ == "__main__":

    N = 100
    k = 10
    omega = 1.6
    recovery_rate = 0.1
    R0 = 10
    t_run_total = 1000

    AM = tc.EdgeActivityModel(N,
                           k/(N-1.),
                           omega,
                           t0 = 2000
                          )
    infection_rate = R0 / k * recovery_rate


    SIS = tc.SIS(N,t_run_total,infection_rate,recovery_rate,
            number_of_initially_infected=N,
            sampling_dt=0.0,
            )

    print(simulate_and_measure_i_inf(AM, SIS,t_equilibrate=900))
