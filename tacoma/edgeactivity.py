# -*- coding: utf-8 -*-
"""
This module provides functions related to the edge activity
temporal network model.
"""
from copy import deepcopy

import numpy as np
from numpy import random

from scipy.integrate import ode
from scipy.integrate import simps
from scipy.special import gamma as Gamma

import tacoma as tc

from tacoma.flockwork import convert_to_varying_edge_activity_parameters


def naive_varying_rate_edge_activity_simulation(N, t, network_densities, activity_rates, tmax):
    r"""
    Do a naive simulation of edge activity model systems where the network density and edge activity rates
    vary over time as step functions. It is called `naive` because the rate change will not
    be considered when evaluating the inter-event time at time points when the rates change.
    I.e. at a rate change at time `t`, the new model will be initiated as if the 
    last event happened at time `t`.

    Parameters
    ----------
    N : int
        number of nodes
    t : numpy.ndarray of float
        time points at which :math:`\alpha` and :math`\beta` change
    network_densities : numpy.ndarray of float
        values of :math:`\rho` for the corresponding time values in ``t``.
    activity_rates : numpy.ndarray of float
        values of :math:`\omega` for the corresponding time values in ``t``.
    tmax : float
        time at which the simulation is supposed to end.

    Returns
    -------
    edge_changes : :class:`_tacoma.edge_changes`
        The simulated edge activity model instance.
    """

    if len(t) != len(network_densities) or len(activity_rates) != len(network_densities):
        raise ValueError('The arrays `t`, `activity_rates` and `network_densities` must have the same length.' )

    assert(t[-1]<tmax)
    t = np.append(t,tmax)

    all_networks = []
    it = 0
    initial_edges = []

    result = None
    for r_, w_ in zip(network_densities, activity_rates):

        dt = t[it+1] - t[it]

        if r_ == 1.0:
            r_ = (0.5*(N-1)*N-1)/(0.5*(N-1)*N)
        elif r_ == 0.0:
            r_ = 1e-100

        #print("rho =", r_, "omega =", w_)

        if w_ > 0.0:

            EAM = tc.EdgeActivityModel( N, r_, w_, save_temporal_network = True)
            if it > 0:
                EAM.set_initial_edgelist(0.0, initial_edges)
            tc.simulate_EdgeActivityModel(EAM, dt, reset_simulation_objects = False)
            this_result = EAM.edge_changes
            if it < len(activity_rates) - 1:
                initial_edges = EAM.get_current_edgelist()

        else:
            this_result = tc.edge_changes()
            this_result.edges_initial = initial_edges
            this_result.N = N
            this_result.t0 = 0.0
            this_result.tmax = dt

        all_networks.append(this_result)

        it += 1

    result = tc.concatenate(all_networks)

    return result


