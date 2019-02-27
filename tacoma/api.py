# -*- coding: utf-8 -*-
"""
This module provides wrapper functions for the core functions provided in the
C++-module `_tacoma`.
"""

import numpy as np

import _tacoma as _tc

from _tacoma import edge_changes as ec
from _tacoma import edge_lists as el
from _tacoma import edge_lists_with_histograms as el_h
from _tacoma import edge_changes_with_histograms as ec_h

edge_changes = _tc.edge_changes
edge_lists = _tc.edge_lists
edge_lists_with_histograms = _tc.edge_lists_with_histograms
edge_changes_with_histograms = _tc.edge_changes_with_histograms
edge_trajectories = _tc.edge_trajectories
edge_trajectory_entry = _tc.edge_trajectory_entry

from tacoma.network_formats import adjacency_matrices
from tacoma.network_formats import sparse_adjacency_matrices

edge_weight = _tc.edge_weight
social_trajectory_entry = _tc.social_trajectory_entry

SI = _tc.SI
SIS = _tc.SIS
SIR = _tc.SIR
SIRS = _tc.SIRS
node_based_SIS = _tc.node_based_SIS
eSIS = _tc.eSIS
EdgeActivityModel = _tc.EdgeActivityModel
QS_SIS = _tc.QS_SIS
coverage_SIS = _tc.coverage_SIS
cluster_size_SIS = _tc.cluster_size_SIS
FlockworkPModel = _tc.FlockworkPModel
MARKOV_SIS = _tc.MARKOV_SIS


def _get_raw_temporal_network(temporal_network):
    """Return an instance of `edge_changes` if the temporal network was an instance of 
    `edge_changes_with_histograms` (return `edge_lists` for `edge_lists_with_histograms`)."""

    if type(temporal_network) == ec_h:
        temporal_network = ec(temporal_network)
    elif type(temporal_network) == el_h:
        temporal_network = el(temporal_network)

    return temporal_network


def measure_group_sizes_and_durations(temporal_network, ignore_size_histogram=False, verbose=False):
    """Measures aggregated group size distribution, group size histogram for each time point, 
    contact durations, inter-contact durations, durations of each group size, and the aggregated social network.

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        An instance of a temporal network.
    ignore_size_histogram : bool, optional
        Don't compute the single time point group size histograms (save time and memory). default : False
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    :class:`_tacoma.group_sizes_and_durations`
        The result of the measurements.
        Note that `inter-contact durations` is just durations of groups of size 1 and hence correspond to the
        first entry of :class:`_tacoma.group_sizes_and_durations.group_durations`.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = _tc.measure_group_sizes_and_durations_for_edge_changes(temporal_network,
                                                                        ignore_size_histogram,
                                                                        verbose,
                                                                        )
    elif type(temporal_network) == el:
        result = _tc.measure_group_sizes_and_durations_for_edge_lists(temporal_network,
                                                                      ignore_size_histogram,
                                                                      verbose,
                                                                      )
    else:
        raise ValueError('Unknown temporal network format: ' +
                         str(type(temporal_network)))

    return result


def aggregated_network(temporal_network):
    """Measures the static network composed.

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        An instance of a temporal network.

    Returns 
    -------
    aggregated_network : :obj:`dict` where keys are pairs of ints and values are floats
        Keys represent edges and values represent the aggregated time spent together (in units
        of time provided by the temporal network

    """

    return measure_group_sizes_and_durations(temporal_network).aggregated_network


def degree_distribution(temporal_network):
    """Measures the time-averaged degree distribution.

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        An instance of a temporal network.

    Returns 
    -------
    P_k : :obj:`list` of `float`
        A list of N entries where the k-th entry is the time-averaged probability of finding a node with degree k.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = _tc.degree_distribution_from_edge_changes(temporal_network)
    elif type(temporal_network) == el:
        result = _tc.degree_distribution_from_edge_lists(temporal_network)
    else:
        raise ValueError('Unknown temporal network format: ' +
                         str(type(temporal_network)))

    return result


def bin(temporal_network, dt=0.0, N_time_steps=0, verbose=False):
    """Bins a temporal network for each `dt` (after each step, respectively, if `N_time_steps` was provided).

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        An instance of a temporal network.
    dt : float
        The demanded bin size. default : 0.0
    N_time_steps : int
        Number of time bins (use either this or dt). default : 0
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    :class:`_tacoma.edge_lists` 
        An edge_lists instance with one edge list describing the network as a list of all edges
        that were present in the last time bin.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = _tc.bin_from_edge_changes(
            temporal_network, dt, N_time_steps, verbose)
    elif type(temporal_network) == el:
        result = _tc.bin_from_edge_lists(
            temporal_network, dt, N_time_steps, verbose)
    else:
        raise ValueError('Unknown temporal network format: ' +
                         str(type(temporal_network)))

    return result


def slice(temporal_network, new_t0, new_tmax, verbose=False):
    """Get a slice of a temporal network for times new_t0 <= t < new_tmax.

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        An instance of a temporal network.
    new_t0 : float
        Where the slice should begin.
    new_t0 : float
        Where the slice should end.
    verbose: bool, optional
        Be chatty.
        default : False

    Returns
    -------
    :class:`_tacoma.edge_lists` or :class:`_tacoma.edge_changes` 
        Sliced network (type depending on input type).
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = _tc.slice_edge_changes(
            temporal_network, new_t0, new_tmax, verbose)
    elif type(temporal_network) == el:
        result = _tc.slice_edge_lists(
            temporal_network, new_t0, new_tmax, verbose)
    else:
        raise ValueError('Unknown temporal network format: ' +
                         str(type(temporal_network)))

    return result


def sample(temporal_network, dt=0.0, N_time_steps=0, sample_aggregates=False, verbose=False):
    """Samples a temporal network after each `dt` (after each step, respectively, if `N_time_steps` was provided).

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        An instance of a temporal network.
    dt : float
        The demanded bin size. default : 0.0
    N_time_steps : int
        Number of time bins (use either this or dt). default : 0
    sample_aggregates : bool, optional
        If this is True, the following happens. If an edge is active during a time bin,
        it will appear in the final graph at the end of this time bin. It will then exist
        until the ende of the coming time bin. (This is different from the binning procedure). 
        This does not seem to make very much sense so maybe just don't use it.
        default : False
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    :class:`_tacoma.edge_lists`
        An edge_lists instance with one edge list describing the network states after each
        discrete time step.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if isinstance(temporal_network, ec):
        result = _tc.sample_from_edge_changes(temporal_network,
                                              dt,
                                              N_time_steps,
                                              sample_aggregates,
                                              verbose,
                                              )
    elif isinstance(temporal_network, el):
        result = _tc.sample_from_edge_lists(temporal_network,
                                            dt,
                                            N_time_steps,
                                            sample_aggregates,
                                            verbose,
                                            )
    else:
        raise ValueError('Unknown temporal network format: ' +
                         str(type(temporal_network)))

    return result


def binned_social_trajectory(temporal_network, node, dt=0.0, N_time_steps=0, verbose=False):
    """Computes the binned social trajectory of node `node` in `temporal_network`, i.e. the groups it was part of and pairs of times (t0, t1) when it was part of a certain group. In this case, binning makes it easier to, e.g. see the days where a group was active.

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        An instance of a temporal network.
    node : int
        The node for which to compute the social trajectory
    dt : float
        The demanded bin size. default : 0.0
    N_time_steps : int
        Number of time bins (use either this or dt). default : 0
    verbose: bool, optional
        Be chatty.

    Returns
    -------
        :obj:`list` of :obj:`list` of int
            A list, one entry for each dt. Each entry is a list of group ids.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = _tc.binned_social_trajectory_from_edge_changes(temporal_network,
                                                                node,
                                                                dt,
                                                                N_time_steps,
                                                                verbose
                                                                )
    elif type(temporal_network) == el:
        result = _tc.binned_social_trajectory_from_edge_lists(temporal_network,
                                                              node,
                                                              dt,
                                                              N_time_steps,
                                                              verbose
                                                              )
    else:
        raise ValueError('Unknown temporal network format: ' +
                         str(type(temporal_network)))

    return result


def social_trajectory(temporal_network, node, verbose=False):
    """Computes the social trajectory of node `node` in `temporal_network`, i.e. the groups it was part of and pairs of times (t0, t1) when it was part of a certain group.

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        An instance of a temporal network.
    node : int
        The node for which to compute the social trajectory
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    :obj:`list` of :class:`_tacoma.social_trajectory_entry`
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = _tc.social_trajectory_from_edge_changes(
            temporal_network, node, verbose)
    elif type(temporal_network) == el:
        result = _tc.social_trajectory_from_edge_lists(
            temporal_network, node, verbose)
    else:
        raise ValueError('Unknown temporal network format: ' +
                         str(type(temporal_network)))

    return result

def markov_epidemics(temporal_network_or_model, epidemic_object, max_dt, is_static=False, reset_simulation_objects=True, verbose=False):
    """Integrates an epidemic markov process on the provided temporal network.

    Parameters
    ----------
    temporal_network_or_model : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, :class:`_tacoma.edge_lists_with_histograms`, :class:`_tacoma.EdgeActivityModel`, or :class:`_tacoma.FlockworkPModel`.
        An instance of a temporal network or network model.
    epidemic_object : :class:`_tacoma.MARKOV_SIS`
        An initialized epidemic object.
    max_dt : float
        The maximally allowed integration step.
    is_static : bool, default : False
        The algorithm works a bit differently if it knows that the network is actually static.
        It works only with instances of :class:`_tacoma.edge_lists`.
    reset_simulation_objects : bool, default : True
        Whether or not the simulation objects should be reset in the beginning.
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    None
        But the observables are saved in the :mod:`_tacoma` epidemic Markov-object.
    """

    tn_or_mdl = temporal_network_or_model

    if type(tn_or_mdl) in [ec, ec_h]:
        temporal_network = _get_raw_temporal_network(tn_or_mdl)
        if type(epidemic_object) == MARKOV_SIS:
            _tc.markov_SIS_on_edge_changes(temporal_network, epidemic_object, max_dt, verbose)
        else:
            raise ValueError('Invalid epidemic object type: ' + str(type(epidemic_object)))
    elif type(tn_or_mdl) in [el, el_h]:
        temporal_network = _get_raw_temporal_network(tn_or_mdl)
        if type(epidemic_object) == MARKOV_SIS:
            _tc.markov_SIS_on_edge_lists(temporal_network, epidemic_object, max_dt, is_static, verbose)
        else:
            raise ValueError('Invalid epidemic object type: ' + str(type(epidemic_object)))
    elif type(tn_or_mdl) in [EdgeActivityModel]:
        temporal_network = tn_or_mdl
        if type(epidemic_object) == MARKOV_SIS:
            _tc.markov_SIS_on_EdgeActivityModel(temporal_network, epidemic_object, max_dt, reset_simulation_objects, verbose)
        else:
            raise ValueError('Invalid epidemic object type: ' + str(type(epidemic_object)))
    elif type(tn_or_mdl) in [FlockworkPModel]:
        temporal_network = tn_or_mdl
        if type(epidemic_object) == MARKOV_SIS:
            _tc.markov_SIS_on_FlockworkPModel(temporal_network, epidemic_object, max_dt, reset_simulation_objects, verbose)
        else:
            raise ValueError('Invalid epidemic object type: ' + str(type(epidemic_object)))
    else:
        raise ValueError('Invalid temporal network/model format: ' + str(type(tn_or_mdl)))

def gillespie_epidemics(temporal_network_or_model, epidemic_object, is_static=False, reset_simulation_objects=True, verbose=False):
    """Simulates an epidemic process on the provided temporal network using the Gillespie stochastic simulation algorithm.

    Parameters
    ----------
    temporal_network_or_model : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, :class:`_tacoma.edge_lists_with_histograms`, :class:`_tacoma.EdgeActivityModel`, or :class:`_tacoma.FlockworkPModel`.
        An instance of a temporal network or network model.
    epidemic_object : :class:`_tacoma.SI`, :class:`_tacoma.SIS`, :class:`_tacoma.SIR`, :class:`_tacoma.SIRS`, :class:`_tacoma.node_based_SIS`, :class:`_tacoma.eSIS`
        An initialized epidemic object.
    is_static : bool, default : False
        The algorithm works a bit differently if it knows that the network is actually static.
        It works only with instances of :class:`_tacoma.edge_lists`.
    reset_simulation_objects : bool, default : True
        Whether or not the simulation objects should be reset in the beginning.
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    None
        But the observables are saved in the :mod:`_tacoma` epidemic object.
    """

    tn_or_mdl = temporal_network_or_model

    if type(tn_or_mdl) in [ec, ec_h]:
        temporal_network = _get_raw_temporal_network(tn_or_mdl)
        if type(epidemic_object) == SI:
            _tc.gillespie_SI_on_edge_changes(temporal_network, epidemic_object, verbose)
        elif type(epidemic_object) == SIS:
            _tc.gillespie_SIS_on_edge_changes(temporal_network, epidemic_object, verbose)
        elif type(epidemic_object) == SIR:
            _tc.gillespie_SIR_on_edge_changes(temporal_network, epidemic_object, verbose)
        elif type(epidemic_object) == SIRS:
            _tc.gillespie_SIRS_on_edge_changes(temporal_network, epidemic_object, verbose)
        elif type(epidemic_object) == node_based_SIS:
            _tc.gillespie_node_based_SIS_on_edge_changes(temporal_network, epidemic_object, verbose)
        elif type(epidemic_object) == eSIS:
            _tc.gillespie_eSIS_on_edge_changes(temporal_network, epidemic_object, verbose)
        elif type(epidemic_object) == coverage_SIS:
            _tc.gillespie_coverage_SIS_on_edge_changes(temporal_network, epidemic_object, verbose)
        elif type(epidemic_object) == cluster_size_SIS:
            _tc.gillespie_cluster_size_SIS_on_edge_changes(temporal_network, epidemic_object, verbose)
        else:
            raise ValueError('Invalid epidemic object type: ' + str(type(epidemic_object)))
    elif type(tn_or_mdl) in [el, el_h]:
        temporal_network = _get_raw_temporal_network(tn_or_mdl)
        if type(epidemic_object) == SI:
            _tc.gillespie_SI_on_edge_lists(temporal_network, epidemic_object, is_static, verbose)
        elif type(epidemic_object) == SIS:
            _tc.gillespie_SIS_on_edge_lists(temporal_network, epidemic_object, is_static, verbose)
        elif type(epidemic_object) == SIR:
            _tc.gillespie_SIR_on_edge_lists(temporal_network, epidemic_object, is_static, verbose)
        elif type(epidemic_object) == SIRS:
            _tc.gillespie_SIRS_on_edge_lists(temporal_network, epidemic_object, is_static, verbose)
        elif type(epidemic_object) == node_based_SIS:
            _tc.gillespie_node_based_SIS_on_edge_lists(temporal_network, epidemic_object, is_static, verbose)
        elif type(epidemic_object) == eSIS:
            _tc.gillespie_eSIS_on_edge_lists(temporal_network, epidemic_object, is_static, verbose)
        elif type(epidemic_object) == coverage_SIS:
            _tc.gillespie_coverage_SIS_on_edge_lists(temporal_network, epidemic_object, is_static, verbose)
        elif type(epidemic_object) == cluster_size_SIS:
            _tc.gillespie_cluster_size_SIS_on_edge_lists(temporal_network, epidemic_object, is_static, verbose)
        else:
            raise ValueError('Invalid epidemic object type: ' + str(type(epidemic_object)))
    elif type(tn_or_mdl) in [EdgeActivityModel]:
        temporal_network = tn_or_mdl
        if type(epidemic_object) == SI:
            _tc.gillespie_SI_on_EdgeActivityModel(temporal_network, epidemic_object, reset_simulation_objects, verbose)
        elif type(epidemic_object) == SIS:
            _tc.gillespie_SIS_on_EdgeActivityModel(temporal_network, epidemic_object, reset_simulation_objects, verbose)
        elif type(epidemic_object) == SIR:
            _tc.gillespie_SIR_on_EdgeActivityModel(temporal_network, epidemic_object, reset_simulation_objects, verbose)
        elif type(epidemic_object) == SIRS:
            _tc.gillespie_SIRS_on_EdgeActivityModel(temporal_network, epidemic_object, reset_simulation_objects, verbose)
        elif type(epidemic_object) == node_based_SIS:
            _tc.gillespie_node_based_SIS_on_EdgeActivityModel(temporal_network, epidemic_object, reset_simulation_objects, verbose)
        elif type(epidemic_object) == eSIS:
            _tc.gillespie_eSIS_on_EdgeActivityModel(temporal_network, epidemic_object, reset_simulation_objects, verbose)
        elif type(epidemic_object) == coverage_SIS:
            _tc.gillespie_coverage_SIS_on_EdgeActivityModel(temporal_network, epidemic_object, reset_simulation_objects, verbose)
        elif type(epidemic_object) == cluster_size_SIS:
            _tc.gillespie_cluster_size_SIS_on_EdgeActivityModel(temporal_network, epidemic_object, reset_simulation_objects, verbose)
        else:
            raise ValueError('Invalid epidemic object type: ' + str(type(epidemic_object)))
    elif type(tn_or_mdl) in [FlockworkPModel]:
        temporal_network = tn_or_mdl
        if type(epidemic_object) == SI:
            _tc.gillespie_SI_on_FlockworkPModel(temporal_network, epidemic_object, reset_simulation_objects, verbose)
        elif type(epidemic_object) == SIS:
            _tc.gillespie_SIS_on_FlockworkPModel(temporal_network, epidemic_object, reset_simulation_objects, verbose)
        elif type(epidemic_object) == SIR:
            _tc.gillespie_SIR_on_FlockworkPModel(temporal_network, epidemic_object, reset_simulation_objects, verbose)
        elif type(epidemic_object) == SIRS:
            _tc.gillespie_SIRS_on_FlockworkPModel(temporal_network, epidemic_object, reset_simulation_objects, verbose)
        elif type(epidemic_object) == node_based_SIS:
            _tc.gillespie_node_based_SIS_on_FlockworkPModel(temporal_network, epidemic_object, reset_simulation_objects, verbose)
        elif type(epidemic_object) == eSIS:
            _tc.gillespie_eSIS_on_FlockworkPModel(temporal_network, epidemic_object, reset_simulation_objects, verbose)
        elif type(epidemic_object) == coverage_SIS:
            _tc.gillespie_coverage_SIS_on_FlockworkPModel(temporal_network, epidemic_object, reset_simulation_objects, verbose)
        elif type(epidemic_object) == cluster_size_SIS:
            _tc.gillespie_cluster_size_SIS_on_FlockworkPModel(temporal_network, epidemic_object, reset_simulation_objects, verbose)
        else:
            raise ValueError('Invalid epidemic object type: ' + str(type(epidemic_object)))
    else:
        raise ValueError('Invalid temporal network/model format: ' + str(type(tn_or_mdl)))

def gillespie_SIS(temporal_network, SIS, is_static=False, verbose=False):
    """Simulates an SIS process on the provided temporal network using the Gillespie stochastic simulation algorithm.

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        An instance of a temporal network.
    SIS : :class:`_tacoma.SIS`
        An initialized SIS object.
    is_static : bool, default : False
        The algorithm works a bit differently if it knows that the network is actually static.
        It works only with instances of :class:`_tacoma.edge_lists`.
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    None
        But the observables are saved in the :class:`_tacoma.SIS` object.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = _tc.gillespie_SIS_on_edge_changes(
            temporal_network, SIS, verbose)
    elif type(temporal_network) == el:
        result = _tc.gillespie_SIS_on_edge_lists(
            temporal_network, SIS, is_static, verbose)
    else:
        raise ValueError('Unknown temporal network format: ' +
                         str(type(temporal_network)))

    return result


def gillespie_node_based_SIS(temporal_network, SIS, is_static=False, verbose=False):
    """Simulates a node-based SIS process on the provided temporal network using the Gillespie stochastic simulation algorithm.

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        An instance of a temporal network.
    SIS : :class:`_tacoma.SIS_node_based`
        An initialized node-based SIS object.
    is_static : bool, default : False
        The algorithm works a bit differently if it knows that the network is actually static.
        It works only with instances of :class:`_tacoma.edge_lists`.
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    None
        But the observables are saved in the :class:`_tacoma.SIS` object.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = _tc.gillespie_node_based_SIS_on_edge_changes(
            temporal_network, SIS, verbose)
    elif type(temporal_network) == el:
        result = _tc.gillespie_node_based_SIS_on_edge_lists(
            temporal_network, SIS, is_static, verbose)
    else:
        raise ValueError('Unknown temporal network format: ' +
                         str(type(temporal_network)))

    return result


def gillespie_SI(temporal_network, SI, is_static=False, verbose=False):
    """Simulates an SI process on the provided temporal network using the Gillespie stochastic simulation algorithm.

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        An instance of a temporal network.
    SI : :class:`_tacoma.SI`
        An initialized SI object.
    is_static : bool, default : False
        The algorithm works a bit differently if it knows that the network is actually static.
        It works only with instances of :class:`_tacoma.edge_lists`.
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    None
        But the observables are saved in the :class:`_tacoma.SI` object.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = _tc.gillespie_SI_on_edge_changes(
            temporal_network, SI, verbose)
    elif type(temporal_network) == el:
        result = _tc.gillespie_SI_on_edge_lists(
            temporal_network, SI, is_static, verbose)
    else:
        raise ValueError('Unknown temporal network format: ' +
                         str(type(temporal_network)))

    return result


def gillespie_SIR(temporal_network, SIR, is_static=False, verbose=False):
    """Simulates an SIR process on the provided temporal network using the Gillespie stochastic simulation algorithm.

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        An instance of a temporal network.
    SIR : :class:`_tacoma.SIR`
        An initialized SIR object.
    is_static : bool, default : False
        The algorithm works a bit differently if it knows that the network is actually static.
        It works only with instances of :class:`_tacoma.edge_lists`.
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    None
        But the observables are saved in the :class:`_tacoma.SIR` object.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = _tc.gillespie_SIR_on_edge_changes(
            temporal_network, SIR, verbose)
    elif type(temporal_network) == el:
        result = _tc.gillespie_SIR_on_edge_lists(
            temporal_network, SIR, is_static, verbose)
    else:
        raise ValueError('Unknown temporal network format: ' +
                         str(type(temporal_network)))

    return result


def gillespie_SIRS(temporal_network, SIRS, is_static=False, verbose=False):
    """Simulates an SIRS process on the provided temporal network using the Gillespie stochastic simulation algorithm.

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        An instance of a temporal network.
    SIRS : :class:`_tacoma.SIRS`
        An initialized SIRS object.
    is_static : bool, default : False
        The algorithm works a bit differently if it knows that the network is actually static.
        It works only with instances of :class:`_tacoma.edge_lists`.
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    None
        But the observables are saved in the :class:`_tacoma.SIRS` object.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        _tc.gillespie_SIRS_on_edge_changes(temporal_network, SIRS, verbose)
    elif type(temporal_network) == el:
        _tc.gillespie_SIRS_on_edge_lists(
            temporal_network, SIRS, is_static, verbose)
    else:
        raise ValueError('Unknown temporal network format: ' +
                         str(type(temporal_network)))


def convert_to_edge_trajectories(temporal_network, return_edge_similarities=False, verbose=False):
    """Converts a temporal network to an instance of :class:`_tacoma.edge_trajectories`.

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        An instance of a temporal network.
    return_edge_similarities : bool, optional
        If this is `True`, compute the similarity between edges, too. default : False
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    traj : :class:`_tacoma.edge_trajectories`
        The converted temporal network.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        if return_edge_similarities:
            raise ValueError(
                'Please convert to `edge_lists` first as the edge similarity algorithm is only implemented for `edge_lists`')

        result = _tc.edge_trajectories_from_edge_changes(
            temporal_network, return_edge_similarities, verbose)
    elif type(temporal_network) == el:
        result = _tc.edge_trajectories_from_edge_lists(
            temporal_network, return_edge_similarities)
    else:
        raise ValueError('Unknown temporal network format: ' +
                         str(type(temporal_network)))

    return result


def convert_edge_trajectories(traj):
    """
    Convert an instance of :class:`_tacoma.edge_trajectories`
    to an instance of :class:`_tacoma.edge_changes`. Note that there's
    no method to convert to :class:`_tacoma.edge_lists` directly.

    Parameters
    ----------
    traj : :class:`_tacoma.edge_trajectories`
        The temporal network which has to be converted.

    Returns
    -------
    temporal_network : :class:`_tacoma.edge_changes`
        The converted temporal network.
    """

    return _tc.convert_edge_trajectories_to_edge_changes(traj)


def get_edge_trajectories(temporal_network, return_edge_similarities=False, verbose=False):
    """
    Computes a list of time intervals in which each edge existed and optionally the
    similarity between edges.

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, :class:`_tacoma.edge_lists_with_histograms` or :class:`_tacoma.edge_trajectories`
        An instance of a temporal network.
    return_edge_similarities : bool, optional
        If this is `True`, return the similarity between edges. default : False
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    trajectories : :class:`_tacoma.edge_trajectories.trajectories`
        a list of `edge_trajectory_entry` objects, which contain
        the `.time_pairs` attributes, a list of time pairs (t0, t1)
        for time t0 <= t <= t1 in which the edge existed and `.edge`
        a pair of node integers.
    similarities : :class:`_tacoma.edge_trajectories.edge_similarities`
        a list of triples ( u, v, similarity ) where
        `u` and `v` refer to the edge indices in
        `edge_trajectories.trajectories` and similarity
        is the integrated time both edges were switched on
        while being connected to the same node.
    """

    if type(temporal_network) == edge_trajectories:
        result = temporal_network
    else:
        result = convert_to_edge_trajectories(
            temporal_network, return_edge_similarities, verbose)

    if return_edge_similarities:
        return result.trajectories, result.edge_similarities
    else:
        return result.trajectories

    return result


def verify(temporal_network, verbose=False):
    """Checks wether the temporal network is compliant with the demanded formats for analyses. Writes remarks on violations to stdout.

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        An instance of a temporal network.
    verbose: bool, optional
        If this is `True`, cout all errors that were found. default : False

    Returns
    -------
    int
        Number of found errors.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = _tc.verify_edge_changes(temporal_network, verbose)
    elif type(temporal_network) == el:
        result = _tc.verify_edge_lists(temporal_network, verbose)
    else:
        raise ValueError('Unknown temporal network format: ' +
                         str(type(temporal_network)))

    return result


def convert(temporal_network, verbose=False):
    """Converts either an instance of :class:`_tacoma.edge_changes` to an instance of :class:`_tacoma.edge_lists` or vice versa.

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        An instance of a temporal network.
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    An instance of the other format.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = _tc.convert_edge_changes(temporal_network, verbose=False)
    elif type(temporal_network) == el:
        result = _tc.convert_edge_lists(temporal_network, verbose=False)
    else:
        raise ValueError('Unknown temporal network format: ' +
                         str(type(temporal_network)))

    return result


def concatenate(list_of_temporal_networks, verbose = False):
    """Concatenates a list of either :class:`_tacoma.edge_changes` or :class:`_tacoma.edge_lists` to a single instance of :class:`_tacoma.edge_changes` or :class:`_tacoma.edge_lists`, respectively.

    Parameters
    ----------
    temporal_network : :obj:`list` of :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        A list of a temporal networks.
    verbose: bool, default : False
        Be chatty.

    Returns
    -------
    A single instance of the format provided in the lists. 
    """

    list_of_temporal_networks = [_get_raw_temporal_network(
        _t) for _t in list_of_temporal_networks]

    _t = list_of_temporal_networks[0]

    if type(_t) == ec:
        result = _tc.concatenate_edge_changes(
            list_of_temporal_networks, verbose)
    elif type(_t) == el:
        result = _tc.concatenate_edge_lists(list_of_temporal_networks, verbose)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(_t)))

    return result


def mean_degree(temporal_network):
    """Computes the time dependent mean degree of the temporal network.

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        A list of a temporal networks.

    Returns
    -------
    numpy array
        `t` containing the time points.
    numpy array
        `k` containing the mean degree at the corresponding times.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = _tc.mean_degree_from_edge_changes(temporal_network)
    elif type(temporal_network) == el:
        result = _tc.mean_degree_from_edge_lists(temporal_network)
    else:
        raise ValueError('Unknown temporal network format: ' +
                         str(type(temporal_network)))

    result = np.array(result)
    t, k = result[:, 0], result[:, 1]

    return t, k


def edge_counts(temporal_network):
    """Returns the number of edge events and total edge count (In this order: `m_in`, `m_out`, and `m`).

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        A list of a temporal networks.

    Returns
    -------
    :obj:`list` of `int`
        Number of edges coming into the network at the corresponding times.
    :obj:`list` of `int`
        Number of edges leaving the network at the corresponding times.
    :obj:`list` of `int`
        Number of edges existent in the network at the corresponding times.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == el:
        temporal_network = convert(temporal_network)

    if type(temporal_network) == ec:
        m_in, m_out, m = _tc.get_edge_counts(temporal_network)
    else:
        raise ValueError('Unknown temporal network format: ' +
                         str(type(temporal_network)))

    return m_in, m_out, m

