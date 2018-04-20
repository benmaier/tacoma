# -*- coding: utf-8 -*-
"""
This module provides wrapper functions for the functions provided in the
C++-module `_tacoma`. The main purpose of the functions here is to check
the type of the provided temporal network and to call the appropriate
function from `_tacoma`.
"""

import numpy as np

from _tacoma import *

from _tacoma import edge_changes as ec
from _tacoma import edge_lists as el
from _tacoma import edge_lists_with_histograms as el_h
from _tacoma import edge_changes_with_histograms as ec_h


def _get_raw_temporal_network(temporal_network):

    if type(temporal_network) == ec_h:
        temporal_network = ec(temporal_network)
    elif type(temporal_network) == el_h:
        temporal_network = el(temporal_network)

    return temporal_network 


def measure_group_sizes_and_durations(temporal_network,*args,**kwargs):
    """Measures aggregated group size distribution, group size histogram for each time point, contact durations, inter-contact durations, durations of each group size, and the aggregated social network.

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        An instance of a temporal network.
    ignore_size_histogram : bool, optional
        Don't compute the single time point group size histograms (save time and memory). default : False
    verbose: bool, optional
        Be chatty.

    Returns 
    -------
    :mod:`group_sizes_and_durations`
        The result of the measurements.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        if "ignore_size_histograms" in kwargs:
            val = kwargs.pop("ignore_size_histograms")
            kwargs["ignore_size_histogram_differences"] = val
        result = measure_group_sizes_and_durations_for_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = measure_group_sizes_and_durations_for_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    return result
    
def degree_distribution(temporal_network,*args,**kwargs):
    """Measures the time-averaged degree distribution.

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        An instance of a temporal network.

    Returns 
    -------
    :obj:`list` of `float`
        A list of N entries where the k-th entry is the time-averaged probability of finding a node with degree k.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = degree_distribution_from_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = degree_distribution_from_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    return result

def bin(temporal_network,*args,**kwargs):
    """Bins a temporal network for each `dt` (after each step, respectively, if `N_time_steps` was provided).

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        An instance of a temporal network.
    dt : float
        The demanded bin size. default : 0.0
    N_time_steps : int
        Number of time bins (use either this or dt). default : 0
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    :mod:`edge_lists` 
        An edge_lists instance with one edge list describing the network as a list of all edges
        that were present in the last time bin.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = bin_from_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = bin_from_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    return result
    
def sample(temporal_network,*args,**kwargs):
    """Samples a temporal network after each `dt` (after each step, respectively, if `N_time_steps` was provided).

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        An instance of a temporal network.
    dt : float
        The demanded bin size. default : 0.0
    N_time_steps : int
        Number of time bins (use either this or dt). default : 0
    sample_aggregates : bool, optional
        If this is True, the following happens. If an edge is active during a time bin,
        it will appear in the final graph at the end of this time bin. It will then exist
        until the ende of the coming time bin. (This is different from the binning procedure). 
        default : False
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    :mod:`edge_lists` 
        An edge_lists instance with one edge list describing the network states after every time bin.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = sample_from_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = sample_from_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    return result
    
def binned_social_trajectory(temporal_network,*args,**kwargs):
    """Computes the binned social trajectory of node `node` in `temporal_network`, i.e. the groups it was part of and pairs of times (t0, t1) when it was part of a certain group. In this case, binning makes it easier to, e.g. see the days where a group was active.

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
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
        result = binned_social_trajectory_from_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = binned_social_trajectory_from_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    return result
    
def social_trajectory(temporal_network,*args,**kwargs):
    """Computes the social trajectory of node `node` in `temporal_network`, i.e. the groups it was part of and pairs of times (t0, t1) when it was part of a certain group.

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        An instance of a temporal network.
    node : int
        The node for which to compute the social trajectory
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    :obj:`list` of :mod:`social_trajectory_entry`
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = social_trajectory_from_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = social_trajectory_from_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    return result
    

def gillespie_SIS(temporal_network,*args,**kwargs):
    """Simulates an SIS process on the provided temporal network using the Gillespie stochastic simulation algorithm.``

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        An instance of a temporal network.
    SIS : :mod:`SIS`
        An initialized SIS object.
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    None
        But the observables are saved in the :mod:`SIS` object.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = gillespie_SIS_on_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = gillespie_SIS_on_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    return result
    
def gillespie_SI(temporal_network,*args,**kwargs):
    """Simulates an SI process on the provided temporal network using the Gillespie stochastic simulation algorithm.``

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        An instance of a temporal network.
    SI : :mod:`SI`
        An initialized SI object.
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    None
        But the observables are saved in the :mod:`SI` object.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = gillespie_SI_on_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = gillespie_SI_on_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    return result
    

def gillespie_SIR(temporal_network,*args,**kwargs):
    """Simulates an SIR process on the provided temporal network using the Gillespie stochastic simulation algorithm.``

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        An instance of a temporal network.
    SIR : :mod:`SIR`
        An initialized SIR object.
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    None
        But the observables are saved in the :mod:`SIR` object.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = gillespie_SIR_on_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = gillespie_SIR_on_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    return result
    
def gillespie_SIRS(temporal_network,*args,**kwargs):
    """Simulates an SIRS process on the provided temporal network using the Gillespie stochastic simulation algorithm.``

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        An instance of a temporal network.
    SIRS : :mod:`SIRS`
        An initialized SIRS object.
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    None
        But the observables are saved in the :mod:`SIRS` object.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = gillespie_SIRS_on_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = gillespie_SIRS_on_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))


def get_edge_trajectories(temporal_network,*args,**kwargs):
    """Computes the time intervals in which each edge existed.

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        An instance of a temporal network.
    return_edge_similarities : bool, optional
        If this is `True`, return the similarity between edges. default : False
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    :mod:`edge_trajectories`
        edge_trajectories.trajectories -- a list of `edge_trajectory_entry` objects, which contain
                                          the `.time_pairs` attributes, a list of time pairs (t0, t1)
                                          for time t0 <= t <= t1 in which the edge existed and `.edge`
                                          a pair of node integers.
        edge_trajectories.edge_similarities -- a list of triples ( u, v, similarity ) where
                                               `u` and `v` refer to the edge indices in
                                               `edge_trajectories.trajectories` and similarity
                                               is the integrated time both edges were switched on
                                               while being connected to the same node.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        if 'return_edge_similarities' in kwargs and\
           kwargs['return_edge_similarities']:
            raise ValueError('Please convert to `edge_lists` first as the edge similarity algorithm is only implemented for `edge_lists`')

        result = edge_trajectories_from_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = edge_trajectories_from_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    return result

def verify(temporal_network,*args,**kwargs):
    """Checks wether the temporal network is compliant with the demanded formats for analyses. Writes remarks on violations to stdout.

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        An instance of a temporal network.
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    int
        Number of found errors.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = verify_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = verify_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    return result

def convert(temporal_network,*args,**kwargs):
    """Converts either an instance of :mod:`edge_changes` to an instance of :mod:`edge_lists` or vice versa.

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        An instance of a temporal network.
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    An instance of the other format.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = convert_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = convert_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    return result

def concatenate(list_of_temporal_networks,*args,**kwargs):
    """Concatenates a list of either :mod:`edge_changes` or :mod:`edge_lists` to a single instance of :mod:`edge_changes` or :mod:`edge_lists`, respectively.

    Parameters
    ----------
    temporal_network : :obj:`list` of :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        A list of a temporal networks.
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    A single instance of the format provided in the lists. 
    """

    list_of_temporal_networks = [ _get_raw_temporal_network(_t) for _t in list_of_temporal_networks ]

    _t = list_of_temporal_networks[0]

    if type(_t) == ec:
        result = concatenate_edge_changes(list_of_temporal_networks,*args,**kwargs)
    elif type(_t) == el:
        result = concatenate_edge_lists(list_of_temporal_networks,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(_t)))

    return result

def mean_degree(temporal_network,*args,**kwargs):
    """Computes the time dependent mean degree of the temporal network.

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
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
        result = mean_degree_from_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = mean_degree_from_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    result = np.array(result)
    t, k = result[:,0], result[:,1]

    return t, k

def edge_counts(temporal_network):
    """Returns the number of edge events and total edge count (In this order: `m_in`, `m_out`, and `m`).

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
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
        m_in, m_out, m = get_edge_counts(temporal_network)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    return m_in, m_out, m

