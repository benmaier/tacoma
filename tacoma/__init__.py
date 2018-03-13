# -*- coding: utf-8 -*-
"""
TemporAl COntact Modeling and Analysis. Provides fast tools to analyze temporal contact networks, produce surrogate networks using qualitative models and simulate Gillespie processes on them.
"""

__version__ = "0.0.17"

import numpy as np

from _tacoma import *

from _tacoma import edge_changes as ec
from _tacoma import edge_lists as el
from _tacoma import edge_lists_with_histograms as el_h
from _tacoma import edge_changes_with_histograms as ec_h

color_sequence = [ u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf' ]
marker_sequence = ['s','d','o','X','v','<','^','.','>','h','p','P','*','8','H']

def complete_graph(N):
    """
        Get a single frame of a complete network.

        N -- number of nodes

        Returns an instance of `tacoma.edge_lists` with t = [0.0], tmax = 1.0
    """

    edge_list = []
    for i in range(N-1):
        for j in range(i+1,N):
            edge_list.append((i,j))

    this = el()
    this.t = [0.]
    this.tmax = 1.
    this.edges = [edge_list]
    this.N = N

    return this

def _get_raw_temporal_network(temporal_network):

    if type(temporal_network) == ec_h:
        temporal_network = ec(temporal_network)
    elif type(temporal_network) == el_h:
        temporal_network = el(temporal_network)

    return temporal_network 


def measure_group_sizes_and_durations(temporal_network,*args,**kwargs):
    """
        Measures aggregated group size distribution, group size histogram for each time point,
        contact durations, inter-contact durations, durations of each group size.

        py::arg("temporal_network"),
        py::arg("ignore_size_histograms") = false, -- don't compute the single time point group size histograms (save time and memory)
        py::arg("verbose") = false

        Returns an instance of `group_sizes_and_durations`.
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
    
def bin(temporal_network,*args,**kwargs):
    """
        Bins a temporal network for each `dt` (after each step, respectively,
        if `N_time_steps` was provided).

        py::arg("temporal_network"), -- the temporal network
        py::arg("dt") = 0.0,         -- dt of the time bin
        py::arg("N_time_steps") = 0, -- number of time bins (use either this or dt).
        py::arg("sample_aggregates") = false, -- if this is True, behavior is equal to `bin_temporal_network()`
        py::arg("verbose") = false

        Returns and instance of `edge_lists` with one edge list describing the network as a list of all edges
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
    """
        Samples a temporal network after each `dt` (after each step, respectively,
        if `N_time_steps` was provided).

        py::arg("temporal_network"), -- the temporal network
        py::arg("dt") = 0.0,         -- dt of the time bin
        py::arg("N_time_steps") = 0, -- number of time bins (use either this or dt).
        py::arg("sample_aggregates") = false, -- if this is True, behavior is equal to `bin_temporal_network()`
        py::arg("verbose") = false

        Returns and instance of `edge_lists` with one edge list describing the network states after every time bin.
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
    """
        Computes the binned social trajectory of node `node` in `temporal_network`, i.e.
        the groups it was part of and pairs of times (t0, t1) when it was part of
        a certain group.

        py::arg("temporal_network"), -- the temporal network
        py::arg("node"),             -- the node for which to compute the trajectory
        py::arg("dt") = 0.0,         -- dt of the time bin
        py::arg("N_time_steps") = 0, -- number of time bins (use either this or dt).
        py::arg("verbose") = false

        Returns a list, one entry for each dt. Each entry is a list of group ids.
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
    """
        Computes the social trajectory of node `node` in `temporal_network`, i.e.
        the groups it was part of and pairs of times (t0, t1) when it was part of
        a certain group.

        py::arg("temporal_network"), -- the temporal network
        py::arg("node"),             -- the node for which to compute the trajectory
        py::arg("verbose") = false

        Returns a list of `social_trajectory_entry`.
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
    """
        Simulates an SIS process on the provided temporal network

        py::arg("temporal_network"),
        py::arg("SIS"),
        py::arg("verbose") = false

        Returns `None`, but the observables are saved in the `SIS` object.
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
    """
        Simulates an SI process on the provided temporal network

        py::arg("temporal_network"),
        py::arg("SI"),
        py::arg("verbose") = false

        Returns `None`, but the observables are saved in the `SI` object.
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
    """
        Simulates an SIR process on the provided temporal network

        py::arg("temporal_network"),
        py::arg("SIR"),
        py::arg("verbose") = false

        Returns `None`, but the observables are saved in the `SIR` object.
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
    """
        Simulates an SIRS process on the provided temporal network

        py::arg("temporal_network"),
        py::arg("SIRS"),
        py::arg("verbose") = false

        Returns `None`, but the observables are saved in the `SIRS` object.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = gillespie_SIRS_on_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = gillespie_SIRS_on_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))


def get_edge_trajectories(temporal_network,*args,**kwargs):
    """
        Computes the times in which each edge existed.

        py::arg("temporal_network"),
        py::arg("return_edge_similarities"),
        py::arg("verbose") = false

        Returns an instance of `edge_trajectories`
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
        result = edge_trajectories_from_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = edge_trajectories_from_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    return result

def verify(temporal_network,*args,**kwargs):
    """
        Checks wether the temporal network is compliant with the demanded
        formats for analyses. Writes remarks on violations in the console.

        py::arg("temporal_network"),
        py::arg("verbose") = false
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
    """
        Converts either an instance of `edge_changes` to an instance of
        `edge_lists` or vice versa.

        py::arg("temporal_network"),
        py::arg("verbose") = false

        Returns an instance of the other format.
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
    """
        Concatenates a list of either `edge_changes` or `edge_lists`
        to a single instance of `edge_changes` or`edge_lists`, respectively.

        py::arg("list_of_temporal_networks"),
        py::arg("verbose") = false

        Returns a single instance of the format provided in the lists.
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
    """
        Computes the time dependent mean degree of the temporal network.

        py::arg("temporal_network"), -- the temporal network

        Returns, two numpy arrays, `t` containing the time points and `k` containing the mean degree.
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

from .data_io import *
from .tools import *
