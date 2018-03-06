# -*- coding: utf-8 -*-
"""
TemporAl COntact Modeling and Analysis. Provides fast tools to analyze temporal contact networks, produce surrogate networks using qualitative models and simulate Gillespie processes on them.
"""

__version__ = "0.0.17"

from _tacoma import *

from _tacoma import edge_changes as ec
from _tacoma import edge_lists as el
from _tacoma import edge_lists_with_histograms as el_h
from _tacoma import edge_changes_with_histograms as ec_h

color_sequence = [ u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf' ]
marker_sequence = ['s','d','o','X','v','<','^','.','>','h','p','P','*','8','H']

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
        result = tc.measure_group_sizes_and_durations_for_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = tc.measure_group_sizes_and_durations_for_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    return result
    
def bin_temporal_network(temporal_network,*args,**kwargs):
    """
        Bins a temporal network for each `dt` (after each step, respectively,
        if `N_time_steps` was provided).

        py::arg("temporal_network"), -- the temporal network
        py::arg("node"),             -- the node for which to compute the trajectory
        py::arg("dt") = 0.0,         -- dt of the time bin
        py::arg("N_time_steps") = 0, -- number of time bins (use either this or dt).
        py::arg("sample_aggregates") = false, -- if this is True, behavior is equal to `bin_temporal_network()`
        py::arg("verbose") = false

        Returns and instance of `edge_lists` with one edge list describing the network as a list of all edges
        that were present in the last time bin.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = tc.bin_from_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = tc.bin_from_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    return result
    
def sample_temporal_network(temporal_network,*args,**kwargs):
    """
        Samples a temporal network after each `dt` (after each step, respectively,
        if `N_time_steps` was provided).

        py::arg("temporal_network"), -- the temporal network
        py::arg("node"),             -- the node for which to compute the trajectory
        py::arg("dt") = 0.0,         -- dt of the time bin
        py::arg("N_time_steps") = 0, -- number of time bins (use either this or dt).
        py::arg("sample_aggregates") = false, -- if this is True, behavior is equal to `bin_temporal_network()`
        py::arg("verbose") = false

        Returns and instance of `edge_lists` with one edge list describing the network states after every time bin.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = tc.sample_from_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = tc.sample_from_edge_lists(temporal_network,*args,**kwargs)
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
        result = tc.binned_social_trajectory_from_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = tc.binned_social_trajectory_from_edge_lists(temporal_network,*args,**kwargs)
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
        result = tc.social_trajectory_from_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = tc.social_trajectory_from_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    return result
    

def gillespie_SIS(temporal_network,*args,**kwargs):
    """
        Simulates an SIS process on the provided temporal network

        py::arg("temporal_network"),
        py::arg("Dyn_SIS"),
        py::arg("verbose") = false

        Returns `None`, but the observables are saved in the `Dyn_SIS` object.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == ec:
        result = tc.gillespie_SIS_on_edge_changes(temporal_network,*args,**kwargs)
    elif type(temporal_network) == el:
        result = tc.gillespie_SIS_on_edge_lists(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    return result
    



