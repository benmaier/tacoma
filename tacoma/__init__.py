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

def flockwork_P_args(temporal_network,*args,**kwargs):
    """
        Bins an `edge_changes` instance for each `dt` (after each step, respectively,
        if `N_time_steps` was provided) and computes the rewiring rate and P 
        from the binned `edges_in` and `edges_out`. For DTU data use dt = 3600,
        for sociopatterns use dt = 600.

        py::arg("temporal_network"),             -- the temporal network
        py::arg("t_run_total") = None,           -- this is just plainly copied to the returned 
                                                    kwargs. If it is set to `None`, t_run_total
                                                    will be set to `temporal_network.tmax`
        py::arg("dt") = 0.0,                     -- dt of the time bin
        py::arg("N_time_steps") = 0,             -- number of time bins (use either this or dt).
        py::arg("aggregated_network") = {}       -- dict[edge] -> similarity, if this is given,
                                                    the kwargs are supposed to be for the function
                                                    `flockwork_P_varying_rates_neighbor_affinity`,
                                                    you can get this network from the `aggregated_network`
                                                    property from the results returned from
                                                    `measure_group_sizes_and_durations`
        py::arg("ensure_empty_network") = false, -- if this is True, bins where the original network
                                                    is empty (n_edges = 0) will be an artificially 
                                                    set high gamma with P = 0, such that nodes will
                                                    lose all edges.
        py::arg("use_preferential_node_selection") = false -- this is just plainly copied to 
                                                              the returned kwargs if `aggregated_network`
                                                              is not empty
        py::arg("verbose") = false

        Returns a dictionary of kwargs for the functions `flockwork_P_varying_rates` or 
        `flockwork_P_varying_rates_neighbor_affinity`.
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    new_kwargs = {}
    with_affinity = 'aggregated_network' in kwargs and len(kwargs['aggregated_network']) > 0
    if with_affinity and\
       'use_preferential_node_selection' in kwargs:
            new_kwargs['use_preferential_node_selection'] = kwargs.pop('use_preferential_node_selection')

    if 't_run_total' in kwargs and kwargs['t_run_total'] is not None:
        new_kwargs['t_run_total'] = kwargs.pop('t_run_total')
    elif 't_run_total' in kwargs and kwargs['t_run_total'] is None:
        new_kwargs['t_run_total'] = temporal_network.tmax
        kwargs.pop('t_run_total')
    else:
        new_kwargs['t_run_total'] = temporal_network.tmax

    if type(temporal_network) == el:
        temporal_network = convert(temporal_network)

    if type(temporal_network) == ec:
        kw = get_flockwork_P_args(temporal_network,*args,**kwargs)
        new_kwargs['E'] = kw.E
        new_kwargs['N'] = kw.N
        new_kwargs['P'] = kw.P
        new_kwargs['rewiring_rate'] = kw.rewiring_rate
        new_kwargs['tmax'] = kw.tmax
        if with_affinity:
            new_kwargs['neighbor_affinity'] = kw.neighbor_affinity
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(_t)))

    return new_kwargs


from .data_io import *
