# -*- coding: utf-8 -*-

from __future__ import print_function
import numpy as np
import numpy.polynomial.polynomial as poly
from numpy.polynomial import Polynomial

from scipy.optimize import fsolve
from scipy.optimize import minimize

from tacoma import edge_changes as ec
from tacoma import edge_lists as el
from tacoma import edge_lists_with_histograms as el_h
from tacoma import edge_changes_with_histograms as ec_h
from tacoma import get_logarithmic_histogram
from tacoma.power_law_fitting import fit_power_law_clauset
from tacoma import _get_raw_temporal_network
from tacoma import mean_coordination_number
from tacoma import convert
from tacoma import get_flockwork_P_args
from tacoma import get_flockwork_P_node_parameters_gamma_and_P
from tacoma import get_flockwork_P_node_parameters_alpha_and_beta_from_gamma_and_P
from tacoma import get_flockwork_alpha_beta_args
from tacoma import get_flockwork_node_parameters_alpha_and_beta
from tacoma import get_flockwork_node_rates_alpha_and_beta
import tacoma as tc

from scipy.stats import lognorm

# ========================================================= ZSBB ==================================================

def ZSBB_mean_coordination_number(b0,lam,N,b1):

    s = 0.0
    for n in range(1,N):
        s += (n+1.0) / ((n+1.0)*b1 - 1.0) * ((1.0-lam)/lam)**(n-1)

    pi_10 = ( 0.5 / ( b0 - (2*lam - 1.0)/(3*lam - 1.0)) + 0.5/lam * s )**(-1.0)

    return pi_10 / (2.0*lam) * s

def ZSBB_b0_func(b0,lam):
    return b0 * (2.0+(1.0-lam/(2*lam-1))) + 1

def estimate_ZSBB_args(temporal_network,
                       group_sizes_and_durations = None,
                       fit_discrete = False,
                       dt = None,
                       ):

    if fit_discrete and dt is None:
        raise ValueError('If the data is supposed to be treated as discrete, a value for `dt` must be provided in order to norm the group size durations.')
    elif not fit_discrete:
        dt = 1.

    result = group_sizes_and_durations
    if result is None:
        result = tc.measure_group_sizes_and_durations(temporal_network)

    m_in, m_out, m = tc.edge_counts(temporal_network)

    max_edge_events = sum(m_in) + sum(m_out)

    N = temporal_network.N

    # let's estimate b1 with the distribution of group_durations (pairs)
    #group_size = 2
    #b_ones = []
    #while group_size <= N: 
    #    values = result.group_durations[group_size] 
    #    alpha_1, err, xmin = fit_power_law_clauset(values)
    #    if (err/alpha_1)>0.05:
    #        break
    #    b_ones.append((alpha_1-1)/group_size)
    #    print b_ones
    #    group_size += 1
    #b1 = np.mean(b_ones)
    #print alpha_1

    values = result.group_durations[2] 
    alpha_1, err, xmin = fit_power_law_clauset(values)
    b1 = (alpha_1 - 1) / 2.0
    if b1 < 0.5:
        b1 = 0.51
    elif b1 > 1:
        b1 = 1.0

    # let's estimate b0 with the distribution of inter-contact durations

    if fit_discrete:
        values = np.array(result.group_durations[1]) / dt
        values = np.array(values,dtype=int)
        alpha_0, err, xmin = fit_power_law_clauset(values+1,discrete=True)
    else:
        alpha_0, err, xmin = fit_power_law_clauset(result.group_durations[1])

    mean_n = mean_coordination_number(result)

    def equations(p):
        b0, lam = p
        n = ZSBB_mean_coordination_number(b0,lam,N,b1)
        al0 = ZSBB_b0_func(b0,lam)
        return (n - mean_n, al0 - alpha_0)

    def cost(p):
        b0, lam = p
        n = ZSBB_mean_coordination_number(b0,lam,N,b1)
        al0 = ZSBB_b0_func(b0,lam)
        return np.abs((n - mean_n)/mean_n) + np.abs((al0 - alpha_0)/alpha_0)

    #b0, lam = fsolve(equations,(0.7,0.7))
    res = minimize(cost,(0.7,0.7))
    b0, lam = res.x

    if lam < 0.5:
        lam = 0.51
    elif lam > 1:
        lam = 1.0

    if b0 < 0.5:
        b0 = 0.51
    elif b0 > 1:
        b0 = 1.0

    if b0 <= (2*lam-1) / (3*lam-1.0):
        b0 = (2*lam-1) / (3*lam-1.0) + 0.01


    kwargs = {}
    kwargs['b0'] = b0
    kwargs['b1'] = b1
    kwargs['lambda'] = lam

    #temporal_network = _get_raw_temporal_network(temporal_network)

    kwargs['E'] = []
    kwargs['N'] = N
    kwargs['return_after_equilibration_only'] = True
    kwargs['t_equilibration'] = float(10000*N)
    kwargs['max_edge_events_to_end_simulation'] = max_edge_events

    return kwargs

# ==================================================== FLOCKWORK P-MODEL ===============================================

def estimate_flockwork_P_args_for_single_nodes(temporal_network,use_event_rate_method=False,*args,**kwargs):
    """Bins an `edge_changes` instance for each `dt` (after each step, respectively,
    if `N_time_steps` was provided) and computes the rewiring rate gamma and probability
    to stay alone P  from the binned `edges_in` and `edges_out`.
    Additionally

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        An instance of a temporal network.
    t_run_total : float
        this is just plainly copied to the returned kwargs. If it is set to `None`, t_run_total will be set to `temporal_network.tmax`
    dt : float
        The demanded bin size. default : 0.0
    N_time_steps : int
        Number of time bins (use either this or dt). default : 0
    aggregated_network : :obj:`dict` of :obj:`tuple` of int -> float, optional
        dict(edge -> similarity), if this is given,
        the kwargs are supposed to be for the function
        :mod:`flockwork_P_varying_rates_neighbor_affinity`,
        you can get this network from the `aggregated_network`
        property from the results returned by
        :mod:`measure_group_sizes_and_durations`. default : `{}`
    ensure_empty_network : bool, optional
        if this is True, bins where the original network
        is empty (n_edges = 0) will be an artificially 
        set high gamma with P = 0, such that nodes will
        lose all edges. default : False
    use_preferential_node_selection : bool, optional
        this is just plainly copied to 
        the returned kwargs if `aggregated_network`
        is not empty. default : False
    verbose: bool, optional
        Be chatty.


    Returns
    -------
    :obj:`dict`
        kwargs for the functions :mod:`flockwork_P_varying_rates` or 
        :mod:`flockwork_P_varying_rates_neighbor_affinity`, if `aggregated_network` was provided
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == el:
        temporal_network = convert(temporal_network)

    if type(temporal_network) == ec:
        # kwargs.pop('use_event_rate_method')
        kw = estimate_flockwork_P_args(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(_t)))

    if not use_event_rate_method:

        # get rewiring rate and P
        gamma = np.array(kw['rewiring_rate'])
        t, gamma_t = gamma[:,0], gamma[:,1]
        P_t = np.array(kw['P'])

        # convert to active reconnection event rate and active disconnection event rate
        alpha = gamma_t * P_t
        beta = gamma_t * (1 - P_t)

        # compute single node rewiring rate and P factors
        g_node, P_node = get_flockwork_P_node_parameters_gamma_and_P(temporal_network, kw['rewiring_rate'], kw['P'])
        g_node = np.array(g_node)
        P_node = np.array(P_node)

        new_g_node = []
        new_P_node = []

        # for each time bin
        for t, g, P, a, b in zip(t, gamma_t, P_t, alpha, beta):

            # convert to active reconnection event rates and active disconnection event rates
            a_node = np.array(g_node*g) * np.array(P_node*P)
            b_node = np.array(g_node*g) * (1-np.array(P_node*P))

            # set negative rate factors to zero and norm factors
            b_node[b_node<0] = 0.0
            a_node[a_node<0] = 0.0
            if not np.all(a_node == 0.0):
                a_node /= a_node.mean()
            if not np.all(b_node == 0.0):
                b_node /= b_node.mean()

            # multiply activity rates with normed factors s.t. the mean of a and b is conserved
            a_node *= a
            b_node *= b

            # convert back to gamma and P
            new_g = a_node + b_node
            _temp_new_g = new_g.copy()
            _temp_new_g[new_g==0.0] = 1.0
            new_P = a_node / _temp_new_g

            new_g_node.append((t, new_g.tolist()))
            new_P_node.append(new_P.tolist())
    else:
        # get rewiring rate and P
        gamma = np.array(kw['rewiring_rate'])
        t, gamma_t = gamma[:,0], gamma[:,1]
        P_t = np.array(kw['P'])

        # convert to active reconnection event rate and active disconnection event rate
        alpha = gamma_t * P_t
        beta = gamma_t * (1 - P_t)

        # compute single node rewiring rate and P factors
        a_node, b_node = get_flockwork_P_node_parameters_alpha_and_beta_from_gamma_and_P(temporal_network, kw['rewiring_rate'], kw['P'])
        a_node = np.array(a_node)
        b_node = np.array(b_node)

        print(gamma_t)

        new_g_node = []
        new_P_node = []

        # for each time bin
        for t, g, P, a, b in zip(t, gamma_t, P_t, alpha, beta):

            this_a_node = a_node.copy()
            this_b_node = b_node.copy()

            if not np.all(this_a_node == 0.0):
                this_a_node /= this_a_node.mean()
            if not np.all(this_b_node == 0.0):
                this_b_node /= this_b_node.mean()

            # multiply activity rates with normed factors s.t. the mean of a and b is conserved
            this_a_node *= a
            this_b_node *= b

            # convert back to gamma and P
            new_g = this_a_node + this_b_node
            _temp_new_g = new_g.copy()
            _temp_new_g[new_g==0.0] = 1.0
            new_P = this_a_node / _temp_new_g

            new_g_node.append((t, new_g.tolist()))
            new_P_node.append(new_P.tolist())
        

    kw['P'] = new_P_node
    kw['rewiring_rates'] = new_g_node

    kw.pop('rewiring_rate')

    return kw

    
def estimate_flockwork_P_args(temporal_network,*args,**kwargs):
    """Bins an `edge_changes` instance for each `dt` (after each step, respectively,
    if `N_time_steps` was provided) and computes the rewiring rate gamma and probability
    to stay alone P  from the binned `edges_in` and `edges_out`. For DTU data use 
    dt = 3600, for sociopatterns use dt = 600.

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        An instance of a temporal network.
    t_run_total : float
        this is just plainly copied to the returned kwargs. If it is set to `None`, t_run_total will be set to `temporal_network.tmax`
    dt : float
        The demanded bin size. default : 0.0
    N_time_steps : int
        Number of time bins (use either this or dt). default : 0
    aggregated_network : :obj:`dict` of :obj:`tuple` of int -> float, optional
        dict(edge -> similarity), if this is given,
        the kwargs are supposed to be for the function
        :mod:`flockwork_P_varying_rates_neighbor_affinity`,
        you can get this network from the `aggregated_network`
        property from the results returned by
        :mod:`measure_group_sizes_and_durations`. default : `{}`
    ensure_empty_network : bool, optional
        if this is True, bins where the original network
        is empty (n_edges = 0) will be an artificially 
        set high gamma with P = 0, such that nodes will
        lose all edges. default : False
    use_preferential_node_selection : bool, optional
        this is just plainly copied to 
        the returned kwargs if `aggregated_network`
        is not empty. default : False
    verbose: bool, optional
        Be chatty.


    Returns
    -------
    :obj:`dict`
        kwargs for the functions :mod:`flockwork_P_varying_rates` or 
        :mod:`flockwork_P_varying_rates_neighbor_affinity`, if `aggregated_network` was provided
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

def estimate_flockwork_alpha_beta_args_for_single_nodes(temporal_network,
                                                        estimate_for_each_time_bin=False,
                                                        apply_mean_correction=False,*args,**kwargs):
    """Bins an `edge_changes` instance for each `dt` (after each step, respectively,
    if `N_time_steps` was provided) and computes the reconnection rate alpha and disconnection
    rate beta from the binned `edges_in` and `edges_out`.
    Additionally

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        An instance of a temporal network.
    t_run_total : float
        this is just plainly copied to the returned kwargs. If it is set to `None`, t_run_total will be set to `temporal_network.tmax`
    dt : float
        The demanded bin size. default : 0.0
    N_time_steps : int
        Number of time bins (use either this or dt). default : 0
    aggregated_network : :obj:`dict` of :obj:`tuple` of int -> float, optional
        dict(edge -> similarity), if this is given,
        the kwargs are supposed to be for the function
        :mod:`flockwork_P_varying_rates_neighbor_affinity`,
        you can get this network from the `aggregated_network`
        property from the results returned by
        :mod:`measure_group_sizes_and_durations`. default : `{}`
    ensure_empty_network : bool, optional
        if this is True, bins where the original network
        is empty (n_edges = 0) will be an artificially 
        set high gamma with P = 0, such that nodes will
        lose all edges. default : False
    verbose: bool, optional
        Be chatty.


    Returns
    -------
    :obj:`dict`
        kwargs for the functions :mod:`flockwork_P_varying_rates` or 
        :mod:`flockwork_P_varying_rates_neighbor_affinity`, if `aggregated_network` was provided
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    verbose = 'verbose' in kwargs and kwargs['verbose']

    ensure_empty_network = 'ensure_empty_network' in kwargs and kwargs['ensure_empty_network']
    k_over_k_real_scaling = 1.0
    if 'k_over_k_real_scaling' in kwargs:
        k_over_k_real_scaling = kwargs['k_over_k_real_scaling']

    if type(temporal_network) == el:
        temporal_network = convert(temporal_network)

    if type(temporal_network) == ec:
        # kwargs.pop('use_event_rate_method')
        kw = estimate_flockwork_alpha_beta_args(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(_t)))

    if not estimate_for_each_time_bin:

        # get rewiring rate and P
        alpha = np.array(kw['reconnection_rate'])
        t, alpha = alpha[:,0], alpha[:,1]
        beta = np.array(kw['disconnection_rate'])

        # compute single node rewiring rate and P factors
        #print("estimating node specific parameters")
        a_node, b_node = get_flockwork_node_parameters_alpha_and_beta(temporal_network, 
                                                                      kw['reconnection_rate'], 
                                                                      kw['disconnection_rate'],
                                                                      k_over_k_real_scaling,
                                                                      apply_mean_correction,
                                                                      verbose = verbose,
                                                                     )
        #print("done")
        a_node = np.array(a_node)
        b_node = np.array(b_node)

        new_a_node = []
        new_b_node = []

        if ensure_empty_network:
            m_in, m_out, m = tc.edge_counts(temporal_network)
            t_m = np.array([temporal_network.t0] + temporal_network.t + [temporal_network.tmax])
            m = np.array(m)

            last_bin_was_for_emptying = False

        # for each time bin
        for it, (t_, a, b) in enumerate(zip(t, alpha, beta)):

            # multiply activity rates with normed factors s.t. the mean of a and b is conserved
            this_a_node = a_node * a
            this_b_node = b_node * b

            if ensure_empty_network:
                
                if it+2 < len(t):

                    # compute the number of total edges in this time bin and the next one
                    ndx1 = np.where(np.logical_and(t_m >= t[it], t_m < t[it+1]))
                    ndx2 = np.where(np.logical_and(t_m >= t[it+1], t_m < t[it+2]))

                    m1 = np.sum(m[ndx1])
                    m2 = np.sum(m[ndx2])

                    # if there's no edges in this and the next time bin, ensure that all edges are removed in here
                    if m1 == 0 and m2 == 0:
                        if not last_bin_was_for_emptying:
                            this_a_node = np.zeros_like(this_a_node)
                            this_b_node = np.ones_like(this_b_node) * np.log(temporal_network.N) / (t[it+1] - t[it])
                        last_bin_was_for_emptying = True
                    else:
                        last_bin_was_for_emptying = False


            new_a_node.append((t_, this_a_node.tolist()))
            new_b_node.append(this_b_node.tolist())

        kw['reconnection_rates'] = new_a_node
        kw['disconnection_rates'] = new_b_node

        kw.pop('disconnection_rate')
        kw.pop('reconnection_rate')
    else:
        alpha_t, beta_t = get_flockwork_node_rates_alpha_and_beta(temporal_network, kw['reconnection_rate'], kw['disconnection_rate'],apply_mean_correction)

        mean_alpha = np.array(kw['reconnection_rate'])
        t, mean_alpha = mean_alpha[:,0], mean_alpha[:,1]

        new_alpha = list(zip(t, alpha_t))

        kw['reconnection_rates'] = new_alpha
        kw['disconnection_rates'] = beta_t

        kw.pop('disconnection_rate')
        kw.pop('reconnection_rate')

    return kw

def estimate_flockwork_alpha_beta_args_for_single_nodes_randomly(temporal_network,sigma,
                                                                 *args,**kwargs):
    """Bins an `edge_changes` instance for each `dt` (after each step, respectively,
    if `N_time_steps` was provided) and computes the reconnection rate alpha and disconnection
    rate beta from the binned `edges_in` and `edges_out`.
    Additionally

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        An instance of a temporal network.
    t_run_total : float
        this is just plainly copied to the returned kwargs. If it is set to `None`, t_run_total will be set to `temporal_network.tmax`
    dt : float
        The demanded bin size. default : 0.0
    N_time_steps : int
        Number of time bins (use either this or dt). default : 0
    aggregated_network : :obj:`dict` of :obj:`tuple` of int -> float, optional
        dict(edge -> similarity), if this is given,
        the kwargs are supposed to be for the function
        :mod:`flockwork_P_varying_rates_neighbor_affinity`,
        you can get this network from the `aggregated_network`
        property from the results returned by
        :mod:`measure_group_sizes_and_durations`. default : `{}`
    ensure_empty_network : bool, optional
        if this is True, bins where the original network
        is empty (n_edges = 0) will be an artificially 
        set high gamma with P = 0, such that nodes will
        lose all edges. default : False
    verbose: bool, optional
        Be chatty.

    Returns
    -------
    :obj:`dict`
        kwargs for the functions :mod:`flockwork_P_varying_rates` or 
        :mod:`flockwork_P_varying_rates_neighbor_affinity`, if `aggregated_network` was provided
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if type(temporal_network) == el:
        temporal_network = convert(temporal_network)

    if type(temporal_network) == ec:
        # kwargs.pop('use_event_rate_method')
        kw = estimate_flockwork_alpha_beta_args(temporal_network,*args,**kwargs)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(_t)))

    s = sigma
    scale = 1.0
    # get rewiring rate and P
    alpha = np.array(kw['reconnection_rate'])
    t, alpha = alpha[:,0], alpha[:,1]
    beta = np.array(kw['disconnection_rate'])

    # compute single node rewiring rate and P factors
    a_node = lognorm.rvs(s, scale=scale, size=temporal_network.N)
    b_node = lognorm.rvs(s, scale=scale, size=temporal_network.N)

    a_node /= a_node.mean()
    b_node /= b_node.mean()

    new_a_node = []
    new_b_node = []

    # for each time bin
    for t, a, b in zip(t, alpha, beta):

        # multiply activity rates with normed factors s.t. the mean of a and b is conserved
        this_a_node = a_node * a
        this_b_node = b_node * b

        new_a_node.append((t, this_a_node.tolist()))
        new_b_node.append(this_b_node.tolist())

    kw['reconnection_rates'] = new_a_node
    kw['disconnection_rates'] = new_b_node

    kw.pop('disconnection_rate')
    kw.pop('reconnection_rate')

    return kw

def estimate_flockwork_alpha_beta_args(temporal_network,*args,**kwargs):
    """Bins an `edge_changes` instance for each `dt` (after each step, respectively,
    if `N_time_steps` was provided) and computes the reconnection rate alpha and 
    disconnection rate beta from the binned `edges_in` and `edges_out`. For DTU data use 
    dt = 3600, for sociopatterns use dt = 600.

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        An instance of a temporal network.
    t_run_total : float
        this is just plainly copied to the returned kwargs. If it is set to `None`, t_run_total will be set to `temporal_network.tmax`
    dt : float
        The demanded bin size. default : 0.0
    N_time_steps : int
        Number of time bins (use either this or dt). default : 0
    aggregated_network : :obj:`dict` of :obj:`tuple` of int -> float, optional
        dict(edge -> similarity), if this is given,
        the kwargs are supposed to be for the function
        :mod:`flockwork_P_varying_rates_neighbor_affinity`,
        you can get this network from the `aggregated_network`
        property from the results returned by
        :mod:`measure_group_sizes_and_durations`. default : `{}`
    ensure_empty_network : bool, optional
        if this is True, bins where the original network
        is empty (n_edges = 0) will be an artificially 
        set high gamma with P = 0, such that nodes will
        lose all edges. default : False
    use_preferential_node_selection : bool, optional
        this is just plainly copied to 
        the returned kwargs if `aggregated_network`
        is not empty. default : False
    verbose: bool, optional
        Be chatty.


    Returns
    -------
    :obj:`dict`
        kwargs for the functions :mod:`flockwork_P_varying_rates` or 
        :mod:`flockwork_P_varying_rates_neighbor_affinity`, if `aggregated_network` was provided
    """

    temporal_network = _get_raw_temporal_network(temporal_network)

    if 'aggregated_network' in kwargs and (kwargs['aggregated_network'] is None):
        kwargs['aggregated_network'] = {}

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
        kw = get_flockwork_alpha_beta_args(temporal_network,*args,**kwargs)
        new_kwargs['E'] = kw.E
        new_kwargs['N'] = kw.N
        new_kwargs['disconnection_rate'] = kw.disconnection_rate
        new_kwargs['reconnection_rate'] = kw.reconnection_rate
        new_kwargs['tmax'] = kw.tmax
        if with_affinity:
            new_kwargs['neighbor_affinity'] = kw.neighbor_affinity
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(_t)))

    return new_kwargs 

# =============================================== DYNAMIC RGG ===============================================

def estimate_dynamic_RGG_args(sampled_or_binned_temporal_network,
                              dt = None,
                              periodic_boundary_conditions_for_link_building = True,
                              group_sizes_and_durations = None,
                              ):

    if not periodic_boundary_conditions_for_link_building:
        raise ValueError('So far, only parameter estimation for periodic_boundary_conditions_for_link_building = True has been implemented')

    if group_sizes_and_durations is not None:
        result = group_sizes_and_durations
    else:
        result = tc.measure_group_sizes_and_durations(sampled_or_binned_temporal_network)

    mean_m = tc.mean_group_size(result)
    log_y = np.log(mean_m)

    params = [0.7895411,  1.28318048, 0.0]
    params[-1] -= log_y

    def get_root(p):
        a, b, c = params
        r = b**2 - 4*a*c

        if r >= 0:
            r1 = (-b + np.sqrt(r)) / (2.0*a)     
            r2 = (-b - np.sqrt(r)) / (2.0*a)
            return max(r1,r2)
        else:
            return None
        """
        elif r == 0:
            r = -b/2.0*a
            return r
        """

    density = get_root(params)

    """
    print params
    this_poly = Polynomial(params)
    print this_poly
    roots = poly.polyroots(this_poly)
    print roots
    density = roots.max()
    """

    if dt is None:
        dt = sampled_or_binned_temporal_network.t[1] - sampled_or_binned_temporal_network.t[0]

    mean_link_duration = np.array(result.contact_durations).mean() / dt

    kwargs = {}
    kwargs['N'] = sampled_or_binned_temporal_network.N
    kwargs['t_run_total'] = len(sampled_or_binned_temporal_network.t)
    kwargs['mean_link_duration'] = mean_link_duration
    kwargs['critical_density'] = density

    return kwargs


if __name__ == "__main__":
    import tacoma as tc
    import matplotlib.pyplot as pl
    from tacoma.analysis import temporal_network_group_analysis
    

    test = tc.dtu_week()
    rewiring_rate = test.gamma
    P = test.P

    fw = tc.flockwork_P_varying_rates([],100,P,24*3600,rewiring_rate,tmax=24*3600*7)
    fw_binned = tc.sample(fw,dt=300)
    fw_binned_result = tc.measure_group_sizes_and_durations(fw_binned)

    kwargs = get_ZSBB_parameters(fw_binned,fw_binned_result,fit_discrete=True,dt=300.)
    print("lambda =", kwargs['lambda'])
    print("b0 =", kwargs['b0'])
    print("b1 =", kwargs['b1'])
    kwargs['t_run_total'] = (len(fw_binned.t) + 1)*kwargs['N']
    zsbb = tc.ZSBB_model(**kwargs)
    zsbb_binned = tc.sample(zsbb,dt=kwargs['N'])
    zsbb_binned_result = tc.measure_group_sizes_and_durations(zsbb_binned)


    fig, ax, data = temporal_network_group_analysis(fw_binned_result)
    temporal_network_group_analysis(zsbb_binned_result,
                                    time_normalization_factor = 300./kwargs['N'],
                                    ax=ax)

    pl.show()
