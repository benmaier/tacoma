# -*- coding: utf-8 -*-
"""
This module provides functions related to the flockwork
temporal network model.
"""
from copy import deepcopy

import numpy as np
from numpy import random

from scipy.integrate import ode
from scipy.integrate import simps
from scipy.special import gamma as Gamma

from tacoma.power_law_fitting import fit_power_law_clauset

from _tacoma import flockwork_P_varying_rates
from tacoma import _get_raw_temporal_network

import tacoma as tc



def flockwork_P_equilibrium_group_size_distribution(N, P):
    """Get the equilibrium group size distribution of a Flockwork-P model
    given node number N and probability to reconnect P.

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

    N = int(N)
    P = float(P)

    assert N > 2
    assert P >= 0.
    assert P <= 1.

    if P == 0.:
        dist = np.zeros((N+1,))
        dist[1] = N
    elif P == 1.:
        dist = np.zeros((N+1,))
        dist[-1] = 1.
    else:
        dist = [N*(1.-P)]

        N_fak = N - np.arange(1, N-1)
        j_fak = ((N-1)-P * np.arange(1, N-1))
        div = N_fak / j_fak
        cum_product_div = np.cumprod(div)
        for m in range(2, N):
            #dist.append( (-1)**(m%2) * float(N)/float(m) * (P-1.) * np.prod(N_fak[1:m]/j_fak[1:m]) * P**(m-1) )
            dist.append(float(N)/float(m) *
                        (1-P) * cum_product_div[m-2] * P**(m-1))

        value = P
        for j in range(1, N-1):
            #value *= float(N-j-1) / ((P-N+1.) / P + (j-1))
            value *= float(N-j-1) * P / ((N-1)-P*j)
        #value *= P**(N-1)
        dist.append(value)

        dist = [0.] + dist

    return np.array(dist)


def group_size_distribution_asymptotics(N, P, mmax=None, simple_pochhammer_approximation=True):
    """Get the asymptotic equilibrium group size distribution of a Flockwork-P model
    given node number N and probability to reconnect P.

    Parameters
    ----------
    N : int
        Number of nodes
    P : float
        Probability to reconnect
    mmax : int, default = None
        The maximum group size for which to calculate the asymptotics
    simple_pochhammer_approximation : bool, default = True
        Whether to use a simple first order Pochhammer approximation
        or the Gamma-function approximation.

    Returns
    -------
    ms : numpy.ndarray
        group size vector
    dist : numpy.ndarray
        Asymptotic group size distribution
    """

    N = int(N)
    P = float(P)

    assert N > 2
    assert P >= 0.
    assert P < 1.

    if mmax is None:
        mmax = int(N) // 2

    ms = np.arange(1,mmax+1)

    if P == 0.:
        ms = np.arange(N+1)
        dist = np.zeros((N+1,))
        dist[1] = N
    else:
        dist = [(1.-P)]

        for m in ms[1:]:
            if simple_pochhammer_approximation:
                factor = ( 1 - P/(N-1)*(m-1)*(m-2)/2 )**(-1)
            else:
                factor = ((-1)**(m % 2) * (P*m/(N-1)/np.exp(1))**m * 2*np.pi / Gamma(-(N-1)/P) * m**(-(N-1)/P - 0.5))**(-1)
            this = 1/m * (P/np.exp(1))**(m-1) * ((N-1)/(N-m))**(N-m+0.5) * factor

            dist.append(this)

    return ms, N*np.array(dist)

def flockwork_P_equilibrium_configuration(N, 
                                          P, 
                                          shuffle_nodes=True, 
                                          return_histogram=False, 
                                          seed=0, 
                                          shuffle_group_sizes=False,
                                          dist=None,
                                          use_binomial=True,
                                          construct_edges=True,
                                          ):
    """Get an equilibrium configuration of a Flockwork-P model
    given node number N and probability to reconnect P.

    Parameters
    ----------
    N : int
        Number of nodes
    P : float
        Probability to reconnect
    shuffle_nodes : bool, default : True
        Shuffle the node order in which nodes are distributed
        to groups. 'True' is recommended.
    return_group_size_histogram : bool, default : False
        Return a numpy array containing the counts of groups
        of size :math:`g`.
    seed : int, default : 0
        The random seed. RNG is initialized randomly if ``seed = 0``.
    shuffle_group_sizes : bool, default : False
        Shuffle the order of group sizes in which nodes are distributed
        to groups. 'True' is recommended.
    dist : np.array, default : None
        Is computed from N and P if `None` is provided
    use_binomial : bool, default : True
        whether to draw from a Binomial instead of a Poisson
        (conceptually more sound)
    construct_edgelist : bool, default : True
        if to actually construct the edgelist.
        
    
    Returns
    -------
    :obj:`list` of :obj:`tuple` of int
        edge list of equilibrium configuration
    numpy.ndarray
        group size counter of this configuration (only
        if ``return_group_size_histogram`` is `True`)
    """

    if seed > 0:
        random.seed(seed)

    if dist is None:
        dist = flockwork_P_equilibrium_group_size_distribution(N, P)

        dist = dist[1:]

    total_group_number = sum(dist)
    size_dist = dist / total_group_number
    num_components = []
    edges = []

    if shuffle_nodes:
        node_ints = random.permutation(N)
    else:
        node_ints = np.arange(N)

    # in the beginning, theres N nodes left
    nodes_left = N

    # ... and no group counter contains any groups yet
    C_m = np.zeros(N)

    # while there are still nodes left to distribute
    while nodes_left > 0:

        # loop through group sizes in descending order.
        # start with the smallest group size that
        # may contain all of the nodes left
        if shuffle_group_sizes:
            group_sizes = 1 + random.permutation(nodes_left)
        else:
            group_sizes = np.arange(1,nodes_left+1)[::-1]
        for m in group_sizes:

            if (nodes_left-m) < 0:
                continue

            # if the expected number of groups of this size is not zero
            if dist[m-1] > 0. and nodes_left >= m:

                #new_N_m = random.binomial(total_group_number,size_dist[m-1])
                # dist carries the expected number of m-groups in equilibrium,
                # so we draw a number of m groups from a Poisson distribution
                # with mean dist[m-1]
                if use_binomial:
                    max_draws = N // m
                    p = dist[m-1] / max_draws
                    new_C_m = random.binomial(max_draws,p)
                else:
                    new_C_m = random.poisson(dist[m-1])

                # if the new number of groups of size m is larger than the previously drawn number
                if new_C_m > C_m[m-1]:
                    # accept the new number and add the difference in group count of this group size
                    delta_C_m = int(new_C_m - C_m[m-1])
                elif nodes_left == 1 and m == 1:
                    # if there's only one node left, just add it to the loners
                    delta_C_m = 1
                else:
                    delta_C_m = 0

                # add the additional number of groups of this size
                for group_instance in range(delta_C_m):
                    # for u in xrange(nodes_left-1,nodes_left-m,-1):
                    #    for v in xrange(u-1,nodes_left-m-1,-1):
                    #        edges.append((node_ints[u],node_ints[v]))
                    if nodes_left - m < 0:
                        break

                    # add fully connected clusters to the edge set
                    if construct_edges:
                        if m > 1 and nodes_left-m >= 0:
                            edges.extend([tuple(sorted((int(node_ints[u]), int(node_ints[v]))))
                                          for u in range(nodes_left-m, nodes_left-1)
                                          for v in range(u+1, nodes_left)])
                        elif nodes_left - m < 0:
                            break

                    # remove the grouped nodes from the pool of remaining nodes
                    nodes_left -= m

                    # save this group as an accepted group
                    # and increment the group count
                    C_m[m-1] += 1

                    """
                    # if there's no nodes left to distribute,
                    # N_m was chosen too large and we abandon
                    if nodes_left <= 0.:
                        break
                    """

            # go to next smaller group size

    if return_histogram:
        return edges, np.append(np.array([0.]), C_m)
    else:
        return edges


def flockwork_P_mean_degree_for_varying_rates(flockwork_P_params, N=None):
    """Compute the mean group size distribution for a Flockwork-P system with varying rates.

    Parameters
    ----------
    flockwork_P_params : :obj:`dict`
        Contains all parameters necessary for a Flockwork-P simulation, especially the 
        time-dependent rewiring rate and time-dependent reconnection probability
    N : int, default : None
        If given, compute everything for `N` nodes, where `N` is different from 
        `N` in `flockwork_P_params`.

    Returns
    -------
    t : numpy.array
        An array of times at which the mean degree was evaluated
    k : numpy.array
        An array of mean degree values corresponding to the times in t.
    """

    data = flockwork_P_params

    if N is None:
        N = data['N']

    # compute the initial mean degree
    k0 = 2*len(data['E']) / float(N)

    # get arrays for rewiring
    gamma = np.array(data['rewiring_rate'])
    T, gamma = gamma[:, 0], gamma[:, 1]
    T = np.append(T, data['tmax'])
    gamma = np.append(gamma, gamma[-1])

    P = np.array(data['P'])
    P = np.append(P, P[-1])

    # define the linear ODE describing the mean degree evolution
    def dkdt(t, k, g_, P_):
        return 2 * g_ * P_ - 2*k * (g_*(1-P_))

    # intialize the integrator
    r = ode(dkdt)
    r.set_integrator('dopri5')

    k = [k0]
    new_t = [T[0]]

    i = 0
    for t_, g_, P_ in zip(T[:-1], gamma[:-1], P[:-1]):

        # for every interval set the initial mean degree
        # and new parameters
        r.set_initial_value(k[-1], t_)
        r.set_f_params(g_, P_)

        # for 10 time points within this interval,
        # compute the expected mean degree
        this_t = np.linspace(t_, T[i+1], 10)
        for t__ in this_t[1:]:
            new_t.append(t__)
            this_k = r.integrate(t__)
            k.append(this_k)

        # increase interval
        i += 1

    return np.array(new_t), np.array(k)

def flockwork_P_mean_degree_for_varying_alpha_beta(N,k_initial,t,alpha,beta,tmax):
    r"""Compute the mean group size distribution for a Flockwork-P system with varying rates.

    Parameters
    ----------
    N : int
        Number of nodes
    k_initial : float
        initial mean degree
    t : numpy.ndarray of float
        time points at which :math:`\alpha(t)` and 
        :math:`\beta(t)` change
    alpha : numpy.ndarray of float
        active reconnection rate associated with
        the time points in ``t``
    beta : numpy.ndarray of float
        active disconnection rate associated with
        the time points in ``t``
    tmax : float
        final time

    Returns
    -------
    t : numpy.array
        An array of times at which the mean degree was evaluated
    k : numpy.array
        An array of mean degree values corresponding to the times in t.
    """

    # compute the initial mean degree
    k0 = k_initial

    # get arrays for rewiring
    t = np.append(t, tmax)
    alphas = np.append(alpha, alpha[-1])
    betas = np.append(beta,beta[-1])

    # define the linear ODE describing the mean degree evolution
    def dkdt(t, k, a_, b_):
        return 2 * a_ - 2*k * (b_)

    # intialize the integrator
    r = ode(dkdt)
    r.set_integrator('dopri5')

    k = [k0]
    new_t = [t[0]]

    i = 0
    for t_, a_, b_ in zip(t[:-1], alphas[:-1], betas[:-1]):

        # for every interval set the initial mean degree
        # and new parameters
        r.set_initial_value(k[-1], t_)
        r.set_f_params(a_, b_)

        # for 10 time points within this interval,
        # compute the expected mean degree
        this_t = np.linspace(t_, t[i+1], 10)
        for t__ in this_t[1:]:
            new_t.append(t__)
            this_k = r.integrate(t__)
            k.append(this_k)

        # increase interval
        i += 1

    return np.array(new_t), np.array(k)



def flockwork_P_mean_group_size_distribution_for_varying_rates(flockwork_P_params, N=None):
    """Compute the mean group size distribution for a Flockwork-P system with varying rates.

    Parameters
    ----------
    flockwork_P_params : :obj:`dict`
        Contains all parameters necessary for a Flockwork-P simulation, especially the 
        time-dependent rewiring rate and time-dependent reconnection probability
    N : int, default : None
        If given, compute everything for `N` nodes, where `N` is different from 
        `N` in `flockwork_P_params`.

    Returns
    -------
    mean_distribution : numpy.array
        An array of length `N` with its `i`-th entry containing the mean number of
        groups of size `m = i + 1`.
    """

    if N is None:
        N = flockwork_P_params['N']

    # estimate mean degree from integrating ODE
    new_t, k = flockwork_P_mean_degree_for_varying_rates(flockwork_P_params, N)

    # from equilibrium assumption k = P/(1-P) compute adjusted P
    new_P = k / (k+1)

    distro = []

    # for every time point and adjusted P, compute the equilibrium group size distribution
    for P_ in new_P:
        this_distro = flockwork_P_equilibrium_group_size_distribution(N, P_)
        distro.append(this_distro[1:])

    # compute the mean group size distribution as a time integral over the
    # group size distribution
    distro = np.array(distro)
    mean_distro = np.trapz(distro, x=new_t, axis=0) / (new_t[-1] - new_t[0])

    return mean_distro


def estimated_mean_group_size_distribution(temporal_network):
    """Compute the mean group size distribution for a temporal network under the assumption
    that it can be described reasonably by a flockwork-P model.

    Parameters
    ----------
    temporal_network : :mod:`edge_changes` or :mod:`edge_lists`
        A temporal network.

    Returns
    -------
    mean_distribution : numpy.array
        The average group size distribution of the temporal network which is closer to
        to the _true_ group size distribution than measuring over the binned system.
        The result is an array of length `N` with its `i`-th entry containing the mean number of
        groups of size `m = i + 1`.
    """

    new_t, k = tc.mean_degree(temporal_network)
    N = temporal_network.N

    # from equilibrium assumption k = P/(1-P) compute adjusted P
    new_P = k / (k+1)

    distro = []

    # for every time point and adjusted P, compute the equilibrium group size distribution
    for P_ in new_P:
        this_distro = flockwork_P_equilibrium_group_size_distribution(N, P_)
        distro.append(this_distro[1:])

    # compute the mean group size distribution as a time integral over the
    # group size distribution
    distro = np.array(distro)
    mean_distro = np.trapz(distro, x=new_t, axis=0) / (new_t[-1] - new_t[0])

    return mean_distro


def flockwork_P_mean_group_size_distribution_from_mean_degree_distribution(flockwork_P_params, dk, N=None):
    r"""Compute the mean group size distribution for a Flockwork-P system with varying rates from the 
    mean degree distribution which is fitted as :math:`\left\langle k\right\rangle^{-\alpha}`,
    hence this returns

    .. math::

        \left\langle N_m \right\rangle = \int dk P(k)  \times N_m( k/(k+1) )

    Parameters
    ----------
    flockwork_P_params : :obj:`dict`
        Contains all parameters necessary for a Flockwork-P simulation, especially the 
        time-dependent rewiring rate and time-dependent reconnection probability
    dk : float
        resolution of <k>-space for the solution of the integral.
    N : int, default : None
        If given, compute everything for `N` nodes, where `N` is different from 
        `N` in `flockwork_P_params`.

    Returns
    -------
    mean_distribution : numpy.array
        An array of length `N` with its `i`-th entry containing the mean number of
        groups of size `m = i + 1`.
    """

    if N is None:
        N = flockwork_P_params['N']

    # estimate mean degree from integrating ODE
    new_t, k = flockwork_P_mean_degree_for_varying_rates(flockwork_P_params, N)

    kmin = 2.0 / N
    ind = np.where(k >= kmin)
    new_t = new_t[ind]
    k = k[ind]

    alpha, err, xmin = fit_power_law_clauset(k)
    kmin = k.min()
    kmax = k.max()

    norm = 1/(1-alpha) * (kmax**(-alpha+1) - kmin**(-alpha+1))

    def dist(k_): return k_**(-alpha) / norm

    k = np.linspace(kmin, kmax, int((kmax-kmin) / dk) + 1)

    # from equilibrium assumption k = P/(1-P) compute adjusted P
    new_P = k / (k+1)

    distro = []

    # for every time point and adjusted P, compute the equilibrium group size distribution
    for P_ in new_P:
        this_distro = flockwork_P_equilibrium_group_size_distribution(N, P_)
        distro.append(this_distro[1:])

    distro = np.array(distro)

    mean_distro = simps(dist(k)[:, None] * distro, x=k, axis=0)

    return mean_distro


def naive_varying_rate_flockwork_simulation(N, t, reconnection_rates, disconnection_rates, tmax):
    r"""
    Do a naive simulation of a Flockwork systems where the reconnection and disconnection rate
    vary over time as step functions. It is called `naive` because the rate change will not
    be considered when evaluating the inter-event time at time points when the rates change.
    I.e. at a rate change at time `t`, the new Flockwork model will be initiated as if the 
    last event happened at time `t`. This error is neglibile if the upper bounded mean inter-event time
    at rate changes is 
    :math:`\tau=(N(\alpha_{\mathrm{min}}+\beta_{\mathrm{min}})^{-1}\ll \Deltat_\mathrm{min}`
    where :math:`\Deltat_\mathrm{min}` is the minimal time between two rate changes.

    Parameters
    ----------
    N : int
        number of nodes
    t : numpy.ndarray of float
        time points at which :math:`\alpha` and :math`\beta` change
    reconnection_rates : numpy.ndarray of float
        values of :math:`\alpha` for the corresponding time values in ``t``.
    disconnection_rates : numpy.ndarray of float
        values of :math:`\beta` for the corresponding time values in ``t``.
    tmax : float
        time at which the simulation is supposed to end.

    Returns
    -------
    edge_changes : :class:`_tacoma.edge_changes`
        The simulated Flockwork instance.
    """

    if len(t) != len(reconnection_rates) or len(reconnection_rates) != len(disconnection_rates):
        raise ValueError('The arrays `t`, `reconnection_rates` and `disconnection_rates` must have the same length' )

    assert(t[-1]<tmax)
    t = np.append(t,tmax)

    all_networks = []
    it = 0

    result = None
    for a_, b_ in zip(reconnection_rates, disconnection_rates):

        dt = t[it+1] - t[it]
        if a_ == 0.0 and b_ == 0.0:
            P = 0
        else:
            P = a_ / (a_+b_)
        g = (a_+b_)

        if it == 0:
            initial_edges = flockwork_P_equilibrium_configuration(N,P)

        FW = tc.FlockworkPModel(initial_edges, N, g, P, save_temporal_network = True)
        FW.simulate(dt)
        this_result = FW.edge_changes

        if it < len(reconnection_rates) - 1:
            initial_edges = FW.get_current_edgelist()

        all_networks.append(this_result)

        it += 1

    result = tc.concatenate(all_networks)

    return result
        


def flockwork_P(N, P, t_run_total, rewiring_rate = 1.0, initial_edges = None, seed = 0, return_edge_changes_with_histograms=False):
    r"""
    Simulate a flockwork P-model where the disconnection rate per 
    node is :math:`\gamma` and reconnection probability is :math:`P`.
    In order to start with an equilibrated initial state, use
    :func:`tacoma.flockwork.flockwork_P_equilibrium_configuration`
    or just pass ``initial_edges = None`` to this function.
    The total event rate is :math:`N\gamma` such that in
    one unit of time there's an expected number of
    :math:`N` events.

    Parameters
    ----------
    N : int
        number of nodes
    P : float
        The reconnection probability. Has to be :math:`0\leq P\leq 1`
    t_run_total : float
        The total run time in units of :math:`\gamma^{-1}`.
    rewiring_rate : float, default : 1.0
        Event rate :math:`\gamma` per node per unit of time.
    initial_edges : list of tuple of int, default : None
        The initial state of the network as an edge list.
        If `None` is provided, the initial state will be taken
        from an equilibrium configuration generated with
        :func:`tacoma.flockwork.flockwork_P_equilibrium_configuration`
    seed : int, default : 0
        The random seed.
    return_edge_changes_with_histograms : bool, default : False
        Instead of the converted :class:`_tacoma.edge_changes`,
        return the original instance of 
        :class:`_tacoma.edge_changes_with_histograms`.

    Returns
    -------
    :class:`_tacoma.edge_changes`
        The simulated network. if return_edge_changes_with_histograms is `True`,
        returns an instance of :class:`_tacoma.edge_changes_with_histograms` instead.

    """

    if initial_edges is None:
        initial_edges = flockwork_P_equilibrium_configuration(N, P)

    FW = tc.FlockworkPModel(initial_edges, N, rewiring_rate, P, save_temporal_network = True)
    FW.simulate(dt)
    this_result = FW.edge_changes

    return this_result

def degree_distribution(N,P):
    """Get the equilibrium degree distribution of a Flockwork-P model
    given node number N and probability to reconnect P.

    Parameters
    ----------
    N : int
        Number of nodes
    P : float
        Probability to reconnect

    Returns
    -------
    numpy.ndarray
        Degree distribution of this configuration. The :math:`k`-th entry
        of this array is the probability that a node has degree :math:`k`.
    """

    C_m = flockwork_P_equilibrium_group_size_distribution(N, P)
    P_k = np.arange(1,N+1) * C_m[1:] / N

    return P_k

def degree_moment(N,P,m):
    r"""Get the :math:`m`-th moment of the degree of an Flockwork-P model
    equilibrium configuration
    given node number N and probability to reconnect P.

    Parameters
    ----------
    N : int
        Number of nodes
    P : float
        Probability to reconnect

    Returns
    -------
    <k^m> : float
        The :math:`m`-th moment of the degree distribution.
    """

    P_k = degree_distribution(N, P)

    return (np.arange(N)**m).dot(P_k)

def mean_degree(N,P,m):
    r"""
    Get the exact theoretical mean degree :math:`\left< k\right>`
    for a Flockwork-P model.

    Parameters
    ----------
    N : int
        Number of nodes
    P : float
        Probability to reconnect

    Returns
    -------
    <k> : float
        The mean degree.
    """

    return degree_moment(N, P, 1)

def convert_to_edge_activity_parameters_plus_minus(N, P, gamma=1.0):
    r"""
    Convert the Flockwork-P parameters math:`P` and :math:`\gamma` 
    to the corresponding parameters in the edge activity model
    :math:`\omega^{+}` and :math:`\omega^{-}`.

    Parameters
    ----------
    N : int
        Number of nodes
    P : float
        Probability to reconnect
    gamma : float, default = 1
        Rewiring rate.

    Returns
    -------
    omega_plus : float
        The rate with which inactive edges are activated.
    omega_minus : float
        The rate with which active edges are deactivated.
    """

    k = degree_moment(N, P, 1)
    k2 = degree_moment(N, P, 2)

    omega_minus = 2*gamma*(1-P/(N-1)*k2/k)
    rho = k / (N-1)
    omega = omega_minus * rho
    omega_plus = omega / (1-rho)

    return omega_plus, omega_minus
    
def convert_to_edge_activity_parameters(N, P, gamma=1.0):
    r"""
    Convert the Flockwork-P parameters math:`P` and :math:`\gamma` 
    to the corresponding parameters in the edge activity model
    :math:`\rho` and :math:`\omega`.

    Parameters
    ----------
    N : int
        Number of nodes
    P : float
        Probability to reconnect
    gamma : float, default = 1
        Rewiring rate.

    Returns
    -------
    rho : float
        The network density.
    omega : float
        The rate with which either 
    """

    k = degree_moment(N, P, 1)
    k2 = degree_moment(N, P, 2)

    omega_minus = 2*gamma*(1-P/(N-1)*k2/k)
    rho = k / (N-1)
    omega = omega_minus * rho

    return rho, omega
    

if __name__ == "__main__":

    N = 10
    P = 0.5

    dist = flockwork_P_equilibrium_group_size_distribution(N, P)

    print(dist, sum([m * h for m, h in enumerate(dist)]))

