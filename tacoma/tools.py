from math import log10
import json

import numpy as np

import tacoma as tc

from _tacoma import edge_changes as ec
from _tacoma import edge_lists as el
from _tacoma import edge_lists_with_histograms as el_h
from _tacoma import edge_changes_with_histograms as ec_h


from scipy.optimize import curve_fit
from scipy.stats import lognorm
from scipy.special import gamma as Gamma
from scipy.stats import weibull_min

from scipy.integrate import quad

from lmfit import minimize, Parameters

def complete_graph(N,tmax=1.0):
    """Get a single frame which consists of a complete network.

    Parameters
    ----------
    N : int
        Number of nodes.
    tmax : float, default : 1.0
        Maximum time of the static network.
        
    Returns 
    -------
    :class:`_tacoma.edge_lists`
        An instance of `tacoma.edge_lists` with t = [0.0], tmax = ``tmax``.
    """

    edge_list = []
    for i in range(N-1):
        for j in range(i+1,N):
            edge_list.append((i,j))

    this = tc.edge_lists()
    this.t = [0.]
    this.tmax = tmax
    this.edges = [edge_list]
    this.N = N

    return this

def convert_static_network(N,edge_list,tmax=1.0):
    """Get a single frame which consists of the static network.

    Parameters
    ----------
    N : int
        Number of nodes.
    edge_list : :obj:`list` of :obj:`tuple` of int
        The edges of the static graph
    tmax : double, default : 1.0
        The maximum time until the network is looped.
        
    Returns 
    -------
    :class:`_tacoma.edge_lists`
        An instance of `tacoma.edge_lists` with t = [0.0], tmax = ``tmax``.
    """

    new_edge_list = []

    for e in edge_list:
        if e[0] > e[1]:
            new_edge_list.append((e[1],e[0]))
        elif e[1] > e[0]:
            new_edge_list.append((e[0],e[1]))

    this = tc.edge_lists()
    this.t = [0.]
    this.tmax = tmax
    this.edges = [ new_edge_list ]
    this.N = N

    return this

def get_logarithmic_histogram(data,
                              bins, #number of bins
                              return_bin_means = True,
                              density = True):
    """Get a logarithmic histogram.

    Parameters
    ----------
    data : array-like
        The data to bin.
    bins : int
        The number of bins.
    return_bin_means : bool, default : True
        return the geometric means of binning intervals,
        otherwise return bin edges
    density : bool, default : True
        return probability density, if `False`, return counts

    Returns
    -------
    x : numpy.ndarray
        Either geometric intervals or interval edges
    y : numpy.ndarray
        Either probability density or counts.
    """
    
    data = np.array(data)
    data = data[data>0]
    MIN = min(data)
    MAX = max(data)

    # check if bins is an array, if not, make an array
    try:
        bins[0]
    except:
        # if this fails, we assume `bins` is an integer
        bins = np.logspace(log10(MIN), log10(MAX), bins+1,base=10.)

    y, x = np.histogram(data,
                        bins = bins,
                        density = density,
                       )
    
    new_x = np.sqrt(x[1:]*x[:-1])

    if return_bin_means:
        return new_x, y
    else:
        return x, y

def group_size_histogram(group_sizes_and_durations):
    """Returns the mean number of groups of size `g`
    (mean over both ensemble and time).

    Parameters
    ==========
    group_sizes_and_durations : :mod:`group_sizes_and_durations`
        Result from the function :mod:`measure_group_sizes_and_durations`

    Returns
    =======
    g : :obj:`numpy.ndarray` of float
        Group sizes
    N : :obj:`numpy.ndarray` of float
        Mean number of groups of the corresponding size in `g`.
    """

    group_size_histogram = np.array([
                            (size, val)\
                            for size, val in enumerate(group_sizes_and_durations.aggregated_size_histogram)\
                            if val > 0.
                        ],dtype=float)

    x_group, y_group = group_size_histogram[:,0], group_size_histogram[:,1]
 
    return x_group, y_group


def mean_coordination_number(group_sizes_and_durations):
    """Returns the mean coordination number (mean over both ensemble and time).
    Following the definition by Zhao, Stehle, Bianconi, Barrat,
    the coordination number of node `i` is equal to the size of the group
    it is part of.

    Parameters
    ==========
    group_sizes_and_durations : :mod:`group_sizes_and_durations`
        Result from the function :mod:`measure_group_sizes_and_durations`

    Returns
    =======
    mean_coordination_number : float
        Temporal and ensemble mean of a node's group size.
    """

    m, N_m = group_size_histogram(group_sizes_and_durations)

    N = len( group_sizes_and_durations.aggregated_size_histogram) - 1

    P_k = m * N_m / N # degree distribution for completely connected groups
    
    return m.dot(P_k) - 1.0


def mean_number_of_groups(group_sizes_and_durations):
    """Returns the mean number of groups (mean over both ensemble and time).

    Parameters
    ==========
    group_sizes_and_durations : :mod:`group_sizes_and_durations`
        Result from the function :mod:`measure_group_sizes_and_durations`

    Returns
    =======
    mean_number_of_groups : float
        Temporal and ensemble mean of the total number of groups
        a network consists of at a certain time.
    """

    m, N_m = group_size_histogram(group_sizes_and_durations)

    return N_m.sum()


def mean_group_size(group_sizes_and_durations):
    """Returns the mean group size (mean over both ensemble and time).

    Parameters
    ==========
    group_sizes_and_durations : :mod:`group_sizes_and_durations`
        Result from the function :mod:`measure_group_sizes_and_durations`

    Returns
    =======
    mean_group_size : float
        Temporal and ensemble mean of the group size of all groups
        a network consists of at a certain time.
    """


    N = len( group_sizes_and_durations.aggregated_size_histogram) - 1
    c = mean_number_of_groups(group_sizes_and_durations)

    return float(N) / c

def slow_mean_degree(temporal_network):
    """Returns the mean degree (mean over ensemble) but it takes
    ages to compute. You should instead use :mod:`mean_degree`.

    Parameters
    ==========
    temporal_network : :mod:`edge_lists` or :mod:`edge_changes`

    Returns
    =======
    t : float
        Temporal and ensemble mean of the node degree.
    mean_degree : float
        Temporal and ensemble mean of the node degree.
    """


    temporal_network = tc._get_raw_temporal_network(temporal_network)

    if type(temporal_network) == tc.ec:
        t = np.array([temporal_network.t0])
        t = np.append(t,temporal_network.t)
        m = [ len(temporal_network.edges_initial) ]

        for i in range(len(temporal_network.edges_in)):
            m_in = len(temporal_network.edges_in[i])
            m_out = len(temporal_network.edges_out[i])
            dm = m_in - m_out
            m.append( m[-1] + dm )

        m.append(m[-1])
        t = np.append(t,temporal_network.tmax)

        m = np.array(m, dtype=float)
        k = 2.0*m / float(temporal_network.N)

    elif type(temporal_network) == tc.el:
        t = np.array(temporal_network.t)
        m = [ ]

        for i in range(len(temporal_network.edges)):
            _m = len(temporal_network.edges[i])
            m.append( _m )

        m.append(m[-1])
        t = np.append(t,temporal_network.tmax)

        m = np.array(m, dtype=float)
        k = 2.0*m / float(temporal_network.N)

    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))


    return t, k

def time_average(t,x,tmax=None):
    """
    Return a temporal average of an observable `x(t)`, 

    .. math::

        \\overline x = \\frac{1}{t_\\mathrm{max}-t_0} 
                       \\int\\limits_{t_0}^{t_\\mathrm{max}} dt\, x(t).

    Where it's assumed that `x` changes as a step function, such that

    .. math::

        \\overline x = \\frac{1}{t_\\mathrm{max}-t_0} 
                       \\left( (t_\\mathrm{max} - t_{N_t-1})  x(t_{N_t-1}) +
                                \\sum_{i=1}^{N_t-1} (t_{i}-t_{i-1})x(t_{i-1}) \\right).

    Parameters
    ----------
    t : numpy.ndarray
        Times at which `x` changes
    x : numpy.ndarray
        Value of the observable at the corresponding times.
    tmax : float, default : None
        If this is None, the integral is computed until 
        time t[-1], otherwise it's evaluated until tmax.

    Returns
    -------
    average_x : float
        The temporal average of `x`.
    """
    
    if len(t) != len(x):
        raise ValueError("t and x must have the same shape")

    if tmax is not None:
        t = np.append(t,tmax)
        x = np.append(x,x[-1])


    dt = t[1:] - t[:-1]

    sum_dt = dt.sum()
    
    return dt.dot(x[:-1]) / sum_dt

def time_RMSE(t,x1,x2,tmax=None):
    r"""
    Get the root mean squared error (RMSE) of two observables
    over the same time. Combine with :func:`tacoma.tools.sample_a_function`.

    .. math::

        \mathrm{RMSE} = \left( 
                            \frac{1}{t_\mathrm{max}-t_0}
                            \int\limits_{t_0}^{t_\mathrm{max}} dt\, 
                                \left( x_1(t) - x_2(t) \right)^2.
                        \right)^{1/2}

    Parameters
    ----------
    t : numpy.ndarray
        times at which ``x1`` and ``x2`` change
    x1 : numpy.ndarray
        The first observable to compute the RMSE.
    x2 : numpy.ndarray
        The second observable to compute the RMSE.
    tmax : float, default : None
        If this is None, the integral is computed until 
        time t[-1], otherwise it's evaluated until tmax.

    Returns
    -------
    RMSE : float
        The root mean squared error.
    """
    
    if len(t) != len(x1):
        raise ValueError("t and x1 must have the same shape")

    if len(t) != len(x2):
        raise ValueError("t and x2 must have the same shape")

    if tmax is not None:
        t = np.append(t,tmax)
        x1 = np.append(x1,x1[-1])
        x2 = np.append(x2,x2[-1])

    return np.sqrt(time_average(t,(x1-x2)**2))

def bin_a_function(x,y,bins,mode='mean'):
    r"""
    Bin a step function :math:`y(t)` over time bins. Binning can be done using
    either ``numpy.sum`` or ``numpy.mean``.

    Parameters
    ----------
    x : numpy.ndarray
        domain values at which `y` changes
    y : numpy.ndarray
        values of the cunt
    bins : numpy.ndarray
        Bin to those bin edges
    mode : string, default : mean
        Build either a ``mean`` over the bins or a ``sum``.

    Returns
    -------
    new_y : numpy.ndarray
        binned observable
    """
    if mode=='mean':
        cumfunc = np.mean
    elif mode=='sum':
        cumfunc = np.sum
    indices = np.searchsorted(x,bins)
    new_y = [0]
    for i in range(1,len(indices)):
        if indices[i-1] == indices[i]:
            _y = 0
        else:
            _y = cumfunc(y[indices[i-1]:indices[i]])
        new_y.append( _y )

    return new_y

def sample_a_function(x,y,time_points,sample_width=0):
    r"""
    Sample an observable :math:`y` which is a step function
    and changes at corresponding values :math:`t`.

    Parameters
    ----------
    t : numpy.ndarray
        domain values at which `y` changes
    y : numpy.ndarray
        observable to sample
    time_points : numpy.ndarray
        time points for which to sample.
    sample_width : float, default : 0.0
        sample as a temporal average over [t-sample_width/2, t+sample_width/2].

    Returns
    -------
    new_y : numpy.ndarray
        values of `y` sampled at `time_points`.
    """

    x = t
    new_y = []
    for bin in time_points:
        if sample_width > 0:
            indices = np.searchsorted(x,[bin-sample_width/2., bin+sample_width/2.])
            dt = x[(indices[0]+1):indices[1]] - x[indices[0]:(indices[1]-1)]
            y_ = y[(indices[0]+1):indices[1]]
            val = np.dot(dt,y_)
            norm = np.sum(dt)
            #x_ = x[(indices[0]+1):indices[1]]
        #y_ = y[indices[0]:indices[1]]
        #val = np.trapz(y=y_,x=x_)
        #norm = np.trapz(y=np.ones_like(x_),x=x_)
            new_y.append( val/norm)
        else:
            indices = np.searchsorted(x,[bin])
            if indices[0] == 0:
                this_index = 0
            else:
                this_index = indices[0]-1
            new_y.append(y[this_index])

    return np.array(new_y)


def rescale_time(temporal_network, new_t0, new_tmax):
    """Rescale the time in this temporal network (inplace).
    
    Parameters
    ==========
    temporal_network : :mod:`edge_lists` or :mod:`edge_changes`
    new_t0 : float
        The new value of t0.
    new_tmax : float
        The new value of tmax.

    Returns
    =======
    temporal_network : :mod:`edge_lists` or :mod:`edge_changes`
        Same instance as input.
    """

    if hasattr(temporal_network,'t0'):
        this_t0 = temporal_network.t0
        temporal_network.t0 = new_t0
    else:
        this_t0 = temporal_network.t[0]

    this_T = temporal_network.tmax - this_t0
    new_T = new_tmax - new_t0
    
    temporal_network.t = [ (t - this_t0) / this_T * new_T + new_t0 for t in temporal_network.t ]
    temporal_network.tmax = new_tmax

    return temporal_network

def number_of_discovered_edges(temporal_network):
    """Get the total number of discovered unique edges C(t), i.e. the contact coverage.
    
    Parameters
    ==========
    temporal_network : :class:`_tacoma.edge_trajectories`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes` or :obj:`list` of :class:`_tacoma.edge_trajectory_entry`

    Returns
    =======
    t : numpy.ndarray
        Time points at which new edges have been discovered
    C : numpy.ndarray
        total number of edges discovered up to time t.
    """

    if type(temporal_network) in [ ec, el, el_h, ec_h ]:
        traj = tc.get_edge_trajectories(temporal_network)
    elif type(temporal_network) == list and type(temporal_network[0]) == tc.edge_trajectory_entry:
        traj = temporal_network
    elif type(temporal_network) == edge_trajectories:
        traj = temporal_network.trajectories
    else:
        raise ValueError("Unknown type for temporal network:", type(temporal_network))

    t = []
    count = []
    for iedge, entry in enumerate(traj):
        t0 = entry.time_pairs[0][0]
        if len(t)>0:
            if t[-1] == t0:
                count[-1] = iedge+1
            else:
                count.append(iedge+1)
                t.append(t0)
        else:
            count.append(iedge+1)
            t.append(t0)

    return np.array(t), np.array(count,dtype=float)

def get_edge_probability_and_rate(temporal_network):
    """
    For each edge compute the probability that it is active and the rate
    with which it is activated.
    
    Parameters
    ==========
    temporal_network : :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes` or :class:`_tacoma.edge_trajectories`

    Returns
    =======
    p : numpy.ndarray
        The probability to be switched on for each observed edge of the network 
        (the remaining un-observed edges have probability p = 0).
    omega : numpy.ndarray
        The rate with which the observed edges are switched on 
        :math:`\\omega = \\left(\\frac{1}{\\tau^+} + \\frac{1}{\\tau^{-}}\\right)^{-1}`
        (the remaining un-observed edges have rate :math:`\\omega = 0`).
    """

    if type(temporal_network) in [ ec, el, el_h, ec_h ]:
        traj = tc.get_edge_trajectories(temporal_network)
        tmax = temporal_network.tmax
        if type(temporal_network) in [ ec, ec_h ]:
            t0 = temporal_network.t0
        else:
            t0 = temporal_network.t[0]
    elif type(temporal_network) == edge_trajectories:
        traj = temporal_network.trajectories
        tmax = temporal_network.tmax
        t0 = temporal_network.t0
    else:
        raise ValueError("Unknown type for temporal network:", type(temporal_network))

    T = tmax - t0

    connection_probability = np.empty(len(traj))
    activity_rate = np.empty(len(traj))

    for iedge, entry in enumerate(traj):

        t_on = 0.0
        for interval in entry.time_pairs:
            t_on += interval[1] - interval[0]

        activity_rate[iedge] = len(entry.time_pairs) / T
        connection_probability[iedge] = t_on / T

    return connection_probability, activity_rate

def get_reduced_time(t, intervals_to_discard):
    """
    Reduce the provided time in a way that intervals
    are cut out.

    Parameters
    ----------
    t : numpy.ndarray
        time array to reduce
    intervals_to_discard : list of tuple of float
        The time intervals which have to be discarded.

    Returns
    -------
    new_t : numpy.ndarray
        reduced time with equal shape to ``t``. Each time in
        intervals which have to be discarded are remapped
        to the beginning of the intervals.
    """

    x_ = t.copy()

    offset = 0.0

    for interval in intervals_to_discard_for_fit:
        t0, t1 = interval
        x_[np.logical_and(x>=t0, x<t1)] = t0 - offset

        x_[x>=t1] -= t1 - t0
        offset += t1 - t0

    return x_

def fit_contact_coverage(N, time, contact_coverage, intervals_to_discard_for_fit=[], kind='gamma'):
    """
    Fit the contact coverage :math:`C(t)` to an edge activity rate distribution of kind ``kind``.

    Parameters
    ----------
    N : int
        Number of nodes.
    time : numpy.ndarray
        Time points at which the contact converage changed.
    contact_coverage : numpy.ndarray
        Contact coverage at the corresponding times in ``time``.
    intervals_to_discard_for_fit : list of tuple of float
        a list of time intervals which have to be discarded for the fit
    kind : string, default : 'gamma'
        Rate distribution model for fit. Can be ``gamma``, ``exponential``,
        ``uniform``, ``normal``, ``delta``, ``power-law``.

    Returns
    -------
    fit_function : function
        The fit function
    popt : tuple of float
        The found optimal parameters to pass to ``fit_function``
    pstd : tuple of float
        standard uncertainty of the parameters
    """


    fac = N*(N-1)/2.
    edge_count = contact_coverage


    if kind == 'gamma':
        fit = lambda x, alpha, scale: fac * (1 - (scale/(scale + get_reduced_time(x, intervals_to_discard_for_fit)) )**alpha)
        popt, pcov = curve_fit(fit, time, edge_count,[0.5,10.0],maxfev=10000)
    elif kind == 'exponential':
        alpha = 1.0
        fit = lambda x, scale: fac * (1 - (scale/(scale + get_reduced_time(x, intervals_to_discard_for_fit)) )**alpha)
        popt, pcov = curve_fit(fit, time, edge_count,[10.0],maxfev=10000)
    elif kind == 'uniform':

        def fit(x, alpha):
            t = get_reduced_time(x, intervals_to_discard_for_fit)
            result = fac * (1 - (1-np.exp(-alpha*t)) / (alpha * t))
            result[t==0] = 0
            return result

        popt, pcov = curve_fit(fit, time, edge_count,[0.1],maxfev=10000)
    elif kind == 'normal':
        fit = lambda x, mean, sigma: fac * (1 - np.exp(-mean* get_reduced_time(x, intervals_to_discard_for_fit)+sigma**2/2.0*get_reduced_time(x, intervals_to_discard_for_fit)**2))
        popt, pcov = curve_fit(fit, time, edge_count,[0.5,0.01],maxfev=10000)
    elif kind == 'delta':
        fit = lambda x, scale: fac * (1 - np.exp(-scale * get_reduced_time(x, intervals_to_discard_for_fit)))
        popt, pcov = curve_fit(fit, time, edge_count,[0.5],maxfev=10000)
    elif kind == 'power-law':
        def fit(x,alpha, xmin,xmax,time_is_reduced = False):
            #print(alpha, xmin, xmax)
            A = (1-alpha) / (xmax**(1-alpha) - xmin**(1-alpha))
            if time_is_reduced:
                t = x
            else:
                t = get_reduced_time(x, intervals_to_discard_for_fit)
            x_ = np.logspace(np.math.log10(xmin),np.math.log10(xmax),100)

            #observable = x_[None,:]**(-alpha)*np.exp(-t[:,None]*x_[None,:])
            #print(observable.shape)
            #M = np.trapz(observable, x = x_[None,:], axis=0,)
            #print(M.shape)
            M = np.array([ np.trapz(x_**(-alpha)*np.exp(-t_*x_), x = x_) for t_ in t])
            #M = np.array([ quad(lambda x: A * x**(-alpha)*np.exp(-t_*x), xmin,xmax)[0] for t_ in t ])
            return fac * (1-M)

        def residual(params, x, time_is_reduced=True):
            alpha, xmin, xmax = params['alpha'].value, params['xmin'].value, params['xmax'].value
            #print(alpha, xmin, xmax)
            return fit(x, alpha, xmin, xmax,time_is_reduced) - edge_count
        reduced_time = get_reduced_time(time, intervals_to_discard_for_fit)

        params = Parameters()
        params.add('alpha',value = 1.0, min=1e-16)
        params.add('xmin',value = 1e-5, min=1e-16)
        params.add('xmax',value = 10, min=1e-16)
        out = minimize(residual, params, args = (reduced_time,))

        popt = (out.params['alpha'].value, out.params['xmin'].value, out.params['xmax'].value)

        #print(help(out))
        #print(out.pretty_print())

        # TODO: once lmfit is fixed, this has to change to out.covar
        pcov = np.zeros(len(popt))

        #popt, pcov = curve_fit(fit, time, edge_count,[1,1e-10,10],maxfev=10000)
#    elif kind == 'weibull':
#        #def fit(x, k, omega_c):
#        #    nmax = 18
#        #    moments = np.array([ omega_c**n * Gamma(1+n/k) for n in range(nmax+1) ] )
#        #    def M(t_,nmax):
#        #        return np.sum(np.array([ t_**n * moments[n] / np.math.factorial(n) for n in range(nmax+1)]), axis=0)
#        #        
#        #    return fac * (1-M(-t,nmax))
#        def fit(t_, k, omega_c):
#            eps = 1e-16
#            integrand = lambda x, t : weibull_min(k,scale=omega_c).pdf(x) * np.exp(-t*x)
#            t = get_reduced_time(t_, intervals_to_discard_for_fit)
#            return np.array([quad(integrand,eps,np.inf,args=(t_,)) for this_t in t])
#
#        popt, pcov = curve_fit(fit, time, edge_count,[0.5, 0.01],maxfev=10000)
#    elif kind == 'lognormal':
#
#
#        def fit(x, mu, sigma):
#            nmax = 15
#            moments = np.array([ np.exp(n*mu + n**2*sigma**2 / 2.0 for n in range(nmax+1) ))
#            def M(t_,nmax):
#                pass
#                
#
#            t = get_reduced_time(x, intervals_to_discard_for_fit)
#
#            pdf = lognorm(sigma,scale=np.exp(mu)).pdf
#            #weights = lognorm.rvs(sigma,scale=np.exp(mu),size=N*(N-1)//2)
#            #print(weights)
#            return fac * ( 1.0 - np.array([np.mean(np.exp(-weights*x_)) for x_ in x]))
#
#        popt, pcov = curve_fit(fit, time, edge_count,[0.5,0.5],maxfev=10000)
    else:
        raise ValueError('Unknown fit function:', kind)
    #popt, pcov = curve_fit(fit, fit_x, fit_y,[1./fac,fac,10.0],maxfev=10000)
    #popt, pcov = curve_fit(fit, fit_x, fit_y,[2,fac,10.0],maxfev=10000)


    return fit, popt, np.sqrt(np.diag(pcov))


