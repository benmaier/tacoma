from math import log10
import json

import numpy as np

import tacoma as tc

from scipy.optimize import curve_fit
from scipy.stats import lognorm

def complete_graph(N):
    """Get a single fram which consists of a complete network.

    Parameters
    ----------
    N : int
        Number of nodes.
        
    Returns 
    -------
    :mod:`edge_lists`
        An instance of `tacoma.edge_lists` with t = [0.0], tmax = 1.0.
    """

    edge_list = []
    for i in range(N-1):
        for j in range(i+1,N):
            edge_list.append((i,j))

    this = tc.edge_lists()
    this.t = [0.]
    this.tmax = 1.
    this.edges = [edge_list]
    this.N = N

    return this

def convert_static_network(N,edge_list,tmax):
    """Get a single frame which consists of the static network.

    Parameters
    ----------
    N : int
        Number of nodes.
    edge_list : :obj:`list` of :obj:`tuple` of int
        The edges of the static graph
    tmax : double
        The maximum time until the network is looped.
        
    Returns 
    -------
    :mod:`edge_lists`
        An instance of `tacoma.edge_lists` with t = [0.0], tmax = 1.0.
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

    group_size_histogram = np.array([
                            (size, val)\
                            for size, val in enumerate(group_sizes_and_durations.aggregated_size_histogram)\
                            if val > 0.
                        ],dtype=float)

    x_group, y_group = group_size_histogram[:,0], group_size_histogram[:,1]
 
    return x_group, y_group


def mean_coordination_number(group_sizes_and_durations):


    m, N_m = group_size_histogram(group_sizes_and_durations)

    N = len( group_sizes_and_durations.aggregated_size_histogram) - 1

    P_k = m * N_m / N # degree distribution for completely connected groups
    
    return m.dot(P_k) - 1.0


def mean_number_of_groups(group_sizes_and_durations):

    m, N_m = group_size_histogram(group_sizes_and_durations)

    return N_m.sum()


def mean_group_size(group_sizes_and_durations):

    N = len( group_sizes_and_durations.aggregated_size_histogram) - 1
    c = mean_number_of_groups(group_sizes_and_durations)

    return float(N) / c

def slow_mean_degree(temporal_network):

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
    
    if len(t) != len(x):
        raise ValueError("t and x must have the same shape")

    if tmax is not None:
        t = np.append(t,tmax)
        x = np.append(x,x[-1])


    dt = t[1:] - t[:-1]

    sum_dt = dt.sum()
    
    return dt.dot(x[:-1]) / sum_dt

def time_RMSE(t,x1,x2,tmax=None):
    
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
    result = tc.get_edge_trajectories(temporal_network)
    traj = result.trajectories
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

def get_reduced_time(x, intervals_to_discard_for_fit):

    x_ = x.copy()

    offset = 0.0

    for interval in intervals_to_discard_for_fit:
        t0, t1 = interval
        x_[np.logical_and(x>=t0, x<t1)] = t0 - offset

        x_[x>=t1] -= t1 - t0
        offset += t1 - t0

    return x_

def fit_number_of_discovered_edges(N, time, edge_count, intervals_to_discard_for_fit=[],kind='gamma'):


    fac = N*(N-1)/2.

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
    elif kind == 'lognormal':
        def fit(x, mu, sigma):

            t = get_reduced_time(x, intervals_to_discard_for_fit)
            weights = lognorm.rvs(sigma,scale=np.exp(mu),size=N*(N-1)//2)
            print(weights)
            return fac * ( 1.0 - np.array([np.mean(np.exp(-weights*x_)) for x_ in x]))

        popt, pcov = curve_fit(fit, time, edge_count,[0.5,0.5],maxfev=10000)
    else:
        raise ValueError('Unknown fit function:', kind)
    #popt, pcov = curve_fit(fit, fit_x, fit_y,[1./fac,fac,10.0],maxfev=10000)
    #popt, pcov = curve_fit(fit, fit_x, fit_y,[2,fac,10.0],maxfev=10000)


    return fit, popt, np.sqrt(np.diag(pcov))


def load_json_dict(fn):
    with open(fn,'r') as f:
         this_dict = json.load(f)

    return this_dict
