from math import log10

import numpy as np

import tacoma as tc

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
