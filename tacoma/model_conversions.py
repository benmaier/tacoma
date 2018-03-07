# -*- coding: utf-8 -*-

import numpy as np

from scipy.optimize import fsolve
from scipy.optimize import minimize

from tacoma import edge_changes as ec
from tacoma import edge_lists as el
from tacoma import edge_lists_with_histograms as el_h
from tacoma import edge_changes_with_histograms as ec_h
from tacoma.analysis import get_logarithmic_histogram
from tacoma.power_law_fitting import fit_power_law_clauset
from tacoma import _get_raw_temporal_network

def ZSBB_mean_coordination_number(b0,lam,N,b1):

    s = 0.0
    for n in range(1,N):
        s += (n+1.0) / ((n+1.0)*b1 - 1.0) * ((1.0-lam)/lam)**(n-1)

    pi_10 = ( 0.5 / ( b0 - (2*lam - 1.0)/(3*lam - 1.0)) + 0.5/lam * s )**(-1.0)

    return pi_10 / (2.0*lam) * s

def b0_func(b0,lam):
    return b0 * (2.0+(1.0-lam/(2*lam-1))) + 1


def get_mean_coordination_number(group_sizes_and_durations):


    group_size_histogram = np.array([
                            (size, val)\
                            for size, val in enumerate(group_sizes_and_durations.aggregated_size_histogram)\
                            if size > 0
                        ],dtype=float)

    N = group_size_histogram.shape[0]
    m = group_size_histogram[:,0]   # group sizes
    N_m = group_size_histogram[:,1] # mean number of groups existing

    P_k = m * N_m / N # degree distribution
    
    return m.dot(P_k) - 1.0


def estimate_ZSBB_parameters(temporal_network,
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

    mean_n = get_mean_coordination_number(result)

    def equations(p):
        b0, lam = p
        n = ZSBB_mean_coordination_number(b0,lam,N,b1)
        al0 = b0_func(b0,lam)
        return (n - mean_n, al0 - alpha_0)

    def cost(p):
        b0, lam = p
        n = ZSBB_mean_coordination_number(b0,lam,N,b1)
        al0 = b0_func(b0,lam)
        return np.abs((n - mean_n)/mean_n) + np.abs((al0 - alpha_0)/alpha_0)

    #b0, lam = fsolve(equations,(0.7,0.7))
    res = minimize(cost,(0.7,0.7))
    b0, lam = res.x

    if lam < 0.5:
        lam = 0.51
    elif lam > 1:
        lam = 1

    if b0 < 0.5:
        b0 = 0.51
    elif b0 > 1:
        b0 = 1

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
    kwargs['t_equilibration'] = 10000*N

    return kwargs

if __name__ == "__main__":
    import tacoma as tc
    import matplotlib.pyplot as pl
    from tacoma.analysis import temporal_network_group_analysis
    

    test = tc.dtu_week()
    rewiring_rate = test.gamma
    P = test.P

    fw = tc.flockwork_P_varying_rates([],100,P,24*3600,rewiring_rate,tmax=24*3600*7)
    fw_binned = tc.sample_temporal_network(fw,dt=300)
    fw_binned_result = tc.measure_group_sizes_and_durations(fw_binned)

    kwargs = estimate_ZSBB_parameters(fw_binned,fw_binned_result,fit_discrete=True,dt=300.)
    print "lambda =", kwargs['lambda']
    print "b0 =", kwargs['b0']
    print "b1 =", kwargs['b1']
    kwargs['t_run_total'] = (len(fw_binned.t) + 1)*kwargs['N']
    zsbb = tc.ZSBB_model(**kwargs)
    zsbb_binned = tc.sample_temporal_network(zsbb,dt=kwargs['N'])
    zsbb_binned_result = tc.measure_group_sizes_and_durations(zsbb_binned)


    fig, ax, data = temporal_network_group_analysis(fw_binned_result)
    temporal_network_group_analysis(zsbb_binned_result,
                                    time_normalization_factor = 300./kwargs['N'],
                                    ax=ax)
    pl.show()
