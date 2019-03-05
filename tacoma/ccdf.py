
import tacoma as tc

import numpy as np
from numpy import log10

def get_bin_means(bin_edges,logarithmic_bins=False):
    if logarithmic_bins:
        return np.sqrt(bin_edges[1:] * bin_edges[:-1])
    else:
        return 0.5 * (bin_edges[1:] + bin_edges[:-1])

def get_ccdf(data):

    x = np.append([data.min()-1.0], np.sort(data))
    y = 1 - np.arange(0,len(data)+1) / len(data)

    return x, y

def get_ccdf_from_distribution(x,y):

    new_x_min = x.min() - 1
    max_ndx = np.where(y>0)[0][-1] + 1

    CDF = np.append([0.],np.cumsum(y))
    CDF = (1 - CDF/CDF.max())[:max_ndx+1]

    x = np.append([new_x_min], x)[:max_ndx+1]

    return x, CDF


def compute_ccdfs(binned_temporal_network,max_group,time_normalization_factor=1./3600.,n_bins=50,logarithmic_bins=False):

    t_fw, k_fw = tc.mean_degree(binned_temporal_network)

    if logarithmic_bins:
        bins = np.append([0.],np.logspace(log10(k_fw[k_fw>0.0].min())-0.1,log10(k_fw.max()),n_bins) )
    else:
        bins = np.append([0.],np.linspace(k_fw[k_fw>0.0].min(), k_fw.max(), n_bins) ) 

    x_k, y_k = get_ccdf(k_fw)
    y_k = tc.sample_a_function(x_k, y_k, bins)
    x_k = bins

    result = tc.measure_group_sizes_and_durations(binned_temporal_network)

    grp_sizes = np.array(result.aggregated_size_histogram[1:])
    m = np.arange(1,len(grp_sizes)+1)
    m, grp_sizes = get_ccdf_from_distribution(m, grp_sizes)

    durations = np.array(result.contact_durations) * time_normalization_factor

    if logarithmic_bins:
        bins = np.append([0.],np.logspace(log10(durations.min())-0.1,log10(durations.max()),n_bins) )
    else:
        bins = np.append([0.],np.linspace(durations.min(), durations.max(), n_bins) )

    x_contact, y_contact = get_ccdf(durations)
    y_contact = tc.sample_a_function(x_contact, y_contact, bins)
    x_contact = bins

    y_groups = []
    x_groups = []

    for group_size in range(1,max_group+1):
        durations = np.array(result.group_durations[group_size]) * time_normalization_factor


        if len(durations) <= 2:
            x = []
            y = []
        else:
            if logarithmic_bins:
                bins = np.append([0.],np.logspace(log10(durations.min())-0.1,log10(durations.max()),n_bins) )
            else:
                bins = np.append([0.],np.linspace(durations.min(), durations.max(), n_bins) )

            x, y = get_ccdf(durations)
            y = tc.sample_a_function(x_contact, y_contact, bins)
            x = bins

        #if group_size == 1:
        #    print('\n',alpha,'\n')
        x_groups.append(x)
        y_groups.append(y)

    xs = [x_k, [], x_contact ] + x_groups
    ys = [y_k, grp_sizes, y_contact ] + y_groups

    return xs, ys

def compute_all_bins(binned_temporal_network,max_group,time_normalization_factor=1./3600.,n_bins=50,logarithmic_bins=False):

    t_fw, k_fw = tc.mean_degree(binned_temporal_network)

    if logarithmic_bins:
        k_fw = k_fw[k_fw>0]
        bins = np.logspace(log10(k_fw.min()),log10(k_fw.max()),n_bins+1)
    else:
        bins = n_bins

    y_k, x_k = np.histogram(k_fw,bins=bins,density=True)
    x_k = get_bin_means(x_k,logarithmic_bins)

    result = tc.measure_group_sizes_and_durations(binned_temporal_network)

    grp_sizes = np.array(result.aggregated_size_histogram[1:])
    max_ndx = np.where(grp_sizes>0)[0][-1]
    grp_sizes = grp_sizes[:max_ndx+1]

    durations = np.array(result.contact_durations) * time_normalization_factor
    if logarithmic_bins:
        bins = np.logspace(log10(durations.min()),log10(durations.max()),n_bins+1)
    else:
        bins = n_bins
    y_contact, x_contact = np.histogram(durations,bins=n_bins,density=True)
    x_contact = get_bin_means(x_contact,logarithmic_bins)

    y_groups = []
    x_groups = []

    for group_size in range(1,max_group+1):
        durations = np.array(result.group_durations[group_size]) * time_normalization_factor

        n = int(min([np.sqrt(len(durations)),n_bins]))

        if len(durations) <= 6:
            x = []
            y = []
        else:
            if logarithmic_bins:
                bins = np.logspace(log10(durations.min()),log10(durations.max()),n)
            else:
                bins = n_bins
            y, x = np.histogram(durations,bins=bins,density=True)
            x = get_bin_means(x,logarithmic_bins)

        #if group_size == 1:
        #    print('\n',alpha,'\n')
        x_groups.append(x)
        y_groups.append(y)

    xs = [x_k, [], x_contact ] + x_groups
    ys = [y_k, grp_sizes, y_contact ] + y_groups

    return xs, ys


if __name__ == "__main__":

    import matplotlib.pyplot as pl

    orig = tc.load_json_taco('~/.tacoma/ht09.taco')
    orig_binned = tc.bin(orig,20.)
    result = tc.measure_group_sizes_and_durations(orig_binned)

    n_bins = 100

    durations = np.array(result.group_durations[1]) / 3600.

    bins = np.append([0.],np.logspace(log10(durations.min())-1,log10(durations.max()),n_bins) )

    x, y = get_ccdf(durations)
    y_sampled = tc.sample_a_function(x,y,bins)

    print("====== HEAD ======")

    print("original", x[:4], y[:4])
    print("sampled", bins[:4], y_sampled[:4])

    print("====== TAIL ======")
    print("original", x[-4:], y[-4:])
    print("sampled", bins[-4:], y_sampled[-4:])

    fig, ax = pl.subplots(1,2)

    ax[0].step(x,y,where='post')
    ax[0].plot(bins, y_sampled)

    ax[0].set_xscale('log')
    ax[0].set_yscale('log')

    P = np.array(result.aggregated_size_histogram)[1:]
    m = np.arange(1,len(P)+1)
    print(len(P), len(m))

    x, y = get_ccdf_from_distribution(m,P)
    print(x,y)

    ax[1].step(x,y,where='post')
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')

    pl.show()
