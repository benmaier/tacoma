# -*- coding: utf-8 -*-
"""
This module provides functions to analyze and plot
the results from temporal network analyses, especially
from the :func:`tacoma.api.measure_group_sizes_and_durations`
routine. These functions are not imported by default
as they require `matplotlib` which cannot be easily installed
on some systems like compute clusters due to a missing
X server. If you want to use the routines of this
submodule, please make sure `matplotlib` is installed.
"""

from math import log10

from collections import Counter

import numpy as np

try:
    import matplotlib.pyplot as pl
except ImportError as e:
    print("\033[1m tacoma does not install `matplotlib` as a dependency. Please install it manually. \033[0m")
    raise e


from tacoma import marker_sequence
from tacoma import color_sequence
from tacoma import get_logarithmic_histogram

import _tacoma


def temporal_network_group_analysis(result,
                                    max_group=5,
                                    time_normalization_factor=1.,
                                    time_unit=None,
                                    plot_step=False,
                                    fit_power_law=False,
                                    ax=None,
                                    bins=100,
                                    bin_dt=None,
                                    use_logarithmic_histogram=True,
                                    markersize=4,

                                    ):
    """Analyze the result of 
    :func:`measure_group_sizes_and_durations`
    and plot it.

    Parameters
    ----------
    result : :mod:`tacoma.group_sizes_and_durations`
        The result of the temporal network analysis provided by
        :mod:`tacoma.measure_group_sizes_and_durations`
    max_group : int, optional
        The maximum group size for which to plot the duration distribution.
        default: 5
    time_normalization_factor : float, optional
        Factor with which the time in the duration lists are rescaled
        default: 1.0
    time_unit : :obj:`str`, optional
        Time unit to put on the axes. default : None
    plot_step : bool, optional
        If True, plot step functions instead of markers. default : False
    fit_power_law : bool, optional
        If True, fit and plot a power law to the distributions. default: False
    ax : :obj:`list` of matplotlib axes objects
        The axes to plot to (have to be a list of length 3 at least). If 
        set to None, a figure and axes will be created and returned. 
        default : None
    bins : int, default : 100
        number of bins for the histogram
    bin_dt : float, default : None
        if given, do discrete binning to bins of this time-length
    use_logarithmic_histogram : bool, default : True
        if True, use logarithmicly growing bin sizes, otherwise use
        constant bin size
    markersize : int, default : 4,
        markersize for the plots

    Returns
    -------
    matplotlib figure
        A matplotlib figure object.
    array_like of matplotlib axes
        As the name says.
    :obj:`dict`
        The results of the analysis (the distributions).
    """

    if ax is None:
        fig, ax = pl.subplots(2, 2)
        ax = ax.flatten()
    else:
        fig = None

    res = {}
    res_sizes = plot_group_size_histogram(result,
                                          ax[0],
                                          plot_step=plot_step,
                                          fit_power_law=fit_power_law,
                                          markersize=markersize,
                                          )
    res_contacts = plot_contact_durations(result,
                                          ax[1],
                                          time_normalization_factor=time_normalization_factor,
                                          time_unit=time_unit,
                                          plot_step=plot_step,
                                          fit_power_law=fit_power_law,
                                          use_logarithmic_histogram=use_logarithmic_histogram,
                                          bins=bins,
                                          bin_dt=bin_dt,
                                          markersize=markersize,
                                          )
    res_dur = plot_group_durations(result,
                                   ax[2],
                                   time_normalization_factor=time_normalization_factor,
                                   max_group=max_group,
                                   time_unit=time_unit,
                                   plot_step=plot_step,
                                   fit_power_law=fit_power_law,
                                   use_logarithmic_histogram=use_logarithmic_histogram,
                                   bins=bins,
                                   bin_dt=bin_dt,
                                          markersize=markersize,
                                   )

    res.update(res_sizes)
    res.update(res_contacts)
    res.update(res_dur)

    return fig, ax, res

def detailed_temporal_network_group_analysis(result,
                                    P_k,
                                    max_group=5,
                                    time_normalization_factor=1.,
                                    time_unit=None,
                                    plot_step=False,
                                    fit_power_law=False,
                                    ax=None,
                                    bins=100,
                                    bin_dt=None,
                                    use_logarithmic_histogram=True,
                                    marker=None,
                                    markersize=4,
                                    label=None,
                                    ):
    """Analyze the result of 
    :func:`measure_group_sizes_and_durations`
    and plot it.

    Parameters
    ----------
    result : :mod:`tacoma.group_sizes_and_durations`
        The result of the temporal network analysis provided by
        :mod:`tacoma.measure_group_sizes_and_durations`
    P_k : list or numpy.ndarray of float
        average degree distribution
        entry the :math:`k`-th entry of ``P_k`` contains the average
        probability that a node has degree :math:`k`
    max_group : int, optional
        The maximum group size for which to plot the duration distribution.
        default: 5
    time_normalization_factor : float, optional
        Factor with which the time in the duration lists are rescaled
        default: 1.0
    time_unit : :obj:`str`, optional
        Time unit to put on the axes. default : None
    plot_step : bool, optional
        If True, plot step functions instead of markers. default : False
    fit_power_law : bool, optional
        If True, fit and plot a power law to the distributions. default: False
    ax : :obj:`list` of matplotlib axes objects
        The axes to plot to (have to be a list of length 3 at least). If 
        set to None, a figure and axes will be created and returned. 
        default : None
    bins : int, default : 100
        number of bins for the histogram
    bin_dt : float, default : None
        if given, do discrete binning to bins of this time-length
    use_logarithmic_histogram : bool, default : True
        if True, use logarithmicly growing bin sizes, otherwise use
        constant bin size
    marker : str, default : None,
        if set, all curves will be drawn using this marker
    markersize : int, default : 4,
        markersize for the plots
    label : str, default : None,
        if set, all curves will be associated with this label

    Returns
    -------
    matplotlib figure
        A matplotlib figure object.
    array_like of matplotlib axes
        As the name says.
    :obj:`dict`
        The results of the analysis (the distributions).
    """

    if ax is None:
        fig, ax = pl.subplots(int(np.ceil((max_group-1)/4))+1, 4)
        ax = ax.flatten()
    else:
        fig = None

    res = {}

    res_degree = plot_degree_distribution(P_k,
                                          ax[0],
                                          plot_step=plot_step,
                                          marker=marker,
                                          markersize=markersize,
                                          label=label,
                                          )
    res_sizes = plot_group_size_histogram(result,
                                          ax[1],
                                          plot_step=plot_step,
                                          fit_power_law=fit_power_law,
                                          marker=marker,
                                          markersize=markersize,
                                          label=label,
                                          )
    res_contacts = plot_contact_durations(result,
                                          ax[2:4],
                                          time_normalization_factor=time_normalization_factor,
                                          time_unit=time_unit,
                                          plot_step=plot_step,
                                          fit_power_law=fit_power_law,
                                          use_logarithmic_histogram=use_logarithmic_histogram,
                                          bins=bins,
                                          bin_dt=bin_dt,
                                          marker=marker,
                                          markersize=markersize,
                                          label=label,
                                          )

    for g in range(2, max_group+1):
        res_dur = plot_group_durations(result,
                                       ax[4+g-2],
                                       time_normalization_factor=time_normalization_factor,
                                       min_group=g,
                                       max_group=g,
                                       time_unit=time_unit,
                                       plot_step=plot_step,
                                       fit_power_law=fit_power_law,
                                       use_logarithmic_histogram=use_logarithmic_histogram,
                                       bins=bins,
                                       bin_dt=bin_dt,
                                       markersize=markersize,
                                       marker=marker,
                                       label=label,
                                       )
        res.update(res_dur)

    res.update(res_sizes)
    res.update(res_contacts)
    res.update(res_dur)
    res.update(res_degree)

    return fig, ax, res

def plot_degree_distribution(P_k,
                              ax,
                              xlabel='node degree $k$',
                              plot_step=False,
                              markersize=4,
                              marker=None,
                              label=None,
                              ):

    P = np.array(P_k)
    k = np.arange(np.max(np.where(P>0)[0])+1)
    P = P[k]

    if marker is None:
        marker = 'o'

    if plot_step:
        ax.step(k, P,
                where='mid',
                label=label,
                )
    else:
        ax.plot(k, P, marker,
                ms=markersize,
                mew=1,
                mfc='None',
                label=label,
                )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r'avg. probability $\overline{P_k}$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    res = {'degree_distribution': (k, P)}

    return res

def plot_group_size_histogram(result,
                              ax,
                              xlabel='group size $g$',
                              plot_step=False,
                              fit_power_law=False,
                              markersize=4,
                              marker = None,
                              label=None,
                              ):

    if marker is None:
        marker = 'o'

    group_size_histogram = np.array([
        (size, val)
        for size, val in enumerate(result.aggregated_size_histogram)
        if val > 0.
    ], dtype=float)

    x_group, y_group = group_size_histogram[:, 0], group_size_histogram[:, 1]

    if plot_step:
        ax.step(x_group, y_group,
                where='mid',
                label=label,
                )
    else:
        ax.plot(x_group, y_group, marker,
                ms=markersize,
                mew=1,
                mfc='None',
                label=label,
                )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r'avg. group-count $\overline{n_g}$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    res = {'size_histogram': (x_group, y_group)}

    return res


def plot_contact_durations(result,
                           ax,
                           marker=None,
                           xlabel='duration',
                           bins=100,  # number of bins
                           time_normalization_factor=1.,
                           time_unit=None,
                           bin_dt=None,
                           plot_step=False,
                           fit_power_law=False,
                           use_logarithmic_histogram=True,
                           markersize=4,
                           label=None,
                           ):
    if marker is None:
        markers = ['o', 'd']
    else:
        markers = [marker] * 2

    if label is None:
        labels = ['contact', 'inter-contact']
    else:
        labels = [label] * 2


    if bin_dt is not None:
        use_discrete_dt = True
    else:
        use_discrete_dt = False

    if not hasattr(ax,'__len__'):
        a = [ax, ax]
    else:
        a = ax

    durs = [np.array(result.contact_durations, dtype=float),
            np.array(result.group_durations[1], dtype=float)]
    res = {}

    for i, dur in enumerate(durs):
        if not plot_step:
            if use_logarithmic_histogram:
                x, y = get_logarithmic_histogram(
                    time_normalization_factor*dur, bins)
            elif use_discrete_dt:
                c = Counter(dur / bin_dt)
                total = sum(c.values())
                x = []
                y = []
                for x_, y_ in c.items():
                    x.append(x_* bin_dt)
                    y.append(y_/total / bin_dt)
            else:
                y, x = np.histogram(dur*time_normalization_factor,bins=bins,density=True)
                x = 0.5*(x[1:]+x[:-1])
                print(x.shape,y.shape)

            a[i].plot(x, y, ls='', marker=markers[i], 
                    label=labels[i],
                    ms=markersize,
                    mew=1,
                    mfc='None'
                    )
        else:
            if use_logarithmic_histogram:
                x, y = get_logarithmic_histogram(
                    time_normalization_factor*dur, bins, return_bin_means=False)
            elif use_discrete_dt:
                c = Counter(dur / bin_dt)
                total = sum(c.values())
                x = []
                y = []
                for x_, y_ in c.items():
                    x.append(x_* bin_dt)
                    y.append(y_/total / bin_dt)
                x.append(x[-1]+1)
            else:
                y, x = np.histogram(dur*time_normalization_factor,bins=bins,density=True)
            a[i].step(x, np.append(y, y[-1]),
                    where='post',
                    label=labels[i],
                    )

        res[labels[i]] = (x, y)

    if time_unit is not None:
        xlabel += ' [' + time_unit + ']'
    ylabel = 'probability density'
    if time_unit is not None:
        ylabel += ' [1/' + time_unit + ']'

    for ax in a: 
 
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend()

    return res


def plot_group_durations(result,
                         ax,
                         min_group=2,
                         max_group=5,
                         xlabel='duration',
                         bins=100,
                         time_normalization_factor=1.,
                         bin_dt=None,
                         time_unit=None,
                         plot_step=False,
                         fit_power_law=False,
                         use_logarithmic_histogram=True,
                         markersize=4,
                         marker=None,
                         label = None,
                         ):
    if label is None:
        labels = [ '$g=%d$' % size for size in range(min_group, max_group+1) ]
    else:
        labels = [label] * (max_group+1-min_group)

    if marker is None:
        this_marker_sequence = marker_sequence
    else:
        this_marker_sequence = [marker] * max_group


    if bin_dt is not None:
        use_discrete_dt = True
    else:
        use_discrete_dt = False

    res = {}

    for size in range(min_group, max_group+1):
        if len(result.group_durations[size]) > 6:
            data = time_normalization_factor * \
                np.array(result.group_durations[size], dtype=float)
            if not plot_step:
                if use_logarithmic_histogram:
                    x, y = get_logarithmic_histogram(data, bins)
                elif use_discrete_dt:
                    c = Counter(data / bin_dt)
                    total = sum(c.values())
                    x = []
                    y = []
                    for x_, y_ in c.items():
                        x.append(x_* bin_dt)
                        y.append(y_/total / bin_dt)
                else:
                    y, x = np.histogram(data*time_normalization_factor,bins=bins,density=True)
                    x = 0.5*(x[1:]+x[:-1])



                ax.plot(x, y,
                        ls='',
                        marker=this_marker_sequence[size % len(this_marker_sequence)],
                        label=labels[size-min_group],
                        ms=markersize,
                        mew=1,
                        mfc='None'
                        )
            else:
                if use_logarithmic_histogram:
                    x, y = get_logarithmic_histogram(
                        data, bins, return_bin_means=False)
                elif use_discrete_dt:
                    c = Counter(data / bin_dt)
                    total = sum(c.values())
                    x = []
                    y = []
                    for x_, y_ in c.items():
                        x.append(x_* bin_dt)
                        y.append(y_/total / bin_dt)
                    x.append(x[-1]+1)
                else:
                    y, x = np.histogram(data*time_normalization_factor,bins=bins,density=True)

                ax.step(x, np.append(y, y[-1]),
                        where='post',
                        label=labels[size-min_group],
                        )
            res[size] = (x, y)

    ax.set_xscale('log')
    ax.set_yscale('log')
    if time_unit is not None:
        xlabel += ' [' + time_unit + ']'

    ylabel = 'probability density'
    if time_unit is not None:
        ylabel += ' [1/' + time_unit + ']'

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend()

    return res


def plot_binned_social_trajectory(traj,
                                  ax,
                                  time_unit='d',
                                  ):
    for t, groups in enumerate(traj):
        for group in groups:
            ax.plot([t, t+1], [group, group], '-k')

    ax.set_xlabel('time ['+time_unit+']')
    ax.set_ylabel('group id')


def plot_social_trajectory(traj,
                           ax,
                           time_normalization_factor=1.,
                           time_unit=None
                           ):
    for group, entry in enumerate(traj):
        for t in entry.time_pairs:
            t = np.array(t)
            ax.plot(t*time_normalization_factor, [group, group], '-k')

    xlabel = 'time'
    if time_unit is not None:
        xlabel += '['+time_unit+']'
    ax.set_xlabel(xlabel)
    ax.set_ylabel('group id')


if __name__ == "__main__":

    import tacoma as tc

    # let's simulate a similar dynamic RGG
    RGG_edge_lists = _tacoma.dynamic_RGG(N=412,
                                         t_run_total=int(
                                             24*2*3600 / 300),  # 2 days
                                         mean_link_duration=5.,
                                         periodic_boundary_conditions_for_link_building=False,
                                         record_sizes_and_durations=False,
                                         # verbose = True)
                                         seed=2335
                                         )

    RGG_result = tc.measure_group_sizes_and_durations(RGG_edge_lists)

    fig, ax, data = temporal_network_group_analysis(
        RGG_result, time_normalization_factor=300./3600., time_unit='h')
    fig.tight_layout()

    pl.show()
