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
                                          )
    res_contacts = plot_contact_durations(result,
                                          ax[1],
                                          time_normalization_factor=time_normalization_factor,
                                          time_unit=time_unit,
                                          plot_step=plot_step,
                                          fit_power_law=fit_power_law,
                                          )
    res_dur = plot_group_durations(result,
                                   ax[2],
                                   time_normalization_factor=time_normalization_factor,
                                   max_group=max_group,
                                   time_unit=time_unit,
                                   plot_step=plot_step,
                                   fit_power_law=fit_power_law,
                                   )

    res.update(res_sizes)
    res.update(res_contacts)
    res.update(res_dur)

    return fig, ax, res


def plot_group_size_histogram(result,
                              ax,
                              marker='.',
                              xlabel='group size $m$',
                              plot_step=False,
                              fit_power_law=False,
                              ):

    group_size_histogram = np.array([
        (size, val)
        for size, val in enumerate(result.aggregated_size_histogram)
        if val > 0.
    ], dtype=float)

    x_group, y_group = group_size_histogram[:, 0], group_size_histogram[:, 1]

    if plot_step:
        ax.step(x_group, y_group,
                where='mid',
                )
    else:
        ax.plot(x_group, y_group, 'o',
                ms=4,
                mew=1,
                mfc='None'
                )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r'mean number of $m$-sized groups $\overline{N_m}$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    res = {'size_histogram': (x_group, y_group)}

    return res


def plot_contact_durations(result,
                           ax,
                           markers=['o', 'd'],
                           xlabel='duration',
                           bins=100,  # number of bins
                           time_normalization_factor=1.,
                           time_unit=None,
                           plot_step=False,
                           fit_power_law=False,
                           ):

    durs = [np.array(result.contact_durations, dtype=float),
            np.array(result.group_durations[1], dtype=float)]
    labels = ['contact', 'inter-contact']
    res = {}

    for i, dur in enumerate(durs):
        if not plot_step:
            x, y = get_logarithmic_histogram(
                time_normalization_factor*dur, bins)
            ax.plot(x, y, ls='', marker=markers[i], label=labels[i],
                    ms=4,
                    mew=1,
                    mfc='None'
                    )
        else:
            x, y = get_logarithmic_histogram(
                time_normalization_factor*dur, bins, return_bin_means=False)
            ax.step(x, np.append(y, y[-1]),
                    where='post'
                    )

        res[labels[i]] = (x, y)

    ax.set_xscale('log')
    ax.set_yscale('log')
    if time_unit is not None:
        xlabel += ' [' + time_unit + ']'
    ax.set_xlabel(xlabel)

    ylabel = 'probability density'
    if time_unit is not None:
        ylabel += ' [1/' + time_unit + ']'

    ax.set_ylabel(ylabel)
    ax.legend()

    return res


def plot_group_durations(result,
                         ax,
                         max_group=5,
                         xlabel='duration',
                         bins=100,
                         time_normalization_factor=1.,
                         time_unit=None,
                         plot_step=False,
                         fit_power_law=False,
                         ):

    res = {}

    for size in range(2, max_group+1):
        if len(result.group_durations[size]) > 6:
            data = time_normalization_factor * \
                np.array(result.group_durations[size], dtype=float)
            if not plot_step:
                x, y = get_logarithmic_histogram(data, bins)
                ax.plot(x, y,
                        ls='',
                        marker=marker_sequence[size % len(marker_sequence)],
                        label='$m=%d$' % size,
                        ms=4,
                        mew=1,
                        mfc='None'
                        )
            else:
                x, y = get_logarithmic_histogram(
                    data, bins, return_bin_means=False)
                ax.step(x, np.append(y, y[-1]),
                        where='post',
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
            ax.plot(t*time_normalization_factor, [group, group], '-')

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
