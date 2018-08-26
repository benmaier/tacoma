from __future__ import print_function
import matplotlib.pyplot as pl
from matplotlib.collections import LineCollection
import networkx as nx
import numpy as np

from scipy.optimize import curve_fit

import tacoma as tc

import community

_layout_function = 'graphviz'

def draw_edge_lists(L):
    """this draws a force-directed layout for each snapshot of a temporal network and hence should be used with caution"""
    from rocsNWL.drawing import draw
    from rocsNWL.drawing import get_pos

    G = nx.Graph()
    G.add_nodes_from(range(L.N))
    G.add_edges_from(L.edges[0])
    
    pos, _ = get_pos(G,layout_function=_layout_function)

    fig,ax = pl.subplots(1,len(L.edges),figsize=(len(L.edges)*3,4))
    ax = ax.flatten()

    draw(G,pos=pos,ax=ax[0],layout_function=_layout_function)

    for i in range(1,L.N):
        G = nx.Graph()
        G.add_nodes_from(range(L.N))
        G.add_edges_from(L.edges[i])
        
        draw(G,pos=pos,ax=ax[i],layout_function=_layout_function)

def fit_function(x, alpha, scale, fac, intervals_to_discard_for_fit ):

    x_ = x.copy()

    offset = 0.0

    for interval in intervals_to_discard_for_fit:
        t0, t1 = interval
        x_[np.logical_and(x>=t0, x<t1)] = t0 - offset

        x_[x>=t1] -= t1 - t0
        offset += t1 - t0

    #print(x_)

    #pl.figure()
    #pl.plot(x, x_)
    #pl.show()

    #return fac * (1.0 - np.exp(-alpha*(x_-x_.min())))
    return fac * (1.0 - (scale/(scale+x_))**(alpha) )

    

def draw_edges(traj,
               time_normalization_factor = 1.,
               time_unit = None,
               ax = None,
               fit = False,
               edge_order = None,
               color = None,
               intervals_to_discard_for_fit = [],
               fit_color = 'k',
               return_fit_params = False,
               ):

    if ax is None:
        fig, ax = pl.subplots(1,1)
    else:
        fig = None

    if color is None:
        color = 'k'


    lines = []
    max_i = len(traj)
    all_t_max = []
    all_t_min = []
    max_node = []
    for i, entry in enumerate(traj):
        all_t_max.append(entry.time_pairs[-1][-1] * time_normalization_factor)
        all_t_min.append(entry.time_pairs[0][0] * time_normalization_factor)
        max_node.extend(entry.edge)
        for t_ in entry.time_pairs:
            t_ = np.array(t_) * time_normalization_factor

            if edge_order is not None:
                y = edge_order[i]
            else:
                y = i

            lines.append([ t_, [y,y]])

    #if intervals_to_discard_for_fit is not None:
    #    fit_x = []
    #
    #    for t_ in all_t_min:
    #
    #else:
    fit_x = np.array(all_t_min)

    fit_y = np.arange(len(traj),dtype=float)
    

    lines = [list(zip(x, y)) for x, y in lines]
    colors = [ color for _ in range(len(lines)) ] 

    ax.add_collection(LineCollection(lines,colors=colors,alpha=0.5,linewidth=1))

    t0 = min(all_t_min)
    ax.set_ylim(-1,max_i)
    ax.set_xlim(t0,max(all_t_max))

    xlabel = 'time'
    if time_unit is not None:
        xlabel += ' ['+time_unit+']'

    ylabel = 'edge id'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)



    if fit:
        N = max(max_node) + 1
        fac = N*(N-1)/2.
        fit = lambda x, alpha, scale: fit_function(x, alpha, scale, fac, intervals_to_discard_for_fit)
        #popt, pcov = curve_fit(fit, fit_x, fit_y,[1./fac,fac,10.0],maxfev=10000)
        #popt, pcov = curve_fit(fit, fit_x, fit_y,[2,fac,10.0],maxfev=10000)
        popt, pcov = curve_fit(fit, fit_x, fit_y,[0.5,10.0],maxfev=10000)

        #print (popt)

        ax.plot(fit_x, fit(fit_x,*popt),'r')

        log_y = np.log(fit_y) - 1.
        log_x = np.log(fit_x) - 1.

    if not return_fit_params:
        return fig, ax
    else:
        return fig, ax, popt

def edge_activity_plot(temporal_network,
                       time_normalization_factor = 1.,
                       time_unit = None,
                       ax = None,
                       fit = False,
                       edge_order = None,
                       color = None,
                       intervals_to_discard_for_fit = [],
                       fit_color = None,
                       return_fit_params = False,
                       ):
    result = tc.get_edge_trajectories(temporal_network)
    traj = result.trajectories
    return draw_edges(
                       traj,
                       time_normalization_factor = time_normalization_factor,
                       time_unit = time_unit,
                       ax = ax,
                       fit = fit,
                       edge_order = edge_order,
                       color = color,
                       intervals_to_discard_for_fit = intervals_to_discard_for_fit,
                       fit_color = fit_color,
                       return_fit_params = return_fit_params,
                     )

def get_edge_order(edge_traj,threshold=0.):

    G = get_edge_graph(edge_traj,threshold = 0.)

    partition = community.best_partition(G)
    N_comm = max([v for v in partition.values()]) + 1

    comm = [ [] for i in range(N_comm) ]

    for k, v in partition.items():
        comm[v].append(k)

    order = []
    for module in comm:
        order.extend(module)

    order = np.argsort(order)

    return order


def get_edge_graph(edge_traj,threshold = 0.):

    N_edges = len(edge_traj.trajectories)
    G = nx.Graph()
    G.add_nodes_from(range(N_edges))
    G.add_edges_from([ (u,v) for u,v,val in edge_traj.edge_similarities if val > threshold])

    return G


if __name__ == "__main__":
    import time

    L = tc.edge_lists()

    L.N = 3
    L.t = [0.0,1.0,2.0]
    L.tmax = 3.0
    L.edges = [ 
                [
                  (0,1)
                ],
                [
                  (1,2), (0,2)
                ],
                [
                  (0,1)
                ],
               ]

    L = tc.dynamic_RGG(100,100,mean_link_duration = 10)
    #F = tc.flockwork_P_varying_rates([],100,[0.5],100,[(0.0,1.0)],tmax=100)
    F = L
    FBIN = tc.bin(F,dt=1)
    #draw_rows(FBIN)

    start = time.time()
    result = tc.get_edge_trajectories(FBIN)
    end = time.time()
    traj = result.trajectories
    print(result.edge_similarities)

    print("needed ", end-start, "seconds")
    draw_edges(traj,fit=True)

    start = time.time()
    result = tc.get_edge_trajectories(F)
    end = time.time()
    traj = result.trajectories
    print("needed ", end-start, "seconds")
    draw_edges(traj)



    #draw_edge_lists(L)

    pl.show()
