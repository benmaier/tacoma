import matplotlib.pyplot as pl
from matplotlib.collections import LineCollection
import networkx as nx
import numpy as np

import tacoma as tc

from rocsNWL.drawing import draw
from rocsNWL.drawing import get_pos

_layout_function = 'graphviz'

def draw_edge_lists(L):

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


def draw_rows(L):

    t = L.t + [L.tmax]

    edge_numbers = {}

    fig, ax = pl.subplots(1,1)

    def get_edge_number(edge):
        if edge not in edge_numbers:
            edge_numbers[edge] = len(edge_numbers)
        return edge_numbers[edge]

    for i in range(len(t)-1):
        t_ = t[i:i+2]
        for edge in L.edges[i]:
            num = get_edge_number(edge)
            ax.plot(t_,[num,num],'-k')

def draw_edges(traj,time_normalization_factor=1.,time_unit=None):

    fig, ax = pl.subplots(1,1)
    lines = []
    max_i = len(traj)
    all_t_max = []
    all_t_min = []
    for i, entry in enumerate(traj):
        all_t_max.append(entry.time_pairs[-1][-1] * time_normalization_factor)
        all_t_min.append(entry.time_pairs[0][0] * time_normalization_factor)
        for t_ in entry.time_pairs:
            t_ = np.array(t_) * time_normalization_factor
            lines.append([ t_, [i,i]])

    fit_x = np.array(all_t_min)
    fit_y = np.arange(len(traj),dtype=float)

    

    lines = [zip(x, y) for x, y in lines]

    ax.add_collection(LineCollection(lines))

    ax.set_ylim(-1,max_i)
    ax.set_xlim(min(all_t_min),max(all_t_max))

    xlabel = 'time'
    if time_unit is not None:
        xlabel += ' ['+time_unit+']'

    ylabel = 'edge id'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


    return fig, ax






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

    L = tc.dynamic_RGG(100,100,mean_link_duration = 5)
    F = tc.flockwork_P_varying_rates([],100,[0.5],100,[(0.0,1.0)],tmax=100)
    FBIN = tc.bin_temporal_network(F,dt=1)
    #draw_rows(FBIN)

    start = time.time()
    traj = tc.edge_trajectories(FBIN)
    end = time.time()

    print "needed ", end-start, "seconds"
    draw_edges(traj)

    start = time.time()
    traj = tc.edge_trajectories(F)
    end = time.time()
    print "needed ", end-start, "seconds"
    draw_edges(traj)



    #draw_edge_lists(L)

    pl.show()
