import sys
from collections import Counter

import matplotlib.pyplot as pl
import numpy as np

from netwulf import visualize
import networkx as nx

N = 500

from tacoma.flockwork import flockwork_P_equilibrium_configuration 
from tacoma.flockwork import flockwork_P_equilibrium_group_size_distribution

"""
G = nx.Graph()
G.add_nodes_from(range(N))
G.add_edges_from(flockwork_P_equilibrium_configuration(N, P))

visualize(G,config={'Node size exponent': 0, 'Node size': 2})
"""

N_meas = 1000
Ps = [0.01,0.1,0.3,0.5,0.7,0.9,0.99]

for P in Ps:
    this_hist = np.zeros((N+1,))
    for meas in range(N_meas):

        edges, hist = flockwork_P_equilibrium_configuration(N, P, return_histogram=True)

        if len(edges) != len(set(edges)):
            print("P =", P, " duplicate edges")
            e = Counter(edges)
            print(e.most_common(1))
            sys.exit()
            #G = nx.Graph()
            #G.add_nodes_from(range(N))
            #G.add_edges_from(flockwork_P_equilibrium_configuration(N, P))

            #visualize(G,config={'Node size exponent': 0, 'Node size': 2})

        this_hist += hist

    this_hist /= N_meas

    P_th = flockwork_P_equilibrium_group_size_distribution(N, P)


    p, = pl.plot(np.arange(1,N+1), P_th[1:],'-')
    pl.plot(np.arange(1,N+1), this_hist[1:],'--s',mfc='none',c=p.get_color())

#pl.xscale('log')
pl.yscale('log')

pl.ylim(1/N**2,N)

pl.show()
