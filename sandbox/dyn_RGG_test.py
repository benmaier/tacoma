import tacoma as tc
from rocsNWL import draw
import matplotlib.pyplot as pl
import networkx as nx
from collections import Counter
import numpy as np
import time

N = 400

start = time.time()

edge_lists = tc.dynamic_RGG(N = N,
                            t_run_total = 1000,
                            mean_link_duration = 5.,
                            periodic_boundary_conditions_for_link_building = True,
                            record_sizes_and_durations = True,
                            #verbose = True)
                            seed = 2335
                            )

end = time.time()

print("needed", end-start, "seconds")

G = nx.Graph()
G.add_edges_from(edge_lists.edges[-1])
G.add_nodes_from(list(range(N)))

#draw(G)

#pl.show()

size_counter = Counter()
histograms = edge_lists.size_histograms

duration_counter = Counter(edge_lists.group_durations)

for hist in histograms:

    size_counter += Counter(hist)

counters = [size_counter, duration_counter]

fig, ax = pl.subplots(2,1)

for ia, a in enumerate(ax):

    sizes = np.array(list(counters[ia].items()),dtype=float)
    size = sizes[:,0]
    count = sizes[:,1]

    a.set_xscale('log')
    a.set_yscale('log')

    ax[ia].plot(size,count,'s')

ax[0].set_xlabel('group size')
ax[1].set_xlabel('contact duration [steps]')

print("mean link duration =", np.mean(edge_lists.group_durations))

pl.show()
