import EqFlockwork
import EpiFlockwork
from flockworks import flockwork
import networkx as nx
import pylab as pl
import seaborn as sns
from numpy import *


N = 1000
Q = 0.9
k0 = min([N-1,1./(1.-Q)])


#=========== create and equilibrate =============
G = nx.fast_gnp_random_graph(N,p=k0/N)
edge_list = EqFlockwork.equilibrate(G.edges(),N,Q)
#edge_list = EqFlockwork.equilibrate([],N,Q)

G_eq = nx.Graph()
G_eq.add_edges_from(edge_list)
G_eq.add_nodes_from(xrange(N))

F = flockwork(Q,G_eq)

print F.mean_degree(), 1./(1.-Q)
print "number of groups: ", 
mean_C, mean_L = F.mean_clustering_coeff(with_L=True)
print "mean clustering coeff:",  mean_C, Q/1.24
print "mean number of links between neighbors:",  mean_L, 0.5*Q/(1-Q)**2
N_mod, sizes = F.components_dist()
print "number of components:", N_mod, (1.-Q)/(2.-Q)*N
print "mean size of components:", mean(sizes), (2.-Q)/(1.-Q)


#pos = nx.nx_pydot.graphviz_layout(F.G)
#nx.draw(F.G,pos)

fig,ax = pl.subplots(1,1)

k,hist = F.deg_dist()

ax.plot(k,hist)
ax.set_xscale("log")
ax.set_yscale("log")




pl.show()


