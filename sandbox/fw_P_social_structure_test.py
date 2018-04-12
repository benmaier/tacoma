import tacoma as tc
import networkx as nx
import numpy as np
import pylab as pl
import seaborn as sns

N = 5

sfig, sax = pl.subplots(1,2,figsize=(10,5))

neighbor_affinity = [ 
                        ([ 1, 2 ], [ 0.9, 0.1]),
                        ([ 0, 2 ], [ 0.9, 0.9]),
                        ([ 0, 1, 3 ], [ 0.1, 0.9, 0.5 ]),
                        ([ 2, 4, ], [ 0.5, 0.01 ]),
                        ([ 3 ], [ 0.01 ])
                    ]

A = np.zeros((N,N))
for u in range(N):
    nodes = neighbor_affinity[u][0]
    vals = neighbor_affinity[u][1]
    for v,val in zip(*neighbor_affinity[u]):
        A[u,v] = val


max_A = np.max(np.max(A))
A = A/max_A

sax[0].imshow(A)



edge_changes = tc.flockwork_P_varying_rates_neighbor_affinity([],N,[0.4,0.],600000,[ (0, 1./3600.), (300,2./3600.) ],neighbor_affinity, 600,seed=2545,use_preferential_node_selection=True)



fw_t = np.array([0.]+edge_changes.t)
edges_in = list(edge_changes.edges_in)
edges_out = list(edge_changes.edges_out)
fw_dt = fw_t[1:] - fw_t[:-1]

G = nx.Graph()
G.add_nodes_from(list(range(N)))
A_fw = np.zeros((N,N))

fw_k = []

for it in range(1,len(fw_t)):
    print(it/float(len(fw_t)))
    dt = fw_t[it]-fw_t[it-1]
    for u in G.nodes():
        for v in G.neighbors(u):
            A_fw[u,v] += dt

    fw_k.append(np.mean([v for k,v in G.degree()]))
    G.remove_edges_from(edges_out[it-1])
    G.add_edges_from(edges_in[it-1])

#fw_k.append(np.mean(G.degree().values()))
#ax[5].step(fw_t/24./3600.,fw_k,where="post",alpha=0.3)

#for u in xrange(f2fw.f2f.N):
#    ax[3].plot(A[u,:],A_fw[u,:],'.',alpha=.1)

max_A_fw = np.max(np.max(A_fw))
A_fw = A_fw/max_A_fw

sax[1].imshow(A_fw)

print(A)
print(A_fw)

pl.show()
