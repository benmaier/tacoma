
import tacoma as tc
import numpy as np
import matplotlib.pyplot as pl
import networkx as nx
import time

N = 100
compl = tc.complete_graph(N)


R0 = 2
k = N - 1.0
rho = 1.0
I0 = 4

eta = R0*rho/k

sir = tc.SIR(N,200,eta,rho,I0, save_infection_events = True)
tc.gillespie_SIR(compl,sir,is_static=True)

I = np.array(sir.I)
t = np.array(sir.time)

infection_events = sir.infection_events

print("number saved infection events:", len(infection_events))
print("expected number of saved infection events:", np.max(sir.R)-I0)

T = nx.Graph()
T.add_edges_from(infection_events)
print("is the infection graph a forest?: ", nx.is_forest(T))

pl.figure()
pl.plot(t,I)

pl.show()
