import numpy as np
from tacoma.epidemics import SIR_weighted


import matplotlib.pyplot as pl
import networkx as nx
import tacoma as tc

colors = ['k','r']

N = 20
R0 = 10
recovery_rate = 1.0
I0 = 10
tmax = 1e9

#G = nx.fast_gnp_random_graph(N, k/(N-1))
#G = nx.random_regular_graph(k, N)
G = nx.complete_graph(N)
k = N-1.0
infection_rate = R0 * recovery_rate / k

edge_weight_tuples = [ ( u, v, 1.0) for u, v in G.edges() ]

t_sample = np.linspace(0,10,100)
i_sample = np.zeros_like(t_sample)
i_sample_tc = np.zeros_like(t_sample)
N_measurements = 1000
for meas in range(N_measurements): 

    sir = SIR_weighted(N, edge_weight_tuples, infection_rate, recovery_rate, I0)

    t, I, R = sir.simulation(tmax)

    print("simulation time:", t[-1])

    i_sample += tc.sample_a_function(t,I/N, t_sample) / N_measurements
    #pl.plot(t[[0,-1]], [1-1/R0]*2)
    #pl.plot(t[[0,-1]], [I_mean/N]*2, '--', lw=2,)
    #==============================

    _G = tc.convert_static_network(N, list(G.edges()))
    sir = tc.SIR(N, tmax, infection_rate, recovery_rate, I0)

    tc.gillespie_SIR(_G, sir)
    t = np.array(sir.time)
    I = np.array(sir.I)

    #pl.figure()
    print("simulation time tacoma:", t[-1])

    i_sample_tc += tc.sample_a_function(t,I/N, t_sample) / N_measurements
    #pl.plot(t[[0,-1]], [1-1/R0]*2)
    #pl.plot(t[[0,-1]], [I_mean/N]*2, '--', lw=2,)
    #==============================

pl.step(t_sample, i_sample, where='post', color=colors[1],label='Python method')
pl.step(t_sample, i_sample_tc, where='post', color=colors[0],label='C++ method')
pl.ylim([0,1.0])
pl.legend()
pl.show()


        


            


