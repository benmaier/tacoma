import numpy as np

import matplotlib.pyplot as pl

import tacoma as tc
from tacoma.model_conversions import estimate_flockwork_P_args

sample_aggregates = False
N_time_steps = 2

print("===== edge_lists => edge_lists =====")

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
              (1,2)
            ],
           ]

params = estimate_flockwork_P_args(L,dt=1.1,adjust_last_bin_if_dt_does_not_fit=True)

g = np.array(params['rewiring_rate'])
t, g = g[:,0], g[:,1]

t = np.append(t, params['tmax'])
g = np.append(g, g[-1])
print(t, g)

P = np.array(params['P'])
P = np.append(P, P[-1])
print(t, P)

fig, ax = pl.subplots(1,2)

ax[0].step(t,g,where='post')
ax[1].step(t,P,where='post')

# second estimation

params = estimate_flockwork_P_args(L, dt = 1.0)

g = np.array(params['rewiring_rate'])
t, g = g[:,0], g[:,1]

t = np.append(t, params['tmax'])
g = np.append(g, g[-1])
print(t, g)

P = np.array(params['P'])
P = np.append(P, P[-1])
print(t, P)
ax[0].step(t,g,where='post')
ax[1].step(t,P,where='post')

pl.show()

