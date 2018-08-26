from _tacoma import activity_model
from tacoma import tc
import matplotlib.pyplot as pl

import numpy as np

N = 1000
k = 3
rho = k/(N-1.0)

tau = 4
omega = 1.0/tau
t_run_total = 10

temporal_network = activity_model(N,rho,omega,t_run_total,seed=1)

print("errors in generated temporal network =",tc.verify(temporal_network))

fig, ax = pl.subplots(1,2)

t1, k_ = tc.mean_degree(temporal_network)
ax[0].plot(t1,k_)
ax[0].plot(t1[[0,-1]],[tc.time_average(t1, k_)]*2,'--',lw=3)
ax[0].plot(t1[[0,-1]],[k]*2,'-',lw=2)

t, C = tc.number_of_discovered_edges(temporal_network)
ax[1].plot(t,C,'s',mfc='None')
m = N*(N-1)/2
ax[1].plot(t,m-m*np.exp(-omega*t))

pl.show()
