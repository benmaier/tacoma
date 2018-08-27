from _tacoma import activity_model
from tacoma import tc
import matplotlib.pyplot as pl

import numpy as np

N = 100
k = 1
rho = k/(N-1.0)
m = N*(N-1)/2

tau = 4
omega = 1.0/tau
t_run_total = 10000

temporal_network = activity_model(N,rho,omega,t_run_total,seed=1)

print("errors in generated temporal network =",tc.verify(temporal_network))

fig, ax = pl.subplots(1,3,figsize=(10,3))

t1, k_ = tc.mean_degree(temporal_network)

k_mean = tc.time_average(t1, k_)
print(k_mean)
ax[0].plot(t1,k_,label='simulation')
ax[0].plot(t1[[0,-1]],[k_mean]*2,'--',lw=3,label='temporal average')
ax[0].plot(t1[[0,-1]],[k]*2,'-',lw=2,label='demanded')

t, C = tc.number_of_discovered_edges(temporal_network)
m0 = len(temporal_network.edges_initial)
C *= (m-m0) / m
C += m0
ndx = np.where(t>0)
t = t[ndx]
C = C[ndx]
ax[1].plot(t,C,'s',mfc='None',label='simulation')
ax[1].plot(t,m-m*np.exp(-omega*t),label='theory')
#ax[1].set_xscale('log')
#ax[1].set_yscale('log')

P_k = np.array(tc.degree_distribution(temporal_network))
ks = np.where(P_k>0)[0]

from scipy.stats import binom
n = N-1
p = rho
P_theory = binom(n,p).pmf(ks)

ax[2].plot(ks, P_k[ks],'s',mfc='None',label='simulation')
ax[2].plot(ks, P_theory,label='theory')


for i in range(2):
    ax[i].set_xlabel(r'time $t\times\omega$')
ax[2].set_xlabel(r'degree $k$')

ax[0].set_ylabel(r'mean degree')
ax[1].set_ylabel(r'contact coverage')
ax[2].set_ylabel(r'probability $\overline{P(k)}$')

ax[0].legend()
ax[1].legend()
ax[2].legend()
fig.tight_layout()
pl.show()
