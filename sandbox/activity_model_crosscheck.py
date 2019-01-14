from _tacoma import activity_model
from tacoma import tc
import matplotlib.pyplot as pl

from _tacoma import EdgeActivityModel
from _tacoma import simulate_EdgeActivityModel

import numpy as np

from time import time

N = 150
k = 5
rho = k/(N-1.0)
m = N*(N-1)/2

tau = 4
omega = 1.0/tau
t_run_total = 1000

for testid in range(3):

    if testid == 0:
        start = time()
        temporal_network = activity_model(N,rho,omega,t_run_total,seed=1)
        end = time()
        figtitle = "Traditional simulation"
        print("Traditional simulation took",end-start,"seconds.")

    elif testid == 1:
        start = time()
        AM = EdgeActivityModel(N, rho, omega,use_rejection_sampling_of_non_neighbor=True)
        simulate_EdgeActivityModel(AM, t_run_total)
        temporal_network = AM.edge_changes
        end = time()
        figtitle = "Rejection sampling on model class + template"
        print("Rejection sampling on model class + template",end-start,"seconds.")

    elif testid == 2:
        start = time()
        AM = EdgeActivityModel(N, rho, omega, use_rejection_sampling_of_non_neighbor=False,verbose=False)
        temporal_network = simulate_EdgeActivityModel(AM, t_run_total,verbose=False)
        temporal_network = AM.edge_changes
        end = time()
        figtitle = "Traditional sampling on model class + template"
        print("Gillespie simulation with traditional sampling on model class + template",end-start,"seconds.")


    print("errors in generated temporal network =",tc.verify(temporal_network))

    fig, ax = pl.subplots(1,3,figsize=(10,3))

    t1, k_ = tc.mean_degree(temporal_network)

    k_mean = tc.time_average(t1, k_)
    print(k_mean)
    ax[0].plot(t1,k_,label='simulation')
    ax[0].plot(t1[[0,-1]],[k_mean]*2,'--',lw=3,label='temporal average')
    ax[0].plot(t1[[0,-1]],[k]*2,'-',lw=2,label='demanded')

    t, C = tc.contact_coverage(temporal_network)
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
    fig.suptitle(figtitle)
    fig.tight_layout()
pl.show()
