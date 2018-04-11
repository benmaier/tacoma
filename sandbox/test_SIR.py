
import tacoma as tc
import numpy as np
import matplotlib.pyplot as pl
import GillEpi
import networkx as nx
import time

N = 100
compl = tc.complete_graph(N)
G = nx.complete_graph(N)

R0s = np.logspace(-0.5,1.0,10)

k = N - 1.0
rho = 1.0
I0 = 1

N_meas = 100
Rs = np.zeros((len(R0s),N_meas))
R_gillepi_s = np.zeros((len(R0s),N_meas))
R_s = np.zeros((len(R0s),N_meas))

fig, ax = pl.subplots(4,3,figsize=(18,14))
ax = ax.flatten()


for iR0,R0 in enumerate(R0s):

    for meas in range(N_meas):
        eta = R0*rho/k

        start1 = time.time()
        sir = tc.SIR(N,200,eta,rho,I0)

        tc.gillespie_SIR(compl,sir,is_static=True)
        end1 = time.time()

        Rs[iR0,meas] = sir.R[-1]

        
        start2 = time.time()
        sir2 = GillEpi.SIR(
                          G,
                          infection_rate = eta,
                          recovery_rate = rho,
                          infection_seeds = 1,
                         )

        sir2.simulate()
        end2 = time.time()


        r, tr = sir2.get_r_of_t()
        i, ti = sir2.get_i_of_t()
        R_gillepi_s[iR0,meas] = r[-1]*N

        ax[iR0].plot(tr,r*N,'-',c='b',alpha=0.1,lw=1)
        ax[iR0].plot(ti,i*N,'-',c='c',alpha=0.1,lw=1)
        ax[iR0].plot(sir.time,sir.I,'-',c='r',alpha=0.1,lw=1)
        ax[iR0].plot(sir.time,sir.R,'-',c='g',alpha=0.1,lw=1)

        print(R0, meas, end1-start1, end2-start2)

pl.figure()
pl.plot(R0s,Rs.mean(axis=1)/N)
pl.plot(R0s,R_gillepi_s.mean(axis=1)/N)
pl.plot(R0s[R0s>1],(1-1./R0s[R0s>1]))

pl.xscale('log')


pl.show()
