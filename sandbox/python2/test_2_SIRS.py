import tacoma as tc
import numpy as np
import matplotlib.pyplot as pl

N = 50
compl = tc.complete_graph(N)

R0s = np.logspace(-0.1,1.0,10)

k = N - 1.0
rho = 1.0
I0 = 25

N_meas = 10
#Rs = np.zeros((len(R0s),N_meas))
R0 = 4.0
eta = R0*rho/k
print eta, rho

sirs = tc.SIRS(N,200,eta,rho,1.0,number_of_initially_infected=I0)

tc.gillespie_SIRS(compl,sirs)

pl.plot(sirs.time,sirs.R)
pl.plot(sirs.time,sirs.I)

#Rs[iR0,meas] = sirs.R[-1]

#pl.plot(R0s,Rs.mean(axis=1)/N)
#pl.plot(R0s,(1-1./R0s))

#pl.xscale('log')

pl.show()
