import tacoma as tc
import numpy as np
import matplotlib.pyplot as pl

N = 5
compl = tc.complete_graph(N)

R0s = np.logspace(-0.1,1.0,10)

k = N - 1.0
rho = 1.0
I0 = 1

N_meas = 10
#Rs = np.zeros((len(R0s),N_meas))
R0 = 1.5
eta = R0*rho/k
print eta, rho

sir = tc.Dyn_SIR(N,200,eta,rho,I0,verbose=True)

tc.gillespie_SIR_on_edge_lists(compl,sir,verbose=True)

Rs[iR0,meas] = sir.R[-1]

pl.plot(R0s,Rs.mean(axis=1)/N)
pl.plot(R0s,(1-1./R0s))

pl.xscale('log')

pl.show()
