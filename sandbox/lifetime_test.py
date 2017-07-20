import cFlockwork as cF
import numpy as np
from time import time
from flockworks import get_mean_group_lifetime

import pylab as pl

tau_theory = lambda P: (80 - 80 *P + 10* P**2 - P**3)/((-4 + P)**2* (-2 + P)* (-1 + P))
tau_theory = lambda P:         (4*(-117649 + 268912*P - 211288*P**2 + 58996*P**3 + 3185*P**4 - 3836*P**5 +\
        600*P**6))/((-7 + P)**2*(-1 + P)*(2401 - 4802*P + 3479*P**2 - 1078*P**3 + 120*P**4))

N = 8
Ps = np.linspace(0.01,0.95,20)
tau = np.zeros_like(Ps)
tau_inv = np.zeros_like(Ps)

seed = 345698
for iP,P in enumerate(Ps):
    seed += 1
    a = cF.simulate_P_group_lifetimes(N,P,seed,100000)
    tau[iP] = np.mean(a)
    tau_inv[iP] = get_mean_group_lifetime(N,P)

pl.plot(Ps,tau,label='sim')
pl.plot(Ps,tau_theory(Ps),label='mathematica')
pl.plot(Ps,tau_inv,label='fundamental matrix')
pl.legend()

pl.show()
