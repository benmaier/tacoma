import tacoma as tc
import matplotlib.pyplot as pl
from tacoma.flockwork import flockwork_P_equilibrium_configuration

import time


N = 1000
P = [ 0.5 ] 
rewiring_rate = [ (0.0,1.0) ]
t_run_total = 10
tmax = 10
seed = 7925
infection_rate = 1.0
recovery_rate = 0.001

E = flockwork_P_equilibrium_configuration(N,P[0])

fwP_ec = tc.flockwork_P_varying_rates(E,N,P,t_run_total,rewiring_rate,tmax,seed=seed)

SIS = tc.SIS(N,t_run_total*2,infection_rate,recovery_rate,seed = seed)
SIS.set_node_status([1]+[0 for i in range(1,N)])
tc.gillespie_SIS(fwP_ec,SIS)
t0, I0 = list(SIS.time), list(SIS.I)
node_status = SIS.get_node_status()


pl.plot(t0, I0)

start = time.time()
SIS = tc.SIS(N,t_run_total*2,infection_rate,recovery_rate,seed = seed)
SIS.set_node_status(node_status)
tc.gillespie_SIS(fwP_ec,SIS)
t1, I1 = list(SIS.time), list(SIS.I)

t1 = [ _ + t0[-1] for _ in t1]

pl.plot(t1, I1)




pl.show()
