import tacoma as tc
import matplotlib.pyplot as pl
from tacoma.flockwork import flockwork_P_equilibrium_configuration

import time


N = 100
P = [ 0.5 ] 
rewiring_rate = [ (0.0,1.0) ]
t_run_total = 500.0
tmax = 1000.0
seed = 7925
infection_rate = 1.0
recovery_rate = 0.1

E, _ = flockwork_P_equilibrium_configuration(N,P[0])

fwP_ec = tc.flockwork_P_varying_rates(E,N,P,t_run_total,rewiring_rate,tmax,seed=seed)
fwP_el = tc.convert(fwP_ec)

start = time.time()
SIS = tc.node_based_SIS(N,t_run_total*2,infection_rate,recovery_rate,N//2,seed = seed,verbose=True)
tc.gillespie_node_based_SIS(fwP_ec,SIS,verbose=True)
end = time.time()

print("node-based simulation on edge_changes took", end-start,"seconds")

pl.plot(SIS.time, SIS.I)

start = time.time()
SIS = tc.node_based_SIS(N,t_run_total*2,infection_rate,recovery_rate,N//2,seed = seed)
tc.gillespie_node_based_SIS(fwP_el,SIS)
end = time.time()

print("node-based simulation on edge_lists took", end-start,"seconds")

pl.plot(SIS.time, SIS.I)


start = time.time()
SIS = tc.SIS(N,t_run_total*2,infection_rate,recovery_rate,N//2,seed = seed,verbose=True)
tc.gillespie_SIS(fwP_ec,SIS,verbose=True)
end = time.time()

print("edge-based simulation on edge_changes took", end-start,"seconds")

pl.plot(SIS.time, SIS.I)

start = time.time()
SIS = tc.SIS(N,t_run_total*2,infection_rate,recovery_rate,N//2,seed = seed)
tc.gillespie_SIS(fwP_el,SIS)
end = time.time()

print("edge-based simulation on edge_lists took", end-start,"seconds")

pl.plot(SIS.time, SIS.I)





pl.show()
