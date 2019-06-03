import tacoma as tc
import matplotlib.pyplot as pl
from tacoma.flockwork import flockwork_P_equilibrium_configuration
import _tacoma as _tc

import time


N = 200
t_run_total = 500.0
tmax = 10.0
seed = 7925
R0 = 1.2
recovery_rate = 1.0
infection_rate = R0/(N-1)*recovery_rate

complete_graph = tc.complete_graph(N)

start = time.time()
SIS = tc.SIS(N,t_run_total*2,infection_rate,recovery_rate,N//2,seed = seed)
tc.gillespie_SIS(complete_graph,SIS)
end = time.time()

print("simulation on edge_lists took", end-start,"seconds")

pl.plot(SIS.time, SIS.I)

mv_SIS = tc.MARKOV_SIS(N,t_run_total*2,infection_rate,recovery_rate,0.01,N//2,seed = seed)

start = time.time()
tc.markov_epidemics(complete_graph,mv_SIS,0.01)
#print(mv_SIS.time, mv_SIS.I)
end = time.time()
print("Markov on edge_changes took", end-start,"seconds")

pl.plot(mv_SIS.time, mv_SIS.I)



pl.ylim([0,N])



pl.show()
