import tacoma as tc
import matplotlib.pyplot as pl
from tacoma.flockwork import flockwork_P_equilibrium_configuration, degree_moment
import _tacoma as _tc

import time


N = 200
P = [ 0.03 ] 
rewiring_rate = [ (0.0,1.0) ]
t_run_total = 500.0
tmax = 2000.0
seed = 7925
R0 = 1.0
infection_rate = 1.0
recovery_rate = 0.03

infection_rate = R0/degree_moment(N,P[0],1)*recovery_rate


E = flockwork_P_equilibrium_configuration(N,P[0])

fwP_ec = tc.flockwork_P_varying_rates(E,N,P,t_run_total,rewiring_rate,tmax,seed=seed)
fwP_el = tc.convert(fwP_ec)

start = time.time()
SIS = tc.SIS(N,t_run_total*2,infection_rate,recovery_rate,N//2,seed = seed)
tc.gillespie_SIS(fwP_ec,SIS)
end = time.time()
print("simulation on edge_changes took", end-start,"seconds")
pl.plot(SIS.time, SIS.I,label='sim on ec')

mv_SIS = tc.MARKOV_SIS(N,t_run_total*2,infection_rate,recovery_rate,0.01,N//2,seed = seed)

start = time.time()
tc.markov_epidemics(fwP_ec,mv_SIS,1.0)
#print(mv_SIS.time, mv_SIS.I)
end = time.time()
print("Markov on edge_changes took", end-start,"seconds")

pl.plot(mv_SIS.time, mv_SIS.I,label='markov on ec')


start = time.time()
SIS = tc.SIS(N,t_run_total*2,infection_rate,recovery_rate,N//2,seed = seed)
tc.gillespie_SIS(fwP_el,SIS)
end = time.time()

print("simulation on edge_lists took", end-start,"seconds")

pl.plot(SIS.time, SIS.I,label='sim on el')

mv_SIS = tc.MARKOV_SIS(N,t_run_total*2,infection_rate,recovery_rate,0.01,N//2,seed = seed)

start = time.time()
tc.markov_epidemics(fwP_el,mv_SIS,1.0)
#print(mv_SIS.time, mv_SIS.I)
end = time.time()
print("Markov on edge_lists took", end-start,"seconds")

pl.plot(mv_SIS.time, mv_SIS.I,label='markov on el',ls='--',lw=0.5)

FW = tc.FlockworkPModel(E, N, rewiring_rate[0][1], P[0])
mv_SIS = tc.MARKOV_SIS(N,t_run_total*2,infection_rate,recovery_rate,0.01,N//2,seed = seed)

start = time.time()
tc.markov_epidemics(FW, mv_SIS,max_dt = 0.01)
#_tc.markov_SIS_on_FlockworkPModel(FW, mv_SIS,max_dt = 0.01,verbose=False)
end = time.time()
print("Markov on Model took", end-start,"seconds")


pl.plot(mv_SIS.time, mv_SIS.I,label='markov on model')

pl.ylim([0,N])

pl.legend()



pl.show()
