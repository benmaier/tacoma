import tacoma as tc
import pdb
from tacoma.epidemics import simulate_quasi_stationary_SIS_on_static_network 

N = 100
R0 = 1.1
rho = 1.0
eta = rho * R0 / (N-1.0)

tmax = 100000

M = 100

G = tc.complete_graph(N,tmax=tmax)


QS = tc.QS_SIS(N,
               tmax,
               eta,
               rho,
               M,
               sampling_rate=2*rho,
               number_of_initially_infected=N,
               sample_network_state=False,
               )


I, I2 = simulate_quasi_stationary_SIS_on_static_network(G, QS,verbose=False)


#pdb.set_trace()

print(I, I2)


