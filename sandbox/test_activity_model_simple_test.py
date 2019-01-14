from _tacoma import activity_model
from tacoma import tc
import matplotlib.pyplot as pl

print(" ")
print("=============================== STARTING ===============================")
print("========================================================================")
print(" ")

N = 3
k = 1
rho = k/(N-1.0)

tau = 1
omega = 1.0/tau
t_run_total = 3

temporal_network = activity_model(N,rho,omega,t_run_total,verbose=True,seed=1)

tn = temporal_network
print(tn.edges_initial)
print(tn.edges_in)
print(tn.edges_out)
print(tn.t0)
print(tn.t)

print(tc.verify(tn))

t, k = tc.mean_degree(temporal_network)
#print(t)

tC, C = tc.contact_coverage(temporal_network)
print(tC, C)

fig, ax = pl.subplots(1,2)

ax[0].plot(t,k)
ax[1].plot(tC,C)
pl.show()
