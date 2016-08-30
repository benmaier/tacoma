import EqFlockwork
import EpiFlockwork
from flockworks import flockwork
import networkx as nx
import pylab as pl
import seaborn as sns
from numpy import *

N = 100
Q = 0.1

#=========== create and equilibrate =============
G = nx.complete_graph(N)
#edge_list = EqFlockwork.equilibrate(G.edges(),N,Q)
edge_list = G.edges()

init_seed = 3833
N_measurements = 50

kw = {}
kw['E'] = edge_list
kw['N'] = N
kw['Q'] = Q
kw['t_run_total'] = 100
kw['rewiring_rate'] = 0.0
kw['recovery_rate'] = 1.0

kw["number_of_vaccinated"] = 0
kw["number_of_infected"] = 10

R0s = linspace(0.5,5.,10)
print R0s
epi_curve = zeros_like(R0s)
epi_all = zeros((len(R0s),N_measurements))

seed = init_seed
for iR0, R0 in enumerate(R0s):
    current_i_infs = []
    for meas in range(N_measurements):
        print R0, kw['recovery_rate'], N-1
        kw['infection_rate'] = R0 * kw['recovery_rate'] / (N-1.)
        kw["seed"] = seed

        result = EpiFlockwork.SIS(**kw)
        I = result.I_of_t

        print R0, I[-1]
        I = array(I)
        seed += 1
        current_i_infs.append(I[-1,1]/float(N))
    epi_curve[iR0] = mean(current_i_infs)
    epi_all[iR0,:] = array(current_i_infs)
    print iR0



R0s_violin = array([ [R0]*N_measurements for R0 in R0s] ).flatten()
print array([ [R0]*N_measurements for R0 in R0s] ).flatten(), epi_all.flatten()
ax = sns.violinplot(R0s_violin, epi_all.flatten())
#pl.plot(R0s,epi_curve)
ax.plot(R0s[R0s>=1],1-1./R0s[R0s>=1])
print R0s[R0s>=1],1-1./R0s[R0s>=1]
pl.show()
    



