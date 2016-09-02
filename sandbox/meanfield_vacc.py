import EqFlockwork
import cFlockwork
from flockworks import flockwork
import networkx as nx
import pylab as pl
import seaborn as sns
from numpy import *
import os

N = 200
Q = 0.1

#=========== create and equilibrate =============
G = nx.complete_graph(N)
#edge_list = EqFlockwork.equilibrate(G.edges(),N,Q)
edge_list = G.edges()

init_seed = 383
N_measurements = 50

only_nontrivial_eq = True

kw = {}
kw['E'] = edge_list
kw['N'] = N
kw['Q'] = Q
kw['t_run_total'] = 110
kw['rewiring_rate'] = 0.0
kw['recovery_rate'] = 1.0

kw["number_of_vaccinated"] = 0
kw["number_of_infected"] = 10

R0s = linspace(0.5,5.,10)
vrs = linspace(0,0.9,5)
print R0s
if os.path.exists("./epi_curve.npy"):
    epi_curve = load("./epi_curve.npy")
else:
    epi_curve = zeros((len(vrs),len(R0s)))


    seed = init_seed
    for ivr, vr in enumerate(vrs):
        kw["number_of_vaccinated"] = int(N*vr)
        for iR0, R0 in enumerate(R0s):
            current_i_infs = []
            for meas in range(N_measurements):
                print R0, kw['recovery_rate'], N-1
                kw['infection_rate'] = R0 * kw['recovery_rate'] / (N-1.)
                kw["seed"] = seed
                I, SI, R0_, new_edge_list = cFlockwork.SIS(**kw)
                print R0, I[-1]
                I = array(I)
                seed += 1
                if only_nontrivial_eq:
                    if I[-1,1] > 0:
                        current_i_infs.append(I[-1,1]/float(N))
                else:
                    current_i_infs.append(I[-1,1]/float(N))

            curr_res = current_i_infs
            
            if len(curr_res)==0:
                val = 0.
                err = 0.
            elif len(curr_res)==1:
                val = curr_res[0]
                err = curr_res[0]
            else:
                val = mean(curr_res)
                err = std(curr_res) / sqrt(len(curr_res)-1)


            epi_curve[ivr,iR0] = val


    save("epi_curve.npy",epi_curve)

fig,ax = pl.subplots(1,2)

for ivr,vr in enumerate(vrs):
    p, = ax[0].plot(R0s,epi_curve[ivr],'o')
    ax[0].plot(R0s[R0s>=1],(1-1./R0s[R0s>=1])-vr,c=p.get_color())
    print R0s[R0s>=1],1-1./R0s[R0s>=1]

for iR0, R0 in enumerate(R0s):
    p, = ax[1].plot(vrs,epi_curve[:,iR0],'o')
    ax[1].plot(vrs,epi_curve[0,iR0]-vrs,c=p.get_color())
    

pl.show()
    



