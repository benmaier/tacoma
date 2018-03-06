import cFlockwork as cF
import matplotlib.pyplot as pl
import DynGillEpi as gill
import numpy as np


N = 50
t_run_total = 10
complete_graph = []

for i in range(N-1):
    for j in range(i+1,N):
        complete_graph.append((i,j))

new_list = [ complete_graph for g in range(t_run_total) ]

N_meas = 10
L_ = cF.dynamic_RGG(50,t_run_total=t_run_total,mean_link_duration=2)
L = cF.edge_lists()
L.copy_from(L_)
L.t = [0., 0.1,4.0, 4.1,5.2,6.9,7.5,8.1,8.5,9.0]
#print len(L.edges)

L.edges = new_list

etas = np.logspace(-0.1,0.5,10)
R0 = np.logspace(-0.1,1,10)
rho = 1.0
etas = R0*rho/(N-1.)

#R0 = etas * (N-1.)/rho
cF_vals = np.zeros((len(etas),N_meas))
dyn_vals = np.zeros((len(etas),N_meas))

seed = 12

import time

for ieta,eta in enumerate(etas):
    for meas in range(N_meas):
        start = time.time()
        L_sis = cF.Dyn_SIS(L.N, 100, eta, rho,number_of_initially_infected = 25,seed=seed)
        cF.gillespie_SIS_on_edge_lists(L,L_sis)
        end = time.time()
        #print "cFlockwork needed", end-start, "seconds"

        t = np.array(L_sis.time)
        dt = t[1:] - t[:-1]
        I_mean = np.array(dt).dot(L_sis.I[:-1])/sum(dt)

        if L_sis.I[-1] == 0:
            I_mean = 0

        start = time.time()
        result = gill.SIS_Poisson_homogeneous(L.N,L.edges,eta,rho,100,seed=seed,initial_number_of_infected=25,verbose=False)
        end = time.time()
        #print "dynGillEpi needed", end-start, "seconds"

        I_2_mean = np.mean(np.array(result.I[0][:-1],dtype=float))
        I_2_mean = result.I[0][-1]
        if result.I[0][-1] == 0:
            I_2_mean = 0

        t = np.array(result.true_t)
        dt = t[1:] - t[:-1]
        I_2_mean = np.array(dt).dot(result.true_I[:-1])/sum(dt)
        if result.I[0][-1] == 0:
            I_2_mean = 0


        cF_vals[ieta,meas] = I_mean
        dyn_vals[ieta,meas] = I_2_mean

        #pl.plot(L_sis.time,L_sis.I)
        #pl.plot(np.arange(len(result.I[0])),result.I[0])
        #pl.title("$R_0 = %4.2f$ meas = %d" %(R0[ieta], meas))

        #pl.show()
        

        seed += 1
        print meas


cF_means = cF_vals.mean(axis=1)
dyn_means = dyn_vals.mean(axis=1)

print R0
print etas

pl.plot(R0,cF_means/N,label='cF')
pl.plot(R0,dyn_means/N,label='dyn')
pl.plot(R0[R0>1],1-1/R0[R0>1])
pl.xscale("log")
#pl.yscale("log")
pl.legend()


pl.show()

