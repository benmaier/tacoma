import tacoma as tc
import matplotlib.pyplot as pl
import DynGillEpi as gill
import numpy as np


N_meas = 10
L_ = tc.dynamic_RGG(50,t_run_total=10,mean_link_duration=2)
L = tc.edge_lists()
L.copy_from(L_)
etas = np.logspace(0.1,1,10)
eta = 10.
rho = 0.5
tc_vals = np.zeros((len(etas),N_meas))
dyn_vals = np.zeros((len(etas),N_meas))

seed = 123434543

for ieta,eta in enumerate(etas):
    for meas in range(N_meas):
        L_sis = tc.Dyn_SIS(L.N, 100, eta, rho,number_of_initially_infected = 25,seed=seed)
        tc.gillespie_SIS_on_edge_lists(L,L_sis)

        t = np.array(L_sis.time)
        dt = t[1:] - t[:-1]
        I_mean = np.array(dt).dot(L_sis.I[:-1])/sum(dt)

        result = gill.SIS_Poisson_homogeneous(L.N,L.edges,eta,rho,100,seed=seed,initial_number_of_infected=25)
        tc_vals[ieta,meas] = I_mean
        dyn_vals[ieta,meas] = np.mean(np.array(result.I[0][:-1],dtype=float))

        pl.plot(L_sis.time,L_sis.I)
        pl.plot(np.arange(len(result.I[0])),result.I[0])

        pl.show()
        

        seed += 1
        print(meas)


tc_means = tc_vals.mean(axis=1)
dyn_means = dyn_vals.mean(axis=1)

pl.plot(etas,tc_means)
pl.plot(etas,dyn_means)

pl.show()

