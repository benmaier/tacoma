import tacoma as tc
from tacoma.drawing import draw_edges
from tacoma.analysis import temporal_network_group_analysis
import matplotlib.pyplot as pl

from tacoma.model_conversions import estimate_ZSBB_args
from tacoma.model_conversions import estimate_flockwork_P_args
from tacoma.model_conversions import estimate_dynamic_RGG_args

import time
import numpy as np


# ======== get original network =============
socio = tc.load_json_taco("~/.tacoma/dtu_1_weeks.taco")
socio_binned = tc.bin(socio,dt=300)

socio_result = tc.measure_group_sizes_and_durations(socio)

# ============== generate surrogate network from flockwork_P model ==============
fwP_params = estimate_flockwork_P_args(socio_binned,dt=3600.,aggregated_network=socio_result.aggregated_network)
fwP = tc.flockwork_P_varying_rates_neighbor_affinity(**fwP_params)

fwP_params = estimate_flockwork_P_args(socio_binned,dt=3600.)
fwP = tc.flockwork_P_varying_rates(**fwP_params)

fwP_binned = tc.bin(fwP,dt=300)

N = fwP.N

R0 = 1.0
rho = 3.0 / (24*3600.)
dt = 3600.
t_simulation = 2*fwP.tmax
t_sample = np.arange(int(t_simulation / dt)+1,dtype=float) * dt

N_meas = 30

fig, ax = pl.subplots(1,3,figsize=(12,4))

for tn in [socio, fwP, fwP_binned]:
    start = time.time()
    t, k = np.array(tc.mean_degree(tn))
    end = time.time()

    print "took", end-start, "seconds"

    line, = ax[0].step(t,k,where='post',lw=1)

    mean_k = tc.time_average(t, k)
    print mean_k
    eta = R0 * rho / mean_k

    
    i_sample = np.zeros_like(t_sample)

    for meas in range(N_meas):

        sis = tc.SIS(N, t_simulation, eta, rho, number_of_initially_infected = 10)

        tc.gillespie_SIS(tn, sis)

        t = np.array(sis.time)
        i = np.array(sis.I,dtype=float) / N

        this_sample = tc.sample_a_function(t,i,t_sample)

        i_sample += this_sample / N_meas
        
        ax[2].plot(t_sample,this_sample,c=line.get_color(),alpha=0.1)

    ax[1].plot(t_sample,i_sample)




pl.show()
