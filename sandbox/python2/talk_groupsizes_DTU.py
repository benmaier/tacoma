import matplotlib
import matplotlib.pyplot as pl

import tacoma as tc
from tacoma.drawing import draw_edges
from tacoma.analysis import temporal_network_group_analysis
import matplotlib.pyplot as pl

from tacoma.model_conversions import estimate_ZSBB_args
from tacoma.model_conversions import estimate_flockwork_P_args
from tacoma.model_conversions import estimate_dynamic_RGG_args

import numpy as np

SIZE = 10 
matplotlib.rc('font', size=SIZE,family='Arial')
matplotlib.rc('axes', titlesize=SIZE)
matplotlib.rc('mathtext', default='regular') 

# ======== get original network =============
socio = tc.load_json_taco("~/.tacoma/dtu_1_weeks.taco")
socio_binned = tc.bin(socio,dt=300)

# ========= plot properties ==============
socio_result = tc.measure_group_sizes_and_durations(socio)
fig, ax, data = temporal_network_group_analysis(socio_result,time_normalization_factor = 1/3600.,time_unit='h')

new_fig, new_ax = pl.subplots(1,1,figsize=(4.5,3))

x_group, y_group = data['size_histogram']
new_ax.plot(x_group,y_group,'o',
                ms=4,
                mew=1,
                mfc='None',
                label = 'DTU one week',
            )
new_ax.legend()
new_ax.set_xlim([0.8,200])
new_ax.set_ylim([1e-4,1e3])
new_ax.set_xscale("log")
new_ax.set_yscale("log")
new_ax.set_xlabel('group size m')
new_ax.set_ylabel('mean number of m-sized groups')
new_fig.tight_layout()
new_fig.savefig("talk/dtu_1.pdf")


# ============== generate surrogate network from flockwork_P model ==============
#fwP_params = estimate_flockwork_P_args(socio_binned,dt=1800.,aggregated_network=socio_result.aggregated_network)
fwP_params = estimate_flockwork_P_args(socio_binned,dt=3600.)
#tc.write_fwP_args(fwP_params,"./fwP_args_dtu_1_weeks.json")
#fwP = tc.flockwork_P_varying_rates_neighbor_affinity(**fwP_params)
seed = np.random.randint(100000)
print seed
seed = 88566
fwP = tc.flockwork_P_varying_rates(seed=seed,**fwP_params)
fwP_result = tc.measure_group_sizes_and_durations(fwP)
fig, ax, data = temporal_network_group_analysis(fwP_result,time_unit='h',time_normalization_factor = 1/3600.)
x_group, y_group = data['size_histogram']
new_ax.plot(x_group,y_group,'*',
                c = tc.color_sequence[4], 
                ms=5,
                mew=1,
                mfc='None',
                label = 'flockwork original',
            )
new_ax.legend()
new_ax.set_xlim([0.8,200])
new_ax.set_ylim([1e-4,1e3])
new_ax.set_xscale("log")
new_ax.set_yscale("log")
new_fig.tight_layout()
new_fig.savefig("talk/dtu_2.pdf")

fwP_binned = tc.bin(fwP,dt=300)

fwP_result = tc.measure_group_sizes_and_durations(fwP_binned)
fig, ax, data = temporal_network_group_analysis(fwP_result,time_unit='h',time_normalization_factor = 1/3600.)
x_group, y_group = data['size_histogram']
new_ax.plot(x_group,y_group,'d',
                c = tc.color_sequence[3], 
                ms=4,
                mew=1,
                mfc='None',
                label = r'flockwork binned $\Delta t=5$min',
            )
new_ax.legend()
new_ax.set_xlim([0.8,200])
new_ax.set_ylim([1e-4,1e3])
new_ax.set_xscale("log")
new_ax.set_yscale("log")
new_fig.tight_layout()
new_fig.savefig("talk/dtu_3.pdf")


pl.show()
