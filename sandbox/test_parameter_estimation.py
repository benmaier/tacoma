import tacoma as tc
from tacoma.drawing import draw_edges
from tacoma.analysis import temporal_network_group_analysis
import matplotlib.pyplot as pl

from tacoma.model_conversions import estimate_ZSBB_args
from tacoma.model_conversions import estimate_flockwork_P_args
from tacoma.model_conversions import estimate_dynamic_RGG_args


# ======== get original network =============
socio = tc.load_sociopatterns_hypertext_2009()
socio_binned = tc.bin(socio,dt=20.)

# ========= plot properties ==============
socio_result = tc.measure_group_sizes_and_durations(socio)
fig, ax, data = temporal_network_group_analysis(socio_result,time_unit=socio.time_unit)
fig.tight_layout()
traj = tc.get_edge_trajectories(socio)
draw_edges(traj.trajectories,ax=ax[3])

# ========== generate surrogate from ZSBB model ======
ZSBB_params = estimate_ZSBB_args(socio,group_sizes_and_durations=socio_result)
ZSBB_params['t_run_total'] = len(socio_binned.t) * socio_binned.N / 5.0
zsbb = tc.ZSBB_model(**ZSBB_params)
this_t0 = zsbb.t0
zsbb.t0 = 0.
zsbb.tmax -= this_t0
zsbb.t = [ t-this_t0 for t in zsbb.t]

# bin the result because the original network is binned, too
zsbb_binned = tc.bin(zsbb,N_time_steps=int((zsbb.tmax-zsbb.t0) / (socio_binned.N / 5.0)))

# plot properties
result = tc.measure_group_sizes_and_durations(zsbb_binned)
#print result.group_durations[4]
fig, ax, data = temporal_network_group_analysis(result,max_group=3)
fig.tight_layout()
traj = tc.get_edge_trajectories(zsbb_binned)
draw_edges(traj.trajectories,ax=ax[3])

# ============== generate surrogate network from flockwork_P model ==============
fwP_params = estimate_flockwork_P_args(socio_binned,dt=(socio_binned.tmax-socio_binned.t[0])/120.)
fwP = tc.flockwork_P_varying_rates(**fwP_params)
fwP_binned = tc.bin(fwP,dt=20.)

result = tc.measure_group_sizes_and_durations(fwP_binned)
fig, ax, data = temporal_network_group_analysis(result,time_unit=socio.time_unit)
fig.tight_layout()
traj = tc.get_edge_trajectories(fwP_binned)
draw_edges(traj.trajectories,ax=ax[3],time_unit=socio.time_unit)

# ============== generate surrogate network from dynamic RGG model ==============
dRGG_params = estimate_dynamic_RGG_args(socio_binned,group_sizes_and_durations=socio_result)
dRGG = tc.dynamic_RGG(**dRGG_params)

result = tc.measure_group_sizes_and_durations(dRGG)
fig, ax, data = temporal_network_group_analysis(result,time_normalization_factor=20.,time_unit=socio.time_unit)
fig.tight_layout()
traj = tc.get_edge_trajectories(dRGG)
draw_edges(traj.trajectories,ax=ax[3],time_normalization_factor=20.,time_unit=socio.time_unit)


pl.show()
