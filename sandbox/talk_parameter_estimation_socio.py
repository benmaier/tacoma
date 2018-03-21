import matplotlib
import matplotlib.pyplot as pl

import tacoma as tc
from tacoma.drawing import draw_edges
from tacoma.analysis import temporal_network_group_analysis
from tacoma.model_conversions import estimate_ZSBB_args
from tacoma.model_conversions import estimate_flockwork_P_args
from tacoma.model_conversions import estimate_dynamic_RGG_args

import numpy as np

SIZE = 10 
matplotlib.rc('font', size=SIZE,family='Arial')
matplotlib.rc('axes', titlesize=SIZE)
matplotlib.rc('mathtext', default='regular') 


# ======== get original network =============
socio = tc.load_sociopatterns_hypertext_2009()
socio_binned = tc.bin(socio,dt=20.)

fig_k, ax_k = pl.subplots(1,1,figsize=(4,3))
t, k = tc.mean_degree(socio_binned)
ax_k.plot(t/3600.,k,'-',lw=2)
ax_k.set_xlabel("time t [h]")
ax_k.set_ylabel(r"mean degree $\left\langle k\right\rangle$")
fig_k.tight_layout()
fig_k.tight_layout()
fig_k.savefig('talk/mean_degree_a.pdf')

# ========= plot properties ==============
socio_result = tc.measure_group_sizes_and_durations(socio)

axs = []
figs = []
for i in range(4):
    fig, ax = pl.subplots(1,1,figsize=(4,3))
    axs.append(ax)
    figs.append(fig)

__, _, data_s = temporal_network_group_analysis(socio_result,time_unit='h',
        time_normalization_factor = 1/3600.,
        ax = axs,
        )
#figs.tight_layout()
traj = tc.get_edge_trajectories(socio)
draw_edges(traj.trajectories,ax=axs[3],time_unit='h',
        time_normalization_factor = 1/3600.,
        )
axs[0].set_ylabel(r'number of $m$-sized groups')

for i, fig in enumerate(figs):
    fig.tight_layout()
    fig.tight_layout()
    fig.savefig('talk/socio_'+str(i)+'a.pdf',dpi=300)

# ============== generate surrogate network from flockwork_P model ==============
fwP_params = estimate_flockwork_P_args(socio_binned,dt=(socio_binned.tmax-socio_binned.t[0])/120.,
                                       aggregated_network = socio_result.aggregated_network
                                      )
seed = np.random.randint(100000)
#seed = 53220
#seed = 64346
seed = 17407
print seed
fwP = tc.flockwork_P_varying_rates_neighbor_affinity(seed=seed,**fwP_params)
#fwP = tc.flockwork_P_varying_rates(seed=seed,**fwP_params)
fwP_binned = tc.bin(fwP,dt=20.)

tc.write_json_taco(fwP_binned,"~/.tacoma/fw_binned_ht09.taco")
tc.write_json_taco(socio_binned,"~/.tacoma/binned_ht09.taco")
#fwP_binned = fwP

result = tc.measure_group_sizes_and_durations(fwP_binned)
fig, ax, data = temporal_network_group_analysis(result,time_unit='h',
        time_normalization_factor = 1/3600.,
        )
#fig.tight_layout()
fig, ax = pl.subplots(1,1,figsize=(4,3))
traj = tc.get_edge_trajectories(fwP_binned)
draw_edges(traj.trajectories,ax=ax,time_unit='h',
        time_normalization_factor = 1/3600.,
        )
ax.set_ylim([0,2200])


fig.tight_layout()
fig.savefig('talk/socio_edges_FW.pdf',dpi=300)

fig_dt, ax_dt = pl.subplots(1,1,figsize=(4,3))
ts = np.array([fwP.t0] + fwP.t + [fwP.tmax])
dt = ts[1:] - ts[:-1]
x, y = tc.get_logarithmic_histogram(dt,100)


ax_dt.plot(x,y,'s',mfc='None')
ax_dt.set_xscale('log')
ax_dt.set_yscale('log')
ax_dt.set_xlabel(r'inter-event times $\Delta t$ [s]')
ax_dt.set_ylabel(r'probability density [1/s]')

fig_dt.tight_layout()
fig_dt.savefig("talk/distrubution_inter_event_times.pdf")


print data_s.keys()
x_size, y_size = data['size_histogram']
xc, yc = data['contact']
xic, yic = data['inter-contact']


axs[0].plot(x_size, y_size,'s',mfc='None',c=tc.color_sequence[2],ms=6)
axs[1].plot(xc, yc,'s',mfc='None',ms=3)
axs[1].plot(xic, yic,'s',mfc='None',ms=3)

for size in range(2,6):
    xgr, ygr = data[size]
    #axs[2].plot(xgr,ygr,ls='',marker=tc.marker_sequence[size-2],c=tc.color_sequence[size-2],mfc='None',ms=3)
    axs[2].plot(xgr,ygr,ls='',marker=tc.marker_sequence[size],c=tc.color_sequence[size+2],mfc='None',ms=6)

for i, fig in enumerate(figs):
    fig.tight_layout()
    fig.savefig('talk/socio_'+str(i)+'b.pdf',dpi=300)

t, k = tc.mean_degree(fwP_binned)
ax_k.plot(t/3600.,k,'-',lw=2)
fig_k.tight_layout()
fig_k.tight_layout()
fig_k.savefig('talk/mean_degree_b.pdf')

#figs.savefig('socio_2.pdf',dpi=300)

gamma = np.array(fwP_params['rewiring_rate'] + [(fwP_params['tmax'], fwP_params['rewiring_rate'][-1][-1])] )
t, gamma = gamma[:,0], gamma[:,1]
P = np.array(fwP_params['P'] + [fwP_params['P'][-1]])

fig_p, ax_p = pl.subplots(1,3,figsize=(10,3))

ax_p[0].step(t/3600.,gamma*3600,where='post')
ax_p[1].step(t/3600.,P,where='post')
ax_p[2].step(t/3600.,gamma*P*3600,where='post')

ax_p[0].set_xlabel('time [h]')
ax_p[1].set_xlabel('time [h]')
ax_p[2].set_xlabel('time [h]')
ax_p[0].set_ylabel('event rate per node $\gamma$ [1/h]')
ax_p[1].set_ylabel('probability to reconnect P')
ax_p[2].set_ylabel(r'reconnection rate per node $\gamma\times P$ [1/h]')

fig_p.tight_layout()
fig_p.savefig('talk/gamma_P.pdf')

pl.show()
