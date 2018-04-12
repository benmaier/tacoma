import tacoma as tc
import numpy as np
import matplotlib.pyplot as pl


N = 100

dtu_week = tc.dtu_week()
E = []
P = dtu_week.P
gamma = dtu_week.gamma
tmax = dtu_week.tmax
t_run_total = tmax
seed = 45
plot_size = False

print("simulating")
result = tc.flockwork_P_varying_rates(E,N,P,t_run_total,gamma,tmax,seed=seed)
print("done")


this = tc.edge_changes()
this.t = result.t
this.N = result.N
this.edges_initial = result.edges_initial
this.t0 = result.t0
this.tmax = result.tmax
this.edges_in = result.edges_in
this.edges_out = result.edges_out
print(len(this.t), len(this.edges_in))
print("first time point: ", this.t[0])

second_result = tc.measure_group_sizes_and_durations_for_edge_changes(
                                                     this,
                                                     #verbose=True,
                                                     )
print("intial edges", result.edges_initial)

third_measurement_durations = []
initial_edges = set(this.edges_initial)
current_edges = {}
for edge in initial_edges:
    current_edges[edge] = this.t0

for t_, e_in, e_out in zip(this.t,this.edges_in,this.edges_out):
    
    for edge in e_in:
        current_edges[edge] = t_

    for edge in e_out:
        if edge in initial_edges:
            initial_edges.erase(edge)
        else:
            third_measurement_durations.append(t_ - current_edges[edge])

new_fig, ax = pl.subplots(1,2)
ax = ax.flatten()


print(second_result.contact_durations[:10], third_measurement_durations[:10])

x_ = lambda _: 0.5*(_[1:]+_[:-1])
y, x = np.histogram(second_result.contact_durations,bins=200)
ax[0].plot(x_(x),y,'.')
y, x = np.histogram(third_measurement_durations,bins=200)
ax[0].plot(x_(x),y,'.')
ax[0].set_xscale('log')
ax[0].set_yscale('log')



pl.show()
