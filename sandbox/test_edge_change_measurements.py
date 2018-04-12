import tacoma as tc
import numpy as np
import matplotlib.pyplot as pl


N = 50
t_sim = 1000*N
t_eq = 0*N
t_total = t_sim + t_eq
lambda_ = 0.7
b0 = 0.6
b1 = 0.8

plot_size = False

print("simulating")
result = tc.ZSBB_model([],N,lambda_,b0,b1,t_sim,t_equilibration=t_eq,seed=1346,record_sizes_and_durations=True)
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


#print result.durations
#print second_result.contact_durations
#print result.size_histograms
#print second_result.size_histograms

#print "all size histograms are the same:", all([ hist1 == hist2 for hist1, hist2 in zip(result.size_histograms, second_result.size_histograms)])
print("all durations are the same:", all([ float(dur1) == float(dur2) for dur1, dur2 in zip(result.contact_durations,second_result.contact_durations)]))

fig, ax = pl.subplots(1,2,figsize=(10,5))

for size, group_dur in enumerate(second_result.group_durations):
    if len(group_dur) > 0:
        unique, counts = np.unique(group_dur, return_counts=True)
        #print counts
        ax[1].plot(unique, counts,'-.',label='$N=%d$' % (size,))

        unique2, counts2 = np.unique(result.group_durations[size], return_counts=True)
        ax[1].plot(unique2, counts2,'-.',label='$N=%d$' % (size,))

        #if size == 1:
        #    print [ (u1,u2) for u1,u2 in zip(unique,unique2)]
        #    print [ c1 == c2 for c1,c2 in zip(counts,counts2)]
        #    print [ (c1,c2) for c1,c2 in zip(counts,counts2)]
        print("size =", size, "; all counts equal: ", all([ c1 == c2 for c1,c2 in zip(counts,counts2)]))


ax[1].set_xscale('log')
ax[1].set_yscale('log')

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



x_ = lambda _: 0.5*(_[1:]+_[:-1])
y, x = np.histogram(result.contact_durations,bins=200)
ax[0].plot(x_(x),y,'.')
y, x = np.histogram(second_result.contact_durations,bins=200)
ax[0].plot(x_(x),y,'.')
y, x = np.histogram(third_measurement_durations,bins=200)
ax[0].plot(x_(x),y,'.')
ax[0].set_xscale('log')
ax[0].set_yscale('log')



pl.show()
