import cFlockwork as cF
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

print "simulating"
result = cF.ZSBB_model([],N,lambda_,b0,b1,t_sim,t_equilibration=t_eq,seed=1346,record_sizes_and_durations=True)
print "done"

this = cF.edge_changes()
this.t = result.t
this.edges_in = result.edges_in
this.edges_out = result.edges_out
print len(this.t), len(this.edges_in)
print "first time point: ", this.t[0]

second_result = cF.measure_group_sizes_and_durations_for_edge_changes(
                                                     [],
                                                     N,
                                                     0.,
                                                     this,
                                                     #verbose=True,
                                                     )


#print result.durations
#print second_result.contact_durations
#print result.size_histograms
#print second_result.size_histograms

#print "all size histograms are the same:", all([ hist1 == hist2 for hist1, hist2 in zip(result.size_histograms, second_result.size_histograms)])
print "all durations are the same:", all([ float(dur1) == float(dur2) for dur1, dur2 in zip(result.contact_durations,second_result.contact_durations)])


for size, group_dur in enumerate(second_result.group_durations):
    if len(group_dur) > 0:
        unique, counts = np.unique(group_dur, return_counts=True)
        #print counts
        pl.plot(unique, counts,'-.',label='$N=%d$' % (size,))

        unique2, counts2 = np.unique(result.group_durations[size], return_counts=True)
        pl.plot(unique2, counts2,'-.',label='$N=%d$' % (size,))

        #if size == 1:
        #    print [ (u1,u2) for u1,u2 in zip(unique,unique2)]
        #    print [ c1 == c2 for c1,c2 in zip(counts,counts2)]
        #    print [ (c1,c2) for c1,c2 in zip(counts,counts2)]
        print "size =", size, "; all counts equal: ", all([ c1 == c2 for c1,c2 in zip(counts,counts2)])



pl.xscale('log')
pl.yscale('log')


pl.show()
