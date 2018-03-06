import tacoma as tc
import numpy as np
import matplotlib.pyplot as pl


N = 10
t_run = 10000

result = tc.dynamic_RGG(N=N,t_run_total=t_run,mean_link_duration=10.,
                            periodic_boundary_conditions_for_link_building = False,
                            record_sizes_and_durations = True,
                            critical_density = 1.0,
                            #verbose = True)
                            seed = 2335
                            )

this = tc.edge_lists(result)

second_result = tc.measure_group_sizes_and_durations_for_edge_lists(this,
                                                     #verbose=True,
                                                     )

#print result.durations
#print second_result.contact_durations
#print result.size_histograms
#print second_result.size_histograms

print "all size histograms are the same:", all([ hist1 == hist2 for hist1, hist2 in zip(result.size_histograms, second_result.size_histograms)])
print "all durations are the same:", all([ float(dur1) == float(dur2) for dur1, dur2 in zip(result.group_durations,second_result.contact_durations)])


for size, group_dur in enumerate(second_result.group_durations):
    if len(group_dur) > 0:
        unique, counts = np.unique(group_dur, return_counts=True)
        print counts
        pl.plot(unique, counts,'-.',label='$N=%d$' % (size,))


pl.xscale('log')
pl.yscale('log')


pl.show()
