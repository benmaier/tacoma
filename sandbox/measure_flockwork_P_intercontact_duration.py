import numpy as np
import cFlockwork as cF
import matplotlib.pyplot as pl
import scipy.sparse as sprs
import time

dtu_week = cF.dtu_week()
E = []
P = dtu_week.P
gamma = dtu_week.gamma
N = dtu_week.N
tmax = dtu_week.tmax
t_run_total = tmax
seed = 45

edge_changes = cF.flockwork_P_varying_rates(E,N,P,t_run_total,gamma,tmax,seed=seed)

print "measuring groups and durations, ignore histograms differences =", True
start = time.time()
result = cF.measure_group_sizes_and_durations_for_edge_changes(edge_changes,ignore_size_histogram_differences=True)
end = time.time()
print "took", end - start, "seconds"


new_fig, axes = pl.subplots(2,2)
axes = axes.flatten()

for group_size in range(1,6):
    y, x = np.histogram(result.group_durations[group_size],bins=100)
    x = 0.5*(x[1:] + x[:-1])
    axes[0].plot(x,y,'.')
    axes[0].set_xscale('log')
    axes[0].set_yscale('log')

    
size_histogram = np.array([ (size,sum(durations)) \
                            for size, durations in enumerate(result.group_durations)\
                            if durations > 0 ])
axes[1].plot(size_histogram[:,0],size_histogram[:,1],'.')
axes[1].set_xscale('log')
axes[1].set_yscale('log')

social_network = result.aggregated_network

i, j , data = [], [], []
[ ( i.append(k[0]), j.append(k[1]), data.append(np.log(v))) for k, v in social_network.iteritems() ]

A = sprs.csr_matrix((data,(i,j)),shape=(N,N))
A += A.T

axes[2].imshow(A.todense())

y, x = np.histogram(np.array(result.contact_durations)/3600.,bins=200)
x = 0.5*(x[1:] + x[:-1])
axes[3].plot(x,y,'.')
axes[3].set_xscale('log')
axes[3].set_yscale('log')

pl.show()

