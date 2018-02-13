import cFlockwork as cF
import matplotlib.pyplot as pl
from collections import Counter
from itertools import izip
import numpy as np

def get_hist_from_counter(c):

    data = np.array(c.items(),dtype=float)

    x = data[:,0]
    y = data[:,1] / data[:,1].sum()

    return x,y

def get_hist_from_list(l,bins=400):
    y,x = np.histogram(l,bins=bins)
    y = np.array(y,dtype=float)
    y /= y.sum()
    x = 0.5*(x[:1]+x[:-1])
    return x,y

N = 1000
t_sim = 10000*N
t_eq = 10000*N
t_total = t_sim + t_eq
lambda_ = 1.0
b0 = 0.6
b1 = 0.8

plot_size = False

print "simulating"
result = cF.ZSBB_model([],N,lambda_,b0,b1,t_sim,t_equilibration=t_eq,seed=1346,record_sizes_and_durations=True)
print "done"

fig, ax = pl.subplots(1,3,figsize=(12,4))

x,y = get_hist_from_list(np.array(result.contact_durations,dtype=float)/float(N))
#x,y = get_hist_from_counter(Counter(result.contact_durations))
ax[1].plot(x,y,'s')
y2 = (1+x)**(-2*b1-1)
ax[1].plot(x,y2/y2.sum())
x,y = get_hist_from_list(np.array(result.inter_contact_durations,dtype=float)/float(N))
ax[1].plot(x,y,'.')
y2 = (1+x)**(-2*b0-1)
ax[1].plot(x,y2/y2.sum())
ax[1].set_yscale('log')
ax[1].set_xscale('log')


if plot_size:
    size_count = Counter(result.initial_size_histogram)
    this_count = Counter(size_count)

    for ch in result.group_changes:
        this_count += Counter(ch)
        size_count += this_count

    x,y = get_hist_from_counter(size_count)
    ax[0].plot(x,y,'s')
    ax[0].set_yscale('log')
    ax[0].set_xscale('log')

m_edges = 0
ks = []
print len(result.edges_in[0])
print result.edges_in[0]
edges_in = result.edges_in
edges_out = result.edges_out

for e_in,e_out in izip(edges_in,edges_out):
    #print ch, e_in, e_out
    m_in = len(e_in)
    m_out = len(e_out)
    m_edges += m_in
    m_edges -= m_out
    mean_k= 2*m_edges /float(N)
    ks.append(mean_k)


t = np.array(result.t,dtype=float)
t /= N

print len(ks), len(result.t)


ax[2].plot(t,ks,'-')
ax[2].plot(np.array([t_eq,t_eq],dtype=float)/N,[0,0.25],'-')
ax[2].set_xscale('log')




pl.show()
