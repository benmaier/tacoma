import tacoma as tc
import matplotlib.pyplot as pl
from collections import Counter

import numpy as np

def get_hist_from_counter(c):

    data = np.array(list(c.items()),dtype=float)

    x = data[:,0]
    y = data[:,1] / data[:,1].sum()

    return x,y

def get_hist_from_list(l,bins=400):
    y,x = np.histogram(l,bins=bins)
    y = np.array(y,dtype=float)
    y /= y.sum()
    x = 0.5*(x[:1]+x[:-1])
    return x,y

N = 412
t_sim = 1000*N
t_eq = 1000*N
t_total = t_sim + t_eq
lambda_ = 0.52
b0 = 0.5
b1 = 0.5

plot_size = True

print("simulating")
result = tc.ZSBB_model([],N,lambda_,b0,b1,t_sim,t_equilibration=t_eq,seed=1346,record_sizes_and_durations=True)
print("done")

fig, ax = pl.subplots(2,2,figsize=(8,8))
ax = ax.flatten()

x,y = get_hist_from_list(np.array(result.contact_durations,dtype=float)/float(N))
#x,y = get_hist_from_counter(Counter(result.contact_durations))
ax[1].plot(x,y,'s',label='contacts')
x,y = get_hist_from_list(np.array(result.inter_contact_durations,dtype=float)/float(N))
ax[1].plot(x,y,'s',label='intercontacts')
ax[1].legend()
ax[1].set_xlabel('duration [steps/N]')
ax[1].set_ylabel('count')
#ax[1].plot(x,(1+x)**(-2*b0-1))
ax[1].set_yscale('log')
ax[1].set_xscale('log')

from plot_dtu import plot_dtu



if plot_size:
    size_count = np.array([ 
                            ( size, sum(durations)) for \
                            size, durations in enumerate(result.group_durations)\
                            if len(durations) > 0
                          ],dtype=float)
    x = size_count[:,0]
    y = size_count[:,1]
    norm = y.sum()

    ax[0].plot(x,y/norm,'s')
    ax[0].set_yscale('log')
    ax[0].set_xscale('log')
    ax[0].set_xlabel('group size')
    ax[0].set_ylabel('probability')

plot_dtu(ax[0])

m_edges = 0
ks = []
print(len(result.edges_in[0]))
print(result.edges_in[0])
edges_in = result.edges_in
edges_out = result.edges_out

for e_in,e_out in zip(edges_in,edges_out):
    #print ch, e_in, e_out
    m_in = len(e_in)
    m_out = len(e_out)
    m_edges += m_in
    m_edges -= m_out
    mean_k = 2*m_edges /float(N)
    ks.append(mean_k)


t = np.array(result.t,dtype=float)
t /= N

print(len(ks), len(result.t))


ax[2].plot(t,ks,'-')
ax[2].plot(np.array([t_eq,t_eq],dtype=float)/N,[0,max(ks)],'-')
ax[2].set_xscale('log')
ax[2].set_xlabel('t')
ax[2].set_ylabel(r'$\left\langle k\right\rangle$')

group_durations = result.group_durations

size = 1
max_size = 5
for size in range(size,max_size+1):
    x,y = get_hist_from_list(np.array(group_durations[size],dtype=float)/N,bins=100)
    ax[3].plot(x,y,'.')

ax[3].set_xscale('log')
ax[3].set_yscale('log')

ax[3].set_xlabel('duration [steps/N]')
ax[3].set_ylabel('count')

fig.tight_layout()



pl.show()
