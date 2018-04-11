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

N = 100

print "simulating"
result = cF.ZSBB_model([],N,1.0,0.7,0.7,1000000,seed=1346235)
print "done"

fig, ax = pl.subplots(1,3,figsize=(12,4))

size_count = Counter()
size_count[1] = N
this_count = Counter(size_count)

m_edges = 0
edges_in = result.edges_in
edges_out = result.edges_out
ks = []
for e_in,e_out in izip(edges_in,edges_out):
    m_in = len(e_in)
    m_out = len(e_out)
    m_edges += m_in
    m_edges -= m_out
    mean_k = m_edges / (2.*N)
    ks.append(mean_k)


t = np.array(result.t,dtype=float)
t /= N

new_t = t[::N]
new_k = ks[::N]



ax[2].plot(new_t,new_k,'-')




pl.show()
