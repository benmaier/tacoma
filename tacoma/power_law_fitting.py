# -*- coding: utf-8 -*-


from __future__ import print_function
import numpy as np

from scipy.special import zeta
from scipy.optimize import fmin
from scipy.optimize import newton

def sum_until_value_threshold(f,start_at,thresh):
    
    root_function = lambda x: np.abs(f(x)) - thresh

    n_max = newton(root_function,start_at+10)
    n = np.arange(start_at,int(n_max)+1,dtype=float)
    return np.sum(f(n))
    

def fit_power_law_clauset(data,
                          x_min = None,
                          x_max = None,
                          discrete = False,
                         ):
    """after https://arxiv.org/pdf/0706.1062.pdf"""
    
    x = np.array(data,dtype=float)
    if x_max is not None:
        x = x[x<=x_max]
    if x_min is None:
        x_min = x.min()
    else:
        x = x[x>=x_min]
    n = float(len(x))

    if not discrete:
        alpha = 1.0 + n/(np.log(x/x_min)).sum()
        err = (alpha-1.0) / np.sqrt(n)
    else:
        sum_ln_x = np.log(x).sum()
        negative_likelihood = lambda alpha: n*np.log(zeta(alpha,x_min)) + alpha*sum_ln_x
        minimum = fmin(negative_likelihood,10,disp=False)
        alpha = minimum[0]
        z = zeta(alpha,x_min)
        z_1 = sum_until_value_threshold(lambda x: -np.log(x+x_min)/(x+x_min)**alpha,
                                        start_at=x_min,
                                        thresh=1e-6)
        z_2 = sum_until_value_threshold(lambda x: +np.log(x+x_min)**2/(x+x_min)**alpha,
                                        start_at=x_min,
                                        thresh=1e-6)
        err = (n * (z_2/z - (z_1/z)**2) )**(-1.0)
 
    return alpha, err, x_min

def power_law(alpha,x_min,discrete=False):
    if discrete:
        return lambda x: x**(-alpha)/zeta(alpha,x_min)
    else:
        return lambda x: (x/xmin)**(-alpha)/x_min*(alpha-1.0)

if __name__ == "__main__":
    import networkx as nx
    from collections import Counter
    import matplotlib.pyplot as pl

    G = nx.barabasi_albert_graph(100000,2)
    k = [ G.degree(n) for n in G.nodes() ] 
    c = Counter(k)
    dat = np.array([(a,b) for a,b in c.items()])
    x = dat[:,0]
    y = dat[:,1]/float(dat[:,1].sum())
    alpha, err, x_min = fit_power_law_clauset(k,discrete=True)
    print(alpha, err)
    pl.figure()
    pl.plot(x,y,'.')
    pl.xscale('log')
    pl.yscale('log')
    pl.xlabel('degree $k$')
    pl.ylabel('probability')
    pl.plot(x,power_law(alpha,x_min,discrete=True)(x))

    pl.show()
