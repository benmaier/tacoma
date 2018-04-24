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

def power_law_log_likelihood_continuous(
                         data,
                         alpha,
                         x_min = None,
                         x_max = None,
                         ):
    """ Get Log-Likelihood of a power-law fit to the given continuous data
    after https://arxiv.org/pdf/0706.1062.pdf .

    Parameters
    ----------
    data : `numpy.array` or :obj:`list` of `float`
        values measured in an experiment
    alpha : `float`
        Fitted exponent of the power law P(x) ~ x^{-alpha}
    x_min : `float`, optional, default: None
        If given, throws away all data below `x_min`.
        Otherwise it's inferred from the data.
    x_max : `float`, optional, default: None
        if given, throws away all data above `x_max`
        Otherwise it's inferred from the data.

    Returns
    -------
    L : `float`
        Log-Likelihood of alpha given the data.
    """

    x = np.array(data,dtype=float)

    if x_max is not None:
        x = x[x<=x_max]
    if x_min is None:
        x_min = x.min()
    else:
        x = x[x>=x_min]

    n = float(len(x))

    L = n * np.log(alpha-1) - np.log(x_min) - alpha * np.log(x/x_min).sum()

    return L

def power_law_log_likelihood_discrete(
                         data,
                         alpha,
                         x_min = None,
                         x_max = None
                         ):
    """ Get Log-Likelihood of a power-law fit to the given discrete data
    after https://arxiv.org/pdf/0706.1062.pdf .

    Parameters
    ----------
    data : `numpy.array` or :obj:`list` of `float`
        values measured in an experiment
    alpha : `float`
        Fitted exponent of the power law P(x) ~ x^{-alpha}
    x_min : `float`, optional, default: None
        If given, throws away all data below `x_min`.
        Otherwise it's inferred from the data.
    x_max : `float`, optional, default: None
        if given, throws away all data above `x_max`
        Otherwise it's inferred from the data.

    Returns
    -------
    L : `float`
        Log-Likelihood of alpha given the data.
    """


    x = np.array(data,dtype=float)

    if x_max is not None:
        x = x[x<=x_max]
    if x_min is None:
        x_min = x.min()
    else:
        x = x[x>=x_min]

    n = float(len(x))

    L = - n*np.log(zeta(alpha,x_min)) - alpha * np.log(x).sum()

    return L

def power_law_log_likelihood_discrete_distribution(
                         x,
                         P,
                         alpha,
                         ):
    """ Get Log-Likelihood of a power-law fit to the given 
    discrete distribution
    after https://arxiv.org/pdf/0706.1062.pdf .

    Parameters
    ----------
    x : `numpy.array` or :obj:`list` of `float`
        mean values of the bins (for discrete data 
        this is just the values of the discrete data).
    P : `numpy.array` or :obj:`list` of `float` 
        The corresponding weights for bin means given 
        in `x` (`P` does not have to be normalized)
    alpha : `float`
        Fitted exponent of the power law P(x) ~ x^{-alpha}

    Returns
    -------
    L : `float`
        Log-Likelihood of alpha given the distribution P(x).
    """

    x = np.array(x,dtype=float)
    P = np.array(P,dtype=float)

    x_min = x.min()

    n = P.sum()
    sum_ln_x = np.log(x).dot(P)

    L = - n * np.log(zeta(alpha,x_min)) - alpha * sum_ln_x

    return L

def fit_power_law_clauset(data,
                          x_min = None,
                          x_max = None,
                          discrete = False,
                         ):
    """ Fits a power-law to the distribution that the given data follows
    after https://arxiv.org/pdf/0706.1062.pdf .

    Parameters
    ----------
    data : `numpy.array` or :obj:`list` of `float`
        values measured in an experiment
    x_min : `float`, optional, default: None
        If given, throws away all data below `x_min`.
        Otherwise it's inferred from the data.
    x_max : `float`, optional, default: None
        if given, throws away all data above `x_max`
        Otherwise it's inferred from the data.
    discrete: `bool`, optional, default: False
        Fit the discrete version of the law using Riemann's
        zeta function.

    Returns
    -------
    alpha : `float`
        The exponent of the power law P(x) ~ x^{-alpha}.
    error : `float`
        The standard error of alpha.
    xmin : `float`
        The inferred minimum value of the data.
    """
    
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

def fit_discrete_power_law_from_distribution(
                          x,
                          P,
                         ):
    """ Fits a power-law to the given distribution
    by a fake redrawing from the given distribution
    after https://arxiv.org/pdf/0706.1062.pdf .

    Parameters
    ----------
    x : `numpy.array` or :obj:`list` of `float`
        mean values of the bins (for discrete data 
        this is just the values of the discrete data).
    P : `numpy.array` or :obj:`list` of `float` 
        The corresponding weights for bin means given 
        in `x` (`P` does not have to be normalized)

    Returns
    -------
    alpha : `float`
        The exponent of the power law P(x) ~ x^{-alpha}.
    error : `float`
        The standard error of alpha.
    xmin : `float`
        The inferred minimum value of the data.
    """
    
   
    x = np.array(x,dtype=float)
    P = np.array(P,dtype=float) 

    x_min = x.min()
    n = P.sum()
    sum_ln_x = np.log(x).dot(P)

    negative_likelihood = lambda alpha: n * np.log(zeta(alpha,x_min)) + alpha * sum_ln_x
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
    dat = np.array(sorted([(a,b) for a,b in c.items()],key=lambda x:x[0]))
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

    alpha2, err2, x_min2 = fit_discrete_power_law_from_distribution(x,y)
    pl.plot(x,power_law(alpha2,x_min2,discrete=True)(x),'--')

    pl.show()
