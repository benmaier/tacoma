import numpy as np
import matplotlib.pyplot as pl
import pickle as pickle
from scipy.special import kl_div
from rocsNWL import marker_sequence


def get_distribution_from_counter(c,m=None):
    if m is None:
        m = max(c)
    s = float(sum(c.values()))

    keys = np.arange(1,m+1,dtype=float)
    hist = np.zeros_like(keys,dtype=float)
    for key in keys:
        k = int(key)
        hist[k-1] = c[k] / s

    return keys, hist


def get_size_distribution_kl_div(counter_test,counter_orig,symmetric=True):
    test_max = max(counter_test)
    orig_max = max(counter_orig)
    sum_test = float(sum(counter_test.values()))
    sum_orig = float(sum(counter_orig.values()))

    keys = np.arange(1,max(test_max,orig_max)+1,dtype=int)

    x = []
    y = []

    for key in keys:
        x.append(counter_orig[key] / sum_orig)
        y.append(counter_test[key] / sum_test)

    result = kl_div(x,y)

    if symmetric:
        result += kl_div(y,x)
        result /= 2.

    return result


def log_least_squares(y1,y2,offset_punishment_for_zero = 1):
    y1 = np.array(y1)
    y2 = np.array(y2)

    global_min = min(
                        np.min(y1[y1>0.]),
                        np.min(y2[y2>0.])
                    )

    result = 0.
    indices = np.where(np.logical_and(y1>0.,y2>0.))

    result += np.sum((np.log(y1[indices])-np.log(y2[indices]))**2)
    indices = np.where(np.logical_and(y1>0.,y2==0.))
    result += np.sum((np.log(y1[indices])-(np.log(global_min)-offset_punishment_for_zero))**2)
    indices = np.where(np.logical_and(y1==0.,y2>0.))
    result += np.sum((np.log(y2[indices])-(np.log(global_min)-offset_punishment_for_zero))**2)

    result /= len(y1)

    return result

def plot_dtu(ax):

    orig_fn = 'dtu_group_size_distribution_1_weeks.pickle'

    with open(orig_fn,'r') as f:
        orig_counter = pickle.load(f)

    orig_x, orig_hist = get_distribution_from_counter(orig_counter,m=412)
    orig_mean = orig_x.dot(orig_hist)
    orig_std = np.sqrt(orig_hist.dot( orig_x**2 - orig_mean**2 ))

    ax.plot(orig_x, orig_hist,
                    ls = 'None',
                    marker = '.',
                    ms = 3,
                    mfc = 'None',
                    mew = 0.8,
                    label='DTU'
                    )

