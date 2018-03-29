"""
This module provides functions related to the flockwork
temporal network model.
"""

import numpy as np
from numpy import random

def flockwork_P_equilibrium_group_size_distribution(N,P):
    """Get the equilibrium group size distribution of a Flockwork-P model
    given node number N and probability to reconnect P.

    Parameters
    ----------
    N : int
        Number of nodes
    P : float
        Probability to reconnect

    Returns
    -------
    numpy.ndarray
        Group size distribution of this configuration. The m-th entry
        of this array contains the expected number of nodes of groups
        of size m.
    """

    N = int(N)
    P = float(P)

    assert N > 2
    assert P >= 0.
    assert P <= 1.

    if P == 0.:
        dist = np.zeros((N+1,))
        dist[1] = N
    elif P == 1.:
        dist = np.zeros((N+1,))
        dist[-1] = 1.
    else:
        dist = [ N*(1.-P) ]

        N_fak = N - np.arange(1,N-1)
        j_fak = ( P * np.arange(1,N-1) - N+1.)
        div = N_fak / j_fak
        cum_product_div = np.cumprod(div)
        for m in range(2,N):
            #dist.append( (-1)**(m%2) * float(N)/float(m) * (P-1.) * np.prod(N_fak[1:m]/j_fak[1:m]) * P**(m-1) )
            dist.append( (-1)**(m%2) * float(N)/float(m) * (P-1.) * cum_product_div[m-2] * P**(m-1) )


        value = (-1)**( N % 2 ) * P
        for j in range(1,N-1):
            value *= float(N-j-1) / ( (P-N+1.) / P + (j-1))
        #value *= P**(N-1)
        dist.append(value)
        
        dist = [0.] + dist

    return np.array(dist)


def flockwork_P_equilibrium_configuration(N,P,shuffle_nodes = True):
    """Get an equilibrium configuration of a Flockwork-P model
    given node number N and probability to reconnect P.

    Parameters
    ----------
    N : int
        Number of nodes
    P : float
        Probability to reconnect
    shuffle_nodes : bool, optional
        Shuffle the node order in which nodes are distributed
        to groups. 'True' is recommended. default : True

    Returns
    -------
    :obj:`list` of :obj:`tuple` of int
        edge list of equilibrium configuration
    numpy.ndarray
        group size counter of this configuration

    """

    dist = flockwork_P_equilibrium_group_size_distribution(N,P)

    dist = dist[1:]

    total_group_number = sum(dist)
    size_dist = dist / total_group_number
    num_components = []
    edges = []

    if shuffle_nodes:
        node_ints = random.permutation(N)
    else:
        node_ints = np.arange(N)

    # in the beginning, theres N nodes left
    nodes_left = N

    # ... and no group counter contains any groups yet
    C_m = np.zeros(N)

    # while there are still nodes left to distribute
    while nodes_left > 0:

        # loop through group sizes in descending order.
        # start with the smallest group size that
        # may contain all of the nodes left
        for m in xrange(nodes_left,0,-1):

            # if the expected number of groups of this size is not zero
            if dist[m-1] > 0. and nodes_left >= m:

                #new_N_m = random.binomial(total_group_number,size_dist[m-1])
                # dist carries the expected number of m-groups in equilibrium,
                # so we draw a number of m groups from a Poisson distribution
                # with mean dist[m-1]
                new_C_m = random.poisson(dist[m-1])

                # if the new number of groups of size m is larger than the previously drawn number
                if new_C_m>C_m[m-1]:
                    # accept the new number and add the difference in group count of this group size
                    delta_C_m = int(new_C_m - C_m[m-1])
                elif nodes_left == 1 and m == 1:
                    # if there's only one node left, just add it to the loners
                    delta_C_m = 1
                else:
                    delta_C_m = 0

                # add the additional number of groups of this size
                for group_instance in xrange(delta_C_m):
                    #for u in xrange(nodes_left-1,nodes_left-m,-1):
                    #    for v in xrange(u-1,nodes_left-m-1,-1):
                    #        edges.append((node_ints[u],node_ints[v]))

                    # add fully connected clusters to the edge set
                    if m > 1:
                        edges.extend([ (node_ints[u],node_ints[v])\
                                        for u in xrange(nodes_left-1,nodes_left-m,-1) \
                                            for v in xrange(u-1,nodes_left-m-1,-1) ])

                    # remove the grouped nodes from the pool of remaining nodes
                    nodes_left -= m

                    # save this group as an accepted group
                    # and increment the group count
                    C_m[m-1] += 1

                    # if there's no nodes left to distribute,
                    # N_m was chosen too large and we abandon
                    if nodes_left <= 0.:
                        break

            # go to next smaller group size

    return edges, np.append(np.array([0.]), C_m)

if __name__ == "__main__":

    N = 10
    P = 0.5

    dist = flockwork_P_equilibrium_group_size_distribution(N,P)

    print dist, sum([ m * h for m,h in enumerate(dist)])
