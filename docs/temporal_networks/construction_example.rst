A construction example
======================

Networks generated in pure python can be easily
described as a `tacoma` temporal network. In a very simple model,
the status of a temporal network changes to a new random network
every :math:`\left\langle \tau \right\rangle` where :math:`\tau`
follows an exponential distribution.

.. code:: python

    import numpy as np
    import tacoma as tc
    import networkx as nx

    # static structure parameters
    N = 100
    mean_degree = 1.5
    p = mean_degree / (N-1.0)

    # temporal parameters
    edge_lists = []
    mean_tau = 1.0
    t0 = 0.0
    tmax = 100.0
    t = []
    this_time = t0

    # Generate a new random network 
    while this_time < tmax:
        G = nx.fast_gnp_random_graph(N, p)
        these_edges = list(G.edges())
        t.append(this_time)
        edge_lists.append(these_edges)

        this_time += np.random.exponential(scale=1/mean_tau)

    # save to _tacoma-object
    el = tc.edge_lists()

    el.N = N
    el.t = t
    el.edges = edge_lists
    el.tmax = tmax

    print("Number of mistakes:", tc.verify(el))

