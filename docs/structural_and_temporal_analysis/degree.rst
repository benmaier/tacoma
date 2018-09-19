Edge counts and degree
======================

The first and simplest network observables to analyze are the
network density (or mean degree or total number of edges)
and degree distribution.

Average degree distribution
---------------------------

To compute the average degree distribution

.. math::

    \overline{P_k} = \frac{1}{t_\mathrm{max}-t_0} 
                     \int\limits_{t_0}^{t_\mathrm{max}}dt\ P_k(t)

you can use the function :func:`tacoma.api.degree_distribution`.
Here's an example

.. code:: python

    >>> P_k = tc.degree_distribution(temporal_network)
    >>> P_k 
    [ 0.0, 0.21, 0.45, 0.1 ]

The `k`-th entry of the list ``P_k`` is the average probability
that a node has degree `k`.

Mean degree
-----------

The time-dependent mean degree :math:`\left\langle k \right\rangle(t)`
is one way of computing the density of
a temporal network and is computed as

.. math::

    \left\langle k \right\rangle(t) = 2m(t)/N.

To compute it, use the function :func:`tacoma.api.mean_degree`

.. code:: python

    >>> t, mean_k = tc.mean_degree(temporal_network)
    >>> t
    [ 0.0, 0.1, 0.25 ]
    >>> mean_k
    [ 1.2, 0.1, 4.7 ]

The mean degree is constant in between changing times. In order
to compute the average mean degree, use :func:`tacoma.tools.time_average`

    >>> average_mean_k = tc.time_average(t, 
                                         mean_k, 
                                         tmax=temporal_network.tmax
                                        )

Edge counts
-----------

Measuring the mean degree 
does not tell you anything about the network activity, 
i.e. the exchange of edges in the network. To find out about 
activity, 
one can measure the total amount of edges being switched off 
and the total
amount of edges being created at events times :math:`t`.
The observables :math:`m(t)`, :math:`m_\mathrm{in}(t)`, 
and :math:`m_\mathrm{out}(t)` can be computed using
the function :func:`tacoma.api.edge_counts` as

.. code:: python

    >>> m, m_in, m_out = tc.edge_counts(temporal_network)

While ``m`` contains the number of edges present in the temporal
network at :math:`t_0` and all event times, ``m_in`` is a list
of counts of all edges being created at all event times
(``m_out`` being a list of counts of all edges leaving the network),
so the length of the list ``m`` is :math:`N_e+1`, while the length 
of the lists ``m_in`` and ``m_out`` is :math:`N_e`.

In data where the inter-event time is large but ``m_in`` and ``m_out``
are typically large, as well, one could interpret both to reflect the
number of edges being created/deleted in the time interval
:math:`(t_i,t_{i+1}]`.

Contact coverage
----------------

An observable to analyze the distribution of edge activity rates in
a temporal network is the total number of uniquely observed edges up
to time :math:`t`, the contact coverage :math:`C(t)`. It can be
computed using the function :func:`tacoma.tools.contact_coverage` as

.. code:: python

    >>> t, C = tc.contact_coverage(temporal_network)

``t`` contains the time points at which ``C`` changes and ``C`` contains
the corresponding count of uniquely observed edges. :math:`C(t)` is 
hence a strict monotonically increasing function.


