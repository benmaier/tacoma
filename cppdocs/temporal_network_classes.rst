Temporal network classes
========================

Undirected and unweighted temporal networks are composed of :math:`N` nodes
and up to :math:`m_{\mathrm{max}}=N(N+1)/2` edges, where each edge :math:`(i,j)` can be
described as a series of events where the edge is either switched on 
or switched off. One way of expressing that is to define the temporal
adjacency matrix

.. math::
    A_{ij}(t) = \begin{cases} 1, & (i,j)\ \mathrm{connected\ at\ time\ } t\\
                              0, & \mathrm{else}.
                \end{cases}

In `tacoma`, we will interpret temporal networks as if they were recorded in an experiment.
We expect that over the course of time :math:`t_0\leq t < t_\mathrm{max}` in which we
record activity, we will encounter :math:`N` nodes from the node set 
:math:`V={0,1,\dots,N-1}` (nodes posses an integer label).

The experiment begins at time :math:`t_0`, where the network consists of an 
edge set :math:`E_0 \subseteq \{i,j: V\times V, i<j\}`. Then, each time the network
changes, we denote that time by an entry in a time vector :math:`t`. Each entry
in the time vector corresponds to a network change event and thus to a change in the edge set.
We call the total number of change events :math:`N_e`, such that the vector :math:`t` has
:math:`N_e` entries.
In between consecutive 
times, the network is constant. After the last recorded event, we kept the experiment running
until the maximum time :math:`t_\mathrm{max}` without observing any change and stopped recording
at :math:`t_\mathrm{max}`.

There's three data structures implemented in this package, all of which capture the situation
described above in different ways and are useful in different situations.

Edge lists
~~~~~~~~~~

The class :class:`_tacoma.edge_lists` consists of a collection of complete edge lists,
each time the network changes, a complete edge list of the network after the change is saved.
It has the following attributes.

- :math:`N` : The total number of nodes
- :math:`t` : A vector of length :math:`N_e+1`. The 0-th entry contains the time of the beginning of the
  experiment :math:`t_0`
- `edges` : A vector of length :math:`N_e+1` where each entry contains an edge list, describing the 
  network after the change which occured at the corresponding time in :math:`t`. 
  The 0-th entry contains the edge list of the beginning of the experiment :math:`t_0`
- :math:`t_\mathrm{max}` : The time at which the experiment ended.

Additionally, 


Edge changes
~~~~~~~~~~~~


The class :class:`_tacoma.edge_changes` consists of a collection of both edges being created
and edges being deleted.
It has the following attributes.

- :math:`N` : The total number of nodes.
- :math:`t_0` : The time of the beginning of the experimen.
- `edges_initial` : The edge list of the beginning of the experiment at :math:`t_0`.
- :math:`t` : A vector of length :math:`N_e`, each time corresponding to a change in the network.
- :math:`t_\mathrm{max}` : The time at which the experiment ended.

Additionally, 


Edge trajectories
~~~~~~~~~~~~~~~~~


