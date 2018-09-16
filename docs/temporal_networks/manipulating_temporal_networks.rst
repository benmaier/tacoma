Manipulating temporal networks
==============================

Crucial functionalities needed for working with temporal networks are any which
can also be found in typical video editing software. We want to be able
to cut slices out of them, concatenate different networks, and to rescale time.

Slicing
~~~~~~~

Say we have an experiment conducted between :math:`8\mathrm{h}\leq t < 48\mathrm{h}`.
We are only interested in the activity during the first night, though, i.e. for
:math:`20\mathrm{h}\leq t < 32\mathrm{h}`.

We can do

.. code:: python

    night = tc.slice(temporal_network, new_t0 = 20, new_tmax = 32)

The corresponding function can be found at :class:`tacoma.api.slice`.

Concatenating
~~~~~~~~~~~~~

We might have data of two consecutive weeks, but in two files. We can load both files,
construct ``temporal_network_A`` and ``temporal_network_B`` and then concatenate them
by

.. code:: python

    temporal_network = tc.concatenate([temporal_network_A, temporal_network_B])

.. attention::
    
    - All temporal networks in the list must be of equal type.
    - Temporal networks do not have to have the same number of nodes. However,
      if :math:`N_A` and :math:`N_B` differ, it is useful to provide 
      `int_to_node`-dictionaries containing unique identifiers for nodes.
      When concatenating, these dictionaries will be merged and node integers
      will be remapped accordingly.
    - All times :math:`t` of a temporal network will be remapped to have the
      beginning of the experiment at :math:`t_\mathrm{max}` of the predecessing
      network in the list.
      
The corresponding function can be found at :class:`tacoma.api.concatenate`.

Rescaling time
~~~~~~~~~~~~~~

We can rescale time to either speed up or slow down all events, or to change
the unit of time. For example, to change from seconds to hours, do

.. code:: python

    ec_hour = tc.rescale_time(ec_seconds,
                              new_t0 = ec_seconds.t0 / 3600.,
                              new_tmax = ec_seconds.tmax / 3600.)

To speed up the time by a factor 2, you could do

.. code:: python

    ec_speedy = tc.rescale_time(ec,
                                new_t0 = ec.t0,
                                new_tmax = ec.tmax / 2.)

The corresponding function can be found at :class:`tacoma.api.rescale_time`.          

