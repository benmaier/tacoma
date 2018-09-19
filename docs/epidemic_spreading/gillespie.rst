Dynamic Gillespie SSA
=====================

Gillespie's stochastic simulation algorithm is a powerful tool with a successful history
for the simulation of reaction processes in homogeneous media and static networks.
In the following we explain `Vestergaard's and Génois's adaption for temporal networks`_.

How it works
------------

You should check out the original paper linked above for a detailed explanation
of the alorithm, the study is well-presented and accessible.

To sum up, the algorithm depends on event rates, which in turn depend
on the structure of the temporal network and the constitutents of the
stochastic simulation model depending on transition events.

After each event, transition rates :math:`\lambda_i` are updated according to the 
rules of the model and the structure of the temporal network; subsequently
the total rate :math:`\Lambda=\sum_i \lambda_i` is computed. Then, 
a dimensionless time :math:`\tau` is drawn from an exponential distribution.
A dimension is attributed to :math:`\tau` using the total rate :math:`\Lambda`.
If within the expected change of time :math:`\tau` the network structure
changes, :math:`\lambda_i` and :math:`\Lambda` are updated accordingly 
and the remaining time :math:`\Delta \tau` is updated until
:math:`\tau` is reached.

The adpated Gillespie algorithm hence is a neat modular algorithm
in a sense that it solely depends on rates
and how they change when the network structure changes. This implies
that in an implementation we can focus on two parts: A generalized
implementation of the adapted Gillespie algorithm which solely manages
the temporal network and determining the next event. This module can then 
be coupled to a stochastic simulation model which defines the 
transition events and interprets changes
in the network structure to update transition rates.

How it's implemented
--------------------

In the C++-core :mod:`_tacoma`, a Gillespie function is implemented to
be supplied with a ``temporal_network`` and a stochastic simulation
``model``.

The interface between the Gillespie algorithm and a hypothetical 
stochastic model ``model`` is defined with the following functions.

- ``model.update_network`` : This is used by the Gillespie function
  to tell the model that the network has been updated. The model
  then has to update its states internally. This has to be an
  overloaded function with a second version 
  where not only a new edge list is given but lists of edges leaving
  the network and edges being created can be passed, too.
- ``model.get_rates_and_Lambda`` : This function computes a vector
  of event rates and the total rate, based on its internal state.
  The Gillepsie function then uses these objects to determine
  the new time and the event happening at this time.
- ``model.make_event`` : The Gillespie function passes an integer
  to this function which is the index of the corresponding
  event in the rate vector. The model then has to update its
  internal state according to its rules.

Further necessary functions include

- ``model.simulation_ended`` : ask the model whether the simulation
  is essentially over (e.g. if there's no infected anymore)
- ``model.update_observables`` : ask the model to update its
  internal observals because the simulation is about to end
- ``model.print`` : a function to offer a status check for the
  Gillespie function's verbose argument
- ``model.reset`` : the possibility to wind back the model, e.g.
  for using the same model instance for a second simulation.


How to develop an own model
---------------------------

# TODO (this chapter is a bit complicated)


.. _Vestergaard's and Génois's adaption for temporal networks: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004579
