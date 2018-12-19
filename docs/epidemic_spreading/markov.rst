Dynamic Markov integration
==========================

In contrast to Gillespie's simulation one can write down an approximate equation
for the temporal evolution of the nodes' probability to be infected, (see Eq. (1) in
`Epidemic Threshold in Continuous-Time Evolving Networks`_) which can then be
integrated. This is slower than the Gillespie simulation but leads to smaller 
fluctuations which makes it easier to obtain the epidemic threshold.

How it works
------------

Let's denote the probability that node `u` is infected at time `t` as :math:`p_u(t)`.
Further, we assume that these probabilities are uncorrelated and that the maximum
integration time step :math:`\Delta t_\mathrm{max}` is small. For an 
SIS-process with infection rate :math:`\eta` and recovery rate :math:`\varrho`, the temporal evolution is then
governed by

.. math::

    p_u(t+\Delta t) = p_u(t) (1-\varrho\Delta t) + (1-p_u(t))
                  \left[
                        1-\prod_{v=1}^N(1-\eta\Delta t A_{vu}(t)p_v(t))
                  \right]

Here, :math:`\Delta t=\mathrm{min}\{\Delta t_{\mathrm{network\ change}}, \Delta t_{\mathrm{max}}\}`
is either the time until the network changes next, or the maximally allowed 
integration time step :math:`\Delta t_\mathrm{max}`, whichever is smaller.

How it's implemented
--------------------

In the C++-core :mod:`_tacoma`, the Markov function is implemented to
be supplied with a ``temporal_network`` and a markov integration
``model``.

The interface between the Markov integration function and a hypothetical 
Markov model ``model`` is defined with the following functions.

- ``model.update_network`` : This is used by the Markov integration function
  to tell the model that the network has been updated. The model
  then has to update its states internally. 
- ``model.step`` : The Markov integration function passes the current time
  and the desired :math:`\Delta t` for the integration 
  to this function. The model then has perform the integration step
  according to its rules.

Further necessary functions include

- ``model.simulation_ended`` : ask the model whether the simulation
  is essentially over (if the accumulated probability is smaller
  than some defined threshold ``model.minimum_I``)
- ``model.update_observables`` : ask the model to update its
  internal observals because the simulation is about to end
- ``model.print`` : a function to offer a status check for the
  Gillespie function's verbose argument
- ``model.reset`` : the possibility to wind back the model, e.g.
  for using the same model instance for a second simulation.

How to integrate for an SIS-Model
---------------------------------

First, generate a Markov SIS-object.

.. code:: python

    mv_SIS = tc.MARKOV_SIS(
                        N, #numbder of nodes
                        t_run_total, # run time of the integration
                        infection_rate,
                        recovery_rate,
                        minimum_I, # minimal accumulated probability that nodes are infected
                                   #(below that, the integration is stopped)
                        number_of_initally_infected=N//2,
                        seed = seed)

Then, we can integrate and plot the results

.. code:: python

    tc.markov_epidemics(temporal_network, mv_SIS, max_dt=0.01)
    pl.plot(mv_SIS.time, mv_SIS.I)


We can also integrate the equation on a model, on the fly, e.g. the edge activity model

.. code:: python

    AM = tc.EdgeActivityModel(N, # number of nodes
                          k/(N-1.), # network density
                          omega, # edge activity rate
                          save_temporal_network=False)

    tc.markov_epidemics(AM, mv_SIS, max_dt=0.01)
    pl.plot(mv_SIS.time, mv_SIS.I)


How to develop an own model
---------------------------

# TODO (this chapter is a bit complicated)


.. _Epidemic Threshold in Continuous-Time Evolving Networks: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.068302
