Epidemic classes
================

Simulations of epidemic spreading in `tacoma` work by providing
epidemic compartmental classes to the adapted Gillespie SSA class.
State transitions rates are defined within the epidemic classes,
the simulation is then performed in the Gillespie class, where
the observables are written back to the epidemic class instances.

SI
--
In a susceptible-infected dynamic, nodes can be in two compartments,
the susceptible compartment `S` and the infected compartment `I`.
While only nodes are in compartments, it is a link-based reaction
process where links between an `S`-node and an `I` node transition
to links between an `I`-node and another `I`-node with infection rate
:math:`\eta` by turning a susceptible to an infected. The corresponding
reaction equation is

.. math::

    S + I \stackrel{\eta}{\longrightarrow} I + I.

A susceptible-infected dynamic is initialized in the following way, using
the class :class:`_tacoma.SI`

.. code:: python

    SI = tc.SI(N, #number of nodes
               t_simulation, # maximum time of the simulation
               infection_rate,
               number_of_initially_infected = int(N), # optional, default: 1
               number_of_initially_vaccinated = 0, # optional, default: 0
               seed = 792, # optional, default: randomly initiated
               save_infection_events, # optional, default: false
              )

The infection rate is the expected number of infection events between a
susceptible-infected pair per unit of time of the temporal network.

The instance of the `SI` class is then passed to the corresponding
Gillespie function :func:`tacoma.api.gillespie_SI` for the simulation.

.. code:: python

    tc.gillespie_SI(temporal_network, SI)

During the simulation, the following observables are written to
the `SI` object.

- ``SI.time`` : A time-ordered list of floats where each entry is a time
  point at which one of the observable changed. In between these
  times the observables are constant.
- ``SI.I``: A list of ints where each entry is the total number of infected
  at the corresponding time in ``SI.time``
- ``SI.SI``: A list of ints where each entry is the total number of
  susceptible-infected contacts at the corresponding time in ``SI.time``
- ``SI.infection_events``: A list of pairs of ints, where each entry is
  the edge along which the infection event occurred at the corresponding time
  in ``SI.time``. Each edge is given in the form ``(infection_source,
  infection_target)``. Only saved if the flag ``save_infection_events``
  is set to `True`.

Plot the results as

.. code:: python

    import matplotlib.pyplot as pl

    pl.step(SI.time, SI.I)

.. note::

    - If the time of the epidemic spreading simulation is larger than
      the duration of the temporal network the network is automatically
      looped.
    - The simulation works on both :class:`_tacoma.edge_lists` and
      :class:`_tacoma.edge_changes`.

SIS
---
In a susceptible-infected-susceptible dynamic, nodes can be in two compartments,
the susceptible compartment `S` and the infected compartment `I`.

Infection reactions are a link-based reaction
process where links between an `S`-node and an `I` node transition
to links between an `I`-node and another `I`-node with infection rate
:math:`\eta` by turning a susceptible to an infected. The corresponding
reaction equation is

.. math::

    S + I \stackrel{\eta}{\longrightarrow} I + I.

Furthermore, nodes can recover with recovery rate :math:`\rho` to
become susceptible again. The corresponding reaction equation is

.. math::

    I \stackrel{\rho}{\longrightarrow} S

An SIS dynamic is initialized in the following way, using
the class :class:`_tacoma.SIS`

.. code:: python

    SIS = tc.SIS(N, #number of nodes
                 t_simulation, # maximum time of the simulation
                 infection_rate,
                 recovery_rate,
                 number_of_initially_infected = int(N), # optional, default: 1
                 number_of_initially_vaccinated = 0, # optional, default: 0
                 seed = 792, # optional, default: randomly initiated
                )

The infection rate is the expected number of infection events between a
single susceptible-infected pair per unit of time of the temporal network.
The recovery rate is the expected number of recovery events of a single node
per unit of time of the temporal network.

The instance of the `SIS` class is then passed to the corresponding
Gillespie function :func:`tacoma.api.gillespie_SIS` for the simulation.

.. code:: python

    tc.gillespie_SIS(temporal_network, SIS)

During the simulation, the following observables are written to
the `SIS` object.

- ``SIS.time`` : A time-ordered list of floats where each entry is a time
  point at which one of the observable changed. In between these
  times the observables are constant.
- ``SIS.I``: A list of ints where each entry is the total number of infected
  at the corresponding time in ``SIS.time``
- ``SIS.R0``: A list of floats where each entry is the basic
  reproduction number at the corresponding time in ``SIS.time``. The basic
  reproduction number is computed as
  :math:`R_0 = \left\langle k\right\rangle (t) \eta / \rho`.
- ``SIS.SI``: A list of ints where each entry is the total number of
  susceptible-infected contacts at the corresponding time in ``SIS.time``

Plot the results as

.. code:: python

    import matplotlib.pyplot as pl

    pl.step(SIS.time, SIS.I)

.. note::

    - If the time of the epidemic spreading simulation is larger than
      the duration of the temporal network the network is automatically
      looped.
    - The simulation works on both :class:`_tacoma.edge_lists` and
      :class:`_tacoma.edge_changes`.

SIR
---
In a susceptible-infected-recovered dynamic,
nodes can be in three compartments,
the susceptible compartment `S`, the infected compartment `I`,
and the recovered compartment `R`. Recovered notes cannot
take part in any reaction anymore.

Links between an `S`-node and an `I` node transition
to links between an `I`-node and another `I`-node with infection rate
:math:`\eta` by turning a susceptible to an infected. The corresponding
reaction equation is

.. math::

    S + I \stackrel{\eta}{\longrightarrow} I + I.

Furthermore, nodes can recover with recovery rate :math:`\rho` to
become recovered (or removed). The corresponding reaction equation is

.. math::

    I \stackrel{\rho}{\longrightarrow} R

An SIR dynamic is initialized in the following way, using
the class :class:`_tacoma.SIR`

.. code:: python

    SIR = tc.SIR(N, #number of nodes
                 t_simulation, # maximum time of the simulation
                 infection_rate,
                 recovery_rate,
                 number_of_initially_infected = int(N), # optional, default: 1
                 number_of_initially_vaccinated = 0, # optional, default: 0
                 seed = 792, # optional, default: randomly initiated
                )

The infection rate is the expected number of infection events between a
single susceptible-infected pair per unit of time of the temporal network.
The recovery rate is the expected number of recovery events of a single node
per unit of time of the temporal network.

The instance of the `SIR` class is then passed to the corresponding
Gillespie function :func:`tacoma.api.gillespie_SIR` for the simulation.

.. code:: python

    tc.gillespie_SIR(temporal_network, SIR)

During the simulation, the following observables are written to
the `SIR` object.

- ``SIR.time`` : A time-ordered list of floats where each entry is a time
  point at which one of the observable changed. In between these
  times the observables are constant.
- ``SIR.I``: A list of ints where each entry is the total number of infected
  at the corresponding time in ``SIR.time``
- ``SIR.R``: A list of ints where each entry is the total number of recovered
  at the corresponding time in ``SIR.time``
- ``SIR.R0``: A list of floats where each entry is the basic
  reproduction number at the corresponding time in ``SIR.time``. The basic
  reproduction number is computed asR
  :math:`R_0 = \left\langle k\right\rangle (t) \eta / \rho`.
- ``SIR.SI``: A list of ints where each entry is the total number of
  susceptible-infected contacts at the corresponding time in ``SIR.time``

Plot the results as

.. code:: python

    import matplotlib.pyplot as pl

    pl.step(SIR.time, SIR.I)
    pl.step(SIR.time, SIR.R)

.. note::

    - If the time of the epidemic spreading simulation is larger than
      the duration of the temporal network the network is automatically
      looped.
    - The simulation works on both :class:`_tacoma.edge_lists` and
      :class:`_tacoma.edge_changes`.


SIRS
----
In a susceptible-infected-recovered-susceptible dynamic,
nodes can be in three compartments,
the susceptible compartment `S`, the infected compartment `I`,
and the recovered compartment `R`. Recovered notes can now lose
their immunity with waning immunity rate :math:`\omega`.
The reaction equation is

.. math::

    R \stackrel{\omega}{\longrightarrow} S

Links between an `S`-node and an `I` node transition
to links between an `I`-node and another `I`-node with infection rate
:math:`\eta` by turning a susceptible to an infected. The corresponding
reaction equation is

.. math::

    S + I \stackrel{\eta}{\longrightarrow} I + I.

Furthermore, nodes can recover with recovery rate :math:`\rho` to
become recovered (or removed). The corresponding reaction equation is

.. math::

    I \stackrel{\rho}{\longrightarrow} R

An SIR dynamic is initialized in the following way, using
the class :class:`_tacoma.SIRS`

.. code:: python

    SIRS = tc.SIRS(N, #number of nodes
                   t_simulation, # maximum time of the simulation
                   infection_rate,
                   recovery_rate,
                   waning_immunity_rate,
                   number_of_initially_infected = int(N), # optional, default: 1
                   number_of_initially_vaccinated = 0, # optional, default: 0
                   seed = 792, # optional, default: randomly initiated
                  )

The infection rate is the expected number of infection events between a
single susceptible-infected pair per unit of time of the temporal network.
The recovery rate is the expected number of recovery events of a single node
per unit of time of the temporal network.
The waning immunity is the expected number of events of a single recovered
becoming susceptible per unit of time of the temporal network.

The instance of the `SIRS` class is then passed to the corresponding
Gillespie function :func:`tacoma.api.gillespie_SIRS` for the simulation.

.. code:: python

    tc.gillespie_SIRS(temporal_network, SIRS)

During the simulation, the following observables are written to
the `SIRS` object.

- ``SIRS.time`` : A time-ordered list of floats where each entry is a time
  point at which one of the observable changed. In between these
  times the observables are constant.
- ``SIRS.I``: A list of ints where each entry is the total number of infected
  at the corresponding time in ``SIRS.time``
- ``SIRS.R``: A list of ints where each entry is the total number of recovered
  at the corresponding time in ``SIRS.time``
- ``SIRS.R0``: A list of floats where each entry is the basic
  reproduction number at the corresponding time in ``SIRS.time``. The basic
  reproduction number is computed as
  :math:`R_0 = \left\langle k\right\rangle (t) \eta / \rho`.
- ``SIRS.SI``: A list of ints where each entry is the total number of
  susceptible-infected contacts at the corresponding time in ``SIRS.time``

Plot the results as

.. code:: python

    import matplotlib.pyplot as pl

    pl.step(SIRS.time, SIRS.I)
    pl.step(SIRS.time, SIRS.R)

.. note::

    - If the time of the epidemic spreading simulation is larger than
      the duration of the temporal network the network is automatically
      looped.
    - The simulation works on both :class:`_tacoma.edge_lists` and
      :class:`_tacoma.edge_changes`.
