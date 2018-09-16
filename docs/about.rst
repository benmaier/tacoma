
About this project
==================

Tacoma is an acronym for *(T)empor(A)l (CO)ntact (M)odeling and (A)nalysis*. 
It is a joint C++/Python-package for the modeling and analysis
of undirected and unweighted temporal networks, with a focus on (but not
limited to) human face-to-face contact networks.

Quick example
-------------

In order to download the SocioPatterns `‘Hypertext 2009’-dataset`_ and
visualize it interactively, do the following.

.. code:: python

   >>> import tacoma as tc
   >>> from tacoma.interactive import visualize
   >>> temporal_network = tc.download_and_convert_sociopatterns_hypertext_2009()
   100% [........................................................] 67463 / 67463
   >>> visualize(temporal_network, frame_dt = 20)

This is the result:

.. figure:: https://github.com/benmaier/tacoma/raw/master/img/ht09_extensive_example.gif
   :alt: visualization example

Why should I use tacoma?
------------------------

Pros
~~~~

- networks are natively described in continuous time (which includes a description in
  discrete time)
- two main native formats to describe temporal networks
  (:class:`_tacoma.edge_lists` and :class:`_tacoma.edge_changes`), a third way, a sorted
  list of ``on``-intervals for each edge called
  ``tc.edge_trajectories`` is available, but algorithms work on the two
  native formats only
- the simple portable file-format ``.taco`` as a standardized way to
  share temporal network data (which is just the data dumped to a
  ``.json``-file, a simple file format readable from a variety of
  languages)
- easy functions to produce surrogate temporal networks from four
  different models
- easy way to simulate Gillespie (here, epidemic spreading) processes
  on temporal networks
- easy framework to develop new Gillespie-simulations algorithms on
  temporal networks
- multiple and simple ways to interactively visualize temporal networks
- simple functions to manipulate temporal networks (slice, concatenate,
  rescale time, sample, bin, convert)
- simple functions to analyze structural and statistical properties of
  temporal networks (mean degree, degree distribution, group size
  distribution, group life time distributions, etc.)
- fast algorithms due to C++-core (*fast* as in *faster than pure
  Python*)
- relatively fast and easy to compile since it only depends on the
  C++11-stdlib and `pybind11`_ without the large overhead of ``Boost``

Cons
~~~~

-  no support for directed temporal networks yet
-  no support for weighted temporal networks yet

Install
-------

If you get compiling errors, make sure that `pybind11`_ is installed.

::

   $ git clone https://github.com/benmaier/tacoma
   $ pip install ./tacoma

Note that a C++11-compiler has to be installed on the system before
installing ``tacoma``.

.. _‘Hypertext 2009’-dataset: http://www.sociopatterns.org/datasets/hypertext-2009-dynamic-contact-network/
.. _pybind11: https://github.com/pybind/pybind11
