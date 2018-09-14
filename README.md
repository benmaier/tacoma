![logo](logo/new_logo_grey.png)

TemporAl COntact Modeling and Analysis. Provides fast tools to analyze temporal contact networks, 
produce surrogate networks using qualitative models and simulate Gillespie processes on them.

## Quick example

In order to download the SocioPatterns 
['Hypertext 2009'-dataset](http://www.sociopatterns.org/datasets/hypertext-2009-dynamic-contact-network/) 
and visualize it interactively, do the following.

```python
>>> import tacoma as tc
>>> from tacoma.interactive import visualize
>>> temporal_network = tc.download_and_convert_sociopatterns_hypertext_2009()
100% [..............................................................................] 67463 / 67463
>>> visualize(temporal_network, frame_dt = 20)
```

![visualization example](https://github.com/benmaier/tacoma/raw/master/img/ht09_example.gif)

## What is tacoma?

`tacoma` is a joint C++/Python-package for the modeling and analysis of undirected and 
unweighted temporal networks, with a focus on (but not limited to) human face-to-face contact networks.

### Pros of using tacoma

* networks are natively described in continuous time
* two main native formats to describe temporal networks (`tc.edge_lists` and `tc.edge_changes`),
  a third way, a sorted list of `on`-intervals for each list called `tc.edge_trajectories` is
  available, but algorithms work on the two native formats only
* the simple portable file-format `.taco` as a standardized way to share temporal network data
  (which is just the data dumped to a `.json`-file, a simple file format readable from a
  variety of languages)
* easy way to simulate Gillespie (here, epidemic spreading) processes on temporal networks
* easy framework to develop new Gillespie-simulations algorithms on temporal networks
* multiple and simple ways to interactively visualize temporal networks
* simple functions to manipulate temporal networks (slice, concatenate, rescale time, sample, bin, convert)
* simple functions to analyze structural and statistical properties of temporal networks
  (mean degree, degree distribution, group size distribution, group life time distributions, etc.)
* fast algorithms due to C++-core (_fast_ as in _faster than pure Python_)
* relatively fast and easy to compile since it only depends on the C++11-stdlib 
  and [https://github.com/pybind/pybind11](pybind11) without the large overhead of `Boost`

### Cons of using tacoma

* no support for directed temporal networks yet
* no support for weighted temporal networks yet

## Install

If you get compiling errors, make sure that [https://github.com/pybind/pybind11](pybind11) is installed.

    $ git clone https://github.com/benmaier/tacoma
    $ pip install pybind11
    $ pip install ./tacoma

## Examples

Check out the [sandbox directory](https://github.com/benmaier/tacoma/tree/master/sandbox])
