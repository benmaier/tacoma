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

## Install

If you get compiling errors, make sure that [https://github.com/pybind/pybind11](pybind11) is installed.

    $ git clone https://github.com/benmaier/tacoma
    $ pip install pybind11
    $ pip install ./tacoma

## Some functionalities
