![logo](logo/new_logo_grey.png)

TemporAl COntact Modeling and Analysis. Provides fast tools to analyze temporal contact networks, produce surrogate networks using qualitative models and simulate Gillespie processes on them.

## Quick example

```ipython
In [1]: import tacoma as tc

In [2]: from tacoma.interactive import visualize

In [3]: temporal_network = tc.download_and_convert_sociopatterns_hypertext_2009()
100% [..............................................................................] 67463 / 67463
In [4]: visualize(temporal_network, frame_dt = 20)
```

![visualization example](https://github.com/benmaier/tacoma/raw/master/img/ht09_example.gif)

## Install

If you get compiling errors, make sure that [https://github.com/pybind/pybind11](pybind11) is installed.

    $ git clone https://github.com/benmaier/tacoma
    $ pip install pybind11
    $ pip install ./tacoma

## Some Explanations

## Examples

### Python

    $ python sandbox/meanfieldtest.py
