![logo](logo/new_logo_grey.png)

TemporAl COntact Modeling and Analysis. Provides fast tools to analyze temporal contact networks, produce surrogate networks using qualitative models and simulate Gillespie processes on them.

## Quick example

```python
import tacoma as tc
from tacoma.interactive import visualize

# define temporal network as a list of edge changes
temporal_network = tc.edge_changes()
temporal_network.N = 10
temporal_network.edges_initial = [ (0,1), (2,3), (1,7), (3,5), (1,9), (7,2) ]
temporal_network.t0 = 0.0
temporal_network.t = [ 0.8, 2.4 ]
temporal_network.tmax = 3.1
temporal_network.edges_in = [ 
                              [ (0, 5), (3, 6) ], 
                              [ (3, 7), (4, 9), (7, 8) ],
                            ]
temporal_network.edges_out = [ 
                                [ (0, 1) ],
                                [ (2, 3), (3, 6) ],
                             ]

visualize(temporal_network, frame_dt = 0.05)
```

![visualization example](https://github.com/benmaier/tacoma/raw/master/img/tacoma_example.gif)

## Install

If you get compiling errors, make sure that [https://github.com/pybind/pybind11](pybind11) is installed.

    $ git clone https://github.com/benmaier/tacoma
    $ pip install pybind11
    $ pip install ./tacoma

## Some Explanations

## Examples

### Python

    $ python sandbox/meanfieldtest.py
