# -*- coding: utf-8 -*-
"""
Tacoma stands for TemporAl COntact Modeling and Analysis. 
It provides fast tools to analyze temporal contact networks, 
produce surrogate networks using qualitative models 
and simulate Gillespie processes on them. Additionally 
provides some visualization tools.
"""

__version__ = "0.0.26"

from _tacoma import *

color_sequence = [ '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf' ]
marker_sequence = ['s','d','o','X','v','<','^','.','>','h','p','P','*','8','H']

from .api import _get_raw_temporal_network
from .api import *
from .data_io import *
from .tools import *

# from .load_model_parameters import *
