# -*- coding: utf-8 -*-
"""
Tacoma stands for TemporAl COntact Modeling and Analysis. 
It provides fast tools to analyze temporal contact networks, 
produce surrogate networks using qualitative models 
and simulate Gillespie processes on them. Additionally 
provides some visualization tools.
"""

__version__ = "0.0.20"

from _tacoma import *

color_sequence = [ u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf' ]
marker_sequence = ['s','d','o','X','v','<','^','.','>','h','p','P','*','8','H']

from .api import _get_raw_temporal_network
from .api import *
from .data_io import *
from .tools import *
from .load_model_parameters import *
