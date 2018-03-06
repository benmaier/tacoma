from _tacoma import *

from _tacoma import edge_changes as ec
from _tacoma import edge_lists as el
from _tacoma import edge_lists_with_histograms as el_h
from _tacoma import edge_changes_with_histograms as ec_h

color_sequence = [ u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf' ]
marker_sequence = ['s','d','o','X','v','<','^','.','>','h','p','P','*','8','H']

def measure_group_sizes_and_durations(temporal_network):

    if type(temporal_network) == ec_h:
        temporal_network = ec(temporal_network)
        is_ec_format = True
    elif type(temporal_network) == el_h:
        temporal_network = el(temporal_network)
        is_ec_format = False

    if type(temporal_network) == ec:
        result = tc.measure_group_sizes_and_durations_for_edge_changes(temporal_network)
    elif type(temporal_network) == el:
        result = tc.measure_group_sizes_and_durations_for_edge_lists(temporal_network)
    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    return result
    
