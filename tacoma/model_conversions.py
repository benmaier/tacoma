from tacoma import edge_changes as ec
from tacoma import edge_lists as el
from tacoma import edge_lists_with_histograms as el_h
from tacoma import edge_changes_with_histograms as ec_h
import tacoma as tc
from tacoma.analysis import get_logarithmic_histogram

def estimate_ZSBB_parameters(temporal_network,result=None):


    if type(temporal_network) == ec_h:
        temporal_network = ec(temporal_network)
        is_ec_format = True
    elif type(temporal_network) == el_h:
        temporal_network = el(temporal_network)
        is_ec_format = False

    if result is None:
        if type(temporal_network) == ec:
            result = tc.measure_group_sizes_and_durations_from_edge_changes(temporal_network)
        elif type(temporal_network) == el:
            result = tc.measure_group_sizes_and_durations_from_edge_lists(temporal_network)

    
    get_logarithmic_histogram()

    

