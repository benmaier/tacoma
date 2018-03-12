import json

import tacoma as tc
from tacoma import _get_raw_temporal_network

def write_json_taco(temporal_network,fp):
    """
        Writes the provided temporal network to a .taco-file 
        (which is actually json)

        py::arg("temporal_network") -- an instance of `edge_changes` or `edge_lists`
        py::arg("fp")               -- a string containing a path or a file-like object
    """

    temporal_network = _get_raw_temporal_network(temporal_network)
    this_data = {}

    if type(temporal_network) == tc.edge_changes:
        this_data['type'] = 'edge_changes'
        this_data['t'] = temporal_network.t
        this_data['t0'] = temporal_network.t0
        this_data['tmax'] = temporal_network.tmax
        this_data['N'] = temporal_network.N
        this_data['edges_initial'] = temporal_network.edges_initial
        this_data['edges_in'] = temporal_network.edges_in
        this_data['edges_out'] = temporal_network.edges_out

    elif type(temporal_network) == tc.edge_lists:
        this_data['type'] = 'edge_lists'
        this_data['t'] = temporal_network.t
        this_data['tmax'] = temporal_network.tmax
        this_data['N'] = temporal_network.N
        this_data['edges'] = temporal_network.edges

    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    if isinstance(fp,basestring):
        fp = open(fp,'w')

    json.dump(this_data,fp)

    fp.close()

def load_json_taco(fp):
    """
        Loads a temporal network from a .taco-file 
        (which is actually json)

        py::arg("fp") -- a string containing a path or a file-like object
    """

    file_is_string = isinstance(fp,basestring)

    if file_is_string:
        fp = open(fp,'r')

        this_data = json.load(fp)

    if file_is_string:
        fp.close()

    if this_data['type'] == 'edge_changes':
        temporal_network = tc.edge_changes()
        temporal_network.t = this_data['t']
        temporal_network.t0 = this_data['t0']
        temporal_network.tmax = this_data['tmax']
        temporal_network.N = this_data['N']
        temporal_network.edges_initial = this_data['edges_initial']
        temporal_network.edges_in = this_data['edges_in']
        temporal_network.edges_out = this_data['edges_out']

    elif this_data['type'] == 'edge_lists':
        temporal_network = tc.edge_lists()
        temporal_network.t = this_data['t']
        temporal_network.tmax = this_data['tmax']
        temporal_network.N = this_data['N']
        temporal_network.edges = this_data['edges']

    else:
        raise ValueError('file is corrupted, unknown temporal network format: ' + this_data['type'])

    return temporal_network

def read_json_taco(fp):
    """
        Loads a temporal network from a .taco-file 
        (which is actually json) by simply calling `load_json_taco` because
        I'm too stupid to remember if it's actually 'read' or 'load' smh.

        py::arg("fp") -- a string containing a path or a file-like object
    """

    return load_json_taco(fp)


