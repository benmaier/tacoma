import os 
import json
import gzip # for sociopatterns data
import csv

import tacoma as tc
from tacoma import _get_raw_temporal_network

import wget

def mkdirp_customdir(directory='~/.tacoma/'):
    directory = os.path.abspath( os.path.expanduser(directory) )
    if not os.path.exists(directory):
        os.makedirs(directory)

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
        this_data['int_to_node'] = temporal_network.int_to_node
        this_data['notes'] = temporal_network.notes
        this_data['time_unit'] = temporal_network.time_unit

    elif type(temporal_network) == tc.edge_lists:
        this_data['type'] = 'edge_lists'
        this_data['t'] = temporal_network.t
        this_data['tmax'] = temporal_network.tmax
        this_data['N'] = temporal_network.N
        this_data['edges'] = temporal_network.edges
        this_data['int_to_node'] = temporal_network.int_to_node
        this_data['notes'] = temporal_network.notes
        this_data['time_unit'] = temporal_network.time_unit

    else:
        raise ValueError('Unknown temporal network format: ' + str(type(temporal_network)))

    if isinstance(fp,basestring):
        fp = os.path.abspath( os.path.expanduser(fp) )
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
        fp = os.path.abspath( os.path.expanduser(fp) )
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
        temporal_network.int_to_node = this_data['int_to_node']
        temporal_network.notes = this_data['notes']
        temporal_network.time_unit = this_data['time_unit']

    elif this_data['type'] == 'edge_lists':
        temporal_network = tc.edge_lists()
        temporal_network.t = this_data['t']
        temporal_network.tmax = this_data['tmax']
        temporal_network.N = this_data['N']
        temporal_network.edges = this_data['edges']
        temporal_network.int_to_node = this_data['int_to_node']
        temporal_network.notes = this_data['notes']
        temporal_network.time_unit = this_data['time_unit']

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

def download_and_convert_sociopatterns_hypertext_2009(url="http://www.sociopatterns.org/files/datasets/003/ht09_contact_list.dat.gz",
                                                      filename="~/.tacoma/ht09.taco",
                                                      ):
    """
        Download the SocioPatterns 'Hypertext 2009 dynamic contact network' data,
        extract it and save it as taco. This data is actually binned in intervals
        of [t-20s, t].

        py::arg(url)      -- this is the url to the gzipped data
        py::arg(filename) -- this is the path where the taco should be saved to.
        
        If you use this data, please cite

        L. Isella et al.,  What's in a crowd? Analysis of face-to-face behavioral networks, 
        Journal of Theoretical Biology 271, 166 (2011).
    """

    # get directory name for download

    directory, single = os.path.split(os.path.abspath(os.path.expanduser(filename)))
    mkdirp_customdir(directory)

    # download
    wget.download(url,out=directory)

    # open gzipped file
    gzip_file = os.path.join(directory, 'ht09_contact_list.dat.gz')
    with gzip.open(gzip_file,'rb') as f:
        reader = csv.reader(f,delimiter='\t')

        # mappings of nodes to integers
        node_to_int = {}
        int_to_node = {}

        # get an initial t_old
        # (this is done to detect changes in the tsv
        t_old = None

        # list of edge lists
        edges = []

        # time points
        time = []
        for row in reader:
            t = float( int(row[0]) - 20 ) #this is to account for the interval choice [t-20s, t]

            # if the time changed, we save the new time and 
            # prepare to save new edges
            if t_old != t:
                edges.append([])
                time.append(t)

            # get the edge
            i = int(row[1])
            j = int(row[2])

            # map the edge to integers
            if i not in node_to_int:
                this_int = len(node_to_int)
                node_to_int[i] = len(node_to_int)
                int_to_node[this_int] = i

            if j not in node_to_int:
                this_int = len(node_to_int)
                node_to_int[j] = len(node_to_int)
                int_to_node[this_int] = j

            # save the edge
            edges[-1].append(tuple(sorted([
                                    node_to_int[i],
                                    node_to_int[j]
                                    ])))
            t_old = t

        N = len(node_to_int)
        tmax = time[-1] + 20.0


    # get a new `edge_lists` instance
    el = tc.edge_lists()

    el.N = N
    el.tmax = tmax
    el.edges = edges
    el.t = time
    el.time_unit = 's'
    el.notes = """
        This data is binned.

        In this data, t0 = 0.0 corresponds to 8am on Jun 29th 2009 (UNIX time 1246255200).

        For more info, please visit http://www.sociopatterns.org/datasets/hypertext-2009-dynamic-contact-network/ .

        If you use this data, please cite

        L. Isella et al.,  What's in a crowd? Analysis of face-to-face behavioral networks, 
        Journal of Theoretical Biology 271, 166 (2011).
        """

    # verifying that this is a valid temporal network
    tc.verify(el)

    # save this edge_lists instance
    with open(os.path.abspath(os.path.expanduser(filename)),'w') as f:
        write_json_taco(el,f)

    # remove the downloaded gzipped file
    os.remove(gzip_file)

    return el

def load_sociopatterns_hypertext_2009(filename="~/.tacoma/ht09.taco"):
    """
        Once `download_sociopatterns_hypertext_2009` was called,
        use this function to retrieve an `edge_lists` instance
        of the conference data set 'Hypertext 2009 dynamic contact network'
        (from the SocioPatterns project).

        py::arg(filename) -- this is the path where the taco was saved too
                             (default: ~/.tacoma/ht09.taco)
        
        If you use this data, please cite

        L. Isella et al.,  What's in a crowd? Analysis of face-to-face behavioral networks, 
        Journal of Theoretical Biology 271, 166 (2011).
    """

    filename = os.path.abspath(os.path.expanduser(filename))

    return load_json_taco(filename)

if __name__ == "__main__":
    el = download_and_convert_sociopatterns_hypertext_2009()

    print el.N
    print el.tmax
    print el.t[:10]
    print el.edges[:10]
    print el.notes
    print el.time_unit




