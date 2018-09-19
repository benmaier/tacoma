"""
This module provides routines to load and write
temporal networks to the `taco`-fileformat.
"""
from __future__ import print_function

import os
import json
import gzip  # for sociopatterns data
import csv

import wget

import tacoma as tc
from tacoma import _get_raw_temporal_network


def mkdirp_customdir(directory='~/.tacoma/'):
    """simulate `mkdir -p` functionality"""

    directory = os.path.abspath(os.path.expanduser(directory))
    if not os.path.exists(directory):
        os.makedirs(directory)


def write_json_taco(temporal_network, fp):
    """Writes a temporal network to a .taco-file (which is actually in json format).

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        An instance of a temporal network.
    fp : file-like or :obj:`str`
        write to this file

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
        raise ValueError('Unknown temporal network format: ' +
                         str(type(temporal_network)))

    if isinstance(fp, str):
        fp = os.path.abspath(os.path.expanduser(fp))
        fp = open(fp, 'w')

    json.dump(this_data, fp, separators=(',', ':'))

    fp.close()


def load_json_taco(fp):
    """
    Loads a temporal network from a .taco-file (which is actually in json format).

    Parameters
    ----------
    fp : file-like or :obj:`str`
        read from this file

    Returns
    -------
    temporal network
        type as given in the .taco-file
    """

    file_is_string = isinstance(fp, str)

    if file_is_string:
        fp = os.path.abspath(os.path.expanduser(fp))
        fp = open(fp, 'r')

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
        temporal_network.int_to_node = {
            int(i): s for i, s in this_data['int_to_node'].items()}
        temporal_network.notes = this_data['notes']
        temporal_network.time_unit = this_data['time_unit']

    elif this_data['type'] == 'edge_lists':
        temporal_network = tc.edge_lists()
        temporal_network.t = this_data['t']
        temporal_network.tmax = this_data['tmax']
        temporal_network.N = this_data['N']
        temporal_network.edges = this_data['edges']
        temporal_network.int_to_node = {
            int(i): s for i, s in this_data['int_to_node'].items()}
        temporal_network.notes = this_data['notes']
        temporal_network.time_unit = this_data['time_unit']

    else:
        raise ValueError(
            'file is corrupted, unknown temporal network format: ' + this_data['type'])

    return temporal_network


def read_json_taco(fp):
    """Loads a temporal network from a .taco-file (which is actually in json format) 
    by simply calling :mod:`load_json_taco` because I'm too stupid to remember 
    if it's actually 'read' or 'load' smh.

    Parameters
    ----------
    fp : file-like or :obj:`str`
        read from this file

    Returns
    -------
    temporal_network : `edge_lists` or `edge_changes`
        type as given in the .taco-file
    """

    return load_json_taco(fp)


def download_and_convert_sociopatterns_hypertext_2009(url="http://www.sociopatterns.org/files/datasets/003/ht09_contact_list.dat.gz",
                                                      filename="~/.tacoma/ht09.taco",
                                                      ):
    """Download the SocioPatterns 'Hypertext 2009 dynamic contact network' data,
    extract it and save it as taco. This data is actually binned in intervals
    of `[t-20s, t]`.

    Parameters
    ----------
    url : :obj:`str`, optional
        The url from which the tsv-data should be retrieved
    filename : :obj:`str`, optional
        this is the path where the taco will be saved to. default : "~/.tacoma/ht09.taco"

    Returns
    -------
    edge_lists : :mod:`edge_lists`
        The temporal network of the 'Hypertext 2009 dynamic contact network'.

    Notes
    -----

    If you use this data, please cite

    ::
        L. Isella et al.,  What's in a crowd? Analysis of face-to-face behavioral networks, 
        Journal of Theoretical Biology 271, 166 (2011).
    """

    # get directory name for download
    directory, _ = os.path.split(os.path.abspath(os.path.expanduser(filename)))
    mkdirp_customdir(directory)

    # download
    wget.download(url, out=directory)

    # open gzipped file
    gzip_file = os.path.join(directory, 'ht09_contact_list.dat.gz')
    with gzip.open(gzip_file, mode='rt') as f:
        reader = csv.reader(f, delimiter='\t')

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
            # this is to account for the interval choice [t-20s, t]
            t = float(int(row[0]) - 20)

            # if the time changed, we save the new time and
            # prepare to save new edges
            if t_old != t:
                if (t_old is not None) and (t - t_old > 20):
                    edges.append([])
                    time.append(t_old+20)

                edges.append([])
                time.append(t)

            # get the edge
            i = int(row[1])
            j = int(row[2])

            # map the edge to integers
            if i not in node_to_int:
                this_int = len(node_to_int)
                node_to_int[i] = len(node_to_int)
                int_to_node[this_int] = str(i)

            if j not in node_to_int:
                this_int = len(node_to_int)
                node_to_int[j] = len(node_to_int)
                int_to_node[this_int] = str(j)

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
    el.int_to_node = int_to_node

    # verifying that this is a valid temporal network
    tc.verify(el)

    # save this edge_lists instance
    with open(os.path.abspath(os.path.expanduser(filename)), 'w') as f:
        write_json_taco(el, f)

    # remove the downloaded gzipped file
    os.remove(gzip_file)

    return el


def download_and_convert_sociopatterns_high_school_2013(url="http://www.sociopatterns.org/wp-content/uploads/2015/07/High-School_data_2013.csv.gz",
                                                        filename="~/.tacoma/hs13.taco",
                                                        ):
    """Download the SocioPatterns 'High school 2013 dynamic contact network' data,
    extract it and save it as taco. This data is actually binned in intervals
    of `[t-20s, t]`.

    Parameters
    ----------
    url : :obj:`str`, optional
        The url from which the tsv-data should be retrieved
    filename : :obj:`str`, optional
        this is the path where the taco will be saved to. default : "~/.tacoma/hs13.taco"

    Returns
    -------
    edge_lists : :mod:`edge_lists`
        The temporal network of the 'High school 2013 dynamic contact network'.

    Notes
    -----

    If you use this data, please cite

    :: [HS13]
        R. Mastrandrea, J. Fournet, A. Barrat,
        Contact patterns in a high school: a comparison between data collected
        using wearable sensors, contact diaries and friendship surveys.
        PLoS ONE 10(9): e0136497 (2015)
    """

    # get directory name for download

    directory, _ = os.path.split(os.path.abspath(os.path.expanduser(filename)))
    mkdirp_customdir(directory)

    # download
    wget.download(url, out=directory)

    # open gzipped file
    gzip_file = os.path.join(directory, 'High-School_data_2013.csv.gz')
    with gzip.open(gzip_file, mode='rt') as f:
        reader = csv.reader(f, delimiter=' ')

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

        count = 0
        for row in reader:

            if count == 0:
                t0 = int(row[0]) - 20

            # this is to account for the interval choice [t-20s, t]
            t = float(int(row[0]) - 20 - t0)

            # if the time changed, we save the new time and
            # prepare to save new edges
            if t_old != t:
                if (t_old is not None) and (t - t_old > 20):
                    edges.append([])
                    time.append(t_old+20)

                edges.append([])
                time.append(t)

            # get the edge
            i = int(row[1])
            j = int(row[2])

            # map the edge to integers
            if i not in node_to_int:
                this_int = len(node_to_int)
                node_to_int[i] = len(node_to_int)
                int_to_node[this_int] = str(i)

            if j not in node_to_int:
                this_int = len(node_to_int)
                node_to_int[j] = len(node_to_int)
                int_to_node[this_int] = str(j)

            # save the edge
            edges[-1].append(tuple(sorted([
                node_to_int[i],
                node_to_int[j]
            ])))
            t_old = t

            count += 1

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

        In this data, t0 = 0.0 corresponds to UNIX time """ + str(t0) + """.

        For more info, please visit 
        http://www.sociopatterns.org/datasets/high-school-contact-and-friendship-networks/ .

        If you use this data, please cite

        R. Mastrandrea, J. Fournet, A. Barrat,
        Contact patterns in a high school: a comparison between data collected
        using wearable sensors, contact diaries and friendship surveys.
        PLoS ONE 10(9): e0136497 (2015)
        """
    el.int_to_node = int_to_node

    # verifying that this is a valid temporal network
    tc.verify(el)

    # save this edge_lists instance
    with open(os.path.abspath(os.path.expanduser(filename)), 'w') as f:
        write_json_taco(el, f)

    # remove the downloaded gzipped file
    os.remove(gzip_file)

    return el


def write_fwP_args(args, filename):
    """Dump Flockwork-P arguments to a json-file"""

    filename = os.path.abspath(os.path.expanduser(filename))

    with open(filename, 'w') as f:
        json.dump(args, f)


def load_fwP_args(filename):
    """Load Flockwork-P arguments from a json-file"""

    filename = os.path.abspath(os.path.expanduser(filename))

    with open(filename, 'r') as f:
        args = json.load(f)

    return args


def load_sociopatterns_hypertext_2009(filename="~/.tacoma/ht09.taco"):
    """Once :func:`tacoma.data_io.download_sociopatterns_hypertext_2009` was called,
    use this function to retrieve an :mod:`edge_lists` instance
    of the conference data set 'Hypertext 2009 dynamic contact network'
    (from the SocioPatterns project).

    Parameters
    ----------
    filename : :obj:`str`, optional
        this is the path where the taco was saved to. default : "~/.tacoma/ht09.taco"

    Returns
    -------
    edge_lists : :mod:`edge_lists`
        The temporal network of the 'Hypertext 2009 dynamic contact network'.

    If you use this data, please cite

    ::
        L. Isella et al.,  What's in a crowd? Analysis of face-to-face behavioral networks, 
        Journal of Theoretical Biology 271, 166 (2011).
    """

    filename = os.path.abspath(os.path.expanduser(filename))

    if not os.path.exists(filename):
        raise ValueError(
            "File "+filename+" does not exist. Have you called `tacoma.download_and_convert_sociopatterns_hypertext_2009()` before?")

    return load_json_taco(filename)


def load_sociopatterns_high_school_2013(filename="~/.tacoma/hs13.taco"):
    """Once :func:`tacoma.data_io.download_sociopatterns_high_school_2013` was called,
    use this function to retrieve an :mod:`edge_lists` instance
    of the conference data set 'High school 2013 dynamic contact network'
    (from the SocioPatterns project).

    Parameters
    ----------
    filename : :obj:`str`, optional
        this is the path where the taco was saved to. default : "~/.tacoma/hs13.taco"

    Returns
    -------
    edge_lists : :mod:`edge_lists`
        The temporal network of the 'High school 2013 dynamic contact network'.

    Notes
    -----
    If you use this data, please cite

    ::
        R. Mastrandrea, J. Fournet, A. Barrat,
        Contact patterns in a high school: a comparison between data collected using wearable sensors, contact diaries and friendship surveys.
        PLoS ONE 10(9): e0136497 (2015)
    """

    filename = os.path.abspath(os.path.expanduser(filename))

    if not os.path.exists(filename):
        raise ValueError(
            "File "+filename+" does not exist. Have you called tacoma.download_and_convert_sociopatterns_high_school_2013() before?")

    return load_json_taco(filename)


def write_edge_trajectory_coordinates(temporal_network, filename, filter_for_duration=0.0):
    """Write the coordinates of the edge activation periods to a json-file
    such that each entry corresponds to a line to be drawn.

    Parameters
    ----------
    temporal_network : :mod:`edge_changes`, :mod:`edge_lists`, :mod:`edge_changes_with_histograms`, or :mod:`edge_lists_with_histograms`
        An instance of a temporal network.
    filename : :obj:`str`
        Write to this file.
    """

    traj = tc.get_edge_trajectories(temporal_network)

    try:
        t0 = temporal_network.t[0]
    except:
        t0 = temporal_network.t0
    tmax = temporal_network.tmax

    coords = []
    for i_edge, entry in enumerate(traj):
        for time_pair in entry.time_pairs:
            t_i = time_pair[0]
            t_f = time_pair[1]
            if t_f - t_i > filter_for_duration:
                coords.append((i_edge, t_i, t_f))

    data = {
        'xlim': (t0, tmax),
        'ylim': (0, len(traj)),
        'data': coords,
    }

    filename = os.path.abspath(os.path.expanduser(filename))
    with open(filename, 'w') as f:
        json.dump(data, f)


def load_json_dict(fn):
    """Load a dictionary from a JSON-file"""
    with open(fn, 'r') as f:
        this_dict = json.load(f)

    return this_dict


if __name__ == "__main__":
    el = download_and_convert_sociopatterns_hypertext_2009()

    print(el.N)
    print(el.tmax)
    print(el.t[:10])
    print(el.edges[:10])
    print(el.notes)
    print(el.time_unit)
