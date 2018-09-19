# -*- coding: utf-8 -*-
"""
This module provides the necessary functions to start up a local
HTTP server and open an interactive d3-visualization of one or
several temporal networks.
"""
from __future__ import print_function

import os
import sys

import http.server
import webbrowser
import time
import threading

import wget
import shutil

import json
from distutils.dir_util import copy_tree

import tacoma as tc
from tacoma.data_io import mkdirp_customdir

standard_config = {
    "temporal_network_files": [
    ],
    "edges_coordinate_files": [
    ],
    "plot_width": 320,
    "titles": [
    ],
    "network_plot_height": 250,
    "edges_plot_height": 100,
    "padding": 10,
    "start_it": 0,
    "node_radius": 2.5,
    "link_distance": 10,
    "node_charge": -8,
    "edge_line_width": 1,
    "font_size_in_px": 14,
    "link_width": 1,
    "d3_format_string": ".3f",
}

html_source_path = os.path.join(tc.__path__[0], 'interactive')


def _make_and_get_directory(path):
    """Simulate ``mkdir -p`` and return the path of the repository"""
    directory, _ = os.path.split(
        os.path.abspath(os.path.expanduser(path))
    )
    mkdirp_customdir(directory)
    return directory


def download_d3():
    """Download `d3.v4.min.js` and save it in `~/.tacoma/d3`, if the file does not exist yet."""

    url = "https://cdnjs.cloudflare.com/ajax/libs/d3/4.13.0/d3.min.js"
    filename = "~/.tacoma/web/d3.v4/d3.v4.min.js"
    target_name = "d3.min.js"
    if not os.path.exists(os.path.abspath(os.path.expanduser(filename))):

        # get directory name for download
        directory = _make_and_get_directory(filename)

        # download
        print("downloading d3 ...")
        wget.download(url, out=directory)

        # move to better name
        os.rename(os.path.join(directory, target_name),
                  os.path.abspath(os.path.expanduser(filename)))


def prepare_visualization_directory():
    """Move all files from the tacoma/interactive directory to ~/.tacoma/web"""
    src = html_source_path
    dst = os.path.abspath(os.path.expanduser("~/.tacoma/web/"))
    if not os.path.exists(os.path.join(dst, "index.html")):
        copy_tree(src, dst)


class StoppableHTTPServer(http.server.HTTPServer):
    """Taken from https://stackoverflow.com/questions/268629/how-to-stop-basehttpserver-serve-forever-in-a-basehttprequesthandler-subclass """

    def __init__(self, server_address, handler, subfolder):
        http.server.HTTPServer.__init__(self, server_address, handler)
        self.subfolder = subfolder

        while subfolder.endswith('/'):
            subfolder = subfolder[:-1]

        self.subjson = subfolder + '_config.json'

    def run(self):
        try:
            self.serve_forever()
        except OSError:
            pass

    def stop_this(self):
        # Clean-up server (close socket, etc.)
        print('was asked to stop the server')
        self.server_close()

        # try:
        if os.path.exists(self.subjson):
            os.remove(self.subjson)
        # except
        if os.path.exists(self.subfolder):
            try:
                # os.rmdir(self.subfolder)
                shutil.rmtree(self.subfolder)
            except FileNotFoundError as e:
                raise e

        # os.chdir(self.cwd)
        print('deleted all files')

    # def __del__(self):
    #    self.stop_this()


def _get_prepared_network(tn, dt, time_unit, time_normalization_factor):
    """Prepare the provided network (i.e. bin it in discrete time)"""

    tn_b = tc.bin(tn, dt=dt)  # rebin the network
    # rescale the network's time
    tn_b.t = [t * time_normalization_factor for t in tn_b.t]
    tn_b.tmax *= time_normalization_factor
    if time_unit is None:
        time_unit = ""
    tn_b.time_unit = time_unit  # set time unit

    return tn_b


def visualize(temporal_networks,
              frame_dt,
              time_normalization_factor=1,
              new_time_unit=None,
              titles=None,
              config=None,
              port=8226):
    """
    Visualize a temporal network or a list of temporal networks interactively.
    This routine starts up an HTTP server, bins the networks according to the
    time step ``frame_dt`` and copies them to ``~/.tacoma/web``. Subsequently,
    a the interaction is started in the standard browser.

    The visualization is stopped with KeyboardInterrupt. The temporary
    temporal network files will subsequently be deleted.

    Parameters
    ----------
    temporal_networks : an instance of :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_lists_with_histograms`, :class:`_tacoma.edge_changes_with_histograms` or a list containing those.
        The temporal networks to visualize. If a list is provided, all networks need to have the
        same `t0` and `tmax`.
    frame_dt : float
        The duration of a frame in the visualization.
    time_normalization_factor : float, default : 1.0
        Rescale time with this factor.
    time_unit : string, default : None,
        Unit of time of the visualization.
    titles : string or list of strings, default : None
        Titles to put on the figures of the corresponding temporal networks.
    config : dict
        Configuration values for the JavaScript visualization.
    port : int, default : 8226
        Port of the started HTTP server.

    Notes
    -----

    The configuration dictionary is filled with values to control
    the appearance of the visualizations. The standard configuration is

    .. code:: python

        config = {
            "plot_width" : 320 ,
            "network_plot_height" : 250,
            "edges_plot_height" : 100,
            "padding" : 10,
            "start_it" : 0,
            "node_radius" : 2.5,
            "link_distance" : 10,
            "node_charge": -8,
            "edge_line_width" : 1,
            "font_size_in_px" : 14,
            "link_width" : 1,
            "d3_format_string": ".3f",
        }
    """

    if not hasattr(temporal_networks, '__len__'):
        temporal_networks = [temporal_networks]

    if titles is None:
        titles = ["" for _ in temporal_networks]
    elif type(titles) == str or not hasattr(titles, '__len__'):
        titles = [titles]

    # print(titles)

    # define the server address
    # server_address = ('127.0.0.1', port)

    path = "~/.tacoma/web/"
    web_dir = os.path.abspath(os.path.expanduser(path))

    # download d3 if that did not happen yet
    download_d3()

    # copy the html and js files for the visualizations
    prepare_visualization_directory()

    # create a subfolder based on the current time
    subfolder = "tmp_{:x}".format(int(time.time()*1000))
    mkdirp_customdir(directory=web_dir)
    subfolder_path = os.path.join(web_dir, subfolder)
    mkdirp_customdir(directory=subfolder_path)

    # change directory to this directory
    print("changing directory to", web_dir)
    print("starting server here ...", web_dir)
    cwd = os.getcwd()
    os.chdir(web_dir)

    server = StoppableHTTPServer(("127.0.0.1", 8226),
                                 http.server.SimpleHTTPRequestHandler,
                                 subfolder_path,
                                 )

    for itn, tn in enumerate(temporal_networks):
        print("preparing network", titles[itn])
        tn_b = _get_prepared_network(
            tn, frame_dt, new_time_unit, time_normalization_factor)
        taco_fname = os.path.join(subfolder, subfolder+'_'+str(itn)+'.taco')
        edge_fname = os.path.join(subfolder, subfolder+'_'+str(itn)+'.json')
        tc.write_edge_trajectory_coordinates(tn_b,
                                             os.path.join(web_dir, edge_fname),
                                             filter_for_duration=frame_dt * time_normalization_factor)
        tc.write_json_taco(tn_b, os.path.join(web_dir, taco_fname))

        standard_config['temporal_network_files'].append(taco_fname)
        standard_config['edges_coordinate_files'].append(edge_fname)
        standard_config['titles'].append(titles[itn])

    if type(config) == dict:
        standard_config.update(config)

    with open(os.path.join(web_dir, subfolder+'_config.json'), 'w') as f:
        json.dump(standard_config, f)

    # ========= start server ============
    thread = threading.Thread(None, server.run)
    thread.start()

    webbrowser.open("http://localhost:8226/?data=" + subfolder)

    try:
        while True:
            time.sleep(2)
    except KeyboardInterrupt:
        # thread.join()
        print('stopping server ...')
        server.stop_this()
        # thread.join()

    # time.sleep(1)

    print('changing directory back to', cwd)

    os.chdir(cwd)


if __name__ == "__main__":
    # download_d3()
    a = tc.load_json_taco("~/.tacoma/ht09.taco")
    visualize(a, frame_dt=20, titles='HT09', new_time_unit='h',
              time_normalization_factor=1./3600.)
