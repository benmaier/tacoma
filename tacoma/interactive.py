from __future__ import print_function

import os
import sys

"""
vers = sys.version_info[0]
if vers == 3:
    from http import server
"""


import http.server
import webbrowser
import time
import threading

import wget

import tacoma as tc
from tacoma.data_io import mkdirp_customdir

def _make_and_get_directory(path):
    directory, single = os.path.split(os.path.abspath(os.path.expanduser(path)))
    mkdirp_customdir(directory)
    return directory

def download_d3():
    """Download `d3.v4.min.js` and save it in `~/.tacoma/d3`."""

    url = "https://cdnjs.cloudflare.com/ajax/libs/d3/4.13.0/d3.min.js"
    filename = "~/.tacoma/web/d3/d3.v4.min.js"
    target_name = "d3.min.js"
    if not os.path.exists( os.path.abspath(os.path.expanduser(filename)) ):

        # get directory name for download
        directory = _make_and_get_directory(filename)

        # download
        print("downloading d3 ...")
        wget.download(url,out=directory)

        # move to better name
        os.rename(os.path.join(directory, target_name), os.path.abspath(os.path.expanduser(filename)) )

class StoppableHTTPServer(http.server.HTTPServer):
    """ from https://stackoverflow.com/questions/268629/how-to-stop-basehttpserver-serve-forever-in-a-basehttprequesthandler-subclass """

    def __init__(self, server_address, handler, subfolder):
        http.server.HTTPServer.__init__(self, server_address, handler)
        self.subfolder = subfolder

        while subfolder.endswith('/'):
            subfolder = subfolder[:-1]
            
        self.subjson = subfolder + '.json'

        self.thread = threading.Thread(None, self.run)
        self.thread.start()

    def run(self):
        try:
            self.serve_forever()
        except KeyboardInterrupt:
            self.stop_this()
        finally:
            self.stop_this()

    def stop_this(self):
        # Clean-up server (close socket, etc.)
        self.server_close()
        os.remove(self.subjson)
        shutil.rmtree(self.subfolder)
        self.thread.join()

    def __del__(self):
        self.stop_this()

def visualize(temporal_networks,port = 8226):


    if not hasattr(temporal_networks,'__len__'):
        temporal_networks = [ temporal_networks ]

    # define the server address
    server_address = ('127.0.0.1', port)

    path = "~/.tacoma/web/"
    web_dir = os.path.abspath(os.path.expanduser(path))

    # download d3 if that did not happen yet
    download_d3()
    
    # create a subfolder based on the current time
    subfolder = "{:x}".format(int(time.time()*1000))
    mkdirp_customdir(directory = web_dir)
    mkdirp_customdir(directory = os.path.join(web_dir,subfolder))

    # change directory to this directory
    print("starting server at", web_dir)
    cwd = os.getcwd()
    os.chdir(web_dir)

    server = StoppableHTTPServer(("127.0.0.1", 8226),
                                  http.server.SimpleHTTPRequestHandler,
                                  os.path.join(web_dir,subfolder)
                                  )
    # thread = threading.Thread(None, server.run)
    # thread.start()

    webbrowser.open("http://localhost:8226/")

    try:
        while True:
            time.sleep(2)
    except KeyboardInterrupt:
        server.stop_this()

    #time.sleep(1)
    #os.chdir(cwd)

if __name__ == "__main__":
    #download_d3()
    a = tc.load_json_taco("~/.tacoma/ht09.taco")
    visualize(a)


