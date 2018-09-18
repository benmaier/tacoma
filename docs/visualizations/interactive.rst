Interactive
===========

`tacoma` offers the possibility to interactively visualize temporal networks.
The visualization process works comparable to a video, where the temporal network
is binned to a certain duration ``dt``, called a frame. Then the visualization ticks
through the frames.

The visualization can be started with the function :func:`tacoma.interactive.visualize`.

For a visualization, a new directory is created in `tacoma`'s internal directory, ``~/.tacoma/web``,
and filled with temporary tacos. The program then starts a local HTTP-server
and boots the visualization. In order to stop the server, use the ``KeyboardInterrupt`` signal,
i.e. in the console use ``Ctr+C`` and in Jupyter notebooks use the stop button.

Here's an example on how it works. I think the GUI is pretty self-explanatory.

.. code:: python

   >>> import tacoma as tc
   >>> from tacoma.interactive import visualize
   >>> temporal_network = tc.download_and_convert_sociopatterns_hypertext_2009()
   100% [........................................................] 67463 / 67463
   >>> visualize(temporal_network, frame_dt = 20)

This is the result:

.. figure:: https://github.com/benmaier/tacoma/raw/master/img/ht09_extensive_example.gif
   :alt: visualization example

You can pause the visualization using ``space`` and tick through single frames using the ← and → keys.

The concurrent visualization of multiple temporal networks is possible if you
provide a list of temporal networks and all of them have the same duration.

In order to tailor the visualization to your needs you can pass the argument
`config` to the function, wich is a dictionary of settings for the visualization.
The standard configuration is


.. code:: python

    config = {

        # width of the visualization
        "plot_width" : 320 ,

        # Height of the network plot
        "network_plot_height" : 250,

        # Height of the edge activity plot
        "edges_plot_height" : 100,

        # whitespace between plots
        "padding" : 10,

        # at which frame to start the visualization
        "start_it" : 0,

        # node radius
        "node_radius" : 2.5,

        # mean distance between nodes
        "link_distance" : 10,

        # repelling force between nodes
        "node_charge": -8,

        # width of edges in edge activity plot
        "edge_line_width" : 1,

        # font size of titles and temporal labels
        "font_size_in_px" : 14,

        # width of links in the network plot
        "link_width" : 1,

        # precision of time label
        "d3_format_string": ".3f",
    }
