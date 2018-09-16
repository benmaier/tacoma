Read, write and share
=====================

Storing and sharing data is obviously rather important. Hence we need a file-format which enables
us to easily read, write and share the data present in :class:`_tacoma.edge_lists` and 
:class:`_tacoma.edge_changes`.

The status quo --- csv
----------------------

Typically, data is shared in csv-format (or tsv)
where the first column refers to a discrete time, 
the second column contains a source node and the third column contains a target node,
such as 

.. csv-table:: A csv-example for edge lists
    :header: "time", "source", "target"

    20, "Alice", "Bob"
    20, "Alice", "Clara"
    60, "Darren", "Alice"

These formats usually describe edge lists at discrete times.

For data sets describing edge changes, this format often includes a fourth column
containing a descriptor for whether the edge was activated or deactivated at
the corresponding time in the first columns, such as

.. csv-table:: A csv-example for edge changes
    :header: "time", "source", "target", "in or out"

    20, "Alice", "Bob", "in"
    20, "Alice", "Clara", "in"
    60, "Alice", "Bob", "out"
    60, "Alice", "Clara", "out"
    60, "Darren", "Alice", "in"

The only advantage of these formats is that they're easily readable by humans.

The disadvantages are

- do not contain the total number of nodes
- no description of the data, information about the experiment
- ambiguous dimension of time
- final time of the experiment :math:`t_\mathrm{max}` is not given
- ambiguous about the class of the provided data

A new format --- The taco
-------------------------

The number of cases in which data is opened to be read by a human is significantly
lower than the number of cases in which it is opened to be read by a machine.
My personal opinion is therefore that data should be stored in a format easily readible
by a computer and converted to appropriate formats whenever a human wants to access 
the data directly.

`tacoma` comes with its own data-format called ".taco" which is simply the whole data 
of either a :class:`_tacoma.edge_lists` or a :class:`_tacoma.edge_changes` object
dumped to a JSON
object and then written to a file. We chose JSON since there exists packages to read
and write JSON in nearly every programming environment. Furthermore, JSON-strings are
considerably more light-weight than, e.g. XML-strings.

Hence, not unlike to a real taco, we dump all the good stuff to a single shell --- 
delicious, but
difficult for human consumption.

Writing
~~~~~~~

In `tacoma`, we can save a temporal network to a `taco` using :func:`tacoma.data_io.write_json_taco`.

.. code::
    
    tc.write_json_taco(temporal_network, '/path/to/temporal_network.taco')

The fields in the written JSON object are the same as they are in the formats defined in
":ref:`temp_network_classes`", besides the additional field ``'type'`` which can either be ``'edge_lists'``
or ``'edge_changes'``.

Considering the example temporal network given in ":ref:`temp_network_classes`"
the JSON-string for :class:`_tacoma.edge_lists` is congruent to

.. code:: javascript

    { 
      'type': 'edge_lists', 
      't': [0.0, 1.0, 1.5, 3.0, 4.0, 7.0, 7.31], 
      'tmax': 8.1,
      'N': 8, 
      'edges': [ [[0, 1], [1, 7]], 
                 [[0, 1]], 
                 [[0, 1], [1, 7]],
                 [[2, 5], [1, 7]],
                 [[2, 5]], 
                 [[0, 1], [2, 5]], 
                 [[0, 1]]
               ],
      'int_to_node': { 
            '0': 'Alice', 
            '1': 'Bob', 
            '2': 'Clara', 
            '3': 'Darren',
            '4': 'Elle',
            '5': 'Felicitas',
            '6': 'George',
            '7': 'Harriett'
            }, 
      'notes': 'This experiment was conducted as a test.', 
      'time_unit': 's'
    }

The actual data is, however, minified and looks more like

.. code:: javascript

    {"type":"edge_lists","t":[0.0,1.0,1.5,3.0,4.0,7.0,7.31],"tmax":8.1,"N":8,"edges":[[[0,1],[1,7]],[[0,1]],[[0,1],[1,7]],[[2,5],[1,7]],[[2,5]],[[0,1],[2,5]],[[0,1]]],"int_to_node":{"0":"Alice","1":"Bob","2":"Clara","3":"Darren","4":"Elle","5":"Felicitas","6":"George","7":"Harriett"},"notes":"This experiment was conducted as a test.","time_unit":"s"}

For :class:`_tacoma.edge_changes` it would look like

.. code:: javascript

    {
      'type': 'edge_changes', 
      't': [1.0, 1.5, 3.0, 4.0, 7.0, 7.31], 
      't0': 0.0, 
      'tmax': 8.1, 
      'N': 8, 
      'edges_initial': [[0, 1], [1, 7]], 
      'edges_in': [[], [[1, 7]], [[2, 5]], [], [[0, 1]], []], 
      'edges_out': [[[1, 7]], [], [[0, 1]], [[1, 7]], [], [[2, 5]]], 
      'int_to_node': {
            '0': 'Alice', 
            '1': 'Bob', 
            '2': 'Clara', 
            '3': 'Darren', 
            '4': 'Elle', 
            '5': 'Felicitas', 
            '6': 'George', 
            '7': 'Harriett'
          }, 
      'notes': 'This experiment was conducted as a test.', 
      'time_unit': 's'
    }

The actual data is, however, minified and looks more like

.. code:: javascript

    {"type":"edge_changes","t":[1.0,1.5,3.0,4.0,7.0,7.31],"t0":0.0,"tmax":8.1,"N":8,"edges_initial":[[0,1],[1,7]],"edges_in":[[],[[1,7]],[[2,5]],[],[[0,1]],[]],"edges_out":[[[1,7]],[],[[0,1]],[[1,7]],[],[[2,5]]],"int_to_node":{"0":"Alice","1":"Bob","2":"Clara","3":"Darren","4":"Elle","5":"Felicitas","6":"George","7":"Harriett"},"notes":"This experiment was conducted as a test.","time_unit":"s"}

Reading
~~~~~~~

Reading temporal network data from a taco is as simple as using :func:`tacoma.data_io.load_json_taco`.

.. code:: python

    temporal_network = tc.load_json_taco('temporal_network.taco')

Converting csv to taco
~~~~~~~~~~~~~~~~~~~~~~

As indicated above, converting csv-data to data actually
usable by algorithms can turn out quite tideous.
Below you can find a commented example on how to load a csv-file
(here from the SocioPatterns `‘Hypertext 2009’-dataset`_) and
convert it to a taco, taken directly from :mod:`tacoma.data_io`

.. code:: python

    import gzip
    import csv

    import tacoma as tc

    # open gzipped file
    gzip_file = 'ht09_contact_list.dat.gz'
    with gzip.open(gzip_file,mode='rt') as f:
        reader = csv.reader(f,delimiter='\t')

        # mappings of nodes to integers
        node_to_int = {}
        int_to_node = {}

        # get an initial t_old
        # (this is done to detect changes in the tsv)
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

                # When the time changed more than dt,
                # append an instance of an empty edge list
                # at t = t_old + dt
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
    with open('ht09.taco','w') as f:
        tc.write_json_taco(el,f)


Edge coordinate files
---------------------

:func:`tacoma.data_io.write_edge_trajectory_coordinates`


Loading existing data
---------------------

:func:`tacoma.data_io.download_and_convert_sociopatterns_high_school_2013`

:func:`tacoma.data_io.download_and_convert_sociopatterns_hypertext_2009`

:func:`tacoma.data_io.load_sociopatterns_high_school_2013`

:func:`tacoma.data_io.load_sociopatterns_hypertext_2009`



.. _‘Hypertext 2009’-dataset: http://www.sociopatterns.org/datasets/hypertext-2009-dynamic-contact-network/
