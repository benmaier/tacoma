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
