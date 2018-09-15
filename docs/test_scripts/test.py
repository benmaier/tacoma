import tacoma as tc

temporal_network = tc.edge_lists()

temporal_network.N = 8
temporal_network.t = [ 0., 1., 1.5, 3., 4., 7., 7.31 ]
temporal_network.tmax = 8.1
temporal_network.edges = [
                            [ (0, 1), (1, 7) ],
                            [ (0, 1) ],
                            [ (0, 1), (1, 7) ],
                            [ (2, 5), (1, 7) ],
                            [ (2, 5) ],
                            [ (0, 1), (2, 5) ],
                            [ (0, 1) ]
                        ]
temporal_network.time_unit = 's'
temporal_network.notes = 'This experiment was conducted as a test.'
temporal_network.int_to_node = {
                0 : 'Alice',
                1 : 'Bob',
                2 : 'Clara',
                4 : 'Darren',
                5 : 'Elle',
                5 : 'Felicitas',
                6 : 'George',
                7 : 'Harriett',
              }

C = tc.convert(temporal_network)

import pprint
pp = pprint.PrettyPrinter(indent=4)

pp.pprint(C.N)
pp.pprint(C.t0)
pp.pprint(C.t)
pp.pprint(C.tmax)
pp.pprint(C.edges_initial)
pp.pprint(C.edges_in)
pp.pprint(C.edges_out)

