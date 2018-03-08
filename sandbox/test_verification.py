import tacoma as tc

L = tc.edge_lists()

L.N = 2
L.t = [0.,0.2,0.1]
L.tmax = 0.11
L.edges = [
            [ (0,1), (0,0) ],
            [ (1,2) ],
        ]

tc.verify_edge_lists(L)
