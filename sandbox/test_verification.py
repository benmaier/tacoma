import tacoma as tc

L = tc.edge_lists()

L.N = 2
L.t = [0.,0.2,0.1]
L.tmax = 0.11
L.edges = [
            [ (0,1), (0,0), (1,0) ],
            [ (1,2) ],
        ]

tc.verify(L)

print 

C = tc.edge_changes()

C.N = 2
C.t0 = 0
C.t = [0.2,0.1]
C.tmax = 0.05
C.edges_initial = [
                    (5,4), (0,1), (0,1), 
                ]

C.edges_in = [ [(2,3), (0,0), (2,3)],[(1,0),(0,1)] ]
C.edges_out = [ [(0,0), (1,0), (1,0)],[(1,0),(0,1) ] ]


tc.verify(C)
