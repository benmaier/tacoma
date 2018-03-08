import tacoma as tc

sample_aggregates = False
N_time_steps = 2

print "===== edge_lists => edge_lists ====="

L = tc.edge_lists()

L.N = 3
L.t = [0.0,1.0,2.0]
L.tmax = 3.0
L.edges = [ 
            [
              (0,1)
            ],
            [
              (1,2), (0,2)
            ],
            [
              (0,1)
            ],
           ]


new = tc.concatenate_edge_lists([L,L,L])

print new.N
print new.t
print new.tmax
print new.edges

print "===== edge_changes => edge_changes ====="

C = tc.edge_changes()

C.N = 3
C.edges_initial = [ (0,1) ]
C.t0 = 1.0
C.tmax = 3.0
C.t = [ 2.0, ]
C.edges_in = [
                [
                    (1,2), (0,2)
                ],
             ]
C.edges_out = [
                [
                    (0,1)
                ],
              ]

new = tc.concatenate_edge_changes([C,C,C])

print new.N
print new.t0
print new.t
print new.tmax
print new.edges_initial
print new.edges_in
print new.edges_out

