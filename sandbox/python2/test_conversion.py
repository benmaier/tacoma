import tacoma as tc

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


new = tc.convert(L)

print "N =", new.N
print "t0 =", new.t0
print "t =", new.t
print "tmax =", new.tmax
print "edges_in =", new.edges_in
print "edges_out =", new.edges_out
print "edges_initial =", new.edges_initial

new = tc.convert(new)

print 
print "N =", new.N
print "t =", new.t
print "tmax =", new.tmax
print "edges =", new.edges

print "=========== edge changes =============="
C = tc.edge_changes()

C.N = 3
C.edges_initial = [ (0,1) ]
C.t0 = 0.0
C.tmax = 3.0
C.t = [ 1.0, 2.0 ]
C.edges_in = [
                [
                    (1,2), (0,2)
                ],
                [
                    (0,1),
                ],
             ]
C.edges_out = [
                [
                    (0,1)
                ],
                [
                    (1,2), (0,2)
                ],
              ]

new = tc.convert(C)
print 
print "N =", new.N
print "t =", new.t
print "tmax =", new.tmax
print "edges =", new.edges

print "=========== edge changes =============="
new = tc.convert(new)

print
print "N =", new.N
print "t0 =", new.t0
print "t =", new.t
print "tmax =", new.tmax
print "edges_in =", new.edges_in
print "edges_out =", new.edges_out
print "edges_initial =", new.edges_initial
print new.notes

