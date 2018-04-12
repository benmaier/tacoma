import tacoma as tc
from tacoma.data_io import load_json_taco
from tacoma.data_io import write_json_taco

print "===== edge_lists ====="

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

print "   == old"
print L.N
print L.t
print L.tmax
print L.edges
print L.int_to_node

write_json_taco(L,"test_write.taco")
new = load_json_taco("test_write.taco")

print "   == new"
print new.N
print new.t
print new.tmax
print new.edges

print "===== edge_changes ====="

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

print "   == old"
print C.N
print C.t
print C.t0
print C.tmax
print C.edges_initial
print C.edges_in
print C.edges_out

write_json_taco(C,"test_write.taco")
new = load_json_taco("test_write.taco")

print "   == new"
print new.N
print new.t
print new.t0
print new.tmax
print new.edges_initial
print new.edges_in
print new.edges_out

