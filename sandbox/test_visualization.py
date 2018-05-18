import tacoma as tc
from tacoma.interactive import visualize

print("===== edge_lists => edge_lists =====")

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


print("=========== edge changes ==============")
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


print("========= start visualizing ========")

visualize([L, C], frame_dt = 0.01, titles = ["From edge_lists", "From edge_changes"])
