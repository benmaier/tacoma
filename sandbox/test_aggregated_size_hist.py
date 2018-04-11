import tacoma as tc
import numpy as np

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

result = tc.measure_group_sizes_and_durations_for_edge_lists(L)

print("N_m =", result.aggregated_size_histogram)
print("sum(m*N_m) =", np.dot(np.arange(0,L.N+1),result.aggregated_size_histogram), "should be N =", L.N)


print("===== edge_changes => edge_lists =====")

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

result = tc.measure_group_sizes_and_durations_for_edge_changes(C)

print("N_m =", result.aggregated_size_histogram)
print("sum(m*N_m) =", np.dot(np.arange(0,C.N+1),result.aggregated_size_histogram), "should be N =", C.N)
