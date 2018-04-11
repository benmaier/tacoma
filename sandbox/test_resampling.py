import tacoma as tc

sample_aggregates = False
N_time_steps = 2

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
L.time_unit = 's'

new = tc.bin_from_edge_lists(L,N_time_steps=N_time_steps,verbose=True)

print(new.N)
print(new.t)
print(new.tmax)
print(new.edges)
print(new.notes)
print(new.time_unit)

print("=====binning from edge lists=======")

new = tc.bin_from_edge_lists(L,N_time_steps=N_time_steps,verbose=True)

print(new.N)
print(new.t)
print(new.tmax)
print(new.edges)

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

new = tc.sample_from_edge_changes(C,N_time_steps=N_time_steps,sample_aggregates=sample_aggregates,verbose=True)

print(new.N)
print(new.t)
print(new.tmax)
print(new.edges)

new = tc.bin_from_edge_changes(C,N_time_steps=N_time_steps,verbose=True)

print(new.N)
print(new.t)
print(new.tmax)
print(new.edges)
