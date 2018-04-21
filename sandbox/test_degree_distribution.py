import tacoma as tc

print("===== edge_lists =====")

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

print("errors in temporal network:", tc.verify(L))

P_k = tc.degree_distribution_from_edge_lists(L)
for k, P in enumerate(P_k):
    print(k, P)

P_k = tc.degree_distribution(L)
for k, P in enumerate(P_k):
    print(k, P)

print("===== edge_changes => edge_changes =====")

C = tc.edge_changes()

C.N = 3
C.edges_initial = [ (0,1) ]
C.t0 = 0.0
C.tmax = 3.0
C.t = [ 1.0, 2.0, ]
C.edges_in = [
                [
                    (1,2), (0,2)
                ],
                [
                    (0,1)
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

P_k = tc.degree_distribution_from_edge_changes(C)
for k, P in enumerate(P_k):
    print(k, P)

P_k = tc.degree_distribution(C)
for k, P in enumerate(P_k):
    print(k, P)

