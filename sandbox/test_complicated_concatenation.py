import tacoma as tc

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
L.int_to_node = {
            0: 'a',
            1: 'b',
            2: 'c',
        }

L2 = tc.edge_lists()

L2.N = 4
L2.t = [0.0,1.0,2.0]
L2.tmax = 3.0
L2.edges = [ 
            [
              (1,3)
            ],
            [
              (1,3), (2,3)
            ],
           ]
L2.int_to_node = {
            0: 'a',
            1: 'd',
            2: 'b',
            3: 'c',
            }

new = tc.concatenate([L,L2,L])

print(new.N)
print(new.t)
print(new.tmax)
print(new.edges)

print("===== edge_changes => edge_changes =====")

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
C.int_to_node = {
            0: 'a',
            1: 'b',
            2: 'c',
        }

C2 = tc.edge_changes()

C2.N = 4
C2.edges_initial = [ (1,3) ]
C2.t0 = 1.0
C2.tmax = 3.0
C2.t = [ 2.0, ]
C2.edges_in = [
                [
                    (2,3), (0,2)
                ],
             ]
C2.edges_out = [
                [
                    (1,3)
                ],
              ]
C2.int_to_node = {
            0: 'a',
            1: 'd',
            2: 'b',
            3: 'c',
        }

new = tc.concatenate([C,C2,C])

print(new.N)
print(new.t0)
print(new.t)
print(new.tmax)
print(new.edges_initial)
print(new.edges_in)
print(new.edges_out)
print(new.int_to_node)
