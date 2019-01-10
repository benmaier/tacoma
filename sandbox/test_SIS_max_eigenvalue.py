import tacoma as tc
from tacoma.epidemics import get_SIS_max_eigenvalue, get_SIS_critical_infection_rate, get_SIS_critical_recovery_rate

L = tc.edge_lists()

L.N = 4
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

from scipy.sparse import csc_matrix

A = tc.sparse_adjacency_matrices(L,sparse_generator=csc_matrix)


print(get_SIS_max_eigenvalue(A, 1.03,1.0))

print(get_SIS_critical_infection_rate(A, 1.0))
print(get_SIS_critical_recovery_rate(A, 1.0))
