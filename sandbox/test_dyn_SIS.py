import tacoma as tc
import matplotlib.pyplot as pl
import DynGillEpi as gill
import numpy as np

sample_aggregates = False
N_time_steps = 5

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

L_sis = tc.SIS(L.N, 100, 1.0, 0.1,3,0,12234235,True)

tc.gillespie_SIS_on_edge_lists(L,L_sis,verbose=True)

pl.plot(L_sis.time,L_sis.I)
pl.plot(L_sis.time,L_sis.SI)
pl.plot(L_sis.time,L_sis.R0)

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

C_sis = tc.SIS(C.N, 100, 1.0, 0.1,3,0,12234235,True)

tc.gillespie_SIS_on_edge_changes(C,C_sis,verbose=True)

pl.plot(C_sis.time,C_sis.I,'--')
pl.plot(C_sis.time,C_sis.SI,'--')
pl.plot(C_sis.time,C_sis.R0,'--')

print("================ SIS vestergaard ==============")
result = gill.SIS_Poisson_homogeneous(L.N,L.edges,1.0,0.1,100,seed=12234235,initial_number_of_infected=3,verbose=True)

print(result)

#It = np.array(result[0].I_of_t)
#t = It[:,0]
#I = It[:,1]
I = result.I[0]
t = np.arange(len(result.I[0]))

pl.plot(t,I,'--')
#pl.plot(result.true_t,result.true_I,'--')
#pl.plot(result.t,result.SI,'--')

pl.show()

