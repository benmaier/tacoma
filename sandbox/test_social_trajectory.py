import tacoma as tc

sample_aggregates = True
N_time_steps = 2

print "======================== BINNING =============================="

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


traj = tc.binned_social_trajectory_from_edge_lists(L,2,N_time_steps=N_time_steps,verbose=True)

print traj

traj = tc.social_trajectory_from_edge_lists(L,2,verbose=True)
for entry in traj:
    print 
    print "    hash =", entry.hash
    print "    size =", entry.size
    print "    time =", entry.time_pairs



print "===== edge_changes => edge_lists ====="

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

traj = tc.binned_social_trajectory_from_edge_changes(C,2,N_time_steps=N_time_steps,verbose=True)
print traj

traj = tc.social_trajectory_from_edge_changes(C,2,verbose=True)
for entry in traj:
    print 
    print "    hash =", entry.hash
    print "    size =", entry.size
    print "    time =", entry.time_pairs
