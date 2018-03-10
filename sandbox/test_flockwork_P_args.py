import tacoma as tc

ec = tc.convert( tc.dynamic_RGG(3,3,mean_link_duration=1.) )

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


fw_args = tc.get_flockwork_P_args(C,dt=1.)

print fw_args.rewiring_rate
print fw_args.P
print fw_args.__dict__
