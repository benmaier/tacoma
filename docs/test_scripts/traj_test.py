import tacoma as tc

traj = tc.edge_trajectories()
entry = tc.edge_trajectory_entry

traj.N = 8
traj.t0 = 0.0
traj.tmax = 8.1
traj.trajectories = [
                        entry( (0,1), [(0., 3.), (7.0, 8.1)] ),
                        entry( (2,5), [(3., 7.31)] ),
                        entry( (1,7), [(0., 1.), (1.5, 4.0)] ),
                    ]
traj.time_unit = 's'
traj.notes = 'This experiment was conducted as a test.'
traj.int_to_node = {
                0 : 'Alice',
                1 : 'Bob',
                2 : 'Clara',
                4 : 'Darren',
                5 : 'Elle',
                5 : 'Felicitas',
                6 : 'George',
                7 : 'Harriett',
              }

C = tc.convert_edge_trajectories(traj)

import pprint
pp = pprint.PrettyPrinter(indent=4)

pp.pprint(C.N)
pp.pprint(C.t0)
pp.pprint(C.t)
pp.pprint(C.tmax)
pp.pprint(C.edges_initial)
pp.pprint(C.edges_in)
pp.pprint(C.edges_out)

