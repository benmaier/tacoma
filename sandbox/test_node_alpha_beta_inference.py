import tacoma as tc

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

dt = 1.0
k_scaling = 1.0

from tacoma.model_conversions import estimate_flockwork_alpha_beta_args_for_single_nodes
single_node_params = estimate_flockwork_alpha_beta_args_for_single_nodes(
                                       C,
                                       dt = dt,
                                       adjust_last_bin_if_dt_does_not_fit = True,
                                       k_over_k_real_scaling = k_scaling,
                                       apply_mean_correction = True,
                                       verbose = True,
                                       )

