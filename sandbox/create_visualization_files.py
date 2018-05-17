import tacoma as tc
from tacoma.model_conversions import estimate_flockwork_P_args

def get_prepared_network(tn,dt):
    
    tn_b = tc.bin(tn,dt=dt) # rebin the network
    tn_b.t = [ t / 3600. for t in tn_b.t ] # rescale the network's time
    tn_b.tmax /= 3600.
    tn_b.time_unit = 'h' # set time unit

    return tn_b


# ============ HT 09 ==============

tn = tc.load_json_taco("~/.tacoma/ht09.taco")
tn_b = get_prepared_network(tn,dt=20)

tc.write_edge_trajectory_coordinates(tn_b,
                                     "~/Sites/tacoma/data/ht09_edge_trajectories.json",
                                     filter_for_duration = 0)
tc.write_json_taco(tn_b,"~/Sites/tacoma/data/ht09_binned.taco")


aggregated_network = tc.measure_group_sizes_and_durations(tn).aggregated_network
fw_params = estimate_flockwork_P_args(tn,
                                     dt=120,
                                     k_over_k_real_scaling=2.05,
                                     aggregated_network=aggregated_network,
                                     ensure_empty_network=True,
                                     adjust_last_bin_if_dt_does_not_fit=True)

fw = tc.flockwork_P_varying_rates_neighbor_affinity(**fw_params)
fw_b = get_prepared_network(fw,20)
tc.write_edge_trajectory_coordinates(fw_b,
                                     "~/Sites/tacoma/data/fw_ht09_edge_trajectories.json",
                                     filter_for_duration = 0)
tc.write_json_taco(fw_b,"~/Sites/tacoma/data/fw_ht09_binned.taco")

# ============ DTU ==============

tn = tc.load_json_taco("~/.tacoma/dtu_1_weeks.taco")
tn_b = get_prepared_network(tn,dt=300)

tc.write_edge_trajectory_coordinates(tn_b,
                                     "~/Sites/tacoma/data/dtu_edge_trajectories.json",
                                     filter_for_duration = 1.)
tc.write_json_taco(tn_b,"~/Sites/tacoma/data/dtu_binned.taco")


aggregated_network = tc.measure_group_sizes_and_durations(tn).aggregated_network
fw_params = estimate_flockwork_P_args(tn,
                                     dt=3000,
                                     k_over_k_real_scaling=1.3334181228501472,
                                     aggregated_network=aggregated_network,
                                     ensure_empty_network=True,
                                     adjust_last_bin_if_dt_does_not_fit=True)

fw = tc.flockwork_P_varying_rates_neighbor_affinity(**fw_params)
fw_b = get_prepared_network(fw,300)
tc.write_edge_trajectory_coordinates(fw_b,
                                     "~/Sites/tacoma/data/fw_dtu_edge_trajectories.json",
                                     filter_for_duration = 1.)
tc.write_json_taco(fw_b,"~/Sites/tacoma/data/fw_dtu_binned.taco")

# ============ HS13 ==============

tn = tc.load_json_taco("~/.tacoma/hs13.taco")
tn_b = get_prepared_network(tn,dt=20)

tc.write_edge_trajectory_coordinates(tn_b,
                                     "~/Sites/tacoma/data/hs13_edge_trajectories.json",
                                     filter_for_duration = 40/3600.)
tc.write_json_taco(tn_b,"~/Sites/tacoma/data/hs13_binned.taco")


aggregated_network = tc.measure_group_sizes_and_durations(tn).aggregated_network
fw_params = estimate_flockwork_P_args(tn,
                                     dt=1200,
                                     k_over_k_real_scaling=1.839,
                                     aggregated_network=aggregated_network,
                                     ensure_empty_network=True,
                                     adjust_last_bin_if_dt_does_not_fit=True)

fw = tc.flockwork_P_varying_rates_neighbor_affinity(**fw_params)
fw_b = get_prepared_network(fw,20)
tc.write_edge_trajectory_coordinates(fw_b,
                                     "~/Sites/tacoma/data/fw_hs13_edge_trajectories.json",
                                     filter_for_duration = 40/3600.)
tc.write_json_taco(fw_b,"~/Sites/tacoma/data/fw_hs13_binned.taco")
"""
a = tc.load_json_taco("~/.tacoma/dtu_1_weeks.taco")
tc.write_edge_trajectory_coordinates(a,"~/Sites/tacoma/data/dtu_1_weeks_edge_trajectories.json",filter_for_duration = 1200.)

a = tc.load_json_taco("~/.tacoma/hs13.taco")
"""
