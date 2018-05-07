import tacoma as tc

a = tc.load_json_taco("~/.tacoma/ht09.taco")
tc.write_edge_trajectory_coordinates(a,"~/Sites/tacoma-interactive/ht09_edge_trajectories.json",filter_for_duration = 0)

a = tc.load_json_taco("~/.tacoma/dtu_1_weeks.taco")
tc.write_edge_trajectory_coordinates(a,"~/Sites/tacoma-interactive/dtu_1_weeks_edge_trajectories.json",filter_for_duration = 1200.)
