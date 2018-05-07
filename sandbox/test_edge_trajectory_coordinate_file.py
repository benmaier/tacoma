import tacoma as tc

a = tc.load_json_taco("~/.tacoma/ht09.taco")
tc.write_edge_trajectory_coordinates(a,"ht09_edge_trajectories.json")
