import cFlockwork as cF

dyn_RGG = cF.dynamic_RGG(N=3,t_run_total=10,mean_link_duration=2.0)
test_list = cF.edge_lists()
test_list.copy_from(dyn_RGG)

print test_list.t, dyn_RGG.t
print test_list.tmax, dyn_RGG.tmax
print test_list.edges
print dyn_RGG.edges
print test_list.N, dyn_RGG.N

zsbb = cF.ZSBB_model([],3,0.6,0.6,0.6,t_run_total=100)

test_list = cF.edge_changes()
test_list.copy_from(zsbb)

print test_list.t, zsbb.t
print test_list.t0, zsbb.t0
print test_list.tmax, zsbb.tmax
print test_list.edges_in
print zsbb.edges_in
print test_list.edges_out
print zsbb.edges_out
print test_list.N, zsbb.N
