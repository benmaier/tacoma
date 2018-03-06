import tacoma as tc
ec = tc.flockwork_P_varying_rates([],10,[0.4,0.8], 600,[ (0, 1./3600.), (300,2./3600.) ], 600,seed=25456)

print ec.edges_out[:5]
print ec.edges_in[:5]
