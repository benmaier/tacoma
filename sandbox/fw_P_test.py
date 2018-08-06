import tacoma as tc
ec = tc.flockwork_P_varying_rates([],10,[0.4,0.8], 600,[ (0, 1.), (300,2.) ], 600,seed=25456)

print(ec.edges_out[:5])
print(ec.edges_in[:5])

import numpy as np

ec = tc.flockwork_P_varying_rates_for_each_node([],
                                                10,
                                                [np.random.random(10).tolist(),
                                                 np.random.random(10).tolist()], 
                                                600,
                                                [ (0, np.random.random(10).tolist()),
                                                  (300, np.random.random(10).tolist()) ], 
                                                600,
                                                seed=25456)

print(ec.edges_out[:5])
print(ec.edges_in[:5])
