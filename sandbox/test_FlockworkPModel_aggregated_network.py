import tacoma as tc
from tacoma.flockwork import flockwork_P_equilibrium_configuration as eq_conf

N = 2
P = 0.6
#E = eq_conf(N,P)
k = P/(1-P)
R0 = 5
recovery_rate = 1
gamma = recovery_rate
infection_rate = R0 * recovery_rate / k
I0 = 1
tmax = 200

E = [ (i,j) for i in range(N-1) for j in range(i+1,N) ]
E = []

print("initial edge list", E)

FW = tc.FlockworkPModel(
                      E,  # initial edge list (list of tuple of ints)
                      N, # number of nodes
                      gamma, # rate with which anything happens
                      P,     # probability to reconnect
                      save_temporal_network=True,
                      save_aggregated_network=True,                    
                      verbose=True,
                      )

SIS = tc.SIS(N, #number of nodes
             tmax, # maximum time of the simulation
             infection_rate,
             recovery_rate,
             number_of_initially_infected = I0, # optional, default: 1
             verbose = True,
            )

tc.gillespie_epidemics(FW, SIS, verbose=True)
tsim = SIS.time[-1]

ec = FW.edge_changes
ec.tmax = tsim

print("ec.N)           ,",ec.N)
print("ec.t0)          ,",ec.t0)
print("ec.t)           ,",ec.t)
print("ec.tmax)        ,",ec.tmax)
print("ec.edges_in)    ,",ec.edges_in)
print("ec.edges_out)   ,",ec.edges_out)
print("ec.edges_initial,",ec.edges_initial)

if len(ec.t) > 0:

    aggtc = tc.aggregated_network(ec)
    aggfw = FW.finish_and_get_aggregated_network(tsim)

    print(len(aggtc), len(aggfw))

    for (e_tc, v_tc), (e_fw, v_fw) in zip(sorted(aggtc.items(),key=lambda x: x[1]), 
                                          sorted(aggfw.items(),key=lambda x: x[1])):
        print(e_tc, v_tc, e_fw, v_fw)


