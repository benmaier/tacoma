from tacoma.flockwork import flockwork_P

# Define structural properties
N = 100
k = 2.0
P = k / (k+1)

# The disconnection rate per node is set 
# to gamma = 1, so the total run time is 
# in units of 1/gamma = t_run_total
t_run_total = 100

# In order to start with any initial state,
# pass an edge list to this function.
# If you don't, the function assumes
# you want to start with an equilibrium
# configuration and generates one.
fw = flockwork_P(N, P, t_run_total)

# verify that the Flockwork is in a proper
# format
import tacoma as tc
print("Number of rule violations in generated temporal network:", 
      tc.verify(fw))
print(type(fw))
