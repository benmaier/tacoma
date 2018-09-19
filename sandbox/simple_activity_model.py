import tacoma as tc

# Structural parameters
N = 100
k = 3
rho = k/(N-1.0)

# Temporal parameters
tau = 4
omega = 1.0/tau
t_run_total = 10

# Simulate
temporal_network = tc.activity_model(N, rho, omega, t_run_total)

# Draw
from tacoma.drawing import edge_activity_plot
import matplotlib.pyplot as pl

edge_activity_plot(temporal_network)