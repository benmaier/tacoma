import matplotlib.pyplot as pl

import tacoma as tc
from tacoma.drawing import draw_edges
from tacoma.drawing import get_edge_order

import numpy as np

#tn = tc.load_sociopatterns_hypertext_2009()
tn = tc.load_json_taco('~/.tacoma/dtu_1_weeks.taco')
tn = tc.convert(tn)

edge_traj = tc.get_edge_trajectories(tn,return_edge_similarities=True)

print(edge_traj.edge_similarities[:10])

edge_order = get_edge_order(edge_traj,threshold=3600.)


print(edge_order)
print(np.all(edge_order == np.sort(edge_order)))

draw_edges(edge_traj.trajectories,edge_order=edge_order)
draw_edges(edge_traj.trajectories)


pl.show()
