import os
import tacoma as tc

from tacoma.drawing import edge_activity_plot
import matplotlib.pyplot as pl

hs13 = tc.load_sociopatterns_high_school_2013()
fig, ax = edge_activity_plot(hs13,
                             time_normalization_factor = 1/3600.,
                             time_unit='h',
                             alpha = 1.0,  # opacity
                             linewidth = 1.5,
                             )

fig.tight_layout()
fig.savefig('./hs13_edge_activity.png',dpi=300)

pl.show()

