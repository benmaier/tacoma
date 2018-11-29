import tacoma as tc
import matplotlib.pyplot as pl
import numpy as np

import _tacoma as _tc

N = 100
k = 4

AM = _tc.ActivityModel(N,
                       k/(N-1),
                       2/N/(N-1.)
                      )
SIS = _tc.SIS(N,100,1.0,1.0,number_of_initially_infected=100)

_tc.gillespie_SIS_ActivityModel(AM,SIS)

t = np.array(SIS.time)
I = np.array(SIS.I)

pl.plot(t,I,'--')

pl.show()

