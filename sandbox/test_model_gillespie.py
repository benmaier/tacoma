import sys 

import tacoma as tc
import matplotlib.pyplot as pl
import numpy as np

import _tacoma as _tc

N = 100

t_run_total = 10000/N
recovery_rate = 1.0

seed = 12323423

R0s = [0.5,1.0,1.2,1.5,2,4,6]

for k in range(1,8,3):
    pl.figure()
    pl.title('$k={0:d}$'.format(k))
    for omega in np.logspace(-2,0,4,base=N):
        pl.figure()
        curve_1 = []
        curve_2 = []
        for R0 in R0s:

            print(k, omega, R0)
            AM = _tc.ActivityModel(N,
                                   k/(N-1.),
                                   omega,
                                   verbose = False,
                                  )
            infection_rate = R0 / k * recovery_rate


            SIS = _tc.SIS(N,t_run_total,infection_rate,recovery_rate,
                    number_of_initially_infected=N,
                    verbose=False,
                    seed=seed,
                    sampling_dt=1,
                    )

            _tc.gillespie_SIS_on_ActivityModel(AM,SIS,verbose=False)            

            t = np.array(SIS.time)
            I = np.array(SIS.I)

            this_pl, = pl.plot(t,I,'s',ms=2,alpha=0.5,mfc='None')
            mean_I_1 = tc.time_average(t, I, tmax=t_run_total)



            SIS = _tc.SIS(N,
                          t_run_total,
                          infection_rate,
                          recovery_rate,
                          number_of_initially_infected=N,
                          sampling_dt=1)

            tn = _tc.activity_model(N,k/(N-1.), omega, t_run_total)
            tc.gillespie_SIS(tn, SIS)

            t = np.array(SIS.time)
            I = np.array(SIS.I)

            mean_I_2 = tc.time_average(t, I, tmax=t_run_total)

            curve_1.append(mean_I_1)
            curve_2.append(mean_I_2)
            this_pl, = pl.plot(t,I,'-',lw=1,alpha=0.5,c=this_pl.get_color())

        #this_pl, = pl.plot(R0s,curve_1,'s',ms=2,alpha=0.5,mfc='None')
        #this_pl, = pl.plot(R0s,curve_2,'-',lw=1,alpha=0.5,c=this_pl.get_color())

pl.show()

