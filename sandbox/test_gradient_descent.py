import tacoma as tc


orig = tc.convert( tc.load_json_taco('~/.tacoma/ht09.taco') )
#orig = tc.convert( tc.load_json_taco('~/.tacoma/hs13.taco') )
#orig = tc.load_json_taco('~/.tacoma/dtu_1_weeks.taco') 
dt_for_inference = 120.
dt_binning = 20.

fetch_k_scaling = True

if fetch_k_scaling:
    k_scaling = tc.estimate_k_scaling_gradient_descent( 
                                            orig, 
                                            dt_for_inference = dt_for_inference, 
                                            dt_for_binning = dt_binning, 
                                            measurements_per_configuration = 20,
                                            learning_rate = 0.5,
                                            relative_error = 1e-2,
                                         )
else:
    k_scaling = 5

from tacoma.model_conversions import estimate_flockwork_P_args
import matplotlib.pyplot as pl
import numpy as np

t_orig, k_orig = tc.mean_degree(tc.bin(orig,dt_binning))

fig, ax = pl.subplots(1,2,sharex=True,sharey=True)
ax[0].plot(t_orig, k_orig,label='original')
ax[1].plot(t_orig, k_orig,label='original')

n_k_meas = 6

for iscl, scaling in enumerate([k_scaling, 1.0]):

    these_ks = []

    for meas in range(n_k_meas):


        args = estimate_flockwork_P_args(
                orig,
                dt = dt_for_inference,
                k_over_k_real_scaling = scaling,
                ensure_empty_network = True,
                adjust_last_bin_if_dt_does_not_fit = True
                )

        fwP = tc.flockwork_P_varying_rates(**args)

        fwP_binned = tc.bin(fwP,dt_binning)

        t_fw, k_fw = tc.mean_degree(fwP_binned)
        these_ks.append(k_fw)
        print("scaling =", scaling, " ; measurement =", meas)
        ndx = np.where(np.logical_and(k_orig>0.0, k_fw>0.0))
        print("k_fw / k_orig = ",np.mean(k_fw[ndx] / k_orig[ndx]))

    these_ks = np.array(these_ks).mean(axis=0)


    ax[iscl].plot(t_fw, these_ks,'-',c=tc.color_sequence[iscl+1], alpha=0.8,label='scaling = {0:4.3f}'.format(scaling))
#pl.yscale('log')

pl.show()

