"""
This module provides an API to load parameters for the model functions
in order to produce surrogate networks that have statistics similar to
the original temporal networks SocioPatterns HT09 and one week of DTU
data.
"""

from __future__ import print_function

import os

import tacoma as tc

# get package path
path = os.path.join(tc.__path__[0],'resources')

def load_ht09_dyn_RGG_params():
    """Load the standard parameters to generate SocioPatterns HT09 
    surrogate networks using the dynamic RGG model.

    Returns
    -------
    :obj:`dict`
        The `kwargs` to pass to :func:`tacoma.dynamic_RGG`.
    """

    fn = os.path.join(path,'ht09_dyn_RGG_params.json')
    return tc.load_json_dict(fn)

def load_ht09_ZSBB_params():
    """Load the standard parameters to generate SocioPatterns HT09 
    surrogate networks using the ZSBB model.

    Returns
    -------
    :obj:`dict`
        The `kwargs` to pass to :func:`tacoma.ZSBB_model`.
    """

    fn = os.path.join(path,'ht09_zsbb_params.json')
    return tc.load_json_dict(fn)

def load_ht09_flockwork_params(scaled=False):
    """Load the standard parameters to generate SocioPatterns HT09 
    surrogate networks using the Flockwork-P model with varying rates.

    Parameters
    ----------
    scaled : bool, optional
        If this is `True`, load the rewiring rate `gamma(t)` and proba-
        bility to reconnect `P(t)` rescaled with a corrective factor.
        This factor emerges because in the original network the mean degree
        is overestimated due to binning of edges. This overestimation
        typically leads to an underestimation of `gamma(t)` and an 
        overestimation of `P(t)`.
        
    Returns
    -------
    :obj:`dict`
        The `kwargs` to pass to :func:`tacoma.flockwork_P_varying_rates`.
    """

    if scaled:
        fn = os.path.join(path,'ht09_fwP_params_scaled.json')
    else:
        fn = os.path.join(path,'ht09_fwP_params_unscaled.json')

    return tc.load_json_dict(fn)

def load_ht09_flockwork_unscaled_params():
    """Load the standard parameters to generate SocioPatterns HT09 
    surrogate networks using the Flockwork-P model with varying rates,
    not corrected for overestimation of the mean degree in the original 
    data.

    Returns
    -------
    :obj:`dict`
        The `kwargs` to pass to :func:`tacoma.flockwork_P_varying_rates`.
    """

    return load_ht09_flockwork_params(False)

def load_ht09_flockwork_scaled_params():
    """Load the standard parameters to generate SocioPatterns HT09 
    surrogate networks using the Flockwork-P model with varying rates,
    corrected for overestimation of the mean degree in the original data.

    Returns
    -------
    :obj:`dict`
        The `kwargs` to pass to :func:`tacoma.flockwork_P_varying_rates`.
    """

    return load_ht09_flockwork_params(True)

def load_dtu_dyn_RGG_params():
    """Load the standard parameters to generate surrogate networks 
    for one week of DTU data using the dynamic RGG model.

    Returns
    -------
    :obj:`dict`
        The `kwargs` to pass to :func:`tacoma.dynamic_RGG`.
    """

    fn = os.path.join(path,'dtu_dyn_RGG_params.json')
    return tc.load_json_dict(fn)

def load_dtu_ZSBB_params():
    """Load the standard parameters to generate surrogate networks 
    for one week of DTU data using the ZSBB model.

    Returns
    -------
    :obj:`dict`
        The `kwargs` to pass to :func:`tacoma.ZSBB_model`.
    """

    fn = os.path.join(path,'dtu_zsbb_params.json')
    return tc.load_json_dict(fn)

def load_dtu_flockwork_params(scaled=False):
    """Load the standard parameters to generate surrogate networks 
    for one week of DTU data using the Flockwork-P model with varying rates.

    Parameters
    ----------
    scaled : bool, optional
        If this is `True`, load the rewiring rate `gamma(t)` and proba-
        bility to reconnect `P(t)` rescaled with a corrective factor.
        This factor emerges because in the original network the mean degree
        is overestimated due to binning of edges. This overestimation
        typically leads to an underestimation of `gamma(t)` and an 
        overestimation of `P(t)`.
        
    Returns
    -------
    :obj:`dict`
        The `kwargs` to pass to :func:`tacoma.flockwork_P_varying_rates`.
    """

    if scaled:
        fn = os.path.join(path,'dtu_fwP_params_scaled.json')
    else:
        fn = os.path.join(path,'dtu_fwP_params_unscaled.json')

    return tc.load_json_dict(fn)

def load_dtu_flockwork_unscaled_params():
    """Load the standard parameters to generate surrogate networks 
    for one week of DTU data using the Flockwork-P model with varying rates,
    not corrected for overestimation of the mean degree in the original data.

    Returns
    -------
    :obj:`dict`
        The `kwargs` to pass to :func:`tacoma.flockwork_P_varying_rates`.
    """

    return load_dtu_flockwork_params(False)

def load_dtu_flockwork_scaled_params():
    """Load the standard parameters to generate surrogate networks 
    for one week of DTU data using the Flockwork-P model with varying rates,
    corrected for overestimation of the mean degree in the original data.

    Returns
    -------
    :obj:`dict`
        The `kwargs` to pass to :func:`tacoma.flockwork_P_varying_rates`.
    """

    return load_dtu_flockwork_params(True)


def load_hs13_dyn_RGG_params():
    """Load the standard parameters to generate SocioPatterns HS13
    surrogate networks using the dynamic RGG model.

    Returns
    -------
    :obj:`dict`
        The `kwargs` to pass to :func:`tacoma.dynamic_RGG`.
    """

    fn = os.path.join(path,'ht09_dyn_RGG_params.json')
    return tc.load_json_dict(fn)

def load_hs13_ZSBB_params():
    """Load the standard parameters to generate SocioPatterns HS13 
    surrogate networks using the ZSBB model.

    Returns
    -------
    :obj:`dict`
        The `kwargs` to pass to :func:`tacoma.ZSBB_model`.
    """

    fn = os.path.join(path,'ht09_zsbb_params.json')
    return tc.load_json_dict(fn)

def load_hs13_flockwork_params(scaled=False):
    """Load the standard parameters to generate SocioPatterns HS13
    surrogate networks using the Flockwork-P model with varying rates.

    Parameters
    ----------
    scaled : bool, optional
        If this is `True`, load the rewiring rate `gamma(t)` and proba-
        bility to reconnect `P(t)` rescaled with a corrective factor.
        This factor emerges because in the original network the mean degree
        is overestimated due to binning of edges. This overestimation
        typically leads to an underestimation of `gamma(t)` and an 
        overestimation of `P(t)`.
        
    Returns
    -------
    :obj:`dict`
        The `kwargs` to pass to :func:`tacoma.flockwork_P_varying_rates`.
    """

    if scaled:
        fn = os.path.join(path,'hs13_fwP_params_scaled.json')
    else:
        fn = os.path.join(path,'hs13_fwP_params_unscaled.json')

    return tc.load_json_dict(fn)

def load_hs13_flockwork_unscaled_params():
    """Load the standard parameters to generate SocioPatterns HS13 
    surrogate networks using the Flockwork-P model with varying rates,
    not corrected for overestimation of the mean degree in the original 
    data.

    Returns
    -------
    :obj:`dict`
        The `kwargs` to pass to :func:`tacoma.flockwork_P_varying_rates`.
    """

    return load_hs13_flockwork_params(False)

def load_hs13_flockwork_scaled_params():
    """Load the standard parameters to generate SocioPatterns HS13 
    surrogate networks using the Flockwork-P model with varying rates,
    corrected for overestimation of the mean degree in the original data.

    Returns
    -------
    :obj:`dict`
        The `kwargs` to pass to :func:`tacoma.flockwork_P_varying_rates`.
    """

    return load_hs13_flockwork_params(True)

if __name__ == "__main__":

    print(load_ht09_flockwork_params(scaled=False)),
    print(load_ht09_flockwork_params(scaled=True)),
    print(load_dtu_flockwork_params(scaled=False)),
    print(load_dtu_flockwork_params(scaled=True)),

    for f in [
                load_ht09_dyn_RGG_params,
                load_ht09_ZSBB_params,
                load_dtu_dyn_RGG_params,
                load_dtu_ZSBB_params,
            ]:
        print(f())
