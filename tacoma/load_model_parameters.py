import os
import tacoma as tc

path = os.path.join(tc.__path__[0],'resources')

def load_ht09_dyn_RGG_params():
    fn = os.path.join(path,'ht09_dyn_RGG_params.json')
    return tc.load_json_dict(fn)

def load_ht09_ZSBB_params():
    fn = os.path.join(path,'ht09_zsbb_params.json')
    return tc.load_json_dict(fn)

def load_ht09_flockwork_params(scaled=False):

    if scaled:
        fn = os.path.join(path,'ht09_fwP_params_scaled.json')
    else:
        fn = os.path.join(path,'ht09_fwP_params_unscaled.json')

    return tc.load_json_dict(fn)

def load_dtu_dyn_RGG_params():
    fn = os.path.join(path,'dtu_dyn_RGG_params.json')
    return tc.load_json_dict(fn)

def load_dtu_ZSBB_params():
    fn = os.path.join(path,'dtu_zsbb_params.json')
    return tc.load_json_dict(fn)

def load_dtu_flockwork_params(scaled=False):

    if scaled:
        fn = os.path.join(path,'dtu_fwP_params_scaled.json')
    else:
        fn = os.path.join(path,'dtu_fwP_params_unscaled.json')

    return tc.load_json_dict(fn)


if __name__ == "__main__":

    print load_ht09_flockwork_params(scaled=False),
    print load_ht09_flockwork_params(scaled=True),
    print load_dtu_flockwork_params(scaled=False),
    print load_dtu_flockwork_params(scaled=True),

    for f in [
                load_ht09_dyn_RGG_params,
                load_ht09_ZSBB_params,
                load_dtu_dyn_RGG_params,
                load_dtu_ZSBB_params,
            ]:
        print f()
