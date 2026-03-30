
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.conversion import localconverter
import os.path as osp
import sys
import xarray as xr

def pars2r(pars_dict):
    """Convert a Python dict to an R named list."""
    r_list = ro.ListVector(
        {k: ro.FloatVector([v]) if isinstance(v, (int, float))
            else ro.FloatVector(v)
            for k, v in pars_dict.items()}
    )
    return r_list



def r2dict(varbs):
    result = {}
    for name, val in zip(varbs.names, varbs):
        # nested ListVector — recurse
        if isinstance(val, ro.vectors.ListVector):
            result[name] = r2dict(val)
        else:
            try:
                arr = np.array(val, dtype=float)
                result[name] = arr[0] if arr.size == 1 else arr
            except (ValueError, TypeError):
                result[name] = val
    return result 

if __name__ == '__main__':
    ROOT_DIR = '/projects/schultz/d.sasaki/km_scale_model/' + \
                'mom6cobalt_25th/20240723_zstar/tasks/' + \
                '202603_cbed_R2py'
    FPATH     = '/home/d.sasaki/scratch/mom_experiments/cbed_test_001/outputs_raw'
    FTOPO     = '/home/d.sasaki/schultz/d.sasaki/km_scale_model/mom6cobalt_25th/mom_tools/data/grid/nwa25_interped/netcdf3/ocean_topog.nc'
    CACHE_DIR = osp.join(ROOT_DIR, 'data/cache/scratch_test')


    sys.path.append(osp.join(ROOT_DIR,'scripts/'))  
    import model_reader as mr


    script_path = osp.join(ROOT_DIR,'src/cbed_R/cbed_v1_func.R')

    ds_dict = mr.read_variables(ROOT_DIR,FPATH,FTOPO, CACHE_DIR) 
    print(ds_dict.keys())



    dscob  = ds_dict['dscobalt_btm'].sel(xh=-67, yh=42, method='nearest').mean(dim='time')
    dscob2 = ds_dict['dscobalt_tr'].sel(xh=-67, yh=42, method='nearest').mean(dim='time')
    dsmom  = ds_dict['dsmom'].sel(xh=-67, yh=42, method='nearest')

    dscob.load()
    dscob2.load()
    dsmom.load()

    ocean_depth = xr.open_dataset(FTOPO)
    btm_temp = float(dsmom['temp'].values)
    btm_O2   = float(dscob['btm_o2'].values)
    btm_no3  = float(dscob['btm_no3'].values)
    btm_nh4  = float(dscob2['nh4'].values)
    btm_dic  = float(dscob['btm_dic'].values)
    btm_talk  = float(dscob['btm_alk'].values)
    btm_salt = float(dsmom['salt'].values)

    fntot    = float(dscob['fndet_btm'].values +
                dscob['fndi_btm'].values  +
                dscob['fnsm_btm'].values  +
                dscob['fnmd_btm'].values  +
                dscob['fnlg_btm'].values)

    # need to adjust the rates for the appropriate time period
    w_tem = float(
        (dscob['fcadet_arag_btm'].values * 100/2.71  +
        dscob['fcadet_calc_btm'].values * 100/2.94)  +
        (dscob['fsidet_btm'].values +
        dscob['fsimd_btm'].values  +
        dscob['fsilg_btm'].values) * 60/2.65         +
        (dscob['ffedet_btm'].values +
        dscob['ffedi_btm'].values  +
        dscob['ffesm_btm'].values  +
        dscob['ffemd_btm'].values  +
        dscob['ffelg_btm'].values) * 160/5.24        +
        (dscob['fpdet_btm'].values +
        dscob['fpdi_btm'].values   +
        dscob['fpsm_btm'].values   +
        dscob['fpmd_btm'].values   +
        dscob['fplg_btm'].values)  * 120/2.3         +
        dscob['flithdet_btm'].values / 2.65          +
        fntot * 6.625 * 22.4 / 0.9
        ) / 1e4 * 60*60*24*365  # cm/yr



    source = lambda x: f"source(\"{x}\")"
    r = ro.r
    r(source(script_path))
    r(source('scripts/default_CBED_pars.R'))  # contains get_default_pars function

    # set up parameters
    r_pars = r('get_default_pars()')
    _default_pars = r2dict(r_pars)
    # edit _default_pars as needed

    _default_pars["J.OM"]    = max([0, fntot])    * 6.625*1000*3600*24 *365/10
    _default_pars["O2.w"]    = max([0, btm_O2])   * 1e3
    _default_pars["w"]       = max([0, 0.0359])
    _default_pars["temp"]    = max([0, btm_temp])
    _default_pars["por.0"]   = max([0, 0.8111])
    _default_pars["por.inf"] = max([0, 0.8111])
    _default_pars["NO3.w"]   = max([0, btm_no3])  * 1e3
    _default_pars["NH4.w"]   = max([0, btm_nh4])  * 1e3
    _default_pars["DIC"]     = max([0, btm_dic])
    _default_pars["TAlk.w"]  = max([0, btm_talk])
    _default_pars["S"]       = max([0, btm_salt])
    _default_pars["depth"]   = max([0, 200])



    # _default_pars
    r_pars_updated = pars2r(_default_pars)

    r("""
    manual_control_bioturbation <- F
    manual_control_bioirrigation <- F
    manual_OM_decay_rate <- F
    """)

    ro.globalenv["py_pars"] = r_pars_updated
    r_out_ss = r("cbed_model(py_pars)")  # running the model with py_pars



    test = r2dict(r_out_ss)


    for i in test:
        print(i, test[i].shape)