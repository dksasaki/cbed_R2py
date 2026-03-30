
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.conversion import localconverter
import os.path as osp
import sys
import xarray as xr
import pandas as pd
from multiprocessing import Pool



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

def empty_ds(i):
    ds_empty = xr.Dataset(
        {name: (["level_0", "level_1"], np.full((20, 1), np.nan))
        for name in names},
        coords={
            "level_0": np.arange(20),
            "level_1": [i],
        }
    )
    return ds_empty

def init_worker():
    global r
    import rpy2.robjects as ro
    r = ro.r
    r(f'source("{ROOT_DIR}/src/cbed_R/cbed_v1_func.R")')
    r(f'source("{ROOT_DIR}/scripts/default_CBED_pars.R")')


def cbed_wrapped(dsmom_a, dscob_a, dscob2_a, i,j,cont):
    dsmom_b = dsmom_a.isel(xh=i, yh=j)
    dscob2_b = dscob2_a.isel(xh=i, yh=j)
    dscob_b = dscob_a.isel(xh=i, yh=j)

    ocean_depth = xr.open_dataset(FTOPO)
    btm_temp = float(dsmom_b['temp'].values)
    btm_O2   = float(dscob_b['btm_o2'].values)
    btm_no3  = float(dscob_b['btm_no3'].values)
    btm_nh4  = float(dscob2_b['nh4'].values)
    btm_dic  = float(dscob_b['btm_dic'].values)
    btm_talk  = float(dscob_b['btm_alk'].values)
    btm_salt = float(dsmom_b['salt'].values)

    fntot    = float(dscob_b['fndet_btm'].values +
                dscob_b['fndi_btm'].values  +
                dscob_b['fnsm_btm'].values  +
                dscob_b['fnmd_btm'].values  +
                dscob_b['fnlg_btm'].values)

    # need to adjust the rates for the appropriate time period
    w_tem = float(
        (dscob_b['fcadet_arag_btm'].values * 100/2.71  +
        dscob_b['fcadet_calc_btm'].values * 100/2.94)  +
        (dscob_b['fsidet_btm'].values +
        dscob_b['fsimd_btm'].values  +
        dscob_b['fsilg_btm'].values) * 60/2.65         +
        (dscob_b['ffedet_btm'].values +
        dscob_b['ffedi_btm'].values  +
        dscob_b['ffesm_btm'].values  +
        dscob_b['ffemd_btm'].values  +
        dscob_b['ffelg_btm'].values) * 160/5.24        +
        (dscob_b['fpdet_btm'].values +
        dscob_b['fpdi_btm'].values   +
        dscob_b['fpsm_btm'].values   +
        dscob_b['fpmd_btm'].values   +
        dscob_b['fplg_btm'].values)  * 120/2.3         +
        dscob_b['flithdet_btm'].values / 2.65          +
        fntot * 6.625 * 22.4 / 0.9
        ) / 1e4 * 60*60*24*365  # cm/yr



    # source = lambda x: f"source(\"{x}\")"
    # r = ro.r
    # r(source(script_path))
    # r(source('scripts/default_CBED_pars.R'))  # contains get_default_pars function

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

    # names were obtained from the cbed code
    names = ["OM1","OM2","OM3","O2","NH4","NO3", "ODU","DIC","TAlk"]
    aux = r2dict(r_out_ss)['y']
    midx = pd.MultiIndex.from_product([np.arange(20),[cont]])
    df = pd.DataFrame(aux, columns=names, index=midx)
    ds = df.to_xarray()
    return ds


def run_point(args):
    cont, i, j, valid_points = args
    if (i, j) in valid_points:
        return cont, cbed_wrapped(dsmom_a, dscob_a, dscob2_a, i, j, cont)
    else:
        return cont, empty_ds(cont)


if __name__ == '__main__':
    ROOT_DIR = '/projects/schultz/d.sasaki/km_scale_model/' + \
                'mom6cobalt_25th/20240723_zstar/tasks/' + \
                '202603_cbed_R2py'
    FPATH     = '/home/d.sasaki/scratch/mom_experiments/cbed_test_001/outputs_raw'
    FTOPO     = '/home/d.sasaki/schultz/d.sasaki/km_scale_model/mom6cobalt_25th/mom_tools/data/grid/nwa25_interped/netcdf3/ocean_topog.nc'
    CACHE_DIR = osp.join(ROOT_DIR, 'data/cache/scratch_test')

    nproc = int(sys.argv[1])

    sys.path.append(osp.join(ROOT_DIR,'scripts/'))  
    import model_reader as mr


    script_path = osp.join(ROOT_DIR,'src/cbed_R/cbed_v1_func.R')

    ds_dict = mr.read_variables(ROOT_DIR,FPATH,FTOPO, CACHE_DIR) 
    print(ds_dict.keys())



    dscob  = ds_dict['dscobalt_btm'].mean(dim='time')
    dscob2 = ds_dict['dscobalt_tr'].mean(dim='time')
    dsmom  = ds_dict['dsmom']

    dscob.load()
    dscob2.load()
    dsmom.load()

    dsmom_a  =dsmom.isel(xh=slice(-4,None), yh=slice(0,4))
    dscob2_a =dscob2.isel(xh=slice(-4,None), yh=slice(0,4))
    dscob_a =dscob.isel(xh=slice(-4,None), yh=slice(0,4))
  
    
    # del(dscob, dscob2, dsmom)
    jvec, ivec = np.where(~np.isnan(dscob_a['btm_o2'].values))
    valid_points = set(zip(ivec, jvec))


    ny, nx = dscob_a['btm_o2'].values.shape
    jm, im = np.meshgrid(np.arange(ny), np.arange(nx), indexing='ij')

    # ds = {}
    # cont = 0
    # for i,j in zip(im.ravel(), jm.ravel()):

    #     if (i,j) in valid_points:
    #         ds[cont] = cbed_wrapped(dsmom_a, dscob_a, dscob2_a,  i, j, cont)
    #     else:
    #         ds[cont] = empty_ds(cont)

    #     cont += 1

    # -- parallel implementation --
    args = [(cont, i, j, valid_points)
            for cont, (i, j) in enumerate(zip(im.ravel(), jm.ravel()))]

    with Pool(processes=nproc, initializer=init_worker) as pool:
        results = pool.map(run_point, args)

    ds = dict(results)

    # -- grouping results --


    test_dict = [ds[i] for i in ds]

    # concat flat list along a new dimension
    combined = xr.concat(test_dict, dim="level_1")

    # reshape each variable into (level_0, level_1, level_2) = (20, 2, 2)
    ds_3d = xr.Dataset(
        {name: (["level_0", "level_1", "level_2"],
                combined[name].values.reshape(20, ny, nx))
        for name in combined.data_vars},
        coords={
            "level_0": combined.level_0.values,
            "level_1": np.arange(ny),
            "level_2": np.arange(nx),
        }
    )

    ds_3d.to_netcdf(osp.join(ROOT_DIR,'data/cache/cbed_mom6.nc'))