
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
import logging
import model_reader as mr
from porosity2cbed import porosity_main
import tomllib
import time


def load_config(config_path='config.toml'):
    with open(config_path, 'rb') as f:
        cfg = tomllib.load(f)

    paths  = cfg['paths']
    chunks = cfg['chunks']

    ROOT_DIR  = paths['root_dir']
    FPATH     = paths['fpath']
    FTOPO     = paths['ftopo']
    FGRD      = paths['fgrd']
    FTEMPLATE = paths['ftemplate']
    FOUT      = paths['fout']
    CACHE_DIR = paths['cache_dir']
    nx_chunk  = chunks['nx']
    ny_chunk  = chunks['ny']

    return ROOT_DIR, FPATH, FTOPO, FGRD, FTEMPLATE, FOUT, CACHE_DIR, nx_chunk, ny_chunk


# defaults of the [cbed] section of config.toml
CBED_DEFAULTS = {
    'method': 'steady_state',   # 'steady_state' (original) | 'newton_krylov'
    'model':  'fast',           # cbed2py model used by 'newton_krylov': 'fast' | 'v2'
    'params': {},               # cbed2py parameter overrides ('newton_krylov')
    'spinup': {},               # keywords of cbed2py.spinup.newton_krylov_spinup
}
METHODS = ('steady_state', 'newton_krylov')


def load_cbed_config(config_path='config.toml'):
    """[cbed] section of config.toml, with defaults (see CBED_DEFAULTS)."""
    with open(config_path, 'rb') as f:
        cfg = tomllib.load(f).get('cbed', {})
    out = {**CBED_DEFAULTS, **cfg}
    if out['method'] not in METHODS:
        raise ValueError(f"[cbed] method must be one of {METHODS}, got {out['method']!r}")
    return out


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

    names=("OM1","OM2","OM3","O2","NH4","NO3", "ODU","DIC","TAlk")

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
    if CBED_CFG['method'] == 'newton_krylov':
        import cbed2py
        cbed2py.set_model(CBED_CFG['model'])
        cbed2py.default_params()        # sources the cbed2py R model once per worker
        return
    import rpy2.robjects as ro
    r = ro.r
    r(f'source("{ROOT_DIR}/src/cbed_R/cbed_v1_func.R")')
    r(f'source("{ROOT_DIR}/scripts/default_CBED_pars.R")')


def cbed_wrapped(dsmom_a, dscob_a, dscob2_a,dsporo_a, i,j,cont, start_r=False):

    if start_r is False:
        global r

    print(f"Running point i={i} j={j} cont={cont}", flush=True)

    dsmom_b = dsmom_a.isel(xh=i, yh=j)
    dscob2_b = dscob2_a.isel(xh=i, yh=j)
    dscob_b = dscob_a.isel(xh=i, yh=j)
    dsporo_b = dsporo_a.isel(xh=i, yh=j)

    ocean_depth = xr.open_dataset(FTOPO)
    btm_temp = float(dsmom_b['temp'].values)
    btm_O2   = float(dscob_b['btm_o2'].values)
    btm_no3  = float(dscob_b['btm_no3'].values)
    btm_nh4  = float(dscob2_b['nh4'].values)
    btm_talk  = float(dscob2_b['talk'].values)

    btm_dic  = float(dscob_b['btm_dic'].values)
    # btm_talk  = float(dscob_b['btm_alk'].values)
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

    porosity = float(dsporo_b['porosity'].values)

    if start_r:
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
    _default_pars["por.0"]   = max([0, porosity])
    _default_pars["por.inf"] = max([0, porosity])
    _default_pars["NO3.w"]   = max([0, btm_no3])  * 1e3
    _default_pars["NH4.w"]   = max([0, btm_nh4])  * 1e3
    _default_pars["DIC"]     = max([0, btm_dic])
    _default_pars["TAlk.w"]  = max([0, btm_talk])
    _default_pars["S"]       = max([0, btm_salt])
    _default_pars["depth"]   = max([0, 200])

    pars_items = ["J.OM", "O2.w", "w", "temp", "por.0",
                  "por.inf", "NO3.w", "NH4.w", "DIC",
                  "TAlk.w", "S", "depth"]
    


    # _default_pars
    r_pars_updated = pars2r(_default_pars)

    r("""
    manual_control_bioturbation <- F
    manual_control_bioirrigation <- F
    manual_OM_decay_rate <- F
    """)

    ro.globalenv["py_pars"] = r_pars_updated
    try:
        r_out_ss = r("cbed_model(py_pars)")
    except Exception as e:
        print(f"FAILED i={i} j={j} cont={cont}", flush=True)
        print(f'{i} {j}',', '.join(f"{k}={_default_pars[k]}" for k in pars_items), flush=True)
        return empty_ds(cont)
    # names were obtained from the cbed code
    names = ["OM1","OM2","OM3","O2","NH4","NO3", "ODU","DIC","TAlk"]
    aux = r2dict(r_out_ss)['y']
    midx = pd.MultiIndex.from_product([np.arange(20),[cont]])
    df = pd.DataFrame(aux, columns=names, index=midx)
    ds = df.to_xarray()
    return ds


# ---------------------------------------------------------------------------
# method 'newton_krylov': cbed2py transient model spun up to its periodic
# (cyclostationary) state under a one-year climatology of the forcing
# ---------------------------------------------------------------------------
STATE_NAMES = ["OM1", "OM2", "OM3", "O2", "NH4", "NO3", "ODU", "DIC", "TAlk"]
J_CONV = 6.625 * 1000 * 3600 * 24 * 365 / 10   # mol N m-2 s-1 -> umol C cm-2 yr-1 (as cbed_wrapped)
DIAG_KEYS = ('status', 'wall_time', 'model_years', 'newton_steps', 'imbalance', 'rms_G')
# status: 1 converged, 0 not converged (best state kept), -1 failed, -2 land


def forcing_fields(dscob_t, dscob2_t):
    """
    Time-varying CBED forcing [CBED units] from the COBALT bottom fields:
      btm: J.OM, O2.w, NO3.w, DIC.w   (cobalt_btm)
      tr : NH4.w, TAlk.w              (cobalt_tracers, bottom level)
    Concentrations mol kg-1 -> mM (x1e3), as done for O2/NO3/NH4 in cbed_wrapped.
    """
    fntot = (dscob_t['fndet_btm'] + dscob_t['fndi_btm'] + dscob_t['fnsm_btm'] +
             dscob_t['fnmd_btm'] + dscob_t['fnlg_btm'])
    btm = xr.Dataset({'J.OM':  fntot * J_CONV,
                      'O2.w':  dscob_t['btm_o2'] * 1e3,
                      'NO3.w': dscob_t['btm_no3'] * 1e3,
                      'DIC.w': dscob_t['btm_dic'] * 1e3})
    tr = xr.Dataset({'NH4.w':  dscob2_t['nh4'] * 1e3,
                     'TAlk.w': dscob2_t['talk'] * 1e3})
    return btm, tr


def annual_climatology(ds):
    """
    Mean seasonal cycle of a time series, on coordinate 'phase' [yr] in (0, 1)
    (phase 0 = Jan 1): by day of year (day 366 dropped) for sub-monthly data,
    by month otherwise. Multi-year input is averaged over the years.
    """
    if ds.sizes['time'] < 2:
        raise ValueError('newton_krylov needs time-varying forcing (time size >= 2)')
    dt_days = np.median(np.diff(ds.time.values)) / np.timedelta64(1, 'D')
    if dt_days < 20:
        clim = ds.groupby('time.dayofyear').mean('time')
        clim = clim.sel(dayofyear=clim.dayofyear <= 365)
        phase = (clim.dayofyear.values - 0.5) / 365
        clim = clim.rename(dayofyear='phase')
    else:
        clim = ds.groupby('time.month').mean('time')
        phase = (clim.month.values - 0.5) / 12
        clim = clim.rename(month='phase')
    return clim.assign_coords(phase=phase)


def _clean(x, name):
    """concentrations/fluxes >= 0; missing values -> annual mean"""
    x = np.asarray(x, dtype=float)
    ok = np.isfinite(x)
    if not ok.any():
        raise ValueError(f'forcing {name} is all NaN')
    return np.maximum(np.where(ok, x, np.mean(x[ok])), 0.0)


def point_forcing(clim_btm_b, clim_tr_b):
    """one-year forcing dict for cbed2py at one point (on the cobalt_btm phases)"""
    t = clim_btm_b.phase.values
    frc = {'time': t}
    for v in clim_btm_b.data_vars:
        frc[v] = _clean(clim_btm_b[v].values, v)
    tp = clim_tr_b.phase.values
    for v in clim_tr_b.data_vars:
        frc[v] = np.interp(t, tp, _clean(clim_tr_b[v].values, v), period=1.0)
    return frc


def cbed_nk_wrapped(dsmom_a, clim_btm_a, clim_tr_a, dsporo_a, i, j, cont):
    """
    CBED at one point with method 'newton_krylov'. Returns (ds, diag):
    ds   the periodic state at phase 0 (Jan 1), same layout as cbed_wrapped
    diag spin-up diagnostics (DIAG_KEYS)
    """
    from cbed2py.spinup import newton_krylov_spinup

    print(f"Running point i={i} j={j} cont={cont} (newton_krylov)", flush=True)
    dsmom_b = dsmom_a.isel(xh=i, yh=j)
    frc = point_forcing(clim_btm_a.isel(xh=i, yh=j), clim_tr_a.isel(xh=i, yh=j))
    porosity = float(dsporo_a.isel(xh=i, yh=j)['porosity'].values)

    params = dict(CBED_CFG['params'])
    # constant values: annual means of the forcing (steady-state initial
    # condition, and J for jpoc mode 'longterm_mean')
    params.update({v: float(np.mean(frc[v])) for v in frc if v != 'time'})
    params.update({'temp':    max(0, float(dsmom_b['temp'].values)),
                   'S':       max(0, float(dsmom_b['salt'].values)),
                   'por.0':   max(0, porosity),
                   'por.inf': max(0, porosity),
                   'w':       0.0359,
                   'depth':   200.0})

    spin = {'max_newton': 10, 'verbose': False, **CBED_CFG['spinup']}
    y, hist, _ = newton_krylov_spinup(params, frc, **spin)
    h = hist[-1]
    diag = {'status': 1 if h['converged'] else 0, 'model_years': h['runs'],
            'newton_steps': h['newton'], 'imbalance': h['imbalance'], 'rms_G': h['rms_G']}
    if not h['converged']:
        print(f"NOT CONVERGED i={i} j={j} cont={cont}: imbalance {h['imbalance']:+.4f}, "
              f"|G|rms {h['rms_G']:.2e}", flush=True)

    N = y.size // len(STATE_NAMES)
    midx = pd.MultiIndex.from_product([np.arange(N), [cont]])
    df = pd.DataFrame(y.reshape(len(STATE_NAMES), N).T, columns=STATE_NAMES, index=midx)
    return df.to_xarray(), diag


def run_point(args):
    cont, iv, jv, i,j, valid_points = args
    # iv, jv are global indices that are comparable with valid_points set
    # i,j are local indices used to find the point in the dataset through isel
    diag = dict.fromkeys(DIAG_KEYS, np.nan)
    if (jv,iv) not in valid_points:
        print(i,j,'land')
        diag['status'] = -2
        return cont, empty_ds(cont), diag

    t0 = time.perf_counter()
    if CBED_CFG['method'] == 'newton_krylov':
        try:
            ds, d = cbed_nk_wrapped(dsmom_a, clim_btm_a, clim_tr_a, dsporo_a, i, j, cont)
            diag.update(d)
        except Exception as e:
            print(f"FAILED i={i} j={j} cont={cont}: {e}", flush=True)
            ds = empty_ds(cont)
            diag['status'] = -1
    else:
        ds = cbed_wrapped(dsmom_a, dscob_a, dscob2_a, dsporo_a, i, j, cont)
        diag['status'] = -1 if np.isnan(ds['OM1'].values).all() else 1
    diag['wall_time'] = time.perf_counter() - t0
    return cont, ds, diag


def chunk_bounds(chunk, ny, nx, n_chunks_x=2, n_chunks_y=2):
    """(x_start, x_end, y_start, y_end) of a chunk (0 = full domain)"""
    if chunk == 0:
        return 0, nx, 0, ny
    chunk_idx = chunk - 1
    ci = chunk_idx % n_chunks_x
    cj = chunk_idx // n_chunks_x

    x_edges = np.linspace(0, nx, n_chunks_x + 1, dtype=int)
    y_edges = np.linspace(0, ny, n_chunks_y + 1, dtype=int)
    return x_edges[ci], x_edges[ci + 1], y_edges[cj], y_edges[cj + 1]


def get_chunk(chunk, dsmom, dscob, dscob2, dsporo, ny, nx, n_chunks_x=2, n_chunks_y=2):
    if chunk == 0:
        return dsmom, dscob, dscob2, dsporo, np.arange(nx), np.arange(ny), np.arange(nx), np.arange(ny)

    x_start, x_end, y_start, y_end = chunk_bounds(chunk, ny, nx, n_chunks_x, n_chunks_y)

    ivx = np.arange(x_start, x_end)
    ivy = np.arange(y_start, y_end)
    ix  = np.arange(x_end - x_start)
    iy  = np.arange(y_end - y_start)

    dsmom_a  = dsmom.isel(xh=slice(x_start, x_end), yh=slice(y_start, y_end))
    dscob_a  = dscob.isel(xh=slice(x_start, x_end), yh=slice(y_start, y_end))
    dscob2_a = dscob2.isel(xh=slice(x_start, x_end), yh=slice(y_start, y_end))
    dsporo_a = dsporo.isel(xh=slice(x_start, x_end), yh=slice(y_start, y_end))

    return dsmom_a, dscob_a, dscob2_a, dsporo_a, ivx, ivy, ix, iy

if __name__ == '__main__':

    ROOT_DIR, FPATH, FTOPO, FGRD, FTEMPLATE, FOUT, CACHE_DIR, nx_chunk, ny_chunk = load_config()
    CBED_CFG = load_cbed_config()
    print(f"method={CBED_CFG['method']}", flush=True)

    nproc    = int(sys.argv[1])
    chunk    = int(sys.argv[2])



    script_path = osp.join(ROOT_DIR,'src/cbed_R/cbed_v1_func.R')

    # mr.read_variables will need three datasets containing mom6 and cobalt info
    # the datasets MUST BE NAMED AS: 
    # - cobalt_btm.nc
    # - cobalt_tr.nc
    # - mom6.nc
    # ncdump of the files I used can be found in CBED_R2PY/aux
    ds_dict = mr.read_variables(ROOT_DIR,FPATH,FTOPO, CACHE_DIR) 


    # ncdump of the files I used below can be found in CBED_R2PY/aux
    ds_poro = porosity_main(FGRD=FGRD,        # porosity file (subhadeep prepared it)
                            ROOTDIR=ROOT_DIR, # rood dir we are working at
                            FPATH=FTEMPLATE,  # template file for output
                            FTOPO=FTOPO,      # mom6 topography file
                            FOUT=FOUT)        # output/cache file
    ds_poro.load()
    print(ds_dict.keys(), flush=True)

    print("reading datasets")
    dscob  = ds_dict['dscobalt_btm'].mean(dim='time')
    dscob2 = ds_dict['dscobalt_tr'].mean(dim='time')
    dsmom  = ds_dict['dsmom']

    dscob.load()
    dscob2.load()
    dsmom.load()
    print(f"models loaded: {list(ds_dict.keys())}", flush=True)

    ny, nx = dscob['btm_o2'].values.shape
    print(f"grid shape: ny={ny} nx={nx}", flush=True)

    dsmom_a, dscob_a, dscob2_a, dsporo_a,ivx, ivy, ix, iy = get_chunk(chunk,
                                                   dsmom,
                                                   dscob,
                                                   dscob2,
                                                   ds_poro,
                                                   ny,
                                                   nx,
                                                   n_chunks_x=nx_chunk,
                                                   n_chunks_y=ny_chunk)

    print(f"chunk={chunk} ix={ix[[0,-1]]} iy={iy[[0,-1]]}", flush=True)

    # one-year forcing climatology of the chunk (method 'newton_krylov' only)
    clim_btm_a = clim_tr_a = None
    if CBED_CFG['method'] == 'newton_krylov':
        x0, x1, y0, y1 = chunk_bounds(chunk, ny, nx, n_chunks_x=nx_chunk, n_chunks_y=ny_chunk)
        sub = dict(xh=slice(x0, x1), yh=slice(y0, y1))
        btm, tr = forcing_fields(ds_dict['dscobalt_btm'].isel(**sub),
                                 ds_dict['dscobalt_tr'].isel(**sub))
        clim_btm_a = annual_climatology(btm).astype('float32').load()
        clim_tr_a  = annual_climatology(tr).astype('float32').load()
        print(f"forcing climatology: {clim_btm_a.sizes['phase']} phases (cobalt_btm), "
              f"{clim_tr_a.sizes['phase']} phases (cobalt_tr)", flush=True)

    # identifying indices with land points
    jivalid = np.where(~np.isnan(dscob_a['btm_o2'].values))  # selecting mask in local index
    jvm, ivm = np.meshgrid(ivy, ivx, indexing='ij') #  global indices of the subset
    
    # mapping jivalid onto global indices
    jvec = jvm[jivalid] 
    ivec = ivm[jivalid]
    # local valid water indices (used to identify water points in the loop)
    valid_points = set(zip(jvec, ivec))  


    jm, im = np.meshgrid(iy, ix, indexing='ij')    # local indices


    
    
    # -- parallel implementation --
    args = [(cont, iv, jv, i, j, valid_points)
            for cont, (iv, jv, i, j) in enumerate(zip(ivm.ravel(),  #
                                                      jvm.ravel(),
                                                      im.ravel(),
                                                      jm.ravel()))]
    
    print(f"starting pool: nproc={nproc}", flush=True)
    with Pool(processes=nproc, initializer=init_worker, maxtasksperchild=20) as pool:
        results = pool.map(run_point, args)

    ds    = {cont: d for cont, d, _ in results}
    diags = [g for _, _, g in results]

    # -- grouping results --


    test_dict = [ds[i] for i in ds]

    # concat flat list along a new dimension
    combined = xr.concat(test_dict, dim="level_1")

    # reshape each variable into (level_0, level_1, level_2) = (20, 2, 2)
    ds_3d = xr.Dataset(
        {name: (["level_0", "level_1", "level_2"],
                combined[name].values.reshape(20, iy.size, ix.size))
        for name in combined.data_vars},
        coords={
            "level_0": combined.level_0.values,
            "level_1": iy,
            "level_2": ix
        }
    )

    ds_3d.attrs['cbed_method'] = CBED_CFG['method']
    if CBED_CFG['method'] == 'newton_krylov':
        ds_3d.attrs['cbed2py_model'] = CBED_CFG['model']
        ds_3d.attrs['state_phase'] = 'periodic state at phase 0 of the forcing year (Jan 1)'
    ds_3d.to_netcdf(osp.join(ROOT_DIR,f'data/cache/cbed_mom6_{chunk:03d}.nc'))

    # per-point diagnostics (status, wall time, spin-up convergence). Not named
    # cbed*.nc so that cbed_stitch.py's cbed_pattern does not pick it up.
    ds_diag = xr.Dataset(
        {k: (["level_1", "level_2"],
             np.array([g[k] for g in diags], dtype=float).reshape(iy.size, ix.size))
         for k in DIAG_KEYS},
        coords={"level_1": iy, "level_2": ix},
        attrs={'cbed_method': CBED_CFG['method'],
               'status': '1 converged, 0 not converged, -1 failed, -2 land'})
    ds_diag.to_netcdf(osp.join(ROOT_DIR, f'data/cache/spinup_diag_{chunk:03d}.nc'))
    ok = ds_diag['status'].values
    print(f"points: {np.sum(ok == 1)} ok, {np.sum(ok == 0)} not converged, "
          f"{np.sum(ok == -1)} failed, {np.sum(ok == -2)} land; "
          f"mean wall time {np.nanmean(ds_diag['wall_time'].values):.2f} s/point", flush=True)
    