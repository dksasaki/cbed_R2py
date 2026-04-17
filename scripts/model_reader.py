
import os
import os.path as osp
import glob
import xarray as xr
import numpy as np
import logging
from dask.distributed import Client, LocalCluster
from dask.distributed import progress, wait

log = logging.getLogger(__name__)

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

# if 'client' not in globals():
#     cluster = LocalCluster(n_workers=16,
#                            threads_per_worker=1,
#                            memory_limit='32GB')
#     client = Client(cluster)

_client = None
_cluster = None

def get_client(n_workers=16, threads_per_worker=1, memory_limit='32GB'):
    global _client, _cluster
    if _client is None:
        _cluster = LocalCluster(n_workers=n_workers,
                                threads_per_worker=threads_per_worker,
                                memory_limit=memory_limit)
        _client = Client(_cluster)
    return _client


def julian2npdatetime(ds):
    """Convert time coordinate from Julian calendar (cftime) to numpy datetime64."""
    time = np.array([np.datetime64(i) for i in ds.time.values]).astype('datetime64[ns]')
    return ds.assign_coords(time=time)

def read_mom6cobalt(paths, ftopo, varbs, chunk_dict=None, topog=False):
    _ = get_client()

    def _assert_variables(ds, varbs, path):
        assert_list = []
        for  v in varbs:
            if v in ds:
                pass
            else:
                assert_list.append(v)
        assert len(assert_list)==0, f"{assert_list} not found in {path}"

    def _merge_topo(ftopo, dsmom):
        dstopo = xr.open_dataset(ftopo)
        dstopo = dstopo.rename(nx='xh', ny='yh')
        dstopo = dstopo.assign_coords(xh=dsmom.xh, yh=dsmom.yh)
        
        return xr.merge([dsmom,dstopo])
    

    assert(paths),"must be a string"
    
    f = lambda x: sorted(glob.glob(x))
    paths = f(paths)
    datasets = []
    for path in paths:  # we use a for instead of mfdataset to avoid a memory issue
        print(path)
        ds = xr.open_dataset(path, lock=False, decode_timedelta=True)

        _assert_variables(ds, varbs, path)

        ds = julian2npdatetime(ds)
        aux = ds[varbs]

        # make sure that cobalt3d, mom63d and mom62d work with this function
        if chunk_dict is None:
            if 'z_l' in ds:
                chunk_dict = dict(xh=100, yh=100, time=20, z_l=5)
            elif 'z_l' in ds:
                chunk_dict = dict(xh=100, yh=100, time=20, zl=5)
            else:
                chunk_dict = dict(xh=100, yh=100, time=20)
        else:
            assert type(chunk_dict) is dict,"chunk_dict must either None or type(dict)"
            
            for i in chunk_dict:
                assert i in ds, f"{i} not found i {path}"

        datasets.append(aux.chunk(**chunk_dict))

    dscobalt = xr.concat(datasets,
                   dim='time',
                   coords='minimal',
                   compat='override',
                   data_vars='minimal')

    if topog:
        dscobalt = _merge_topo(ftopo, dscobalt)
    return dscobalt

def _variables_model():
    varbs_cobalt_btm = [
        # Dissolved inorganic
        "btm_o2",           # oxygen
        "btm_no3",          # nitrate
        # "btm_nh4",          # ammonium
        "btm_dic",          # dissolved inorganic carbon
        # "btm_alk",          # alkalinity

        # Nitrogen fluxes
        "fndet_btm",        # detrital N
        "fndi_btm",         # diazotroph N
        "fnlg_btm",         # large phyto N
        "fnmd_btm",         # medium phyto N
        "fnsm_btm",         # small phyto N

        # Silicon fluxes
        "fsidet_btm",       # detrital Si
        "fsimd_btm",        # medium phyto Si
        "fsilg_btm",        # large phyto Si

        # Iron fluxes
        "ffedet_btm",       # detrital Fe
        "ffedi_btm",        # diazotroph Fe
        "ffesm_btm",        # small phyto Fe
        "ffemd_btm",        # medium phyto Fe
        "ffelg_btm",        # large phyto Fe

        # Carbonate fluxes
        "fcadet_arag_btm",  # aragonite detrital Ca
        "fcadet_calc_btm",  # calcite detrital Ca

        # Phosphorus fluxes
        "fpdet_btm",        # detrital P
        "fpdi_btm",         # diazotroph P
        "fpsm_btm",         # small phyto P
        "fpmd_btm",         # medium phyto P
        "fplg_btm",         # large phyto P

        # Lithogenic
        "flithdet_btm",     # lithogenic detritus
    ]
    varbs_mom6   = ['temp', 'salt']
    varbs_cobalt_tr = ['nh4', 'talk']
    return varbs_cobalt_btm, varbs_mom6, varbs_cobalt_tr

def read_variables(root_dir: str, fpath: str, ftopo: str, cache_dir:str) -> tuple:
    """
    Load and cache MOM6/COBALT model output datasets — bottom COBALT 
    tracers, surface COBALT tracers, and MOM6 ocean daily fields — from raw 
    NetCDF files into a local cache directory. On first run, each dataset is 
    computed via Dask, persisted in memory, and saved to disk; on subsequent 
    runs the cached files are loaded directly, skipping the expensive 
    computation. Returns the three datasets as a tuple.

    Variables in each dataset are controlled by `_variables_model()`. You
    can call it from your script directly for explicit information.


    Parameters
    ----------
    root_dir : str
        Root directory for cache output.
    fpath : str
        Directory containing raw model NetCDF files.
    ftopo : str
        Path to ocean topography file.
    cache_dir: str
        Path to cache

    Returns
    -------
    tuple of xarray.Dataset
        (dsmom, dscobalt_btm, dscobalt_tr)

    Raises
    ------
    FileNotFoundError
        If fpath or ftopo do not exist.
    """
    _ = get_client()

    if not osp.exists(fpath):
        raise FileNotFoundError(f"Model output path not found: {fpath}")
    if not osp.exists(ftopo):
        raise FileNotFoundError(f"Topography file not found: {ftopo}")

    os.makedirs(cache_dir, exist_ok=True)
    log.info(f"Cache directory: {cache_dir}")

    varbs_cobalt_btm, varbs_mom6, varbs_cobalt_tr = _variables_model()

    def _persist_and_save(ds, fsave):
        ds1 = ds.persist()
        progress(ds1)
        wait(ds1)
        log.info(f"Saving to {fsave}")
        ds1.to_netcdf(fsave)
        return ds1

    def _cached(fsave, fn):
        if glob.glob(fsave):
            log.info(f"Loading from cache: {fsave}")
            return xr.open_dataset(fsave, chunks={})
        log.info(f"Cache miss, computing: {fsave}")
        return fn(fsave)

    def _read_cobalt_btm(fsave):
        fpaths = osp.join(fpath, '*cobalt_btm.nc')
        ds = read_mom6cobalt(fpaths, ftopo, varbs_cobalt_btm)
        return _persist_and_save(ds, fsave)

    def _read_cobalt_tr(fsave):
        fpaths = osp.join(fpath, '*cobalt_tracers.nc')
        ds = read_mom6cobalt(fpaths, ftopo, varbs_cobalt_tr, chunk_dict={'z_l': 52})
        ds = ds.ffill(dim='z_l').bfill(dim='z_l').isel(z_l=-1)
        return _persist_and_save(ds, fsave)

    def _read_mom(fsave):
        fpaths = osp.join(fpath, '*ocean_daily.nc')
        ds = read_mom6cobalt(fpaths, ftopo, varbs_mom6, chunk_dict={'zl': 1})
        ds = ds.isel(zl=-1).mean(dim='time')
        return _persist_and_save(ds, fsave)

    dscobalt_btm = _cached(osp.join(cache_dir, 'cobalt_btm.nc'), _read_cobalt_btm)
    dscobalt_tr  = _cached(osp.join(cache_dir, 'cobalt_tr.nc'),  _read_cobalt_tr)
    dsmom        = _cached(osp.join(cache_dir, 'mom6.nc'),        _read_mom)

    ds_dict = {}
    ds_dict['dsmom']        = dsmom
    ds_dict['dscobalt_btm'] = dscobalt_btm
    ds_dict['dscobalt_tr']  = dscobalt_tr
    return ds_dict


if __name__ =='__main__':
    
    ROOT_DIR = '/projects/schultz/d.sasaki/km_scale_model/' + \
                'mom6cobalt_25th/20240723_zstar/tasks/' + \
                '202603_cbed_R2py'
    FPATH     = '/home/d.sasaki/scratch/mom_experiments/cbed_test_001/outputs_raw'
    FTOPO     = '/home/d.sasaki/schultz/d.sasaki/km_scale_model/mom6cobalt_25th/mom_tools/data/grid/nwa25_interped/netcdf3/ocean_topog.nc'
    CACHE_DIR = osp.join(ROOT_DIR, 'data/cache/scratch_test')

    os.chdir(ROOT_DIR)

    ds_dict = read_variables(ROOT_DIR, FPATH, FTOPO, CACHE_DIR)



