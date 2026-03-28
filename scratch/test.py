
import os
import os.path as osp
import glob
import xarray as xr
import numpy as np

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

from dask.distributed import Client, LocalCluster
from dask.distributed import progress

if 'client' not in globals():
    cluster = LocalCluster(n_workers=16,
                           threads_per_worker=1,
                           memory_limit='32GB')
    client = Client(cluster)


def julian2npdatetime(ds):
    """Convert time coordinate from Julian calendar (cftime) to numpy datetime64."""
    time = np.array([np.datetime64(i) for i in ds.time.values]).astype('datetime64[ns]')
    return ds.assign_coords(time=time)

def read_mom6cobalt(paths, ftopo, varbs, chunk_dict=None):
    
    def _assert_variables(ds, varbs, path):
        assert_list = []
        for  v in varbs:
            if v in ds:
                pass
            else:
                assert_list.append(v)
        assert len(assert_list)==0, f"{assert_list} not found in {path}"


    assert(paths),"must be a string"
    
    f = lambda x: sorted(glob.glob(x))
    paths = f(paths)
    datasets = []
    for path in paths:  # we use a for instead of mfdataset to avoid a memory issue
        print(path)
        ds = xr.open_dataset(path, lock=False)

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
    return dscobalt

def _variables_model():
    varbs_cobalt_btm = [
        # Dissolved inorganic
        "btm_o2",           # oxygen
        "btm_no3",          # nitrate
        # "btm_nh4",          # ammonium
        "btm_dic",          # dissolved inorganic carbon
        "btm_alk",          # alkalinity

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
    varbs_cobalt_tr = ['nh4']
    return varbs_cobalt_btm, varbs_mom6, varbs_cobalt_tr

def read_variables():
    f = lambda x: sorted(glob.glob(x))

    fpath = '/home/d.sasaki/scratch/mom_experiments/cbed_test_001/outputs_raw'
    ftopo = '/home/d.sasaki/schultz/d.sasaki/km_scale_model/mom6cobalt_25th/mom_tools/data/grid/nwa25_interped/netcdf3/ocean_topog.nc'

    varbs_cobalt_btm, varbs_mom6, varbs_cobalt_tr = _variables_model()


    fpaths_cobalt     = osp.join(fpath,'*cobalt_btm.nc')
    dscobalt_btm          = read_mom6cobalt(fpaths_cobalt, ftopo, varbs_cobalt_btm)


    fpaths_cobalt     = osp.join(fpath,'*cobalt_tracers.nc')
    dscobalt_tr       = read_mom6cobalt(fpaths_cobalt, ftopo, varbs_cobalt_tr, chunk_dict={'z_l':1})
    dscobalt_tr       = dscobalt_tr.ffill(dim='z_l') \
                                .bfill(dim='z_l')
    dscobalt_tr       = dscobalt_tr.isel(z_l=-1)

    fpaths_mom6       = osp.join(fpath,'*ocean_daily.nc')
    dsmom             = read_mom6cobalt(fpaths_mom6, ftopo, varbs_mom6, chunk_dict={'zl':1})
    dsmom             = dsmom.isel(zl=-1)
    dsmom             = dsmom.mean(dim='time')

    return dsmom, dscobalt_btm, dscobalt_tr


if __name__ =='__main__':
    dsmom, dscobalt_btm, dscobalt_tr = read_variables()
    dsmom.load()
    dscobalt_btm.load()
    dscobalt_tr.load()