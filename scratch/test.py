
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

def read_mom6cobalt(paths, ftopo, varbs):
    
    assert(paths),"must be a string"
    
    f = lambda x: sorted(glob.glob(x))
    paths = f(paths)
    datasets = []
    for path in paths:  # we use a for instead of mfdataset to avoid a memory issue
        print(path)
        ds = xr.open_dataset(path, lock=False)
        ds = julian2npdatetime(ds)
        aux = ds[varbs]

        # make sure that cobalt3d, mom63d and mom62d work with this function
        if 'z_l' in ds:
            chunck_dict = dict(xh=100, yh=100, time=20, z_l=5)
        elif 'z_l' in ds:
            dict(xh=100, yh=100, time=20, zl=5)
        else:
            dict(xh=100, yh=100, time=20)
        
        datasets.append(aux.chunk(chunck_dict))

    dscobalt = xr.concat(datasets,
                   dim='time',
                   coords='minimal',
                   compat='override',
                   data_vars='minimal')
    return dscobalt

def _variables_model():
    varbs_cobalt = [
        # Dissolved inorganic
        "btm_o2",           # oxygen
        "btm_no3",          # nitrate
        "btm_nh4",          # ammonium
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
    return varbs_cobalt, varbs_mom6

f = lambda x: sorted(glob.glob(x))

fpath = '/home/d.sasaki/scratch/mom_experiments/cbed_test_001/outputs_raw0'
ftopo = '/home/d.sasaki/schultz/d.sasaki/km_scale_model/mom6cobalt_25th/mom_tools/data/grid/nwa25_interped/netcdf3/ocean_topog.nc'

varbs_cobalt, varbs_mom6 = _variables_model()


# fpaths_cobalt     = osp.join(fpath,'*sediment_cbed.nc')
# dscobalt          = read_mom6cobalt(fpath, ftopo, varbs_cobalt)

fpaths_mom6       = osp.join(fpath,'*ocean_daily.nc')
dsmom             = read_mom6cobalt(fpath, ftopo, varbs_mom6)