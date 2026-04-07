import xarray as xr
import model_reader as mr
import os
import os.path as osp
import xesmf as xe
import glob



class PorosityRegridder:
    def __init__(self, fpath, ftopo, fgrd, fout):
        self.fpath = fpath
        self.ftopo = ftopo
        self.fgrd  = fgrd
        self.fout  = fout

    def load(self, save=False, usecache=False):
        if usecache and glob.glob(self.fout):
            print('Loading cached file.')
            return xr.open_dataset(self.fout)

        ds1    = mr.read_mom6cobalt(self.fpath, self.ftopo, varbs=['temp'])
        has_xy = self.fgrd.endswith('.grd')
        chunks = dict(x=400, y=400) if has_xy else dict(lon=200, lat=200)
        dsgrd  = xr.open_dataset(self.fgrd, chunks=chunks)

        lon_dim = 'x'   if has_xy else 'lon'
        lat_dim = 'y'   if has_xy else 'lat'
        varb    = 'z'   if has_xy else 'porosity'
        factr   = 1e-2  if has_xy else 1

        ds = dsgrd.sel({
            lon_dim: slice(ds1.xh.min() - 1, ds1.xh.max() + 1),
            lat_dim: slice(ds1.yh.min() - 1, ds1.yh.max() + 1),
        })

        if has_xy:
            ds = ds.rename(x='lon', y='lat')

        reg   = xe.Regridder(ds, ds1.rename(xh='lon', yh='lat'), 'bilinear')
        dsout = reg(ds)

        dsout[varb] = (dsout[varb]
                       .ffill(dim='lon').ffill(dim='lat')
                       .bfill(dim='lon').bfill(dim='lat')) * factr

        if varb != 'porosity':
            dsout = dsout.rename({varb: 'porosity'})

        dsout.load()
        dsout = dsout.rename({'lon': 'xh', 'lat': 'yh'})

        if save:
            dsout.to_netcdf(self.fout)

        return dsout


def porosity_main(save=False, usecache=True):
    print('WARNING, paths are hardcoded in porosity_main function [in porosity2cbed.py]')
    FGRD = '/home/d.sasaki/schultz/data/cbed_supporting_data/subhadeep/globalporosity_map.grd'
    ROOTDIR   = '/projects/schultz/d.sasaki/km_scale_model/mom6cobalt_25th/20240723_zstar/tasks/202603_cbed_R2py'
    FPATH     = '/home/d.sasaki/scratch/mom_experiments/cbed_test_001/outputs_raw/19930101.ocean_daily.nc'
    FTOPO     = '/home/d.sasaki/schultz/d.sasaki/km_scale_model/mom6cobalt_25th/mom_tools/data/grid/nwa25_interped/netcdf3/ocean_topog.nc'
    FOUT = osp.join(ROOTDIR,'data/cache/porosity_neus25.nc')
    
    pr    = PorosityRegridder(
        fpath=FPATH,
        ftopo=FTOPO,
        fgrd=FGRD,
        fout=FOUT
    )
    dsout = pr.load(save=False, usecache=True)

    return dsout

if __name__ == '__main__':
    porosity_main()
    
# FGRD = '/home/d.sasaki/schultz/data/cbed_supporting_data/subhadeep/seafloor_porosity.nc'

# def porosity_file( save=False, usecache=False):
#     if usecache:
#         if glob.glob(FOUT):
#             print('Loading cached file.')
#             return xr.open_dataset(FOUT)

#     ds1   = mr.read_mom6cobalt(FPATH, FTOPO, varbs=['temp'])
#     chunks = dict(x=400, y=400) if FGRD.endswith('.grd') else dict(lon=200, lat=200)
#     dsgrd  = xr.open_dataset(FGRD, chunks=chunks)

#     # Detect coordinate names
#     has_xy  = 'x' in dsgrd.coords and 'y' in dsgrd.coords
#     lon_dim = 'x'  if has_xy else 'lon'
#     lat_dim = 'y'  if has_xy else 'lat'
#     varb    = 'z'  if has_xy else 'porosity'
#     factr   = 1 if FGRD.endswith('.nc') else 1e-2

#     ds = dsgrd.sel({
#         lon_dim: slice(ds1.xh.min() - 1, ds1.xh.max() + 1),
#         lat_dim: slice(ds1.yh.min() - 1, ds1.yh.max() + 1),
#     })

#     if has_xy:
#         ds = ds.rename(x='lon', y='lat')

#     reg   = xe.Regridder(ds, ds1.rename(xh='lon', yh='lat'), 'bilinear')
#     dsout = reg(ds)

#     dsout[varb] = (dsout[varb]
#                    .ffill(dim='lon').ffill(dim='lat')
#                    .bfill(dim='lon').bfill(dim='lat')) * factr

#     if varb != 'porosity':
#         dsout = dsout.rename({varb: 'porosity'})

#     dsout.load()
#     dsout = dsout.rename({'lon': 'xh', 'lat':'yh'})

#     if save:
#         dsout.to_netcdf(fout)

#     return dsout

# def porosity_main(save=False, usecache=True):
#     print('WARNING, paths are hardcoded in porosity_main function [in porosity2cbed.py]')
#     FGRD = '/home/d.sasaki/schultz/data/cbed_supporting_data/subhadeep/globalporosity_map.grd'

#     ROOTDIR   = '/projects/schultz/d.sasaki/km_scale_model/mom6cobalt_25th/20240723_zstar/tasks/202603_cbed_R2py'
#     FPATH     = '/home/d.sasaki/scratch/mom_experiments/cbed_test_001/outputs_raw/19930101.ocean_daily.nc'
#     FTOPO     = '/home/d.sasaki/schultz/d.sasaki/km_scale_model/mom6cobalt_25th/mom_tools/data/grid/nwa25_interped/netcdf3/ocean_topog.nc'
#     FOUT = osp.join(ROOTDIR,'data/cache/porosity_neus25.nc')
#     dsout = porosity_file(save=save, usecache=usecache)
#     return dsout

# if __name__ == '__main__':
    
#     dsout = porosity_file(save=True)