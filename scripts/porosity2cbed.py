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


    def reconfigure(self, dsout=None, fout=None, save=False):
        """Reformat regridded dataset to ncdump-compliant NetCDF3."""
        if dsout is None:
            dsout = xr.open_dataset(self.fout)
        if fout is None:
            fout = self.fout.replace('.nc', '_reconfigured.nc')

        FILL = 0.8111111

        # rename xh/yh back to lat/lon if needed
        rename = {k: v for k, v in {'xh': 'lon', 'yh': 'lat'}.items() if k in dsout.dims}
        if rename:
            dsout = dsout.rename(rename)

        if 'time' not in dsout.dims:
            import pandas as pd
            t = pd.date_range('1993-01-01', periods=1, freq='YS')
            dsout = dsout.expand_dims(time=t)

        por = dsout['porosity'].astype('float32')
        por = por.where(por.notnull(), other=FILL)

        dsout = xr.Dataset(
            {'porosity': por},
            coords={'time': dsout['time'], 'lat': dsout['lat'], 'lon': dsout['lon']}
        )

        dsout['lat'].attrs  = {'units': 'degrees_north', 'long_name': 'lat', 'cartesian_axis': 'Y'}
        dsout['lon'].attrs  = {'units': 'degrees_east',  'long_name': 'lon', 'cartesian_axis': 'X'}
        dsout['time'].attrs = {'long_name': 'time'}
        dsout['porosity'].attrs = {
            'units':      'unitless',
            'long_name':  'Seafloor sediment porosity',
            '_FillValue': FILL,
        }

        dsout['porosity'] = dsout['porosity'].transpose('time', 'lat', 'lon')

        if save:
            encoding = {
                'time': {'dtype': 'float64', 'calendar': 'NOLEAP'},
            }
            dsout.to_netcdf(fout, format='NETCDF3_CLASSIC', encoding=encoding, unlimited_dims=['time'])

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
    ds = pr.load(save=False, usecache=True)
    pr.reconfigure(dsout=ds, fout=osp.join(ROOTDIR, 'data/cache/porosity_final.nc'), save=True)
    return ds

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