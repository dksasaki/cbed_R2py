"""Write time-mean fntot_btm (mol N m-2 s-1) on the full regional grid, for data_override.

Usage: python make_fntot_mean.py <cobalt_btm.nc> <out.nc>

fntot_btm = fndet_btm + fndi_btm + fnsm_btm + fnmd_btm + fnlg_btm (as in cbed_wrapper.py).
Layout follows porosity_final.nc: time(1), lat, lon; data_override interpolates it to the model grid.
Land (NaN) is set to 0.
"""
import sys

import numpy as np
import xarray as xr

cob, out = sys.argv[1:3]

comps = ["fndet_btm", "fndi_btm", "fnsm_btm", "fnmd_btm", "fnlg_btm"]
with xr.open_dataset(cob) as c:
    fntot = sum(c[v].astype("f8") for v in comps).load()
    lon = c["xh"].values
    lat = c["yh"].values

print(f"grid {lon.size} x {lat.size}, land/NaN cells set to 0: {int(fntot.isnull().sum())}")
a = fntot.fillna(0.0).values
assert np.all(a >= 0), "negative fntot_btm"

time = xr.DataArray([0.0], dims="time",
                    attrs={"long_name": "time", "units": "days since 1993-01-01 00:00:00", "calendar": "NOLEAP"})
da = xr.DataArray(
    a[np.newaxis].astype("f4"),
    dims=("time", "lat", "lon"),
    coords={
        "time": time,
        "lat": ("lat", lat, {"units": "degrees_north", "long_name": "lat", "cartesian_axis": "Y"}),
        "lon": ("lon", lon, {"units": "degrees_east", "long_name": "lon", "cartesian_axis": "X"}),
    },
    attrs={"units": "mol m-2 s-1",
           "long_name": "time-mean total N sinking flux to bottom (fndet+fndi+fnsm+fnmd+fnlg)"},
)
ds = da.to_dataset(name="fntot_btm")
ds.attrs = {"source": cob, "history": "written by make_fntot_mean.py"}

encoding = {v: {"_FillValue": None} for v in ["time", "lat", "lon", "fntot_btm"]}
ds.to_netcdf(out, format="NETCDF3_64BIT", unlimited_dims=["time"], encoding=encoding)
print(f"wrote {out}")
