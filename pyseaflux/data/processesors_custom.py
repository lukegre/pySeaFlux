"""
Contains processors that are specific to the specific datasets
"""

import numpy as np
import xarray as xr


def era5_gcp_raw_fname_mapper(s):
    import re
    
    pattern = re.compile(r'(\d{4})/(\d{2})/(\d{2})/(\w+)/surface.nc')
    yyyy, mm, dd, var = re.findall(pattern, s)[0]

    fname = f'{yyyy}/era5-{yyyy}{mm}{dd}-{var}.nc'

    return fname


def compute_wind_moments(ds):

    u = ds.u10
    v = ds.v10
    second = (u**2 + v**2).expand_dims(moment=[2])
    third  = (u**3 + v**3).expand_dims(moment=[3])
    first  = (second**0.5).assign_coords(moment=[1])

    moments = (
        xr.concat([first, second, third], dim='moment')
        .assign_attrs(
            long_name='Wind moments',
            units='m^moment/s^moment',
            description='Wind moments computed from u10 and v10',
        )
    )

    return moments


def era5_coarsen_to_1deg(ds_025deg):

    s = 1
    lat = np.arange(-90 + s/2, 90, s)
    lon = np.arange(-180 + s/2, 182, s)

    lat025 = ds_025deg.lat.values
    lon025 = ds_025deg.lon.values

    assert (lat025.size > 600) & (lat025.size < 740), 'Latitude dimension is not correct'
    assert (lon025.size > 1430) & (lon025.size < 1480), 'Longitude dimension is not correct'

    n_ext = 12
    lon180_ext = np.hstack([lon025, lon025[:n_ext]])
    lon360_ext = np.hstack([lon025, lon025[:n_ext] + 360])
    
    ds_025deg = ds_025deg.sel(lon=lon180_ext).assign_coords(lon=lon360_ext)

    ds_1deg_offset = ds_025deg.coarsen(lat=4, lon=4, boundary='pad').mean()
    ds_1deg = ds_1deg_offset.interp(lat=lat, lon=lon, method='linear')
    ds_1deg = ds_1deg.sel(lat=slice(-90, 90), lon=slice(-180, 180))

    return ds_1deg


def ncep_interp(ds):

    s = 1
    lat = np.arange(-90 + s/2, 90, s)
    lon = np.arange(-180 + s/2, 182, s)

    ds_interp = ds.interp(lat=lat, lon=lon, method='linear')
    ds_interp = ds.sel(lat=slice(-90, 90), lon=slice(-180, 180))

    return ds_interp
