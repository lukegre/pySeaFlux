

def area_grid(resolution=1):
    """Calculate the area of each grid cell for a user-provided
    grid cell resolution. Area is in square meters, but resolution
    is given in decimal degrees."""
    import numpy as np
    import xarray as xr

    # Calculations needs to be in radians
    lats = np.deg2rad(np.arange(-90, 90.1, resolution))
    r_sq = 6371000 ** 2
    n_lats = int(360.0 / resolution)
    area = (
        r_sq
        * np.ones(n_lats)[:, None]
        * np.deg2rad(resolution)
        * (np.sin(lats[1:]) - np.sin(lats[:-1]))
    )
    xda = xr.DataArray(
        area.T,
        dims=["lat", "lon"],
        coords={
            "lat": np.arange(-90 + 0.5, 90),
            "lon": np.arange(-180 + 0.5, 180),
        },
        attrs={
            "long_name": "area_per_pixel",
            "description": "area per pixel",
            "units": "m^2",
        },
    )

    return xda


def preserve_xda(func):
    """
    Function wrapper that will return output to xr.DataArray
    if any of the inputs are DataArrays.

    The coordinates from the first input term that is a DataArray will be
    used to define the output DataArray

    If the output is a tuple and the second output is a dictionary, this
    output will be used to create the attributes of the DataArray.
    """
    from functools import wraps
    import xarray as xr

    @wraps(func)
    def wrapper(*args, **kwargs):
        xda = None
        for a in args:
            if isinstance(a, xr.DataArray):
                xda = a
                break

        out = func(*args, **kwargs)
        istuple = isinstance(out, tuple)
        if istuple:
            second_isdict = isinstance(out[1], dict)

        if second_isdict:
            data = out[0]
            attrs = out[1]
        else:
            data = out
            attrs = {}

        out = xr.DataArray(
            data=data,
            dims=xda.dims,
            coords=xda.coords,
            attrs=attrs
        )
        return out

    return wrapper