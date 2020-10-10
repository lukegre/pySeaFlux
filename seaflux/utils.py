def earth_radius(lat):
    from numpy import deg2rad, sin, cos

    lat = deg2rad(lat)
    a = 6378137
    b = 6356752
    r = (
        ((a ** 2 * cos(lat)) ** 2 + (b ** 2 * sin(lat)) ** 2)
        / ((a * cos(lat)) ** 2 + (b * sin(lat)) ** 2)
    ) ** 0.5

    return r


def area_grid(lat, lon, return_dataarray=False):
    """Calculate the area of each grid cell for a user-provided
    grid cell resolution. Area is in square meters, but resolution
    is given in decimal degrees.
    Based on the function in
    https://github.com/chadagreene/CDT/blob/master/cdt/cdtarea.m
    """
    from numpy import meshgrid, deg2rad, gradient, cos

    ylat, xlon = meshgrid(lat, lon)
    R = earth_radius(ylat)

    dlat = deg2rad(gradient(ylat, axis=1))
    dlon = deg2rad(gradient(xlon, axis=0))

    dy = dlat * R
    dx = dlon * R * cos(deg2rad(ylat))

    area = dy * dx

    if not return_dataarray:
        return area
    else:
        from xarray import DataArray

        xda = DataArray(
            area.T,
            dims=["lat", "lon"],
            coords={"lat": lat, "lon": lon},
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
        else:
            return out

        if (xda is None) & second_isdict:
            return out[0]

        attrs = xda.attrs
        if second_isdict:
            data = out[0]
            attrs.update(out[1])
        else:
            data = out
            attrs = xda.attrs

        out = xr.DataArray(data=data, dims=xda.dims, coords=xda.coords, attrs=attrs)
        return out

    return wrapper
