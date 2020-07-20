def earth_radius(lat):
    from numpy import deg2rad, sin, cos
    lat = deg2rad(lat)
    a = 6378137
    b = 6356752
    r = (((a**2 * cos(lat))**2 + (b**2 * sin(lat))**2) /
        ((a * cos(lat))**2 + (b * sin(lat))**2))**0.5

    return r


def area_grid(lat, lon, return_dataarray=False):
    """
    Calculate the area of each grid cell for a user-provided
    grid cell resolution.

    Based on the function in
    https://github.com/chadagreene/CDT/blob/master/cdt/cdtarea.m

    Parameters
    ----------
    lat : array (1D)
        latitudes ranging between -90 and 90
    lon : array (1D)
        longitudes ranging between (-180, 180)
    return_dataarray : bool (False)
        if True, returns a xr.DataArray with lat/lon as coords

    Returns
    -------

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
            coords={
                "lat": lat,
                "lon": lon,
            },
            attrs={
                "long_name": "area_per_pixel",
                "description": "area per pixel",
                "units": "m^2",
            },
        )
        return xda
