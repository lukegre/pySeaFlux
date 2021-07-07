"""
Area calculations
-----------------

Calculates the area of pixels for a give grid input.
"""


def earth_radius(lat):
    """Calculate the radius of the earth for a given latitude

    Args:
        lat (array, float): latitude value (-90 : 90)

    Returns:
        array: radius in metres
    """
    from numpy import cos, deg2rad, sin

    lat = deg2rad(lat)
    a = 6378137
    b = 6356752
    r = (
        ((a ** 2 * cos(lat)) ** 2 + (b ** 2 * sin(lat)) ** 2)
        / ((a * cos(lat)) ** 2 + (b * sin(lat)) ** 2)
    ) ** 0.5

    return r


def area_grid(lat, lon, return_dataarray=False):
    """Calculate the area of each grid cell for given lats and lons

    Args:
        lat (array): latitudes in decimal degrees of length N
        lon (array): longitudes in decimal degrees of length M
        return_dataarray (bool, optional): if True returns xr.DataArray, else array

    Returns:
        array, xr.DataArray: area of each grid cell in meters

    References:
        https://github.com/chadagreene/CDT/blob/master/cdt/cdtarea.m
    """
    from numpy import cos, deg2rad, gradient, meshgrid

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
            attrs=dict(
                long_name="Area per pixel",
                units="m^2",
                description=(
                    "Area per pixel as calculated by SeaFlux. The non-"
                    "spherical shape of Earth is taken into account."
                ),
            ),
        )

        return xda


def get_area_from_dataset(dataarray, lat_name="lat", lon_name="lon"):
    """
    Calculate the grid cell area from a xr.Dataset or xr.DataArray.
    """
    da = dataarray
    x = da.lon.values
    y = da.lat.values

    area = area_grid(y, x, return_dataarray=True)

    return area
