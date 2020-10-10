from xarray import register_dataarray_accessor


@register_dataarray_accessor("area")
class SeafluxUtils:
    def __init__(self, data):
        self._obj = data

    def __call__(self, lat_name="lat", lon_name="lon"):
        """
        Returns the area of the grid cells if lat and lon
        are present. You can adjust the names of lat and lon.
        Output units are in m^2
        """
        from .utils import area_grid

        xda = self._obj

        assert lat_name in xda.coords, f"{lat_name} is not in data array"
        assert lon_name in xda.coords, f"{lon_name} is not in data array"

        lat = xda[lat_name].values
        lon = xda[lon_name].values

        area = area_grid(lat, lon, return_dataarray=True)
        area = area.rename(lat=lat_name, lon=lon_name)

        return area
