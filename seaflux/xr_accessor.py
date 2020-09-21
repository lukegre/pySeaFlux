from xarray import register_dataarray_accessor


@register_dataarray_accessor('area')
class SeafluxUtils:
    def __init__(self, data):
        self._obj = data

    def __call__(self, lat_name='lat', lon_name='lon'):
        from . utils import area_grid

        xda = self._obj

        assert lat_name in xda.coords, f'{lat_name} is not in data array'
        assert lon_name in xda.coords, f'{lon_name} is not in data array'

        lat = xda[lat_name].values
        lon = xda[lon_name].values

        area = area_grid(lat, lon, return_dataarray=True)

        return area
