"""
Utilities for data set homogenisation, metadata, and plotting
"""


def save_seaflux(xds, dest, variable_name):
    """Helper function to save netCDF files with the correct meta data"""
    from datetime import datetime as dt
    from pathlib import Path as path

    import xarray as xr

    from .config import citation, contact, version, year_range

    if isinstance(xds, xr.DataArray):
        xds = xds.to_dataset(name=variable_name)
    elif not isinstance(xds, xr.Dataset):
        raise TypeError("xds must be a dataset or dataarray")

    assert variable_name in xds
    if "time" in xds:
        xds = xds.sel(time=slice(str(year_range[0]), str(year_range[1])))
        t0, t1 = xds.time.dt.year.values[[0, -1]]
        date_range = f"_{t0}-{t1}"
    else:
        date_range = ""

    sname = f"SeaFlux_{version}_{variable_name}{date_range}.nc"
    full_path = path(dest) / sname

    xds.attrs["citation"] = citation
    xds.attrs["prodct_version"] = version
    xds.attrs["date_created"] = dt.now().strftime("Y%-%m-%d")
    xds.attrs["contact"] = contact

    print(f"[SeaFlux] Saving {variable_name} to {full_path}")

    encoding = {k: dict(zlib=True, complevel=4) for k in xds}
    xds.to_netcdf(full_path, encoding=encoding)

    return full_path


def add_history(xds, message):
    """
    Adds history to xr.Datasets with a time stamp and [SeaFlux] prefix
    """
    import pandas as pd

    time = pd.Timestamp.today().strftime("%Y-%m-%d")
    message = message.replace("'", "").replace('"', "")

    history = xds.attrs.get("history", "").strip()
    history = history if history.endswith(";") or history == "" else history + "; "
    history += f" [SeaFlux @ {time}] {message}; "

    xds = xds.assign_attrs(history=history.strip())

    return xds


def get_kwargs():
    """
    gets keyword=value pairs for a function and returns as a dict
    """
    import inspect

    frame = inspect.currentframe().f_back
    keys, _, _, values = inspect.getargvalues(frame)
    kwargs = {}
    for key in keys:
        if key != "self":
            kwargs[key] = values[key]
    return kwargs


def rename_coords(
    xds,
    lat=["latitude", "Lat", "ylat", "y", "Y"],
    lon=["longitude", "Lon", "xlon", "x", "X"],
    time=["Time", "tmnth", "t", "T"],
    region=["reg", "regions"],
    **kwargs,
):
    """renames coordinates to time, lat, lon from given presents"""

    coord_names = dict(
        lat=lat,
        lon=lon,
        time=time,
        region=region,
        **kwargs,
    )

    rename = {}
    for key in xds.coords:
        for new, old in coord_names.items():
            if key in old:
                rename[key] = new

    xds = xds.rename(**rename)

    if rename != {}:
        xds = add_history(xds, f"renamed variables to match protocol {str(rename)}")

    return xds


def order_lon_180(xds):
    """flips longitudes to -180 : 180"""
    if "lon" not in xds.dims:
        return xds

    lon = xds.lon.copy().values
    xds = xds.assign_coords(lon=(lon - 180) % 360 - 180).sortby("lon")

    if any(lon != xds.lon.values):
        xds = add_history(xds, "Flipped longitudes -180:180")

    return xds


def sort_lats(xds):
    """sort lats"""
    if "lat" not in xds.dims:
        return xds

    import numpy as np

    lat = xds.lat.values
    delta = np.diff(lat[[0, -1]])
    if delta < 0:
        xds = xds.sortby("lat")
        xds = add_history(xds, "Sorted lats")

    return xds


def interpolate_coords(xds):
    """interpolates coordinates onto a global grid"""
    import numpy as np

    dims = []
    if "lat" in xds.dims:
        dims += ("lat",)
    if "lon" in xds.dims:
        dims += ("lon",)

    must_interpolate = []
    for dim in dims:
        coords = xds[dim].values
        if ((coords + 0.5) % 1).sum():
            must_interpolate += (dim,)
    if must_interpolate == []:
        return xds

    if "lon" in must_interpolate:
        x = np.arange(-179.5, 180)
        xds = xds.interp(lon=x)
        xds = add_history(xds, "Interpolated `lon` onto grid centers")
    if "lat" in must_interpolate:
        y = np.arange(-89.5, 90)
        xds = xds.interp(lat=y)
        xds = add_history(xds, "Interpolated `lat` onto grid centers")

    return xds


def transpose(xds):
    """transpose data to [other], time, lat, lon"""
    reccap_order = ["time", "lat", "lon"]
    coords = list(xds.dims)

    if reccap_order == coords:
        return xds

    order = []
    for key in reccap_order:
        if key in coords:
            order += (key,)
            coords.remove(key)
    order = coords + order
    xds = xds.transpose(*order)
    xds = add_history(xds, f"Transposed dimensions to ({str(order)[1:-1]})")

    return xds


def center_time_on_15th(xds):
    """center monthly data on the 15th of the month"""
    import numpy as np
    import xarray as xr

    if "time" not in xds.dims:
        return xds

    dt = np.timedelta64(14, "D")
    t = xds.time.values.astype("datetime64[M]") + dt

    xds["time"] = xr.DataArray(data=t, dims=["time"], coords={"time": t})

    xds = add_history(xds, "Time set to 15th of each month")

    return xds


def add_coord_attrs(xds):
    """Ensures that files are Georeferened for Panoply"""
    xds["lat"].attrs = dict(
        standard_name="latitude",
        units="degrees_north",
        axis="Y",
    )
    xds["lon"].attrs = dict(
        standard_name="longitude",
        units="degrees_east",
        axis="X",
    )
    return xds


def preprocess(
    *args,
    rename_coordinates=True,
    sort_latitudes=True,
    center_months=True,
    interpolate_coordinates=True,
    lon_0_180=True,
    transpose_dims=True,
    add_coord_attributes=True,
):
    """
    Args are additional functions that are applied in sequence
    A function generator that can be used as follows:
        reccap2_formatted_xds = xr.open_mfdataset(
            fname, decode_times=False,
            preprocess=data.preprocess(decode_times=True))

    There are several options that can be switched on or off with boolean.
    All changes are documented and appended to xds.attrs.history
    """

    def dataprep(ds):
        "wrapped function"

        for func in args:
            ds = func(ds)

        if rename_coordinates:
            ds = rename_coords(ds)
        if transpose_dims:
            ds = transpose(ds)
        if sort_latitudes:
            ds = sort_lats(ds)
        if lon_0_180:
            ds = order_lon_180(ds)
        if interpolate_coordinates:
            ds = interpolate_coords(ds)
        if center_months:
            ds = center_time_on_15th(ds)
        if add_coord_attributes:
            ds = add_coord_attrs(ds)

        return ds

    return dataprep


def subplot_map(pos=111, proj=None, round=True, land_color="w", **kwargs):
    """
    Makes an axes object with a cartopy projection for the current figure

    Parameters
    ----------
    pos: int/list [111]
        Either a 3-digit integer or three separate integers
        describing the position of the subplot. If the three
        integers are *nrows*, *ncols*, and *index* in order, the
        subplot will take the *index* position on a grid with *nrows*
        rows and *ncols* columns. *index* starts at 1 in the upper left
        corner and increases to the right.
        *pos* is a three digit integer, where the first digit is the
        number of rows, the second the number of columns, and the third
        the index of the subplot. i.e. fig.add_subplot(235) is the same as
        fig.add_subplot(2, 3, 5). Note that all integers must be less than
        10 for this form to work.
    proj: crs.Projection()
        the cartopy coord reference system object to create the projection.
        Defaults to crs.PlateCarree(central_longitude=0) if not given
    round: bool [True]
        If the projection is stereographic, round will cut the corners and
        make the plot round
    land_color: str ['w']
        the color of the land patches
    **kwargs:
        passed to fig.add_subplot(**kwargs)

    Returns
    -------
    An axes object attached to the current figure. If no figure, a figure
    will be created as is default for plt.gcf()
    """
    import matplotlib.path as mpath
    import matplotlib.pyplot as plt
    import numpy as np

    from cartopy import crs, feature

    proj = crs.PlateCarree(0) if proj is None else proj

    fig = plt.gcf()
    ax = fig.add_subplot(pos, projection=proj, **kwargs)

    # makes maps round
    stereo_maps = (
        crs.Stereographic,
        crs.NorthPolarStereo,
        crs.SouthPolarStereo,
    )
    if isinstance(ax.projection, stereo_maps) & round:

        theta = np.linspace(0, 2 * np.pi, 100)
        center, radius = [0.5, 0.5], 0.475
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)

        ax.set_boundary(circle, transform=ax.transAxes)

    # adds features
    ax.add_feature(feature.LAND, color=land_color, zorder=4)
    ax.add_feature(feature.COASTLINE, lw=0.5, zorder=4)
    ax.outline_patch.set_lw(0.5)
    ax.outline_patch.set_zorder(5)

    return ax


def style_line_subplot(ax, add_zero_line=True, xlim=None, y_range=None):
    """Styles line plots to look ready for publication"""
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np

    mpl.rcParams["axes.spines.right"] = False
    mpl.rcParams["axes.spines.top"] = False

    if ax is None:
        ax = plt.gca()
    plt.sca(ax)

    plt.xticks(rotation=0, ha="center")
    plt.xlabel("")
    plt.title("")

    if ax.get_ylim()[0] < 0 < ax.get_ylim()[1]:
        ax.axhline(0, color="k", ls="--", lw=0.5, zorder=0)

    if xlim is None:
        xlim = ax.get_xlim()
    ax.set_xticks(np.arange("1980", "2020", 5, dtype="datetime64[Y]"))
    ax.set_xticklabels(np.arange(1980, 2020, 5))
    ax.set_xlim(*xlim)

    if y_range is not None:
        center = np.mean(ax.get_ylim())
        upp = center + y_range / 2
        low = center - y_range / 2
        ax.set_ylim(low, upp)

    return ax
