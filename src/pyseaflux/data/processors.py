import numpy as np
import xarray as xr
from loguru import logger

# required to make custom processors available in the namespace
from .custom_funcs import *


def run_processors(ds, processes: list):
    global_namespace = globals()
    for p in processes:
        logger.debug(f"Running processor: {p}")
        if isinstance(p, str):
            func = global_namespace[p]
        elif callable(p):
            func = p

        ds = add_history_wrapper(func)(ds)

    return ds


def as_float32(ds):
    return ds.astype("float32")


def subset(ds, **kwargs):
    """
    Subset the dataset
    """
    if kwargs:
        varnames = list(kwargs.keys())
        return ds[varnames]
    else:
        return ds


def rename(ds, **kwargs):
    """
    Rename the variables in the dataset
    """
    ds = subset(ds, **kwargs)
    return ds.rename(kwargs)


def surface(ds: xr.Dataset) -> xr.Dataset:
    """
    Return the surface of the dataset
    """

    check_dim(ds, "depth")

    return ds.sel(depth=0, method="nearest").drop("depth")


def lon_180(ds: xr.Dataset) -> xr.Dataset:
    """
    Convert longitudes to -180:180 format
    """

    check_dim(ds, "lon")
    lon = ds.lon
    ds = ds.assign_coords(lon=_lon180(lon)).sortby("lon")

    return ds


def _lon180(arr):
    return (arr + 180) % 360 - 180


def time_month_start(ds):
    """
    Set the time to the start of the month
    """
    time_m0 = ds.time.astype("datetime64[M]")
    return ds.assign_coords(time=time_m0)


def grid_edge_to_center_lon(ds: xr.Dataset):
    return _grid_edge_to_center(ds, "lon", 180)


def grid_edge_to_center_lat(ds: xr.Dataset):
    return _grid_edge_to_center(ds, "lat", 90)


def _grid_edge_to_center(ds, dim, limit):
    def _extend_coord(c, n_ext=2):
        c = np.array(c)

        dc = np.nanmedian(np.diff(c))

        start = [c[0] - dc * i for i in range(n_ext, 0, -1)]
        center = c
        end = [c[-1] + dc * i for i in range(1, n_ext + 1)]

        out = np.concatenate([start, center, end])

        return out

    def _wrap_coord(c, n_ext=2):
        c = np.array(c)

        start = c[-n_ext:]
        center = c
        end = c[:n_ext]

        out = np.concatenate([start, center, end])

        return out

    coord = ds[dim]

    name = coord.name

    coord_ext_edge_lbl = _extend_coord(coord)
    coord_ext_edge_val = _wrap_coord(coord)

    coord_cntr = np.convolve(coord_ext_edge_lbl, np.ones(2) / 2, mode="valid")
    coord_cntr = coord_cntr[(coord_cntr > -limit) & (coord_cntr < limit)]

    ds_ext_edge = ds.sel(**{name: coord_ext_edge_val}).assign_coords(
        **{name: coord_ext_edge_lbl}
    )

    ds_cntr = ds_ext_edge.interp(**{name: coord_cntr})

    return ds_cntr


def sort_lat(ds: xr.Dataset) -> xr.Dataset:
    """
    Sort the dataset by latitude
    """

    check_dim(ds, "lat")

    return ds.sortby("lat")


def check_dim(ds, dim):
    dims = list(ds.dims)
    if dim not in dims:
        raise ValueError(f"{dim} dimension not found in dataset with dims {dims}")


def add_history_wrapper(func):
    import pandas as pd

    func_path = func.__module__
    func_name = func.__name__
    func_call = f"{func_path}.{func_name}"

    def wrapper(ds, *args, **kwargs):
        time = pd.Timestamp.today().strftime("%Y-%m-%d %H:%M:%S")
        new_history = f" [pySeaFlux @ {time}] {func_call}"
        old_history = ds.attrs.get("history", "")
        list_history = old_history.split("; ") + [new_history]
        list_history = [h for h in list_history if h]

        ds = func(ds, *args, **kwargs)

        ds.attrs["history"] = ";\n".join(list_history)

        return ds

    return wrapper


def resample_to_monthly(ds):
    def check_month_copmlete(ds):
        """
        Ensures that the month has all days of the given month
        """

        time = ds.time
        month = time.dt.month
        months = np.unique(month)

        n_months = months.size
        n_days_in_mon = time.dt.days_in_month[0].item()
        n_days_in_ds = np.unique(time.dt.dayofyear).size

        single_month = n_months == 1
        correct_days = n_days_in_mon == n_days_in_ds

        if not single_month:
            raise ValueError("Dataset has more than one month: ", months)
        if not correct_days:
            current = time.to_index()[0]
            raise ValueError(
                f"{current:%Y-%m} is not complete with {n_days_in_ds} days instead of {n_days_in_mon}"
            )

    groups = ds.time.groupby("time.month")
    for month, group in groups:
        check_month_copmlete(group)

    ds = ds.resample(time="1MS").mean()

    return ds
