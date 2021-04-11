"""
Module to directly download the latest version of SeaFlux data from Zenodo
"""

from pathlib import Path as path

from . import utils


seaflux_data_path = path(utils.__file__).resolve().parent.absolute()
catalog_name = str(seaflux_data_path / "zenodo.yml")


def get_zenodo_catalog():
    """fetches the default catalog and returns as a dictionary. The dictionary
    is presented as a YAML file if you are using IPython/Jupyter"""
    import fetch_data as fd

    return fd.read_catalog(catalog_name)


_default_catalog = get_zenodo_catalog()
_key = list(_default_catalog.keys())[0]
_default_catalog = _default_catalog[_key]
_dest = _default_catalog.get("dest", None)


def get_seaflux_data(catalog_name=catalog_name, dest=_dest, n_jobs=1, verbose=False):
    """Downloads SeaFlux data from Zenodo using the default yaml file containing
    the paths to the latest SeaFlux data. The data is downloaded and then
    combined. You can create your own yaml file to customise the files you want
    to access."""

    from datetime import datetime as dt

    import fetch_data as fd
    import xarray as xr

    from . import config
    from .utils import preprocess

    cat = fd.read_catalog(catalog_name)
    key = list(cat.keys())[0]
    entry = cat[key]
    entry["dest"] = dest

    flist = fd.download(**entry, n_jobs=n_jobs, verbose=verbose)

    xds = xr.open_mfdataset(flist, preprocess=preprocess())
    xds = xds.assign_attrs(
        product_name="SeaFlux",
        product_version=config.version,
        date_accessed=dt.now().strftime("%Y-%m-%d"),
        contact=config.contact,
        **entry["meta"]
    )

    return xds
