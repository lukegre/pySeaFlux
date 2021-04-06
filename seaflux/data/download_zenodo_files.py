"""
Module to directly download the latest version of SeaFlux data from Zenodo
"""


from . import utils
from .zenodo_config import zenodo_url


def download_dict_urls(dest, urls_dict):
    """Downloads urls and assigns them in a dictionary

    Note that URLS have to be HTTP/S links that are not protected by passwords

    Args:
        dest (str): path to where files will be downloaded
        urls_dict (dict): keys are names of entries, values are urls

    Returns:
        dict: keys are the input keys and values are the destination paths
    """
    import logging
    import os

    import pooch

    max_len = max([len(k) for k in urls_dict])
    formatted = "{{k: <{max_len}}} {{v}}".format(max_len=max_len)
    pretty_urls = "\n".join([formatted.format(k=k, v=v) for k, v in urls_dict.items()])
    logging.info(f"The following data will be downloaded:\n{pretty_urls}")

    out = {}
    for key, url in urls_dict.items():
        out[key] = pooch.retrieve(
            url=url,
            known_hash=None,
            fname=url.split("/")[-1].split("?")[0],
            path=os.path.expanduser(dest),
            downloader=pooch.HTTPDownloader(progressbar=True),
        )
    return out


def fetch_data(dest_path="~/Downloads/", keep_attrs=False, **urls):
    """Fetches SeaFlux netCDF data from given URLs

    Downloads data and stores at the specified destination path.
    Note that the key, value pairs are variable names and corresponding
    netCDF files that will be downloaded. All netCDF variables must
    be mergable, i.e., must contain the same dimensions and sizes.

    Args:
        dest_path (str): the destination where files will be stored
        urls (key=value): netCDF_fname and variable_name combinations. values
        can also be set to None for flexibility. These keys will not be downloaded

    Returns:
        xr.Dataset: a merged dataset of the given variables and metadata

    """
    import xarray as xr

    from pandas import Timestamp

    from . import utils

    keep = [k for k in urls if isinstance(urls[k], str)]
    urls = {k: urls[k] for k in keep}

    fnames = download_dict_urls(dest_path, urls)

    ds = xr.Dataset()
    for key, fname in fnames.items():
        xds = xr.open_mfdataset(
            fname,
            preprocess=utils.preprocess(),
            chunks={"time": 12, "lat": 90, "lon": 180},
        )
        ds[key] = xds[key]
        ds[key].attrs["history"] = xds.attrs.get("history", "")

    if not keep_attrs:
        date = Timestamp.today()
        ds.attrs = dict(
            source=zenodo_url,
            doi="https://doi.org/10.5281/zenodo.4133802",  # link to all versions
            contact="luke.gregor@usys.ethz.ch",
            citation=(
                "Fay, A. R., Gregor, L., Landschützer, P., McKinley, G. A., "
                "Gruber, N., Gehlen, M., Iida, Y., Laruelle, G. G., Rödenbeck, C., "
                "and Zeng, J.: Harmonization of global surface ocean pCO2 mapped "
                "products and their flux calculations; an improved estimate of the "
                "ocean carbon sink, Earth Syst. Sci. Data Discuss. [preprint], "
                " https://doi.org/10.5194/essd-2021-16, in review, 2021."
            ),
            history=f"Downloaded with SeaFlux code for Python on {date:%Y-%m-%d} ",
        )

    return ds


def flux_calc_data(
    dest_path="~/Downloads/",
    sol_Weiss74=f"{zenodo_url}/files/SeaFluxV2021.01_solWeis74.nc?download=1",
    ice=f"{zenodo_url}/files/SeaFluxV2021.01_icefrac_1988-2018.nc?download=1",
    kw_scaled=f"{zenodo_url}/files/SeaFluxV2021.01_kwScaled16.5cmhr_1988-2018.nc?download=1",
    pCO2atm=f"{zenodo_url}/files/SeaFluxV2021.01_pCO2atm_NOAAmbl%2BERA5mslp_1988-2018.nc?download=1",
):
    """Fetches SeaFlux data for flux calculations, excluding the filling data

    Downloads data and stores at the specified destination path.
    Note that the key, value pairs are variable names and corresponding
    netCDF files that will be downloaded. All netCDF variables must
    be mergable, i.e., must contain the same dimensions and sizes.

    Args:
        dest_path (str): the destination where files will be stored
        sol_Weiss74 (str): url to the solubility dataset
        ice (str): url to the sea ice fraction dataset
        kw_scaled (str): url to the kw values scaled
        pCO2atm (str): url to the atmospheric pCO2 data

    Returns:
        xr.Dataset: a merged dataset of the given variables and metadata. The
        returned data is shaped: time x lat x lon. time is monthly, centered on
        the 15th of each month. lat ranges from -89.5 : 89.5. lon ranges from
        -179.5 : 179.5.
    """
    from ..area import get_area_from_dataset

    urls = utils.get_kwargs()
    urls.pop("dest_path")

    ds = fetch_data(dest_path, **urls)
    ds["area"] = get_area_from_dataset(ds)

    return ds


def scaled_spco2_for_filling(
    dest_path="~/Downloads",
    spco2_clim_scaled=f"{zenodo_url}/files/SeaFluxV2021.01_spCO2filled_1988-2018.nc?download=1",
    scaling_factor=f"{zenodo_url}/files/SeaFluxV2021.01_spCO2filled_1988-2018.nc?download=1",
):
    """Fetches SeaFlux data for filling pseudo-global pCO2 data

    Downloads data and stores at the specified destination path.
    Note that the key, value pairs are variable names and corresponding
    netCDF files that will be downloaded. All netCDF variables must
    be mergable, i.e., must contain the same dimensions and sizes.

    Args:
        dest_path (str): the destination where files will be stored
        spco2_clim_scaled (str): url to the scaled climatology that is used to
        fill the data.
        scaling_factor (str): url to the scaling factor (in the same file as
        spco2_clim_scaled)

    Returns:
        xr.Dataset: the data used to fill pCO2 data products. The product used
        to fill the missing data is from the MPI-ULB-SOMFFN product by
        Landschuetzer et al. (2020). This product is a climatology and is scaled
        to match the ensemble of data products that contains the following
        products: MPI-SOMFFN, CMEMS-FFNN2, NIES-FFNN, JMA-MLR, CSIR-ML6, and
        JENA-MLS. The ratio between the MPI-ULB-SOMFFN climatology and the
        ensemble for the abovementioned products is represented by the scaling
        factor. The returned data is shaped: time x lat x lon. time is monthly,
        centered on the 15th of each month. lat ranges from -89.5 : 89.5. lon
        ranges from -179.5 : 179.5.
    """

    urls = utils.get_kwargs()
    urls.pop("dest_path")

    ds = fetch_data(dest_path, **urls)

    return ds
