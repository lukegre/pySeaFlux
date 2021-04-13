from pathlib import Path as path


base = str(path(__file__).resolve().parent.parent)


def area(aux_catalog_fname, dest="../data/output/"):
    """Computes the area of the SeaFlux grid cells"""
    import xarray as xr

    from fetch_data import read_catalog

    from ..area import get_area_from_dataset
    from .utils import save_seaflux

    cat = read_catalog(aux_catalog_fname)

    ds = xr.open_mfdataset(download_sst_ice(cat["oisst_v2"])).sst.rename("temp")
    area = get_area_from_dataset(ds)

    sname = save_seaflux(area, dest, "area")
    return sname


def solubility(aux_catalog_fname, dest="../data/output/"):
    """Computes SeaFlux solubility from SST, Salt and Pres"""
    import xarray as xr

    from fetch_data import read_catalog

    from ..solubility import solubility_weiss1974
    from .utils import save_seaflux

    cat = read_catalog(aux_catalog_fname)

    print("[SeaFlux] fetching SST, Salinity, and sea-level pressure")
    ds = xr.merge(
        [
            xr.open_dataset(download_sst_ice(cat["oisst_v2"])).sst.rename("temp")
            + 273.15,
            xr.open_dataset(download_salinity(cat["en4_g10"])).salinity.rename("salt"),
            xr.open_dataset(
                download_era5_slp(download_dest=cat["era5_mslp"]['dest'])
            ).sp.rename("mslp")
            / 101325,
        ]
    )

    time_mask = ds.to_array("tmp").notnull().all("tmp").any(["lat", "lon"])
    ds = ds.where(time_mask, drop=True)

    # unit analysis
    # mol / L / atm --> mol / m3 / uatm
    # mol . L-1 . atm-1 * (1e3L . m-3) * (1e-6 atm . uatm-1) = * 1e-3
    print("[SeaFlux] calculating solubility using Weiss (1974)")
    arr = solubility_weiss1974(ds.salt, ds.temp, press_atm=ds.mslp) * 1e-3

    sol = xr.DataArray(
        data=arr,
        coords=ds.temp.coords,
        dims=ds.temp.dims,
        attrs={
            "description": "CO2 solubility in seawater using the formulation of Weiss 1974",
            "units": "mol/m3/uatm",
            "long_name": "CO2 solubility in seawater",
        },
    )

    sname = save_seaflux(sol, dest, "sol")

    return sname


def sea_ice_cover(aux_catalog_fname, dest="../data/output/"):
    """Calculates SeaFlux sea ice cover as a fraction"""
    import xarray as xr

    from fetch_data import read_catalog

    from .utils import save_seaflux

    cat = read_catalog(aux_catalog_fname)

    fname = download_sst_ice(cat["oisst_v2"])
    variable = "ice"
    ice = xr.open_mfdataset(fname)["icec"].rename(variable) / 100
    ice = ice.sel(time=slice("1982", None))

    sname = save_seaflux(ice, dest, variable)

    return sname


def download_era5_slp(
    year=list(range(1982, 2021)),
    download_dest="../data/raw/ERA5_mslp/",
    process_dest="../data/processed/era5_mslp_monthly.nc",
    **kwargs,
):
    """
    Shortcut for fetching era5 data. requires the `cdsapi` to be
    correctly set up (~/.cdsapi)
    Uses data from `reanalysis-era5-single-levels`
    Fetches data as monthly files to keep requests to a reasonable size
    """
    import logging
    import os

    from pathlib import Path

    import cdsapi
    import xarray as xr

    from joblib import Parallel, delayed
    from numpy import ndarray

    from .utils import preprocess

    if path(process_dest).is_file():
        return process_dest

    if isinstance(year, (list, tuple, ndarray)):
        logging.info(f"downloading to: {download_dest}")
        inputs = [dict(year=y, download_dest=download_dest) for y in year]
        flist = Parallel(n_jobs=8)(
            delayed(download_era5_slp)(**input_dict) for input_dict in inputs
        )
        ds = xr.open_mfdataset(flist, preprocess=preprocess())
        ds.to_netcdf(
            process_dest, encoding={k: dict(zlib=True, complevel=4) for k in ds}
        )
        return ds

    year = str(year)
    sname = os.path.join(download_dest, f"ERA5_surfpress_monthly_{year}.nc")

    if os.path.isfile(sname):
        return sname

    Path(sname).parent.mkdir(exist_ok=True, parents=True)
    cds_client = cdsapi.Client()
    cds_client.retrieve(
        "reanalysis-era5-single-levels-monthly-means",
        {
            "product_type": "monthly_averaged_reanalysis",
            "format": "netcdf",
            "variable": "surface_pressure",
            "year": year,
            "month": [
                "01",
                "02",
                "03",
                "04",
                "05",
                "06",
                "07",
                "08",
                "09",
                "10",
                "11",
                "12",
            ],
            "time": ["00:00"],
        },
        sname,
    )

    return sname


def download_salinity(
    catalog_entry, verbose=True, 
    process_dest="../data/processed/en4_salt_temp.nc",
):
    """Downloads salinity from MetOffice for 1982 until today"""
    import xarray as xr

    from fetch_data import download

    from .utils import preprocess

    if path(process_dest).is_file():
        return process_dest

    flist = download(**catalog_entry, verbose=verbose)
    ds = preprocess()(
        xr.open_mfdataset(paths=flist)[["salinity"]]
        .sel(depth=0, method="nearest")
        .drop("depth")
    )

    encode = {k: dict(zlib=True, complevel=4) for k in ds}
    ds.load().to_netcdf(process_dest, encoding=encode)

    return process_dest


def download_sst_ice(
    catalog_entry,
    process_dest="../data/processed/noaa_oisst_sst_icec.nc",
):
    """Downloads OISSTv2 data from NOAA"""
    import xarray as xr

    from fetch_data import download

    from .utils import preprocess

    if path(process_dest).is_file():
        return process_dest

    flist = download(**catalog_entry)

    ds = (
        xr.open_mfdataset(paths=flist, preprocess=preprocess())
        .where(lambda a: a.icec.notnull())
        .drop("time_bnds")
    )

    ds.to_netcdf(process_dest, encoding={k: dict(zlib=True, complevel=4) for k in ds})

    return process_dest
