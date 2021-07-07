"""
Atmospheric pCO2
----------------
This module is not imported by default and has to be manually imported.
This is because it is used only in production of the SeaFlux atmospheric dataset

Contains functions for the full process:
1. Download the NOAA marine boundary layer product (done)
2. Interpolate the product onto a standard grid (done)
3. Download MSLP (ERA5), SST (AVHRR), salinity (EN4) for ATM calculations (done)
4. Convert atmospheric xCO2 to pCO2 with the data
5. Interpolate the final year (2020)

Note that there are no tests for this code, so it is very likely to break
"""

from pathlib import Path as path


base = str(path(__file__).resolve().parent.parent)


def main(
    noaa_mbl_url,
    download_dest="../data/raw/",
    aux_catalog_name="../data/aux_data.yml",
    processed_dest="../data/processed/",
    output_dest="../data/output/",
):
    """to be called when creating the atmospheric pCO2"""
    import xarray as xr

    from fetch_data import read_catalog
    from pandas import Timestamp

    from .aux_vars import download_era5_slp, download_salinity, download_sst_ice
    from .utils import center_time_on_15th, preprocess, save_seaflux

    if path(output_dest).is_file():
        return output_dest

    cat = read_catalog(aux_catalog_name)

    salt = download_salinity(cat["en4_g10"], f"{processed_dest}/en4_salt_temp.nc")
    temp = download_sst_ice(cat["oisst_v2"], f"{processed_dest}/noaa_oisst_sst_icec.nc")
    pres = download_era5_slp(
        download_dest=cat["era5_mslp"]["dest"],
        process_dest=f"{processed_dest}/era5_mslp_monthly.nc",
    )

    ds = xr.merge(
        [
            xr.open_dataset(salt)["salinity"].rename("saltPSU"),
            xr.open_dataset(temp)["sst"].rename("tempC"),
            xr.open_dataset(pres)["sp"].rename("presPa"),
        ]
    )

    noaa_mbl_xco2 = (
        download_noaa_mbl(
            noaa_mbl_url,
            download_dest=f"{download_dest}/co2_GHGreference_surface.txt",
            target_lat=ds.lat.values,
            target_lon=ds.lon.values,
        )
        .resample(time="1MS")
        .mean()
    )

    t0, t1 = ds.time.values[[0, -1]]
    noaa_mbl_xco2 = center_time_on_15th(noaa_mbl_xco2).sel(time=slice(t0, t1))

    t0, t1 = noaa_mbl_xco2.time.values[[0, -1]]
    ds = ds.sel(time=slice(t0, t1))

    atm_pco2 = atm_xCO2_to_pCO2(
        noaa_mbl_xco2, ds.presPa.where(ds.tempC.notnull()) / 100, ds.tempC, ds.saltPSU
    )

    atm_pco2 = preprocess()(
        xr.DataArray(
            data=atm_pco2,
            dims=ds.tempC.dims,
            coords=ds.tempC.coords,
            name="pco2atm",
            attrs=dict(
                long_name=(
                    "partial_pressure_of_carbon_dioxide_in_the_marine_boundary_layer"
                ),
                short_name="pco2atm",
                units="uatm",
                description=(
                    "Atmospheric pCO2 for the marine boundary layer is calculated "
                    "from the NOAAs marine boundary layer pCO2 with: xCO2 * (Patm "
                    "- pH2O). Where pH2O is calculated using vapour pressure from "
                    "Dickson et al. (2007)"
                ),
                history=(
                    getattr(noaa_mbl_xco2, "history", "").strip(";") + ";\n"
                    f"[SeaFlux @ {Timestamp.today():%Y-%m-%d}] "
                    f"pCO2 calculated from xCO2 * (Patm - pH2O), where "
                    f"pH2O is calculated with Dickson et al. (2007)"
                ),
                citation=(
                    "Ed Dlugokencky and Pieter Tans, NOAA/ESRL "
                    "(www.esrl.noaa.gov/gmd/ccgg/trends/)"
                ),
            ),
        )
    )

    variable = "pco2atm"
    pco2atm = interpolate_year(atm_pco2).to_dataset(name=variable)

    sname = save_seaflux(pco2atm, output_dest, variable)

    return sname


def atm_xCO2_to_pCO2(xCO2_ppm, slp_hPa, tempSW_C, salt):
    """
    Convert atmospheric xCO2 to pCO2 with correction for water vapour pressure
        pCO2atm = xCO2atm * (Press - pH2O)

    Args:
        xCO2_ppm (array): atmospheric, or marine boundary layer mole fraction of CO2 (NOAA MBL)
        slp_hPa (array): sea water temperature in degrees C (ERA5 recommended)
        tempSW_C (array): atmospheric pressure in hecto Pascal (NOAA AVHRR OISSTv2 recommended)
        salt (array): sea surface salinity in PSU (EN4 salinity)

    Returns:
        array: note that output will be an np.ndarray regardless of input
    """
    from numpy import array

    from .. import check_units as check
    from .. import vapour_pressure as vapress

    print("[SeaFlux] Converting xCO2 to pCO2")
    xCO2 = array(xCO2_ppm)
    # check units and mask where outsider of range
    Tsw = check.temp_K(tempSW_C + 273.15)
    Ssw = check.salt(salt)
    Patm = check.pres_atm(slp_hPa / 1013.25)

    pH2O = vapress.dickson2007(Ssw, Tsw)

    pCO2atm = xCO2 * (Patm - pH2O)

    return pCO2atm


def read_noaa_mbl_url(noaa_mbl_url, dest):
    """Downloads url and reads in the MBL surface file

    Args:
        noaa_mbl_url (str): the address for the noaa surface file
        dest (str): the destination to which the raw file will be saved

    Returns:
        pd.Series: multindexed series of xCO2 with (time, lat) as coords.
    """
    import re

    from pathlib import Path

    import numpy as np
    import pandas as pd
    import pooch

    # save to temporary location with pooch
    print(
        f"[SeaFlux] Downloading {noaa_mbl_url} to {dest} and reading in as pd.DataFrame"
    )

    dest = Path(dest)
    fname = pooch.retrieve(
        url=noaa_mbl_url,
        known_hash=None,
        path=str(dest.parent),
        fname=str(dest.name),
    )

    # find start line
    is_mbl_surface = False
    for start_line, line in enumerate(open(fname)):
        if re.findall("MBL.*SURFACE", line):
            is_mbl_surface = True
        if not line.startswith("#"):
            break
    if not is_mbl_surface:
        raise Exception(
            "The file at the provided url is not an MBL SURFACE file. "
            "Please check that you have provided the surface url. "
        )

    # read fixed width file CO2
    df = pd.read_fwf(fname, skiprows=start_line, header=None, index_col=0)
    df.index.name = "date"
    # every second line is uncertainty
    df = df.iloc[:, ::2]
    # latitude is given as sin(lat)
    df.columns = np.rad2deg(np.arcsin(np.linspace(-1, 1, 41)))

    # resolve time properly
    year = (df.index.values - (df.index.values % 1)).astype(int)
    day_of_year = ((df.index.values - year) * 365 + 1).astype(int)
    date_strings = ["{}-{:03d}".format(*a) for a in zip(year, day_of_year)]
    date = pd.to_datetime(date_strings, format="%Y-%j")
    df = df.set_index(date)
    df = df.iloc[:-1]  # remove the last value that is for 2020-01-01

    # renaming indexes (have to stack for that)
    df = df.stack()
    index = df.index.set_names(["time", "lat"])
    df = df.set_axis(index)

    df.source = noaa_mbl_url

    return df


def download_noaa_mbl(
    noaa_mbl_url,
    download_dest="../data/raw/co2_GHGreference_surface.txt",
    target_lat=None,
    target_lon=None,
    interp_method="linear",
):
    """
    Downloads the NOAA marine boundary layer xCO2 and grids it
    to a defined grid if target lat and lon provided.

    Args:
        noaa_mbl_url (str): surface xCO2 file from https://www.esrl.noaa.gov/gmd/ccgg/mbl/index.html
        target_lat (None, array_like): if None, the default lats will be returned,
        if array-like then will return xCO2 interpolated onto the given latitudes
        target_lon (None, array-like): if None, data will not be broadcast (expanded) along
        latitudes, if array-like, then will broadcast to those longitudes
        inter_method (str): can be linear or nearest

    Returns:
        xr.DataArray: the data array contains xCO2 interpolated onto the given lats and lons
    """
    import numpy as np
    import xarray as xr

    from pandas import Timestamp

    history = (
        f"[SeaFlux@{Timestamp.today():%Y-%m-%dT%H:%M}]: "
        f"downloaded NOAA MBL data from {noaa_mbl_url}, "
    )

    df = read_noaa_mbl_url(noaa_mbl_url, download_dest)

    print("[SeaFlux] Converting pd.DataFrame to xr.DataArray")
    xda = df.to_xarray()

    if target_lat is not None:
        history += f"latitude interpolated with {interp_method}, "
        xda = xda.interp(lat=target_lat, method=interp_method)

    if target_lon is not None:
        history += "longitude broadcast"
        lon = xr.DataArray(np.ones_like(target_lon), dims=["lon"], coords=[target_lon])
        xda = xda * lon

    xda.attrs = dict(
        units="ppm",
        product="NOAA Greenhouse Gas Marine Boundary Layer Reference",
        history=history,
        source="https://www.esrl.noaa.gov/gmd/ccgg/mbl/index.html",
        description=(
            "mole fraction of CO2 for the marine boundary layer varying by "
            "latitude and time. Note that values are constant along "
            "longitudes. "
        ),
    )

    return xda


def interpolate_year(co2_dataarray):
    """Interpolates atmospheric pCO2 based on average increases

    The data is not truly interpolated, but rather follows the
    method shown below.

    1. [anom] calculate climatological anomaly (clim - mean)
    2. [diff] calculate average difference from DEC to JAN
    3. [last] get last value of last full year of CO2

    Args:
        co2_dataarray (xr.DataArray): a data array of monthly pCO2atm
        year (int, float): the year to be interpolated, typically, the last

    Returns:
        xr.DataArray: the input array + 'interpolated' data with added meta
    """
    import pandas as pd
    import xarray as xr

    print("[SeaFlux] Interpolating missing year from previous years")

    co2 = co2_dataarray

    t0, t1 = [
        str(s)
        for s in (
            co2.groupby("time.year")
            .count("time")
            .where(lambda a: a != 0, drop=True)
            .year.values[[0, -1]]
            .astype(str)
        )
    ]
    co2 = co2.sel(time=slice(t0, t1))

    t0 = int(t1) + 1
    t1 = t0 + 1
    time = pd.date_range(f"{t0}", f"{t1}", freq="1MS", closed="left")
    time += pd.Timedelta(14, "D")

    anom = co2.groupby("time.month").mean("time") - co2.mean("time")
    diff = co2.diff("time")[::12].mean("time")
    last = co2[-1].drop("time")

    intp = (
        (last + diff + anom)
        .rename(month="time")
        .assign_coords(time=time)
        .transpose("time", "lat", "lon")
    )

    today = pd.Timestamp.today()
    history = (
        f"[SeaFlux @ {today:%Y-%m-%d}] The year {t0} has been "
        "interpolated by adding the climatological monthly values "
        "to the last value of the NOAA MBL time series. An offset "
        "between JAN and DEC has also been added to the "
        "climatological anomaly."
    )

    history = "; ".join([co2.attrs.get("history", ""), history])
    co2 = xr.concat([co2, intp], "time").assign_attrs(history=history)

    return co2
