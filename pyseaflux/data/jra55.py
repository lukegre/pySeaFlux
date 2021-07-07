#!/usr/bin/env python
"""
Fetch JRA-55 data from the UCAR RDAMS server
Script adapted from rdams_client.py
"""
from .utils import parallel


jra_meta = {
    "doi": "https://doi.org/10.2151/jmsj.2015-001",
    "info": "https://jra.kishou.go.jp/JRA-55/leaflet/JRA-55_leaflet_display.pdf",
    "citation": "KOBAYASHI, S., OTA, Y., HARADA, Y., EBITA, A., MORIYA, M., ONODA, H.,  ONOGI, K., KAMAHORI, H., KOBAYASHI, C., ENDO, H., MIYAOKA, K.,  TAKAHASHI, K. (2015). The JRA-55 Reanalysis: General Specifications and  Basic Characteristics. Journal of the Meteorological Society of Japan.  Series II, 93(1), 5-48. https://doi.org/10.2151/jmsj.2015-001",
    "description": "JMA carried out the second reanalysis project (known as the Japanese 55-year  Reanalysis, or JRA-55) using a more sophisticated DA system based on the  operational system as of December 2009, and newly prepared dataset of past observations. The analysis period covers the 55 years from 1958, when regular radiosonde observation began on a global basis.",
}


class RDAMScookies:
    def __init__(self):
        from pathlib import Path as _path

        self.password_file = _path("~/.rdamsrc").expanduser().absolute().resolve()

        self.__version__ = "2.0.1"
        self.__author__ = (
            "Doug Schuster (schuster@ucar.edu), Riley Conroy (rpconroy@ucar.edu)"
        )

    def _get_userinfo(self):
        """Get username and password from the command line."""
        import getpass

        user = input("Enter your RDA username or email: ")
        pasw = getpass.getpass("Enter your RDA password: ")
        self._write_pw_file(user, pasw)
        return (user, pasw)

    def _write_pw_file(self, username, password):
        """Write out file with user information."""
        import codecs

        with open(self.password_file, "w") as file_open:
            npwstring = username + "," + password
            ob_str = codecs.encode(npwstring, "rot_13")
            file_open.write(ob_str)

    def _read_pw_file(self):
        """Read user information from pw file.

        Args:
            pwfile (str): location of password file.

        Returns:
            (tuple): (username,password)
        """
        import codecs

        with open(self.password_file, "r") as f:
            pwstring = codecs.decode(f.read(), "rot_13")
            (username, password) = pwstring.split(",", 2)
        return (username, password)

    def get_authentication(self):
        """Attempts to get authentication.

        Args:
            pwfile (str): location of password file.

        Returns:
            (tuple): username, passord
            (None): If using .netrc file
        """
        import os

        pwfile = self.password_file
        if os.path.isfile(pwfile) and os.path.getsize(pwfile) > 0:
            return self._read_pw_file()
        else:
            return self._get_userinfo()

    def get_cookies(self, username=None, password=None):
        """Authenticates with RDA and returns authentication cookies.

        The user must authenticate with
        authentication cookies per RDA policy.

        Args:
            username (str): RDA username. Typically the user's email.
            password (str): RDA password.

        Returns:
            requests.cookies.RequestsCookieJar: Login request's cookies.
        """
        import requests

        if username is None and password is None:
            username, password = self.get_authentication()

        login_url = "https://rda.ucar.edu/cgi-bin/login"
        values = {"email": username, "passwd": password, "action": "login"}
        ret = requests.post(login_url, data=values)
        if ret.status_code != 200:
            print("Bad Authentication")
            print(ret.text)
            exit(1)
        return ret.cookies


def make_jra_6hrly_urls(
    url_fmt="http://rda.ucar.edu/data/ds628.0/anl_surf/{t0:%Y}/anl_surf.{code_variable}.reg_tl319.{t0:%Y%m%d%H}_{t1:%Y%m%d%H}",
    t0=None,
    t1=None,
    variable_codes=["033_ugrd", "034_vgrd"],
):
    from datetime import datetime

    from pandas import Timedelta, date_range

    today = datetime.now().strftime("%Y-%m-%d")

    period0 = date_range("1982-01-01", "2014-01-01", freq="1AS")
    period1 = date_range("2014-01-01", today, freq="1MS")
    dates = period0.join(period1, how="outer")

    i = dates.notnull()
    if t0 is not None:
        i *= dates >= t0
    if t1 is not None:
        i *= dates <= t1
    dates = dates[i]

    urls = []
    for i in range(dates.size - 1):
        for v in variable_codes:
            t0 = dates[i]
            t1 = dates[i + 1] - Timedelta("6H")
            urls += (url_fmt.format(t0=t0, t1=t1, code_variable=v),)

    return urls


@parallel
def grib_to_netcdf(input_filename, output_filename, load_first=False, overwrite=False):
    """converts a grib file to netcdf and returns the file name

    Args:
        input_filename (str): path to the input file
        output_filename (str): path to the output file
        load_first (bool): will load the netCDF first if True, might be faster
        overwrite (bool): will overwrite files if they exist, otherwise skips

    Returns:
        str: returns the output_filename if successful
    """

    def remove_grib_attrs(attrs):
        return {k.replace("GRIB_", ""): v for k, v in attrs.items()}

    import logging

    from pathlib import Path as path

    from xarray import open_dataset

    output_filename = path(output_filename)

    if output_filename.is_file() and not overwrite:
        return str(output_filename)

    output_filename.parent.mkdir(exist_ok=True, parents=True)

    logging.log(15, f"converting GRIB to netCDF4: {output_filename}")
    xds = open_dataset(input_filename, engine="cfgrib")
    if load_first:
        xds = xds.load()

    drop = [k for k, v in xds.coords.items() if v.size == 1]
    xds = xds.drop(drop)

    for key in xds:
        xds[key].attrs = remove_grib_attrs(xds[key].attrs)
    xds.attrs = remove_grib_attrs(xds.attrs)

    xds.to_netcdf(
        output_filename, encoding={k: {"complevel": 4, "zlib": True} for k in xds}
    )

    return str(output_filename)


def calculate_wind_speed(flist, u="u10", v="v10"):
    from numpy import arange
    from seaflux.data.utils import preprocess
    from xarray import Dataset, open_mfdataset

    prep = preprocess(
        center_months=False, interpolate_coordinates=False, lon_0_180=False
    )

    xds = open_mfdataset(flist, preprocess=prep)

    second_moment = xds["u10"] ** 2 + xds["v10"] ** 2
    wind_speed = Dataset()
    wind_speed["wind_speed_1st"] = second_moment ** 0.5
    wind_speed["wind_speed_2nd"] = second_moment
    wind_speed["wind_speed_3rd"] = wind_speed.wind_speed_1st ** 3

    wind_speed_moments = (
        wind_speed.resample(time="1MS", loffset="14D")
        .mean()  # monthly resolution
        .interp(lat=arange(-89.5, 90), lon=arange(0.5, 361))  # 1deg interp
        .roll(lon=180, roll_coords=False)  # fills the gap in data at 359.5
        .interpolate_na("lon", limit=2)
        .roll(lon=-180, roll_coords=False)
        .sel(lon=slice(0, 360))  # make sure we only have 360 lons
    )

    describe = "Calculated at the native model resolution and then coarsened. "
    wind_speed_moments["wind_speed_1st"].attrs = dict(
        description="JRA 55 wind speed. " + describe,
        units="m.s-1",
        processing=xds.attrs.get("history", "not available"),
    )
    wind_speed_moments["wind_speed_2nd"].attrs = dict(
        description="2nd moment of the JRA 55 wind speed. " + describe,
        units="m2.s-2",
        processing=xds.attrs.get("history", "not available"),
    )
    wind_speed_moments["wind_speed_3rd"].attrs = dict(
        description="3rd moment of the JRA 55 wind speed. " + describe,
        units="m3.s-3",
        processing=xds.attrs.get("history", "not available"),
    )

    return wind_speed_moments


def get_jra55_wind_speed(
    url="leave empty - replaced in function",  # for transparency
    download_dest="../data/raw/jra_55/{file_format}/{year}",
    process_dest="../data/processed/jra55_wind_speed_moments.nc",
    years=range(1982, 2021),
    verbose=False,
    n_jobs=8,
):
    """
    TODO: add readme to netCDF folder.
    TODO: add readme to
    """
    from pathlib import Path as path

    from dask.diagnostics import ProgressBar
    from fetch_data import download
    from fetch_data.core import create_download_readme
    from fetch_data.utils import commong_substring
    from pandas import Timestamp
    from xarray import concat

    years = list(years)
    process_dest = p = path(process_dest)
    process_dest = p.parent / f"{p.stem}_{years[0]}-{years[-1]}{p.suffix}"

    if path(process_dest).is_file():
        return process_dest
    else:
        print(f"File does not exist: {process_dest}")

    cookies = RDAMScookies().get_cookies()
    grib_names = []
    for y in years:
        t0 = Timestamp(f"{y}")
        t1 = Timestamp(f"{y+1}")
        grib_names += download(
            # JRA URLs switch from annual to monthly in 2014
            url=make_jra_6hrly_urls(t0=t0, t1=t1),
            dest=download_dest.format(
                year=y, file_format="grib"
            ),  # store the data per year
            login=dict(cookies=cookies),
            verbose=verbose,
            n_jobs=n_jobs,
            log_name="../downloading.log",
            readme_fname="../README.txt",
            meta=jra_meta,
        )

    # replace the path '/grib/' with netcdf for he conversion
    netcdf_names = [f.replace("/grib/", "/netcdf/") + ".nc" for f in grib_names]
    # the function grib_to_netcdf has been made to run in parallel with the decorator
    flist = grib_to_netcdf(grib_names, netcdf_names, n_jobs=n_jobs)

    jra_meta["processing"] = (
        "Data has been converted from grib file format to netCDF4 using the cfgrib "
        "package. Variables without dimensions have been dropped. "
    )
    jra_meta["variables"] = "u10, v10"
    jra_meta["grib_source"] = download_dest.format(year="YYYY", file_format="grib")
    create_download_readme(
        "README.md",
        url=commong_substring(grib_names) + "...",
        dest=str(path(download_dest.format(year="YYYY", file_format="netcdf")).parent),
        meta=jra_meta,
    )

    # we get the folders for each year
    folders = sorted(list(set([path(f).parent for f in flist])))
    xds = []
    for folder in folders:
        # list nc files - assumes u10 and v10 in the folder
        ylist = list(folder.glob("*.nc"))
        xds += (calculate_wind_speed(ylist),)

    with ProgressBar():
        wind_speed = concat(xds, "time").load()

    process_dest = path(process_dest)
    process_dest.parent.mkdir(exist_ok=True, parents=True)

    jra_meta["processing"] = (
        "Data has been converted from grib file format to netCDF4)."
        "u10 and v10 data was loaded and the wind_speed was calculated with "
        "(u10^2 + v10^2)^0.5. The first, second, and third moments "
        "(wind_speed^n) were calculated from wind_speed. Note that these "
        "variables were calculated at the model resolution and then scaled "
        "up to monthly by 1 degree to preserve the variability that would "
        "otherwise be lost in the squared function. "
    )
    jra_meta["variables"] = "wind_speed, wind_speed^2, wind_speed^3"
    jra_meta["netcdf_source"] = download_dest.format(year="YYYY", file_format="netcdf")

    wind_speed.attrs = jra_meta
    wind_speed.to_netcdf(
        str(process_dest),
        encoding={k: {"complevel": 4, "zlib": True} for k in wind_speed},
    )

    return str(process_dest)
