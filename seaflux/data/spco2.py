"""
Create data for SeaFlux filling data for surface pCO2 products
"""

import numpy as np
import pandas as pd
import xarray as xr

from fetch_data import download

from .utils import preprocess


def main(catalog_fname="../data/spco2_data.yml", dest="../data/output/", verbose=True):
    """
    Runs functions to compute the filling data for SeaFlux
    """
    from .utils import save_seaflux

    ds = SOCOMensemble(catalog_fname, verbose=verbose).data

    # remove JMA for a longer time series (1985-2019)
    drop = ["seamask", "MPI_ULB_SOMFFN", "month", "JMA_MLR"]
    pco2 = ds.drop(drop).to_array("product")
    clim = ds.MPI_ULB_SOMFFN.sel(month=pco2.time.dt.month).drop("month")

    mask = pco2.notnull().all("product")
    pco2 = pco2.where(mask).dropna("time", "all")

    variable = "spco2filler"
    products = pco2.product.values.astype(str)
    reference_info = f"Ensemble of spco2 products: {', '.join(products)}"

    filler = make_filling_data(pco2, clim, "SOCOM ensemble", varname=variable)
    filler.attrs["climatology_data"] = "MPI-ULB-SOMFFN"
    filler.attrs["reference_data"] = reference_info

    sname = save_seaflux(filler, dest, variable)

    return sname


def calculate_fluxes(spco2, apco2, tempC, salt, preshPa, wind, ice):
    """
    TODO: might delete
    """
    import seaflux as sf

    fluxes = sf.flux_bulk(
        tempC,
        salt,
        spco2,
        apco2,
        preshPa,
        wind,
        kw_func=sf.gas_transfer_CO2.k_Wa92,
        kw_scaling=16,
    )

    xds = xr.Dataset()
    xds["fgco2"] = fluxes * (1 - ice)
    xds["apco2"] = apco2
    xds["ice"] = ice
    xds["area"] = sf.utils.area_grid()

    xds.fgco2.attrs["units"] = "gC/m2/day"
    xds.fgco2.attrs["long_name"] = "air-sea flux of CO2"
    xds.fgco2.attrs["description"] = (
        "Air-sea co2 flux where upward is positive: \n"
        "FCO2 = K0 * kw * (spco2 - apco2) * ice \n"
        "kw = Wanninkhof (1992) scaled to 16 cm/hr; \n"
        "K0 = Weiss (1974); \n"
        "apco2 = xCO2mbl * (P - pH2O), using pH2O from Dickson et al 2007; \n"
        "spco2 = data products in spco2; \n"
        "ice = (1 - sea ice fraction)\n"
        "ERA5 winds + MSLP, OSTIA SST, EN4 salinity were used as auxiliary data"
    )

    xds.apco2.attrs.pop("citation")
    xds.apco2.attrs["long_name"] = "atmospheric partial pressure of CO2"
    xds.apco2.attrs["units"] = "uatm"
    xds.apco2.attrs["source"] = "https://www.esrl.noaa.gov/gmd/ccgg/mbl/"
    xds.apco2.attrs["description"] = (
        "Atmospheric pCO2 was calculated with:\n"
        "pCO2 = xCO2mbl * (P - pH2O)\n\n"
        "xCO2mbl = NOAAs xCO2 marine boundary layer product (surface) interpolated linearly along latitudes and extrapolated along longitudes;\n"
        "P = atmospheric pressure from ERA5;\n "
        "pH2O = water vapour pressure calculated using Dickson et al. (2007);"
    )

    xds.ice.attrs["long_name"] = "sea ice fraction"
    xds.ice.attrs["units"] = "/100"
    xds.ice.attrs["description"] = "Sea ice fraction from the OSTIAv2 product"

    return xds


def make_filling_data(reference, climatology, reference_name, varname="spco2filler"):
    """
    Calculates the filling data based on a reference dataset and a climatology.

    The reference name will be used for metadata. Varname is to ensure that
    the correct variable name is assigned to the filling output.
    """
    from datetime import datetime as dt

    from ..area import get_area_from_dataset

    print("[SeaFlux] Calculating filler data from scaled climatology")
    ref = reference.dropna("time", "all")
    clim = climatology.dropna("time", "all")
    area = get_area_from_dataset(clim)

    dims = list(set(ref.dims) - set(["time"]))

    factor = (ref / clim).weighted(area).mean(dims)
    filler = clim * factor

    ds = xr.Dataset(
        attrs=dict(
            date=dt.today().strftime("%Y-%m-%d"),
            software="https://github.com/lukegre/SeaFlux",
            contact="Luke Gregor (luke.gregor@usys.ethz.ch)",
            description=(
                "A surface pCO2 product to fill marginal and high latitude "
                "seas for surface ocean pCO2 products that have pseudo-global "
                "coverage. This should be used for global flux integration "
                "calculations to ensure that the area being compared is truly "
                "global. The filling product is based on the MPI-ULB-SOMFFN "
                "product (Landschuetzer et al., 2020)."
            ),
        )
    )

    ds["scaling_factor"] = factor.assign_attrs(
        description=(
            "A scaling factor for surface ocean pCO2 to scale the MPI-ULB-SOMFFN "
            "climatological product. The scaling factor is calculated as the "
            f"ratio between a reference dataset ({reference_name}) and the "
            "climatology: (reference / climatology).mean([lat, lon]), where the "
            "average is weighted by pixel area. For more info see Fay et al. "
            "(2021)."
        )
    )

    ds[varname] = filler.assign_attrs(
        units="uatm",
        description=(
            "Scaled surface ocean pCO2 as estimated using the MPI-ULB-SOMFFN "
            "product. The climatological pco2 is scaled to a reference dataset "
            f"({reference_name}). The scaling is performed with: "
            "reference / climatology. For more details, see Fay et al. (2021)."
        ),
    )

    return ds


class SOCOMensemble:
    """
    Creates an object that makes it easy to access SOCOM data
    """

    def __init__(self, catalog_fname, verbose=True):
        """
        An object that downloads, reads in, homogenizes and combines
        surface ocean pCO2 products.

        The catalog must contain entries that match the functions in the
        SOCOM ensemble object. The catalog contains the urls to where the data
        can be downloaded. Further, the object contains functions are tailormade
        to homogenise the data sets so that they can easily be worked with.

        A default list of ensemble members is stored under ``self.members``
        which excludes the climatology and sea mask. The climatology name is
        stored under ``self.climatology``.

        The data can be accessed under self.data, which will load the dataset
        from memory if loaded, from a file if it exists, process each individual
        product, or download the products. A full pipeline!

        The catalog is not stored on GitHub as it contains passwords. Contact
        Luke for access to this catalog.
        """
        from fetch_data import read_catalog

        self.catalog_fname = catalog_fname
        self.cat = read_catalog(catalog_fname)

        self.members = [
            "jena_mls",
            "mpi_somffn",
            "cmems_ffnn",
            "csir_ml6",
            "not there",
            "nies_fnn",
            "jma_mlr",
        ]
        self.climatology = 'mpi_ulb_somffn'

        missing = [key for key in self.members if key not in self.cat]
        if any(missing):
            for m in missing:
                self.members.remove(m)

        self._data = None
        self.verbose = verbose
        self.aux_catalog_name = '../data/aux_data.yml'

        print("[SeaFlux] Default ensemble members in catalog:", ", ".join(self.members))

    @property
    def data(self):
        """quickly access the data as an object"""
        if self._data is None:
            self._data = self.__call__()
        return self._data

    def __call__(self, dest="../data/processed/spco2_socom_ensemble.nc"):
        """
        Will download all the data, process it and save it to the given destination
        """
        from pathlib import Path as path

        if path(dest).is_file():
            print(f"[SeaFlux] Returning existing file {dest}")
            return xr.open_dataset(dest)

        # mask
        mask = self.get_data("seamask").assign_attrs(long_name="ocean_mask")

        attrs = dict(units="uatm", long_name="surface_ocean_pco2")
        clim = self.get_data(self.climatology).assign_attrs(**attrs)

        ensemble = [self.get_data(k).assign_attrs(**attrs) for k in self.members]

        ds = xr.merge(ensemble).sel(time=slice("1982", None))
        ds = ds.assign_attrs(
            date_created=pd.Timestamp.today().strftime("%Y-%m-%d"),
            created_by="SeaFlux (Luke Gregor)",
            description=(
                "Surface ocean pCO2 from various data-based products. The information for these "
                "products can be found in each of the DataArrays's meta data. These data are "
                "part of the SOCOM ensemble and are updated annually. We limit the data from 1982 to "
                "present. "
            ),
        )
        ds[self.climatology.upper()] = clim
        ds["seamask"] = mask

        encode = {k: dict(zlib=True, complevel=4) for k in ds}
        ds.to_netcdf(dest, encoding=encode)

        return ds

    def get_data(self, name):
        """Helper function to download the data"""
        import logging

        if name not in self.cat:
            raise KeyError(
                f"{name} is not in the given data catalog {self.catalog_fname}"
            )

        func = getattr(self, f"get_{name}", None)
        if func is None:
            raise KeyError(f"get_{name} is not a function in {str(self)}")

        logging.info(f"[SeaFlux] retrieving {name.upper()}")

        kwargs = self.cat[name]
        meta = kwargs.get("meta", {})
        kwargs['verbose'] = self.verbose
        
        url = kwargs.get("url")
        url = url[0] if isinstance(url, list) else url
        
        data = func(kwargs).rename(name.upper()).load()

        data = data.assign_attrs(source=url, **meta)

        return data

    @staticmethod
    def get_mpi_ulb_somffn(entry):
        """processes data"""

        flist = download(**entry)

        xds = xr.open_mfdataset(flist)
        xda = xds.pco2.where(xds.pco2 > 0).coarsen(lat=4, lon=4).mean()
        xda = xda.rename("mpiulb_somffn").rename(time="month")

        pp = preprocess(rename_coordinates=False, center_months=False)

        xda = pp(xda)
        return xda

    @staticmethod
    def get_mpi_somffn(entry):
        """processes data"""

        flist = download(**entry)

        xda = xr.open_mfdataset(flist, drop_variables="date").spco2_raw
        xda = xda.rename("mpi_somffn").assign_attrs(
            units="uatm", source=entry["url"], **entry["meta"]
        )

        xda = preprocess()(xda)

        return xda

    @staticmethod
    def get_jena_mls(entry):
        """processes data"""

        flist = download(**entry)

        xds = xr.open_mfdataset(flist)
        xda = xds.pCO2.resample(mtime="1MS").mean("mtime")

        xda = xda.rename("jena_mls")
        xda = (
            xda.interp(
                lat=np.arange(-89.5, 90), lon=np.arange(-179.5, 180), method="nearest"
            )
            .roll(lon=180, roll_coords=False)
            .interpolate_na("lon", limit=20)
            .roll(lon=-180, roll_coords=False)
            .rename(mtime="time")
            .assign_attrs(
                units="uatm",
                source=entry["url"],
                **entry["meta"],
                history=(
                    "[Seaflux] resampled from daily to monthly and "
                    "interpolated to 1 degree using nearest neighbour interpolation"
                ),
            )
        )

        xda = preprocess()(xda)

        return xda

    @staticmethod
    def get_cmems_ffnn(entry):
        """processes data"""

        flist = download(**entry, n_jobs=8)

        xds = xr.open_mfdataset(flist, combine="nested", concat_dim="time")
        xda = (
            (xds.spco2 * 9.867)
            .assign_coords(longitude=(xds.longitude - 180) % 360 - 180)
            .rename(latitude="lat", longitude="lon")
            .resample(time="1MS")
            .mean()
            .sortby("lon")
            .assign_attrs(units="uatm", source=entry["url"], **entry["meta"])
        )

        xda = preprocess()(xda)

        return xda

    @staticmethod
    def get_nies_fnn(entry):
        """processes data"""
        from warnings import filterwarnings
        
        from fetch_data import read_catalog
        from ..fco2_pco2_conversion import fCO2_to_pCO2
        from .utils import add_history
        from .aux_vars import download_era5_slp, download_sst_ice

        filterwarnings("ignore", category=RuntimeWarning)

        def decode_time(xds):
            """processes data"""
            import pandas as pd

            from datetime_matcher import DatetimeMatcher

            re_date = DatetimeMatcher()

            fname = xds.encoding["source"]
            datetime = re_date.extract_datetime("flux.%Y.ver", fname)
            year = pd.Timestamp(datetime).year

            y0, y1 = str(year), str(year + 1)
            time = pd.date_range(y0, y1, freq="1MS", closed="left")

            xds = xds.rename(month="time").assign_coords(time=time)

            xds = add_history(xds, "decode times manually")

            return xds

        flist = download(**entry)
        xda = xr.open_mfdataset(flist, preprocess=preprocess(decode_time)).fco2
        
        aux_cat = read_catalog('../data/aux_data.yml')
        
        t0, t1 = [str(s) for s in xda.time.values[[0, -1]]]
        sst = xr.open_dataset(download_sst_ice(aux_cat['oisst_v2']))["sst"].sel(time=slice(t0, t1))
        msl = xr.open_dataset(download_era5_slp())["sp"].sel(time=slice(t0, t1)) / 100

        pco2 = xr.DataArray(
            fCO2_to_pCO2(xda, sst, msl),
            coords=xda.coords,
            dims=xda.dims,
            attrs=dict(units="uatm", source=entry["url"], **entry["meta"]),
        )

        pco2 = add_history(
            pco2, "re-shaped data from [year month lat lon] to [time lat lon]."
        )
        pco2 = add_history(
            pco2, "converted fCO2 to pCO2 using OISST v2.1, and ERA5 MSLP"
        )

        return pco2

    @staticmethod
    def get_jma_mlr(entry):
        """processes data"""

        def decode_time(xds):
            """processes data"""
            import pandas as pd

            from seaflux.data.utils import add_history

            time = xds.time
            unit = time.attrs.get("units")
            year = pd.to_datetime(unit.split()[-1]).year

            y0, y1 = str(year), str(year + 1)
            time = pd.date_range(y0, y1, freq="1MS", closed="left")
            xds = xds.assign_coords(time=time)

            xds = add_history(xds, "decode times manually")

            return xds

        flist = download(**entry, n_jobs=8)

        xda = xr.open_mfdataset(
            flist, decode_times=False, preprocess=preprocess(decode_time)
        ).pCO2s
        xda = xda.assign_attrs(units="uatm")

        return xda

    @staticmethod
    def get_csir_ml6(entry):
        """processes data"""

        flist = download(**entry)
        xds = xr.open_mfdataset(flist, preprocess=preprocess())

        xds = xds["spco2"].assign_attrs(
            units="uatm", source=entry["url"], **entry["meta"]
        )

        return xds

    @staticmethod
    def get_seamask(entry):
        """processes data"""

        flist = download(**entry)
        xds = xr.open_mfdataset(flist, preprocess=preprocess())

        xda = xds.seamask.assign_attrs(**entry["meta"])

        return xda
