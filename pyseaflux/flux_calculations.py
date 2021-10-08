"""
High level functions
--------------------


"""
from . import solubility as sol
from .area import get_area_from_dataset


def flux_bulk(
    temp_C,
    salt,
    pCO2_sea_uatm,
    pCO2_air_uatm,
    pres_hPa,
    kw_cmhr,
):
    """
    Calculates bulk air-sea CO2 fluxes

    .. math::
        FCO_2 = k_w \\cdot K_0 \\cdot \\Delta pCO_2

    Args:
        temp_C (array): temperature from OISST in degCelcius with an allowable
            range of [-2:45]
        salt (array): salinity from EN4 in PSU. Allowable range [5 : 50]
        pCO2_sea_uatm (array): partial pressure of CO2 in the sea in micro-atmospheres.
            Allowable range is [50 : 1000]
        pCO2_air_uatm (array): partial pressure of CO2 in the air in micro-atmospheres.
            Allowable range is [50 : 1000].
        press_hPa (array): atmospheric pressure in hecto-Pascals with allowable range [500:1500]
        kw_cmhr (array): the gas transfer velocity in (cm/hr). Given the careful choices
            involved in estimating kw, we require the user to explicitly provide kw.
            kw can be calculated with pyseaflux.gas_transfer_velocity.<func>. Things to be
            aware of when calculating kw: wind product and scaling coeffient of gas transfer,
            resolution resampling, and the formulation (i.e. quadratic, cubic).

    Returns:
        array:
            Sea-air CO2 flux where positive is out of the ocean and negative is
            into the ocean. Units are gC.m-2.day-1 (grams Carbon per metre squared
            per day). If the input is an xarray.DataArray, then the output will be
            a data array with fluxes, globally integrated flux, and the area used to
            integrate the fluxes.
    """
    import xarray as xr

    pres_atm = pres_hPa / 1013.25
    temp_K = temp_C + 273.15

    pCO2sea = pCO2_sea_uatm * 1e-6  # to atm
    pCO2air = pCO2_air_uatm * 1e-6

    K0 = sol.solubility_weiss1974(salt, temp_K, pres_atm)

    """unit analysis
    kw = (cm . hr-1) * hr . day-1 . cm-1 . m
    kw = m . day-1   """
    kw = kw_cmhr * (24 / 100)

    # molar mas of carbon in g . mmol-1
    mC = 12.0108 * 1000  # (g . mol-1) / (mmol . mol-1)

    """ unit analysis
    flux = (m . day-1) .  (mol . L-1 . atm-1) . atm . (gC . mmol-1)
    flux = (m . day-1) . (mmol . m-3 . atm-1) . atm . (gC . mmol-1)
    flux = gC . m-2 . day-1   """
    CO2flux_bulk = kw * K0 * (pCO2sea - pCO2air) * mC

    if isinstance(CO2flux_bulk, xr.DataArray):
        area = get_area_from_dataset(CO2flux_bulk)
        ds = xr.Dataset()
        ds["fgco2"] = CO2flux_bulk.assign_attrs(
            units="gC/m2/day",
            description="Air sea CO2 fluxes calculated using the bulk formulation.",
        )
        ds["area"] = area
        ds["fgco2_global"] = (
            (CO2flux_bulk * area * 365)
            .sum(["lat", "lon"])
            .assign_attrs(units="gC/Yr", description="integrated fluxes fgco2 * area")
        )
        return ds
    else:
        return CO2flux_bulk
