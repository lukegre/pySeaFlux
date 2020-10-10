"""
Common tools used to correct and calculate pCO2, fCO2 and fluxes

Please cite the use of this code as Gregor et al. (2019) but make sure to
reference Woolf et al. (2016) if the woolf flux method is used. Same applies
for all other functions in this script (Weiss 1974, Dickson et al. 2007 ...)

https://doi.org/10.5194/gmd-2019-46
"""
import warnings

from . import check_units as check
from . import fco2_pco2_conversion as f2p
from . import gas_transfer_velocity
from . import solubility as sol
from .utils import preserve_xda


def flux_woolf2016_rapid(
    temp_bulk_C,
    salt_bulk,
    pCO2_bulk_uatm,
    pCO2_air_uatm,
    press_hPa,
    wind_ms,
    kw_func=gas_transfer_velocity.k_Ni00,
    cool_skin_bias=-0.14,
    salty_skin_bias=0.1,
):
    """
    Calculates air sea CO2 fluxes using the RAPID model as defined by Woolf et
    al. (2016), where the concentration of CO2 in the skin and foundation
    layers are used to calculate the fluxes rather than delta pCO2 (latter is
    called bulk flux).

    We calculate the skin temperature and salinity using a cool and salty skin
    bias as defined in Woolf et al. (2016). The defaults are 0.14 degC and
    0.1 PSU as taken from FluxEngine.

    **Assumptions: ** This function is set up to use AVHRR only OISST which
    reports temperatures at 1m depth based on a buoy correction (Banzon et al.
    2016). We make the assumption that this bulk temperature is equivalent to
    foundation temperature (where nighttime and daytime temperatures are the
    same). We also assume that EN4 salinity is foundation salinity (this is
    probably more accurate than the first assumtion). Lastly we assume that the
    ML estimated fCO2 is bulk fCO2 – we use bulk variable inputs (SSS and SST).

    Parameters
    ----------
    temp_bulk_C : np.array
        temperature from OISST in deg Celcius with an allowable range of
        [-2 : 45]
    salt_bulk : np.array
        salinity from EN4 in PSU. Allowable range [5 : 50]
    pCO2_bulk_uatm : np.array
        partial pressure of CO2 in the sea in micro-atmospheres, assuming that
        it was measured/predicted at the same level as the temperature and
        salinity (See our assumptions above). Allowable range is [50 : 1000]
    pCO2_air_uatm : np.array
        partial pressure of CO2 in the air in micro-atmospheres. Allowable
        range is [50:1000].
    press_hPa : np.array
        atmospheric pressure in hecto-Pascals with an allowable range of
        [500 : 1500] hPa
    wind_ms : np.array
        wind speed in metres per second with an allowable range of [0 : 40]
    kw_func : callable
        a function that returns the gas transfer velocity in cm/hr. The default
        is the gas transfer volicty as calculated by Ho et al. (2006). This
        is the prefered method of Goddijn-Murphy et al. (2016). Other functions
        are available in the `gas_transfer` class. If you'd like to use your
        own inputs must be wind speed (m/s) and temperature (degC) and output
        must be cm/hr
    cool_skin_bias : float
        The temperature difference between the foundation/bulk temperature and
        the skin temperature as suggested by Wolf et al. (2016). The default is
        0.14 degC where this will be subtracted from the bulk temperature, i.e.
        the surface is cooler due to the cooling effect of winds.
    salty_skin_bias : float
        The salinity difference between the foundation and skin layers. This is
        driven by evaporation and defaults to 0.1 (will be added to salinity).

    Reurns
    ------
    FCO2 : np.array
        Sea-air CO2 flux where positive is out of the ocean and negative is
        into the ocean. Units are gC.m-2.day-1 (grams Carbon per metre squared
        per day)
    """
    from numpy import array
    from xarray import DataArray

    warnings.warn("This function has not been tested yet")

    if isinstance(pCO2_bulk_uatm, DataArray):
        var = pCO2_bulk_uatm.copy()  # attribute preservation
    else:
        var = None

    press_atm = array(press_hPa) / 1013.25

    SSTfnd_C = array(temp_bulk_C)
    SSTskn_C = SSTfnd_C - cool_skin_bias  # from default FluxEngine config
    SSTfnd_K = SSTfnd_C + 273.15
    SSTskn_K = SSTskn_C + 273.15
    SSTdelta = SSTfnd_C - SSTskn_C

    SSSfnd = array(salt_bulk)
    SSSskn = SSSfnd + salty_skin_bias  # from default FluxEngine config

    pCO2sea = array(pCO2_bulk_uatm) * 1e-6  # to atm
    pCO2air = array(pCO2_air_uatm) * 1e-6

    # checking units
    press_atm = check.pres_atm(press_atm)
    SSTfnd_K = check.temp_K(SSTfnd_K)
    SSSfnd = check.salt(SSSfnd)
    pCO2sea = check.CO2_mol(pCO2sea)
    pCO2air = check.CO2_mol(pCO2air)
    wind_ms = check.wind_ms(wind_ms)

    fCO2sea = pCO2sea * f2p.virial_coeff(SSTfnd_K, press_atm)
    fCO2air = pCO2air * f2p.virial_coeff(SSTskn_K, press_atm)

    # units in mol . L-1 . atm-1
    K0fnd = sol.solubility_woolf2016(SSSfnd, SSTfnd_K, SSTdelta, press_atm)
    K0skn = sol.solubility_woolf2016(SSSskn, SSTskn_K, SSTdelta, press_atm)

    # molar mass of carbon (gC/mol * kg/g)
    mC = 12.0108 * 1000  # kg . mol-1

    # CONC : UNIT ANALYSIS
    #         solubility         *  pCO2 *  molar mass
    # conc = (mol . L-1 . atm-1) * (atm) * (kg . mol-1)
    # conc = mol. mol-1 . L-1 . atm . atm-1 * kg
    # conc = kg . L-1    |||    gC . m-3
    # Bulk uses skin, equilibrium and rapid use foundation for concSEA
    concSEA = K0fnd * fCO2sea * mC
    concAIR = K0skn * fCO2air * mC

    # KW : UNIT ANALYSIS
    # kw = (cm / 100) / (hr / 24)
    # kw = m . day-1
    kw = kw_func(wind_ms, SSTskn_C) * (24 / 100)

    # FLUX : UNIT ANALYSIS
    # flux = (m . day-1) * (g . m-3)
    # flux = gC . m . m-3 . day-1
    # flux = gC . m-2 . day-1
    CO2flux_woolfe = kw * (concSEA - concAIR)

    if isinstance(var, DataArray):
        kw_name = kw_func.__name__[2:]
        attributes = dict(
            units="gC / m2 / day",
            description=f"sea-air CO2 fluxes calculated with {kw_name}",
            long_name="sea-air CO2 fluxes",
        )

        CO2flux_woolfe = DataArray(
            data=CO2flux_woolfe, coords=var.coords, attrs=attributes
        )

    return CO2flux_woolfe


@preserve_xda
def flux_bulk(
    temp_bulk_C,
    salt_bulk,
    pCO2_bulk_uatm,
    pCO2_air_uatm,
    press_hPa,
    wind_ms,
    kw_func=gas_transfer_velocity.k_Ni00,
):
    """
    Calculates bulk air-sea CO2 fluxes: FCO2 = kw * K0 * dfCO2, without
    defining skin and foundation conc. differences as in the RAPID model.

    Parameters
    ----------
    temp_bulk_C : np.array
        temperature from OISST in degCelcius with an allowable range of [-2:45]
    salt_bulk : np.array
        salinity from EN4 in PSU. Allowable range [5 : 50]
    pCO2_bulk_uatm : np.array
        partial pressure of CO2 in the sea in micro-atmospheres. Allowable
        range is [50 : 1000]
    pCO2_air_uatm : np.array
        partial pressure of CO2 in the air in micro-atmospheres. Allowable
        range is [50 : 1000].
    press_hPa : np.array
        atmospheric pressure in hecto-Pascals with allowable range [500:1500]
    wind_ms : np.array
        wind speed in metres per second with an allowable range of [0 : 40]
    kw_func : callable
        a function that returns the gas transfer velocity in cm/hr. The default
        is the gas transfer volicty as calculated by Ho et al. (2006). This
        is the prefered method of Goddijn-Murphy et al. (2016). Other functions
        are available in the `gas_transfer` class. If you'd like to use custom
        inputs must be wind speed (m/s) and temperature (degC) and output must
        be cm/hr.

    Reurns
    ------
    FCO2 : np.array
        Sea-air CO2 flux where positive is out of the ocean and negative is
        into the ocean. Units are gC.m-2.day-1 (grams Carbon per metre squared
        per day)
    """
    from numpy import array

    press_atm = array(press_hPa) / 1013.25

    SSTfnd_C = array(temp_bulk_C)
    SSTfnd_K = SSTfnd_C + 273.15

    SSSfnd = array(salt_bulk)

    pCO2sea = array(pCO2_bulk_uatm) * 1e-6  # to atm
    pCO2air = array(pCO2_air_uatm) * 1e-6

    # checking units
    SSTfnd_K = check.temp_K(SSTfnd_K)
    SSSfnd = check.salt(SSSfnd)
    press_atm = check.pres_atm(press_atm)
    pCO2sea = check.CO2_mol(pCO2sea)
    pCO2air = check.CO2_mol(pCO2air)
    wind_ms = check.wind_ms(wind_ms)

    fCO2sea = f2p.pCO2_to_fCO2(pCO2sea, SSTfnd_K, press_atm)
    fCO2air = f2p.pCO2_to_fCO2(pCO2air, SSTfnd_K, press_atm)

    K0blk = sol.solubility_weiss1974(SSSfnd, SSTfnd_K, press_atm)

    # molar mas of carbon in g . mmol-1
    mC = 12.0108 * 1000  # (g . mol-1) / (mmol . mol-1)

    # KW : UNIT ANALYSIS
    # kw = (cm . hr-1) * hr . day-1 . cm-1 . m
    # kw = m . day-1
    kw = kw_func(wind_ms, SSTfnd_C) * (24 / 100)

    # flux = (m . day-1) .  (mol . L-1 . atm-1) . atm . (gC . mmol-1)
    # flux = (m . day-1) . (mmol . m-3 . atm-1) . atm . (gC . mmol-1)
    # flux = gC . m-2 . day-1
    CO2flux_bulk = kw * K0blk * (fCO2sea - fCO2air) * mC

    kw_name = kw_func.__name__[2:]
    meta = dict(
        units="gC / m2 / day",
        description=f"sea-air CO2 fluxes calculated with {kw_name}",
        long_name="sea-air CO2 fluxes",
    )

    return CO2flux_bulk, meta
