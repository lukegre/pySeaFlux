"""
High level functions
--------------------


"""
from . import check_units as check
from . import gas_transfer_velocity
from . import solubility as sol


def flux_bulk(
    temp_bulk_C,
    salt_bulk,
    pCO2_bulk_uatm,
    pCO2_air_uatm,
    press_hPa,
    wind_ms,
    kw_func=gas_transfer_velocity.k_Ni00,
    as_dict=False,
):
    """
    Calculates bulk air-sea CO2 fluxes

    .. math::
        FCO_2 = k_w \\cdot K_0 \\cdot \\Delta pCO_2

    .. warning::
        Note that you have to be very aware of the impact on resampled winds
        where there is a loss of variability when averaging. For now, I recommend
        to use the lower level functions

    Args:
        temp_bulk_C (array): temperature from OISST in degCelcius with an allowable
            range of [-2:45]
        salt_bulk (array): salinity from EN4 in PSU. Allowable range [5 : 50]
        pCO2_bulk_uatm (array): partial pressure of CO2 in the sea in micro-atmospheres.
            Allowable range is [50 : 1000]
        pCO2_air_uatm (array): partial pressure of CO2 in the air in micro-atmospheres.
            Allowable range is [50 : 1000].
        press_hPa (array): atmospheric pressure in hecto-Pascals with allowable range [500:1500]
        wind_ms (array): wind speed in metres per second with an allowable range of [0 : 40]
        kw_func (callable): a function that returns the gas transfer velocity in cm/hr.
            The default is the gas transfer volicty as calculated by Ho et al. (2006).
            This is the preferred method of Goddijn-Murphy et al. (2016). Other functions
            are available in the `gas_transfer` class. If you'd like to use custom
            inputs must be wind speed (m/s) and temperature (degC) and output must
            be cm/hr.
        as_dict (bool, optional): If True, will return a dictionary where output is stored
            under ``data`` and meta is stored under ``attrs``. This dictionary can be passed
            to ``xr.DataArray`` for convenience.

    Returns:
        array:
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

    K0blk = sol.solubility_weiss1974(SSSfnd, SSTfnd_K, press_atm)

    # molar mas of carbon in g . mmol-1
    mC = 12.0108 * 1000  # (g . mol-1) / (mmol . mol-1)

    """unit analysis
    kw = (cm . hr-1) * hr . day-1 . cm-1 . m
    kw = m . day-1
    """
    kw = kw_func(wind_ms, SSTfnd_C) * (24 / 100)

    """ unit analysis
    flux = (m . day-1) .  (mol . L-1 . atm-1) . atm . (gC . mmol-1)
    flux = (m . day-1) . (mmol . m-3 . atm-1) . atm . (gC . mmol-1)
    flux = gC . m-2 . day-1
    """
    CO2flux_bulk = kw * K0blk * (pCO2sea - pCO2air) * mC

    if as_dict:
        kw_name = kw_func.__name__[2:]
        return dict(
            data=CO2flux_bulk,
            attrs=dict(
                units="gC . m^-2 . day^-1",
                description=f"sea-air CO2 fluxes calculated with {kw_name}",
                long_name="sea-air CO2 fluxes",
            ),
        )

    return CO2flux_bulk
