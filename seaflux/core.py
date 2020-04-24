"""
Common tools used to correct and calculate pCO2, fCO2 and fluxes

Please cite the use of this code as Gregor et al. (2019) but make sure to
reference Woolf et al. (2016) if the woolf flux method is used. Same applies
for all other functions in this script (Weiss 1974, Dickson et al. 2007 ...)

https://doi.org/10.5194/gmd-2019-46
"""
from . import aux_eqs as eqs
from . import gas_transfer_CO2
from . import unit_checks as check
import warnings


def fCO2_to_pCO2(fCO2SW_uatm, tempSW_C, pres_hPa=1013.25, tempEQ_C=None):
    """
    Convert fCO2 to pCO2 for SOCAT in sea water. A simple version of the
    equation would simply be:
        pCO2sw = fCO2sw / virial_exp
    where the virial expansion is calculated without xCO2

    We get a simple approximate for equilibrator xCO2 with:
        xCO2eq = fCO2sw * deltaTemp(sw - eq) / press_eq

    pCO2sw is then calculated with:
        pCO2sw = fCO2sw / virial_exp(xCO2eq)

    Parameters
    ----------
    fCO2SW_uatm : array
        seawater fugacity of CO2 in micro atmospheres
    tempSW_C : array
        sea water temperature in degrees C
    pres_hPa : array
        equilibrator pressure in kilo Pascals
    tempEQ_C : array
        equilibrator temperature in degrees C

    Returns
    -------
    pCO2SW_uatm : array
        partial pressure of CO2 in seawater

    Note
    ----
    In FluxEngine, they account fully solve for the original xCO2 that is used
    in the calculation of the virial exponent. I use the first estimate of
    xCO2 (based on fCO2 rather than pCO2). The difference between the two
    approaches is so small that it is not significant to be concerned. Their
    correction is more precise, but the difference between their iterative
    correction and our approximation is on the order of 1e-14 atm (1e-8 uatm).

    Examples
    --------
    >>> fCO2_to_pCO2(380, 8)
    381.50806485658234
    >>> fCO2_to_pCO2(380, 8, pres_hPa=985)
    381.4659553134281
    >>> fCO2_to_pCO2(380, 8, pres_hPa=985, tempEQ_C=14)
    381.466027968504
    """

    from numpy import array

    # if equilibrator inputs are None, tempEQ=tempSW
    if tempEQ_C is None:
        tempEQ_was_None = True
        tempEQ_C = tempSW_C

    # standardise the inputs and convert units
    fCO2sw = array(fCO2SW_uatm) * 1e-6
    Tsw = array(tempSW_C) + 273.15
    Teq = array(tempEQ_C) + 273.15
    Peq = array(pres_hPa) / 1013.25

    # check if units make sense
    check.pres_atm(Peq)
    check.CO2_mol(fCO2sw)
    check.temp_K(Tsw)
    check.temp_K(Teq)

    # calculate the CO2 diff due to equilibrator and seawater temperatures
    # if statement is there to save a bit of time
    if tempEQ_was_None:
        dT = 1.
    else:
        dT = eqs.temperature_correction(Tsw, Teq)

    # a best estimate of xCO2 - this is an aproximation
    # one would have to use pCO2 / Peq to get real xCO2
    # Not getting the exact equilibrator xCO2
    xCO2eq = fCO2sw * dT / Peq

    pCO2SW = fCO2sw / eqs.virial_coeff(Tsw, Peq, xCO2eq)
    pCO2SW_uatm = pCO2SW * 1e6

    return pCO2SW_uatm


def pCO2_to_fCO2(pCO2SW_uatm, tempSW_C, pres_hPa=None, tempEQ_C=None):
    """
    Convert fCO2 to pCO2 for SOCAT in sea water. A simple version of the
    equation would simply be:
        fCO2sw = pCO2sw / virial_exp
    where the virial expansion is calculated without xCO2

    We get a simple approximate for equilibrator xCO2 with:
        xCO2eq = pCO2sw * deltaTemp(sw - eq) / press_eq

    fCO2sw is then calculated with:
        fCO2sw = pCO2sw * virial_exp(xCO2eq)

    Parameters
    ----------
    pCO2SW_uatm : array
        seawater fugacity of CO2 in micro atmospheres
    tempSW_C : array
        sea water temperature in degrees C/K
    tempEQ_C : array
        equilibrator temperature in degrees C/K
    pres_hPa : array
        pressure in kilo Pascals

    Returns
    -------
    fCO2SW_uatm : array
        partial pressure of CO2 in seawater

    Note
    ----
    In FluxEngine, they account for the change in xCO2. This error is so small
    that it is not significant to be concerned about it. Their correction is
    more precise, but the difference between their iterative correction and our
    approximation is less than 1e-14 atm (or 1e-8 uatm).

    Examples
    --------
    >>> pCO2_to_fCO2(380, 8)
    378.49789637942064
    >>> pCO2_to_fCO2(380, 8, pres_hPa=985)
    378.53967828231225
    >>> pCO2_to_fCO2(380, 8, pres_hPa=985, tempEQ_C=14)
    378.53960618459695
    """
    from numpy import array

    # if equilibrator inputs are None then make defaults Patm=1, tempEQ=tempSW
    if tempEQ_C is None:
        tempEQ_C = tempSW_C
    if pres_hPa is None:
        pres_hPa = 1013.25

    # standardise the inputs and convert units
    pCO2sw = array(pCO2SW_uatm) * 1e-6
    Tsw = array(tempSW_C) + 273.15
    Teq = array(tempEQ_C) + 273.15
    Peq = array(pres_hPa) / 1013.25

    # check if units make sense
    check.pres_atm(Peq)
    check.CO2_mol(pCO2sw)
    check.temp_K(Tsw)
    check.temp_K(Teq)

    # calculate the CO2 diff due to equilibrator and seawater temperatures
    dT = eqs.temperature_correction(Tsw, Teq)
    # a best estimate of xCO2 - this is an aproximation
    # one would have to use pCO2 / Peq to get real xCO2
    xCO2eq = pCO2sw * dT / Peq

    fCO2sw = pCO2sw * eqs.virial_coeff(Tsw, Peq, xCO2eq)
    fCO2sw_uatm = fCO2sw * 1e6

    return fCO2sw_uatm


def flux_woolf2016_rapid(
    temp_bulk_C,
    salt_bulk,
    pCO2_bulk_uatm,
    pCO2_air_uatm,
    press_hPa,
    wind_ms,
    kw_func=gas_transfer_CO2.k_Ni00,
    kw_scaling=None,
    cool_skin_bias=0.14,
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
    check.temp_K(SSTfnd_K)
    check.salt(SSSfnd)
    check.pres_atm(press_atm)
    check.CO2_mol(pCO2sea)
    check.CO2_mol(pCO2air)
    check.wind_ms(wind_ms)

    fCO2sea = pCO2sea * eqs.virial_coeff(SSTfnd_K, press_atm)
    fCO2air = pCO2air * eqs.virial_coeff(SSTskn_K, press_atm)

    # units in mol . L-1 . atm-1
    K0fnd = eqs.solubility_woolf2016(SSSfnd, SSTfnd_K, SSTdelta, press_atm)
    K0skn = eqs.solubility_woolf2016(SSSskn, SSTskn_K, SSTdelta, press_atm)

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
    kw = kw_func(wind_ms, SSTskn_C, kw_scaling) * (24 / 100)

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


def flux_bulk(
    temp_bulk_C,
    salt_bulk,
    pCO2_bulk_uatm,
    pCO2_air_uatm,
    press_hPa,
    wind_ms,
    kw_func=gas_transfer_CO2.k_Ni00,
    kw_scaling=None,
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
    from xarray import DataArray

    if isinstance(pCO2_bulk_uatm, DataArray):
        var = pCO2_bulk_uatm.copy()  # attribute preservation
    else:
        var = None

    press_atm = array(press_hPa) / 1013.25

    SSTfnd_C = array(temp_bulk_C)
    SSTfnd_K = SSTfnd_C + 273.15

    SSSfnd = array(salt_bulk)

    pCO2sea = array(pCO2_bulk_uatm) * 1e-6  # to atm
    pCO2air = array(pCO2_air_uatm) * 1e-6

    # checking units
    check.temp_K(SSTfnd_K)
    check.salt(SSSfnd)
    check.pres_atm(press_atm)
    check.CO2_mol(pCO2sea)
    check.CO2_mol(pCO2air)
    check.wind_ms(wind_ms)

    fCO2sea = pCO2sea * eqs.virial_coeff(SSTfnd_K, press_atm)
    fCO2air = pCO2air * eqs.virial_coeff(SSTfnd_K, press_atm)

    K0blk = eqs.solubility_weiss1974(SSSfnd, SSTfnd_K, press_atm)

    # molar mas of carbon in g . mmol-1
    mC = 12.0108 * 1000  # (g . mol-1) / (mmol . mol-1)

    # KW : UNIT ANALYSIS
    # kw = (cm . hr-1) * hr . day-1 . cm-1 . m
    # kw = m . day-1
    kw = kw_func(wind_ms, SSTfnd_C, kw_scaling) * (24 / 100)

    # flux = (m . day-1) .  (mol . L-1 . atm-1) . atm . (gC . mmol-1)
    # flux = (m . day-1) . (mmol . m-3 . atm-1) . atm . (gC . mmol-1)
    # flux = gC . m-2 . day-1
    CO2flux_bulk = kw * K0blk * (fCO2sea - fCO2air) * mC

    if isinstance(var, DataArray):
        kw_name = kw_func.__name__[2:]
        attributes = dict(
            units="gC / m2 / day",
            description=f"sea-air CO2 fluxes calculated with {kw_name}",
            long_name="sea-air CO2 fluxes",
        )

        CO2flux_bulk = DataArray(
            data=CO2flux_bulk, coords=var.coords, attrs=attributes
        )

    return CO2flux_bulk
