from . import unit_checks as check


def solubility_weiss1974(salt, temp_K, press_atm=1):
    """
    Calculates the solubility of CO2 in sea water for the calculation of
    air-sea CO2 fluxes. We use the formulation by Weiss (1974) summarised in
    Wanninkhof (2014).

    Parameters
    ----------
    salt : np.array
        salinity in PSU
    temp_K : np.array
        temperature in deg Kelvin
    press_atm : np.array
        pressure in atmospheres. Used in the solubility correction for water
        vapour pressure. If not given, assumed that press_atm is 1atm

    Returns
    -------
    K0 : np.array
        solubility of CO2 in seawater in mol/L/atm

    Examples
    --------
    # from Weiss (1974) Table 2 but with pH2O correction
    >>> solubility_weiss1974(35, 299.15)
    0.029285284543519093
    """

    from numpy import array, exp, log

    T = array(temp_K)
    S = array(salt)
    P = array(press_atm)

    check.temp_K(T)
    check.salt(S)
    check.pres_atm(P)

    # from table in Wanninkhof 2014
    a1 = -58.0931
    a2 = +90.5069
    a3 = +22.2940
    b1 = +0.027766
    b2 = -0.025888
    b3 = +0.0050578

    T100 = T / 100
    K0 = exp(
        a1
        + a2 * (100 / T)
        + a3 * log(T100)
        + S * (b1 + b2 * T100 + b3 * T100 ** 2)
    )

    pH2O = vapress_weiss1980(S, T)
    K0 = K0 / (P - pH2O)

    return K0  # units mol/L/atm


def solubility_woolf2016(salt, temp_K, deltaT, press_atm=1):
    """
    A wrapper around solubility calculated using the Weiss (1974) approach.
    This is taken from the FluxEngine script.

    Parameters
    ----------
    temp_K : np.array
        temperature of sea water at the desired level (e.g. skin)
    salt : np.array
        salinity of seawater
    deltaT : np.array
        SST differences (foundation - skin)

    Returns
    -------
    K0 : np.array
        solubility of CO2 in seawater in mol/L/atm
    """
    K0 = solubility_weiss1974(salt, temp_K, press_atm)

    return K0 * (1 - 0.015 * deltaT)


def schmidt_number(temp_C):
    """
    Calculates the Schmidt number as defined by Jahne et al. (1987) and listed
    in Wanninkhof (2014) Table 1.

    Parameters
    ----------
    temp_C : np.array
        temperature in degrees C

    Returns
    -------
    Sc : np.array
        Schmidt number (dimensionless)

    Examples
    --------
    >>> schmidt_number(20)  # from Wanninkhof (2014)
    668.344

    """

    from numpy import array

    T = array(temp_C)
    check.temp_K(T + 273.15)

    a = +2116.8
    b = -136.25
    c = +4.7353
    d = -0.092307
    e = +0.0007555

    Sc = a + b * T + c * T ** 2 + d * T ** 3 + e * T ** 4

    return Sc


def pressure_height_correction(pres_hPa, tempSW_C, sensor_height=10.0):
    """
    Returns exact sea level pressure if the sensor is measuring at height

    Parameters
    ----------
    pres_hPa : np.array
        Pressure in kiloPascal measured at height
    temp_C : np.array
        Temperature of the seawater in deg C
    sensor_height : float
        the height of the sensor above sea level. Can be negative if you want
        to convert SLP to sensor height pressure

    Return
    ------
    presCor_kPa : np.array
        height corrected pressure
    """
    from numpy import array

    P = array(pres_hPa) * 100  # pressure in Pascal
    T = array(tempSW_C) + 273.15  # temperature in Kelvin

    check.temp_K(T)
    check.pres_atm(P / 101325)

    # Correction for pressure based on sensor height
    R = 8.314  # universal gas constant (J/mol/K)
    M = 0.02897  # molar mass of air in (kg/mol) - Wikipedia
    # Density of air at a given temperature. Here we assume
    # that the air temp is the same as the intake temperature
    d = P / (R / M * T)
    g = 9.8  # gravity in (m/s2)
    h = -sensor_height  # height in (m)
    # correction for atmospheric
    press_height_corr_hpa = (P - (d * g * h)) / 100.0

    return press_height_corr_hpa


def virial_coeff(temp_K, pres_atm, xCO2_mol=None):
    from numpy import exp

    """
    Calculate the ideal gas correction factor for converting pCO2 to fCO2.
    fCO2 = pCO2 * virial_expansion
    pCO2 = fCO2 / virial_expansion

    Based on the Lewis and Wallace 1998 Correction.

    Parameters
    ----------
    press_hPa : np.array
        uncorrected pressure in hPa
    temp_K : np.array
        temperature in degrees Kelvin
    xCO2_mol : np.array
        mole fraction of CO2. Can be pCO2/fCO2 if xCO2 is not defined or can
        leave this as undefined as makes only a small impact on output

    Return
    ------
    virial_expression : np.array
        the factor to multiply with pCO2. Unitless

    Examples
    --------
    The example below is from Dickson et al. (2007)
    >>> 350 * virial_coeff(298.15, 1)  # CO2 [uatm] * correction factor
    348.8836492182758

    References
    ----------
    Weiss, R. (1974). Carbon dioxide in water and seawater: the solubility of a
        non-ideal gas. Marine Chemistry, 2(3), 203–215.
        https://doi.org/10.1016/0304-4203(74)90015-2
    Compared with the Seacarb package in R
    """
    from numpy import array

    T = array(temp_K)
    P = array(pres_atm)
    C = array(xCO2_mol)
    R = 82.057  # gas constant for ATM

    check.temp_K(T)
    check.pres_atm(P)

    # B is the virial coefficient for pure CO2
    B = -1636.75 + 12.0408 * T - 0.0327957 * T ** 2 + 3.16528e-5 * T ** 3
    # d is the virial coefficient for CO2 in air
    d = 57.7 - 0.118 * T

    # "x2" term often neglected (assumed = 1) in applications of Weiss's
    # (1974) equation 9
    if xCO2_mol is not None:
        check.CO2_mol(C)
        x2 = (1 - C) ** 2
    else:
        x2 = 1

    ve = exp(P * (B + 2 * x2 * d) / (R * T))

    return ve


def vapress_weiss1980(salt, temp_K):
    """
    Calculates the water vapour pressure of seawater at a given salinity and
    temerature using the methods defined in Weiss (1974)

    Parameters
    ----------
    salt : np.array
        salinity
    temp_K : np.array
        temperature in deg Kelvin

    Returns
    -------
    sea_vapress : np.array
        sea water vapour pressure in atm

    Examples
    --------
    >>> vapress_weiss1980(35, 25+273.15)  # tempC + 273.15
    0.03065529996317971

    """
    from numpy import array, exp, log

    T = array(temp_K)
    S = array(salt)

    check.temp_K(T)
    check.salt(S)

    # Equation comes straight from Weiss and Price (1980)
    pH2O = exp(
        +24.4543 - 67.4509 * (100 / T) - 4.8489 * log(T / 100) - 0.000544 * S
    )

    return pH2O


def vapress_dickson2007(salt, temp_K):
    """
    Calculates the water vapour pressure of seawater at a given salinity and
    temerature using the methods defined in Dickson et al. (2007; CO2 manual)

    Parameters
    ----------
    salt : np.array
        salinity
    temp_K : np.array
        temperature in deg Kelvin

    Returns
    -------
    sea_vapress : np.array
        sea water vapour pressure in atm

    Examples
    --------
    >>> vapress_dickson2007(35, 298.15)  # from Dickson et al. (2007) Ch 5.3.2
    0.030698866245809465

    """
    from numpy import array, exp

    T = array(temp_K)
    S = array(salt)

    check.temp_K(T)

    ###################################################
    # WATER VAPOUR PRESSURE FOR PURE WATER
    ###################################################
    # alpha coefficients from Wafner and Pruss, (2002)
    a1 = -7.85951783
    a2 = +1.84408259
    a3 = -11.7866497
    a4 = +22.6807411
    a5 = -15.9618719
    a6 = +1.80122502
    # critical points for water
    Pc = 22.064 / 101325.0e-6  # convert to atmosphers
    Tc = 647.096
    # zeta numbers correspond with alpha numbers
    z = 1 - T / Tc
    z1 = z
    z2 = z ** 1.5
    z3 = z ** 3
    z4 = z ** 3.5
    z5 = z ** 4
    z6 = z ** 7.5
    # vapour pressure of pure water
    pure_water = Pc * exp(
        (Tc / T) * (a1 * z1 + a2 * z2 + a3 * z3 + a4 * z4 + a5 * z5 + a6 * z6)
    )

    ###################################################
    # WATER VAPOUR PRESSURE FOR SEA WATER
    ###################################################
    # osmotic coeffcients at 25C - Millero 1974
    c0 = +0.90799
    c1 = -0.08992
    c2 = +0.18458
    c3 = -0.07395
    c4 = -0.00221
    # total molality of dissolved species
    total_molality = 31.998 * S / (1e3 - 1.005 * S)
    B1 = total_molality * 0.5
    B2 = B1 ** 2
    B3 = B1 ** 3
    B4 = B1 ** 4
    osmotic_coeff = c0 + c1 * B1 + c2 * B2 + c3 * B3 + c4 * B4

    seawater = pure_water * exp(-0.018 * osmotic_coeff * total_molality)

    return seawater


def temperature_correction(temp_in, temp_out):
    """
    Calculate a correction factor for the temperature difference between the
    intake and equilibrator. This is based on the empirical relationship used
    in Takahashi et al. 1993.
    pCO2_Tout = pCO2_Tin * T_factor

    Parameters
    ----------
    temp_in : np.array
        temperature at which original pCO2 is measured
    temp_out : np.array
        temperature for which pCO2 should be represented

    Return
    ------
    factor : np.array
        a correction factor to be multiplied to pCO2 (unitless)

    References
    ----------
    Takahashi, Taro et al. (1993). Seasonal variation of CO2 and nutrients in
        the high-latitude surface oceans: A comparative study. Global
        Biogeochemical Cycles, 7(4), 843–878. https://doi.org/10.1029/93GB02263
    """

    from numpy import array, exp

    # see the Takahashi 1993 paper for full description

    Ti = array(temp_in)
    To = array(temp_out)

    factor = exp(0.0433 * (To - Ti) - 4.35e-05 * (To ** 2 - Ti ** 2))

    return factor
