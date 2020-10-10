from . import check_units as check


def weiss1980(salt, temp_K):
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
    from numpy import exp, log

    T = check.temp_K(temp_K)
    S = check.salt(salt)

    # Equation comes straight from Weiss and Price (1980)
    pH2O = exp(+24.4543 - 67.4509 * (100 / T) - 4.8489 * log(T / 100) - 0.000544 * S)

    return pH2O


def dickson2007(salt, temp_K):
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
    from numpy import exp

    T = check.temp_K(temp_K)
    S = check.salt(salt)

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
