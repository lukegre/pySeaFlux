from .aux_eqs import schmidt_number
from . import unit_checks as check


def k_Li86(wind_ms, temp_C, scaled=None):
    """
    Calculates the gas transfer coeffcient for CO2 using the formulation
    of Liss and Merlivat (1986)

    Parameters
    ----------
    wind_ms : np.array
        wind speed in m/s
    temp_C : np.array
        temperature in degrees C

    Returns
    -------
    kw : array
        gas transfer velocity (k600) in cm/hr
    """
    from numpy import array, nanmean, zeros_like

    U = array(wind_ms)
    T = array(temp_C)

    check.temp_K(T + 273.15)
    check.wind_ms(U)

    Sc = schmidt_number(T)
    k = zeros_like(temp_C)

    i1 = U <= 3.6
    i2 = (U > 3.6) & (U < 13.0)
    i3 = U >= 13.0

    k[i1] = (0.17 * U[i1]) * (Sc[i1] / 600) ** (-2.0 / 3.0)
    k[i2] = ((U[i2] - 3.4) * 2.8) * (600 / Sc[i2]) ** 0.5
    k[i3] = ((U[i3] - 8.4) * 5.9) * (600 / Sc[i3]) ** 0.5
    if scaled is not None:
        k *= scaled / nanmean(k)
    return k


def k_Wa92(wind_ms, temp_C, scaled=None):
    """
    Calculates the gas transfer coeffcient for CO2 using the formulation
    of Wanninkhof (1992)
        k660 = 0.39 * u^2

    Parameters
    ----------
    wind_ms : np.array
        wind speed in m/s
    temp_C : np.array
        temperature in degrees C

    Returns
    -------
    kw : array
        gas transfer velocity (k660) in cm/hr
    """
    from numpy import array, nanmean

    U = array(wind_ms)
    T = array(temp_C)

    check.temp_K(T + 273.15)
    check.wind_ms(U)

    Sc = schmidt_number(temp_C)
    k = (0.39 * U ** 2) * (660 / Sc) ** 0.5

    if scaled is not None:
        k *= scaled / nanmean(k)
    return k


def k_Sw07(wind_ms, temp_C, scaled=None):
    """
    Calculates the gas transfer coeffcient for CO2 using the formulation
    of Sweeny et al (2007) who scaled Wanninkhof (1992)
        k660 = 0.27 * u^2

    Parameters
    ----------
    wind_ms : np.array
        wind speed in m/s
    temp_C : np.array
        temperature in degrees C

    Returns
    -------
    kw : array
        gas transfer velocity (k660) in cm/hr
    """
    from numpy import array, nanmean

    U = array(wind_ms)
    T = array(temp_C)

    check.temp_K(T + 273.15)
    check.wind_ms(U)

    Sc = schmidt_number(temp_C)
    k = (0.27 * U ** 2) * (660 / Sc) ** 0.5

    if scaled is not None:
        k *= scaled / nanmean(k)
    return k


def k_Wa99(wind_ms, temp_C, scaled=None):
    """
    Calculates the gas transfer coeffcient for CO2 using the formulation
    of Wanninkhof (1999)
        k600 = 0.0283 * U^3

    Parameters
    ----------
    wind_ms : np.array
        wind speed in m/s
    temp_C : np.array
        temperature in degrees C

    Returns
    -------
    kw : array
        gas transfer velocity (k600) in cm/hr
    """
    from numpy import array, nanmean

    U = array(wind_ms)
    T = array(temp_C)

    check.temp_K(T + 273.15)
    check.wind_ms(U)

    Sc = schmidt_number(temp_C)
    k = (0.0283 * U ** 3) * (600 / Sc) ** 0.5
    if scaled is not None:
        k *= scaled / nanmean(k)
    return k


def k_Ni00(wind_ms, temp_C, scaled=None):
    """
    Calculates the gas transfer coeffcient for CO2 using the formulation
    of Nightingale et al (2000)
        k600 = (0.333 * U) + (0.222 * U^2)

    Parameters
    ----------
    wind_ms : np.array
        wind speed in m/s
    temp_C : np.array
        temperature in degrees C

    Returns
    -------
    kw : array
        gas transfer velocity (k600) in cm/hr
    """
    from numpy import array, nanmean

    U = array(wind_ms)
    T = array(temp_C)

    check.temp_K(T + 273.15)
    check.wind_ms(U)

    Sc = schmidt_number(temp_C)
    k = (0.333 * U + 0.222 * U ** 2) * (600 / Sc) ** 0.5

    if scaled is not None:
        k *= scaled / nanmean(k)
    return k


def k_Ho06(wind_ms, temp_C, scaled=None):
    """
    Calculates the gas transfer coeffcient for CO2 using the formulation
    of Ho et al (2006)
        k600 = 0.266 * U^2

    Parameters
    ----------
    wind_ms : np.array
        wind speed in m/s
    temp_C : np.array
        temperature in degrees C

    Returns
    -------
    kw : array
        gas transfer velocity (k600) in cm/hr
    """
    from numpy import array, nanmean

    U = array(wind_ms)
    T = array(temp_C)

    check.temp_K(T + 273.15)
    check.wind_ms(U)

    Sc = schmidt_number(temp_C)
    k = (0.266 * U ** 2) * (600 / Sc) ** 0.5

    if scaled is not None:
        k *= scaled / nanmean(k)
    return k


def k_Wa09(wind_ms, temp_C, scaled=None):
    """
    Calculates the gas transfer coeffcient for CO2 using the formulation
    of Wanninkhof et al. (2009)
        k660 = 3. + (0.1 * U) + (0.064 * U^2) + (0.011 * U^3)

    Parameters
    ----------
    wind_ms : np.array
        wind speed in m/s
    temp_C : np.array
        temperature in degrees C

    Returns
    -------
    kw : array
        gas transfer velocity (k660) in cm/hr
    """
    from numpy import array, nanmean

    U = array(wind_ms)
    T = array(temp_C)

    check.temp_K(T + 273.15)
    check.wind_ms(U)

    Sc = schmidt_number(temp_C)
    k = (3.0 + 0.1 * U + 0.064 * U ** 2 + 0.011 * U ** 3) * (660 / Sc) ** 0.5

    if scaled is not None:
        k *= scaled / nanmean(k)
    return k


def k_Mc01(wind_ms, temp_C, scaled=None):
    """
    Calculates the gas transfer coeffcient for CO2 using the formulation
    of McGillis et al. (2001)
        k660 = 3.3 + (0.026 * U^3)

    Parameters
    ----------
    wind_ms : np.array
        wind speed in m/s
    temp_C : np.array
        temperature in degrees C

    Returns
    -------
    kw : array
        gas transfer velocity (k660) in cm/hr
    """
    from numpy import array, nanmean

    U = array(wind_ms)
    T = array(temp_C)

    check.temp_K(T + 273.15)
    check.wind_ms(U)

    Sc = schmidt_number(temp_C)
    k = 3.3 + (0.026 * U ** 3) * (660 / Sc) ** 0.5

    if scaled is not None:
        k *= scaled / nanmean(k)
    return k
