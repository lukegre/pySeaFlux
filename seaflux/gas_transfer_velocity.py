from . import check_units as check
from .utils import preserve_xda


@preserve_xda
def kw_scaled(wind_speed, wind_grid_stdev, temp_C, ice_frac, scaling=16):
    """
    Calculates kw using the wind^2 formulation for different wind products and
    scales them according to a given global 14C bomb value.
    kw = a * <U^2> * (Sc / 660)^-0.5     (1)

    If kw for the wind product and resolution has already been calculated,
    then this step does not have to be performed. Else, the downscaling and
    wind product used will have an influence (Sweeney et al. 2007). Unlike in
    Sweeney et al. (2007), we assume that downscaled temperature will not
    have a large impact on the Schmidt number - this assumption breaks down
    when the time period is too long.

    Using formulation 1 we can calculate a:
    a = scaled / (<U^2> * (Sc / 660)^-0.5 * (1 - ice_frac))
    kw = a * <U^2> * (Sc / 660)^-0.5

    Parameters
    ----------
    wind_speed : np.array
        wind speed in m/s
    wind_grid_stdev : np.array
        variance of the wind speed in each grid cell for resampled data
    temp_C : np.array
        temperature in degrees C
    ice_frac : np.array
        the sea ice cover as a fraction
    scaling : int [16 cm/hr]
        the global estimate of 14C bomb flux. The default (16 cm/hr) is from
        Wanninkhof (2013). But Sweeney et al. (2007) suggest 14.6 cm/hr

    Returns
    -------
    kw_scaled : np.array
        gas transfer velocity in cm/hr
    attributes : dict
        a dictionary containing information about the calculation to be used
        for creating netCDF files

    """
    from numpy import array, nansum, nan_to_num, isnan, around

    def calculate_scaled_alpha(U2, Sc, ice, scaling):

        mask = isnan(U2) | isnan(Sc)

        weight = nan_to_num(1 - ice)
        weight[mask] = 0

        a = around(scaling / (nansum(U2 * Sc * weight) / nansum(weight)), 4)

        return a

    Uavg = array(wind_speed)
    Ustd = array(wind_grid_stdev)
    T = array(temp_C)
    ice = array(ice_frac)

    Uavg2 = Uavg ** 2
    U2 = Uavg2 + Ustd ** 2

    Sc = schmidt_number(T)
    Sc660 = (Sc / 660) ** -0.5

    a = calculate_scaled_alpha(U2, Sc660, ice, scaling)

    k = a * Uavg2 * Sc660

    attrs = {
        "long_name": "Gas transfer velocity for CO2",
        "units": "cm/hr",
        "alpha": a,
        "description": (
            f"k = a * <U2> * Sc660 where a is scaled to a "
            f"global 14C bomb estimate of {scaling} cm/hr. "
            f"Note that the global kw is weighted by sea ice "
            f"cover, thus the mean of kw will not be "
            f"{scaling} cm/hr. "
            f"<U2> is the second moment of the wind speed "
            f"calculated as <Uavg^2 + Ustd^2>. "
            f"Sc660 is the Schmidt number for CO2 normalised "
            f"to 660 as in Wanninkhof (1992): (Sc / 660)^-0.5. "
            f"For more details on this calculation see "
            f"Sarmiento and Gruber (2006), Wanninkhof (2014) "
            f"and Sweeney et al. (2007). "
        ),
    }
    return k, attrs


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

    check.temp_K(temp_C + 273.15)
    T = temp_C

    a = +2116.8
    b = -136.25
    c = +4.7353
    d = -0.092307
    e = +0.0007555

    Sc = a + b * T + c * T ** 2 + d * T ** 3 + e * T ** 4

    return Sc


def k_Li86(wind_ms, temp_C):
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
    from numpy import array, zeros_like

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

    return k


def k_Wa92(wind_ms, temp_C):
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
    from numpy import array

    U = array(wind_ms)
    T = array(temp_C)

    check.temp_K(T + 273.15)
    check.wind_ms(U)

    Sc = schmidt_number(temp_C)
    k = (0.39 * U ** 2) * (660 / Sc) ** 0.5

    return k


def k_Wa99(wind_ms, temp_C):
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
    from numpy import array

    U = array(wind_ms)
    T = array(temp_C)

    check.temp_K(T + 273.15)
    check.wind_ms(U)

    Sc = schmidt_number(temp_C)
    k = (0.0283 * U ** 3) * (600 / Sc) ** 0.5

    return k


def k_Ni00(wind_ms, temp_C):
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
    from numpy import array

    U = array(wind_ms)
    T = array(temp_C)

    check.temp_K(T + 273.15)
    check.wind_ms(U)

    Sc = schmidt_number(temp_C)
    k = (0.333 * U + 0.222 * U ** 2) * (600 / Sc) ** 0.5

    return k


def k_Mc01(wind_ms, temp_C):
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
    from numpy import array

    U = array(wind_ms)
    T = array(temp_C)

    check.temp_K(T + 273.15)
    check.wind_ms(U)

    Sc = schmidt_number(temp_C)
    k = 3.3 + (0.026 * U ** 3) * (660 / Sc) ** 0.5

    return k


def k_Ho06(wind_ms, temp_C):
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
    from numpy import array

    U = array(wind_ms)
    T = array(temp_C)

    check.temp_K(T + 273.15)
    check.wind_ms(U)

    Sc = schmidt_number(temp_C)
    k = (0.266 * U ** 2) * (600 / Sc) ** 0.5

    return k


def k_Sw07(wind_ms, temp_C):
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
    from numpy import array

    U = array(wind_ms)
    T = array(temp_C)

    check.temp_K(T + 273.15)
    check.wind_ms(U)

    Sc = schmidt_number(temp_C)
    k = (0.27 * U ** 2) * (660 / Sc) ** 0.5

    return k


def k_Wa09(wind_ms, temp_C):
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
    from numpy import array

    U = array(wind_ms)
    T = array(temp_C)

    check.temp_K(T + 273.15)
    check.wind_ms(U)

    Sc = schmidt_number(temp_C)
    k = (3.0 + 0.1 * U + 0.064 * U ** 2 + 0.011 * U ** 3) * (660 / Sc) ** 0.5

    return k


def k_Wa14(wind_ms, temp_C):
    """
    Calculates the gas transfer coeffcient for CO2 using the formulation
    of Wanninkhof et al. (2014)
        k660 = 0.251 * U^2

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
    from numpy import array

    U = array(wind_ms)
    T = array(temp_C)

    check.temp_K(T + 273.15)
    check.wind_ms(U)

    Sc = schmidt_number(temp_C)
    k = 0.251 * U ** 2 * (660 / Sc) ** 0.5

    return k
