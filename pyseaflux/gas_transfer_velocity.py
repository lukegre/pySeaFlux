"""
Gas Transfer Velocity
---------------------

Modulates the magnitude of the flux between the atmosphere and the ocean.
"""

from . import check_units as check


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
    from numpy import around, array, isnan, nan_to_num, nansum

    def calculate_scaled_alpha(U2, Sc, ice, scaling):
        """scales alpha"""
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

    Args:
        temp_C (array): temperature in degrees C

    Returns:
        array: Schmidt number (dimensionless)

    Examples:
        >>> schmidt_number(20)  # from Wanninkhof (2014)
        668.344

    References:
        Jähne, B., Heinz, G., & Dietrich, W. (1987). Measurement of the
        diffusion coefficients of sparingly soluble gases in water. Journal
        of Geophysical Research: Oceans, 92(C10), 10767–10776.
        https://doi.org/10.1029/JC092iC10p10767
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

    Note:
        This is an old parameterization and we recommend using updated
        parameterisations that are calculated based on the wind product you
        choose to use. We include this parameterisation based purely for
        legacy purposes.

    Args:
        wind_ms (array): wind speed in m/s
        temp_C (array): temperature in degrees C

    Returns:
        kw (array): gas transfer velocity (k600) in cm/hr

    References:
        Liss, P. S., & Merlivat, L. (1986). The Role of Air-Sea Exchange
        in Geochemical Cycling (Vol. 1983, Issue June 1983).
        D. Reidel Publishing Company.
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


def k_Wa92(wind_second_moment, temp_C):
    """
    Calculates the gas transfer coeffcient for CO2 using the formulation
    of Wanninkhof (1992)

    Note:
        This is an old parameterization and we recommend using updated
        parameterisations that are calculated based on the wind product you
        choose to use. We include this parameterisation based purely for
        legacy purposes.

    The gas transfer velocity is scaled from instantaneous wind speeds.
    The study applies a correction to the scaling (0.39) based on instantaneous
    wind speeds to lower it to 0.31. This correction is based on the variability
    of wind.

    .. math::
        k_{660} = 0.31 \\cdot U^2


    Args:
        wind_second_moment (array): wind speed squared in m2/s2. Note that the
        second moment should be calculated at the native resolution of the
        wind to avoid losses of variability when taking the square product.
        temp_C (array): temperature in degrees C

    Returns:
        kw (array): gas transfer velocity (k660) in cm/hr

    References:
        Wanninkhof, R. H. (1992). Relationship between wind speed and gas
        exchange over the ocean. Journal of Geophysical Research, 97(C5),
        7373. https://doi.org/10.1029/92JC00188
    """
    from numpy import array

    U2 = array(wind_second_moment)
    T = array(temp_C)

    check.temp_K(T + 273.15)

    Sc = schmidt_number(temp_C)
    k = (0.31 * U2) * (660 / Sc) ** 0.5

    return k


def k_Wa99(wind_ms, temp_C):
    """
    Calculates the gas transfer coeffcient for CO2 using the formulation
    of Wanninkhof and McGillis (1999)

    The gas transfer velocity has been scaled for in-situ short term wind
    products. Note that using this function for any other wind product is not
    correct.

    .. math::
        k_{600} = 0.0283 \\cdot U^3

    Args:
        wind_ms (array): wind speed in m/s
        temp_C (array): temperature in degrees C

    Returns:
        kw (array): gas transfer velocity (k600) in cm/hr

    References:
        Wanninkhof, R. H., & McGillis, W. R. (1999). A cubic relationship
        between air-sea CO2 exchange and wind speed. Geophysical Research
        Letters, 26(13), 1889–1892. https://doi.org/10.1029/1999GL900363
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

    .. math::
        k_{600} = 0.333 \\cdot U + 0.222 \\cdot U^2

    Args:
        wind_ms (array): wind speed in m/s
        temp_C (array): temperature in degrees C

    Returns:
        kw (array): gas transfer velocity (k600) in cm/hr

    References:
        Nightingale, P. D., Malin, G., Law, C. S., Watson, A. J., Liss, P. S.,
        Liddicoat, M. I., Boutin, J., & Upstill-Goddard, R. C. (2000). In
        situ evaluation of air-sea gas exchange parameterizations using
        novel conservative and volatile tracers. In Global Biogeochemical
        Cycles (Vol. 14, Issue 1, p. 373). https://doi.org/10.1029/1999GB900091
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

    The gas transfer velocity has been scaled for in-situ short term wind
    products. Note that using this function for any other wind product is not
    correct.

    .. math::
        k_{660} = 3.3 + 0.026 \\cdot U^3

    Args:
        wind_ms (array): wind speed in m/s
        temp_C (array): temperature in degrees C

    Returns:
        kw (array): gas transfer velocity (k660) in cm/hr

    References:
        McGillis, W. R., Edson, J. B., Ware, J. D., Dacey, J. W. H., Hare, J. E.,
        Fairall, C. W., & Wanninkhof, R. H. (2001). Carbon dioxide flux
        techniques performed during GasEx-98. Marine Chemistry, 75(4), 267–280.
        https://doi.org/10.1016/S0304-4203(01)00042-1
    """
    from numpy import array

    U = array(wind_ms)
    T = array(temp_C)

    check.temp_K(T + 273.15)
    check.wind_ms(U)

    Sc = schmidt_number(temp_C)
    k = 3.3 + (0.026 * U ** 3) * (660 / Sc) ** 0.5

    return k


def k_Ho06(wind_second_moment, temp_C):
    """
    Calculates the gas transfer coeffcient for CO2 using the formulation
    of Ho et al. (2006)

    The gas transfer velocity is for the QuickSCAT satellite wind product.
    Note that using this function for any other wind product is stricktly
    speaking not correct.

    .. math::
        k_{600} = 0.266 \\cdot U^2

    The parameterization is based on the SOLAS Air-Sea Gas Exchange (SAGE)
    experiment.

    Args:
        wind_ms (array): wind speed in m/s
        temp_C (array): temperature in degrees C

    Returns:
        kw (array): gas transfer velocity (k600) in cm/hr

    References:
        Ho, D. T., Law, C. S., Smith, M. J., Schlosser, P., Harvey, M., & Hill, P.
        (2006). Measurements of air-sea gas exchange at high wind speeds in the Southern
        Ocean: Implications for global parameterizations. Geophysical Research Letters,
        33(16), 1–6. https://doi.org/10.1029/2006GL026817
    """
    from numpy import array

    U2 = array(wind_second_moment)
    T = array(temp_C)

    check.temp_K(T + 273.15)

    Sc = schmidt_number(temp_C)
    k = (0.266 * U2) * (600 / Sc) ** 0.5

    return k


def k_Sw07(wind_second_moment, temp_C):
    """
    Calculates the gas transfer coeffcient for CO2 using the formulation
    Wanninkhof (1992) rescaled by Sweeny et al (2007)

    The gas transfer velocity has been scaled for the NCEP/NCAR reanalysis 1
    product. Note that using this function for any other wind product is not
    correct.

    .. math::
        k_{660} = 0.27 \\cdot U^2

    Args:
        wind_second_moment (array): wind speed squared in m2/s2. Note that the
            second moment should be calculated at the native resolution of the
            wind to avoid losses of variability when taking the square product.
        temp_C (array): temperature in degrees C

    Returns:
        kw (array): gas transfer velocity (k660) in cm/hr

    References:
        Sweeney, C., Gloor, E., Jacobson, A. R., Key, R. M., McKinley, G. A.,
        Sarmiento, J. L., & Wanninkhof, R. H. (2007). Constraining global
        air-sea gas exchange for CO2 with recent bomb 14C measurements.
        Global Biogeochemical Cycles, 21(2). https://doi.org/10.1029/2006GB002784
    """
    from numpy import array

    U2 = array(wind_second_moment)
    T = array(temp_C)

    check.temp_K(T + 273.15)

    Sc = schmidt_number(temp_C)
    k = (0.27 * U2) * (660 / Sc) ** 0.5

    return k


def k_Wa09(wind_ms, temp_C):
    """
    Calculates the gas transfer coeffcient for CO2 using the formulation
    of Wanninkhof et al. (2009)

    The gas transfer velocity has been scaled for the Cross-Calibrated Multi-
    Platform (CCMP) Winds product. Note that using this function for any other
    wind product is not correct.

    .. math::
        k_{660} = 3.0 + 0.1 \\cdot U + 0.064 \\cdot U^2 + 0.011 \\cdot U^3

    Args:
        wind_ms (array): wind speed in m/s
        temp_C (array): temperature in degrees C

    Returns:
        kw (array): gas transfer velocity (k660) in cm/hr

    References:
        Wanninkhof, R. H., Asher, W. E., Ho, D. T., Sweeney, C., & McGillis,
        W. R. (2009). Advances in Quantifying Air-Sea Gas Exchange and
        Environmental Forcing*. Annual Review of Marine Science, 1(1),
        213–244. https://doi.org/10.1146/annurev.marine.010908.163742
    """
    from numpy import array

    U = array(wind_ms)
    T = array(temp_C)

    check.temp_K(T + 273.15)
    check.wind_ms(U)

    Sc = schmidt_number(temp_C)
    k = (3.0 + 0.1 * U + 0.064 * U ** 2 + 0.011 * U ** 3) * (660 / Sc) ** 0.5

    return k


def k_Wa14(wind_second_moment, temp_C):
    """
    Calculates the gas transfer coeffcient for CO2 using the formulation
    of Wanninkhof et al. (2014)

    The gas transfer velocity has been scaled for the Cross-Calibrated Multi-
    Platform (CCMP) Winds product. Note that using this function for any other
    wind product is not correct.

    .. math::
        k_{660} = 0.251 \\cdot U^2

    Args:
        wind_second_moment (array): wind speed squared in m2/s2. Note that the
            second moment should be calculated at the native resolution of the
            wind to avoid losses of variability when taking the square product.
        temp_C (array): temperature in degrees C

    Returns:
        kw (array): gas transfer velocity (k660) in cm/hr

    References:
        Wanninkhof, R. H. (2014). Relationship between wind speed and gas
        exchange over the ocean revisited. Limnology and Oceanography:
        Methods, 12(JUN), 351–362. https://doi.org/10.4319/lom.2014.12.351
    """
    from numpy import array

    U2 = array(wind_second_moment)
    T = array(temp_C)

    check.temp_K(T + 273.15)

    Sc = schmidt_number(temp_C)
    k = 0.251 * U2 * (660 / Sc) ** 0.5

    return k
