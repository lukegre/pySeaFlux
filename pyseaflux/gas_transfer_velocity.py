"""
Gas Transfer Velocity
---------------------

Modulates the magnitude of the flux between the atmosphere and the ocean.
"""


def _add_xarray_attrs(func):
    """A helper function to add attributes to xarray."""
    import re

    from functools import wraps

    from xarray import DataArray

    def get_refs(func):
        """gets the reference from the docs where the reference is formatted
        as \nReferences: ..."""
        found = re.findall("References:(.*)", func.__doc__, flags=re.DOTALL)
        if any(found):
            ref = " ".join([s.strip() for s in found[0].split("\n")]).strip()
            return ref
        else:
            return ""

    def get_code(func):
        """get the formulation of kw from the function code. Requires the line
        to start with k = ..."""
        import inspect

        raw = "".join(inspect.getsource(func))
        found = re.findall("(k = .*)", raw)

        if any(found):
            code = found[0]
            return code
        else:
            return ""

    @wraps(func)
    def wrapper(*args, **kwargs):
        """wrapper that adds the xarray metadata if input is xarray"""
        out = func(*args, **kwargs)
        if isinstance(out, DataArray):
            # full reference based on function name. This is manually added
            names = {
                "k_Li86": "Liss and Merlivat (1986)",
                "k_Wa92": "Wanninkhof (1992)",
                "k_Wa99": "Wanninkhof and McGillis(1999)",
                "k_Ni00": "Nightingale et al. (2000)",
                "k_Mc01": "McGillis et al (2001)",
                "k_Ho06": "Ho et al. (2006)",
                "k_Sw07": "Sweeney et al. (2007)",
                "k_Wa09": "Wanninkhof et al. (2009)",
                "k_Wa14": "Wanninkhof et al. (2014)",
            }
            name = names[func.__name__]
            out = out.assign_attrs(
                units="cm/hr",
                description=f"gas transfer velocity of CO2 in seawater using {name}",
                reference=get_refs(func),
                formulation=get_code(func),
            )
        return out

    return wrapper


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
    from numpy import nanmedian

    if nanmedian(temp_C) > 270:
        raise ValueError("temperature is not in degC")

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
    from numpy import zeros_like

    U = wind_ms
    T = temp_C

    Sc = schmidt_number(T)
    k = zeros_like(temp_C)

    i1 = U <= 3.6
    i2 = (U > 3.6) & (U < 13.0)
    i3 = U >= 13.0

    k[i1] = (0.17 * U[i1]) * (Sc[i1] / 600) ** (-2.0 / 3.0)
    k[i2] = ((U[i2] - 3.4) * 2.8) * (600 / Sc[i2]) ** 0.5
    k[i3] = ((U[i3] - 8.4) * 5.9) * (600 / Sc[i3]) ** 0.5

    return k


@_add_xarray_attrs
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

    U2 = wind_second_moment

    Sc = schmidt_number(temp_C)
    k = (0.31 * U2) * (660 / Sc) ** 0.5

    return k


@_add_xarray_attrs
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

    U = wind_ms

    Sc = schmidt_number(temp_C)
    k = (0.0283 * U ** 3) * (600 / Sc) ** 0.5

    return k


@_add_xarray_attrs
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

    U = wind_ms

    Sc = schmidt_number(temp_C)
    k = (0.333 * U + 0.222 * U ** 2) * (600 / Sc) ** 0.5

    return k


@_add_xarray_attrs
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

    U = wind_ms

    Sc = schmidt_number(temp_C)
    k = 3.3 + (0.026 * U ** 3) * (660 / Sc) ** 0.5

    return k


@_add_xarray_attrs
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

    U2 = wind_second_moment

    Sc = schmidt_number(temp_C)
    k = (0.266 * U2) * (600 / Sc) ** 0.5

    return k


@_add_xarray_attrs
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

    U2 = wind_second_moment

    Sc = schmidt_number(temp_C)
    k = (0.27 * U2) * (660 / Sc) ** 0.5

    return k


@_add_xarray_attrs
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

    U = wind_ms

    Sc = schmidt_number(temp_C)
    k = (3.0 + 0.1 * U + 0.064 * U ** 2 + 0.011 * U ** 3) * (660 / Sc) ** 0.5

    return k


@_add_xarray_attrs
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

    U2 = wind_second_moment

    Sc = schmidt_number(temp_C)
    k = 0.251 * U2 * (660 / Sc) ** 0.5

    return k
