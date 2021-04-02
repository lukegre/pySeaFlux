"""
Conversions of fCO2 - pCO2
--------------------------
"""


def fCO2_to_pCO2(fCO2SW_uatm, tempSW_C, pres_hPa=1013.25, tempEQ_C=None):
    """Convert fCO2 to pCO2 in sea water.

    If equilibrator temperature is provided, we get a simple approximate for
    equilibrator  :math:`xCO_2` that allows for the virial expansion to be
    calculated more accurately. If not, then a simple approximation is good
    enough. See the examples for the differences.

    .. math::
        pCO_2^{sw} = fCO_2^{sw} \\div virial(xCO_2^{eq})

    where :math:`xCO_2^{eq} = fCO_2^{sw} \\times \\Delta T^{(sw - eq)} \\div P^{eq}`


    Args:
        fCO2SW_uatm (array): seawater fugacity of CO2 in micro atmospheres
        tempSW_C (array): sea water temperature in degrees C
        pres_hPa (array, optional): equilibrator pressure in hecto Pascals. Defaults to 1013.25.
        tempEQ_C (array, optional): equilibrator temperature in degrees C. Defaults to None.

    Returns:
        array:
            partial pressure of CO2 in seawater

    Note:
        In FluxEngine, they account fully solve for the original xCO2 that is used
        in the calculation of the virial expansion. I use the first estimate of
        xCO2 (based on fCO2 rather than pCO2). The difference between the two
        approaches is so small that it is not significant to be concerned. Their
        correction is more precise, but the difference between their iterative
        correction and our approximation is on the order of 1e-14 atm (1e-8 uatm).

    Examples:
        >>> fCO2_to_pCO2(380, 8)
        381.50806485658234
        >>> fCO2_to_pCO2(380, 8, pres_hPa=985)
        381.4659553134281
        >>> fCO2_to_pCO2(380, 8, pres_hPa=985, tempEQ_C=14)
        381.466027968504
    """
    from . import auxiliary_equations as eqs
    from . import check_units as check

    # if equilibrator inputs are None, tempEQ=tempSW
    if tempEQ_C is None:
        tempEQ_was_None = True
        tempEQ_C = tempSW_C
    else:
        tempEQ_was_None = False

    # standardise the inputs and convert units
    fCO2sw = check.CO2_mol(fCO2SW_uatm * 1e-6)
    Tsw = check.temp_K(tempSW_C + 273.15)
    Teq = check.temp_K(tempEQ_C + 273.15)
    Peq = check.pres_atm(pres_hPa / 1013.25)

    # calculate the CO2 diff due to equilibrator and seawater temperatures
    # if statement is there to save a bit of time
    if tempEQ_was_None:
        dT = 1.0
    else:
        dT = eqs.temperature_correction(Tsw, Teq)

    # a best estimate of xCO2 - this is an approximation
    # one would have to use pCO2 / Peq to get real xCO2
    # Not getting the exact equilibrator xCO2
    xCO2eq = fCO2sw * dT / Peq

    pCO2SW = fCO2sw / virial_coeff(Tsw, Peq, xCO2eq)
    pCO2SW_uatm = pCO2SW * 1e6

    return pCO2SW_uatm


def pCO2_to_fCO2(pCO2SW_uatm, tempSW_C, pres_hPa=None, tempEQ_C=None):

    """Convert pCO2 to fCO2 in sea water to account for non-ideal behaviour of CO2

    If equilibrator temperature is provided, we get a simple approximate for
    equilibrator  :math:`xCO_2` that allows for the virial expansion to be
    calculated more accurately. If not, then a simple approximation is probably
    good enough. See the examples for the differences.

    .. math::
        pCO_2^{sw} = fCO_2^{sw} \\times virial(xCO_2^{eq})

    where :math:`xCO_2^{eq} = fCO_2^{sw} \\times \\Delta T^{(sw - eq)} \\div P^{eq}`

    Args:
        fCO2SW_uatm (array): seawater fugacity of CO2 in micro atmospheres
        tempSW_C (array): sea water temperature in degrees C
        pres_hPa (array, optional): equilibrator pressure in hecto Pascals. Defaults to 1013.25.
        tempEQ_C (array, optional): equilibrator temperature in degrees C. Defaults to None.

    Returns:
        array: fugacity of CO2 in seawater

    Note:
        In FluxEngine, they account for the change in xCO2. This error is so small
        that it is not significant to be concerned about it. Their correction is
        more precise, but the difference between their iterative correction and our
        approximation is less than 1e-14 atm (or 1e-8 uatm).

    Examples:
        >>> pCO2_to_fCO2(380, 8)
        378.49789637942064
        >>> pCO2_to_fCO2(380, 8, pres_hPa=985)
        378.53967828231225
        >>> pCO2_to_fCO2(380, 8, pres_hPa=985, tempEQ_C=14)
        378.53960618459695
    """
    from . import auxiliary_equations as eqs
    from . import check_units as check

    # if equilibrator inputs are None then make defaults Patm=1, tempEQ=tempSW
    if tempEQ_C is None:
        tempEQ_C = tempSW_C
    if pres_hPa is None:
        pres_hPa = 1013.25

    # standardise the inputs and convert units
    pCO2sw = check.CO2_mol(pCO2SW_uatm * 1e-6)
    Tsw = check.temp_K(tempSW_C + 273.15)
    Teq = check.temp_K(tempEQ_C + 273.15)
    Peq = check.pres_atm(pres_hPa / 1013.25)

    # calculate the CO2 diff due to equilibrator and seawater temperatures
    dT = eqs.temperature_correction(Tsw, Teq)
    # a best estimate of xCO2 - this is an approximation
    # one would have to use pCO2 / Peq to get real xCO2
    xCO2eq = pCO2sw * dT / Peq

    fCO2sw = pCO2sw * virial_coeff(Tsw, Peq, xCO2eq)
    fCO2sw_uatm = fCO2sw * 1e6

    return fCO2sw_uatm


def virial_coeff(temp_K, pres_atm, xCO2_mol=None):
    """
    Calculate the ideal gas correction factor for converting pCO2 to fCO2.

    Based on the Lewis and Wallace 1998 Correction.

    Args:
        temp_K (array): temperature in degrees Kelvin
        pres_atm (array): atmospheric pressure in atmospheres
        xCO2_mol (array, optional): mole fraction of CO2, can
            also be p/fCO2 if xCO2 not available. Can also be None
            which makes a small difference. See examples.

    Returns:
        array: the factor to multiply/divide with pCO2/fCO2. Unitless

        .. math::
            fCO_2 = pCO_2 \\times \\text{virial expansion}

            pCO_2 = fCO_2 \\div \\text{virial expansion}

    Examples:
        From Dickson et al. (2007)

        >>> 350 * virial_coeff(298.15, 1)  # CO2 [uatm] * correction factor
        348.8836492182758

    References:
        Weiss, R. (1974). Carbon dioxide in water and seawater: the solubility of a
        non-ideal gas. Marine Chemistry, 2(3), 203–215.
        https://doi.org/10.1016/0304-4203(74)90015-2

    Compared with the Seacarb package in R
    """
    from numpy import array, exp

    from . import check_units as check

    T = check.temp_K(temp_K)
    P = check.pres_atm(pres_atm)
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
