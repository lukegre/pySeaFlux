"""
Additional equations
--------------------

Not necessarily linked to the marine carbonate system, but are
useful.
"""

import numpy as np


def pressure_height_correction(pres_hPa, tempSW_C, sensor_height=10.0, checks=True):
    """Returns exact sea level pressure if the sensor is measuring at height

    Args:
        pres_hPa (array): Pressure in kiloPascal measured at height
        tempSW_C (array): Temperature of the seawater in deg C
        sensor_height (float, optional): the height of the sensor above sea
            level. Can be negative if you want to convert SLP to sensor height
            pressure. Defaults to 10.0.

    Returns:
        array: height corrected pressure
    """
    from numpy import nanmedian

    if checks:
        if nanmedian(tempSW_C) > 270:
            raise ValueError("Temperature is not in Celsius")
        if nanmedian(pres_hPa) < 10:
            raise ValueError("Pressure is not in hPa")

    T = tempSW_C + 273.15  # temperature in Kelvin
    P = pres_hPa * 100  # pressure in Pascal

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


def temperature_correction(temp_in, temp_out):
    """pCO2 correction factor for temperature changes

    Calculate a correction factor for the temperature difference between the
    intake and equilibrator. This is based on the empirical relationship used
    in Takahashi et al. 1993.

    .. math::

        pCO_2^{Tout} = pCO_2^{Tin} * T^{factor}

    Args:
        temp_in (array): temperature at which original pCO2 is measured (degK or degC)
        temp_out (array): temperature for which pCO2 should be represented

    Returns:
        array: a correction factor to be multiplied to pCO2 (unitless)

    References:
        Takahashi, Taro et al. (1993). Seasonal variation of CO2 and nutrients in
        the high-latitude surface oceans: A comparative study. Global
        Biogeochemical Cycles, 7(4), 843–878. https://doi.org/10.1029/93GB02263
    """
    # see the Takahashi 1993 paper for full description

    Ti = temp_in
    To = temp_out

    factor = np.exp(0.0433 * (To - Ti) - 4.35e-05 * (To ** 2 - Ti ** 2))

    return factor
