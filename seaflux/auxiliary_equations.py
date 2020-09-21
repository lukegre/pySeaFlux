from . import check_units as check


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

    T = check.temp_K(tempSW_C + 273.15)  # temperature in Kelvin
    P = check.pres_atm(pres_hPa / 1013.25)  # pressure in Pascal

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
