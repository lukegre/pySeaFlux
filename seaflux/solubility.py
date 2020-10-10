from . import check_units as check
from . import vapour_pressure as vapress
from .utils import preserve_xda


@preserve_xda
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
    from Weiss (1974) Table 2 but with pH2O correction
    >>> solubility_weiss1974(35, 299.15)
    0.029285284543519093
    """

    from numpy import exp, log

    T = check.temp_K(temp_K)
    S = check.salt(salt)
    P = check.pres_atm(press_atm)

    # from table in Wanninkhof 2014
    a1 = -58.0931
    a2 = +90.5069
    a3 = +22.2940
    b1 = +0.027766
    b2 = -0.025888
    b3 = +0.0050578

    T100 = T / 100
    K0 = exp(
        a1 + a2 * (100 / T) + a3 * log(T100) + S * (b1 + b2 * T100 + b3 * T100 ** 2)
    )

    pH2O = vapress.weiss1980(S, T)
    K0 = K0 / (P - pH2O)
    meta = {
        "description": "CO2 solubility in seawater using the formulation of Weiss 1974",
        "units": "mol/L/atm",
        "long_name": "CO2 solubility in seawater",
    }

    return K0, meta  # units mol/L/atm


def solubility_woolf2016(salt, temp_K, deltaT, press_atm=1):
    """
    A wrapper around solubility calculated using the Weiss (1974) approach.
    From FluxEngine (Shutler et al, 2016; and Holding et al, 2019).

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
