class UnitError(Exception):
    pass


def unit_checker(arr, func, lim, vname):
    """
    Raise an error if the unit_checker are outside the given limit
    """
    from numpy import array, any

    arr = array(arr, ndmin=1)
    outside = func(arr, lim)
    if any(outside):
        fname = func.__name__
        nouts = outside.sum()
        prep = "are" if nouts > 1 else "is"
        plural = "s" if nouts > 1 else ""
        msg = (
            f"{vname} has {nouts} measurement{plural} that "
            f"{prep} {fname} than the limit ({lim}). \nUse "
            "use pd.Series.where/clip to remove/limit this data."
        )
        raise UnitError(msg)


def temp_K(temp_K):
    from numpy import less, greater

    varname = "temperature (K)"
    unit_checker(temp_K, less, 271.15, varname)
    unit_checker(temp_K, greater, 318.15, varname)


def pres_atm(pres_atm):
    from numpy import greater, less

    varname = "Pressure (atm)"
    unit_checker(pres_atm, greater, 1.5, varname)
    unit_checker(pres_atm, less, 0.5, varname)


def CO2_mol(CO2_mol):
    from numpy import greater, less

    varname = "CO2 mole fraction (ppm)"
    unit_checker(CO2_mol, greater, 0.01, varname)
    unit_checker(CO2_mol, less, 50e-6, varname)


def salt(salt):
    from numpy import greater, less

    varname = "Salinity"
    unit_checker(salt, less, 5, varname)
    unit_checker(salt, greater, 50, varname)


def wind_ms(wind_ms):
    from numpy import greater, less

    varname = "Wind speed (m/s)"
    unit_checker(wind_ms, less, 0, varname)
    unit_checker(wind_ms, greater, 40, varname)
