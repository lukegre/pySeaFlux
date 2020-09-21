class UnitError(Exception):
    pass


def check_array_bounds(arr, lims, action='warn', name=''):
    """
    Checks that units are within the given limits. If not, then
    will raise/warn the user. Will always raise an error if more
    than half of the non-nan values are outside the limits.

    Parameters
    ----------
    arr : array-like
        The array that will be checked
    lims : tuple
        lower and upper limits of checks
        note that limits are exclusive (i.e. < and >, and not >=/<=)
    action: string
        raise - will raise an error and not continue
        warn - will throw a warning and mask values with nan
        quiet - same as warn, but without warning
        ignore - nothing will be done, but may result in bad data
    name: string
        if given, will inform the user of the name of the array
        to make debugging easier

    Return
    ------
    arr : array-like
        returns the array, but if warn or quiet, will be masked
        with nans
    """

    from numpy import array, any, nan, isnan
    from warnings import warn

    arr = array(arr, ndmin=1, dtype=float)
    if arr.size <= 2:
        return arr

    outside = (arr < lims[0]) | (arr > lims[1])

    non_nan_count = arr.size - isnan(arr).sum()
    half_outside = outside.sum() > (non_nan_count * 0.5)
    if half_outside:
        raise UnitError(
            f"More than half of the values in {name} are outside the limits "
            f"{str(lims)}. Check that input contains the correct units.")

    msg = (f"There are {outside.sum():d} values that do not fall within "
           f"the given limits {str(lims)}"
           f" of {name}" if name != "" else "")

    if any(outside) & (action == 'raise'):
        raise UnitError(msg)
    elif action == 'warn':
        if any(outside):
            warn(msg, Warning)
        arr[outside] = nan
    elif action == 'quiet':
        arr[outside] = nan
    elif action == 'ignore':
        pass
    else:
        raise Exception("action must have raise/warn/quiet/ignore as inputs")

    return arr


def temp_K(temp_K):
    return check_array_bounds(
        arr=temp_K,
        lims=(270, 318.5),
        action='warn',
        name="temperature (K)"
    )


def pres_atm(pres_atm):
    return check_array_bounds(
        arr=pres_atm,
        lims=(0.5, 1.5),
        action='warn',
        name="Pressure (atm)"
    )


def CO2_mol(CO2_mol):
    return check_array_bounds(
        arr=CO2_mol,
        lims=(5e-6, 0.08),
        action='warn',
        name="CO2 mole fraction (ppm)"
    )


def salt(salt):
    return check_array_bounds(
        arr=salt,
        lims=(0, 50),
        action='warn',
        name="Salinity (PSU)"
    )


def wind_ms(wind_ms):
    return check_array_bounds(
        arr=wind_ms,
        lims=(0, 50),
        action='warn',
        name="Wind speed (m/s)"
    )
