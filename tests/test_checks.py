def test_limits_warn():
    from seaflux.check_units import check_array_bounds
    import numpy as np

    arr = np.random.normal(0, 2, size=100)
    n_outside = (arr < -2).sum() + (arr > 2).sum()

    out = check_array_bounds(arr, lims=(-2, 2), action="quiet", name="test array")
    n_nans = np.isnan(out).sum()

    assert n_nans == n_outside


def test_limits_raise():
    from seaflux.check_units import check_array_bounds, UnitError
    import numpy as np

    arr = np.random.normal(0, 2, size=100)
    try:
        check_array_bounds(arr, lims=(-2, 2), action="raise", name="test array")
    except UnitError:
        pass
