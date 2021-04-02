import seaflux as sf


def test_CO2flux_bulk():
    flux = sf.flux_bulk(25, 35, 300, 400, 1013.25, 20)
    print(flux)

    # seaward is negative
    assert flux < 0
