import seaflux as sf


def test_CO2flux_bulk():
    flux = sf.flux_bulk(25, 35, 300, 400, 1013.25, 20)
    print(flux)

    # seaward is negative
    assert flux < 0


def test_CO2flux_woolf():
    from seaflux import flux_calculations as flux

    flux = flux.flux_woolf2016_rapid(25, 35, 300, 400, 1013.25, 20)

    # seaward is negative
    assert flux < 0
