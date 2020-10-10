import seaflux as sf


def test_fCO2_to_pCO2_full():
    known_val = 381.466027968504
    assert sf.fCO2_to_pCO2(380, 8, pres_hPa=985, tempEQ_C=14) == known_val


def test_fCO2_to_pCO2_wo_press_tempEQ():
    assert sf.fCO2_to_pCO2(380, 8) == 381.50806485658234


def test_fCO2_to_pCO2_wo_tempEQ():
    assert sf.fCO2_to_pCO2(380, 8, pres_hPa=985) == 381.4659553134281


def test_pCO2_to_fCO2_wo_press_tempEQ():
    from seaflux import pCO2_to_fCO2

    known_value = 378.49789637942064
    assert pCO2_to_fCO2(380, 8) == known_value


def test_pCO2_to_fCO2_wo_tempEQ():
    from seaflux import pCO2_to_fCO2

    known_value = 378.53967828231225
    assert pCO2_to_fCO2(380, 8, pres_hPa=985) == known_value


def test_pCO2_to_fCO2_full():
    from seaflux import pCO2_to_fCO2

    known_value = 378.53960618459695
    assert pCO2_to_fCO2(380, 8, pres_hPa=985, tempEQ_C=14) == known_value
