# Tests against PyCO2SYS (originally v1.6.0).
#
# Weiss (1974) CO2 solubility cannot be directly tested because PyCO2SYS evaluates this
# in /kg units while SeaFlux uses /l units.

import numpy as np
import PyCO2SYS as pyco2

import seaflux as sf


# Seed random number generator for reproducibility
rng = np.random.default_rng(7)
npts = 1000

# Set generic test conditions
salinity = np.array([*rng.uniform(low=0, high=50, size=npts - 1), 0])
temperature = np.array([-2, 0, *rng.uniform(low=-2, high=45, size=npts - 2)])  # deg-C
pCO2_or_fCO2 = rng.uniform(low=1, high=1000, size=npts)  # microatm
pres_atm = 1  # atm

# Convert units
temperature_K = temperature + 273.15  # K


def test_vapour_pressure_weiss1980():
    # Calculate vapour pressures of H2O
    vpf_seaflux = sf.vapour_pressure.weiss1980(salinity, temperature_K)
    vpf_pyco2 = 1 - pyco2.gas.vpfactor(temperature, salinity)
    # Compare results - should be virtually identical
    assert np.all(np.isclose(vpf_seaflux, vpf_pyco2, rtol=1e-12, atol=0))


def test_virial_coefficient():
    # Calculate fugacity factors
    fugfac_seaflux = sf.fco2_pco2_conversion.virial_coeff(temperature_K, pres_atm)
    fugfac_pyco2 = pyco2.sys(2100, 8.1, 2, 3, temperature=temperature)[
        "fugacity_factor"
    ]
    # Compare results
    assert np.all(np.isclose(fugfac_seaflux, fugfac_pyco2, rtol=0, atol=1e-7))


def test_pCO2_to_fCO2_conversion():
    # Convert pCO2 to fCO2
    fCO2_pyco2 = pyco2.sys(pCO2_or_fCO2, 8.1, 4, 3, temperature=temperature)["fCO2"]
    fCO2_seaflux = sf.pCO2_to_fCO2(pCO2_or_fCO2, temperature)
    # Compare results - not perfect agreement, but "good enough"(?) i.e. <0.01 microatm
    not_nan = ~np.isnan(fCO2_seaflux)
    assert np.all(
        np.isclose(fCO2_pyco2[not_nan], fCO2_seaflux[not_nan], rtol=0, atol=0.01)
    )


def test_fCO2_to_pCO2_conversion():
    # Convert fCO2 to pCO2
    pCO2_pyco2 = pyco2.sys(pCO2_or_fCO2, 8.1, 5, 3, temperature=temperature)["pCO2"]
    pCO2_seaflux = sf.fCO2_to_pCO2(pCO2_or_fCO2, temperature)
    # Compare results - not perfect agreement, but "good enough"(?) i.e. <0.01 microatm
    not_nan = ~np.isnan(pCO2_seaflux)
    assert np.all(
        np.isclose(pCO2_pyco2[not_nan], pCO2_seaflux[not_nan], rtol=0, atol=0.01)
    )
