# Tests against PyCO2SYS (originally v1.6.0).
#
# Weiss (1974) CO2 solubility cannot be directly tested because PyCO2SYS evaluates this
# in /kg units while SeaFlux uses /l units.

import seaflux as sf, PyCO2SYS as pyco2, numpy as np

# Seed random number generator for reproducibility
rng = np.random.default_rng(7)
npts = 1000


def test_vapour_pressure_weiss1980():
    # Set test conditions
    salinity = np.array([*rng.uniform(low=0, high=50, size=npts), 0, 35])
    temperature = np.array([-2, 0, *rng.uniform(low=-2, high=45, size=npts)])
    # Convert units
    temperature_K = temperature + 273.15
    # Calculate vapour pressures of H2O
    vpf_seaflux = sf.vapour_pressure.weiss1980(salinity, temperature_K)
    vpf_pyco2 = 1 - pyco2.gas.vpfactor(temperature, salinity)
    # Compare results - should be virtually identical
    assert np.all(np.isclose(vpf_seaflux, vpf_pyco2, rtol=1e-12, atol=0))
