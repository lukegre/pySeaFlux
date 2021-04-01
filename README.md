SeaFlux
==============================
[![Build Status](https://github.com/lukegre/seaflux/workflows/Tests/badge.svg)](https://github.com/lukegre/seaflux/actions)
[![Documentation Status](https://readthedocs.org/projects/seaflux/badge/?version=latest)](https://seaflux.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/lukegre/seaflux/branch/master/graph/badge.svg)](https://codecov.io/gh/lukegre/seaflux)
[![License:MIT](https://img.shields.io/badge/License-MIT-lightgray.svg?style=flt-square)](https://opensource.org/licenses/MIT)
[![pypi](https://img.shields.io/pypi/v/seaflux.svg)](https://pypi.org/project/seaflux)
<!-- [![conda-forge](https://img.shields.io/conda/dn/conda-forge/seaflux?label=conda-forge)](https://anaconda.org/conda-forge/seaflux) -->


Calculation of air-sea fluxes


Installing
----------

### GitHub
`pip install git+https://github.com/luke-gregor/SeaFlux.git`

### PyPi
`pip install seaflux`


Overview of functionality
-------------------------

- Calculate sea-air fluxes using the bulk formulation
- Convert pCO2 to fCO2 and *vice versa*.
- Correct or adjust pCO2 for temperature changes
- Scale kw to 14C bomb values for wind products using Wanninkhof's (1992) second moment of the wind speed (requires standard deviation of the wind)
- Download NOAA Marine Boundary Layer xCO2 and related functions for pCO2 conversion
- Calculate the grid cell area (in m^2) for a grid of latitudes and longitudes - also works as an xarray method (`xda.area()`)


To Do
-----
- unit errors should only happen when more than 50% of non-nan values are not valid. Otherwise, raise warning and make output `nan`. Will double up with `pyCO2SYS`.
- tests! Currently there are not any meaningful tests.
- Add contributor file.
- Documentation


--------

<p><small>Project based on the <a target="_blank" href="https://github.com/jbusecke/cookiecutter-science-project">cookiecutter science project template</a>.</small></p>
