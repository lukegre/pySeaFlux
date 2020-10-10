===============================
SeaFlux
===============================

.. image:: https://badgen.net/pypi/v/seaflux
        :target: https://pypi.org/project/seaflux
.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
        :target: https://www.gnu.org/licenses/gpl-3.0


Planning to publish alongside Fay, Gregor, McKinley ... (2020). DOI will be released when code has been checked/tested by various contributors. Contributors will be part of the SeaFlux package citation. Please contact me if you would like to use the package before this. 


Installing
----------

GitHub
......
``pip install https://github.com/luke-gregor/SeaFlux.git``


Overview of functionality
-------------------------

- Calculate sea-air fluxes using the bulk formulation
- Convert pCO2 to fCO2 and *vice versa*.
- Correct or adjust pCO2 for temperature changes
- Scale kw to 14C bomb values for wind products using Wanninkhof's (1992) second moment of the wind speed (requires standard deviation of the wind)
- Download NOAA Marine Boundary Layer xCO2 and related functions for pCO2 conversion
- Calculate the grid cell area (in m^2) for a grid of latitudes and longitudes - also works as an xarray method (``xda.area()``)


To Do
-----
- unit errors should only happen when more than 50% of non-nan values are not valid. Otherwise, raise warning and make output ``nan``. Will double up with ``pyCO2SYS``.
- tests! Currently there are not any meaningful tests.
- Add contributor file.
- Documentation
