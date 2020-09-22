===============================
seaflux
===============================


.. image:: https://travis-ci.com/luke-gregor/SeaFlux.svg?branch=master
   :target: https://travis-ci.com/luke-gregor/SeaFlux
.. image:: https://badgen.net/pypi/v/seaflux
   :target: https://pypi.org/project/seaflux
.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0

**WARNING: STILL IN DEVELOPMENT**

With the following functionality

- Calculate sea-air fluxes using the bulk formulation
- Convert pCO2 to fCO2 and *vice versa*.
- Correct or adjust pCO2 for temperature changes
- Scale :math:`k_w` to 14C bomb values for wind products using Wanninkhof's (1992) second moment of the wind speed (requires standard deviation of the wind)
- Download NOAA Marine Boundary Layer xCO2 and related functions for pCO2 conversion
- Calculate the grid cell area (in :math:`m^2`) for a grid of latitudes and longitudes - also works as an xarray method (`xda.area()`)


To Do
-----
- unit errors should only happen when more than 50% of non-nan values are not valid. Otherwise, raise warning and make output :code:`nan`. Will double up with :code:`pyCO2SYS`.
- tests! Currently there are not any meaningful tests.
- Add contributor file.
- Documentation
