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

- Calculate sea-air fluxes: bulk, rapid-equilibration model (Woolf et al. 2016)
- Convert pCO2 to fCO2 and *vice versa*. 
- Correct or adjust pCO2 for temperature changes
- Download NOAA Marine Boundary Layer xCO2 and convert to pCO2 if pressure, SST and salinity provided :code:`SeaFlux.utils.noaa_mbl_to_pCO2`

To Do
-----
- unit errors should only happen when more than 50% of non-nan values are not valid. Otherwise, raise warning and make output :code:`nan`. Will double up with :code:`pyCO2SYS`.
- tests! Currently there are not any meaningful tests. 
- Add contributor file. 
- Documentation
