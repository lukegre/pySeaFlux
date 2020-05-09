===============================
seaflux
===============================


.. image:: https://img.shields.io/travis/luke-gregor/SeaFlux.svg
        :target: https://travis-ci.org/luke-gregor/SeaFlux
.. image:: https://badgen.net/pypi/v/seaflux
        :target: https://pypi.org/project/seaflux
.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
        :target: https://www.gnu.org/licenses/gpl-3.0


Calculate sea-air fluxes


Now has function to Calculate NOAA Marine Boundary Layer pCO2: :code:`SeaFlux.utils.noaa_mbl_to_pCO2`

To Do
-----
- unit errors should only happen when more than 50% of non-nan values are not valid. Otherwise, raise warning and make output :code:`nan`. Will double up with :code:`pyCO2SYS`.
- tests! Currently there are not any meaningful tests. 
- Add contributor file. 
- Documentation
