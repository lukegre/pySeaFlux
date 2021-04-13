"""
Set up module access for the base package
"""
from . import auxiliary_equations as eqs
from . import gas_transfer_velocity as kw
from .area import get_area_from_dataset
from .data.download_zenodo_files import get_seaflux_data, get_zenodo_catalog
from .data.utils import preprocess as preprocess
from .fco2_pco2_conversion import fCO2_to_pCO2, pCO2_to_fCO2
from .flux_calculations import flux_bulk
