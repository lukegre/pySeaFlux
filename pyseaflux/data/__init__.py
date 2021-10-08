"""
This module is only intended to be by the SeaFlux authors to
download the data required to create the SeaFlux ensemble.
Has links to most data sources (ERA5 might not be included)

Hence, this module is not imported by default and submodules
should be imported on demand.
"""

from .download_zenodo_files import get_seaflux_data, get_zenodo_catalog
