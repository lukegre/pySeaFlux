from pkg_resources import DistributionNotFound, get_distribution

from .core import fCO2_to_pCO2, flux_bulk, flux_woolf2016_rapid, pCO2_to_fCO2
from . import aux_eqs as eqs
from . import gas_transfer_CO2
from . import utils

try:
    __version__ = get_distribution("glidertools").version
except DistributionNotFound:
    __version__ = "version_undefined"
del get_distribution, DistributionNotFound
