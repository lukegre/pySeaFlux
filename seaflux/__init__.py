from pkg_resources import DistributionNotFound, get_distribution

from . import atmospheric as atm
from . import auxiliary_equations as eqs
from . import gas_transfer_velocity
from . import utils

from . fco2_pco2_conversion import fCO2_to_pCO2, pCO2_to_fCO2
from . flux_calculations import flux_bulk
from . utils import area_grid

from . import xr_accessor as _accessor

try:
    __version__ = get_distribution("glidertools").version
except DistributionNotFound:
    __version__ = "version_undefined"
del get_distribution, DistributionNotFound
