from .constants import *
from .time_utils import swe_julday, swe_revjul, swe_deltat
from .planets import (
    swe_calc_ut,
    swe_calc,
    swe_get_ayanamsa_ut,
    swe_get_ayanamsa,
    swe_get_ayanamsa_name,
    swe_set_sid_mode,
)
from .houses import swe_houses, swe_houses_ex, swe_house_name
from .state import (
    set_topo as swe_set_topo,
    set_ephe_path as swe_set_ephe_path,
    set_ephemeris_file as swe_set_ephemeris_file,
)
from .crossing import swe_solcross_ut, swe_mooncross_ut, swe_cross_ut
from .utils import difdeg2n, swe_calc_angles
from .fixed_stars import swe_fixstar_ut
from .context import EphemerisContext  # NEW: Thread-safe context API


# =============================================================================
# PYSWISSEPH-COMPATIBLE FUNCTION ALIASES (without swe_ prefix)
# =============================================================================
# pyswisseph uses function names without the swe_ prefix
# These aliases provide 100% API compatibility with pyswisseph

# Time functions
julday = swe_julday
revjul = swe_revjul
deltat = swe_deltat

# Planet calculation
calc_ut = swe_calc_ut
calc = swe_calc

# Houses
houses = swe_houses
houses_ex = swe_houses_ex
house_name = swe_house_name

# Ayanamsa (sidereal)
get_ayanamsa_ut = swe_get_ayanamsa_ut
get_ayanamsa = swe_get_ayanamsa
get_ayanamsa_name = swe_get_ayanamsa_name
set_sid_mode = swe_set_sid_mode

# Observer location
set_topo = swe_set_topo
set_ephe_path = swe_set_ephe_path
set_ephemeris_file = swe_set_ephemeris_file

# Fixed Stars
fixstar_ut = swe_fixstar_ut

# Crossings
solcross_ut = swe_solcross_ut
mooncross_ut = swe_mooncross_ut




# Helper for Arabic parts
from .arabic_parts import calc_all_arabic_parts

# Constants (planet IDs, flags, sidereal modes)
from .constants import *

__version__ = "0.1.0"
__author__ = "Giacomo Battaglia"
__license__ = "LGPL-3.0"

__all__ = [
    # Thread-safe Context API
    "EphemerisContext",
    # Time functions (both swe_ and non-prefixed aliases)
    "swe_julday",
    "julday",
    "swe_revjul",
    "revjul",
    "swe_deltat",
    "deltat",
    # Planet calculation
    "swe_calc_ut",
    "calc_ut",
    "swe_calc",
    "calc",
    # Houses
    "swe_houses",
    "houses",
    "swe_houses_ex",
    "houses_ex",
    "swe_house_name",
    "house_name",
    # Ayanamsa (sidereal)
    "swe_set_sid_mode",
    "set_sid_mode",
    "swe_get_ayanamsa_ut",
    "get_ayanamsa_ut",
    "swe_get_ayanamsa",
    "get_ayanamsa",
    "swe_get_ayanamsa_name",
    "get_ayanamsa_name",
    # Observer location
    "swe_set_topo",
    "set_topo",
    "swe_set_ephe_path",
    "set_ephe_path",
    "swe_set_ephemeris_file",
    "set_ephemeris_file",
    # Crossings
    "swe_solcross_ut",
    "solcross_ut",
    "swe_mooncross_ut",
    "mooncross_ut",
    "swe_cross_ut",
    # Utilities
    "difdeg2n",
    "swe_calc_angles",
    # Helpers
    "calc_all_arabic_parts",
    # Fixed Stars
    "swe_fixstar_ut",
    "fixstar_ut",
]
