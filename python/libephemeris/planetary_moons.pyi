from collections.abc import Sequence

def register_moon_spk(spk_file: str, moons: Sequence[int] | None = None) -> None: ...
def unregister_moon_spk(spk_file: str) -> None: ...
def list_registered_moons() -> dict[int, str]: ...
def get_moon_name(moon_id: int) -> str: ...
def is_planetary_moon(ipl: int) -> bool: ...
def calc_moon_position(t: float, moon_id: int, iflag: int) -> tuple[float, float, float, float, float, float] | None:
    """Ecliptic-of-date position at t, a TT Julian day; angles degrees, distance AU, rates per day."""
    ...
def get_moon_coverage(moon_id: int) -> tuple[float, float] | None: ...
def close_moon_kernels() -> None: ...

from .constants import (
    FLG_BARYCTR, FLG_EQUATORIAL, FLG_HELCTR, FLG_J2000, FLG_NOABERR, FLG_NONUT, FLG_SIDEREAL, FLG_SPEED, FLG_TOPOCTR, FLG_TRUEPOS, MOON_ARIEL, MOON_CALLISTO, MOON_CHARON, MOON_DEIMOS, MOON_DIONE, MOON_ENCELADUS, MOON_EUROPA, MOON_GANYMEDE, MOON_HYPERION, MOON_IAPETUS, MOON_IO, MOON_MIMAS, MOON_MIRANDA, MOON_OBERON, MOON_OFFSET, MOON_PHOBOS, MOON_RHEA, MOON_TETHYS, MOON_TITAN, MOON_TITANIA, MOON_TRITON, MOON_UMBRIEL, NAIF_ARIEL, NAIF_CALLISTO, NAIF_CHARON, NAIF_DEIMOS, NAIF_DIONE, NAIF_ENCELADUS, NAIF_EUROPA, NAIF_GANYMEDE, NAIF_HYPERION, NAIF_IAPETUS, NAIF_IO, NAIF_MIMAS, NAIF_MIRANDA, NAIF_OBERON, NAIF_PHOBOS, NAIF_RHEA, NAIF_TETHYS, NAIF_TITAN, NAIF_TITANIA, NAIF_TRITON, NAIF_UMBRIEL, PLMOON_OFFSET,
)

from ._core import (
    _module_planetary_moons_NAIF_MARS_BARYCENTER as NAIF_MARS_BARYCENTER,
    _module_planetary_moons_NAIF_JUPITER_BARYCENTER as NAIF_JUPITER_BARYCENTER,
    _module_planetary_moons_NAIF_SATURN_BARYCENTER as NAIF_SATURN_BARYCENTER,
    _module_planetary_moons_NAIF_URANUS_BARYCENTER as NAIF_URANUS_BARYCENTER,
    _module_planetary_moons_NAIF_NEPTUNE_BARYCENTER as NAIF_NEPTUNE_BARYCENTER,
    _module_planetary_moons_NAIF_PLUTO_BARYCENTER as NAIF_PLUTO_BARYCENTER,
    _module_planetary_moons_MOON_NAIF_MAP as MOON_NAIF_MAP,
    _module_planetary_moons_MOON_NAMES as MOON_NAMES,
    _module_planetary_moons_MOON_PARENT_MAP as MOON_PARENT_MAP,
)

from libephemeris._core import (
    EphemerisRangeError,
)

from libephemeris._core import (
    get_loader,
    get_planets,
    get_timescale,
)
