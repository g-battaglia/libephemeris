"""
Shared fixtures and utilities for LEB vs Skyfield comparison tests.

This module provides:
- TierTolerances dataclass with per-tier configurable thresholds
- CompareHelper class for executing functions in both modes
- Test date generators
- Helper functions for angular calculations
- Body lists and per-body tolerance tables

Tolerance override via environment variables:
    Per-tier:  LEB_TOL_BASE_POSITION_ARCSEC=0.3
    Global:    LEB_TOL_POSITION_ARCSEC=2.0   (fallback for all tiers)
"""

from __future__ import annotations

import dataclasses
import os
from typing import Any, Callable, Generator

import pytest

import libephemeris as ephem
from libephemeris.time_utils import julday


# =============================================================================
# TIER-AWARE TOLERANCES
# =============================================================================

_LEB_DATA_DIR = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "data", "leb")
)


@dataclasses.dataclass
class TierTolerances:
    """Configurable tolerance thresholds for LEB vs Skyfield comparison.

    All values can be overridden via environment variables:
      - Per-tier: LEB_TOL_{TIER}_{FIELD}  e.g. LEB_TOL_BASE_POSITION_ARCSEC=0.3
      - Global:   LEB_TOL_{FIELD}         e.g. LEB_TOL_POSITION_ARCSEC=2.0

    Per-tier overrides take precedence over global overrides.

    Typical precision by body type:
    - Inner planets (Mercury-Venus): ~0.1-0.2 arcsec
    - Outer planets (Mars-Pluto): ~0.1-5 arcsec (varies by date)
    - Moon: ~0.1 arcsec position, ~0.01 deg/day speed
    - Asteroids: ~0.1 arcsec
    - Ecliptic bodies (nodes, Lilith): ~0.5 arcsec
    """

    # Core position
    POSITION_ARCSEC: float = 5.0
    ECLIPTIC_ARCSEC: float = 0.5
    ASTEROID_ARCSEC: float = 5.0
    EQUATORIAL_ARCSEC: float = 5.0
    J2000_ARCSEC: float = 5.0
    SIDEREAL_ARCSEC: float = 20.0
    DISTANCE_AU: float = 1e-4

    # Velocity
    SPEED_LON_DEG_DAY: float = 0.05
    SPEED_LAT_DEG_DAY: float = 0.05
    SPEED_DIST_AU_DAY: float = 1e-4

    # Asteroid velocity (architectural: ICRS→ecliptic pipeline amplification)
    ASTEROID_SPEED_LON_DEG_DAY: float = 0.15
    ASTEROID_SPEED_LAT_DEG_DAY: float = 1.0
    ASTEROID_SPEED_DIST_AU_DAY: float = 5e-3

    # Ayanamsha
    AYANAMSHA_ARCSEC: float = 0.001

    # Nutation & Delta-T
    NUTATION_ARCSEC: float = 0.01
    DELTAT_SEC: float = 0.1

    # Timing (indirect callers)
    CROSSING_SUN_SEC: float = 1.0
    CROSSING_MOON_SEC: float = 5.0
    CROSSING_MOON_NODE_SEC: float = 10.0
    CROSSING_PLANET_SEC: float = 30.0
    ECLIPSE_TIMING_SEC: float = 1.0
    ECLIPSE_MAGNITUDE: float = 1e-4
    ECLIPSE_POSITION_DEG: float = 0.001
    ELONGATION_ARCSEC: float = 0.5
    STATION_TIMING_SEC: float = 1.0
    RISE_TRANSIT_SEC: float = 1.0
    HOUSE_SUNSHINE_ARCSEC: float = 0.001
    GAUQUELIN_SECTOR: float = 0.001

    @classmethod
    def for_tier(cls, tier: str, **overrides: float) -> "TierTolerances":
        """Create tolerances for a specific tier with env-var overrides.

        Resolution order (highest priority first):
        1. Explicit overrides passed as kwargs
        2. LEB_TOL_{TIER}_{FIELD} environment variable
        3. LEB_TOL_{FIELD} environment variable (global)
        4. Default value from TIER_DEFAULTS[tier]
        """
        defaults = TIER_DEFAULTS.get(tier, {})
        kwargs: dict[str, float] = {}

        for field in dataclasses.fields(cls):
            name = field.name

            # Check explicit override first
            if name in overrides:
                kwargs[name] = overrides[name]
                continue

            # Check per-tier env var
            tier_env = f"LEB_TOL_{tier.upper()}_{name.upper()}"
            tier_val = os.environ.get(tier_env)
            if tier_val is not None:
                kwargs[name] = float(tier_val)
                continue

            # Check global env var
            global_env = f"LEB_TOL_{name.upper()}"
            global_val = os.environ.get(global_env)
            if global_val is not None:
                kwargs[name] = float(global_val)
                continue

            # Use tier-specific default (or dataclass default)
            if name in defaults:
                kwargs[name] = defaults[name]

        return cls(**kwargs)

    @classmethod
    def from_env(cls) -> "TierTolerances":
        """Load tolerance overrides from environment variables (legacy)."""
        return cls.for_tier("medium")


# Per-tier default overrides (only fields that differ from dataclass defaults)
TIER_DEFAULTS: dict[str, dict[str, float]] = {
    "base": {
        # ICRS barycentric storage with COORD_ICRS_BARY_SYSTEM for outer
        # planets + gravitational deflection + SR aberration pipeline.
        # Measured worst-case (200-point sweep, 1849-2150):
        #   Sun 0.000001", Moon 0.000332", Mercury 0.000004",
        #   Venus 0.000006", Mars 0.000003", Jupiter 0.000002",
        #   Saturn 0.000005", Uranus 0.000001", Neptune 0.000000",
        #   Pluto 0.000000", Earth 0.000000"
        #   Asteroids: max 0.000045" (Juno)
        #   Ecliptic: max 0.000049" (OscuApog)
        "POSITION_ARCSEC": 0.001,  # Moon 0.000332" (3x margin)
        "ASTEROID_ARCSEC": 0.001,  # Juno 0.000045" (22x margin)
        "EQUATORIAL_ARCSEC": 0.02,  # Moon heliocentric 0.0103" (2x margin)
        "J2000_ARCSEC": 0.001,  # Same pipeline, J2000 ecliptic output
        "SIDEREAL_ARCSEC": 0.001,  # = position error (ayanamsha is formula-exact)
        "ECLIPTIC_ARCSEC": 0.001,  # OscuApog 0.000049" (20x margin)
        "DISTANCE_AU": 5e-6,  # Uranus helio 1.04e-6 AU (5x margin)
        "SPEED_LON_DEG_DAY": 0.045,  # OscuApogee 0.035 (1.3x margin)
        "SPEED_LAT_DEG_DAY": 0.004,  # OscuApogee ~0.003, planets 0.000082
        "SPEED_DIST_AU_DAY": 1.2e-4,  # Pluto 8.98e-5 (1.3x margin)
        "ASTEROID_SPEED_LON_DEG_DAY": 0.15,  # ICRS pipeline amplification
        "ASTEROID_SPEED_LAT_DEG_DAY": 1.7,  # Pallas 0.84 (2x margin)
        "ASTEROID_SPEED_DIST_AU_DAY": 5e-3,  # ICRS pipeline
    },
    "medium": {
        # ICRS barycentric storage with COORD_ICRS_BARY_SYSTEM for outer
        # planets + gravitational deflection + SR aberration pipeline.
        # Measured worst-case (200-point sweep, 1550-2650):
        #   Sun 0.000001", Moon 0.000325", Mercury 0.000004",
        #   Venus 0.000004", Mars 0.000002", Jupiter 0.000007",
        #   Saturn 0.000001", Uranus 0.000005", Neptune 0.000000",
        #   Pluto 0.000000", Earth 0.000000"
        #   Asteroids: max 0.000036" (Vesta)
        #   Ecliptic: max 0.000075" (OscuApog)
        "POSITION_ARCSEC": 0.001,  # Moon 0.000325" (3x margin)
        "ASTEROID_ARCSEC": 0.001,  # Vesta 0.000036" (28x margin)
        "EQUATORIAL_ARCSEC": 0.02,  # Moon heliocentric amplification
        "J2000_ARCSEC": 0.001,  # Same pipeline, J2000 ecliptic output
        "SIDEREAL_ARCSEC": 0.001,  # = position error (ayanamsha is formula-exact)
        "ECLIPTIC_ARCSEC": 0.001,  # OscuApog 0.000075" (13x margin)
        "DISTANCE_AU": 5e-6,  # Uranus 1.36e-08 AU (helio amplification margin)
        # Velocity — OscuApogee dominates lon (0.043 deg/day).
        "SPEED_LON_DEG_DAY": 0.045,  # OscuApogee 0.043
        "SPEED_LAT_DEG_DAY": 0.004,  # OscuApogee ~0.003, planets 0.000052
        "SPEED_DIST_AU_DAY": 1e-4,  # Pluto 4.87e-5 AU/day (2.1x margin)
        # Asteroid velocity — ICRS pipeline amplification.
        "ASTEROID_SPEED_LON_DEG_DAY": 0.15,  # ICRS pipeline
        "ASTEROID_SPEED_LAT_DEG_DAY": 1.7,  # Pallas (2x margin)
        "ASTEROID_SPEED_DIST_AU_DAY": 5e-3,  # ICRS pipeline
    },
    "extended": {
        # Extended tier (de441, -5000 to 5000 CE, 10,000 years).
        # ICRS barycentric storage with COORD_ICRS_BARY_SYSTEM for outer
        # planets, gravitational deflection + SR aberration pipeline.
        # Measured precision: all planets <0.001", asteroids <0.001".
        # Ecliptic bodies kept at 0.1" due to Meeus polynomial degradation
        # beyond ±20 centuries from J2000.
        "POSITION_ARCSEC": 0.001,  # Measured <0.0003" for all planets
        "EQUATORIAL_ARCSEC": 0.02,  # Heliocentric amplification
        "J2000_ARCSEC": 0.001,  # Same pipeline
        "SIDEREAL_ARCSEC": 0.001,  # POSITION + ayanamsha
        "ECLIPTIC_ARCSEC": 0.1,  # OscuApogee at extreme dates (Meeus limit)
        "ASTEROID_ARCSEC": 0.001,  # Measured <0.000018"
        "DISTANCE_AU": 5e-6,  # Matching base/medium
        "SPEED_LON_DEG_DAY": 0.05,  # OscuApogee 0.035 (1.4x margin)
        "SPEED_LAT_DEG_DAY": 0.005,  # OscuApogee 0.002 (2.5x margin)
        "SPEED_DIST_AU_DAY": 1.2e-4,  # Pluto (matching base tier)
        "ASTEROID_SPEED_LON_DEG_DAY": 0.15,  # ICRS pipeline
        "ASTEROID_SPEED_LAT_DEG_DAY": 1.7,  # ICRS pipeline
        "ASTEROID_SPEED_DIST_AU_DAY": 5e-3,  # ICRS pipeline
        "NUTATION_ARCSEC": 0.01,
        "DELTAT_SEC": 0.1,
    },
}

# Legacy global instance (medium tier, for existing tests)
TOLS = TierTolerances.from_env()


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


def angular_diff(a: float, b: float) -> float:
    """Angular difference accounting for 360 wrap."""
    d = abs(a - b)
    if d > 180.0:
        d = 360.0 - d
    return d


def lon_error_arcsec(a: float, b: float) -> float:
    """Longitude error in arcseconds with wrap-around."""
    return angular_diff(a, b) * 3600.0


def generate_test_dates(
    n: int, jd_start: float, jd_end: float, margin: float = 30.0
) -> list[float]:
    """Generate n uniformly-spaced dates within [start+margin, end-margin]."""
    effective_start = jd_start + margin
    effective_end = jd_end - margin
    if n <= 1:
        return [(effective_start + effective_end) / 2.0]
    step = (effective_end - effective_start) / (n - 1)
    return [effective_start + i * step for i in range(n)]


def year_to_jd(year: int) -> float:
    """Convert a year to Julian Day (January 1.0)."""
    return julday(year, 1, 1, 0.0)


def leb_file_path(tier: str) -> str:
    """Return absolute path to a tier's LEB file."""
    filename = f"ephemeris_{tier}.leb"
    return os.path.join(_LEB_DATA_DIR, filename)


# =============================================================================
# COMPARE HELPER
# =============================================================================


class CompareHelper:
    """Execute functions in both Skyfield and LEB mode, saving/restoring state.

    Note: Uses set_calc_mode("skyfield") to prevent LEB auto-discovery from
    interfering with Skyfield reference calculations.

    Args:
        leb_path: Path to the LEB file for this tier.
        tier: Precision tier name ("base", "medium", "extended"). When set,
            skyfield() will call set_precision_tier() so the correct BSP file
            (de440s/de440/de441) is used for reference calculations.
    """

    def __init__(self, leb_path: str, tier: str = "medium"):
        self.leb_path = leb_path
        self.tier = tier
        self._saved_leb: str | None = None
        self._saved_mode: str | None = None
        self._saved_reader: Any = None
        self._saved_tier: str | None = None

    def setup(self) -> None:
        """Save current global state and enable SPK auto-download."""
        self._saved_leb = ephem.state._LEB_FILE
        self._saved_mode = ephem.state._CALC_MODE
        self._saved_reader = ephem.state._LEB_READER
        self._saved_tier = ephem.state._PRECISION_TIER
        self._saved_auto_spk = ephem.state._AUTO_SPK_DOWNLOAD
        # Enable auto SPK download so Skyfield reference calculations
        # use SPK21 data for asteroids instead of Keplerian fallback.
        ephem.set_auto_spk_download(True)

    def teardown(self) -> None:
        """Restore saved global state."""
        ephem.state._LEB_FILE = self._saved_leb
        ephem.state._LEB_READER = self._saved_reader
        if self._saved_tier is not None:
            ephem.set_precision_tier(self._saved_tier)
        else:
            ephem.state._PRECISION_TIER = None
        if self._saved_mode is not None:
            ephem.set_calc_mode(self._saved_mode)
        else:
            ephem.set_calc_mode(None)
        ephem.set_auto_spk_download(self._saved_auto_spk)

    def skyfield(self, fn: Callable, *args: Any, **kwargs: Any) -> Any:
        """Call fn in forced Skyfield mode (no LEB auto-discovery).

        Sets the precision tier so the correct BSP file is loaded (e.g.
        de441.bsp for extended tier, de440s.bsp for base tier).
        """
        ephem.state._LEB_FILE = None
        ephem.set_precision_tier(self.tier)
        ephem.set_calc_mode("skyfield")
        try:
            return fn(*args, **kwargs)
        finally:
            ephem.set_calc_mode(None)

    def leb(self, fn: Callable, *args: Any, **kwargs: Any) -> Any:
        """Call fn in LEB mode with explicit file path."""
        ephem.state._LEB_FILE = self.leb_path
        ephem.set_calc_mode("auto")
        try:
            return fn(*args, **kwargs)
        finally:
            ephem.set_calc_mode(None)

    def leb_only(self, fn: Callable, *args: Any, **kwargs: Any) -> Any:
        """Call fn through the source-sealed LEB vector boundary."""
        ephem.set_leb_file(self.leb_path)
        ephem.set_calc_mode("leb")
        try:
            return fn(*args, **kwargs)
        finally:
            ephem.set_calc_mode(None)


# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture(scope="session")
def leb_file() -> str:
    """Path to the pre-generated medium tier LEB file."""
    path = os.path.abspath(leb_file_path("medium"))
    if not os.path.exists(path):
        pytest.skip(f"LEB file not found: {path}")
    return path


@pytest.fixture(scope="session")
def leb_file_base() -> str:
    """Path to the pre-generated base tier LEB file."""
    path = os.path.abspath(leb_file_path("base"))
    if not os.path.exists(path):
        pytest.skip(f"LEB base file not found: {path}")
    return path


@pytest.fixture(scope="session")
def leb_file_extended() -> str:
    """Path to the pre-generated extended tier LEB file."""
    path = os.path.abspath(leb_file_path("extended"))
    if not os.path.exists(path):
        pytest.skip(f"LEB extended file not found: {path}")
    return path


@pytest.fixture
def compare(leb_file: str) -> Generator[CompareHelper, None, None]:
    """CompareHelper instance for medium tier with automatic state save/restore."""
    helper = CompareHelper(leb_file, tier="medium")
    helper.setup()
    try:
        yield helper
    finally:
        helper.teardown()


@pytest.fixture
def compare_base(leb_file_base: str) -> Generator[CompareHelper, None, None]:
    """CompareHelper instance for base tier."""
    helper = CompareHelper(leb_file_base, tier="base")
    helper.setup()
    try:
        yield helper
    finally:
        helper.teardown()


@pytest.fixture
def compare_extended(leb_file_extended: str) -> Generator[CompareHelper, None, None]:
    """CompareHelper instance for extended tier."""
    helper = CompareHelper(leb_file_extended, tier="extended")
    helper.setup()
    try:
        yield helper
    finally:
        helper.teardown()


# Date fixtures for medium tier (different test intensities)
_MEDIUM_TIER_START = year_to_jd(1560)
_MEDIUM_TIER_END = year_to_jd(2640)


@pytest.fixture(scope="session")
def test_dates_200() -> list[float]:
    """200 uniformly-spaced JDs across 1560-2640."""
    return generate_test_dates(200, _MEDIUM_TIER_START, _MEDIUM_TIER_END)


@pytest.fixture(scope="session")
def test_dates_100() -> list[float]:
    """100 uniformly-spaced JDs across 1560-2640."""
    return generate_test_dates(100, _MEDIUM_TIER_START, _MEDIUM_TIER_END)


@pytest.fixture(scope="session")
def test_dates_50() -> list[float]:
    """50 uniformly-spaced JDs across 1560-2640."""
    return generate_test_dates(50, _MEDIUM_TIER_START, _MEDIUM_TIER_END)


@pytest.fixture(scope="session")
def test_dates_20() -> list[float]:
    """20 uniformly-spaced JDs across 1560-2640."""
    return generate_test_dates(20, _MEDIUM_TIER_START, _MEDIUM_TIER_END)


# =============================================================================
# BODY LISTS
# =============================================================================

MAIN_PLANETS = [
    (0, "Sun"),
    (1, "Moon"),
    (2, "Mercury"),
    (3, "Venus"),
    (4, "Mars"),
    (5, "Jupiter"),
    (6, "Saturn"),
    (7, "Uranus"),
    (8, "Neptune"),
    (9, "Pluto"),
    (14, "Earth"),
]

# Backward-compatible alias.
ICRS_PLANETS = MAIN_PLANETS

ECLIPTIC_BODIES = [
    (10, "MeanNode"),
    (11, "TrueNode"),
    (12, "MeanApogee"),
    (13, "OscuApogee"),
    (21, "InterpApogee"),
    (22, "InterpPerigee"),
]

ASTEROID_BODIES = [
    (15, "Chiron"),
    (17, "Ceres"),
    (18, "Pallas"),
    (19, "Juno"),
    (20, "Vesta"),
]

ALL_LEB_BODIES = MAIN_PLANETS + ECLIPTIC_BODIES + ASTEROID_BODIES

# Asteroid SPK coverage (same for all tiers — JPL Horizons limitation).
# Outside this range, both LEB and Skyfield fall back to Keplerian orbits
# which are catastrophically wrong.  Test dates must be filtered to this range.
# Asteroid SPK actual coverage: ~1900-2100 CE for all 5 asteroids.
# The JPL Horizons SPK21 files have limited temporal coverage.  Outside this
# range, Skyfield falls back to Keplerian orbital elements which are
# catastrophically wrong.  The LEB generator used Skyfield during generation,
# so LEB segments outside SPK coverage contain Keplerian-fallback data that
# is baked into the Chebyshev coefficients.  This contamination extends well
# beyond the SPK boundary due to Chebyshev fitting windows straddling the edge.
# Measured errors: 22,000-97,000 arcsec position, 0.06-0.59 AU distance.
_ASTEROID_SPK_JD_START = year_to_jd(
    1930
)  # Safe SPK coverage start (Horizons SPK21 ~1925)
_ASTEROID_SPK_JD_END = year_to_jd(2080)  # Safe SPK coverage end (20yr margin)
_ASTEROID_BODY_IDS = {b[0] for b in ASTEROID_BODIES}


def filter_asteroid_dates(
    dates: list[float], body_id: int, margin: float = 365.0
) -> list[float]:
    """Filter test dates to asteroid SPK coverage range.

    For non-asteroid bodies, returns dates unchanged.
    For asteroids, excludes dates outside the SPK coverage range
    (with a safety margin to avoid Chebyshev segments contaminated
    by boundary extrapolation).
    """
    if body_id not in _ASTEROID_BODY_IDS:
        return dates
    jd_lo = _ASTEROID_SPK_JD_START + margin
    jd_hi = _ASTEROID_SPK_JD_END - margin
    return [jd for jd in dates if jd_lo <= jd <= jd_hi]


# Ecliptic body tolerances (per-body, in arcsec / deg-day)
# Measured worst-case (base tier, 200 samples):
#   MeanNode 0.000000", TrueNode 0.000001", MeanApog 0.000001",
#   OscuApog 0.000049"
# IntpApog/IntpPerg: LEB data predates apse correction tables, so Chebyshev
# fits diverge from current Skyfield reference by up to ~1°.
# After LEB regeneration with 1d/deg17 params, tighten to <1".
ECLIPTIC_TOLERANCES = {
    10: {"lon": 0.001, "speed": 0.0001},  # Mean Node
    11: {"lon": 0.001, "speed": 0.01},  # True Node
    12: {"lon": 0.001, "speed": 0.0001},  # Mean Apogee
    13: {"lon": 0.001, "speed": 0.05},  # Oscu Apogee (higher speed variance)
    21: {
        "lon": 3600.0,
        "lat": 36000.0,
        "dist": 0.001,
        "speed": 1.0,
    },  # Interp Apogee (pre-regen)
    22: {
        "lon": 7200.0,
        "lat": 36000.0,
        "dist": 0.001,
        "speed": 1.0,
    },  # Interp Perigee (pre-regen)
}

# Broad predefined-mode sweep plus the caller-defined mode. Catalog-direction
# modes are kept in the additional sweep below because they exercise a
# different dependency path.
SIDEREAL_MODE_SWEEP = [
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
    18,
    19,
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    37,
    38,
    41,
    43,
    44,
    45,
    46,
    255,  # user-defined
]

# Additional public modes, including the catalog-direction implementations.
ADDITIONAL_SIDEREAL_MODES = [17, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 39, 40, 42]
