"""
pytest configuration and shared fixtures for LibEphemeris tests.

This module provides test infrastructure for libephemeris standalone tests
(no external reference-ephemeris dependency required).
"""

import logging

import pytest
import random
import libephemeris as ephem
from libephemeris.constants import *


# ============================================================================
# GLOBAL SETUP - Suppress debug log noise during tests
# ============================================================================


@pytest.fixture(autouse=True)
def _suppress_debug_logs():
    """Set libephemeris log level to WARNING to reduce noise in test output."""
    logger = logging.getLogger("libephemeris")
    original_level = logger.level
    logger.setLevel(logging.WARNING)
    yield
    logger.setLevel(original_level)


# ============================================================================
# GLOBAL SETUP - Disable strict precision for test suite
# ============================================================================
# By default, strict precision mode is enabled which requires SPK kernels for
# major asteroids (Chiron, Ceres, etc.). Most tests are designed to test the
# Keplerian fallback behavior, so we disable strict mode globally.
# Tests that specifically test strict mode should enable it explicitly.


@pytest.fixture(autouse=True)
def disable_strict_precision_for_tests():
    """Disable strict precision mode for all tests by default."""
    # Save current state
    original = ephem.get_strict_precision()
    # Disable strict mode for tests
    ephem.set_strict_precision(False)
    yield
    # Restore original state
    ephem.set_strict_precision(original)


# ============================================================================
# GLOBAL SETUP - Disable auto SPK download for test suite
# ============================================================================
# By default, auto SPK download is enabled. We disable it in tests to avoid
# network calls during test runs. Tests that specifically test auto-download
# behavior should enable it explicitly.


@pytest.fixture(autouse=True)
def disable_auto_spk_download_for_tests():
    """Disable auto SPK download for all tests by default."""
    original = ephem.get_auto_spk_download()
    ephem.set_auto_spk_download(False)
    yield
    ephem.set_auto_spk_download(original)


# NOTE on calc-mode isolation: since issue #30, close() preserves an
# explicitly set calculation mode (process-level configuration). The
# reset_ephemeris_state fixture below restores _CALC_MODE per test; test
# modules that must run without LEB (regardless of earlier explicit
# set_calc_mode calls in the worker) should monkeypatch
# state._CALC_MODE to None themselves rather than relying on close().


# ============================================================================
# TEST DATA FIXTURES
# ============================================================================


@pytest.fixture
def standard_jd():
    """Standard Julian Day for testing (J2000.0)."""
    return 2451545.0  # 2000-01-01 12:00:00 TT


@pytest.fixture
def j2000_jd():
    """J2000.0 epoch Julian Day."""
    return 2451545.0


@pytest.fixture
def j1900_jd():
    """J1900.0 epoch Julian Day."""
    return 2415020.0


@pytest.fixture
def unix_epoch_jd():
    """Unix epoch (1970-01-01 00:00) Julian Day."""
    return 2440587.5


@pytest.fixture
def test_dates():
    """Collection of test dates spanning different eras."""
    return [
        (2000, 1, 1, 12.0, "J2000"),
        (1980, 5, 20, 0.0, "Past"),
        (2024, 11, 5, 18.0, "Recent"),
        (1950, 10, 15, 6.0, "Mid-century"),
        (1550, 1, 1, 0.0, "DE440 Start"),
        (2650, 1, 1, 0.0, "DE440 End"),
        (2000, 2, 29, 12.0, "Leap year"),
        (1999, 12, 31, 23.99, "Y2K eve"),
        (2000, 1, 1, 0.01, "Y2K"),
    ]


@pytest.fixture
def test_locations():
    """Collection of test locations with various latitudes."""
    return [
        ("Rome", 41.9028, 12.4964, 0),
        ("London", 51.5074, -0.1278, 0),
        ("New York", 40.7128, -74.0060, 0),
        ("Sydney", -33.8688, 151.2093, 0),
        ("Tromso", 69.6492, 18.9553, 0),  # Arctic
        ("McMurdo", -77.8419, 166.6863, 0),  # Antarctic
        ("Equator", 0.0, 0.0, 0),  # Equator
        ("Tokyo", 35.6762, 139.6503, 0),
        ("Cape Town", -33.9249, 18.4241, 0),
        ("Reykjavik", 64.1466, -21.9426, 0),  # Near Arctic Circle
    ]


@pytest.fixture
def polar_locations():
    """Locations at polar latitudes for edge case testing."""
    return [
        ("Arctic 65N", 65.0, 0.0, 0),
        ("Arctic 67N", 67.0, 0.0, 0),  # Above Arctic Circle
        ("Arctic 70N", 70.0, 0.0, 0),
        ("Arctic 80N", 80.0, 0.0, 0),
        ("Arctic 85N", 85.0, 0.0, 0),
        ("North Pole", 89.9, 0.0, 0),
        ("Antarctic 65S", -65.0, 0.0, 0),
        ("Antarctic 70S", -70.0, 0.0, 0),
        ("Antarctic 85S", -85.0, 0.0, 0),
        ("South Pole", -89.9, 0.0, 0),
    ]


@pytest.fixture
def all_planets():
    """All major planets for testing."""
    return [
        (SUN, "Sun"),
        (MOON, "Moon"),
        (MERCURY, "Mercury"),
        (VENUS, "Venus"),
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
        (URANUS, "Uranus"),
        (NEPTUNE, "Neptune"),
        (PLUTO, "Pluto"),
    ]


@pytest.fixture
def inner_planets():
    """Inner planets (faster moving)."""
    return [
        (SUN, "Sun"),
        (MOON, "Moon"),
        (MERCURY, "Mercury"),
        (VENUS, "Venus"),
        (MARS, "Mars"),
    ]


@pytest.fixture
def outer_planets():
    """Outer planets (slower moving)."""
    return [
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
        (URANUS, "Uranus"),
        (NEPTUNE, "Neptune"),
        (PLUTO, "Pluto"),
    ]


@pytest.fixture
def lunar_points():
    """Lunar points (nodes and Lilith)."""
    return [
        (MEAN_NODE, "Mean Node"),
        (TRUE_NODE, "True Node"),
        (MEAN_APOG, "Mean Lilith"),
        (OSCU_APOG, "True Lilith"),
    ]


@pytest.fixture
def asteroids():
    """Main belt asteroids."""
    return [
        (CHIRON, "Chiron"),
        (PHOLUS, "Pholus"),
        (CERES, "Ceres"),
        (PALLAS, "Pallas"),
        (JUNO, "Juno"),
        (VESTA, "Vesta"),
    ]


@pytest.fixture
def all_house_systems():
    """All 19 supported house systems."""
    return [
        (ord("P"), "Placidus"),
        (ord("K"), "Koch"),
        (ord("O"), "Porphyry"),
        (ord("R"), "Regiomontanus"),
        (ord("C"), "Campanus"),
        (ord("E"), "Equal (Asc)"),
        (ord("A"), "Equal (Asc) alt"),
        (ord("W"), "Whole Sign"),
        (ord("M"), "Morinus"),
        (ord("B"), "Alcabitius"),
        (ord("T"), "Polich/Page (Topocentric)"),
        (ord("X"), "Meridian (Axial)"),
        (ord("V"), "Vehlow Equal"),
        (ord("H"), "Horizontal (Azimuthal)"),
        (ord("G"), "Gauquelin"),
        (ord("U"), "Krusinski-Pisa"),
        (ord("F"), "Carter Poli-Equatorial"),
        (ord("Y"), "APC Houses"),
        (ord("N"), "Natural Graduation (Polich/Page var)"),
    ]


@pytest.fixture
def common_house_systems():
    """Most commonly used house systems."""
    return [
        (ord("P"), "Placidus"),
        (ord("K"), "Koch"),
        (ord("O"), "Porphyry"),
        (ord("R"), "Regiomontanus"),
        (ord("C"), "Campanus"),
        (ord("E"), "Equal"),
        (ord("W"), "Whole Sign"),
    ]


@pytest.fixture
def all_ayanamshas():
    """All 47 predefined ayanamsha mode IDs."""
    return [
        (SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
        (SIDM_LAHIRI, "Lahiri"),
        (SIDM_DELUCE, "De Luce"),
        (SIDM_RAMAN, "Raman"),
        (SIDM_USHASHASHI, "Ushashashi"),
        (SIDM_KRISHNAMURTI, "Krishnamurti"),
        (SIDM_DJWHAL_KHUL, "Djwhal Khul"),
        (SIDM_YUKTESHWAR, "Yukteshwar"),
        (SIDM_JN_BHASIN, "JN Bhasin"),
        (SIDM_BABYL_KUGLER1, "Babylonian Kugler 1"),
        (SIDM_BABYL_KUGLER2, "Babylonian Kugler 2"),
        (SIDM_BABYL_KUGLER3, "Babylonian Kugler 3"),
        (SIDM_BABYL_HUBER, "Babylonian Huber"),
        (SIDM_BABYL_ETPSC, "Babylonian ETPSC"),
        (SIDM_ALDEBARAN_15TAU, "Aldebaran 15 Tau"),
        (SIDM_HIPPARCHOS, "Hipparchos"),
        (SIDM_SASSANIAN, "Sassanian"),
        (SIDM_GALCENT_0SAG, "Galactic Center 0 Sag"),
        (SIDM_J2000, "J2000"),
        (SIDM_J1900, "J1900"),
        (SIDM_B1950, "B1950"),
        (SIDM_SURYASIDDHANTA, "Surya Siddhanta"),
        (SIDM_SURYASIDDHANTA_MSUN, "Surya Siddhanta Mean Sun"),
        (SIDM_ARYABHATA, "Aryabhata"),
        (SIDM_ARYABHATA_MSUN, "Aryabhata Mean Sun"),
        (SIDM_SS_REVATI, "SS Revati"),
        (SIDM_SS_CITRA, "SS Citra"),
        (SIDM_TRUE_CITRA, "True Citra"),
        (SIDM_TRUE_REVATI, "True Revati"),
        (SIDM_TRUE_PUSHYA, "True Pushya"),
        (SIDM_GALCENT_RGILBRAND, "Galactic Center Rgilbrand"),
        (SIDM_GALEQU_IAU1958, "Galactic Equator IAU 1958"),
        (SIDM_GALEQU_TRUE, "Galactic Equator True"),
        (SIDM_GALEQU_MULA, "Galactic Equator Mula"),
        (SIDM_GALALIGN_MARDYKS, "Galactic Alignment Mardyks"),
        (SIDM_TRUE_MULA, "True Mula"),
        (SIDM_GALCENT_MULA_WILHELM, "Galactic Center Mula Wilhelm"),
        (SIDM_ARYABHATA_522, "Aryabhata 522"),
        (SIDM_BABYL_BRITTON, "Babylonian Britton"),
        (SIDM_TRUE_SHEORAN, "True Sheoran"),
        (SIDM_GALCENT_COCHRANE, "Galactic Center Cochrane"),
        (SIDM_GALEQU_FIORENZA, "Galactic Equator Fiorenza"),
        (SIDM_VALENS_MOON, "Valens Moon"),
        (SIDM_LAHIRI_1940, "Lahiri 1940"),
        (SIDM_LAHIRI_VP285, "Lahiri VP285"),
        (SIDM_KRISHNAMURTI_VP291, "Krishnamurti VP291"),
        (SIDM_LAHIRI_ICRC, "Lahiri ICRC"),
    ]


@pytest.fixture
def major_ayanamshas():
    """Major ayanamsha modes for testing."""
    return [
        (SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
        (SIDM_LAHIRI, "Lahiri"),
        (SIDM_RAMAN, "Raman"),
        (SIDM_TRUE_CITRA, "True Citra"),
        (SIDM_GALCENT_0SAG, "Galactic Center 0 Sag"),
    ]


@pytest.fixture
def star_based_ayanamshas():
    """Native catalog-direction ayanamshas (require relaxed tolerance)."""
    return [
        SIDM_TRUE_CITRA,
        SIDM_GALCENT_0SAG,
    ]


# ============================================================================
# TOLERANCE FIXTURES
# ============================================================================


@pytest.fixture
def default_tolerances():
    """Default tolerance values for comparisons."""
    return {
        "longitude": 0.001,  # degrees
        "latitude": 0.001,  # degrees
        "distance": 0.0001,  # AU
        "velocity": 0.01,  # degrees/day
        "ayanamsha": 0.06,  # degrees (relaxed for star-based)
        "house_cusp": 0.001,  # degrees
        "julian_day": 1e-10,  # days
        "time_seconds": 0.5,  # seconds (for Delta T)
        "crossing_time": 120,  # seconds
    }


@pytest.fixture
def tight_tolerances():
    """Tight tolerances for high-precision tests."""
    return {
        "longitude": 1 / 3600,  # 1 arcsecond
        "latitude": 1 / 3600,  # 1 arcsecond
        "distance": 1e-6,  # AU
        "velocity": 0.001,  # degrees/day
        "ayanamsha": 0.001,  # degrees
        "house_cusp": 1 / 3600,  # 1 arcsecond
        "julian_day": 1e-12,  # days
        "crossing_time": 1,  # second
    }


@pytest.fixture
def relaxed_tolerances():
    """Relaxed tolerances for approximate calculations."""
    return {
        "longitude": 1.0,  # degrees
        "latitude": 1.0,  # degrees
        "distance": 0.01,  # AU
        "velocity": 0.1,  # degrees/day
        "ayanamsha": 1.0,  # degrees
        "house_cusp": 1.0,  # degrees
        "asteroid": 5.0,  # degrees (Keplerian approx)
        "tno": 10.0,  # degrees (Keplerian approx)
        "lilith": 5.0,  # degrees
    }


# ============================================================================
# RANDOM DATA GENERATORS
# ============================================================================


@pytest.fixture
def random_dates_in_de440_range():
    """Generate random dates within DE440 valid range (1550-2650)."""

    def _generate(n=100, seed=42):
        random.seed(seed)
        dates = []
        for _ in range(n):
            year = random.randint(1550, 2650)
            month = random.randint(1, 12)
            day = random.randint(1, 28)  # Safe for all months
            hour = random.uniform(0, 24)
            jd = ephem.julday(year, month, day, hour)
            dates.append((year, month, day, hour, jd))
        return dates

    return _generate


@pytest.fixture
def random_dates_in_de421_range():
    """Generate random dates within DE421 valid range (1900-2050) - legacy."""

    def _generate(n=100, seed=42):
        random.seed(seed)
        dates = []
        for _ in range(n):
            year = random.randint(1900, 2050)
            month = random.randint(1, 12)
            day = random.randint(1, 28)  # Safe for all months
            hour = random.uniform(0, 24)
            jd = ephem.julday(year, month, day, hour)
            dates.append((year, month, day, hour, jd))
        return dates

    return _generate


@pytest.fixture
def random_locations():
    """Generate random geographic locations."""

    def _generate(n=50, seed=42):
        random.seed(seed)
        locations = []
        for i in range(n):
            lat = random.uniform(-85, 85)
            lon = random.uniform(-180, 180)
            alt = random.uniform(0, 1000)
            locations.append((f"Random_{i}", lat, lon, alt))
        return locations

    return _generate


@pytest.fixture
def random_longitudes():
    """Generate random celestial longitudes."""

    def _generate(n=100, seed=42):
        random.seed(seed)
        return [random.uniform(0, 360) for _ in range(n)]

    return _generate


# ============================================================================
# FAMOUS PEOPLE DATA (for integration tests)
# ============================================================================


@pytest.fixture
def famous_people():
    """Birth data of famous people for testing."""
    return [
        {
            "name": "Albert Einstein",
            "year": 1879,
            "month": 3,
            "day": 14,
            "hour": 11.5,
            "lat": 48.4010,
            "lon": 9.9876,  # Ulm, Germany
        },
        {
            "name": "Marilyn Monroe",
            "year": 1926,
            "month": 6,
            "day": 1,
            "hour": 9.5,
            "lat": 34.0522,
            "lon": -118.2437,  # Los Angeles
        },
        {
            "name": "Mahatma Gandhi",
            "year": 1869,
            "month": 10,
            "day": 2,
            "hour": 7.2,
            "lat": 21.6417,
            "lon": 69.6293,  # Porbandar, India
        },
        {
            "name": "Pablo Picasso",
            "year": 1881,
            "month": 10,
            "day": 25,
            "hour": 23.25,
            "lat": 36.7213,
            "lon": -4.4214,  # Malaga, Spain
        },
        {
            "name": "Nelson Mandela",
            "year": 1918,
            "month": 7,
            "day": 18,
            "hour": 14.0,
            "lat": -31.5875,
            "lon": 28.7833,  # Mvezo, South Africa
        },
    ]


# ============================================================================
# KNOWN ASTRONOMICAL EVENTS
# ============================================================================


@pytest.fixture
def known_equinoxes():
    """Known vernal equinox dates (Sun at 0 Aries)."""
    return [
        (2020, 3, 20, 3.83),  # 2020 vernal equinox ~03:50 UT
        (2021, 3, 20, 9.63),  # 2021 vernal equinox ~09:37 UT
        (2022, 3, 20, 15.53),  # 2022 vernal equinox ~15:33 UT
        (2023, 3, 20, 21.42),  # 2023 vernal equinox ~21:25 UT
        (2024, 3, 20, 3.07),  # 2024 vernal equinox ~03:06 UT
    ]


@pytest.fixture
def known_solstices():
    """Known summer solstice dates (Sun at 90 degrees)."""
    return [
        (2020, 6, 20, 21.73),  # 2020 summer solstice
        (2021, 6, 21, 3.47),  # 2021 summer solstice
        (2022, 6, 21, 9.25),  # 2022 summer solstice
        (2023, 6, 21, 14.95),  # 2023 summer solstice
        (2024, 6, 20, 20.82),  # 2024 summer solstice
    ]


# ============================================================================
# SETUP/TEARDOWN
# ============================================================================


@pytest.fixture(autouse=True)
def reset_ephemeris_state():
    """Reset ephemeris state and clear caches before each test.

    Saves and restores ALL mutable global state that could cause cross-test
    contamination in parallel runs (pytest-xdist). Without this, tests that
    call close(), set_calc_mode(), or pop env vars will corrupt state for
    subsequent tests in the same worker process.
    """
    import os

    from libephemeris import state

    # Save critical global state before the test
    saved_calc_mode = state._CALC_MODE
    saved_leb_file = state._LEB_FILE
    saved_leb_reader = state._LEB_READER
    saved_precision_tier = state._PRECISION_TIER
    saved_ephemeris_file = state._EPHEMERIS_FILE
    saved_ephemeris_file_explicit = state._EPHEMERIS_FILE_EXPLICIT
    saved_topo = state._TOPO

    # Save env vars that some tests pop/modify
    saved_env_mode = os.environ.get("LIBEPHEMERIS_MODE")
    saved_env_leb = os.environ.get("LIBEPHEMERIS_LEB")

    ephem.set_sid_mode(SIDM_FAGAN_BRADLEY)

    yield

    state._CALC_MODE = saved_calc_mode
    state._PRECISION_TIER = saved_precision_tier
    state._EPHEMERIS_FILE = saved_ephemeris_file
    state._EPHEMERIS_FILE_EXPLICIT = saved_ephemeris_file_explicit
    state._TOPO = saved_topo

    state._LEB_FILE = saved_leb_file
    state._LEB_READER = None

    # The fictitious-body trust caches are keyed by a (size, mtime) stat
    # fingerprint that does not change when a monkeypatched DATA_FILES pin is
    # undone at teardown, so a verdict computed under a temporary pin could
    # otherwise leak into the next test in this worker.
    state._FICTITIOUS_SOURCE_CACHE.clear()
    state._FICTITIOUS_TRUST_BY_READER.clear()

    if saved_env_mode is not None:
        os.environ["LIBEPHEMERIS_MODE"] = saved_env_mode
    else:
        os.environ.pop("LIBEPHEMERIS_MODE", None)

    if saved_env_leb is not None:
        os.environ["LIBEPHEMERIS_LEB"] = saved_env_leb
    else:
        os.environ.pop("LIBEPHEMERIS_LEB", None)

    ephem.set_sid_mode(SIDM_FAGAN_BRADLEY)


# ============================================================================
# ASSERTION HELPERS
# ============================================================================


@pytest.fixture
def assert_position_close():
    """Assert two positions are close within tolerance."""

    def _assert(pos1, pos2, tolerance=0.001, msg=""):
        lon_diff = abs(pos1[0] - pos2[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        lat_diff = abs(pos1[1] - pos2[1])
        dist_diff = abs(pos1[2] - pos2[2])

        assert lon_diff < tolerance, f"{msg} Longitude diff {lon_diff} >= {tolerance}"
        assert lat_diff < tolerance, f"{msg} Latitude diff {lat_diff} >= {tolerance}"
        assert dist_diff < tolerance * 10, (
            f"{msg} Distance diff {dist_diff} >= {tolerance * 10}"
        )

    return _assert


@pytest.fixture
def assert_angle_close():
    """Assert two angles are close (handles wrapping)."""

    def _assert(angle1, angle2, tolerance=0.001, msg=""):
        diff = abs(angle1 - angle2)
        if diff > 180:
            diff = 360 - diff
        assert diff < tolerance, f"{msg} Angle diff {diff} >= {tolerance}"

    return _assert


@pytest.fixture
def normalize_angle():
    """Normalize angle to 0-360 range."""

    def _normalize(angle):
        while angle < 0:
            angle += 360
        while angle >= 360:
            angle -= 360
        return angle

    return _normalize


# ============================================================================
# PROGRESS REPORTING
# ============================================================================


class ProgressReporter:
    """
    Progress reporter for long-running test loops.

    Usage:
        def test_long_loop(progress_reporter):
            items = range(100)
            progress = progress_reporter("Testing items", len(items))
            for i, item in enumerate(items):
                # ... do work ...
                progress.update(i, f"item={item}")
            progress.done()

    The reporter will print progress at regular intervals (default 10%),
    showing current iteration, percentage, and optional context.
    """

    def __init__(self, description="Progress", total=100, report_every=10):
        """
        Initialize progress reporter.

        Args:
            description: Description of the operation being tracked
            total: Total number of iterations
            report_every: Report progress every N percent (default 10%)
        """
        self.description = description
        self.total = total
        self.report_every = report_every
        self._last_reported_pct = -1

    def update(self, current, context=""):
        """
        Update progress and print if at a reporting interval.

        Args:
            current: Current iteration (0-based)
            context: Optional context string to display
        """
        if self.total <= 0:
            return

        pct = int((current + 1) * 100 / self.total)

        # Report at every report_every percent
        report_pct = (pct // self.report_every) * self.report_every

        if report_pct > self._last_reported_pct:
            self._last_reported_pct = report_pct
            ctx_str = f" [{context}]" if context else ""
            print(f"  {self.description}: {current + 1}/{self.total} ({pct}%){ctx_str}")

    def done(self, summary=""):
        """Print completion message."""
        if summary:
            print(f"  {self.description}: DONE - {summary}")
        else:
            print(f"  {self.description}: DONE ({self.total} iterations)")


@pytest.fixture
def progress_reporter():
    """
    Fixture to report progress in long-running test loops.

    Usage:
        def test_many_iterations(progress_reporter):
            data = list(range(100))
            pr = progress_reporter("Testing values", len(data))
            for i, val in enumerate(data):
                # ... do test work ...
                pr.update(i, f"val={val}")
            pr.done()
    """

    def _create_reporter(description="Progress", total=100, report_every=10):
        return ProgressReporter(description, total, report_every)

    return _create_reporter


# ============================================================================
# MARKERS
# ============================================================================


def pytest_configure(config):
    """Configure custom pytest markers."""
    config.addinivalue_line("markers", "unit: mark test as a unit test")
    config.addinivalue_line("markers", "integration: mark test as an integration test")
    config.addinivalue_line("markers", "slow: mark test as slow running")
    config.addinivalue_line(
        "markers", "network: marks tests that require network access"
    )
    config.addinivalue_line("markers", "precision: mark test as high-precision")
    config.addinivalue_line("markers", "edge_case: mark test as edge case")
