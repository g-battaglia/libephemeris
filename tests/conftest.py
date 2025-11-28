"""
pytest configuration and shared fixtures for LibEphemeris tests.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


# ============================================================================
# TEST DATA FIXTURES
# ============================================================================


@pytest.fixture
def standard_jd():
    """Standard Julian Day for testing (J2000.0)."""
    return 2451545.0  # 2000-01-01 12:00:00 TT


@pytest.fixture
def test_dates():
    """Collection of test dates spanning different eras."""
    return [
        (2000, 1, 1, 12.0, "J2000"),
        (1980, 5, 20, 0.0, "Past"),
        (2024, 11, 5, 18.0, "Recent"),
        (1950, 10, 15, 6.0, "Mid-century"),
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
    ]


@pytest.fixture
def all_planets():
    """All major planets for testing."""
    return [
        (SE_SUN, "Sun"),
        (SE_MOON, "Moon"),
        (SE_MERCURY, "Mercury"),
        (SE_VENUS, "Venus"),
        (SE_MARS, "Mars"),
        (SE_JUPITER, "Jupiter"),
        (SE_SATURN, "Saturn"),
        (SE_URANUS, "Uranus"),
        (SE_NEPTUNE, "Neptune"),
        (SE_PLUTO, "Pluto"),
    ]


@pytest.fixture
def all_house_systems():
    """All supported house systems."""
    return [
        (b"P", "Placidus"),
        (b"K", "Koch"),
        (b"R", "Regiomontanus"),
        (b"C", "Campanus"),
        (b"E", "Equal"),
        (b"W", "Whole Sign"),
        (b"O", "Porphyry"),
        (b"B", "Alcabitius"),
        (b"T", "Polich/Page"),
        (b"M", "Morinus"),
        (b"X", "Meridian"),
        (b"V", "Vehlow"),
        (b"H", "Horizontal"),
        (b"F", "Carter"),
        (b"U", "Krusinski"),
        (b"N", "Natural Gradient"),
        (b"G", "Gauquelin"),
    ]


@pytest.fixture
def major_ayanamshas():
    """Major ayanamsha modes for testing."""
    return [
        (SE_SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
        (SE_SIDM_LAHIRI, "Lahiri"),
        (SE_SIDM_RAMAN, "Raman"),
        (SE_SIDM_TRUE_CITRA, "True Citra"),
        (SE_SIDM_GALCENT_0SAG, "Galactic Center 0 Sag"),
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
        "velocity": 0.01,  # degrees/day or AU/day
        "ayanamsha": 0.06,  # degrees (relaxed for star-based)
        "house_cusp": 0.001,  # degrees
    }


# ============================================================================
# COMPARISON FIXTURES
# ============================================================================


@pytest.fixture
def compare_with_swisseph():
    """Helper function to compare libephemeris results with SwissEph."""

    def _compare(jd, planet_id, flags, tolerance=0.001):
        """
        Compare planetary positions between implementations.

        Returns:
            (bool, dict): (passed, differences)
        """
        # Calculate with SwissEph
        res_swe, _ = swe.calc_ut(jd, planet_id, flags)

        # Calculate with libephemeris
        res_py, _ = ephem.swe_calc_ut(jd, planet_id, flags)

        # Compare
        diffs = {
            "lon": abs(res_swe[0] - res_py[0]),
            "lat": abs(res_swe[1] - res_py[1]),
            "dist": abs(res_swe[2] - res_py[2]),
        }

        # Handle wrap-around for longitude
        if diffs["lon"] > 180:
            diffs["lon"] = 360 - diffs["lon"]

        passed = all(d < tolerance for d in diffs.values())

        return passed, diffs

    return _compare


# ============================================================================
# SETUP/TEARDOWN
# ============================================================================


@pytest.fixture(autouse=True)
def reset_ephemeris_state():
    """Reset ephemeris state before each test."""
    # Reset to default sidereal mode
    ephem.swe_set_sid_mode(SE_SIDM_FAGAN_BRADLEY)
    swe.set_sid_mode(swe.SIDM_FAGAN_BRADLEY)

    yield

    # Cleanup after test (if needed)
    pass


# ============================================================================
# MARKERS
# ============================================================================


def pytest_configure(config):
    """Configure custom pytest markers."""
    config.addinivalue_line("markers", "unit: mark test as a unit test")
    config.addinivalue_line("markers", "integration: mark test as an integration test")
    config.addinivalue_line("markers", "slow: mark test as slow running")
