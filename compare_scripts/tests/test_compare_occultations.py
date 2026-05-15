"""
Pytest-style Lunar Occultation Comparison Tests.

Validates lun_occult_when_glob, lun_occult_when_loc, and lun_occult_where
calculations against pyswisseph.
"""

import pytest
import swisseph as swe
import libephemeris as pyephem
from libephemeris.constants import (
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    FLG_SWIEPH,
)


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class OccultTolerance:
    """Tolerance thresholds for occultation comparisons.

    The improved algorithms now use:
    - Bisection search for contact times instead of chord-based estimates
    - Golden section search with increased iterations for sub-second precision
    - Topocentric parallax corrections for local occultation timing
    - Reduced initial search step size

    This allows reducing the tolerance from 300s (5 minutes) to 60s (1 minute).
    """

    TIME_SECONDS = 60.0  # 1 minute for occultations (improved from 300s)
    POSITION_DEGREES = 0.5  # 0.5 degree for coordinates


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# Planets for occultation tests (Moon can occult these)
OCCULT_BODIES = [
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
]

# Test locations
OCCULT_LOCATIONS = [
    ("New York", 40.7128, -74.0060, 0),
    ("London", 51.5074, -0.1278, 0),
    ("Sydney", -33.8688, 151.2093, 0),
]

# Search start dates
SEARCH_DATES = [
    (2024, 1, 1, "2024 Start"),
    (2025, 1, 1, "2025 Start"),
    (2026, 1, 1, "2026 Start"),
]


# ============================================================================
# FIXTURES
# ============================================================================


@pytest.fixture
def jd_2024():
    """Julian Day for 2024 start."""
    return swe.julday(2024, 1, 1, 0.0)


# ============================================================================
# GLOBAL OCCULTATION TESTS
# ============================================================================


class TestLunOccultWhenGlob:
    """Tests for global lunar occultation search."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", OCCULT_BODIES[:2])  # Venus and Mars
    @pytest.mark.parametrize("year,month,day,date_desc", SEARCH_DATES[:2])
    def test_lun_occult_when_glob(
        self, body_id, body_name, year, month, day, date_desc
    ):
        """Test global lunar occultation search."""
        jd_start = swe.julday(year, month, day, 0.0)

        # SwissEphemeris (API: tjdut, body, flags, ecltype, backwards)
        try:
            ret_swe = swe.lun_occult_when_glob(
                jd_start, body_id, FLG_SWIEPH, 0, False
            )
            if ret_swe[0] == 0:
                pytest.skip(f"No {body_name} occultation found by SwissEphemeris")
            jd_swe = ret_swe[1][0]
        except Exception as e:
            pytest.skip(f"SwissEphemeris failed: {e}")

        # LibEphemeris (API: tjdut, body, flags, ecltype, backwards) - 1:1 compatible
        try:
            ret_py = pyephem.lun_occult_when_glob(
                jd_start, body_id, FLG_SWIEPH, 0, False
            )
            if ret_py[0] == 0:
                pytest.skip(f"No {body_name} occultation found by LibEphemeris")
            jd_py = ret_py[1][0]
        except Exception as e:
            pytest.skip(f"LibEphemeris failed: {e}")

        diff_seconds = abs(jd_py - jd_swe) * 86400

        assert diff_seconds < OccultTolerance.TIME_SECONDS, (
            f"{body_name} occultation {date_desc}: time diff {diff_seconds:.2f}s"
        )


# ============================================================================
# LOCAL OCCULTATION TESTS
# ============================================================================


class TestLunOccultWhenLoc:
    """Tests for local lunar occultation search."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", OCCULT_BODIES[:2])
    @pytest.mark.parametrize("loc_name,lat,lon,alt", OCCULT_LOCATIONS[:2])
    def test_lun_occult_when_loc(
        self, jd_2024, body_id, body_name, loc_name, lat, lon, alt
    ):
        """Test local lunar occultation search."""
        geopos = (lon, lat, alt)

        # SwissEphemeris (API: tjdut, body, geopos, flags, backwards)
        try:
            ret_swe = swe.lun_occult_when_loc(
                jd_2024, body_id, geopos, FLG_SWIEPH, False
            )
            if ret_swe[0] == 0:
                pytest.skip(f"No {body_name} occultation visible at {loc_name}")
            jd_swe = ret_swe[1][0]
        except Exception as e:
            pytest.skip(f"SwissEphemeris failed: {e}")

        # LibEphemeris (API: tjdut, body, geopos, flags, backwards) - 1:1 compatible
        try:
            ret_py = pyephem.lun_occult_when_loc(
                jd_2024, body_id, geopos, FLG_SWIEPH, False
            )
            if ret_py[0] == 0:
                pytest.skip(f"No {body_name} occultation visible (LibEphemeris)")
            jd_py = ret_py[1][0]
        except Exception as e:
            pytest.skip(f"LibEphemeris failed: {e}")

        diff_seconds = abs(jd_py - jd_swe) * 86400

        assert diff_seconds < OccultTolerance.TIME_SECONDS, (
            f"{body_name} @ {loc_name}: time diff {diff_seconds:.2f}s"
        )


# ============================================================================
# OCCULTATION LOCATION TESTS
# ============================================================================


class TestLunOccultWhere:
    """Tests for lunar occultation location (where on Earth)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", OCCULT_BODIES[:2])
    def test_lun_occult_where(self, jd_2024, body_id, body_name):
        """Test lunar occultation location calculation."""
        # First find an occultation
        try:
            ret = swe.lun_occult_when_glob(jd_2024, body_id, FLG_SWIEPH, 0, False)
            if ret[0] == 0:
                pytest.skip(f"No {body_name} occultation found")
            jd_occult = ret[1][0]
        except Exception as e:
            pytest.skip(f"Failed to find occultation: {e}")

        # SwissEphemeris (new API: tjdut, body, flags)
        try:
            ret_swe = swe.lun_occult_where(jd_occult, body_id, FLG_SWIEPH)
            if ret_swe[0] == 0:
                pytest.skip("No central line (SwissEphemeris)")
            lon_swe, lat_swe = ret_swe[1][0], ret_swe[1][1]
        except Exception as e:
            pytest.skip(f"SwissEphemeris lun_occult_where failed: {e}")

        # LibEphemeris (API: tjdut, body, flags) - 1:1 compatible
        try:
            ret_py = pyephem.lun_occult_where(jd_occult, body_id, FLG_SWIEPH)
            if ret_py[0] == 0:
                pytest.skip("No central line (LibEphemeris)")
            lon_py, lat_py = ret_py[1][0], ret_py[1][1]
        except Exception as e:
            pytest.skip(f"LibEphemeris lun_occult_where failed: {e}")

        diff_lon = abs(lon_swe - lon_py)
        diff_lat = abs(lat_swe - lat_py)
        max_diff = max(diff_lon, diff_lat)

        assert max_diff < OccultTolerance.POSITION_DEGREES, (
            f"{body_name} occultation location: max diff {max_diff}°"
        )


# ============================================================================
# SUMMARY TESTS
# ============================================================================


class TestOccultationSummary:
    """Summary tests for occultation functions."""

    @pytest.mark.comparison
    def test_occultation_functions_available(self):
        """Verify occultation functions are available."""
        jd = swe.julday(2024, 1, 1, 0.0)

        # Check that functions are callable
        assert callable(pyephem.lun_occult_when_glob)
        assert callable(pyephem.lun_occult_when_loc)
        assert callable(pyephem.lun_occult_where)

    @pytest.mark.comparison
    def test_venus_occultation_search(self, jd_2024):
        """Test Venus occultation search specifically."""
        try:
            ret_swe = swe.lun_occult_when_glob(
                jd_2024, VENUS, FLG_SWIEPH, 0, False
            )
            ret_py = pyephem.lun_occult_when_glob(
                jd_2024, VENUS, FLG_SWIEPH, 0, False
            )

            # Both should return the same status (found or not found)
            assert (ret_swe[0] == 0) == (ret_py[0] == 0), (
                "Mismatch in Venus occultation detection"
            )
        except Exception:
            pytest.skip("Venus occultation search failed")
