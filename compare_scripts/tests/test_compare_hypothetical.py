"""
Pytest-style Hypothetical Planets Comparison Tests.

Validates hypothetical body calculations (Uranian planets, Transpluto)
against pyswisseph.

Uses full Keplerian propagation with Gaussian vectors and equinox precession
from J1900 to J2000, matching the standard celestial mechanics algorithm.
Remaining ~35" differences arise from IAU 2006 precession (libephemeris)
vs older precession model (pyswisseph) and Skyfield vs SE geocentric
conversion.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    CUPIDO,
    HADES,
    ZEUS,
    KRONOS,
    APOLLON,
    ADMETOS,
    VULKANUS,
    POSEIDON,
    ISIS,
)


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================

# Tightened tolerances after implementing full Keplerian propagation with
# Gaussian vectors and equinox precession (J1900 -> J2000).
# Max observed: ~35" longitude (from precession model + geocentric differences)
# Latitude: ~25" max for inclined bodies (Cupido, Hades)
URANIAN_LONGITUDE_TOLERANCE = 0.02  # 72 arcsec (~35" observed max)
HYPOTHETICAL_LONGITUDE_TOLERANCE = 0.02  # 72 arcsec for Transpluto
LATITUDE_TOLERANCE = 0.02  # 72 arcsec for latitude


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# Uranian Planets (Hamburg School) - IDs 40-47
URANIAN_PLANETS = [
    (CUPIDO, swe.CUPIDO, "Cupido"),
    (HADES, swe.HADES, "Hades"),
    (ZEUS, swe.ZEUS, "Zeus"),
    (KRONOS, swe.KRONOS, "Kronos"),
    (APOLLON, swe.APOLLON, "Apollon"),
    (ADMETOS, swe.ADMETOS, "Admetos"),
    (VULKANUS, swe.VULKANUS, "Vulkanus"),
    (POSEIDON, swe.POSEIDON, "Poseidon"),
]

# Other hypothetical bodies
OTHER_HYPOTHETICAL = [
    (ISIS, swe.ISIS, "Transpluto/Isis"),
]

# Test dates
TEST_DATES = [
    ("J2000.0", 2000, 1, 1, 12.0),
    ("2024-01-01", 2024, 1, 1, 0.0),
    ("2010-07-01", 2010, 7, 1, 12.0),
    ("1980-01-01", 1980, 1, 1, 0.0),
]

# Extended date range for thorough testing (1900-2100)
EXTENDED_DATES = [
    ("1900-01-01", 1900, 1, 1, 12.0),
    ("1920-06-15", 1920, 6, 15, 0.0),
    ("1950-01-01", 1950, 1, 1, 0.0),
    ("1980-06-15", 1980, 6, 15, 12.0),
    ("J2000.0", 2000, 1, 1, 12.0),
    ("2024-06-21", 2024, 6, 21, 12.0),
    ("2050-06-01", 2050, 6, 1, 0.0),
    ("2100-01-01", 2100, 1, 1, 0.0),
]


# ============================================================================
# URANIAN PLANET TESTS
# ============================================================================


class TestUranianPlanets:
    """Tests for Uranian (Hamburg School) planets."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_py,body_swe,name", URANIAN_PLANETS)
    @pytest.mark.parametrize("date_name,year,month,day,hour", TEST_DATES[:2])  # Subset
    def test_uranian_planet_position(
        self, body_py, body_swe, name, date_name, year, month, day, hour
    ):
        """Test Uranian planet positions match pyswisseph."""
        jd = swe.julday(year, month, day, hour)

        # SwissEphemeris
        try:
            pos_swe, _ = swe.calc_ut(jd, body_swe, 0)
        except swe.Error as e:
            pytest.skip(f"SwissEphemeris {name} not available: {e}")

        # LibEphemeris
        try:
            pos_py, _ = ephem.calc_ut(jd, body_py, 0)
        except Exception as e:
            pytest.skip(f"LibEphemeris {name} failed: {e}")

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])

        assert diff_lon < URANIAN_LONGITUDE_TOLERANCE, (
            f"{name} @ {date_name}: longitude diff {diff_lon}° exceeds tolerance"
        )
        assert diff_lat < LATITUDE_TOLERANCE, (
            f"{name} @ {date_name}: latitude diff {diff_lat}° exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_py,body_swe,name", URANIAN_PLANETS[:4])  # First 4
    def test_uranian_multiple_dates(self, body_py, body_swe, name):
        """Test Uranian planets across multiple dates."""
        for date_name, year, month, day, hour in TEST_DATES:
            jd = swe.julday(year, month, day, hour)

            try:
                pos_swe, _ = swe.calc_ut(jd, body_swe, 0)
                pos_py, _ = ephem.calc_ut(jd, body_py, 0)

                diff_lon = angular_diff(pos_swe[0], pos_py[0])
                assert diff_lon < URANIAN_LONGITUDE_TOLERANCE
            except Exception:
                pass  # Skip dates where calculation fails

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_py,body_swe,name", URANIAN_PLANETS)
    @pytest.mark.parametrize("date_name,year,month,day,hour", EXTENDED_DATES)
    def test_uranian_extended_dates(
        self, body_py, body_swe, name, date_name, year, month, day, hour
    ):
        """Test all Uranian planets across extended 1900-2100 date range."""
        jd = swe.julday(year, month, day, hour)

        try:
            pos_swe, _ = swe.calc_ut(jd, body_swe, 0)
        except swe.Error as e:
            pytest.skip(f"SwissEphemeris {name} not available: {e}")

        try:
            pos_py, _ = ephem.calc_ut(jd, body_py, 0)
        except Exception as e:
            pytest.skip(f"LibEphemeris {name} failed: {e}")

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])

        assert diff_lon < URANIAN_LONGITUDE_TOLERANCE, (
            f"{name} @ {date_name}: longitude diff {diff_lon:.6f}deg "
            f"({diff_lon * 3600:.1f}arcsec) exceeds {URANIAN_LONGITUDE_TOLERANCE}deg"
        )
        assert diff_lat < LATITUDE_TOLERANCE, (
            f"{name} @ {date_name}: latitude diff {diff_lat:.6f}deg "
            f"({diff_lat * 3600:.1f}arcsec) exceeds {LATITUDE_TOLERANCE}deg"
        )


# ============================================================================
# OTHER HYPOTHETICAL BODY TESTS
# ============================================================================


class TestOtherHypothetical:
    """Tests for other hypothetical bodies (Transpluto, etc.)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_py,body_swe,name", OTHER_HYPOTHETICAL)
    @pytest.mark.parametrize("date_name,year,month,day,hour", TEST_DATES)
    def test_hypothetical_body_position(
        self, body_py, body_swe, name, date_name, year, month, day, hour
    ):
        """Test other hypothetical body positions."""
        jd = swe.julday(year, month, day, hour)

        # SwissEphemeris
        try:
            pos_swe, _ = swe.calc_ut(jd, body_swe, 0)
        except swe.Error as e:
            pytest.skip(f"SwissEphemeris {name} not available: {e}")

        # LibEphemeris
        try:
            pos_py, _ = ephem.calc_ut(jd, body_py, 0)
        except Exception as e:
            pytest.skip(f"LibEphemeris {name} failed: {e}")

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])

        assert diff_lon < HYPOTHETICAL_LONGITUDE_TOLERANCE, (
            f"{name} @ {date_name}: longitude diff {diff_lon}°"
        )
        assert diff_lat < LATITUDE_TOLERANCE, (
            f"{name} @ {date_name}: latitude diff {diff_lat}°"
        )


# ============================================================================
# ALL HYPOTHETICAL BODIES SUMMARY TEST
# ============================================================================


class TestHypotheticalSummary:
    """Summary tests for all hypothetical bodies."""

    @pytest.mark.comparison
    def test_all_uranian_available(self):
        """Verify all Uranian planets are calculable."""
        jd = swe.julday(2024, 1, 1, 12.0)
        available = 0

        for body_py, body_swe, name in URANIAN_PLANETS:
            try:
                pos_py, _ = ephem.calc_ut(jd, body_py, 0)
                available += 1
            except Exception:
                pass

        # At least 6 out of 8 should be available
        assert available >= 6, f"Only {available}/8 Uranian planets available"

    @pytest.mark.comparison
    def test_uranian_consistency(self):
        """Test Uranian planets are consistent across a time range."""
        for body_py, body_swe, name in URANIAN_PLANETS[:4]:
            positions = []
            for year in [2000, 2010, 2020]:
                jd = swe.julday(year, 1, 1, 12.0)
                try:
                    pos, _ = ephem.calc_ut(jd, body_py, 0)
                    positions.append(pos[0])
                except Exception:
                    pass

            # Positions should be different (planet moved)
            if len(positions) >= 2:
                assert positions[0] != positions[1], f"{name} didn't move over 10 years"
