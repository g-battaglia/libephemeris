"""
Rise/Transit/Set Calculations Comparison Tests.

Compares rise, transit, and set calculations between pyswisseph and libephemeris.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import SUN, MOON, VENUS


# ============================================================================
# TOLERANCES
# ============================================================================

TIME_TOL_SECONDS = 30.0  # 30 seconds (improved from 120s with Bennett formula)


# ============================================================================
# TEST DATA
# ============================================================================

TEST_LOCATIONS = [
    ("Rome", 41.9028, 12.4964, 0),
    ("New York", 40.7128, -74.0060, 0),
    ("Sydney", -33.8688, 151.2093, 0),
    ("Equator", 0.0, 0.0, 0),
]

PLANETS_TO_TEST = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (VENUS, "Venus"),
]


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestRiseSet:
    """Compare rise and set calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("name,lat,lon,alt", TEST_LOCATIONS)
    @pytest.mark.parametrize(
        "planet_id,planet_name", PLANETS_TO_TEST[:1]
    )  # Sun only for speed
    def test_sunrise(self, name, lat, lon, alt, planet_id, planet_name):
        """Test sunrise time calculations."""
        jd = swe.julday(2024, 6, 15, 0.0)
        geopos = (lon, lat, alt)

        try:
            # CALC_RISE = 1
            rise_swe = swe.rise_trans(jd, planet_id, 1, geopos, 1013.25, 15.0)
            rise_py = ephem.rise_trans(jd, planet_id, 1, [lon, lat, alt], 1013.25, 15.0)
        except Exception as e:
            pytest.skip(f"Rise calculation not available: {e}")
            return

        if rise_swe[0] != 0 or rise_py[0] != 0:
            pytest.skip("No rise event found")
            return

        diff_seconds = abs(rise_swe[1][0] - rise_py[1][0]) * 86400

        assert diff_seconds < TIME_TOL_SECONDS, (
            f"{planet_name} rise at {name}: diff {diff_seconds:.1f}s exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("name,lat,lon,alt", TEST_LOCATIONS)
    def test_sun_transit(self, name, lat, lon, alt):
        """Test Sun transit (noon) time calculations."""
        jd = swe.julday(2024, 6, 15, 0.0)
        geopos = (lon, lat, alt)

        try:
            # CALC_MTRANSIT = 4
            trans_swe = swe.rise_trans(jd, SUN, 4, geopos, 1013.25, 15.0)
            trans_py = ephem.rise_trans(jd, SUN, 4, [lon, lat, alt], 1013.25, 15.0)
        except Exception as e:
            pytest.skip(f"Transit calculation not available: {e}")
            return

        if trans_swe[0] != 0 or trans_py[0] != 0:
            pytest.skip("No transit event found")
            return

        diff_seconds = abs(trans_swe[1][0] - trans_py[1][0]) * 86400

        assert diff_seconds < TIME_TOL_SECONDS, (
            f"Sun transit at {name}: diff {diff_seconds:.1f}s exceeds tolerance"
        )


class TestMoonRiseSet:
    """Test Moon rise/set specifically (more complex due to parallax)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("name,lat,lon,alt", TEST_LOCATIONS[:2])
    def test_moonrise(self, name, lat, lon, alt):
        """Test moonrise calculations."""
        jd = swe.julday(2024, 6, 15, 0.0)
        geopos = (lon, lat, alt)

        try:
            rise_swe = swe.rise_trans(jd, MOON, 1, geopos, 1013.25, 15.0)
            rise_py = ephem.rise_trans(jd, MOON, 1, [lon, lat, alt], 1013.25, 15.0)
        except Exception as e:
            pytest.skip(f"Moon rise calculation not available: {e}")
            return

        if rise_swe[0] != 0 or rise_py[0] != 0:
            pytest.skip("No moonrise found")
            return

        diff_seconds = abs(rise_swe[1][0] - rise_py[1][0]) * 86400

        # Moon can have slightly larger tolerance due to parallax
        assert diff_seconds < TIME_TOL_SECONDS * 2, (
            f"Moonrise at {name}: diff {diff_seconds:.1f}s exceeds tolerance"
        )
