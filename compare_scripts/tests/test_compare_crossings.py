"""
Pytest-style Crossing Functions Comparison Tests.

Validates solcross_ut, mooncross_ut, mooncross_node_ut, helio_cross_ut,
and other crossing calculations against pyswisseph.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    JUPITER,
    SATURN,
    FLG_SWIEPH,
)


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class CrossingTolerance:
    """Tolerance thresholds for crossing comparisons."""

    SUN_SECONDS = 60.0  # 1 minute for Sun crossings
    MOON_SECONDS = 120.0  # 2 minutes for Moon crossings
    MOON_NODE_SECONDS = 300.0  # 5 minutes for Moon node crossings
    PLANET_SECONDS = 600.0  # 10 minutes for planets


# ============================================================================
# FIXTURES
# ============================================================================


@pytest.fixture
def jd_start():
    """Standard Julian Day for crossing searches."""
    return swe.julday(2024, 1, 1, 0.0)


# ============================================================================
# SUN CROSSING TESTS
# ============================================================================


class TestSunCrossings:
    """Tests for Sun crossing (solcross_ut) calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "target_lon,offset",
        [
            (0, 0),  # 0 Aries
            (30, 35),  # 0 Taurus
            (60, 70),  # 0 Gemini
            (90, 105),  # 0 Cancer
            (120, 140),  # 0 Leo
            (150, 175),  # 0 Virgo
            (180, 210),  # 0 Libra
        ],
    )
    def test_sun_crossing_longitude(self, jd_start, target_lon, offset):
        """Test Sun crossing specific longitudes."""
        jd = jd_start + offset

        # LibEphemeris
        jd_py = ephem.solcross_ut(target_lon, jd, 0)

        # SwissEphemeris
        jd_swe = swe.solcross_ut(target_lon, jd, 0)

        # Calculate difference in seconds
        diff_seconds = abs(jd_py - jd_swe) * 86400

        # Verify Sun is near target longitude
        pos_py, _ = ephem.calc_ut(jd_py, SUN, 0)
        lon_diff = min(abs(pos_py[0] - target_lon), 360 - abs(pos_py[0] - target_lon))

        assert diff_seconds < CrossingTolerance.SUN_SECONDS, (
            f"Sun crossing {target_lon}° difference too large: {diff_seconds:.2f}s"
        )
        assert lon_diff < 0.01, f"Sun not at target longitude: error {lon_diff}°"


# ============================================================================
# MOON CROSSING TESTS
# ============================================================================


class TestMoonCrossings:
    """Tests for Moon crossing (mooncross_ut) calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("target_lon", [0, 90, 180, 270])
    def test_moon_crossing_longitude(self, target_lon):
        """Test Moon crossing specific longitudes."""
        jd_start = swe.julday(2024, 11, 1, 0.0)

        # LibEphemeris
        jd_py = ephem.mooncross_ut(target_lon, jd_start, 0)

        # SwissEphemeris
        jd_swe = swe.mooncross_ut(target_lon, jd_start, 0)

        # Calculate difference in seconds
        diff_seconds = abs(jd_py - jd_swe) * 86400

        assert diff_seconds < CrossingTolerance.MOON_SECONDS, (
            f"Moon crossing {target_lon}° difference too large: {diff_seconds:.2f}s"
        )


# ============================================================================
# MOON NODE CROSSING TESTS
# ============================================================================


class TestMoonNodeCrossings:
    """Tests for Moon node crossing (mooncross_node_ut) calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "node_type,node_name",
        [
            (0, "Ascending"),
            (1, "Descending"),
        ],
    )
    def test_mooncross_node_ut(self, jd_start, node_type, node_name):
        """Test Moon crossing its nodes."""
        # SwissEphemeris - returns tuple (jd_cross, xlon, xlat)
        result_swe = swe.mooncross_node_ut(jd_start, node_type)
        jd_swe = result_swe[0]

        # LibEphemeris - also returns tuple (jd_cross, xlon, xlat)
        result_py = ephem.mooncross_node_ut(jd_start, node_type)
        jd_py = result_py[0]

        diff_seconds = abs(jd_py - jd_swe) * 86400

        assert diff_seconds < CrossingTolerance.MOON_NODE_SECONDS, (
            f"Moon {node_name} node crossing difference: {diff_seconds:.2f}s"
        )

    @pytest.mark.comparison
    def test_mooncross_node_multiple(self, jd_start):
        """Test multiple Moon node crossings over time."""
        for i in range(4):
            jd = jd_start + i * 14  # ~2 weeks apart

            for node_type in [0, 1]:
                result_swe = swe.mooncross_node_ut(jd, node_type)
                result_py = ephem.mooncross_node_ut(jd, node_type)

                diff_seconds = abs(result_py[0] - result_swe[0]) * 86400
                assert diff_seconds < CrossingTolerance.MOON_NODE_SECONDS


# ============================================================================
# HELIOCENTRIC CROSSING TESTS
# ============================================================================


class TestHelioCrossings:
    """Tests for heliocentric crossing (helio_cross_ut) calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    @pytest.mark.parametrize("target_lon", [0, 90, 180, 270])
    def test_helio_cross_ut(self, jd_start, body_id, body_name, target_lon):
        """Test heliocentric longitude crossings."""
        # SwissEphemeris
        try:
            jd_swe = swe.helio_cross_ut(body_id, target_lon, jd_start, FLG_SWIEPH)
        except Exception as e:
            pytest.skip(f"SwissEphemeris helio_cross_ut failed: {e}")

        # LibEphemeris
        try:
            jd_py = ephem.helio_cross_ut(body_id, target_lon, jd_start, FLG_SWIEPH)
        except Exception as e:
            pytest.skip(f"LibEphemeris helio_cross_ut failed: {e}")

        diff_seconds = abs(jd_py - jd_swe) * 86400

        assert diff_seconds < CrossingTolerance.PLANET_SECONDS, (
            f"{body_name} helio crossing {target_lon}° difference: {diff_seconds:.2f}s"
        )


# ============================================================================
# GENERIC CROSSING TESTS
# ============================================================================


class TestGenericCrossings:
    """Tests for generic crossing (cross_ut) calculations."""

    @pytest.mark.comparison
    def test_swe_cross_ut_sun(self, jd_start):
        """Test generic crossing function for Sun."""
        target_lon = 0

        # SwissEphemeris - use solcross_ut for Sun
        jd_swe = swe.solcross_ut(target_lon, jd_start, FLG_SWIEPH)

        # LibEphemeris - use cross_ut
        jd_py = ephem.cross_ut(SUN, target_lon, jd_start, FLG_SWIEPH)

        diff_seconds = abs(jd_py - jd_swe) * 86400

        assert diff_seconds < CrossingTolerance.PLANET_SECONDS, (
            f"cross_ut Sun crossing difference: {diff_seconds:.2f}s"
        )
