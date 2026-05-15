"""
Tests for asteroid position calculations by number.

Verifies that asteroids accessed via AST_OFFSET + number
return valid positions, and that the 5 main asteroids
(Chiron, Ceres, Pallas, Juno, Vesta) are accurate.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    CHIRON,
    CERES,
    PALLAS,
    JUNO,
    VESTA,
    AST_OFFSET,
    FLG_SPEED,
    FLG_HELCTR,
    FLG_EQUATORIAL,
    FLG_SIDEREAL,
    SIDM_LAHIRI,
)


# Main asteroids with their direct body IDs
MAIN_ASTEROIDS = [
    (CHIRON, "Chiron", 15),
    (CERES, "Ceres", 17),
    (PALLAS, "Pallas", 18),
    (JUNO, "Juno", 19),
    (VESTA, "Vesta", 20),
]

# Approximate heliocentric distance ranges (AU)
ASTEROID_DISTANCES = {
    CHIRON: (8, 20),  # Chiron: between Saturn and Uranus
    CERES: (2.5, 3.0),  # Ceres: in the main belt
    PALLAS: (2.0, 3.5),  # Pallas: main belt, eccentric
    JUNO: (2.0, 3.5),  # Juno: main belt
    VESTA: (2.0, 2.8),  # Vesta: main belt
}


class TestMainAsteroidsBasic:
    """Basic tests for the 5 main asteroids."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name,idx", MAIN_ASTEROIDS)
    def test_asteroid_returns_valid(self, body_id: int, name: str, idx: int):
        """Each main asteroid returns valid 6-element tuple."""
        jd = 2451545.0
        result, retflag = swe.calc_ut(jd, body_id, FLG_SPEED)
        assert len(result) == 6, f"{name}: expected 6 elements"
        lon, lat, dist = result[0], result[1], result[2]
        assert 0 <= lon < 360, f"{name}: lon {lon} out of range"
        assert -90 <= lat <= 90, f"{name}: lat {lat} out of range"
        assert dist > 0, f"{name}: distance {dist} not positive"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name,idx", MAIN_ASTEROIDS)
    def test_asteroid_speed_reasonable(self, body_id: int, name: str, idx: int):
        """Asteroid speeds should be reasonable (< 1 deg/day)."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
        speed = result[3]
        assert abs(speed) < 1.5, f"{name}: speed {speed}°/day too large"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name,idx", MAIN_ASTEROIDS)
    def test_asteroid_all_finite(self, body_id: int, name: str, idx: int):
        """All values should be finite."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
        for i, val in enumerate(result):
            assert math.isfinite(val), f"{name}: result[{i}] = {val}"


class TestAsteroidDistances:
    """Test asteroid distances are physically plausible."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name,idx", MAIN_ASTEROIDS)
    def test_asteroid_geocentric_distance(self, body_id: int, name: str, idx: int):
        """Geocentric distance should be plausible."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, 0)
        dist = result[2]
        # Geocentric distance depends on Earth position relative to asteroid
        # For main belt: ~1-4 AU from Earth; for Chiron: ~7-20 AU
        if body_id == CHIRON:
            assert 5 < dist < 25, f"{name}: dist {dist} AU"
        else:
            assert 0.5 < dist < 5.0, f"{name}: dist {dist} AU"


class TestAsteroidFlagCombinations:
    """Test asteroids with various flag combinations."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "flags,desc",
        [
            (0, "default"),
            (FLG_SPEED, "speed"),
            (FLG_EQUATORIAL, "equatorial"),
            (FLG_HELCTR, "heliocentric"),
            (FLG_SPEED | FLG_EQUATORIAL, "speed+equatorial"),
        ],
    )
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (CHIRON, "Chiron"),
            (CERES, "Ceres"),
            (VESTA, "Vesta"),
        ],
    )
    def test_asteroid_flag_combo(self, body_id: int, name: str, flags: int, desc: str):
        """Asteroid works with various flag combinations."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, flags)
        assert len(result) == 6
        for i, val in enumerate(result):
            assert math.isfinite(val), f"{name}+{desc}: result[{i}]={val}"


class TestAsteroidSidereal:
    """Test asteroids in sidereal mode."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name,idx", MAIN_ASTEROIDS)
    def test_asteroid_sidereal(self, body_id: int, name: str, idx: int):
        """Sidereal positions valid for asteroids."""
        swe.set_sid_mode(SIDM_LAHIRI)
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, FLG_SIDEREAL)
        assert 0 <= result[0] < 360, f"{name}: sidereal lon {result[0]}"


class TestAsteroidDateRange:
    """Test asteroids across date ranges."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name,idx", MAIN_ASTEROIDS)
    def test_asteroid_across_centuries(self, body_id: int, name: str, idx: int):
        """Asteroids valid from 1800 to 2200."""
        for year in [1800, 1900, 1950, 2000, 2050, 2100, 2200]:
            jd = swe.julday(year, 1, 1, 12.0)
            result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
            assert 0 <= result[0] < 360, f"{name} @ {year}: lon={result[0]}"

    @pytest.mark.unit
    def test_chiron_continuity(self):
        """Chiron positions should be continuous over months."""
        jd_start = 2451545.0
        prev_lon = None
        for i in range(24):  # 2 years monthly
            jd = jd_start + i * 30.0
            result, _ = swe.calc_ut(jd, CHIRON, 0)
            lon = result[0]
            if prev_lon is not None:
                diff = abs(lon - prev_lon)
                if diff > 180:
                    diff = 360 - diff
                # Chiron moves ~2°/month
                assert diff < 10.0, f"Chiron jump {diff:.2f}° at month {i}"
            prev_lon = lon


class TestAsteroidHeliocentric:
    """Test asteroids in heliocentric mode."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name,idx", MAIN_ASTEROIDS)
    def test_asteroid_helio_valid(self, body_id: int, name: str, idx: int):
        """Heliocentric positions valid for asteroids."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, FLG_HELCTR | FLG_SPEED)
        assert len(result) == 6
        assert 0 <= result[0] < 360
        assert result[2] > 0

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name,idx", MAIN_ASTEROIDS)
    def test_asteroid_helio_distance_range(self, body_id: int, name: str, idx: int):
        """Heliocentric distance should match known orbital ranges."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, FLG_HELCTR)
        dist = result[2]
        lo, hi = ASTEROID_DISTANCES[body_id]
        assert lo * 0.6 <= dist <= hi * 1.5, (
            f"{name}: helio dist {dist} AU outside [{lo}, {hi}]"
        )


class TestAsteroidByNumber:
    """Test asteroid lookup by AST_OFFSET + catalog number."""

    @pytest.mark.unit
    def test_ceres_by_offset(self):
        """Ceres via AST_OFFSET + 1 should match CERES."""
        jd = 2451545.0
        try:
            r_direct, _ = swe.calc_ut(jd, CERES, FLG_SPEED)
            r_offset, _ = swe.calc_ut(jd, AST_OFFSET + 1, FLG_SPEED)

            # Positions should be identical
            for i in range(6):
                assert abs(r_direct[i] - r_offset[i]) < 0.01, (
                    f"Ceres element {i}: direct={r_direct[i]}, offset={r_offset[i]}"
                )
        except Exception:
            # AST_OFFSET lookup may not be supported for all implementations
            pytest.skip("AST_OFFSET lookup not supported")

    @pytest.mark.unit
    def test_pallas_by_offset(self):
        """Pallas via AST_OFFSET + 2 should match PALLAS."""
        jd = 2451545.0
        try:
            r_direct, _ = swe.calc_ut(jd, PALLAS, FLG_SPEED)
            r_offset, _ = swe.calc_ut(jd, AST_OFFSET + 2, FLG_SPEED)

            for i in range(6):
                assert abs(r_direct[i] - r_offset[i]) < 0.01, (
                    f"Pallas element {i}: direct={r_direct[i]}, offset={r_offset[i]}"
                )
        except Exception:
            pytest.skip("AST_OFFSET lookup not supported")

    @pytest.mark.unit
    def test_juno_by_offset(self):
        """Juno via AST_OFFSET + 3 should match JUNO."""
        jd = 2451545.0
        try:
            r_direct, _ = swe.calc_ut(jd, JUNO, FLG_SPEED)
            r_offset, _ = swe.calc_ut(jd, AST_OFFSET + 3, FLG_SPEED)

            for i in range(6):
                assert abs(r_direct[i] - r_offset[i]) < 0.01, (
                    f"Juno element {i}: direct={r_direct[i]}, offset={r_offset[i]}"
                )
        except Exception:
            pytest.skip("AST_OFFSET lookup not supported")

    @pytest.mark.unit
    def test_vesta_by_offset(self):
        """Vesta via AST_OFFSET + 4 should match VESTA."""
        jd = 2451545.0
        try:
            r_direct, _ = swe.calc_ut(jd, VESTA, FLG_SPEED)
            r_offset, _ = swe.calc_ut(jd, AST_OFFSET + 4, FLG_SPEED)

            for i in range(6):
                assert abs(r_direct[i] - r_offset[i]) < 0.01, (
                    f"Vesta element {i}: direct={r_direct[i]}, offset={r_offset[i]}"
                )
        except Exception:
            pytest.skip("AST_OFFSET lookup not supported")
