"""
Comprehensive tests for Uranian hypothetical bodies (IDs 40-48).

Verifies that all 9 Uranian bodies (Cupido through Transpluto/Isis)
return valid positions across multiple dates, flag combinations,
and in both geocentric and heliocentric modes.
"""

from __future__ import annotations

import math
import random

import pytest

import libephemeris as swe
from libephemeris.constants import (
    FLG_SPEED,
    FLG_HELCTR,
    FLG_EQUATORIAL,
    FLG_J2000,
    FLG_NOABERR,
    FLG_SIDEREAL,
    SIDM_LAHIRI,
    SIDM_FAGAN_BRADLEY,
)


URANIAN_BODIES = [
    (40, "Cupido"),
    (41, "Hades"),
    (42, "Zeus"),
    (43, "Kronos"),
    (44, "Apollon"),
    (45, "Admetos"),
    (46, "Vulkanus"),
    (47, "Poseidon"),
    (48, "Transpluto"),
]

# Approximate heliocentric distances (AU) for Uranian bodies
URANIAN_DISTANCES = {
    40: (35, 55),
    41: (40, 60),
    42: (50, 70),
    43: (55, 75),
    44: (60, 80),
    45: (65, 85),
    46: (68, 88),
    47: (75, 95),
    48: (80, 105),
}


def _random_jds(n: int, seed: int = 42) -> list[float]:
    """Generate n random JDs in 1900-2100."""
    rng = random.Random(seed)
    return [rng.uniform(2415020.0, 2488070.0) for _ in range(n)]


class TestUranianBasic:
    """Basic functionality tests for Uranian bodies."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", URANIAN_BODIES)
    def test_uranian_returns_valid_position(self, body_id: int, name: str):
        """Each Uranian body returns valid 6-element tuple."""
        jd = 2451545.0
        result, retflag = swe.calc_ut(jd, body_id, FLG_SPEED)
        assert len(result) == 6, f"{name}: expected 6 elements"
        lon, lat, dist = result[0], result[1], result[2]
        assert 0 <= lon < 360, f"{name}: lon {lon} out of range"
        assert -90 <= lat <= 90, f"{name}: lat {lat} out of range"
        assert dist > 0, f"{name}: distance {dist} not positive"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", URANIAN_BODIES)
    def test_uranian_speed_nonzero(self, body_id: int, name: str):
        """Uranian body speeds should be non-zero."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
        speed = result[3]
        assert abs(speed) > 1e-8, f"{name}: speed {speed} is zero"
        # Uranian bodies are very distant, speeds should be small
        assert abs(speed) < 0.1, f"{name}: speed {speed} too large"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", URANIAN_BODIES)
    def test_uranian_latitude_small(self, body_id: int, name: str):
        """Uranian body latitude should be small (hypothetical near-ecliptic)."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, 0)
        lat = result[1]
        # Transpluto has exactly 0 latitude; others have small inclinations
        assert abs(lat) < 5.0, f"{name}: latitude {lat}° too large"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", URANIAN_BODIES)
    def test_uranian_all_finite(self, body_id: int, name: str):
        """All output values should be finite."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
        for i, val in enumerate(result):
            assert math.isfinite(val), f"{name}: result[{i}] = {val} not finite"


class TestUranianDistances:
    """Test Uranian body distances are physically plausible."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", URANIAN_BODIES)
    def test_uranian_distance_range(self, body_id: int, name: str):
        """Uranian body geocentric distance should be in expected range."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, 0)
        dist = result[2]
        lo, hi = URANIAN_DISTANCES[body_id]
        assert lo * 0.7 <= dist <= hi * 1.3, (
            f"{name}: distance {dist} AU outside [{lo}, {hi}]"
        )


class TestUranianHeliocentric:
    """Test Uranian bodies in heliocentric mode."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", URANIAN_BODIES)
    def test_uranian_heliocentric_valid(self, body_id: int, name: str):
        """Heliocentric positions valid for all Uranian bodies."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, FLG_HELCTR | FLG_SPEED)
        assert len(result) == 6
        assert 0 <= result[0] < 360
        assert math.isfinite(result[2])

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", URANIAN_BODIES)
    def test_uranian_helio_vs_geo_distance(self, body_id: int, name: str):
        """Helio and geo distances should differ (parallax from Earth)."""
        jd = 2451545.0
        r_geo, _ = swe.calc_ut(jd, body_id, 0)
        r_helio, _ = swe.calc_ut(jd, body_id, FLG_HELCTR)
        # For very distant bodies, the difference is small but nonzero
        assert abs(r_geo[2] - r_helio[2]) > 0.01, (
            f"{name}: geo dist {r_geo[2]} == helio dist {r_helio[2]}"
        )


class TestUranianFlagCombinations:
    """Test Uranian bodies with various flag combinations."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "flags,desc",
        [
            (0, "default"),
            (FLG_SPEED, "speed"),
            (FLG_EQUATORIAL, "equatorial"),
            (FLG_J2000, "J2000"),
            (FLG_NOABERR, "no aberration"),
            (FLG_HELCTR, "heliocentric"),
            (FLG_SPEED | FLG_EQUATORIAL, "speed+equatorial"),
            (FLG_HELCTR | FLG_SPEED, "helio+speed"),
        ],
    )
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (40, "Cupido"),
            (44, "Apollon"),
            (48, "Transpluto"),
        ],
    )
    def test_uranian_flag_combo(self, body_id: int, name: str, flags: int, desc: str):
        """Uranian body works with various flag combinations."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, flags)
        assert len(result) == 6
        for i, val in enumerate(result):
            assert math.isfinite(val), f"{name}+{desc}: result[{i}] = {val}"


class TestUranianSidereal:
    """Test Uranian bodies in sidereal mode."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", URANIAN_BODIES)
    def test_uranian_sidereal_valid(self, body_id: int, name: str):
        """Sidereal positions valid for all Uranian bodies."""
        swe.set_sid_mode(SIDM_LAHIRI)
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, FLG_SIDEREAL)
        assert 0 <= result[0] < 360, f"{name}: sidereal lon {result[0]}"


class TestUranianDateRange:
    """Test Uranian bodies across date ranges."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", URANIAN_BODIES)
    def test_uranian_across_centuries(self, body_id: int, name: str):
        """Uranian bodies valid from 1600 to 2500."""
        for year in [1600, 1800, 1900, 2000, 2100, 2300, 2500]:
            jd = swe.julday(year, 1, 1, 12.0)
            result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
            assert 0 <= result[0] < 360, f"{name} @ {year}: lon={result[0]}"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", URANIAN_BODIES)
    def test_uranian_50_random_dates(self, body_id: int, name: str):
        """Uranian bodies valid at 50 random dates."""
        jds = _random_jds(50, seed=body_id * 37)
        for jd in jds:
            result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
            assert 0 <= result[0] < 360
            assert math.isfinite(result[3])

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (40, "Cupido"),
            (44, "Apollon"),
            (48, "Transpluto"),
        ],
    )
    def test_uranian_continuity(self, body_id: int, name: str):
        """Uranian body positions should be continuous."""
        jd_start = 2451545.0
        prev_lon = None
        for i in range(100):
            jd = jd_start + i * 30.0  # monthly
            result, _ = swe.calc_ut(jd, body_id, 0)
            lon = result[0]
            if prev_lon is not None:
                diff = abs(lon - prev_lon)
                if diff > 180:
                    diff = 360 - diff
                # Uranian bodies are slow: < 3°/month
                assert diff < 5.0, f"{name}: jump {diff:.2f}° at month {i}"
            prev_lon = lon


class TestUranianHighVolume:
    """High-volume Uranian body tests."""

    @pytest.mark.unit
    @pytest.mark.parametrize("jd", _random_jds(100, seed=4444))
    def test_cupido_100_dates(self, jd: float):
        """Cupido valid at 100 random dates."""
        result, _ = swe.calc_ut(jd, 40, FLG_SPEED)
        assert 0 <= result[0] < 360
        assert abs(result[3]) < 0.1

    @pytest.mark.unit
    @pytest.mark.parametrize("jd", _random_jds(100, seed=5555))
    def test_transpluto_100_dates(self, jd: float):
        """Transpluto valid at 100 random dates."""
        result, _ = swe.calc_ut(jd, 48, FLG_SPEED)
        assert 0 <= result[0] < 360
        # Transpluto latitude is near 0 but not always exactly 0
        assert abs(result[1]) < 0.02
