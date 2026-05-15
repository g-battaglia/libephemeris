"""
Tests for heliocentric planet positions.

Verifies that FLG_HELCTR produces valid heliocentric coordinates
for all supported bodies across a wide date range. Checks output format,
value ranges, physical plausibility, and consistency between heliocentric
and geocentric results.
"""

from __future__ import annotations

import math
import random

import pytest

import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
    EARTH,
    CHIRON,
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    FLG_HELCTR,
    FLG_SPEED,
    FLG_EQUATORIAL,
    FLG_J2000,
    FLG_NOABERR,
    FLG_NOGDEFL,
)


def _random_jds(n: int, seed: int = 42) -> list[float]:
    """Generate n random JDs in the DE440 range (1800-2200)."""
    rng = random.Random(seed)
    jd_min = 2378497.0  # ~1800
    jd_max = 2524594.0  # ~2200
    return [rng.uniform(jd_min, jd_max) for _ in range(n)]


# Known heliocentric distance ranges (AU)
HELIO_DISTANCE_RANGES = {
    MERCURY: (0.30, 0.47),
    VENUS: (0.71, 0.73),
    EARTH: (0.98, 1.02),
    MARS: (1.38, 1.67),
    JUPITER: (4.95, 5.46),
    SATURN: (9.02, 10.05),
    URANUS: (18.3, 20.1),
    NEPTUNE: (29.8, 30.4),
    PLUTO: (29.6, 49.3),
}

# Approximate heliocentric orbital periods in days
HELIO_PERIODS = {
    MERCURY: 87.97,
    VENUS: 224.7,
    EARTH: 365.25,
    MARS: 687.0,
    JUPITER: 4332.6,
    SATURN: 10759.2,
    URANUS: 30688.5,
    NEPTUNE: 60182.0,
    PLUTO: 90560.0,
}

PLANETS = [
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (EARTH, "Earth"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
    (URANUS, "Uranus"),
    (NEPTUNE, "Neptune"),
    (PLUTO, "Pluto"),
]


class TestHeliocentricBasic:
    """Basic heliocentric position tests."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", PLANETS)
    def test_heliocentric_returns_valid_tuple(self, body_id: int, name: str):
        """Heliocentric calculation returns valid 6-element tuple."""
        jd = 2451545.0
        result, retflag = swe.calc_ut(jd, body_id, FLG_HELCTR | FLG_SPEED)
        assert len(result) == 6, f"{name}: expected 6 elements, got {len(result)}"
        for i, val in enumerate(result):
            assert math.isfinite(val), f"{name}: non-finite at index {i}"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", PLANETS)
    def test_heliocentric_longitude_range(self, body_id: int, name: str):
        """Heliocentric longitude should be in [0, 360)."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, FLG_HELCTR)
        assert 0 <= result[0] < 360, f"{name}: longitude {result[0]} out of range"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", PLANETS)
    def test_heliocentric_latitude_small(self, body_id: int, name: str):
        """Heliocentric latitude should be small (planets orbit near ecliptic)."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, FLG_HELCTR)
        # Pluto has high inclination (17°), others < 7°
        max_lat = 20.0 if body_id == PLUTO else 10.0
        assert abs(result[1]) < max_lat, f"{name}: latitude {result[1]}° too large"

    @pytest.mark.unit
    def test_sun_heliocentric_is_zero(self):
        """Sun heliocentric position should be (0,0,0,0,0,0)."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, SUN, FLG_HELCTR)
        for i, val in enumerate(result):
            assert val == pytest.approx(0.0, abs=1e-10), (
                f"Sun heliocentric result[{i}] = {val}, expected 0"
            )


class TestHeliocentricDistances:
    """Tests for heliocentric distance plausibility."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name,dist_min,dist_max",
        [
            (body_id, name, *HELIO_DISTANCE_RANGES[body_id])
            for body_id, name in PLANETS
            if body_id in HELIO_DISTANCE_RANGES
        ],
    )
    def test_heliocentric_distance_in_range(
        self, body_id: int, name: str, dist_min: float, dist_max: float
    ):
        """Heliocentric distance should be within known orbital range."""
        jds = _random_jds(20, seed=body_id)
        for jd in jds:
            result, _ = swe.calc_ut(jd, body_id, FLG_HELCTR)
            dist = result[2]
            assert dist_min * 0.9 <= dist <= dist_max * 1.1, (
                f"{name}: distance {dist} AU outside range "
                f"[{dist_min}, {dist_max}] at JD {jd}"
            )


class TestHeliocentricHighVolume:
    """High-volume heliocentric tests across many dates."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", PLANETS)
    def test_heliocentric_100_dates(self, body_id: int, name: str):
        """Heliocentric calculation succeeds for 100 random dates per planet."""
        jds = _random_jds(100, seed=body_id * 100)
        for jd in jds:
            result, retflag = swe.calc_ut(jd, body_id, FLG_HELCTR | FLG_SPEED)
            assert len(result) == 6
            assert 0 <= result[0] < 360
            assert math.isfinite(result[2])

    @pytest.mark.unit
    @pytest.mark.parametrize("jd", _random_jds(100, seed=9999))
    def test_earth_heliocentric_100_dates(self, jd: float):
        """Earth heliocentric position at 100 dates — distance ~1 AU."""
        result, _ = swe.calc_ut(jd, EARTH, FLG_HELCTR)
        assert 0.98 < result[2] < 1.02, (
            f"Earth heliocentric distance {result[2]} at JD {jd}"
        )


class TestHeliocentricFlagCombinations:
    """Tests for heliocentric with various flag combinations."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "extra_flags,desc",
        [
            (0, "helctr only"),
            (FLG_SPEED, "helctr+speed"),
            (FLG_EQUATORIAL, "helctr+equatorial"),
            (FLG_J2000, "helctr+J2000"),
            (FLG_NOABERR, "helctr+noaberr"),
            (FLG_NOGDEFL, "helctr+nogdefl"),
            (FLG_SPEED | FLG_J2000, "helctr+speed+J2000"),
            (FLG_EQUATORIAL | FLG_SPEED, "helctr+equatorial+speed"),
        ],
    )
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (EARTH, "Earth"),
        ],
    )
    def test_heliocentric_flag_combo(
        self, body_id: int, name: str, extra_flags: int, desc: str
    ):
        """Heliocentric works with various flag combinations."""
        jd = 2451545.0
        flags = FLG_HELCTR | extra_flags
        result, retflag = swe.calc_ut(jd, body_id, flags)
        assert len(result) == 6
        for i, val in enumerate(result):
            assert math.isfinite(val), f"{name} {desc}: non-finite at index {i}"


class TestHeliocentricVsGeocentric:
    """Tests comparing heliocentric and geocentric outputs."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_heliocentric_distance_differs_from_geocentric(
        self, body_id: int, name: str
    ):
        """Heliocentric and geocentric distances should differ."""
        jd = 2451545.0
        r_helio, _ = swe.calc_ut(jd, body_id, FLG_HELCTR)
        r_geo, _ = swe.calc_ut(jd, body_id, 0)
        assert abs(r_helio[2] - r_geo[2]) > 0.1, (
            f"{name}: helio dist {r_helio[2]} == geo dist {r_geo[2]}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
        ],
    )
    def test_heliocentric_longitude_differs_from_geocentric(
        self, body_id: int, name: str
    ):
        """Heliocentric and geocentric longitudes should generally differ."""
        jd = 2451545.0
        r_helio, _ = swe.calc_ut(jd, body_id, FLG_HELCTR)
        r_geo, _ = swe.calc_ut(jd, body_id, 0)
        # They might coincidentally be close, but generally differ
        diff = abs(r_helio[0] - r_geo[0])
        if diff > 180:
            diff = 360 - diff
        # Just verify they're different (not identical)
        assert diff > 0.001, f"{name}: helio and geo longitudes suspiciously identical"

    @pytest.mark.unit
    def test_earth_heliocentric_is_sun_geocentric_opposite(self):
        """Earth helio longitude should be ~180° from Sun geocentric."""
        # Force Skyfield mode to avoid LEB state leakage from other tests
        from libephemeris import state

        saved_leb = state._LEB_FILE
        saved_reader = state._LEB_READER
        state._LEB_FILE = None
        state._LEB_READER = None
        try:
            jd = 2451545.0
            r_earth, _ = swe.calc_ut(jd, EARTH, FLG_HELCTR)
            r_sun, _ = swe.calc_ut(jd, SUN, 0)
        finally:
            state._LEB_FILE = saved_leb
            state._LEB_READER = saved_reader
        diff = abs(r_earth[0] - r_sun[0])
        if diff > 180:
            diff = 360 - diff
        # Earth helio and Sun geo should be roughly opposite
        assert abs(diff - 180) < 1.0, (
            f"Earth helio {r_earth[0]:.2f} vs Sun geo {r_sun[0]:.2f} diff {diff:.2f}°"
        )


class TestHeliocentricSpeed:
    """Tests for heliocentric speed values."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", PLANETS)
    def test_heliocentric_speed_nonzero(self, body_id: int, name: str):
        """Heliocentric speed should be non-zero for all planets."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, FLG_HELCTR | FLG_SPEED)
        speed = result[3]
        assert abs(speed) > 1e-6, (
            f"{name}: heliocentric speed {speed} is essentially zero"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name,expected_speed",
        [
            (MERCURY, "Mercury", 360 / 87.97),
            (VENUS, "Venus", 360 / 224.7),
            (EARTH, "Earth", 360 / 365.25),
            (MARS, "Mars", 360 / 687.0),
            (JUPITER, "Jupiter", 360 / 4332.6),
            (SATURN, "Saturn", 360 / 10759.2),
        ],
    )
    def test_heliocentric_speed_order_of_magnitude(
        self, body_id: int, name: str, expected_speed: float
    ):
        """Heliocentric speed should be roughly 360°/period."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, FLG_HELCTR | FLG_SPEED)
        speed = result[3]
        # Allow factor of 2 variation (eccentricity)
        assert expected_speed * 0.5 < abs(speed) < expected_speed * 2.0, (
            f"{name}: speed {speed}°/day vs expected ~{expected_speed:.4f}"
        )
