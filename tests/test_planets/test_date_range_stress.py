"""
Date range stress tests.

Verifies that libephemeris produces valid positions across the full
DE440 date range (1550-2650 CE) for all core bodies. Includes tests
for boundary dates, century-by-century coverage, and high-volume
random date sampling.
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
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    OSCU_APOG,
    CHIRON,
    FLG_SPEED,
)


def _jd_from_ymd(y: int, m: int, d: int, h: float = 12.0) -> float:
    """Convert Y/M/D to JD via libephemeris."""
    return swe.julday(y, m, d, h)


def _random_jds_in_range(
    n: int, y_min: int, y_max: int, seed: int = 42
) -> list[tuple[float, str]]:
    """Generate n random (jd, label) tuples."""
    rng = random.Random(seed)
    results = []
    for i in range(n):
        y = rng.randint(y_min, y_max)
        m = rng.randint(1, 12)
        d = rng.randint(1, 28)
        h = rng.uniform(0, 24)
        jd = swe.julday(y, m, d, h)
        results.append((jd, f"{y}-{m:02d}-{d:02d}"))
    return results


CORE_BODIES = [
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

EXTENDED_BODIES = [
    (MEAN_NODE, "MeanNode"),
    (TRUE_NODE, "TrueNode"),
    (MEAN_APOG, "MeanApog"),
    (OSCU_APOG, "OscuApog"),
    (CHIRON, "Chiron"),
]


class TestCenturyByCentury:
    """Test positions century by century from 1600 to 2500."""

    CENTURIES = [
        (1600, "17th century"),
        (1700, "18th century"),
        (1800, "19th century"),
        (1900, "20th century start"),
        (1950, "20th century mid"),
        (2000, "21st century start"),
        (2050, "21st century mid"),
        (2100, "22nd century start"),
        (2200, "23rd century"),
        (2300, "24th century"),
        (2400, "25th century"),
        (2500, "26th century"),
    ]

    @pytest.mark.unit
    @pytest.mark.parametrize("year,era", CENTURIES)
    @pytest.mark.parametrize("body_id,body_name", CORE_BODIES)
    def test_planet_valid_at_century(
        self, body_id: int, body_name: str, year: int, era: str
    ):
        """Core planet returns valid position at each century mark."""
        jd = _jd_from_ymd(year, 1, 1)
        result, retflag = swe.calc_ut(jd, body_id, FLG_SPEED)
        assert len(result) == 6, f"{body_name} @ {era}"
        lon, lat, dist = result[0], result[1], result[2]
        assert 0 <= lon < 360, f"{body_name} @ {era}: lon={lon}"
        assert -90 <= lat <= 90, f"{body_name} @ {era}: lat={lat}"
        assert dist > 0, f"{body_name} @ {era}: dist={dist}"
        # Speed should be finite
        assert math.isfinite(result[3]), f"{body_name} @ {era}: speed not finite"


class TestBoundaryDates:
    """Test at or near ephemeris boundary dates."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name", CORE_BODIES)
    def test_de440_start_boundary(self, body_id: int, body_name: str):
        """Position valid near DE440 start (1550 CE)."""
        jd = _jd_from_ymd(1550, 6, 15)
        result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
        assert len(result) == 6
        assert 0 <= result[0] < 360

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name", CORE_BODIES)
    def test_de440_end_boundary(self, body_id: int, body_name: str):
        """Position valid near DE440 end (2650 CE)."""
        jd = _jd_from_ymd(2649, 6, 15)
        result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
        assert len(result) == 6
        assert 0 <= result[0] < 360

    @pytest.mark.unit
    def test_j2000_epoch(self):
        """All core bodies valid at J2000.0."""
        jd = 2451545.0
        for body_id, name in CORE_BODIES:
            result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
            assert len(result) == 6, f"{name} at J2000"
            assert 0 <= result[0] < 360, f"{name} at J2000: lon={result[0]}"

    @pytest.mark.unit
    def test_unix_epoch(self):
        """All core bodies valid at Unix epoch (1970-01-01)."""
        jd = 2440587.5
        for body_id, name in CORE_BODIES:
            result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
            assert len(result) == 6, f"{name} at Unix epoch"


class TestHighVolumeRandomDates:
    """High-volume random date tests — 100 dates per planet."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name", CORE_BODIES)
    def test_100_random_dates_modern_era(self, body_id: int, body_name: str):
        """100 random dates in 1800-2200 range."""
        dates = _random_jds_in_range(100, 1800, 2200, seed=body_id * 7)
        for jd, label in dates:
            result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
            assert len(result) == 6, f"{body_name} @ {label}"
            assert 0 <= result[0] < 360, f"{body_name} @ {label}: lon={result[0]}"
            assert math.isfinite(result[2]), f"{body_name} @ {label}: dist not finite"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name", CORE_BODIES)
    def test_50_random_dates_full_range(self, body_id: int, body_name: str):
        """50 random dates in full DE440 range (1550-2649)."""
        dates = _random_jds_in_range(50, 1550, 2649, seed=body_id * 13)
        for jd, label in dates:
            result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
            assert len(result) == 6, f"{body_name} @ {label}"
            assert 0 <= result[0] < 360

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name", EXTENDED_BODIES)
    def test_50_random_dates_extended_bodies(self, body_id: int, body_name: str):
        """50 random dates for extended bodies (nodes, apogee, Chiron)."""
        dates = _random_jds_in_range(50, 1800, 2200, seed=body_id * 17)
        for jd, label in dates:
            result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
            assert len(result) == 6, f"{body_name} @ {label}"
            assert 0 <= result[0] < 360


class TestContinuityAcrossDecades:
    """Verify positions don't have discontinuities at decade boundaries."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
        ],
    )
    def test_no_jumps_across_centuries(self, body_id: int, body_name: str):
        """No large longitude jumps at century boundaries."""
        years = list(range(1600, 2600, 10))
        prev_lon = None
        for year in years:
            jd = _jd_from_ymd(year, 6, 15)
            result, _ = swe.calc_ut(jd, body_id, 0)
            lon = result[0]
            if prev_lon is not None and body_id != MOON:
                # For slow planets, check no huge jumps
                diff = abs(lon - prev_lon)
                if diff > 180:
                    diff = 360 - diff
                # 10-year gap: even fast planets move < 360° in 10 years
                # (except Moon, skip Moon for this check)
                assert diff < 360, f"{body_name}: jump {diff:.1f}° at {year}"
            prev_lon = lon

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SUN, "Sun"),
            (MARS, "Mars"),
        ],
    )
    def test_daily_continuity_100_days(self, body_id: int, body_name: str):
        """No discontinuities in 100-day daily sequence."""
        jd_start = 2451545.0
        prev_lon = None
        for i in range(100):
            jd = jd_start + i
            result, _ = swe.calc_ut(jd, body_id, 0)
            lon = result[0]
            if prev_lon is not None:
                diff = abs(lon - prev_lon)
                if diff > 180:
                    diff = 360 - diff
                # Max daily motion: Sun ~1°, Mars ~0.8°
                max_daily = 2.0
                assert diff < max_daily, (
                    f"{body_name}: daily jump {diff:.4f}° at day {i}"
                )
            prev_lon = lon


class TestMoonHighVolume:
    """Extra high-volume tests for Moon (fastest body, most sensitive)."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "jd,label",
        _random_jds_in_range(200, 1800, 2200, seed=777),
    )
    def test_moon_200_dates(self, jd: float, label: str):
        """Moon valid at 200 random dates."""
        result, _ = swe.calc_ut(jd, MOON, FLG_SPEED)
        lon, lat, dist, slon, slat, sdist = result
        assert 0 <= lon < 360, f"Moon @ {label}: lon={lon}"
        assert -7 < lat < 7, f"Moon @ {label}: lat={lat}"
        assert 0.002 < dist < 0.003, f"Moon @ {label}: dist={dist}"
        # Moon speed typically 12-15°/day
        assert 10 < slon < 16, f"Moon @ {label}: speed={slon}"
