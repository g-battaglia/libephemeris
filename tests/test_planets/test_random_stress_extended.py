"""
Random date stress tests for planetary positions.

Verifies that calc_ut produces valid, finite results at
many random dates across the DE440 range.
"""

from __future__ import annotations

import math

import numpy as np
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
    FLG_SPEED,
    FLG_SWIEPH,
    FLG_HELCTR,
    FLG_SIDEREAL,
    SIDM_LAHIRI,
)

# DE440 base range: ~1849 to ~2150
JD_1900 = 2415020.5
JD_2100 = 2488069.5


class TestRandomDateStress:
    """Stress test with random dates."""

    @pytest.mark.unit
    def test_sun_100_random_dates(self):
        """Sun position valid at 100 random dates."""
        rng = np.random.default_rng(42)
        jds = rng.uniform(JD_1900, JD_2100, 100)
        for jd in jds:
            result, _ = swe.calc_ut(float(jd), SUN, FLG_SPEED)
            assert 0 <= result[0] < 360, f"JD {jd}: lon={result[0]}"
            assert abs(result[1]) < 1, f"JD {jd}: lat={result[1]}"
            assert 0.98 < result[2] < 1.02, f"JD {jd}: dist={result[2]}"
            assert math.isfinite(result[3])

    @pytest.mark.unit
    def test_moon_100_random_dates(self):
        """Moon position valid at 100 random dates."""
        rng = np.random.default_rng(43)
        jds = rng.uniform(JD_1900, JD_2100, 100)
        for jd in jds:
            result, _ = swe.calc_ut(float(jd), MOON, FLG_SPEED)
            assert 0 <= result[0] < 360
            assert abs(result[1]) < 6  # Moon lat < ~5.3 deg
            assert 0.002 < result[2] < 0.003  # Moon dist ~0.0024-0.0028 AU
            assert math.isfinite(result[3])

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_planet_50_random_dates(self, body_id: int, name: str):
        """Planet positions valid at 50 random dates."""
        rng = np.random.default_rng(body_id * 100 + 44)
        jds = rng.uniform(JD_1900, JD_2100, 50)
        for jd in jds:
            result, _ = swe.calc_ut(float(jd), body_id, FLG_SPEED)
            assert 0 <= result[0] < 360, f"{name} JD {jd}: lon={result[0]}"
            assert abs(result[1]) < 10, f"{name} JD {jd}: lat={result[1]}"
            assert result[2] > 0, f"{name} JD {jd}: dist={result[2]}"


class TestRandomHeliocentricStress:
    """Stress test heliocentric at random dates."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_helio_30_random_dates(self, body_id: int, name: str):
        """Heliocentric positions valid at 30 random dates."""
        rng = np.random.default_rng(body_id * 200 + 55)
        jds = rng.uniform(JD_1900, JD_2100, 30)
        for jd in jds:
            result, _ = swe.calc_ut(float(jd), body_id, FLG_SPEED | FLG_HELCTR)
            assert 0 <= result[0] < 360, f"{name} helio lon={result[0]}"
            assert result[2] > 0, f"{name} helio dist={result[2]}"


class TestRandomSiderealStress:
    """Stress test sidereal positions at random dates."""

    @pytest.mark.unit
    def test_sidereal_sun_50_random_dates(self):
        """Sidereal Sun positions valid at 50 random dates."""
        swe.set_sid_mode(SIDM_LAHIRI)
        rng = np.random.default_rng(66)
        jds = rng.uniform(JD_1900, JD_2100, 50)
        for jd in jds:
            result, _ = swe.calc_ut(float(jd), SUN, FLG_SPEED | FLG_SIDEREAL)
            assert 0 <= result[0] < 360, f"Sid Sun lon={result[0]}"
            assert math.isfinite(result[3])


class TestRandomNodesStress:
    """Stress test lunar nodes at random dates."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (MEAN_NODE, "Mean Node"),
            (TRUE_NODE, "True Node"),
        ],
    )
    def test_nodes_50_random_dates(self, body_id: int, name: str):
        """Lunar node positions valid at 50 random dates."""
        rng = np.random.default_rng(body_id * 300 + 77)
        jds = rng.uniform(JD_1900, JD_2100, 50)
        for jd in jds:
            result, _ = swe.calc_ut(float(jd), body_id, FLG_SPEED)
            assert 0 <= result[0] < 360, f"{name} lon={result[0]}"
            # Node latitude is 0 by definition
            assert abs(result[1]) < 0.01, f"{name} lat={result[1]}"


class TestRandomHousesStress:
    """Stress test houses at random dates and locations."""

    @pytest.mark.unit
    def test_placidus_50_random_configs(self):
        """Placidus houses valid at 50 random date/location combos."""
        rng = np.random.default_rng(88)
        jds = rng.uniform(JD_1900, JD_2100, 50)
        lats = rng.uniform(-60, 60, 50)
        lons = rng.uniform(-180, 180, 50)

        for jd, lat, lon in zip(jds, lats, lons):
            cusps, ascmc = swe.houses(float(jd), float(lat), float(lon), ord("P"))
            assert len(cusps) >= 12
            for i, c in enumerate(cusps[:12]):
                assert 0 <= c < 360, (
                    f"JD={jd:.1f} lat={lat:.1f} lon={lon:.1f}: cusp {i + 1}={c}"
                )
            # ASC and MC should be valid
            assert 0 <= ascmc[0] < 360  # ASC
            assert 0 <= ascmc[1] < 360  # MC
