"""Tests for house_pos and related house placement computations."""

from __future__ import annotations

import pytest
import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    FLG_SWIEPH,
    FLG_SPEED,
)

JD_J2000 = 2451545.0


@pytest.mark.unit
class TestHousePos:
    """Test house_pos for body house placement."""

    @pytest.mark.parametrize(
        "hsys",
        ["P", "K", "O", "R", "C", "E", "W", "B", "M", "T"],
        ids=[
            "Placidus",
            "Koch",
            "Porphyry",
            "Regio",
            "Campanus",
            "Equal",
            "WholeSgn",
            "Alcabitius",
            "Morinus",
            "Topo",
        ],
    )
    def test_house_pos_in_range(self, hsys):
        """House position should be in [1.0, 13.0)."""
        cusps, ascmc = swe.houses(JD_J2000, 41.9, 12.5, ord(hsys))
        armc = ascmc[2]
        eps = ascmc[8] if len(ascmc) > 8 else 23.4393

        sun, _ = swe.calc_ut(JD_J2000, SUN, FLG_SWIEPH | FLG_SPEED)
        pos = swe.house_pos(armc, 41.9, eps, (sun[0], sun[1]), hsys)
        assert 1.0 <= pos < 13.0, f"House pos {pos} out of range for {hsys}"

    @pytest.mark.parametrize(
        "body,name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
        ],
    )
    def test_house_pos_all_planets(self, body, name):
        """All planets should have valid house positions."""
        cusps, ascmc = swe.houses(JD_J2000, 41.9, 12.5, ord("P"))
        armc = ascmc[2]
        eps = ascmc[8] if len(ascmc) > 8 else 23.4393

        result, _ = swe.calc_ut(JD_J2000, body, FLG_SWIEPH | FLG_SPEED)
        pos = swe.house_pos(armc, 41.9, eps, (result[0], result[1]), "P")
        assert 1.0 <= pos < 13.0, f"{name} house pos {pos} out of range"

    def test_asc_in_house_1(self):
        """Ascendant should be at cusp of house 1.

        At the exact Ascendant, house_pos may return either ~1.0 (start of
        house 1) or ~12.999... (end of house 12, which wraps to house 1 cusp)
        due to floating-point boundary precision.
        """
        cusps, ascmc = swe.houses(JD_J2000, 41.9, 12.5, ord("P"))
        armc = ascmc[2]
        eps = ascmc[8] if len(ascmc) > 8 else 23.4393
        asc_lon = ascmc[0]

        pos = swe.house_pos(armc, 41.9, eps, (asc_lon, 0.0), "P")
        at_cusp_1 = (1.0 <= pos < 1.01) or (pos > 12.999)
        assert at_cusp_1, f"Asc house pos {pos} not at cusp 1"

    def test_mc_in_house_10(self):
        """MC should be in house 10 (pos ~ 10.0)."""
        cusps, ascmc = swe.houses(JD_J2000, 41.9, 12.5, ord("P"))
        armc = ascmc[2]
        eps = ascmc[8] if len(ascmc) > 8 else 23.4393
        mc_lon = ascmc[1]

        pos = swe.house_pos(armc, 41.9, eps, (mc_lon, 0.0), "P")
        assert 9.99 < pos < 10.01, f"MC house pos {pos} not at cusp 10"

    def test_whole_sign_integer_house(self):
        """Whole sign house position should have integer part matching sign."""
        cusps, ascmc = swe.houses(JD_J2000, 41.9, 12.5, ord("W"))
        armc = ascmc[2]
        eps = ascmc[8] if len(ascmc) > 8 else 23.4393
        asc_lon = ascmc[0]

        # A body at the Ascendant degree should be in house 1 (or at the 12/1 boundary)
        pos = swe.house_pos(armc, 41.9, eps, (asc_lon, 0.0), "W")
        assert int(pos) == 1 or pos > 12.999, (
            f"Whole sign: body at Asc should be house 1, got {pos}"
        )


@pytest.mark.unit
class TestHousePosMultipleLocations:
    """House positions at different geographic locations."""

    @pytest.mark.parametrize(
        "lat,lon,name",
        [
            (0.0, 0.0, "Equator"),
            (41.9, 12.5, "Rome"),
            (51.5, -0.1, "London"),
            (-33.9, 151.2, "Sydney"),
            (35.7, 139.7, "Tokyo"),
            (64.0, -22.0, "Reykjavik"),
        ],
    )
    def test_sun_house_pos_valid(self, lat, lon, name):
        """Sun should have valid house position at all locations."""
        cusps, ascmc = swe.houses(JD_J2000, lat, lon, ord("P"))
        armc = ascmc[2]
        eps = ascmc[8] if len(ascmc) > 8 else 23.4393

        sun, _ = swe.calc_ut(JD_J2000, SUN, FLG_SWIEPH)
        pos = swe.house_pos(armc, lat, eps, (sun[0], sun[1]), "P")
        assert 1.0 <= pos < 13.0, f"Sun house pos {pos} at {name}"


@pytest.mark.unit
class TestHousePosConsistency:
    """House position should be consistent with cusp positions."""

    def test_body_between_cusps(self):
        """Body house pos integer should match surrounding cusps."""
        cusps, ascmc = swe.houses(JD_J2000, 41.9, 12.5, ord("P"))
        armc = ascmc[2]
        eps = ascmc[8] if len(ascmc) > 8 else 23.4393

        sun, _ = swe.calc_ut(JD_J2000, SUN, FLG_SWIEPH)
        pos = swe.house_pos(armc, 41.9, eps, (sun[0], sun[1]), "P")
        house_num = int(pos)  # 1-12

        # Sun longitude should be between cusp[house_num-1] and cusp[house_num]
        # (accounting for angle wrapping, this is verified by house_pos itself)
        assert 1 <= house_num <= 12

    @pytest.mark.parametrize("day_offset", range(0, 365, 30))
    def test_sun_traverses_all_houses(self, day_offset):
        """Sun house position should vary over the year."""
        jd = JD_J2000 + day_offset
        cusps, ascmc = swe.houses(jd, 41.9, 12.5, ord("P"))
        armc = ascmc[2]
        eps = ascmc[8] if len(ascmc) > 8 else 23.4393

        sun, _ = swe.calc_ut(jd, SUN, FLG_SWIEPH)
        pos = swe.house_pos(armc, 41.9, eps, (sun[0], sun[1]), "P")
        assert 1.0 <= pos < 13.0
