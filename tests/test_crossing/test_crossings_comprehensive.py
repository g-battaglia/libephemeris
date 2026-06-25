"""Tests for crossing functions and station finding with more edge cases."""

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
    SATURN,
    FLG_SWIEPH,
    FLG_SPEED,
)

JD_J2000 = 2451545.0


@pytest.mark.unit
class TestSolcrossAllSigns:
    """Test Sun crossing all zodiac sign boundaries."""

    @pytest.mark.parametrize("degree", list(range(0, 360, 30)))
    def test_solcross_sign_boundary(self, degree):
        """Sun should cross each sign boundary once per year."""
        jd = swe.solcross_ut(float(degree), JD_J2000, FLG_SWIEPH)
        assert jd > JD_J2000
        # Should be within 1 year
        assert jd < JD_J2000 + 370, f"Cross at {degree}° too far: {jd - JD_J2000} days"
        # Verify Sun is actually at the crossing degree
        sun, _ = swe.calc_ut(jd, SUN, FLG_SWIEPH)
        diff = abs(sun[0] - degree)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.001, f"Sun at {sun[0]}° but expected {degree}°"

    def test_solcross_equinox_date(self):
        """Vernal equinox (0° Aries) should occur around March 20."""
        jd = swe.solcross_ut(0.0, JD_J2000, FLG_SWIEPH)
        y, m, d, h = swe.revjul(jd)
        assert m == 3, f"Equinox month: {m}"
        assert 19 <= d <= 21, f"Equinox day: {d}"

    def test_solcross_solstices(self):
        """Summer and winter solstices should be around June/December 21."""
        jd_summer = swe.solcross_ut(90.0, JD_J2000, FLG_SWIEPH)
        y, m, d, h = swe.revjul(jd_summer)
        assert m == 6, f"Summer solstice month: {m}"
        assert 20 <= d <= 22

        jd_winter = swe.solcross_ut(270.0, JD_J2000, FLG_SWIEPH)
        y, m, d, h = swe.revjul(jd_winter)
        assert m == 12, f"Winter solstice month: {m}"
        assert 20 <= d <= 22


@pytest.mark.unit
class TestMooncrossExtended:
    """Extended Moon crossing tests."""

    def test_mooncross_full_zodiac(self):
        """Moon should cross all 12 sign boundaries within 1 month."""
        jd = JD_J2000
        for degree in range(0, 360, 30):
            jd_cross = swe.mooncross_ut(float(degree), JD_J2000, FLG_SWIEPH)
            assert jd_cross > JD_J2000
            # Within 1 month
            assert jd_cross < JD_J2000 + 30, f"Moon cross {degree}° too far"
            # Verify
            moon, _ = swe.calc_ut(jd_cross, MOON, FLG_SWIEPH)
            diff = abs(moon[0] - degree)
            if diff > 180:
                diff = 360 - diff
            assert diff < 0.01, f"Moon at {moon[0]}° expected {degree}°"

    def test_mooncross_consecutive_same_degree(self):
        """Two consecutive Moon crossings of the same degree should be ~27.3 days apart."""
        degree = 0.0
        jd1 = swe.mooncross_ut(degree, JD_J2000, FLG_SWIEPH)
        jd2 = swe.mooncross_ut(degree, jd1 + 1.0, FLG_SWIEPH)
        gap = jd2 - jd1
        # Sidereal month ~27.3 days
        assert 26 < gap < 29, f"Moon return gap: {gap} days"

    def test_mooncross_node_returns_tuple(self):
        """mooncross_node_ut returns (jd, lon, lat) tuple."""
        result = swe.mooncross_node_ut(JD_J2000, FLG_SWIEPH)
        assert len(result) == 3
        jd, lon, lat = result
        assert jd > JD_J2000
        assert abs(lat) < 0.01, "Moon lat should be ~0 at node crossing"

    def test_mooncross_node_sequential(self):
        """Node crossings should be ~13.6 days apart (half draconic month)."""
        result1 = swe.mooncross_node_ut(JD_J2000, FLG_SWIEPH)
        result2 = swe.mooncross_node_ut(result1[0] + 1.0, FLG_SWIEPH)
        gap = result2[0] - result1[0]
        # Half draconic month ~13.6 days
        assert 12 < gap < 15, f"Node crossing gap: {gap} days"


@pytest.mark.unit
class TestStationsExtended:
    """Extended tests for planetary station finding."""

    @pytest.mark.parametrize(
        "body,name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_station_found(self, body, name):
        """Should find a station for each planet."""
        jd, stype = swe.find_station_ut(body, JD_J2000, "any", FLG_SWIEPH)
        assert jd > JD_J2000
        assert stype in ("SR", "SD"), f"{name} station type: {stype}"

    @pytest.mark.parametrize(
        "body,name,max_wait",
        [
            (MERCURY, "Mercury", 120),
            (MARS, "Mars", 800),
            (JUPITER, "Jupiter", 400),
            (SATURN, "Saturn", 400),
        ],
    )
    def test_station_within_reasonable_time(self, body, name, max_wait):
        """Station should be found within a reasonable time."""
        jd, stype = swe.find_station_ut(body, JD_J2000, "any", FLG_SWIEPH)
        wait = jd - JD_J2000
        assert wait < max_wait, f"{name} station in {wait} days (max {max_wait})"

    def test_mercury_station_speed_near_zero(self):
        """At station, planet speed should be near zero."""
        jd, stype = swe.find_station_ut(MERCURY, JD_J2000, "any", FLG_SWIEPH)
        result, _ = swe.calc_ut(jd, MERCURY, FLG_SWIEPH | FLG_SPEED)
        assert abs(result[3]) < 0.01, f"Mercury speed at station: {result[3]}°/day"

    def test_retrograde_station_pair(self):
        """SR followed by SD should define a retrograde period."""
        jd_sr, stype_sr = swe.find_station_ut(MERCURY, JD_J2000, "SR", FLG_SWIEPH)
        assert stype_sr == "SR"
        jd_sd, stype_sd = swe.find_station_ut(
            MERCURY, jd_sr + 1.0, "SD", FLG_SWIEPH
        )
        assert stype_sd == "SD"
        assert jd_sd > jd_sr
        # Mercury retrogrades last ~21 days
        retro_days = jd_sd - jd_sr
        assert 15 < retro_days < 30, f"Mercury retrograde: {retro_days} days"

    def test_sr_sd_alternate(self):
        """Stations should alternate SR/SD."""
        jd = JD_J2000
        prev_type = None
        for _ in range(4):
            jd_s, stype = swe.find_station_ut(MERCURY, jd, "any", FLG_SWIEPH)
            if prev_type is not None:
                assert stype != prev_type, f"Two {stype} in a row"
            prev_type = stype
            jd = jd_s + 1.0
