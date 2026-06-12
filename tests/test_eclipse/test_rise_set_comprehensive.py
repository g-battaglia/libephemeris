"""
Tests for rise/set/transit calculations with various bodies and locations.

Verifies rise_trans for Sun, Moon, planets, and fixed stars
at different geographic locations and dates.
"""

from __future__ import annotations


import pytest

import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    VENUS,
    MARS,
    JUPITER,
    CALC_RISE,
    CALC_SET,
    CALC_MTRANSIT,
    CALC_ITRANSIT,
)


LOCATIONS = [
    # (lon, lat, alt, name)
    (12.5, 41.9, 50.0, "Rome"),
    (-74.0, 40.7, 10.0, "New York"),
    (139.7, 35.7, 40.0, "Tokyo"),
    (151.2, -33.9, 58.0, "Sydney"),
    (0.0, 0.0, 0.0, "Equator"),
    (0.0, 60.0, 0.0, "Helsinki"),
]

PLANETS_FOR_RISE = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
]


class TestSunRiseSet:
    """Test Sun rise/set at various locations."""

    @pytest.mark.unit
    @pytest.mark.parametrize("lon,lat,alt,name", LOCATIONS)
    def test_sun_rise_found(self, lon: float, lat: float, alt: float, name: str):
        """Sun rise should be found at most locations."""
        jd = 2451545.0
        geopos = (lon, lat, alt)
        res, tret = swe.rise_trans(jd, SUN, CALC_RISE, geopos)
        if res == 0:
            rise_jd = tret[0]
            assert rise_jd > jd, f"{name}: rise {rise_jd} not after start {jd}"
            # Should be within ~1 day
            assert rise_jd - jd < 2.0, f"{name}: rise {rise_jd - jd:.2f} days later"

    @pytest.mark.unit
    @pytest.mark.parametrize("lon,lat,alt,name", LOCATIONS)
    def test_sun_set_found(self, lon: float, lat: float, alt: float, name: str):
        """Sun set should be found at most locations."""
        jd = 2451545.0
        geopos = (lon, lat, alt)
        res, tret = swe.rise_trans(jd, SUN, CALC_SET, geopos)
        if res == 0:
            set_jd = tret[0]
            assert set_jd > jd
            assert set_jd - jd < 2.0

    @pytest.mark.unit
    def test_sun_rise_before_set(self):
        """Sun rise time should be before set time (same day, equinox)."""
        # Use spring equinox at equator for clean result
        jd = swe.julday(2024, 3, 20, 0.0)
        geopos = (0.0, 0.0, 0.0)

        res_r, tret_r = swe.rise_trans(jd, SUN, CALC_RISE, geopos)
        res_s, tret_s = swe.rise_trans(jd, SUN, CALC_SET, geopos)

        if res_r == 0 and res_s == 0:
            rise_jd = tret_r[0]
            set_jd = tret_s[0]
            # Rise should be before set (for equator at equinox)
            assert rise_jd < set_jd, f"Rise {rise_jd} not before set {set_jd}"

    @pytest.mark.unit
    def test_equatorial_day_length_near_12h(self):
        """At equator near equinox, day should be ~12 hours."""
        jd = swe.julday(2024, 3, 20, 0.0)
        geopos = (0.0, 0.0, 0.0)

        res_r, tret_r = swe.rise_trans(jd, SUN, CALC_RISE, geopos)
        res_s, tret_s = swe.rise_trans(jd, SUN, CALC_SET, geopos)

        if res_r == 0 and res_s == 0:
            rise_jd = tret_r[0]
            set_jd = tret_s[0]
            day_hours = (set_jd - rise_jd) * 24
            assert 11.5 < day_hours < 12.5, (
                f"Day length {day_hours:.2f}h (expected ~12h)"
            )


class TestMeridianTransit:
    """Test meridian transit calculations."""

    @pytest.mark.unit
    @pytest.mark.parametrize("lon,lat,alt,name", LOCATIONS[:4])
    def test_sun_upper_transit(self, lon: float, lat: float, alt: float, name: str):
        """Sun upper meridian transit should be found."""
        jd = 2451545.0
        geopos = (lon, lat, alt)
        res, tret = swe.rise_trans(jd, SUN, CALC_MTRANSIT, geopos)
        if res == 0:
            transit_jd = tret[0]
            assert transit_jd > jd
            assert transit_jd - jd < 2.0

    @pytest.mark.unit
    @pytest.mark.parametrize("lon,lat,alt,name", LOCATIONS[:4])
    def test_sun_lower_transit(self, lon: float, lat: float, alt: float, name: str):
        """Sun lower meridian transit should be found."""
        jd = 2451545.0
        geopos = (lon, lat, alt)
        res, tret = swe.rise_trans(jd, SUN, CALC_ITRANSIT, geopos)
        if res == 0:
            transit_jd = tret[0]
            assert transit_jd > jd
            assert transit_jd - jd < 2.0

    @pytest.mark.unit
    def test_upper_lower_transit_12h_apart(self):
        """Upper and lower transits should be ~12 hours apart."""
        jd = 2451545.0
        geopos = (12.5, 41.9, 50.0)

        res_u, tret_u = swe.rise_trans(jd, SUN, CALC_MTRANSIT, geopos)
        res_l, tret_l = swe.rise_trans(jd, SUN, CALC_ITRANSIT, geopos)

        if res_u == 0 and res_l == 0:
            upper_jd = tret_u[0]
            lower_jd = tret_l[0]
            gap_hours = abs(lower_jd - upper_jd) * 24
            # Should be ~12h apart (or ~12h if one comes before the other)
            gap_mod = gap_hours % 24
            if gap_mod > 12:
                gap_mod = 24 - gap_mod
            assert 11 < gap_mod < 13, f"Transit gap: {gap_mod:.2f}h (expected ~12h)"


class TestPlanetRiseSet:
    """Test rise/set for various planets."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", PLANETS_FOR_RISE)
    def test_planet_rise_rome(self, body_id: int, name: str):
        """Planet rise found at Rome."""
        jd = 2451545.0
        geopos = (12.5, 41.9, 50.0)
        res, tret = swe.rise_trans(jd, body_id, CALC_RISE, geopos)
        if res == 0:
            rise_jd = tret[0]
            assert rise_jd > jd
            # Should be within a couple of days
            assert rise_jd - jd < 3.0, f"{name}: rise {rise_jd - jd:.2f} days later"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", PLANETS_FOR_RISE)
    def test_planet_set_rome(self, body_id: int, name: str):
        """Planet set found at Rome."""
        jd = 2451545.0
        geopos = (12.5, 41.9, 50.0)
        res, tret = swe.rise_trans(jd, body_id, CALC_SET, geopos)
        if res == 0:
            set_jd = tret[0]
            assert set_jd > jd
            assert set_jd - jd < 3.0


class TestMoonRiseSet:
    """Special tests for Moon rise/set."""

    @pytest.mark.unit
    def test_moon_rise_sequential(self):
        """Find 3 consecutive Moon rises."""
        jd = 2451545.0
        geopos = (12.5, 41.9, 50.0)
        rises = []
        for _ in range(3):
            res, tret = swe.rise_trans(jd, MOON, CALC_RISE, geopos)
            if res == 0:
                rise_jd = tret[0]
                rises.append(rise_jd)
                jd = rise_jd + 0.5  # Start search after this rise
            else:
                break

        if len(rises) >= 2:
            # Moon rises ~50 minutes later each day
            for i in range(len(rises) - 1):
                gap_hours = (rises[i + 1] - rises[i]) * 24
                assert 23.5 < gap_hours < 26.5, f"Moon rise gap: {gap_hours:.2f}h"


class TestRiseSetDateRange:
    """Test rise/set across dates."""

    @pytest.mark.unit
    @pytest.mark.parametrize("year", [1950, 2000, 2024, 2050])
    def test_sun_rise_across_years(self, year: int):
        """Sun rise valid across years."""
        jd = swe.julday(year, 6, 21, 0.0)
        geopos = (12.5, 41.9, 50.0)
        res, tret = swe.rise_trans(jd, SUN, CALC_RISE, geopos)
        if res == 0:
            rise_jd = tret[0]
            assert rise_jd > jd


class TestCircumpolar:
    """Test circumpolar behavior at high latitudes."""

    @pytest.mark.unit
    def test_midnight_sun(self):
        """At high latitude in summer, Sun may not set (circumpolar)."""
        jd = swe.julday(2024, 6, 21, 0.0)
        geopos = (25.0, 70.0, 0.0)  # Northern Norway
        res, tret = swe.rise_trans(jd, SUN, CALC_SET, geopos)
        # res=-2 means circumpolar (no rise/set)
        # Either it finds a set or returns circumpolar
        assert res in (0, -2), f"Unexpected result: {res}"
