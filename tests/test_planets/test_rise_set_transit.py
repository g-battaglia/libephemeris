"""
Tests for rise/set/transit calculations.

Verifies that rise_trans returns valid results for various
bodies, locations, and edge cases including polar regions.
"""

from __future__ import annotations


import pytest

import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    JUPITER,
    VENUS,
    CALC_RISE,
    CALC_SET,
    CALC_MTRANSIT,
    CALC_ITRANSIT,
)


STANDARD_BODIES = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (VENUS, "Venus"),
]

STANDARD_LOCATIONS = [
    (41.9, 12.5, 0, "Rome"),
    (40.7, -74.0, 0, "New York"),
    (35.7, 139.7, 0, "Tokyo"),
    (-33.9, 151.2, 0, "Sydney"),
    (51.5, -0.1, 0, "London"),
    (0.0, 0.0, 0, "Equator"),
]


class TestRiseTransBasic:
    """Basic rise/set/transit tests."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", STANDARD_BODIES)
    def test_sun_rise_returns_valid(self, body_id: int, name: str):
        """Rise calculation returns valid result."""
        jd = 2451545.0
        geopos = (41.9, 12.5, 0)
        res, tret = swe.rise_trans(jd, body_id, CALC_RISE, geopos)
        # res=0 means found, res=-2 means circumpolar
        assert res in (0, -2), f"{name}: unexpected result code {res}"
        if res == 0:
            rise_jd = tret[0]
            assert rise_jd > jd - 1, f"{name}: rise time {rise_jd} too early"
            assert rise_jd < jd + 2, f"{name}: rise time {rise_jd} too late"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", STANDARD_BODIES)
    def test_sun_set_returns_valid(self, body_id: int, name: str):
        """Set calculation returns valid result."""
        jd = 2451545.0
        geopos = (41.9, 12.5, 0)
        res, tret = swe.rise_trans(jd, body_id, CALC_SET, geopos)
        assert res in (0, -2), f"{name}: unexpected result code {res}"
        if res == 0:
            set_jd = tret[0]
            assert set_jd > jd - 1
            assert set_jd < jd + 2

    @pytest.mark.unit
    def test_sun_rise_before_set(self):
        """Sun rise should occur before sun set on same day."""
        jd = 2451545.0  # Noon
        geopos = (41.9, 12.5, 0)
        res_r, tret_r = swe.rise_trans(jd, SUN, CALC_RISE, geopos)
        res_s, tret_s = swe.rise_trans(jd, SUN, CALC_SET, geopos)
        if res_r == 0 and res_s == 0:
            rise_jd = tret_r[0]
            set_jd = tret_s[0]
            # Both should be on the same day or close
            assert abs(rise_jd - set_jd) < 1.0

    @pytest.mark.unit
    def test_transit_returns_valid(self):
        """Upper transit returns valid time."""
        jd = 2451545.0
        geopos = (41.9, 12.5, 0)
        res, tret = swe.rise_trans(jd, SUN, CALC_MTRANSIT, geopos)
        assert res == 0, f"Transit not found, result={res}"
        transit_jd = tret[0]
        assert transit_jd > jd - 1
        assert transit_jd < jd + 2

    @pytest.mark.unit
    def test_lower_transit_returns_valid(self):
        """Lower transit returns valid time."""
        jd = 2451545.0
        geopos = (41.9, 12.5, 0)
        res, tret = swe.rise_trans(jd, SUN, CALC_ITRANSIT, geopos)
        assert res == 0, f"Lower transit not found, result={res}"


class TestRiseTransLocations:
    """Test rise/set at various geographic locations."""

    @pytest.mark.unit
    @pytest.mark.parametrize("lat,lon,alt,name", STANDARD_LOCATIONS)
    def test_sun_rise_at_location(self, lat: float, lon: float, alt: int, name: str):
        """Sun rise is found at standard locations."""
        jd = 2451545.0
        # geopos is (longitude, latitude, altitude) per the public contract:
        # the swapped order passed Tokyo/Sydney longitudes as latitudes >90
        # and now trips the geodetic-domain validation.
        geopos = (lon, lat, alt)
        res, tret = swe.rise_trans(jd, SUN, CALC_RISE, geopos)
        assert res in (0, -2), f"{name}: result={res}"

    @pytest.mark.unit
    @pytest.mark.parametrize("lat,lon,alt,name", STANDARD_LOCATIONS)
    def test_moon_rise_at_location(self, lat: float, lon: float, alt: int, name: str):
        """Moon rise is found at standard locations."""
        jd = 2451545.0
        geopos = (lon, lat, alt)
        res, tret = swe.rise_trans(jd, MOON, CALC_RISE, geopos)
        assert res in (0, -2), f"{name}: result={res}"

    @pytest.mark.unit
    @pytest.mark.parametrize("lat", [30, 45, 55, 60])
    def test_sun_daylight_duration(self, lat: int):
        """Sun is above horizon for reasonable duration at J2000."""
        jd = 2451545.0  # Jan 1
        geopos = (float(lat), 0.0, 0)
        res_r, tret_r = swe.rise_trans(jd, SUN, CALC_RISE, geopos)
        res_s, tret_s = swe.rise_trans(jd, SUN, CALC_SET, geopos)
        if res_r == 0 and res_s == 0:
            daylight_hours = (tret_s[0] - tret_r[0]) * 24
            if daylight_hours < 0:
                daylight_hours += 24
            # In January at 30-60N, daylight should be 7-11 hours
            assert 5 < daylight_hours < 18, f"Lat {lat}: daylight {daylight_hours:.1f}h"


class TestRiseTransPolar:
    """Test rise/set at polar latitudes."""

    @pytest.mark.unit
    @pytest.mark.parametrize("lat", [70.0, 75.0, 80.0, 85.0, 89.0])
    def test_sun_rise_polar_summer(self, lat: float):
        """Sun may be circumpolar in polar summer."""
        # June solstice
        jd = swe.julday(2000, 6, 21, 12.0)
        geopos = (lat, 0.0, 0)
        res, tret = swe.rise_trans(jd, SUN, CALC_RISE, geopos)
        # At extreme latitudes in summer, res=-2 (circumpolar)
        assert res in (0, -2), f"Lat {lat}: unexpected result {res}"

    @pytest.mark.unit
    @pytest.mark.parametrize("lat", [70.0, 75.0, 80.0, 85.0, 89.0])
    def test_sun_rise_polar_winter(self, lat: float):
        """Sun may not rise in polar winter."""
        # December solstice
        jd = swe.julday(2000, 12, 21, 12.0)
        geopos = (lat, 0.0, 0)
        res, tret = swe.rise_trans(jd, SUN, CALC_RISE, geopos)
        assert res in (0, -2), f"Lat {lat}: unexpected result {res}"


class TestRiseTransConsecutiveDays:
    """Test rise/set across consecutive days."""

    @pytest.mark.unit
    def test_sun_rise_daily_for_week(self):
        """Sun rise times should advance ~4 min/day at mid-latitudes."""
        jd_start = 2451545.0
        geopos = (41.9, 12.5, 0)
        rise_times = []
        for i in range(7):
            jd = jd_start + i
            res, tret = swe.rise_trans(jd, SUN, CALC_RISE, geopos)
            if res == 0:
                rise_times.append(tret[0])

        # All rise times should be found (Rome, January)
        assert len(rise_times) == 7, "Not all rise times found"
        # Rise times should be roughly 1 day apart
        for i in range(6):
            gap = rise_times[i + 1] - rise_times[i]
            assert 0.95 < gap < 1.05, f"Day {i}: rise gap {gap:.4f} days"

    @pytest.mark.unit
    def test_moon_rise_daily_for_week(self):
        """Moon rise times should advance ~50 min/day."""
        jd_start = 2451545.0
        geopos = (41.9, 12.5, 0)
        rise_times = []
        for i in range(7):
            jd = jd_start + i
            res, tret = swe.rise_trans(jd, MOON, CALC_RISE, geopos)
            if res == 0:
                rise_times.append(tret[0])

        if len(rise_times) >= 2:
            for i in range(len(rise_times) - 1):
                gap = rise_times[i + 1] - rise_times[i]
                # Moon rises ~50 min later each day
                assert 0.95 < gap < 1.15, f"Day {i}: Moon rise gap {gap:.4f} days"
