"""Tests for rise/set at extreme latitudes and edge cases."""

from __future__ import annotations

import pytest
import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    VENUS,
    JUPITER,
    CALC_RISE,
    CALC_SET,
    CALC_MTRANSIT,
    CALC_ITRANSIT,
)

JD_J2000 = 2451545.0
# June 21 (summer solstice) and December 21 (winter solstice) around 2000
JD_SUMMER = 2451716.0  # ~2000-06-20
JD_WINTER = 2451900.0  # ~2000-12-21


@pytest.mark.unit
class TestRiseSetBasic:
    """Basic rise/set functionality."""

    @pytest.mark.parametrize(
        "lon,lat,name",
        [
            (12.5, 41.9, "Rome"),
            (-0.1, 51.5, "London"),
            (-74.0, 40.7, "New York"),
            (139.7, 35.7, "Tokyo"),
            (0.0, 0.0, "Equator"),
        ],
    )
    def test_sun_rise_found(self, lon, lat, name):
        """Sun rise should be found at temperate latitudes."""
        res, tret = swe.rise_trans(JD_J2000, SUN, CALC_RISE, (lon, lat, 0.0))
        assert res == 0, f"Rise not found at {name}"
        assert tret[0] > JD_J2000

    @pytest.mark.parametrize(
        "lon,lat,name",
        [
            (12.5, 41.9, "Rome"),
            (-0.1, 51.5, "London"),
            (-74.0, 40.7, "New York"),
        ],
    )
    def test_sun_set_found(self, lon, lat, name):
        """Sun set should be found at temperate latitudes."""
        res, tret = swe.rise_trans(JD_J2000, SUN, CALC_SET, (lon, lat, 0.0))
        assert res == 0, f"Set not found at {name}"
        assert tret[0] > JD_J2000

    def test_rise_before_set(self):
        """Sunrise should occur before sunset on the same day."""
        geopos = (12.5, 41.9, 0.0)
        res_r, tret_r = swe.rise_trans(JD_J2000, SUN, CALC_RISE, geopos)
        res_s, tret_s = swe.rise_trans(JD_J2000, SUN, CALC_SET, geopos)
        if res_r == 0 and res_s == 0:
            # Both rise and set on same day
            rise_jd = tret_r[0]
            set_jd = tret_s[0]
            # If rise is found first, set should be after
            if abs(rise_jd - set_jd) < 0.5:
                assert rise_jd < set_jd

    def test_transit_found(self):
        """Upper transit should be found."""
        geopos = (12.5, 41.9, 0.0)
        res, tret = swe.rise_trans(JD_J2000, SUN, CALC_MTRANSIT, geopos)
        assert res == 0
        assert tret[0] > JD_J2000

    def test_lower_transit_found(self):
        """Lower transit should be found."""
        geopos = (12.5, 41.9, 0.0)
        res, tret = swe.rise_trans(JD_J2000, SUN, CALC_ITRANSIT, geopos)
        assert res == 0
        assert tret[0] > JD_J2000

    def test_transit_upper_lower_12h_apart(self):
        """Upper and lower transit should be ~12 hours apart."""
        geopos = (0.0, 0.0, 0.0)
        res_u, tret_u = swe.rise_trans(JD_J2000, SUN, CALC_MTRANSIT, geopos)
        res_l, tret_l = swe.rise_trans(JD_J2000, SUN, CALC_ITRANSIT, geopos)
        if res_u == 0 and res_l == 0:
            diff = abs(tret_u[0] - tret_l[0])
            # Should be ~0.5 day (12h)
            assert 0.45 < diff < 0.55, f"Transit diff {diff * 24}h not ~12h"


@pytest.mark.unit
class TestRiseSetPolar:
    """Rise/set at polar latitudes (midnight sun / polar night)."""

    def test_midnight_sun_arctic_summer(self):
        """At 70°N in summer, Sun should not set (circumpolar)."""
        geopos = (0.0, 70.0, 0.0)
        res, tret = swe.rise_trans(JD_SUMMER, SUN, CALC_SET, geopos)
        # res=-2 means circumpolar (never sets)
        assert res == -2, f"Expected circumpolar at 70N summer, got res={res}"

    def test_polar_night_arctic_winter(self):
        """At 70°N in winter, Sun should not rise (polar night)."""
        geopos = (0.0, 70.0, 0.0)
        res, tret = swe.rise_trans(JD_WINTER, SUN, CALC_RISE, geopos)
        assert res == -2, f"Expected polar night at 70N winter, got res={res}"

    def test_equatorial_always_rises_and_sets(self):
        """At equator, Sun should always rise and set."""
        geopos = (0.0, 0.0, 0.0)
        for jd_offset in range(0, 365, 30):
            jd = JD_J2000 + jd_offset
            res_r, _ = swe.rise_trans(jd, SUN, CALC_RISE, geopos)
            res_s, _ = swe.rise_trans(jd, SUN, CALC_SET, geopos)
            assert res_r == 0, f"No sunrise at equator, JD offset {jd_offset}"
            assert res_s == 0, f"No sunset at equator, JD offset {jd_offset}"


@pytest.mark.unit
class TestRiseSetPlanets:
    """Rise/set for planets."""

    @pytest.mark.parametrize(
        "body,name",
        [
            (MOON, "Moon"),
            (VENUS, "Venus"),
            (JUPITER, "Jupiter"),
        ],
    )
    def test_planet_rise_rome(self, body, name):
        """Planets should rise at Rome."""
        geopos = (12.5, 41.9, 0.0)
        res, tret = swe.rise_trans(JD_J2000, body, CALC_RISE, geopos)
        assert res == 0, f"{name} rise not found at Rome"
        assert tret[0] > JD_J2000

    def test_moon_rise_daily_for_week(self):
        """Moon should rise each day over a week (delays ~50min/day on average)."""
        geopos = (0.0, 45.0, 0.0)
        rise_jds = []
        for day in range(7):
            jd = JD_J2000 + day
            res, tret = swe.rise_trans(jd, MOON, CALC_RISE, geopos)
            if res == 0:
                rise_jds.append(tret[0])
        # Should find at least 5 rises in 7 days
        assert len(rise_jds) >= 5
        # Average delay between consecutive rises should be ~24h50m
        if len(rise_jds) >= 2:
            total_span = (rise_jds[-1] - rise_jds[0]) * 24 * 60  # minutes
            avg_delay = total_span / (len(rise_jds) - 1)
            assert 1400 < avg_delay < 1600, f"Avg moon rise period {avg_delay} min"

    def test_equatorial_day_length_near_12h(self):
        """At equator, day length should be ~12 hours year-round."""
        geopos = (0.0, 0.0, 0.0)
        for offset in [0, 90, 180, 270]:
            jd = JD_J2000 + offset
            res_r, tret_r = swe.rise_trans(jd, SUN, CALC_RISE, geopos)
            res_s, tret_s = swe.rise_trans(jd, SUN, CALC_SET, geopos)
            if res_r == 0 and res_s == 0:
                rise = tret_r[0]
                sset = tret_s[0]
                if sset < rise:
                    sset += 1.0
                daylight = (sset - rise) * 24
                assert 11.5 < daylight < 12.5, (
                    f"Day length {daylight}h at equator, offset={offset}"
                )
