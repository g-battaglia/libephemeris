"""
Tests for edge cases: wrap-around at 0/360 with sidereal,
Earth heliocentric+sidereal, and other boundary conditions.
"""

from __future__ import annotations

import math

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
    EARTH,
    MEAN_NODE,
    TRUE_NODE,
    FLG_SWIEPH,
    FLG_SPEED,
    FLG_SIDEREAL,
    FLG_HELCTR,
    SIDM_LAHIRI,
    SIDM_FAGAN_BRADLEY,
    SIDM_RAMAN,
)


@pytest.fixture(autouse=True)
def _reset_state():
    yield
    swe.close()


JD_J2000 = 2451545.0


class TestSiderealWrapAround:
    """Test that sidereal positions wrap correctly at 0/360 boundary."""

    @pytest.mark.unit
    def test_sidereal_positions_in_range(self):
        """All sidereal positions should be in [0, 360)."""
        swe.set_sid_mode(SIDM_LAHIRI)
        flags = FLG_SWIEPH | FLG_SPEED | FLG_SIDEREAL
        bodies = [SUN, MOON, MERCURY, VENUS, MARS, JUPITER, SATURN]
        for body in bodies:
            pos, _ = swe.calc_ut(JD_J2000, body, flags)
            assert 0.0 <= pos[0] < 360.0, (
                f"Body {body} sidereal lon {pos[0]} out of range"
            )

    @pytest.mark.unit
    def test_sidereal_wrap_across_dates(self):
        """Sidereal positions stay in [0, 360) across many dates."""
        swe.set_sid_mode(SIDM_LAHIRI)
        flags = FLG_SWIEPH | FLG_SPEED | FLG_SIDEREAL
        for offset in range(0, 365, 10):
            jd = JD_J2000 + offset
            pos, _ = swe.calc_ut(jd, SUN, flags)
            assert 0.0 <= pos[0] < 360.0, f"Sun sidereal lon at JD {jd}: {pos[0]}"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "sid_mode",
        [
            SIDM_LAHIRI,
            SIDM_FAGAN_BRADLEY,
            SIDM_RAMAN,
        ],
    )
    def test_sidereal_modes_all_in_range(self, sid_mode):
        """All sidereal modes produce positions in [0, 360)."""
        swe.set_sid_mode(sid_mode)
        flags = FLG_SWIEPH | FLG_SPEED | FLG_SIDEREAL
        for body in [SUN, MOON, MARS]:
            pos, _ = swe.calc_ut(JD_J2000, body, flags)
            assert 0.0 <= pos[0] < 360.0

    @pytest.mark.unit
    def test_sidereal_near_zero_crossing(self):
        """Find a date where Sun sidereal lon is near 0/360 and verify no glitch."""
        swe.set_sid_mode(SIDM_LAHIRI)
        flags = FLG_SWIEPH | FLG_SPEED | FLG_SIDEREAL
        # Scan for 0-crossing (Sun enters sidereal Aries ~April 14)
        # Start from Jan 1 2000 and scan
        prev_lon = None
        for day in range(0, 366):
            jd = JD_J2000 + day
            pos, _ = swe.calc_ut(jd, SUN, flags)
            lon = pos[0]
            assert 0.0 <= lon < 360.0
            if prev_lon is not None and prev_lon > 350.0 and lon < 10.0:
                # Crossing detected — speed should still be positive
                assert pos[3] > 0.0, f"Sun speed negative at 0-crossing: {pos[3]}"
            prev_lon = lon


class TestEarthHeliocentric:
    """Test Earth with heliocentric flag."""

    @pytest.mark.unit
    def test_earth_heliocentric_valid(self):
        """Earth heliocentric returns a valid position (may be non-zero).

        Unlike geocentric Earth (which is always zeros), heliocentric Earth
        returns the Earth's heliocentric ecliptic position (opposite of Sun).
        """
        flags = FLG_SWIEPH | FLG_SPEED | FLG_HELCTR
        pos, _ = swe.calc_ut(JD_J2000, EARTH, flags)
        for i in range(6):
            assert math.isfinite(pos[i]), f"Earth helio pos[{i}] not finite"
        # Heliocentric Earth longitude should be ~180° from Sun
        assert 0.0 <= pos[0] < 360.0

    @pytest.mark.unit
    def test_earth_heliocentric_sidereal_valid(self):
        """Earth heliocentric + sidereal returns valid values."""
        swe.set_sid_mode(SIDM_LAHIRI)
        flags = FLG_SWIEPH | FLG_SPEED | FLG_HELCTR | FLG_SIDEREAL
        pos, _ = swe.calc_ut(JD_J2000, EARTH, flags)
        for i in range(6):
            assert math.isfinite(pos[i])


class TestHeliocentricSidereal:
    """Test heliocentric + sidereal combination."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body", [MERCURY, VENUS, MARS, JUPITER])
    def test_helio_sidereal_in_range(self, body):
        """Heliocentric sidereal positions are in [0, 360)."""
        swe.set_sid_mode(SIDM_LAHIRI)
        flags = FLG_SWIEPH | FLG_SPEED | FLG_HELCTR | FLG_SIDEREAL
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        assert 0.0 <= pos[0] < 360.0

    @pytest.mark.unit
    def test_helio_sidereal_differs_from_tropical(self):
        """Heliocentric sidereal differs from heliocentric tropical by ayanamsa."""
        swe.set_sid_mode(SIDM_LAHIRI)
        ayan = swe.get_ayanamsa_ut(JD_J2000)

        flags_trop = FLG_SWIEPH | FLG_SPEED | FLG_HELCTR
        flags_sid = flags_trop | FLG_SIDEREAL

        pos_t, _ = swe.calc_ut(JD_J2000, MARS, flags_trop)
        pos_s, _ = swe.calc_ut(JD_J2000, MARS, flags_sid)

        diff = (pos_t[0] - pos_s[0]) % 360.0
        if diff > 180.0:
            diff = 360.0 - diff
        assert diff == pytest.approx(ayan, abs=0.5), (
            f"Difference {diff} doesn't match ayanamsa {ayan}"
        )


class TestNodeWrapAround:
    """Test Mean/True Node near 0/360 boundary."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body", [MEAN_NODE, TRUE_NODE])
    def test_node_in_range(self, body):
        """Node longitude always in [0, 360)."""
        flags = FLG_SWIEPH | FLG_SPEED
        for day in range(0, 365, 30):
            jd = JD_J2000 + day
            pos, _ = swe.calc_ut(jd, body, flags)
            assert 0.0 <= pos[0] < 360.0

    @pytest.mark.unit
    def test_mean_node_retrograde(self):
        """Mean node moves retrograde — later date has smaller longitude."""
        flags = FLG_SWIEPH | FLG_SPEED
        pos1, _ = swe.calc_ut(JD_J2000, MEAN_NODE, flags)
        pos2, _ = swe.calc_ut(JD_J2000 + 30, MEAN_NODE, flags)
        # Mean node moves ~0.053 deg/day retrograde
        # After 30 days, it should have moved ~1.6 degrees backward
        # But must account for 0/360 wrap
        diff = pos1[0] - pos2[0]
        if diff < -180:
            diff += 360
        assert diff > 0, "Mean node should have moved retrograde"


class TestLatitudeBoundary:
    """Test latitude near extreme values."""

    @pytest.mark.unit
    def test_moon_latitude_range(self):
        """Moon latitude should be within [-5.3, +5.3] degrees."""
        flags = FLG_SWIEPH | FLG_SPEED
        for day in range(0, 30):
            jd = JD_J2000 + day
            pos, _ = swe.calc_ut(jd, MOON, flags)
            assert -6.0 < pos[1] < 6.0, f"Moon lat {pos[1]} at JD {jd}"

    @pytest.mark.unit
    @pytest.mark.parametrize("body", [MERCURY, VENUS, MARS])
    def test_planet_latitude_reasonable(self, body):
        """Planet latitudes should be within reasonable bounds."""
        flags = FLG_SWIEPH | FLG_SPEED
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        # Planet ecliptic latitudes are small (< 10° typically)
        assert -15.0 < pos[1] < 15.0
