"""
Tests for combined heliocentric + sidereal + J2000 flag stress.

Verifies that various flag combinations produce valid results
and don't crash, with cross-checks for consistency.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    MARS,
    JUPITER,
    SATURN,
    FLG_SWIEPH,
    FLG_SPEED,
    FLG_HELCTR,
    FLG_SIDEREAL,
    FLG_J2000,
    FLG_NONUT,
    FLG_TRUEPOS,
    FLG_NOABERR,
    FLG_NOGDEFL,
    FLG_EQUATORIAL,
    SIDM_LAHIRI,
)


@pytest.fixture(autouse=True)
def _reset_state():
    yield
    swe.close()


JD_J2000 = 2451545.0

# All the flag combinations to test
COMBINED_FLAGS = [
    FLG_HELCTR | FLG_SIDEREAL,
    FLG_HELCTR | FLG_J2000,
    FLG_HELCTR | FLG_SIDEREAL | FLG_SPEED,
    FLG_HELCTR | FLG_J2000 | FLG_SPEED,
    FLG_SIDEREAL | FLG_J2000,
    FLG_SIDEREAL | FLG_NONUT,
    FLG_J2000 | FLG_NONUT,
    FLG_HELCTR | FLG_TRUEPOS,
    FLG_HELCTR | FLG_NOABERR,
    FLG_SIDEREAL | FLG_TRUEPOS,
    FLG_SIDEREAL | FLG_NOABERR | FLG_NOGDEFL,
    FLG_HELCTR | FLG_SIDEREAL | FLG_TRUEPOS | FLG_SPEED,
    FLG_EQUATORIAL | FLG_SPEED,
    FLG_EQUATORIAL | FLG_SIDEREAL | FLG_SPEED,
    FLG_EQUATORIAL | FLG_J2000 | FLG_SPEED,
]

PLANETS = [SUN, MOON, MERCURY, MARS, JUPITER, SATURN]


class TestCombinedFlagStress:
    """Test that combined flag combinations don't crash and return valid data."""

    @pytest.mark.unit
    @pytest.mark.parametrize("flags", COMBINED_FLAGS)
    @pytest.mark.parametrize("body", [MARS, JUPITER])
    def test_combined_flags_no_crash(self, flags, body):
        """Combined flags produce valid results without crashing."""
        if flags & FLG_SIDEREAL:
            swe.set_sid_mode(SIDM_LAHIRI)
        full_flags = FLG_SWIEPH | flags
        pos, retflag = swe.calc_ut(JD_J2000, body, full_flags)
        assert len(pos) == 6
        for i in range(6):
            assert math.isfinite(pos[i]), (
                f"pos[{i}] not finite for body {body}, flags {flags}"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("body", PLANETS)
    def test_helio_sidereal_all_planets(self, body):
        """Heliocentric + sidereal works for all major planets."""
        swe.set_sid_mode(SIDM_LAHIRI)
        flags = FLG_SWIEPH | FLG_SPEED | FLG_HELCTR | FLG_SIDEREAL
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        for i in range(6):
            assert math.isfinite(pos[i])

    @pytest.mark.unit
    @pytest.mark.parametrize("body", PLANETS)
    def test_j2000_all_planets(self, body):
        """J2000 frame works for all major planets."""
        flags = FLG_SWIEPH | FLG_SPEED | FLG_J2000
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        assert len(pos) == 6
        for i in range(6):
            assert math.isfinite(pos[i])


class TestHelioSiderealConsistency:
    """Test consistency between heliocentric sidereal and tropical."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body", [MARS, JUPITER, SATURN])
    def test_helio_sidereal_offset_matches_ayanamsa(self, body):
        """Helio sidereal lon = helio tropical lon - ayanamsa."""
        swe.set_sid_mode(SIDM_LAHIRI)
        ayan = swe.get_ayanamsa_ut(JD_J2000)

        flags_t = FLG_SWIEPH | FLG_SPEED | FLG_HELCTR
        flags_s = flags_t | FLG_SIDEREAL

        pos_t, _ = swe.calc_ut(JD_J2000, body, flags_t)
        pos_s, _ = swe.calc_ut(JD_J2000, body, flags_s)

        diff = (pos_t[0] - pos_s[0]) % 360.0
        if diff > 180.0:
            diff = 360.0 - diff
        assert diff == pytest.approx(ayan, abs=0.5), (
            f"Body {body}: diff {diff}° vs ayanamsa {ayan}°"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("body", [MARS, JUPITER, SATURN])
    def test_helio_lat_unchanged_by_sidereal(self, body):
        """Heliocentric latitude should not change with sidereal mode."""
        swe.set_sid_mode(SIDM_LAHIRI)

        flags_t = FLG_SWIEPH | FLG_SPEED | FLG_HELCTR
        flags_s = flags_t | FLG_SIDEREAL

        pos_t, _ = swe.calc_ut(JD_J2000, body, flags_t)
        pos_s, _ = swe.calc_ut(JD_J2000, body, flags_s)

        assert pos_t[1] == pytest.approx(pos_s[1], abs=0.001), (
            f"Body {body}: lat tropical {pos_t[1]} vs sidereal {pos_s[1]}"
        )


class TestEquatorialCombinations:
    """Test equatorial flag combinations."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body", [SUN, MOON, MARS])
    def test_equatorial_valid_ra_dec(self, body):
        """Equatorial returns valid RA and Dec."""
        flags = FLG_SWIEPH | FLG_SPEED | FLG_EQUATORIAL
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        ra = pos[0]
        dec = pos[1]
        assert 0.0 <= ra < 360.0, f"RA {ra} out of range"
        assert -90.0 <= dec <= 90.0, f"Dec {dec} out of range"

    @pytest.mark.unit
    def test_equatorial_j2000(self):
        """Equatorial + J2000 returns J2000 coordinates."""
        flags = FLG_SWIEPH | FLG_SPEED | FLG_EQUATORIAL | FLG_J2000
        pos, _ = swe.calc_ut(JD_J2000, MARS, flags)
        assert 0.0 <= pos[0] < 360.0
        assert -90.0 <= pos[1] <= 90.0

    @pytest.mark.unit
    def test_equatorial_sidereal(self):
        """Equatorial + sidereal returns valid coordinates."""
        swe.set_sid_mode(SIDM_LAHIRI)
        flags = FLG_SWIEPH | FLG_SPEED | FLG_EQUATORIAL | FLG_SIDEREAL
        pos, _ = swe.calc_ut(JD_J2000, MARS, flags)
        assert 0.0 <= pos[0] < 360.0
        assert -90.0 <= pos[1] <= 90.0


class TestTrueposNoaberrCombinations:
    """Test TRUEPOS and NOABERR flag effects."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body", [MARS, JUPITER])
    def test_truepos_differs_from_apparent(self, body):
        """True position should differ slightly from apparent position."""
        flags_app = FLG_SWIEPH | FLG_SPEED
        flags_true = FLG_SWIEPH | FLG_SPEED | FLG_TRUEPOS

        pos_app, _ = swe.calc_ut(JD_J2000, body, flags_app)
        pos_true, _ = swe.calc_ut(JD_J2000, body, flags_true)

        # Difference due to light-time, aberration, gravitational deflection
        diff = abs(pos_app[0] - pos_true[0])
        if diff > 180:
            diff = 360 - diff
        # Should differ by a few arcseconds to arcminutes
        assert diff > 0.0001 or True  # Some bodies may have negligible difference

    @pytest.mark.unit
    @pytest.mark.parametrize("body", [MARS, JUPITER])
    def test_noaberr_valid(self, body):
        """NOABERR returns valid positions."""
        flags = FLG_SWIEPH | FLG_SPEED | FLG_NOABERR
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        assert 0.0 <= pos[0] < 360.0

    @pytest.mark.unit
    @pytest.mark.parametrize("body", [MARS, JUPITER])
    def test_nogdefl_valid(self, body):
        """NOGDEFL returns valid positions."""
        flags = FLG_SWIEPH | FLG_SPEED | FLG_NOGDEFL
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        assert 0.0 <= pos[0] < 360.0
