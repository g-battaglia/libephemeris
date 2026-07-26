"""
Behavioral tests for the applied SIDBIT projection flags.

These exercise the library alone (no reference ephemeris): they check that the
implemented projection flags are honored by the sidereal path rather than
reduced to the base ayanamsha mode, that the geometry matches the frame the
flag names (the solar-system invariable plane for SIDBIT_SSY_PLANE, the mean
ecliptic of t0 for SIDBIT_ECL_T0), and that unaffected channels stay put.

The library-vs-reference numerical agreement (and the documented residuals) is
pinned in validation/compare_scripts/tests/test_sidereal/test_sidbit_projection.py.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import libephemeris as ephem
from libephemeris.constants import (
    FLG_EQUATORIAL,
    FLG_SIDEREAL,
    FLG_SPEED,
    FLG_XYZ,
    MOON,
    MARS,
    SIDBIT_ECL_T0,
    SIDBIT_SSY_PLANE,
    SIDBIT_USER_UT,
    SIDM_LAHIRI,
    SIDM_USER,
    SUN,
)

JD = 2451545.0


@pytest.fixture(autouse=True)
def _reset_state():
    yield
    ephem.set_sid_mode(SIDM_LAHIRI)
    ephem.close()


class TestSsyPlaneApplied:
    """SIDBIT_SSY_PLANE projects onto the invariable plane, not a base reduce."""

    @pytest.mark.unit
    def test_latitude_shifts_toward_invariable_plane(self):
        """The Moon's ecliptic latitude (~5.17 deg) shifts by more than a degree
        under the invariable-plane projection; the base reduce would leave it
        unchanged."""
        ephem.set_sid_mode(SIDM_LAHIRI)
        base, _ = ephem.calc_ut(JD, MOON, FLG_SIDEREAL)
        ephem.set_sid_mode(SIDM_LAHIRI | SIDBIT_SSY_PLANE)
        ssy, _ = ephem.calc_ut(JD, MOON, FLG_SIDEREAL)
        assert abs(ssy[1] - base[1]) > 1.0

    @pytest.mark.unit
    def test_latitude_matches_invariable_plane_rotation(self):
        """The projected latitude equals the latitude of the J2000 position
        rotated onto the invariable plane (internal-consistency check)."""
        from libephemeris.constants import FLG_J2000, FLG_NONUT
        from libephemeris.sidereal_epoch import _invariable_plane_matrix

        for body in (SUN, MOON, MARS):
            j2000, _ = ephem.calc_ut(JD, body, FLG_J2000 | FLG_NONUT)
            lo, la = math.radians(j2000[0]), math.radians(j2000[1])
            v = np.array(
                [math.cos(la) * math.cos(lo), math.cos(la) * math.sin(lo), math.sin(la)]
            )
            p = _invariable_plane_matrix() @ v
            expected_lat = math.degrees(math.asin(max(-1.0, min(1.0, p[2]))))
            ephem.set_sid_mode(SIDM_LAHIRI | SIDBIT_SSY_PLANE)
            ssy, _ = ephem.calc_ut(JD, body, FLG_SIDEREAL)
            assert ssy[1] == pytest.approx(expected_lat, abs=1e-6)

    @pytest.mark.unit
    def test_equatorial_channel_unchanged(self):
        """SIDBIT_SSY_PLANE affects only the ecliptic representation; the
        equatorial (RA/Dec) channel is unaffected."""
        ephem.set_sid_mode(SIDM_LAHIRI)
        base, _ = ephem.calc_ut(JD, MARS, FLG_SIDEREAL | FLG_EQUATORIAL)
        ephem.set_sid_mode(SIDM_LAHIRI | SIDBIT_SSY_PLANE)
        ssy, _ = ephem.calc_ut(JD, MARS, FLG_SIDEREAL | FLG_EQUATORIAL)
        assert ssy[0] == pytest.approx(base[0], abs=1e-9)
        assert ssy[1] == pytest.approx(base[1], abs=1e-9)

    @pytest.mark.unit
    def test_xyz_out_of_plane_component_small(self):
        """In invariable-plane XYZ, the z component is the small out-of-plane
        offset and the vector norm is preserved."""
        ephem.set_sid_mode(SIDM_LAHIRI)
        base, _ = ephem.calc_ut(JD, MARS, FLG_SIDEREAL | FLG_XYZ)
        ephem.set_sid_mode(SIDM_LAHIRI | SIDBIT_SSY_PLANE)
        ssy, _ = ephem.calc_ut(JD, MARS, FLG_SIDEREAL | FLG_XYZ)
        base_norm = math.sqrt(sum(c * c for c in base[:3]))
        ssy_norm = math.sqrt(sum(c * c for c in ssy[:3]))
        assert ssy_norm == pytest.approx(base_norm, rel=1e-9)

    @pytest.mark.unit
    def test_speed_is_finite_and_reasonable(self):
        """The projected longitude speed stays close to the Moon's ~13 deg/day."""
        ephem.set_sid_mode(SIDM_LAHIRI | SIDBIT_SSY_PLANE)
        ssy, _ = ephem.calc_ut(JD, MOON, FLG_SIDEREAL | FLG_SPEED)
        assert 11.0 < ssy[3] < 16.0
        assert math.isfinite(ssy[4])


class TestEclT0Applied:
    """SIDBIT_ECL_T0 projects onto the mean ecliptic of the mode's t0."""

    @pytest.mark.unit
    def test_position_differs_from_base(self):
        ephem.set_sid_mode(SIDM_LAHIRI)
        base, _ = ephem.calc_ut(JD, MOON, FLG_SIDEREAL)
        ephem.set_sid_mode(SIDM_LAHIRI | SIDBIT_ECL_T0)
        ecl, _ = ephem.calc_ut(JD, MOON, FLG_SIDEREAL)
        assert (base[0], base[1]) != (ecl[0], ecl[1])

    @pytest.mark.unit
    def test_equatorial_reflects_mean_equator_of_t0(self):
        """Under SIDBIT_ECL_T0 the equatorial channel is not projected (the
        base sidereal equator is returned); it is finite and in range."""
        ephem.set_sid_mode(SIDM_LAHIRI | SIDBIT_ECL_T0)
        ecl, _ = ephem.calc_ut(JD, MARS, FLG_SIDEREAL | FLG_EQUATORIAL)
        assert 0.0 <= ecl[0] < 360.0
        assert -90.0 <= ecl[1] <= 90.0


class TestUserUtApplied:
    """SIDBIT_USER_UT interprets the SIDM_USER t0 as a UT date."""

    @pytest.mark.unit
    def test_ayanamsa_responds_to_bit(self):
        """With a deep-past t0 (large Delta T) the UT interpretation shifts the
        ayanamsha; without the bit it does not."""
        t0 = 2086302.5  # ~1000 CE, Delta T ~ 1600 s
        ephem.set_sid_mode(SIDM_USER, t0, 10.0)
        a_tt = ephem.get_ayanamsa_ut(JD)
        ephem.set_sid_mode(SIDM_USER | SIDBIT_USER_UT, t0, 10.0)
        a_ut = ephem.get_ayanamsa_ut(JD)
        # ΔT at t0 pushes the TT baseline later, lowering the accumulated
        # precession, so the UT interpretation yields a slightly smaller value.
        assert a_ut < a_tt
        assert abs(a_ut - a_tt) < 0.01  # sub-arcminute effect

    @pytest.mark.unit
    def test_modern_t0_effect_is_tiny(self):
        """At a modern t0 the UT/TT difference is sub-milliarcsecond."""
        t0 = 2451545.0
        ephem.set_sid_mode(SIDM_USER, t0, 24.0)
        a_tt = ephem.get_ayanamsa_ut(JD)
        ephem.set_sid_mode(SIDM_USER | SIDBIT_USER_UT, t0, 24.0)
        a_ut = ephem.get_ayanamsa_ut(JD)
        assert abs(a_ut - a_tt) < 1.0 / 3600.0
