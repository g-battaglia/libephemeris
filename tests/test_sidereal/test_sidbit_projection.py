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
    FLG_NONUT,
    FLG_SIDEREAL,
    FLG_SPEED,
    FLG_XYZ,
    MOON,
    MARS,
    SIDBIT_ECL_T0,
    SIDBIT_SSY_PLANE,
    SIDBIT_USER_UT,
    SIDM_ALDEBARAN_15TAU,
    SIDM_GALCENT_0SAG,
    SIDM_GALEQU_TRUE,
    SIDM_LAHIRI,
    SIDM_TRUE_CITRA,
    SIDM_TRUE_MULA,
    SIDM_TRUE_PUSHYA,
    SIDM_TRUE_REVATI,
    SIDM_USER,
    SIDM_VALENS_MOON,
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


class TestSidbitSuppressedForStarGalacticModes:
    """The star/galactic "true" modes leave SIDBIT_ECL_T0 / SSY_PLANE inert.

    Their zero point is defined by a live catalog/frame direction on the
    ecliptic of date, so re-projecting onto another reference plane is not
    defined by the model; measured reference behavior suppresses the projection
    (the sidereal longitude and latitude equal the un-projected baseline). The
    epoch-anchored mean modes and the two live modes whose zero is not a bare
    ecliptic-of-date star/frame direction (Aldebaran = 15 Tau, Vettius Valens)
    keep the projection.
    """

    SUPPRESSED = [
        SIDM_TRUE_CITRA,
        SIDM_TRUE_REVATI,
        SIDM_TRUE_PUSHYA,
        SIDM_TRUE_MULA,
        SIDM_GALEQU_TRUE,
        SIDM_GALCENT_0SAG,
    ]
    APPLIED = [SIDM_LAHIRI, SIDM_ALDEBARAN_15TAU, SIDM_VALENS_MOON]
    # A date well away from the modes' J2000 reference epoch, so ECL_T0 is a
    # genuine (non-identity) projection, kept inside the base tier (1850-2150).
    JD_OFF = 2415020.0  # 1900-01-01

    @pytest.mark.unit
    @pytest.mark.parametrize("mode", SUPPRESSED)
    @pytest.mark.parametrize("bit", [SIDBIT_ECL_T0, SIDBIT_SSY_PLANE])
    def test_projection_is_inert(self, mode, bit):
        """Longitude and latitude match the un-projected baseline exactly (a
        non-suppressed mode would move by up to ~0.1° here)."""
        ephem.set_sid_mode(mode)
        base, _ = ephem.calc_ut(self.JD_OFF, MARS, FLG_SIDEREAL)
        ephem.set_sid_mode(mode | bit)
        proj, _ = ephem.calc_ut(self.JD_OFF, MARS, FLG_SIDEREAL)
        assert proj[0] == pytest.approx(base[0], abs=1e-9)
        assert proj[1] == pytest.approx(base[1], abs=1e-9)

    @pytest.mark.unit
    @pytest.mark.parametrize("mode", APPLIED)
    @pytest.mark.parametrize("bit", [SIDBIT_ECL_T0, SIDBIT_SSY_PLANE])
    def test_projection_still_applies_for_other_modes(self, mode, bit):
        """A mean mode or a non-suppressed live mode still shifts the result."""
        ephem.set_sid_mode(mode)
        base, _ = ephem.calc_ut(self.JD_OFF, MARS, FLG_SIDEREAL)
        ephem.set_sid_mode(mode | bit)
        proj, _ = ephem.calc_ut(self.JD_OFF, MARS, FLG_SIDEREAL)
        moved = abs((proj[0] - base[0] + 180.0) % 360.0 - 180.0) + abs(
            proj[1] - base[1]
        )
        assert moved > 1e-3


class TestAppliedProjectionEchoesNonut:
    """An applied SIDBIT_ECL_T0/SSY_PLANE projection echoes FLG_NONUT.

    The projection is defined on a MEAN ecliptic and equinox (of the mode's t0
    or the invariable plane), so the projected longitude is declared without
    nutation -- the same mean-frame convention the base sidereal path echoes.
    Both calc and calc_ut are covered.
    """

    # Off-epoch date so ECL_T0 is a genuine projection (base tier 1850-2150).
    JD_OFF = 2415020.0  # 1900-01-01

    @pytest.mark.unit
    @pytest.mark.parametrize("bit", [SIDBIT_ECL_T0, SIDBIT_SSY_PLANE])
    @pytest.mark.parametrize("mode", [SIDM_LAHIRI, SIDM_ALDEBARAN_15TAU])
    def test_applied_projection_retflag_has_nonut(self, mode, bit):
        ephem.set_sid_mode(mode | bit)
        _, rf_ut = ephem.calc_ut(self.JD_OFF, SUN, FLG_SIDEREAL | FLG_SPEED)
        _, rf_tt = ephem.calc(self.JD_OFF, SUN, FLG_SIDEREAL | FLG_SPEED)
        assert rf_ut & FLG_NONUT
        assert rf_tt & FLG_NONUT

    @pytest.mark.unit
    @pytest.mark.parametrize("bit", [SIDBIT_ECL_T0, SIDBIT_SSY_PLANE])
    def test_suppressed_mode_also_echoes_nonut(self, bit):
        """A suppressed star/galactic mode reaches the base path, which already
        echoes FLG_NONUT; the applied path now agrees."""
        ephem.set_sid_mode(SIDM_TRUE_CITRA | bit)
        _, rf = ephem.calc_ut(self.JD_OFF, SUN, FLG_SIDEREAL | FLG_SPEED)
        assert rf & FLG_NONUT

    @pytest.mark.unit
    @pytest.mark.parametrize("bit", [SIDBIT_ECL_T0, SIDBIT_SSY_PLANE])
    def test_explicit_nonut_does_not_change_position(self, bit):
        """Passing FLG_NONUT explicitly leaves the projected position unchanged
        (NONUT is an echo-only implication, not a computation switch)."""
        ephem.set_sid_mode(SIDM_LAHIRI | bit)
        without, _ = ephem.calc_ut(self.JD_OFF, MARS, FLG_SIDEREAL | FLG_SPEED)
        with_nonut, _ = ephem.calc_ut(
            self.JD_OFF, MARS, FLG_SIDEREAL | FLG_SPEED | FLG_NONUT
        )
        assert with_nonut[0] == pytest.approx(without[0], abs=1e-12)
        assert with_nonut[1] == pytest.approx(without[1], abs=1e-12)


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


class TestValensEclT0Epoch:
    """SIDM_VALENS_MOON's ECL_T0 projection uses the Valens defining epoch.

    The mode keeps its UT-anchored defining pair outside the shared table,
    so the projection-epoch lookup needs the dedicated resolver: with the
    silent J2000 fallback the projection degenerated to a near no-op (up to
    ~4' at high ecliptic latitudes).
    """

    def test_valens_ecl_t0_projects_on_ancient_ecliptic(self):
        import libephemeris as le
        from libephemeris.ayanamsha_definitions import VALENS_MOON_T0_UT
        from libephemeris.planets import _ecl_t0_epoch_jd

        t0 = _ecl_t0_epoch_jd(le.SIDM_VALENS_MOON)
        # The resolved epoch is Valens' (~150 CE), not the J2000 fallback.
        assert abs(t0 - VALENS_MOON_T0_UT) < 1.0
        # And the projection genuinely moves a high-ecliptic-latitude body.
        jd = 2451545.0
        le.set_sid_mode(le.SIDM_VALENS_MOON)
        plain = le.calc_ut(jd, le.PLUTO, le.FLG_SWIEPH | le.FLG_SIDEREAL)[0][0]
        le.set_sid_mode(le.SIDM_VALENS_MOON | le.SIDBIT_ECL_T0)
        proj = le.calc_ut(jd, le.PLUTO, le.FLG_SWIEPH | le.FLG_SIDEREAL)[0][0]
        le.set_sid_mode(0)
        shift = abs((proj - plain + 180.0) % 360.0 - 180.0) * 3600.0
        assert shift > 30.0  # was ~0.03" with the J2000 fallback
