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
        """Under SIDBIT_ECL_T0 the equatorial channel is the mean frame of
        the mode's t0 (see TestEclT0EquatorialFrame); it is finite and in
        range."""
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


class TestUserEclT0Epoch:
    """SIDM_USER's ECL_T0 projection uses the user-supplied t0, not J2000.

    Measured reference behavior: the projection plane is the mean ecliptic
    of the t0 passed to set_sid_mode(SIDM_USER, t0, ayan_t0), taken
    literally. The ECL_T0 latitude shift therefore vanishes when t0 equals
    the computation date and grows by ~0.36 arcsec per year of |date - t0|
    (the secular inclination motion of the ecliptic; Vondrak et al. 2011,
    A&A 534, A22). Before the fix the epoch resolver silently fell back to
    J2000, freezing the plane regardless of t0.
    """

    JD_1990 = 2448076.5 + 0.5  # 1990-07-04 12:00 UT

    def _lat_shift(self, t0: float, body: int = ephem.PLUTO) -> float:
        ephem.set_sid_mode(SIDM_USER, t0, 24.0)
        plain = ephem.calc_ut(self.JD_1990, body, FLG_SIDEREAL)[0][1]
        ephem.set_sid_mode(SIDM_USER | SIDBIT_ECL_T0, t0, 24.0)
        proj = ephem.calc_ut(self.JD_1990, body, FLG_SIDEREAL)[0][1]
        return (proj - plain) * 3600.0

    @pytest.mark.unit
    def test_epoch_resolver_returns_user_t0(self):
        from libephemeris.planets import _ecl_t0_epoch_jd

        t0 = 2430000.25
        ephem.set_sid_mode(SIDM_USER, t0, 0.0)
        assert _ecl_t0_epoch_jd(SIDM_USER) == pytest.approx(t0, abs=1e-9)

    @pytest.mark.unit
    def test_epoch_resolver_converts_user_ut_t0(self):
        from libephemeris.planets import _ecl_t0_epoch_jd
        from libephemeris.time_utils import deltat

        t0 = 2430000.25
        ephem.set_sid_mode(SIDM_USER | SIDBIT_USER_UT, t0, 0.0)
        assert _ecl_t0_epoch_jd(SIDM_USER) == pytest.approx(t0 + deltat(t0), abs=1e-9)

    @pytest.mark.unit
    def test_shift_vanishes_at_t0_equal_date(self):
        # Plane of t0 == plane of date: the projection is a longitude-only
        # rezeroing, so the latitude shift collapses to the numerical floor.
        assert abs(self._lat_shift(self.JD_1990)) < 0.1

    @pytest.mark.unit
    def test_shift_scales_with_epoch_separation(self):
        # ~90 years of secular inclination motion: tens of arcseconds,
        # strictly larger than the ~36-year J2000 baseline. A frozen-J2000
        # plane would make both values identical.
        shift_1900 = self._lat_shift(2415021.0)  # 1900-01-01 12:00
        shift_j2000 = self._lat_shift(JD)
        assert abs(shift_1900) > 20.0
        assert abs(shift_1900) > 2.0 * abs(shift_j2000)


class TestEclT0EquatorialFrame:
    """Equatorial output under SIDBIT_ECL_T0 is the t0 mean-frame position.

    Measured reference behavior (calc and fixstar paths alike): RA/Dec and
    their speeds are the FLG_J2000|FLG_NONUT request reduced to the mean
    equator and equinox of the mode's t0, with no ayanamsha subtraction and
    the plain sidereal retflag echo. SIDBIT_SSY_PLANE alone leaves the calc
    equatorial channel unchanged.
    """

    @pytest.mark.unit
    def test_user_t0_j2000_matches_j2000_request(self):
        from libephemeris.constants import FLG_J2000

        # t0 = J2000: the t0 reduction is the identity, so the output must
        # be bit-identical to the plain FLG_J2000|FLG_NONUT request.
        ephem.set_sid_mode(SIDM_USER | SIDBIT_ECL_T0, JD, 24.0)
        proj, rf = ephem.calc_ut(JD + 6700.0, MARS, FLG_SIDEREAL | FLG_EQUATORIAL)
        base, _ = ephem.calc_ut(
            JD + 6700.0, MARS, FLG_J2000 | FLG_NONUT | FLG_EQUATORIAL
        )
        assert proj[0] == pytest.approx(base[0], abs=1e-9)
        assert proj[1] == pytest.approx(base[1], abs=1e-9)
        assert rf & FLG_SIDEREAL and rf & FLG_NONUT

    @pytest.mark.unit
    def test_user_t0_1900_moves_the_frame(self):
        from libephemeris.constants import FLG_J2000

        # t0 = 1900: the frame differs from J2000 by ~a century of general
        # precession (~46"/yr in RA), and from the plain sidereal equator by
        # the ayanamsha-free reduction: both deltas are large.
        ephem.set_sid_mode(SIDM_USER | SIDBIT_ECL_T0, 2415021.0, 24.0)
        proj, _ = ephem.calc_ut(JD + 6700.0, MARS, FLG_SIDEREAL | FLG_EQUATORIAL)
        j2k, _ = ephem.calc_ut(
            JD + 6700.0, MARS, FLG_J2000 | FLG_NONUT | FLG_EQUATORIAL
        )
        d_ra = abs((proj[0] - j2k[0] + 180.0) % 360.0 - 180.0) * 3600.0
        assert 3000.0 < d_ra < 8000.0

    @pytest.mark.unit
    def test_ssy_plane_alone_keeps_calc_equator(self):
        ephem.set_sid_mode(SIDM_LAHIRI)
        base, _ = ephem.calc_ut(JD, MARS, FLG_SIDEREAL | FLG_EQUATORIAL)
        ephem.set_sid_mode(SIDM_LAHIRI | SIDBIT_SSY_PLANE)
        ssy, _ = ephem.calc_ut(JD, MARS, FLG_SIDEREAL | FLG_EQUATORIAL)
        assert ssy[0] == pytest.approx(base[0], abs=1e-9)
        assert ssy[1] == pytest.approx(base[1], abs=1e-9)

    @pytest.mark.unit
    def test_lahiri_classical_projection_epoch(self):
        from libephemeris.planets import _ecl_t0_epoch_jd

        # The Lahiri VALUE realization is anchored at J2000 (IAE), but the
        # ECL_T0 projection uses the Calendar Reform Committee defining
        # epoch: 23°15' at the vernal equinox of 1956 (CSIR 1955).
        assert _ecl_t0_epoch_jd(SIDM_LAHIRI) == pytest.approx(2435553.5)


class TestSidbitFixstar:
    """SIDBIT_ECL_T0 / SIDBIT_SSY_PLANE are applied on the fixed-star path.

    Measured reference behavior: the star pipeline projects ecliptic output
    onto the plane (same machinery as calc) and reduces equatorial output
    to the mean frame of the projection epoch; before the fix both bits
    were silent no-ops on fixstar (0.000" shift) while the reference moved
    the star by up to ~400" (ecliptic SSY) / ~24° (equatorial).
    """

    @pytest.mark.unit
    def test_ecl_t0_lat_shift_vanishes_at_t0_equal_date(self):
        # SIDM_USER with t0 = date: the projection plane IS the plane of
        # date, so the star's latitude shift collapses; a distant t0 moves
        # it by the secular inclination drift.
        jd = 2448077.0
        ephem.set_sid_mode(SIDM_USER, jd, 24.0)
        plain = ephem.fixstar_ut("Aldebaran", jd, FLG_SIDEREAL)[0]
        ephem.set_sid_mode(SIDM_USER | SIDBIT_ECL_T0, jd, 24.0)
        proj = ephem.fixstar_ut("Aldebaran", jd, FLG_SIDEREAL)[0]
        assert abs(proj[1] - plain[1]) * 3600.0 < 0.1
        ephem.set_sid_mode(SIDM_USER | SIDBIT_ECL_T0, 2415021.0, 24.0)
        far = ephem.fixstar_ut("Aldebaran", jd, FLG_SIDEREAL)[0]
        assert abs(far[1] - plain[1]) * 3600.0 > 10.0

    @pytest.mark.unit
    def test_ssy_plane_moves_ecliptic_latitude(self):
        # Aldebaran sits ~5.5 deg below the ecliptic; the invariable plane
        # (Souami & Souchay 2012) is inclined ~1.58 deg, so the projected
        # latitude moves by a large fraction of a degree.
        ephem.set_sid_mode(SIDM_LAHIRI)
        plain = ephem.fixstar_ut("Aldebaran", JD, FLG_SIDEREAL)[0]
        ephem.set_sid_mode(SIDM_LAHIRI | SIDBIT_SSY_PLANE)
        ssy = ephem.fixstar_ut("Aldebaran", JD, FLG_SIDEREAL)[0]
        assert abs(ssy[1] - plain[1]) * 3600.0 > 300.0

    @pytest.mark.unit
    def test_ssy_equatorial_is_j2000_frame(self):
        from libephemeris.constants import FLG_J2000

        ephem.set_sid_mode(SIDM_LAHIRI | SIDBIT_SSY_PLANE)
        proj = ephem.fixstar_ut("Aldebaran", JD, FLG_SIDEREAL | FLG_EQUATORIAL)[0]
        base = ephem.fixstar_ut(
            "Aldebaran", JD, FLG_J2000 | FLG_NONUT | FLG_EQUATORIAL
        )[0]
        assert proj[0] == pytest.approx(base[0], abs=1e-9)
        assert proj[1] == pytest.approx(base[1], abs=1e-9)

    @pytest.mark.unit
    def test_suppressed_mode_stays_inert_on_fixstar(self):
        ephem.set_sid_mode(SIDM_TRUE_CITRA)
        plain = ephem.fixstar_ut("Aldebaran", JD, FLG_SIDEREAL)[0]
        ephem.set_sid_mode(SIDM_TRUE_CITRA | SIDBIT_ECL_T0)
        proj = ephem.fixstar_ut("Aldebaran", JD, FLG_SIDEREAL)[0]
        assert proj[0] == pytest.approx(plain[0], abs=1e-9)
        assert proj[1] == pytest.approx(plain[1], abs=1e-9)


class TestClassicalEclT0Epochs:
    """Pin the classical SIDBIT_ECL_T0 projection-plane epochs (finding AI-F1).

    The SIDBIT_ECL_T0 projection must refer each mode to the mean ecliptic and
    equinox of its CLASSICAL DEFINING epoch, which for these modes differs from
    the epoch of the value anchor stored in ``AYANAMSHA_DEFINING``. Resolving
    the plane from the value anchor mis-orients the mean equator/ecliptic by up
    to ~29° in RA and tens of arcsec in ecliptic latitude. Each epoch pinned
    below is the round conventional epoch named by the mode's published
    defining statement:

      * DELUCE -> year 0.0 (beginning of the Christian era). Robert De Luce,
        "Constellational Astrology According to the Hindu System" (De Luce
        Publishing Co., Los Angeles, 1963): the constellational and sign
        zodiacs coincide (ayanamsha zero) at the start of the Christian era.
        JD 1721057.5 = 0000 Jan 1.0 (Julian).
      * DJWHAL_KHUL -> 1900.0. The Bailey/Djwhal Khul Aquarian-age doctrine
        (channelled 1940; astronomical realization in Phillip Lindsay, "The
        Beginning of the Age of Aquarius", 2006) sets the Age of Aquarius at
        2117 = 30° and specifies the value via the Synetic Vernal Point at the
        standard epoch 1900.0. JD 2415020.0 = 1900 Jan 0.5.
      * LAHIRI_1940 -> 1900.0. N.C. Lahiri, "Indian Ephemeris of Planets'
        Positions" (1st ed. 1939/1940) — the early Lahiri tradition, tabulated
        at the standard epoch 1900.0 (as are the Raman and Krishnamurti Indian
        ayanamshas). JD 2415020.0 = 1900 Jan 0.5.
      * BABYL_BRITTON -> year 0.0 (beginning of the common era). J.P. Britton,
        "Studies in Babylonian lunar theory: Part III. The introduction of the
        uniform zodiac" (Archive for History of Exact Sciences 64, 2010)
        writes the displacement as Δλ* = C − 1.3828°·Y with C = 3.20° and Y in
        centuries from the common era, so the constant term is anchored at
        Y = 0 (year 0.0). JD 1721057.5 = 0000 Jan 1.0 (Julian).

    One of the five modes in the finding is intentionally NOT overridden here;
    ``test_aldebaran_15tau_pending_projection_range_support`` documents why and
    pins the current (un-overridden) fallback so the exclusion is explicit and
    cannot be silently reverted.
    """

    @pytest.mark.unit
    def test_deluce_plane_is_christian_era(self):
        from libephemeris.constants import SIDM_DELUCE
        from libephemeris.planets import _ecl_t0_epoch_jd

        assert _ecl_t0_epoch_jd(SIDM_DELUCE) == pytest.approx(1721057.5)

    @pytest.mark.unit
    def test_djwhal_khul_plane_is_1900(self):
        from libephemeris.constants import SIDM_DJWHAL_KHUL
        from libephemeris.planets import _ecl_t0_epoch_jd

        assert _ecl_t0_epoch_jd(SIDM_DJWHAL_KHUL) == pytest.approx(2415020.0)

    @pytest.mark.unit
    def test_lahiri_1940_plane_is_1900(self):
        from libephemeris.constants import SIDM_LAHIRI_1940
        from libephemeris.planets import _ecl_t0_epoch_jd

        assert _ecl_t0_epoch_jd(SIDM_LAHIRI_1940) == pytest.approx(2415020.0)

    @pytest.mark.unit
    def test_britton_plane_is_common_era(self):
        from libephemeris.constants import SIDM_BABYL_BRITTON
        from libephemeris.planets import _ecl_t0_epoch_jd

        assert _ecl_t0_epoch_jd(SIDM_BABYL_BRITTON) == pytest.approx(1721057.5)

    @pytest.mark.unit
    def test_aldebaran_15tau_keeps_the_tier_independent_j2000_plane(self):
        """ALDEBARAN_15TAU (14) stays on the J2000 fallback plane.

        Its classical defining plane is the year -100 mean ecliptic (JD
        1684532.5) of the Babylonian normal-star zodiac the mode belongs to
        (C. Fagan, "Zodiacs Old and New", Llewellyn, 1950; P. Huber,
        Centaurus 5, 1958), the same plane the Kugler family uses. Applying
        it is blocked by a data-provisioning problem, not a sourcing one:
        this is a live star-anchored mode, so its ECL_T0 zero point is
        evaluated AT the plane epoch, and year -100 lies outside the
        base/medium DE440 coverage. A revision that fell back to a
        catalog-only anchor there made the same public call return two
        longitudes 0.79 deg apart on the medium and extended tiers — a
        worse defect than the divergence it closed. The epoch therefore
        stays at the J2000 fallback until the zero point can be evaluated
        tier-independently; the divergence is documented in
        docs/comparison/known-differences.md.
        """
        from libephemeris.constants import SIDM_ALDEBARAN_15TAU
        from libephemeris.planets import _ECL_T0_CLASSICAL_EPOCHS, _ecl_t0_epoch_jd

        assert SIDM_ALDEBARAN_15TAU not in _ECL_T0_CLASSICAL_EPOCHS
        assert _ecl_t0_epoch_jd(SIDM_ALDEBARAN_15TAU) == pytest.approx(2451545.0)


class TestFixedEpochJ2000Echo:
    """An explicit caller FLG_J2000 bit is echoed back by the fixed-epoch
    sidereal modes (18/19/20/34): only the INTERNAL J2000 rewrite bit is
    hidden (measured reference retflag keeps the caller's bit)."""

    @pytest.mark.unit
    @pytest.mark.parametrize("mode", [18, 19, 20, 34])
    def test_caller_j2000_bit_echoed(self, mode):
        from libephemeris.constants import FLG_J2000, FLG_SPEED, FLG_SWIEPH

        ephem.set_sid_mode(mode, 0.0, 0.0)
        base = FLG_SWIEPH | FLG_SPEED | FLG_SIDEREAL
        rf_j = ephem.calc_ut(JD, SUN, base | FLG_J2000)[1]
        rf_p = ephem.calc_ut(JD, SUN, base)[1]
        assert rf_j & FLG_J2000
        assert rf_j == rf_p | FLG_J2000


class TestBabylonianFamilyEclT0Epochs:
    """The Kugler triplet shares one classical projection plane (the
    Babylonian norm epoch, year -100; Huber, Centaurus 5, 1958) and
    Hipparchos uses Hipparchus's own era (~-128, Mercier's Hipparchan
    norm): the three Kugler ECL_T0 increments must be identical and the
    epochs pinned."""

    @pytest.mark.unit
    def test_epochs_pinned(self):
        from libephemeris.planets import _ecl_t0_epoch_jd

        for mode in (9, 10, 11):
            assert _ecl_t0_epoch_jd(mode) == pytest.approx(1684532.5)
        assert _ecl_t0_epoch_jd(15) == pytest.approx(1674484.0)

    @pytest.mark.unit
    def test_kugler_triplet_shares_one_plane(self):
        jd = 2450614.5
        incrs = []
        for mode in (9, 10, 11):
            ephem.set_sid_mode(mode, 0.0, 0.0)
            plain = ephem.calc_ut(jd, SUN, FLG_SIDEREAL)[0][1]
            ephem.set_sid_mode(mode | SIDBIT_ECL_T0, 0.0, 0.0)
            proj = ephem.calc_ut(jd, SUN, FLG_SIDEREAL)[0][1]
            incrs.append((proj - plain) * 3600.0)
        assert incrs[0] == pytest.approx(incrs[1], abs=0.01)
        assert incrs[1] == pytest.approx(incrs[2], abs=0.01)
        assert abs(incrs[0]) > 500.0  # a genuinely ancient plane


class TestEclT0EquatorialCalcStarCoherence:
    """The calc and fixstar paths share one ECL_T0 equatorial reduction.

    Both `planets._sidbit_projection_calc` and
    `fixed_stars._sidbit_star_call` must reduce equatorial ECL_T0 output
    with the same transform (`transform_equatorial_epoch_result`) and the
    same projection epoch. This pins the compatibility contract whose
    docstrings previously contradicted each other: each path's projected
    output must equal its own FLG_J2000|FLG_NONUT baseline pushed through
    the shared transform at the mode's t0.
    """

    T0_1900 = 2415021.0

    def _expected(self, base, flags, t0_jd):
        from libephemeris.sidereal_epoch import transform_equatorial_epoch_result

        return transform_equatorial_epoch_result(tuple(base), flags, t0_jd)

    @pytest.mark.unit
    def test_calc_path_uses_shared_transform_at_t0(self):
        from libephemeris.constants import FLG_J2000
        from libephemeris.planets import _ecl_t0_epoch_jd

        ephem.set_sid_mode(SIDM_USER | SIDBIT_ECL_T0, self.T0_1900, 24.0)
        t0 = _ecl_t0_epoch_jd(SIDM_USER)
        assert t0 == pytest.approx(self.T0_1900)
        flags = FLG_SIDEREAL | FLG_EQUATORIAL | FLG_SPEED
        proj, _ = ephem.calc_ut(JD, MARS, flags)
        base, _ = ephem.calc_ut(
            JD, MARS, FLG_J2000 | FLG_NONUT | FLG_EQUATORIAL | FLG_SPEED
        )
        expected = self._expected(base, flags, t0)
        for i in range(6):
            assert proj[i] == pytest.approx(expected[i], abs=1e-9), f"slot {i}"

    @pytest.mark.unit
    def test_star_path_uses_shared_transform_at_t0(self):
        from libephemeris.constants import FLG_J2000
        from libephemeris.planets import _ecl_t0_epoch_jd

        ephem.set_sid_mode(SIDM_USER | SIDBIT_ECL_T0, self.T0_1900, 24.0)
        t0 = _ecl_t0_epoch_jd(SIDM_USER)
        flags = FLG_SIDEREAL | FLG_EQUATORIAL
        proj = ephem.fixstar_ut("Aldebaran", JD, flags)[0]
        base = ephem.fixstar_ut(
            "Aldebaran", JD, FLG_J2000 | FLG_NONUT | FLG_EQUATORIAL
        )[0]
        expected = self._expected(base, flags, t0)
        for i in range(3):
            assert proj[i] == pytest.approx(expected[i], abs=1e-9), f"slot {i}"

    @pytest.mark.unit
    def test_paths_agree_module_and_context(self):
        """The context entry point applies the same reduction as the module."""
        from libephemeris.context import EphemerisContext

        ephem.set_sid_mode(SIDM_USER | SIDBIT_ECL_T0, self.T0_1900, 24.0)
        flags = FLG_SIDEREAL | FLG_EQUATORIAL
        module_xx, module_rf = ephem.calc_ut(JD, MARS, flags)
        ctx = EphemerisContext()
        ctx.set_sid_mode(SIDM_USER | SIDBIT_ECL_T0, self.T0_1900, 24.0)
        ctx_xx, ctx_rf = ctx.calc_ut(JD, MARS, flags)
        assert ctx_rf == module_rf
        for i in range(3):
            assert ctx_xx[i] == pytest.approx(module_xx[i], abs=1e-9), f"slot {i}"

    @pytest.mark.unit
    def test_ssy_plane_divergence_is_the_documented_one(self):
        """Where the two paths differ (SSY_PLANE equatorial): the star path
        reduces to J2000 while the calc path leaves the equator unchanged."""
        from libephemeris.constants import FLG_J2000

        ephem.set_sid_mode(SIDM_LAHIRI | SIDBIT_SSY_PLANE)
        flags = FLG_SIDEREAL | FLG_EQUATORIAL
        star_proj = ephem.fixstar_ut("Aldebaran", JD, flags)[0]
        star_j2k = ephem.fixstar_ut(
            "Aldebaran", JD, FLG_J2000 | FLG_NONUT | FLG_EQUATORIAL
        )[0]
        assert star_proj[0] == pytest.approx(star_j2k[0], abs=1e-9)
        calc_proj, _ = ephem.calc_ut(JD, MARS, flags)
        ephem.set_sid_mode(SIDM_LAHIRI)
        calc_base, _ = ephem.calc_ut(JD, MARS, flags)
        assert calc_proj[0] == pytest.approx(calc_base[0], abs=1e-9)
