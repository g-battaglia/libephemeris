"""Property tests for SIDBIT projection support in sidereal house cusps.

These verify the *structural* behaviour of ``houses_ex`` / ``houses_ex2`` under
the ``SIDBIT_ECL_T0`` (256) and ``SIDBIT_SSY_PLANE`` (512) projection flags,
without any reference-API cusp values:

* On an epoch-anchored ayanamsha mode (Lahiri) the projection bits are NOT a
  no-op -- they move the cusps and the reported ARMC.
* ``SIDBIT_SSY_PLANE`` (a strongly tilted plane) shifts the cusps by degrees;
  ``SIDBIT_ECL_T0`` (a plane near the ecliptic of date) shifts them by only a
  small amount.
* On the star/galactic "true" modes the projection is inert (the same
  suppression the position path applies), so the bit changes nothing.
* The ascendant-anchored systems keep their defining geometry (Equal stays
  30 deg-spaced, Aries stays pinned to 0 deg) in the projected zodiac.
* ``houses_ex2`` speeds stay the coherent finite difference of the projected
  ``houses_ex`` positions, and ``FLG_RADIANS`` is honoured.

Provenance: structural/self-consistency assertions only. No cusp table or
angle value from another implementation is used as an input or an oracle.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    FLG_RADIANS,
    FLG_SIDEREAL,
    SIDBIT_ECL_T0,
    SIDBIT_SSY_PLANE,
    SIDM_LAHIRI,
    SIDM_TRUE_CITRA,
)

JD = 2451545.0  # J2000.0
JD_1972 = swe.julday(1972, 2, 29, 15.0)
JD_2028 = swe.julday(2028, 10, 12, 3.0)
LAT = 41.9028
LON = 12.4964


def _adiff(a: float, b: float) -> float:
    """Signed smallest angular difference a-b in (-180, 180]."""
    return (a - b + 180.0) % 360.0 - 180.0


@pytest.fixture(autouse=True)
def _reset_state():
    yield
    swe.close()


class TestSidbitNotNoOp:
    """The projection bits actually change epoch-mode sidereal houses."""

    @pytest.mark.unit
    @pytest.mark.parametrize("bit", [SIDBIT_ECL_T0, SIDBIT_SSY_PLANE])
    def test_projection_is_applied(self, bit):
        # The projection rebuilds the construction on the reference plane, which
        # re-bases the reported ARMC (ECL_T0 by a few arcsec, SSY_PLANE by
        # degrees) and moves the cusps -- previously the bits were a silent
        # no-op. JD_1972 is an epoch where both signatures are clearly nonzero
        # for both planes (the ECL_T0 cusp shift is frame-dependent and can
        # shrink to sub-arcsecond at other dates, e.g. near its t0 = J2000).
        swe.set_sid_mode(SIDM_LAHIRI)
        plain, plain_a = swe.houses_ex(JD_1972, LAT, LON, ord("P"), FLG_SIDEREAL)
        swe.set_sid_mode(SIDM_LAHIRI | bit)
        proj, proj_a = swe.houses_ex(JD_1972, LAT, LON, ord("P"), FLG_SIDEREAL)
        # The reconstruction re-bases the reported ARMC.
        assert abs(_adiff(proj_a[2], plain_a[2])) > 1.0 / 3600.0
        # At least one cusp moves measurably.
        assert max(abs(_adiff(p, q)) for p, q in zip(proj, plain)) > 1.0 / 3600.0
        # Every cusp is a valid longitude.
        for c in proj:
            assert 0.0 <= c < 360.0

    @pytest.mark.unit
    def test_ssy_plane_moves_more_than_ecl_t0(self):
        """The tilted invariable plane shifts the Ascendant by degrees; the
        near-ecliptic t0 plane shifts it only slightly."""
        swe.set_sid_mode(SIDM_LAHIRI)
        _, plain_a = swe.houses_ex(JD, LAT, LON, ord("P"), FLG_SIDEREAL)

        swe.set_sid_mode(SIDM_LAHIRI | SIDBIT_ECL_T0)
        _, ecl_a = swe.houses_ex(JD, LAT, LON, ord("P"), FLG_SIDEREAL)
        swe.set_sid_mode(SIDM_LAHIRI | SIDBIT_SSY_PLANE)
        _, ssy_a = swe.houses_ex(JD, LAT, LON, ord("P"), FLG_SIDEREAL)

        d_ecl = abs(_adiff(ecl_a[0], plain_a[0]))
        d_ssy = abs(_adiff(ssy_a[0], plain_a[0]))
        assert d_ecl < 0.1  # < 6 arcmin
        assert d_ssy > 1.0  # > 1 degree
        assert d_ssy > d_ecl

    @pytest.mark.unit
    def test_ssy_plane_shifts_reported_armc(self):
        """A full reconstruction on the tilted plane moves the reported ARMC
        (a per-point projection would leave the equatorial ARMC untouched)."""
        swe.set_sid_mode(SIDM_LAHIRI)
        _, plain_a = swe.houses_ex(JD, LAT, LON, ord("P"), FLG_SIDEREAL)
        swe.set_sid_mode(SIDM_LAHIRI | SIDBIT_SSY_PLANE)
        _, ssy_a = swe.houses_ex(JD, LAT, LON, ord("P"), FLG_SIDEREAL)
        assert abs(_adiff(ssy_a[2], plain_a[2])) > 1.0  # ARMC shifts > 1 deg


class TestSidbitSuppressed:
    """On the star/galactic 'true' modes the projection is inert."""

    @pytest.mark.unit
    @pytest.mark.parametrize("bit", [SIDBIT_ECL_T0, SIDBIT_SSY_PLANE])
    def test_true_mode_projection_inert(self, bit):
        swe.set_sid_mode(SIDM_TRUE_CITRA)
        base_c, base_a = swe.houses_ex(JD_2028, LAT, LON, ord("P"), FLG_SIDEREAL)
        swe.set_sid_mode(SIDM_TRUE_CITRA | bit)
        proj_c, proj_a = swe.houses_ex(JD_2028, LAT, LON, ord("P"), FLG_SIDEREAL)
        for a, b in zip(base_c, proj_c):
            assert abs(_adiff(a, b)) < 1e-9
        for a, b in zip(base_a, proj_a):
            assert abs(_adiff(a, b)) < 1e-9


class TestSidbitAscAnchoredSystems:
    """Ascendant-anchored systems keep their geometry in the projected zodiac."""

    @pytest.mark.unit
    @pytest.mark.parametrize("bit", [SIDBIT_ECL_T0, SIDBIT_SSY_PLANE])
    def test_equal_stays_30_spaced(self, bit):
        swe.set_sid_mode(SIDM_LAHIRI | bit)
        cusps, ascmc = swe.houses_ex(JD_2028, LAT, LON, ord("E"), FLG_SIDEREAL)
        assert abs(_adiff(cusps[0], ascmc[0])) < 1e-6  # cusp 1 == Asc
        for i in range(12):
            assert abs(_adiff(cusps[i], (ascmc[0] + 30.0 * i))) < 1e-6

    @pytest.mark.unit
    @pytest.mark.parametrize("bit", [SIDBIT_ECL_T0, SIDBIT_SSY_PLANE])
    def test_aries_stays_pinned(self, bit):
        swe.set_sid_mode(SIDM_LAHIRI | bit)
        cusps, _ = swe.houses_ex(JD_2028, LAT, LON, ord("N"), FLG_SIDEREAL)
        for i in range(12):
            assert abs(_adiff(cusps[i], 30.0 * i)) < 1e-9


class TestSidbitSpeedsAndRadians:
    """houses_ex2 speeds and FLG_RADIANS behave under the projection."""

    @pytest.mark.unit
    @pytest.mark.parametrize("bit", [SIDBIT_ECL_T0, SIDBIT_SSY_PLANE])
    def test_speeds_are_finite_difference(self, bit):
        swe.set_sid_mode(SIDM_LAHIRI | bit)
        _, _, cusps_speed, ascmc_speed = swe.houses_ex2(
            JD_2028, LAT, LON, ord("P"), FLG_SIDEREAL
        )
        dt = 1.0 / 4096.0  # the implementation's exact-binary half-step
        cm, am = swe.houses_ex(JD_2028 - dt, LAT, LON, ord("P"), FLG_SIDEREAL)
        cp, ap = swe.houses_ex(JD_2028 + dt, LAT, LON, ord("P"), FLG_SIDEREAL)
        for i in range(12):
            fd = _adiff(cp[i], cm[i]) / (2.0 * dt)
            assert abs(fd - cusps_speed[i]) < 1e-6
        # ARMC advances at the sidereal rate (~360.99 deg/day), well > 0.
        assert ascmc_speed[2] > 300.0

    @pytest.mark.unit
    @pytest.mark.parametrize("bit", [SIDBIT_ECL_T0, SIDBIT_SSY_PLANE])
    def test_radians_matches_degrees(self, bit):
        swe.set_sid_mode(SIDM_LAHIRI | bit)
        cd, ad = swe.houses_ex(JD_2028, LAT, LON, ord("P"), FLG_SIDEREAL)
        cr, ar = swe.houses_ex(JD_2028, LAT, LON, ord("P"), FLG_SIDEREAL | FLG_RADIANS)
        for i in range(12):
            assert abs(_adiff(cd[i], math.degrees(cr[i]))) < 1e-9
        for i in range(8):
            assert abs(_adiff(ad[i], math.degrees(ar[i]))) < 1e-9
