"""FLG_TRUEPOS / FLG_NOABERR remove the anchor's annual aberration.

The star / galactic-center ayanamsha modes anchor the sidereal zero point to a
live light source. Under a FLG_TRUEPOS or FLG_NOABERR request the anchor's
annual aberration (Bradley 1728; IAU aberration constant kappa ~ 20.49552") is
removed from the zero-point reduction, so the subtracted ayanamsha matches the
aberration-free planet place the same request produces. Modes with a fixed
catalogue anchor (Aldebaran), a geometric galactic-pole direction (GALEQU) or a
precession-formula definition (Mardyks, Valens) are left unchanged.

These are library-internal consistency invariants; no reference output is read
or persisted.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import (
    FLG_NOABERR,
    FLG_SIDEREAL,
    FLG_SWIEPH,
    FLG_TRUEPOS,
    MARS,
    SIDM_ALDEBARAN_15TAU,
    SIDM_GALALIGN_MARDYKS,
    SIDM_GALCENT_0SAG,
    SIDM_GALCENT_COCHRANE,
    SIDM_GALCENT_MULA_WILHELM,
    SIDM_GALCENT_RGILBRAND,
    SIDM_GALEQU_IAU1958,
    SIDM_GALEQU_MULA,
    SIDM_GALEQU_TRUE,
    SIDM_TRUE_CITRA,
    SIDM_TRUE_MULA,
    SIDM_TRUE_PUSHYA,
    SIDM_TRUE_REVATI,
    SIDM_TRUE_SHEORAN,
    SIDM_VALENS_MOON,
)

# Annual aberration is bounded by the IAU aberration constant; the projection
# onto a single ecliptic longitude never exceeds it. A small margin covers the
# ecliptic-latitude foreshortening and light deflection.
ABERRATION_CONSTANT_ARCSEC = 20.49552
SHIFT_CEILING_ARCSEC = ABERRATION_CONSTANT_ARCSEC + 1.0

# A spread of epochs so at least one lands near the anchor's maximum aberration
# projection (the shift is a slow annual sinusoid in ecliptic longitude).
DATES = (2451545.0, 2455000.5, 2458849.5, 2460000.5, 2462000.5)

ABERRANT_MODES = (
    SIDM_TRUE_CITRA,
    SIDM_TRUE_REVATI,
    SIDM_TRUE_PUSHYA,
    SIDM_TRUE_MULA,
    SIDM_TRUE_SHEORAN,
    SIDM_GALCENT_0SAG,
    SIDM_GALCENT_RGILBRAND,
    SIDM_GALCENT_COCHRANE,
    SIDM_GALCENT_MULA_WILHELM,
)

UNAFFECTED_MODES = (
    SIDM_ALDEBARAN_15TAU,  # fixed catalogue anchor
    SIDM_GALEQU_IAU1958,  # geometric galactic pole (already aberration-free)
    SIDM_GALEQU_TRUE,
    SIDM_GALEQU_MULA,
    SIDM_GALALIGN_MARDYKS,  # precession-formula frame
    SIDM_VALENS_MOON,  # epoch-pair Method-B
)


def _signed_arcsec(delta_deg: float) -> float:
    return ((delta_deg + 180.0) % 360.0 - 180.0) * 3600.0


def _shift_arcsec(jd: float, extra_flag: int) -> float:
    _, plain = ephem.get_ayanamsa_ex_ut(jd, FLG_SWIEPH)
    _, shifted = ephem.get_ayanamsa_ex_ut(jd, FLG_SWIEPH | extra_flag)
    return _signed_arcsec(shifted - plain)


@pytest.mark.parametrize("sid_mode", ABERRANT_MODES)
def test_aberrant_anchor_shifts_under_truepos_and_noaberr(sid_mode: int) -> None:
    """The anchor moves by up to the aberration constant, and TP == NB in lib."""
    ephem.set_sid_mode(sid_mode)
    max_shift = 0.0
    for jd in DATES:
        tp = _shift_arcsec(jd, FLG_TRUEPOS)
        nb = _shift_arcsec(jd, FLG_NOABERR)
        # The library uses a single aberration-free toggle for both flags
        # (they differ only by sub-arcsecond light deflection in the reference).
        assert tp == pytest.approx(nb, abs=1e-9)
        assert abs(tp) <= SHIFT_CEILING_ARCSEC
        max_shift = max(max_shift, abs(tp))
    # Over the epoch spread the shift is a real, annual-aberration-scale effect.
    assert max_shift > 5.0


@pytest.mark.parametrize("sid_mode", UNAFFECTED_MODES)
def test_unaffected_modes_ignore_truepos_and_noaberr(sid_mode: int) -> None:
    ephem.set_sid_mode(sid_mode)
    for jd in DATES:
        assert _shift_arcsec(jd, FLG_TRUEPOS) == pytest.approx(0.0, abs=1e-6)
        assert _shift_arcsec(jd, FLG_NOABERR) == pytest.approx(0.0, abs=1e-6)


@pytest.mark.parametrize("jd", DATES)
def test_sidereal_calc_subtracts_the_aberration_free_ayanamsha(jd: float) -> None:
    """calc_ut(SIDEREAL|TRUEPOS) removes the aberration-free anchor ayanamsha.

    The of-date sidereal longitude is the tropical (geometric, TRUEPOS) longitude
    minus the TRUE ayanamsha (mean + nutation) for the same flags, which under
    TRUEPOS is aberration-free. Reconstructing it from the public ayanamsha getter
    proves the sidereal path uses the aberration-free value rather than the plain
    (aberrant) one.
    """
    ephem.set_sid_mode(SIDM_TRUE_CITRA)

    tropical = ephem.calc_ut(jd, MARS, FLG_SWIEPH | FLG_TRUEPOS)[0][0]
    _, aya_truepos = ephem.get_ayanamsa_ex_ut(jd, FLG_SWIEPH | FLG_TRUEPOS)
    _, aya_plain = ephem.get_ayanamsa_ex_ut(jd, FLG_SWIEPH)

    sid_actual = ephem.calc_ut(jd, MARS, FLG_SWIEPH | FLG_SIDEREAL | FLG_TRUEPOS)[0][0]
    expected = (tropical - aya_truepos) % 360.0

    assert _signed_arcsec(sid_actual - expected) == pytest.approx(0.0, abs=0.05)

    # And it is genuinely different from the buggy plain-ayanamsha subtraction
    # whenever the anchor's aberration is non-negligible at this epoch.
    buggy = (tropical - aya_plain) % 360.0
    if abs(_signed_arcsec(aya_truepos - aya_plain)) > 1.0:
        assert abs(_signed_arcsec(sid_actual - buggy)) > 0.5


class TestAnchoredAberrationOnHousesAndFixstar:
    """The houses and fixed-star sidereal paths consume the same anchored
    ayanamsha aberration toggle as the calc chain (measured reference
    behavior: every cusp/angle and star longitude shifts uniformly under
    FLG_TRUEPOS/FLG_NOABERR for the anchored modes, and not at all for the
    epoch-anchored ones)."""

    def test_houses_and_fixstar_shift_for_anchored_mode(self):
        import libephemeris as le

        jd = 2437500.5
        f = le.FLG_SWIEPH | le.FLG_SIDEREAL
        le.set_sid_mode(27, 0.0, 0.0)  # True Citra (anchored)
        a0 = le.houses_ex(jd, 45.0, 12.5, b"P", f)[1][0]
        at = le.houses_ex(jd, 45.0, 12.5, b"P", f | le.FLG_TRUEPOS)[1][0]
        assert abs((at - a0 + 180.0) % 360.0 - 180.0) * 3600.0 > 0.5
        s0 = le.fixstar2_ut("Aldebaran", jd, f)[0][0]
        st = le.fixstar2_ut("Aldebaran", jd, f | le.FLG_TRUEPOS)[0][0]
        # Star shift = star's own aberration + anchor toggle; differs from
        # the non-anchored baseline below by the anchor's response.
        anchored = ((st - s0 + 180.0) % 360.0 - 180.0) * 3600.0
        le.set_sid_mode(1, 0.0, 0.0)  # Lahiri (epoch-anchored)
        b0 = le.fixstar2_ut("Aldebaran", jd, f)[0][0]
        bt = le.fixstar2_ut("Aldebaran", jd, f | le.FLG_TRUEPOS)[0][0]
        plainmode = ((bt - b0 + 180.0) % 360.0 - 180.0) * 3600.0
        assert abs(anchored - plainmode) > 0.5
        # Epoch-anchored houses stay put.
        h0 = le.houses_ex(jd, 45.0, 12.5, b"P", f)[1][0]
        ht = le.houses_ex(jd, 45.0, 12.5, b"P", f | le.FLG_TRUEPOS)[1][0]
        assert abs((ht - h0 + 180.0) % 360.0 - 180.0) * 3600.0 < 0.001
