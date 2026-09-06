# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Frame-consistency regression tests for the Vondrák 2011 long-term model.

These pin two extreme-date fixes that route ``azalt``/``azalt_rev`` and the
public ``sidtime``/``sidtime0`` onto the same long-term (Vondrák 2011) obliquity
and Greenwich mean sidereal time the rest of the library uses for
ecliptic-of-date positions and house cusps. Previously both bypassed it through
the IAU 2006 models (``erfa.obl06`` / ``erfa.gmst06``), which are internally
inconsistent with the library's own frame at deep-time epochs (~5.7" of
obliquity and ~0.66° of GMST at year -3000).

The IAU 2006 models are used here only as the pre-fix baseline: inside the modern
window (1850-2050) the fixes are, by design, a no-op, so the modern result must
still match them to < 1e-6.

The last group pins how :mod:`libephemeris.sidereal_longterm` evaluates the
Vondrák model — the poles from ``erfa.ltpecl``/``erfa.ltpequ``, the of-date mean
obliquity being the angle between them, rebuilt here from ERFA directly — and
the properties the sidereal time has to have as a curve: one continuous,
steadily advancing turn of the Earth over the whole span, on the IAU 2006
expression in the era that expression was fitted for and off it outside.
"""

from __future__ import annotations

import math

import erfa
import numpy as np
import pytest

from libephemeris import sidereal_longterm as sl
from libephemeris import time_utils as tu
from libephemeris.cache import get_true_obliquity
from libephemeris.constants import FLG_EQUATORIAL, FLG_SWIEPH, SUN
from libephemeris.exceptions import EphemerisRangeError
from libephemeris.planets import calc_ut
from libephemeris.precession_vondrak import vondrak_mean_obliquity_deg
from libephemeris.state import get_precision_tier, get_timescale, set_precision_tier
from libephemeris.utils import (
    ECL2HOR,
    EQU2HOR,
    HOR2ECL,
    azalt,
    azalt_rev,
    cotrans,
)

_J2000 = 2451545.0
# Observer near Rome: (longitude E, latitude N, altitude m).
_GEOPOS = (12.5, 41.9, 0.0)
# An arbitrary ecliptic-of-date test position (lon, lat, dist); no ephemeris
# needed, so the self-consistency checks work at any epoch.
_ECL = (137.0, 3.2, 1.0)


def _jd_deep_bce() -> float:
    """Year -3000: well outside the IAU 2006 polynomial fit window."""
    return tu.julday(-3000, 1, 1, 12.0)


def _jd_modern() -> float:
    return tu.julday(2024, 6, 15, 12.0)


def _jd_tt(jd_ut1: float) -> float:
    return get_timescale().ut1_jd(jd_ut1).tt


def _obl06_deg(jd_tt: float) -> float:
    """IAU 2006 mean obliquity (deg) — the pre-fix obliquity source."""
    return math.degrees(erfa.obl06(_J2000, jd_tt - _J2000))


def _gmst06_deg(jd_ut1: float) -> float:
    """IAU 2006 GMST (deg, [0, 360)) — the pre-fix sidereal-time source."""
    jd_tt = jd_ut1 + tu.deltat(jd_ut1)
    return math.degrees(erfa.gmst06(jd_ut1, 0.0, jd_tt, 0.0)) % 360.0


def _wrap180(deg: float) -> float:
    return (deg + 180.0) % 360.0 - 180.0


# --------------------------------------------------------------------------
# Fix 2: sidtime routes through the Vondrák long-term GMST source
# --------------------------------------------------------------------------
def test_sidtime_matches_longterm_gmst_at_3000bce():
    """At year -3000 the public sidtime agrees with the house-engine GMST."""
    jd = _jd_deep_bce()
    # sidtime with nutation=0 is pure GMST (no equation of equinoxes).
    gmst_from_sidtime = (tu.sidtime(jd, 0.0, 23.44, 0.0) * 15.0) % 360.0
    gmst_house_engine = sl.mean_sidereal_time_deg(jd) % 360.0
    assert abs(_wrap180(gmst_from_sidtime - gmst_house_engine)) < 1e-6


def test_sidtime_diverged_from_iau2006_at_3000bce():
    """Guard: the IAU 2006 GMST is the ~0.66° outlier the fix moved off of."""
    jd = _jd_deep_bce()
    gmst_house_engine = sl.mean_sidereal_time_deg(jd) % 360.0
    gmst_iau2006 = _gmst06_deg(jd)
    # Measured divergence is ~0.664°; assert it is a large, real difference.
    assert abs(_wrap180(gmst_iau2006 - gmst_house_engine)) > 0.5


def test_sidtime_modern_unchanged():
    """Inside 1850-2050 sidtime still equals the IAU 2006 GMST (fix is a no-op)."""
    for jd in (_J2000, _jd_modern()):
        gmst_from_sidtime = (tu.sidtime(jd, 0.0, 23.44, 0.0) * 15.0) % 360.0
        assert abs(_wrap180(gmst_from_sidtime - _gmst06_deg(jd))) < 1e-6, jd


# --------------------------------------------------------------------------
# Fix 1: azalt/azalt_rev use the Vondrák mean obliquity
# --------------------------------------------------------------------------
def test_azalt_obliquity_is_vondrak_not_iau2006_at_deep_time():
    """Deep-time: azalt's mean obliquity is Vondrák, several arcsec off obl06."""
    jd_tt = _jd_tt(_jd_deep_bce())
    d_arcsec = (vondrak_mean_obliquity_deg(jd_tt) - _obl06_deg(jd_tt)) * 3600.0
    # Measured ~12.2" at -3000 (Vondrák pole-angle obliquity vs the obl06
    # extrapolation); the pre-fix code would have used the obl06 value here.
    assert abs(d_arcsec) > 3.0


def test_azalt_ecl_and_equ_paths_agree_at_3000bce():
    """ECL2HOR and EQU2HOR agree when both come from the library's own frame.

    The ecliptic-of-date and equatorial-of-date coordinates are related by the
    library's true obliquity (Vondrák mean + nutation). Feeding each through
    azalt yields the same horizon coordinates only if azalt shares that
    obliquity — which it now does. With the pre-fix obl06 obliquity the ECL2HOR
    path would sit ~0.6" off in RA at this epoch.
    """
    jd = _jd_deep_bce()
    eps_true = get_true_obliquity(_jd_tt(jd))
    ra, dec, _ = cotrans(_ECL, -eps_true)
    az_ecl, alt_ecl, _ = azalt(jd, ECL2HOR, _GEOPOS, 1013.25, 15.0, _ECL)
    az_equ, alt_equ, _ = azalt(jd, EQU2HOR, _GEOPOS, 1013.25, 15.0, (ra, dec, 1.0))
    assert abs(_wrap180(az_ecl - az_equ)) * 3600.0 < 1e-3
    assert abs(alt_ecl - alt_equ) * 3600.0 < 1e-3


def test_azalt_roundtrip_self_consistent_at_3000bce():
    """ECL2HOR then HOR2ECL recovers the input ecliptic coordinates."""
    jd = _jd_deep_bce()
    az, alt_true, _ = azalt(jd, ECL2HOR, _GEOPOS, 1013.25, 15.0, _ECL)
    lon2, lat2 = azalt_rev(jd, HOR2ECL, _GEOPOS, az, alt_true)
    assert abs(_wrap180(lon2 - _ECL[0])) * 3600.0 < 1e-3
    assert abs(lat2 - _ECL[1]) * 3600.0 < 1e-3


def _azalt_iau2006_baseline(jd_ut1, jd_tt, ecl, geopos):
    """Pre-fix azalt (true altitude, azimuth) via IAU 2006 obliquity + GMST."""
    lon, lat, _ = geopos
    eps0 = _obl06_deg(jd_tt)
    dpsi_rad, deps_rad = erfa.nut06a(_J2000, jd_tt - _J2000)
    eps = eps0 + math.degrees(deps_rad)
    dpsi_deg = math.degrees(dpsi_rad)
    ra, dec, _ = cotrans(ecl, -eps)
    eoe_deg = dpsi_deg * math.cos(math.radians(eps))
    lst_deg = (_gmst06_deg(jd_ut1) + eoe_deg + lon) % 360.0
    ha = math.radians((lst_deg - ra) % 360.0)
    dec_r = math.radians(dec)
    lat_r = math.radians(lat)
    sin_alt = math.sin(lat_r) * math.sin(dec_r) + math.cos(lat_r) * math.cos(
        dec_r
    ) * math.cos(ha)
    sin_alt = max(-1.0, min(1.0, sin_alt))
    alt = math.degrees(math.asin(sin_alt))
    y = math.sin(ha) * math.cos(dec_r)
    x = math.cos(ha) * math.cos(dec_r) * math.sin(lat_r) - math.sin(dec_r) * math.cos(
        lat_r
    )
    az = math.degrees(math.atan2(y, x)) % 360.0
    return az, alt


def test_azalt_modern_unchanged():
    """At a modern date azalt matches a from-scratch IAU 2006 (pre-fix) baseline.

    Near J2000 Vondrák ≡ IAU 2006 to sub-mas, so the fix must not move the
    modern result. The baseline reproduces the pre-fix algorithm (obl06 mean
    obliquity + IAU 2006 GMST) for the true altitude and azimuth.
    """
    jd = _jd_modern()
    az, alt_true, _ = azalt(jd, ECL2HOR, _GEOPOS, 1013.25, 15.0, _ECL)
    az_ref, alt_ref = _azalt_iau2006_baseline(jd, _jd_tt(jd), _ECL, _GEOPOS)
    assert abs(_wrap180(az - az_ref)) < 1e-6
    assert abs(alt_true - alt_ref) < 1e-6


# --------------------------------------------------------------------------
# Full calc_ut ecliptic frame at deep time (needs the extended DE441
# ephemeris; skipped when it is not available).
# --------------------------------------------------------------------------
def test_azalt_matches_calc_ut_frame_at_3000bce_extended():
    """azalt's ecliptic path matches calc_ut's own ecliptic/equatorial frame."""
    jd = _jd_deep_bce()
    prev_tier = get_precision_tier()
    try:
        set_precision_tier("extended")
        ecl, _ = calc_ut(jd, SUN, FLG_SWIEPH)
        equ, _ = calc_ut(jd, SUN, FLG_SWIEPH | FLG_EQUATORIAL)
    except (EphemerisRangeError, FileNotFoundError, OSError) as exc:
        pytest.skip(f"extended DE441 ephemeris unavailable: {exc}")
    finally:
        set_precision_tier(prev_tier)
    az_ecl, alt_ecl, _ = azalt(
        jd, ECL2HOR, _GEOPOS, 1013.25, 15.0, (ecl[0], ecl[1], ecl[2])
    )
    az_equ, alt_equ, _ = azalt(
        jd, EQU2HOR, _GEOPOS, 1013.25, 15.0, (equ[0], equ[1], equ[2])
    )
    assert abs(_wrap180(az_ecl - az_equ)) * 3600.0 < 1e-2
    assert abs(alt_ecl - alt_equ) * 3600.0 < 1e-2


# --------------------------------------------------------------------------
# The Vondrák model is evaluated by ERFA: poles, obliquity, precession matrix
# --------------------------------------------------------------------------
#: The two epochs where the sidereal time is tied to the IAU 2006 expression,
#: 1850 Jan 1.0 and 2050 Jan 1.0, written out rather than imported so that the
#: tests below depend on behaviour and not on the module's internals.
_ADOPTION_FIRST = 2396758.5
_ADOPTION_LAST = 2469807.5

#: Julian Dates (TT) from the deep past to the far future, both ends of the
#: adoption interval included.
_ERFA_GRID = tuple(
    _J2000 + (year - 2000) * 365.25
    for year in (-13000, -8000, -3000, -1000, 0, 1000, 1850, 2000, 2050, 2500, 5000)
) + (_ADOPTION_FIRST, _ADOPTION_LAST, 2451545.0 + 0.25)


def _epj(jd_tt: float) -> float:
    """Julian epoch, written out here rather than imported from the module."""
    return 2000.0 + (jd_tt - _J2000) / 365.25


def _rot_x(v: np.ndarray, angle: float) -> np.ndarray:
    """Rotate ``v`` about the x-axis by ``angle`` (radians), the module's way."""
    c, s = math.cos(angle), math.sin(angle)
    return np.array([v[0], v[1] * c + v[2] * s, -v[1] * s + v[2] * c])


@pytest.mark.parametrize("jd_tt", _ERFA_GRID)
def test_mean_obliquity_is_the_angle_between_the_erfa_poles(jd_tt):
    """``mean_obliquity_rad`` is acos(ltpecl · ltpequ) at the same epoch."""
    pecl = np.asarray(erfa.ltpecl(_epj(jd_tt)))
    pequ = np.asarray(erfa.ltpequ(_epj(jd_tt)))
    expected = math.acos(max(-1.0, min(1.0, float(np.dot(pecl, pequ)))))
    assert abs(sl.mean_obliquity_rad(jd_tt) - expected) < 1e-15
    assert sl.mean_obliquity_deg(jd_tt) == math.degrees(sl.mean_obliquity_rad(jd_tt))


def _gmst_iau2006_deg(jd_ut1: float) -> float:
    """The IAU 2006 Greenwich mean sidereal time, through ERFA, in degrees."""
    jd_tt = jd_ut1 + tu.deltat(jd_ut1)
    return math.degrees(erfa.gmst06(jd_ut1, 0.0, jd_tt, 0.0))


@pytest.mark.parametrize(
    "jd_ut1", (_ADOPTION_FIRST, 2415020.0, _J2000, 2440000.5, _ADOPTION_LAST)
)
def test_modern_era_is_the_iau_2006_expression(jd_ut1):
    """In the era it was fitted for, the answer is the IAU 2006 GMST.

    The two evaluations are of one expression, so what is left between them is
    the rounding of that expression, not a difference of model.
    """
    got = sl.mean_sidereal_time_deg(jd_ut1)
    assert abs(_wrap180(got - _gmst_iau2006_deg(jd_ut1))) * 3600.0 < 1e-4


@pytest.mark.parametrize(
    ("jd_ut1", "floor_deg"),
    (
        (_J2000 - 4000 * 365.25, 0.05),
        (_J2000 - 8000 * 365.25, 1.0),
        (_J2000 + 8000 * 365.25, 1.0),
    ),
)
def test_remote_epochs_leave_the_fitted_polynomial(jd_ut1, floor_deg):
    """Outside its fit interval the polynomial diverges and is not followed.

    The precession in right ascension is a truncated series around J2000.0; a
    sidereal time that stayed on it would be wrong by degrees at the ends of
    the span, and every house cusp with it.
    """
    departure = abs(
        _wrap180(sl.mean_sidereal_time_deg(jd_ut1) - _gmst_iau2006_deg(jd_ut1))
    )
    assert departure > floor_deg


@pytest.mark.parametrize("epoch", (_ADOPTION_FIRST, _ADOPTION_LAST))
def test_the_curve_has_no_step_at_the_adoption_epochs(epoch):
    """One curve, not two: the sidereal time does not jump where the
    description of it changes.

    Sampled at half-day steps on the far side of an adoption epoch, the value
    detrended by the sidereal rate lies on a straight line; extrapolating that
    line across the epoch reproduces the value there. A step would show up as
    a constant offset between the extrapolation and the value, and the
    curvature of the real curve over one day is orders of magnitude smaller
    than the steps this rules out.
    """
    rate = 360.98564736629  # degrees of sidereal time per day, near enough to detrend

    def detrended(jd: float) -> float:
        return _wrap180(sl.mean_sidereal_time_deg(jd) - rate * (jd - epoch))

    outside = -1.0 if epoch == _ADOPTION_FIRST else 1.0
    a = detrended(epoch + 3.0 * outside)
    b = detrended(epoch + 2.0 * outside)
    c = detrended(epoch + 1.0 * outside)
    quadratic = 3.0 * c - 3.0 * b + a  # extrapolated one step further
    assert abs(_wrap180(quadratic - detrended(epoch))) * 3600.0 < 1e-3


def test_the_sidereal_time_advances_over_the_whole_span():
    """Every step forward in time is a step forward in sidereal time.

    Walked from the first day of the DE441 kernel to its last, the value
    advances by the sidereal rotation of each interval: never a reversal, and a
    rate that stays within a few arcseconds a day of its J2000.0 value even at
    the ends, where the drift of ΔT and the change of the precession rate are
    real and large.
    """
    step = 40.0
    previous = None
    worst = 0.0
    jd = _J2000 - 13000 * 365.25
    last = _J2000 + 15000 * 365.25
    while jd <= last:
        value = sl.mean_sidereal_time_deg(jd)
        if previous is not None:
            advance = (value - previous) % 360.0
            expected = (360.98564736629 * step) % 360.0
            worst = max(worst, abs(_wrap180(advance - expected)) / step)
        previous = value
        jd += step
    assert worst * 3600.0 < 10.0


def test_the_apparent_seam_is_the_mean_one_plus_the_equation_of_the_equinoxes():
    """Nothing but the equation of the equinoxes and the longitude is added."""
    jd = 2451545.0
    dpsi, eps = 0.0047, 23.442600000000002
    mean = sl.mean_sidereal_time_deg(jd)
    for longitude in (-180.0, -37.5, 0.0, 37.5, 180.0):
        expected = (mean + dpsi * math.cos(math.radians(eps)) + longitude) % 360.0
        got = sl.apparent_sidereal_time_deg(
            jd, longitude, dpsi_deg=dpsi, eps_true_deg=eps
        )
        assert got == expected


def test_every_finite_epoch_is_answered():
    """No supported range and no coverage refusal: the unit is a model."""
    for jd in (-1e9, 0.0, 1e9, 1e12):
        value = sl.mean_sidereal_time_deg(jd)
        assert 0.0 <= value < 360.0


def test_infinite_epoch_is_refused_and_nan_propagates():
    """The domain contract of the ERFA-backed evaluators.

    ``math.sin(inf)`` used to raise ``ValueError`` out of the series; ERFA
    would answer NaN. The module keeps the refusal explicit. NaN, which the
    series propagated, still propagates.
    """
    for bad in (math.inf, -math.inf):
        with pytest.raises(ValueError):
            sl.mean_obliquity_rad(bad)
        with pytest.raises(ValueError):
            sl.mean_sidereal_time_deg(bad)
        with pytest.raises(ValueError):
            vondrak_mean_obliquity_deg(bad)
    assert math.isnan(sl.mean_obliquity_rad(math.nan))
    assert math.isnan(sl.mean_obliquity_deg(math.nan))
    assert math.isnan(sl.mean_sidereal_time_deg(math.nan))
    assert math.isnan(
        sl.apparent_sidereal_time_deg(math.nan, dpsi_deg=0.0, eps_true_deg=23.44)
    )
    for bad in (math.inf, -math.inf):
        with pytest.raises(ValueError):
            sl.apparent_sidereal_time_deg(bad, dpsi_deg=0.0, eps_true_deg=23.44)
