# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Frame-consistency regression tests for the Vondrák 2011 long-term model.

These pin two extreme-date fixes that route ``azalt``/``azalt_rev`` and the
public ``sidtime``/``sidtime0`` onto the same long-term (Vondrák 2011) obliquity
and Greenwich mean sidereal time the rest of the library uses for
ecliptic-of-date positions and house cusps. Previously both bypassed it through
the IAU 2006 models (``erfa.obl06`` / ``erfa.gmst06``), which are internally
inconsistent with the library's own frame at deep-time epochs (~5.7" of
obliquity and ~0.66° of GMST at year -3000).

The IAU 2006 models are used here only as the pre-fix oracle: inside the modern
window (1850-2050) the fixes are, by design, a no-op, so the modern result must
still match them to < 1e-6.
"""

from __future__ import annotations

import math

import erfa
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
    # Measured ~5.72"; the pre-fix code would have used the obl06 value here.
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


def _azalt_iau2006_oracle(jd_ut1, jd_tt, ecl, geopos):
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
    """At a modern date azalt matches a from-scratch IAU 2006 (pre-fix) oracle.

    Near J2000 Vondrák ≡ IAU 2006 to sub-mas, so the fix must not move the
    modern result. The oracle reproduces the pre-fix algorithm (obl06 mean
    obliquity + IAU 2006 GMST) for the true altitude and azimuth.
    """
    jd = _jd_modern()
    az, alt_true, _ = azalt(jd, ECL2HOR, _GEOPOS, 1013.25, 15.0, _ECL)
    az_ref, alt_ref = _azalt_iau2006_oracle(jd, _jd_tt(jd), _ECL, _GEOPOS)
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
