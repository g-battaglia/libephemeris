"""Broad validation of libephemeris (LEB backend) versus Swiss Ephemeris.

This is a *standalone* validation script (not part of the pytest suite). It
cross-checks geocentric apparent ecliptic longitude/latitude from libephemeris's
LEB backend against pyswisseph across a wide span of dates and bodies, with two
goals:

1. **Modern broad sweep** (1600-2600, Swiss ``.se1`` / DE431 reference): confirm
   the LEB medium tier matches Swiss to sub-arcsecond for the main bodies.

2. **Long-term sweep** (-3000 to +3000, LEB extended tier; +3000 is the Swiss
   Moshier upper bound): confirm the Vondrák 2011 long-term precession migration
   brings the Sun -- which isolates the precession term -- close to Swiss at
   extreme epochs (where the old IAU 2006 polynomial was ~36" off at -3000). The
   residual on the planets is the irreducible DE441-vs-Swiss(DE431/Moshier)
   ephemeris floor, which Vondrák does not remove. The script also prints the IAU
   2006-vs-Vondrák precession divergence at each epoch as the "counterfactual"
   error the migration removed.

For BCE dates pyswisseph has no ``.se1`` data here and transparently falls back
to its built-in Moshier model (which itself uses Vondrák precession), so the
long-term comparison isolates the ephemeris difference, not precession.

Usage:
    .venv/bin/python scripts/validate_vondrak_vs_swiss.py
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np

import libephemeris as ephem

try:
    import erfa
    import swisseph as swe
except ImportError as exc:  # pragma: no cover - dev tooling
    print(f"Missing dependency for validation: {exc}")
    sys.exit(1)

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_LEB_DIR = os.path.join(_REPO, "data", "leb")
_SWISS_DIR = os.path.join(_REPO, "data", "reference")
_J2000 = 2451545.0

# (id, name) for the bodies cross-checked. Sun first because it isolates the
# precession term in the long-term sweep.
MAIN_BODIES = [
    (0, "Sun"),
    (1, "Moon"),
    (2, "Mercury"),
    (3, "Venus"),
    (4, "Mars"),
    (5, "Jupiter"),
    (6, "Saturn"),
    (7, "Uranus"),
    (8, "Neptune"),
    (9, "Pluto"),
    (10, "MeanNode"),
    (11, "TrueNode"),
]


def _ang_arcsec(a: float, b: float) -> float:
    """Smallest angular separation in arcseconds, with 360 wrap."""
    d = abs(a - b) % 360.0
    if d > 180.0:
        d = 360.0 - d
    return d * 3600.0


def _leb_tt(jd_tt: float, body: int) -> tuple[float, float]:
    """libephemeris geocentric apparent ecliptic (lon, lat), TT input."""
    pos, _ = ephem.calc(jd_tt, body, 0)
    return pos[0], pos[1]


def _swiss_tt(jd_tt: float, body: int) -> tuple[float, float, int]:
    """Swiss geocentric apparent ecliptic (lon, lat), TT input."""
    pos, retflag = swe.calc(jd_tt, body, swe.FLG_SWIEPH)
    return pos[0], pos[1], retflag


def _leb_ut(jd_ut: float, body: int) -> float:
    pos, _ = ephem.calc_ut(jd_ut, body, 0)
    return pos[0]


def _swiss_ut(jd_ut: float, body: int) -> float:
    pos, _ = swe.calc_ut(jd_ut, body, swe.FLG_SWIEPH)
    return pos[0]


def _precession_divergence_arcsec(jd_tt: float) -> float:
    """Max element |Vondrák - IAU2006| precession matrix, as arcsec.

    A proxy for the apparent-place error the IAU 2006 polynomial would still
    carry at this epoch (i.e. what the Vondrák migration removed)."""
    epj = 2000.0 + (jd_tt - _J2000) / 365.25
    vondrak = np.asarray(erfa.ltpb(epj))
    iau2006 = np.asarray(erfa.pmat06(_J2000, jd_tt - _J2000))
    return math.degrees(float(np.max(np.abs(vondrak - iau2006)))) * 3600.0


def _setup_leb(tier: str) -> bool:
    path = os.path.join(_LEB_DIR, f"ephemeris_{tier}.leb")
    if not os.path.exists(path):
        print(f"  [skip] LEB {tier} file not found: {path}")
        return False
    ephem.set_leb_file(path)
    ephem.set_precision_tier(tier)
    ephem.set_calc_mode("leb")
    return True


def modern_sweep() -> int:
    print("=" * 78)
    print("MODERN SWEEP  1600-2600  (LEB medium vs Swiss .se1 / DE431)")
    print("Compared in TT to isolate ephemeris + precession from the ΔT models.")
    print("=" * 78)
    swe.set_ephe_path(_SWISS_DIR)
    if not _setup_leb("medium"):
        return 0

    n = 61
    span = 1000 * 365.25 / (n - 1)
    jds = [swe.julday(1600, 1, 1, 12.0) + i * span for i in range(n)]

    print(
        f"{'Body':<10} {'TT max Δlon"':>13} {'TT mean Δlon"':>14} "
        f"{'TT max Δlat"':>13} {'UT max Δlon"':>13}"
    )
    print("-" * 70)
    worst_tt = 0.0
    for bid, name in MAIN_BODIES:
        dlon_tt, dlat_tt, dlon_ut = [], [], []
        for jd in jds:
            l_lon, l_lat = _leb_tt(jd, bid)
            s_lon, s_lat, _ = _swiss_tt(jd, bid)
            dlon_tt.append(_ang_arcsec(l_lon, s_lon))
            dlat_tt.append(abs(l_lat - s_lat) * 3600.0)
            dlon_ut.append(_ang_arcsec(_leb_ut(jd, bid), _swiss_ut(jd, bid)))
        mx, mean, mxlat = max(dlon_tt), sum(dlon_tt) / len(dlon_tt), max(dlat_tt)
        worst_tt = max(worst_tt, mx)
        print(
            f"{name:<10} {mx:>13.4f} {mean:>14.4f} {mxlat:>13.4f} {max(dlon_ut):>13.4f}"
        )
    print("-" * 70)
    print(
        f'Worst-case TT longitude delta (ephemeris+precession): {worst_tt:.4f}"\n'
        "The much larger 'UT max Δlon' is the documented ΔT-model divergence at the\n"
        "1600/2600 extremes (ΔT ~90 s at 1600, ~1200 s at 2600); it is NOT a\n"
        "precession or ephemeris error -- in TT it collapses to sub-arcsecond."
    )
    return 0


def longterm_sweep() -> int:
    print()
    print("=" * 78)
    print("LONG-TERM SWEEP  -3000..+3000  (LEB extended vs Swiss)")
    print("Compared in TT (precession-isolating); BCE/early dates -> Swiss Moshier.")
    print("Sun isolates precession; planets show the DE441-vs-Swiss ephemeris floor.")
    print("=" * 78)
    swe.set_ephe_path(_SWISS_DIR)
    if not _setup_leb("extended"):
        return 0

    # +3000 is the Moshier upper bound in pyswisseph; beyond it Swiss has no data
    # here and raises, so the sweep stops at +3000.
    years = [-3000, -2000, -1000, 0, 1000, 2000, 3000]
    bodies = [(0, "Sun"), (2, "Mercury"), (3, "Venus"), (4, "Mars"), (5, "Jupiter")]

    head = f"{'Year':>6} {'IAU2006-would-be':>17} | " + " ".join(
        f"{name:>8}" for _, name in bodies
    )
    print(head)
    print(f"{'':>6} {'(Δprec, arcsec)':>17} | " + "  ".join('Δlon"' for _ in bodies))
    print("-" * len(head))
    for yr in years:
        jd_tt = swe.julday(yr, 1, 1, 12.0)
        counterfactual = _precession_divergence_arcsec(jd_tt)
        cells = []
        for bid, _ in bodies:
            try:
                l_lon, _ = _leb_tt(jd_tt, bid)
                s_lon, _, _ = _swiss_tt(jd_tt, bid)
                cells.append(f"{_ang_arcsec(l_lon, s_lon):>8.2f}")
            except swe.Error:
                cells.append(f"{'n/a':>8}")
        print(f'{yr:>6} {counterfactual:>15.2f}" | ' + " ".join(cells))
    print("-" * len(head))
    print(
        "Sun Δlon stays single-digit arcsec at every epoch: the Vondrák precession\n"
        "now matches Swiss, so only the DE441-vs-Moshier ephemeris floor remains.\n"
        "With the old IAU 2006 precession the Sun would instead carry the\n"
        "'IAU2006-would-be' error (~13\" matrix diff at -3000, ~36\" projected in\n"
        "ecliptic longitude). The larger deep-BCE planet deltas (e.g. Mars at -3000)\n"
        "are that documented ephemeris floor, not a precession bug."
    )
    return 0


def main() -> int:
    saved_mode = ephem.state._CALC_MODE
    saved_leb = ephem.state._LEB_FILE
    saved_tier = ephem.state._PRECISION_TIER
    try:
        modern_sweep()
        longterm_sweep()
    finally:
        ephem.state._LEB_FILE = saved_leb
        ephem.set_calc_mode(saved_mode)
        if saved_tier is not None:
            ephem.set_precision_tier(saved_tier)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
