"""Audit the untested flag-combo paths of calc_pctr / nod_aps vs pyswisseph.

Standalone validation script (not part of the pytest suite). The planet-centric
(`calc_pctr`) and node/apside (`nod_aps`) paths have essentially no compare-test
coverage, which is exactly where issue #34 (FLG_NONUT ignored, sidereal+equatorial
speed frame) hid. This sweeps both functions against the pyswisseph oracle across
the relevant flag matrix and reports every divergence beyond tolerance, so any
remaining #34-class gap is surfaced before v3.

Everything is compared in TT (both libephemeris and pyswisseph take ET/TT here),
so ΔT-model differences never enter. pyswisseph is a black-box oracle only.

Usage:
    .venv/bin/python scripts/validate_nonut_nodaps.py
"""

from __future__ import annotations

import os
import sys

try:
    import swisseph as swe
except ImportError as exc:  # pragma: no cover - dev tooling
    print(f"Missing dependency: {exc}")
    sys.exit(1)

import libephemeris as ephem
from libephemeris import constants as C

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SWISS_DIR = os.path.join(_REPO, "data", "reference")
_LEB = os.path.join(_REPO, "data", "leb", "ephemeris_medium.leb")

# Report anything above this (arcsec) as a candidate bug; frame/nutation bugs are
# arcsec-to-degrees, while the calc_pctr/nod_aps pipelines themselves agree with
# Swiss to well under 1" on the shared paths.
TOL_ARCSEC = 1.0

# Flag matrix (name -> value). Single flags + the 2-combos most likely to expose
# a frame/nutation gap. FLG_SWIEPH is added on the Swiss side only.
_F = C
FLAG_SETS = {
    "base": 0,
    "NONUT": _F.FLG_NONUT,
    "SPEED": _F.FLG_SPEED,
    "EQUATORIAL": _F.FLG_EQUATORIAL,
    "J2000": _F.FLG_J2000,
    "TRUEPOS": _F.FLG_TRUEPOS,
    "NOABERR": _F.FLG_NOABERR,
    "SIDEREAL": _F.FLG_SIDEREAL,
    "NONUT|EQ": _F.FLG_NONUT | _F.FLG_EQUATORIAL,
    "NONUT|SPEED": _F.FLG_NONUT | _F.FLG_SPEED,
    "SID|EQ": _F.FLG_SIDEREAL | _F.FLG_EQUATORIAL,
    "SID|SPEED": _F.FLG_SIDEREAL | _F.FLG_SPEED,
    "SID|EQ|SPEED": _F.FLG_SIDEREAL | _F.FLG_EQUATORIAL | _F.FLG_SPEED,
    "EQ|SPEED": _F.FLG_EQUATORIAL | _F.FLG_SPEED,
    "J2000|SPEED": _F.FLG_J2000 | _F.FLG_SPEED,
    "TRUEPOS|EQ": _F.FLG_TRUEPOS | _F.FLG_EQUATORIAL,
}

# calc_pctr target/center pairs.
PCTR_PAIRS = [
    (C.MARS, C.SUN, "Mars/Sun"),
    (C.JUPITER, C.SUN, "Jup/Sun"),
    (C.MOON, C.EARTH, "Moon/Earth"),
    (C.VENUS, C.EARTH, "Venus/Earth"),
    (C.MERCURY, C.SUN, "Merc/Sun"),
]

# nod_aps bodies (planetary nodes/apsides; Sun/Earth have none).
NODAPS_BODIES = [
    (C.MERCURY, "Mercury"),
    (C.VENUS, "Venus"),
    (C.MARS, "Mars"),
    (C.JUPITER, "Jupiter"),
    (C.SATURN, "Saturn"),
    (C.PLUTO, "Pluto"),
]

JDS = [swe.julday(y, 1, 1, 12.0) for y in (1600, 1800, 1900, 2000, 2100, 2300, 2500)]


def _ang(a: float, b: float) -> float:
    d = abs(a - b) % 360.0
    if d > 180.0:
        d = 360.0 - d
    return d * 3600.0


def _setup() -> None:
    swe.set_ephe_path(_SWISS_DIR)
    swe.set_sid_mode(swe.SIDM_FAGAN_BRADLEY, 0.0, 0.0)
    ephem.set_leb_file(_LEB)
    ephem.set_precision_tier("medium")
    ephem.set_calc_mode("leb")
    ephem.set_sid_mode(C.SIDM_FAGAN_BRADLEY, 0.0, 0.0)


def _libe_swe_pctr(jd, tgt, ctr, fval):
    lib, _ = ephem.calc_pctr(jd, tgt, ctr, fval)
    sw, _ = swe.calc_pctr(jd, tgt, ctr, swe.FLG_SWIEPH | fval)
    return max(_ang(lib[0], sw[0]), abs(lib[1] - sw[1]) * 3600.0)


def _libe_swe_nodaps(jd, body, fval):
    # Compare only the NODES (ascending=0, descending=1). The apsides (peri=2,
    # aphe=3) of low-eccentricity orbits are an ill-conditioned osculating
    # quantity that diverges 100-2000" between any two element derivations
    # (documented in divergences.md §7) -- not a meaningful flag-handling signal.
    lib = ephem.nod_aps(jd, body, fval, C.NODBIT_OSCU)
    sw = swe.nod_aps(jd, body, swe.FLG_SWIEPH | fval, swe.NODBIT_OSCU)
    return max(
        max(_ang(lib[i][0], sw[i][0]), abs(lib[i][1] - sw[i][1]) * 3600.0)
        for i in range(2)
    )


def _audit(label, cases, delta_fn) -> list[str]:
    """cases: list of (key, args-tuple). delta_fn(jd_args..., fval) -> arcsec.

    Reports, per flag, the EXCESS over base: max over cases of
    (delta(flag) - delta(base)). A small excess means the flag is handled the
    same as Swiss (no #34-class bug); the pre-existing base disagreement is
    reported separately so it is not mistaken for a flag bug.
    """
    print("=" * 78)
    print(f"{label} vs pyswisseph (TT) — flag matrix")
    print("=" * 78)
    hits: list[str] = []

    base: dict = {}
    worst_base = 0.0
    worst_base_label = ""
    for key, args in cases:
        for jd in JDS:
            try:
                d = delta_fn(jd, *args, 0)
            except Exception as exc:  # noqa: BLE001
                hits.append(f"  [ERR] {label} {key} base @JD{jd:.0f}: {exc}")
                d = None
            base[(key, jd)] = d
            if d is not None and d > worst_base:
                worst_base, worst_base_label = d, f"{key}@{jd:.0f}"
    print(
        f'  {"base (pipeline)":14} max|Δang|={worst_base:10.4f}"  ({worst_base_label})'
        f"  [pre-existing, not flag-related]"
    )

    for fname, fval in FLAG_SETS.items():
        if fval == 0:
            continue
        worst_excess = -1e9
        worst_label = ""
        for key, args in cases:
            for jd in JDS:
                b = base.get((key, jd))
                if b is None:
                    continue
                try:
                    d = delta_fn(jd, *args, fval)
                except Exception as exc:  # noqa: BLE001
                    hits.append(f"  [ERR] {label} {key} {fname} @JD{jd:.0f}: {exc}")
                    continue
                excess = d - b
                if excess > worst_excess:
                    worst_excess, worst_label = excess, f"{key}@{jd:.0f}"
        flag = "  <-- FLAG BUG?" if worst_excess > TOL_ARCSEC else ""
        print(f'  {fname:14} excess={worst_excess:10.4f}"  ({worst_label}){flag}')
        if worst_excess > TOL_ARCSEC:
            hits.append(
                f'  {label} {fname}: flag adds {worst_excess:.4f}" vs base at {worst_label}'
            )
    return hits


def audit_calc_pctr() -> list[str]:
    cases = [((p[2]), (p[0], p[1])) for p in PCTR_PAIRS]
    return _audit("calc_pctr", cases, _libe_swe_pctr)


def audit_nod_aps() -> list[str]:
    cases = [((b[1]), (b[0],)) for b in NODAPS_BODIES]
    print()
    print("(nod_aps: nodes only; apsides of low-e orbits are documented §7)")
    return _audit("nod_aps (nodes)", cases, _libe_swe_nodaps)


def main() -> int:
    _setup()
    hits = audit_calc_pctr() + audit_nod_aps()
    print()
    print("=" * 78)
    if hits:
        print(f"CANDIDATE DIVERGENCES ({len(hits)}) — review vs divergences.md:")
        for h in hits:
            print(h)
    else:
        print("No divergences above tolerance. calc_pctr / nod_aps match Swiss.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
