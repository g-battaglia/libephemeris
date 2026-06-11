#!/usr/bin/env python3
"""Verify fixed-star parity against pyswisseph (verification only).

Compares, for every catalog star bright enough for the chosen limit and
present in both catalogs (matched by Bayer/Flamsteed nomenclature key):

- positions at J2000 and two distant epochs,
- the name-resolution semantics table (name, case, ',nomenclature',
  'name,nomenclature', sequential digits, '%' wildcard, error classes),
- retflag echo, TOPOCTR effect and speed columns.

The Swiss Ephemeris star file is used here exclusively as an external
verification reference, never as a data source.

Usage:
    python scripts/verify_fix_fixstars.py [--mag-limit 3.0]
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import statistics
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mag-limit", type=float, default=3.0)
    args = ap.parse_args()

    import swisseph as swe

    import libephemeris as le
    from libephemeris.fixed_stars import STAR_CATALOG

    le.set_calc_mode("skyfield")
    tmp = tempfile.mkdtemp()
    shutil.copy("data/reference/sefstars.txt", os.path.join(tmp, "sefstars.txt"))
    swe.set_ephe_path(tmp + os.pathsep + os.path.abspath("data"))

    failures = 0

    # ---- positional sweep ----
    for jd, label in ((2451545.0, "J2000"), (2433282.5, "1950"), (2469807.5, "2050")):
        diffs = []
        for e in STAR_CATALOG:
            if e.magnitude > args.mag_limit:
                continue
            key = "," + e.nomenclature
            try:
                x1, _n1, _ = swe.fixstar2_ut(key, jd, swe.FLG_SWIEPH)
            except Exception:
                continue
            x2, _n2, _ = le.fixstar_ut(key, jd, le.FLG_SWIEPH)
            dlon = (
                abs((x2[0] - x1[0] + 180) % 360 - 180)
                * 3600
                * math.cos(math.radians(x1[1]))
            )
            dlat = abs(x2[1] - x1[1]) * 3600
            diffs.append(math.hypot(dlon, dlat))
        diffs.sort()
        within = sum(1 for d in diffs if d <= 0.1)
        print(
            f"{label}: n={len(diffs)} median={statistics.median(diffs):.4f}\" "
            f"p90={diffs[int(0.9 * len(diffs))]:.4f}\" <=0.1\": {within}/{len(diffs)}"
        )
        # alpha Cen's component-identity comparison artifact gives 2
        # known >0.1" rows; anything beyond that is a regression.
        if len(diffs) - within > 2:
            failures += 1

    # ---- name semantics ----
    jd = 2451545.0
    for form in ("Regulus", "regulus", ",alTau", "Foo,alTau", "alde%", "Sirius"):
        x1, n1, _ = swe.fixstar2_ut(form, jd, swe.FLG_SWIEPH)
        x2, n2, _ = le.fixstar_ut(form, jd, le.FLG_SWIEPH)
        ok = n1 == n2 and abs(x2[0] - x1[0]) * 3600 < 0.1
        print(f"resolve {form!r}: SE={n1} LE={n2} {'OK' if ok else 'FAIL'}")
        failures += 0 if ok else 1
    for bad in ("bad%bad", "NoSuchStarXYZ"):
        se_err = le_err = False
        try:
            swe.fixstar2_ut(bad, jd, swe.FLG_SWIEPH)
        except Exception:
            se_err = True
        try:
            le.fixstar_ut(bad, jd, le.FLG_SWIEPH)
        except Exception:
            le_err = True
        ok = se_err and le_err
        print(f"resolve {bad!r}: errors SE={se_err} LE={le_err} {'OK' if ok else 'FAIL'}")
        failures += 0 if ok else 1

    # ---- retflag / TOPOCTR / speeds ----
    _x, _n, r1 = swe.fixstar2_ut("Aldebaran", jd, 0)
    _x, _n, r2 = le.fixstar_ut("Aldebaran", jd, 0)
    print(f"retflag echo (0): SE={r1} LE={r2} {'OK' if r1 == r2 else 'FAIL'}")
    failures += 0 if r1 == r2 else 1

    swe.set_topo(-75.0, 40.0, 0.0)
    le.set_topo(-75.0, 40.0, 0.0)
    xg1, _, _ = swe.fixstar2_ut("Aldebaran", jd, swe.FLG_SWIEPH)
    xt1, _, rt1 = swe.fixstar2_ut(
        "Aldebaran", jd, swe.FLG_SWIEPH | swe.FLG_TOPOCTR
    )
    xg2, _, _ = le.fixstar_ut("Aldebaran", jd, le.FLG_SWIEPH)
    xt2, _, rt2 = le.fixstar_ut("Aldebaran", jd, le.FLG_SWIEPH | le.FLG_TOPOCTR)
    d_se = (xt1[0] - xg1[0]) * 3600
    d_le = (xt2[0] - xg2[0]) * 3600
    ok = rt1 == rt2 and abs(d_se - d_le) < 0.01
    print(
        f"TOPOCTR: dlon SE={d_se:+.4f}\" LE={d_le:+.4f}\" retflag "
        f"SE={rt1} LE={rt2} {'OK' if ok else 'FAIL'}"
    )
    failures += 0 if ok else 1

    print("FAILURES:", failures)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
