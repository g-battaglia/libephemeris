#!/usr/bin/env python3
"""Generate libephemeris/lunar_apse_corrections.py (interpolated-apse residuals).

Reconstructs the generator for the live apse correction tables (the
original lived in a task scratch directory and was not preserved).  The
runtime pipeline in ``lunar.py`` is::

    apogee  = calc_mean_lilith(jd) + _calc_elp2000_apogee_perturbations(jd)
              + apogee_correction(jd)
    perigee = calc_mean_lilith(jd) + 180 + _calc_elp2000_perigee_perturbations(jd)
              + perigee_correction(jd)

so each table entry is the residual against the reference ephemeris::

    correction(jd) = wrap180(reference_lon(jd) - analytic_lon(jd))   [arcsec]

with pyswisseph (Swiss Ephemeris, bodies INTP_APOG=21 / INTP_PERG=22,
TT scale via ``swe.calc``) as the black-box oracle, sampled on a fixed
JD grid (apogee: 10-day step; perigee: 2-day step; both spanning the
DE440 medium range ~1549-2651 CE).

After any recalibration of the trigonometric series
(``leph calibrate perigee`` -> paste coefficients), regenerate the
residual tables with this script, then run ``leph test lunar``.

Provenance / calibration disclosure:
    pyswisseph is used strictly as a black-box oracle. The stored values
    are numeric program *output* (computed apse positions), not Swiss
    Ephemeris source expression. INTP_APOG / INTP_PERG are constructs
    defined by the reference API, so 1:1 behavioral parity requires
    fitting to reference output; the analytic series and this fitting
    pipeline are original. See NOTICE.md, "Calibration Data Disclosure".
    The generated module carries the project dual-license header (the
    tables are project-owned output data).

Usage:
    python scripts/generate_lunar_apse_corrections.py --verify 2000
        Sample N random grid points, recompute the residuals, and report
        agreement with the stored tables (no file written).

    python scripts/generate_lunar_apse_corrections.py --write
        Recompute both tables on the full grid and rewrite
        libephemeris/lunar_apse_corrections.py (~10 min: 241k oracle
        calls).
"""

from __future__ import annotations

import argparse
import datetime as _dt
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _wrap180(deg: float) -> float:
    return (deg + 180.0) % 360.0 - 180.0


def _analytic_apogee(jd_tt: float) -> float:
    from libephemeris.lunar import (
        _calc_elp2000_apogee_perturbations,
        calc_mean_lilith,
    )

    return calc_mean_lilith(jd_tt) + _calc_elp2000_apogee_perturbations(jd_tt)


def _analytic_perigee(jd_tt: float) -> float:
    from libephemeris.lunar import (
        _calc_elp2000_perigee_perturbations,
        calc_mean_lilith,
    )

    return (
        (calc_mean_lilith(jd_tt) + 180.0) % 360.0
        + _calc_elp2000_perigee_perturbations(jd_tt)
    )


def _residual_arcsec(swe, body: int, jd_tt: float, analytic_fn) -> float:
    ref_lon = swe.calc(jd_tt, body, swe.FLG_SWIEPH)[0][0]
    return _wrap180(ref_lon - analytic_fn(jd_tt) % 360.0) * 3600.0


def _grids():
    from libephemeris.lunar_apse_corrections import (
        APOGEE_CT_COUNT,
        APOGEE_CT_JD_START,
        APOGEE_CT_STEP_DAYS,
        PERIGEE_CT_COUNT,
        PERIGEE_CT_JD_START,
        PERIGEE_CT_STEP_DAYS,
    )

    return (
        ("apogee", 21, APOGEE_CT_JD_START, APOGEE_CT_STEP_DAYS, APOGEE_CT_COUNT),
        ("perigee", 22, PERIGEE_CT_JD_START, PERIGEE_CT_STEP_DAYS, PERIGEE_CT_COUNT),
    )


def verify(n_samples: int) -> int:
    import swisseph as swe

    from libephemeris.lunar_apse_corrections import (
        APOGEE_CORRECTIONS,
        PERIGEE_CORRECTIONS,
    )

    swe.set_ephe_path("data")
    rng = random.Random(20260612)
    tables = {"apogee": APOGEE_CORRECTIONS, "perigee": PERIGEE_CORRECTIONS}
    fns = {"apogee": _analytic_apogee, "perigee": _analytic_perigee}
    worst_overall = 0.0
    for name, body, jd0, step, count in _grids():
        stored = tables[name]
        diffs = []
        for _ in range(n_samples):
            i = rng.randrange(count)
            jd = jd0 + i * step
            recomputed = _residual_arcsec(swe, body, jd, fns[name])
            diffs.append(abs(recomputed - float(stored[i])))
        diffs.sort()
        print(
            f"{name}: n={n_samples} median={diffs[len(diffs) // 2]:.3f}\" "
            f"p95={diffs[int(0.95 * len(diffs))]:.3f}\" max={diffs[-1]:.3f}\""
        )
        worst_overall = max(worst_overall, diffs[-1])
    # The stored values are rounded to 0.01"; anything beyond ~0.5" means
    # the trig series changed since generation -> regenerate with --write.
    return 0 if worst_overall < 0.5 else 1


def write_module() -> None:
    import swisseph as swe

    swe.set_ephe_path("data")
    fns = {"apogee": _analytic_apogee, "perigee": _analytic_perigee}
    results = {}
    for name, body, jd0, step, count in _grids():
        vals = []
        for i in range(count):
            jd = jd0 + i * step
            vals.append(_residual_arcsec(swe, body, jd, fns[name]))
            if i % 20000 == 0:
                print(f"{name}: {i}/{count}", flush=True)
        results[name] = (jd0, step, count, vals)

    out = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "libephemeris",
        "lunar_apse_corrections.py",
    )
    (ajd0, astep, acount, avals) = results["apogee"]
    (pjd0, pstep, pcount, pvals) = results["perigee"]
    with open(out, "w") as f:
        f.write(
            f'''# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
AUTO-GENERATED FILE - DO NOT EDIT MANUALLY

Correction tables for interpolated lunar apogee and perigee longitude.
Residual corrections to be added AFTER the trigonometric perturbation series.

Generated by: scripts/generate_lunar_apse_corrections.py
Generated on: {_dt.datetime.now():%Y-%m-%d %H:%M:%S}
Reference: pyswisseph (Swiss Ephemeris) as black-box oracle
Date range: ~1549 to ~2651 CE (DE440 medium ephemeris range)

Method:
  1. Compute mean apogee/perigee from analytical polynomial (calc_mean_lilith)
  2. Compute true perturbation = pyswisseph_position - mean_position
  3. Fit Delaunay-based trigonometric series via least squares
  4. Store residuals (true_perturbation - trig_fit) as correction tables
  5. At runtime: longitude = mean + trig_series + interpolated_correction

Provenance / calibration disclosure:
  The stored values are numeric program *output* — computed positions of
  the interpolated lunar apsides — obtained with pyswisseph used strictly
  as a black-box oracle. They are not Swiss Ephemeris source expression.
  The interpolated-apse bodies (INTP_APOG / INTP_PERG) are constructs
  defined by the reference API, so 1:1 behavioral parity requires fitting
  to reference output; the analytic series and this fitting pipeline are
  original. Disclosed in NOTICE.md, section "Calibration Data Disclosure".
"""

from __future__ import annotations

# =============================================================================
# APOGEE CORRECTION TABLE METADATA
# =============================================================================
APOGEE_CT_JD_START: float = {ajd0!r}
APOGEE_CT_STEP_DAYS: float = {astep!r}
APOGEE_CT_COUNT: int = {acount}

# =============================================================================
# PERIGEE CORRECTION TABLE METADATA
# =============================================================================
PERIGEE_CT_JD_START: float = {pjd0!r}
PERIGEE_CT_STEP_DAYS: float = {pstep!r}
PERIGEE_CT_COUNT: int = {pcount}

# =============================================================================
# APOGEE CORRECTIONS (arcseconds)
# JD {ajd0:.1f} to {ajd0 + (acount - 1) * astep:.1f}
# =============================================================================
APOGEE_CORRECTIONS: tuple[float, ...] = (
'''
        )
        for v in avals:
            f.write(f"    {v:.2f},\n")
        f.write(
            f""")

# =============================================================================
# PERIGEE CORRECTIONS (arcseconds)
# JD {pjd0:.1f} to {pjd0 + (pcount - 1) * pstep:.1f}
# =============================================================================
PERIGEE_CORRECTIONS: tuple[float, ...] = (
"""
        )
        for v in pvals:
            f.write(f"    {v:.2f},\n")
        f.write(")\n")
    print(f"wrote {out}")


def main() -> int:
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--verify", type=int, metavar="N", help="sample-check N points")
    g.add_argument("--write", action="store_true", help="regenerate the module")
    args = ap.parse_args()
    if args.verify:
        return verify(args.verify)
    write_module()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
