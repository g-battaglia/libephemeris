# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Verification harness for the Galilean satellite module.

Records the numeric behavior of ``galilean_moon_positions`` (and the Jupiter
COB offset derived from it through ``moon_theories.constants``) over a dense
JD grid, so that any change to
``libephemeris/moon_theories/galilean.py`` can be proven output-preserving.

Modes:
    --dump FILE        Record every output as JSON (run BEFORE the rewrite).
    --check FILE       Re-run the grid and fail (exit 1) if any moon component
                       drifted by more than --tol-km (default 1e-3 km) or any
                       COB component by more than --tol-cob-au (default 1e-12).
    compare            Informational physical-accuracy check vs a JPL Jupiter
                       satellite SPK (jup365.bsp) when one can be found;
                       skips gracefully otherwise.

Usage:
    .venv/bin/python scripts/verify_galilean.py --dump tasks/results/galilean_baseline.json
    .venv/bin/python scripts/verify_galilean.py --check tasks/results/galilean_baseline.json
    .venv/bin/python scripts/verify_galilean.py compare [--spk path/to/jup365.bsp]

Grid design:
    * 1800-01-01 .. ~2200-01 every 16 days (9132 epochs), wide enough that any
      frame/precession mistake that grows with |t - J2000| is amplified at the
      edges.
    * Golden epochs: the E5 theory epoch (JDE 2443000.5), J2000, and the
      fixed dates used by the unit/compare suites.
    * +/- 1 s pairs around two epochs to pin the central-difference velocity
      path used by ``planets._cob_velocity_correction``.

Gate rationale: the smallest published series terms are ~1e-5 deg
(~0.3 km on Callisto's orbit), while pure floating-point re-association in
an independently structured implementation stays at the centimeter level;
1e-3 km therefore detects any single mistranscribed coefficient without
tripping on legitimate FP reordering. Through the GM-weighted COB dilution
(~2.07e-4) and Jupiter's geocentric distance, even 1 km of moon drift moves
Jupiter's geocentric position by only ~7e-11 arcsec, so the gate is far
tighter than astronomical use requires: it is a transcription-correctness
detector, not an external-output bound.

This script is standalone (project convention: verification via standalone
scripts, not the pytest suite).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
MODULE_PATH = REPO_ROOT / "libephemeris" / "moon_theories" / "galilean.py"

MOON_NAMES = ("Io", "Europa", "Ganymede", "Callisto")

# 16-day grid: 1800-01-01 00:00 TT .. 2200-01-06 (9132 epochs).
GRID_START_JD = 2378496.5
GRID_STEP_DAYS = 16.0
GRID_COUNT = 9132

ONE_SECOND = 1.0 / 86400.0

GOLDEN_JDS = [
    2443000.5,  # E5 theory epoch (t = 0)
    2451545.0,  # J2000.0
    2455000.0,
    2460000.0,
    2415020.5,  # 1900-01-00.5 (B1900-era compare date)
    2433454.0,
    2460477.0,
    2460629.5,
    2460676.5,
    2488069.49995833,  # 2100-12-31 (compare-suite upper edge)
]

VELOCITY_PAIR_JDS = [
    2451545.0 - ONE_SECOND,
    2451545.0 + ONE_SECOND,
    2460000.0 - ONE_SECOND,
    2460000.0 + ONE_SECOND,
]


def iter_jds() -> list[float]:
    jds = [GRID_START_JD + GRID_STEP_DAYS * k for k in range(GRID_COUNT)]
    jds.extend(GOLDEN_JDS)
    jds.extend(VELOCITY_PAIR_JDS)
    return jds


def evaluate(jd: float, ts) -> list[float]:
    """Return the 15 recorded floats for one epoch.

    Twelve moon components (km, Jupiter-centric J2000 ecliptic, order
    Io/Europa/Ganymede/Callisto) followed by the three Jupiter COB offset
    components (AU, equatorial ICRF) — the latter pins the GM weighting and
    the ecliptic->equatorial rotation in ``moon_theories.constants`` too.
    """
    from libephemeris.moon_theories.constants import get_cob_offset
    from libephemeris.moon_theories.galilean import galilean_moon_positions

    values: list[float] = []
    for xyz in galilean_moon_positions(jd):
        values.extend(float(c) for c in xyz)
    values.extend(float(c) for c in get_cob_offset("jupiter barycenter", ts.tt_jd(jd)))
    return values


def _timescale():
    from skyfield.api import load

    return load.timescale()


def _git_head() -> str:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        return out.stdout.strip()
    except Exception:  # noqa: BLE001 - metadata only
        return "unknown"


def do_dump(path: str) -> int:
    ts = _timescale()
    data: dict[str, list[float]] = {}
    for jd in iter_jds():
        data[repr(jd)] = evaluate(jd, ts)
    payload = {
        "meta": {
            "git_head": _git_head(),
            "module_sha256": hashlib.sha256(MODULE_PATH.read_bytes()).hexdigest(),
            "n_epochs": len(data),
            "layout": "12 moon floats (km, J2000 ecliptic) + 3 COB floats (AU)",
        },
        "cases": data,
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh)
    print(f"dumped {len(data)} epochs -> {path}")
    print(f"  module sha256: {payload['meta']['module_sha256']}")
    return 0


def do_check(path: str, tol_km: float, tol_cob_au: float) -> int:
    with open(path, encoding="utf-8") as fh:
        payload = json.load(fh)
    baseline = payload["cases"]
    ts = _timescale()
    checked = drifted = missing = 0
    worst_moon: tuple[float, str, str] = (0.0, "", "")
    worst_cob: tuple[float, str] = (0.0, "")
    for jd in iter_jds():
        key = repr(jd)
        if key not in baseline:
            missing += 1
            continue
        checked += 1
        expected = baseline[key]
        got = evaluate(jd, ts)
        moon_delta = 0.0
        moon_which = ""
        for i in range(12):
            d = abs(got[i] - expected[i])
            if d > moon_delta:
                moon_delta = d
                moon_which = f"{MOON_NAMES[i // 3]}.{'xyz'[i % 3]}"
        cob_delta = max(abs(got[12 + i] - expected[12 + i]) for i in range(3))
        if moon_delta > worst_moon[0]:
            worst_moon = (moon_delta, key, moon_which)
        if cob_delta > worst_cob[0]:
            worst_cob = (cob_delta, key)
        if moon_delta > tol_km or cob_delta > tol_cob_au:
            drifted += 1
            print(
                f"DRIFT jd={key}: moon {moon_delta:.3e} km ({moon_which}), "
                f"cob {cob_delta:.3e} AU"
            )
    print(
        f"check: {checked} epochs, {drifted} drifted, {missing} not in baseline\n"
        f"  worst moon delta {worst_moon[0]:.3e} km"
        f" ({worst_moon[2] or 'n/a'} @ jd {worst_moon[1] or 'n/a'})\n"
        f"  worst COB delta  {worst_cob[0]:.3e} AU @ jd {worst_cob[1] or 'n/a'}"
    )
    return 1 if drifted else 0


def _find_spk(explicit: str | None) -> Path | None:
    candidates = []
    if explicit:
        candidates.append(Path(explicit))
    env_dir = os.environ.get("LIBEPHEMERIS_SPK_DIR")
    if env_dir:
        candidates.append(Path(env_dir) / "jup365.bsp")
    candidates.append(REPO_ROOT / "jup365.bsp")
    candidates.append(Path.home() / ".libephemeris" / "spk" / "jup365.bsp")
    for cand in candidates:
        if cand.is_file():
            return cand
    return None


def do_compare(spk_arg: str | None) -> int:
    """Informational physical-accuracy check vs a JPL satellite SPK."""
    spk_path = _find_spk(spk_arg)
    if spk_path is None:
        print("compare: no jup365.bsp found (use --spk PATH); skipping")
        return 0
    from jplephem.spk import SPK

    from libephemeris.moon_theories.constants import _ecliptic_to_equatorial_j2000
    from libephemeris.moon_theories.galilean import galilean_moon_positions

    kernel = SPK.open(str(spk_path))
    pairs = {(seg.center, seg.target) for seg in kernel.segments}
    moon_ids = (501, 502, 503, 504)
    needed = {(5, 599), *((5, mid) for mid in moon_ids)}
    if not needed.issubset(pairs):
        print(f"compare: {spk_path.name} lacks segments {needed - pairs}; skipping")
        return 0

    worst = {name: 0.0 for name in MOON_NAMES}
    # SPK time argument is TDB; TT is within 1.7 ms — irrelevant at the
    # ~100 km level this informational mode operates at.
    jds = [GRID_START_JD + 64.0 * k for k in range(2283)]  # 64-day sub-grid
    used = 0
    for jd in jds:
        try:
            jup = kernel[5, 599].compute(jd)
        except ValueError:
            continue  # outside SPK coverage
        used += 1
        theory = galilean_moon_positions(jd)
        for idx, mid in enumerate(moon_ids):
            spk_xyz = kernel[5, mid].compute(jd) - jup  # km, equatorial ICRF
            th_eq = _ecliptic_to_equatorial_j2000(theory[idx])
            err = math.dist(th_eq, tuple(spk_xyz))
            worst[MOON_NAMES[idx]] = max(worst[MOON_NAMES[idx]], err)
    kernel.close()
    print(f"compare vs {spk_path.name} ({used} epochs in coverage):")
    for name in MOON_NAMES:
        print(f"  {name:9s} worst |delta| = {worst[name]:9.1f} km")
    print("  (E5-from-Meeus inherent accuracy is ~100-200 km; informational)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", nargs="?", default=None, choices=["compare"])
    ap.add_argument("--dump", metavar="FILE")
    ap.add_argument("--check", metavar="FILE")
    ap.add_argument("--tol-km", type=float, default=1e-3)
    ap.add_argument("--tol-cob-au", type=float, default=1e-12)
    ap.add_argument("--spk", metavar="FILE", help="explicit jup365.bsp path")
    args = ap.parse_args()

    if args.dump:
        return do_dump(args.dump)
    if args.check:
        return do_check(args.check, args.tol_km, args.tol_cob_au)
    if args.mode == "compare":
        return do_compare(args.spk)
    ap.error("choose --dump FILE, --check FILE, or compare")
    return 2


if __name__ == "__main__":
    sys.exit(main())
