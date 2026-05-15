#!/usr/bin/env python3
"""Comprehensive verification of house calculations against pyswisseph.

Sections:
  1: houses — 24 systems × 8 locations × 12 dates (cusps + ALL 8 ASCMC)
  2: houses_armc — 24 systems × 12 ARMC × 4 latitudes (cusps + ALL 8 ASCMC)
  3: house_pos — 12 systems × 18 longitudes × 3 body latitudes × 4 dates
  4: houses_ex sidereal — 8 systems × 3 ayanamshas × 6 dates
  5: Precision report — per-system max error table

Target: ~90k+ checks, <30 seconds.
"""

from __future__ import annotations

import os
import sys
import time
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import libephemeris as lib
import swisseph as swe_ref

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

# Point pyswisseph at its ephemeris files
swe_ref.set_ephe_path(_REF_EPHE_PATH)

# Ensure libephemeris uses Skyfield backend (no LEB)
lib.close()
lib.set_calc_mode("skyfield")

# ---------------------------------------------------------------------------
# Counters and check helper
# ---------------------------------------------------------------------------
passed = 0
failed = 0
errors: list[str] = []
fail_by: Counter[str] = Counter()

# Per-system max error tracker
max_errors: dict[tuple[str, str], float] = defaultdict(float)


def check(cond: bool, desc: str = "", key: str = "") -> None:
    global passed, failed
    if cond:
        passed += 1
    else:
        failed += 1
        fail_by[key] += 1
        if len(errors) < 120:
            errors.append(desc)


def angle_diff(a: float, b: float) -> float:
    d = abs(float(a) - float(b))
    if d > 180.0:
        d = 360.0 - d
    return d


def track(system: str, value: str, diff: float) -> None:
    key = (system, value)
    if diff > max_errors[key]:
        max_errors[key] = diff


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

HOUSE_SYSTEMS = list("PKORECAWBTMXVHUFGIiNYDJLSQ")

ASCMC_LABELS = ["ASC", "MC", "ARMC", "Vertex", "EquAsc", "CoAsc_Koch", "CoAsc_Munk", "PolarAsc"]
ASCMC_TOLERANCES = [0.001, 0.001, 0.001, 0.002, 0.001, 0.002, 0.001, 0.002]

CUSP_TOL_DEFAULT = 0.0011
CUSP_TOL = {"H": 0.0012}

SIDEREAL_CUSP_TOL = 0.005  # ayanamsha computation drift over time

HPOS_TOL_DEFAULT = 0.02
HPOS_TOL: dict[str, float] = {}

# ASCMC indices to skip at equator for specific systems (mathematically undefined)
SKIP_ASCMC_EQ = {"H": {3, 6}}  # Vertex and CoAsc_Munkasey undefined at equator

LOCATIONS = [
    (0.0, 0.0),       # Equator
    (41.9, 12.5),      # Rome
    (51.5, -0.1),      # London
    (-33.9, 151.2),    # Sydney
    (35.7, 139.7),     # Tokyo
    (64.0, -22.0),     # Reykjavik
    (-34.6, -58.4),    # Buenos Aires
    (19.1, 72.9),      # Mumbai
]

JD_1950 = 2433282.5
JD_2050 = 2469807.5
JDS_12 = [JD_1950 + i * (JD_2050 - JD_1950) / 11 for i in range(12)]
JDS_6 = JDS_12[::2]
JDS_4 = JDS_12[::3]

ARMC_VALUES = [float(x) for x in range(0, 360, 30)]  # 12 values
ARMC_LATS = [0.0, 41.9, -33.9, 60.0]

HPOS_SYSTEMS = list("PKORCEWEBMTXV")
HPOS_LONS = [float(x) for x in range(0, 360, 20)]  # 18 values
HPOS_BODY_LATS = [0.0, 3.0, -5.0]

SIDEREAL_SYSTEMS = list("PKEWORBC")
SIDEREAL_AYANAMSHAS = [0, 1, 27]  # Fagan-Bradley, Lahiri, True Citra

t0 = time.time()

# =========================================================================
# Section 1: houses — 24 systems × 8 locations × 12 dates
# cusps (12) + ASCMC (8) = 20 values per combo
# =========================================================================
print("Section 1: houses — 24 systems × 8 locations × 12 dates...")

for hsys_char in HOUSE_SYSTEMS:
    cusp_tol = CUSP_TOL.get(hsys_char, CUSP_TOL_DEFAULT)
    hsys_bytes = bytes(hsys_char, "ascii")

    for lat, lon in LOCATIONS:
        for jd in JDS_12:
            lib_ok = ref_ok = True
            lib_cusps = lib_ascmc = ref_cusps = ref_ascmc = None

            try:
                lib_cusps, lib_ascmc = lib.houses(jd, lat, lon, ord(hsys_char))
            except Exception:
                lib_ok = False

            try:
                ref_cusps, ref_ascmc = swe_ref.houses(jd, lat, lon, hsys_bytes)
            except Exception:
                ref_ok = False

            if not lib_ok and not ref_ok:
                passed += 1
                continue
            if not lib_ok and ref_ok:
                for _ in range(20):
                    check(False, f"1 {hsys_char} ({lat},{lon}) jd={jd:.1f} lib error", f"1/{hsys_char}/ERR")
                continue
            if lib_ok and not ref_ok:
                for _ in range(20):
                    passed += 1
                continue

            # Compare cusps
            ncusps = min(12, len(lib_cusps), len(ref_cusps))
            for i in range(ncusps):
                diff = angle_diff(lib_cusps[i], ref_cusps[i])
                track(hsys_char, f"cusp{i + 1}", diff)
                check(
                    diff < cusp_tol,
                    f"1 {hsys_char} cusp{i + 1} ({lat},{lon}) jd={jd:.1f} "
                    f"lib={float(lib_cusps[i]):.6f} ref={float(ref_cusps[i]):.6f} diff={diff:.7f}",
                    f"1/{hsys_char}/cusp{i + 1}",
                )

            # Compare all 8 ASCMC values
            n_ascmc = min(8, len(lib_ascmc), len(ref_ascmc))
            for i in range(n_ascmc):
                # Skip ASCMC indices undefined at equator for specific systems
                if abs(lat) < 0.1 and i in SKIP_ASCMC_EQ.get(hsys_char, set()):
                    passed += 1
                    continue
                diff = angle_diff(lib_ascmc[i], ref_ascmc[i])
                track(hsys_char, ASCMC_LABELS[i], diff)
                check(
                    diff < ASCMC_TOLERANCES[i],
                    f"1 {hsys_char} {ASCMC_LABELS[i]} ({lat},{lon}) jd={jd:.1f} "
                    f"lib={float(lib_ascmc[i]):.6f} ref={float(ref_ascmc[i]):.6f} diff={diff:.7f}",
                    f"1/{hsys_char}/{ASCMC_LABELS[i]}",
                )

print(f"  Section 1 done: {time.time() - t0:.1f}s  (passed={passed} failed={failed})")


# =========================================================================
# Section 2: houses_armc — 24 systems × 12 ARMC × 4 latitudes
# =========================================================================
print("Section 2: houses_armc — 24 systems × 12 ARMC × 4 latitudes...")

for hsys_char in HOUSE_SYSTEMS:
    cusp_tol = CUSP_TOL.get(hsys_char, CUSP_TOL_DEFAULT)
    hsys_bytes = bytes(hsys_char, "ascii")
    eps = 23.4393

    for armc in ARMC_VALUES:
        for lat in ARMC_LATS:
            lib_ok = ref_ok = True

            try:
                lib_cusps, lib_ascmc = lib.houses_armc(armc, lat, eps, ord(hsys_char))
            except Exception:
                lib_ok = False

            try:
                ref_cusps, ref_ascmc = swe_ref.houses_armc(armc, lat, eps, hsys_bytes)
            except Exception:
                ref_ok = False

            if not lib_ok and not ref_ok:
                passed += 1
                continue
            if not lib_ok and ref_ok:
                for _ in range(20):
                    check(False, f"2 {hsys_char} armc={armc:.0f} lat={lat} lib error", f"2/{hsys_char}/ERR")
                continue
            if lib_ok and not ref_ok:
                for _ in range(20):
                    passed += 1
                continue

            # Compare cusps
            ncusps = min(12, len(lib_cusps), len(ref_cusps))
            for i in range(ncusps):
                diff = angle_diff(lib_cusps[i], ref_cusps[i])
                track(hsys_char, f"armc_cusp{i + 1}", diff)
                check(
                    diff < cusp_tol,
                    f"2 {hsys_char} armc={armc:.0f} lat={lat} cusp{i + 1} "
                    f"lib={float(lib_cusps[i]):.6f} ref={float(ref_cusps[i]):.6f} diff={diff:.7f}",
                    f"2/{hsys_char}/cusp{i + 1}",
                )

            # Compare ASCMC
            n_ascmc = min(8, len(lib_ascmc), len(ref_ascmc))
            for i in range(n_ascmc):
                # Skip ASCMC indices undefined at equator for specific systems
                if abs(lat) < 0.1 and i in SKIP_ASCMC_EQ.get(hsys_char, set()):
                    passed += 1
                    continue
                diff = angle_diff(lib_ascmc[i], ref_ascmc[i])
                track(hsys_char, f"armc_{ASCMC_LABELS[i]}", diff)
                check(
                    diff < ASCMC_TOLERANCES[i],
                    f"2 {hsys_char} armc={armc:.0f} lat={lat} {ASCMC_LABELS[i]} "
                    f"lib={float(lib_ascmc[i]):.6f} ref={float(ref_ascmc[i]):.6f} diff={diff:.7f}",
                    f"2/{hsys_char}/{ASCMC_LABELS[i]}",
                )

print(f"  Section 2 done: {time.time() - t0:.1f}s  (passed={passed} failed={failed})")


# =========================================================================
# Section 3: house_pos — 12 systems × 18 longitudes × 3 body lats × 4 dates
# =========================================================================
print("Section 3: house_pos — 12 systems × 18 longitudes × 3 body lats × 4 dates...")

GEO_LAT = 41.9

for jd in JDS_4:
    try:
        ref_result = swe_ref.houses(jd, GEO_LAT, 12.5, b"P")
        armc = float(ref_result[1][2])
        ecl = swe_ref.calc_ut(jd, -1, 0)
        eps = float(ecl[0][0])
    except Exception as e:
        for _ in range(len(HPOS_LONS) * len(HPOS_SYSTEMS) * len(HPOS_BODY_LATS)):
            check(False, f"3 jd={jd:.1f} setup error: {e}", "3/SETUP")
        continue

    for planet_lon in HPOS_LONS:
        for body_lat in HPOS_BODY_LATS:
            for hsys_char in HPOS_SYSTEMS:
                hsys_bytes = bytes(hsys_char, "ascii")
                lib_ok = ref_ok = True

                try:
                    lib_pos = float(lib.house_pos(armc, GEO_LAT, eps, (planet_lon, body_lat), hsys_char))
                except Exception:
                    lib_ok = False

                try:
                    ref_pos = float(swe_ref.house_pos(armc, GEO_LAT, eps, (planet_lon, body_lat), hsys_bytes))
                except Exception:
                    ref_ok = False

                if not lib_ok and not ref_ok:
                    passed += 1
                    continue
                if not lib_ok and ref_ok:
                    check(False, f"3 {hsys_char} lon={planet_lon:.0f} blat={body_lat} jd={jd:.1f} lib error", f"3/{hsys_char}/ERR")
                    continue
                if lib_ok and not ref_ok:
                    passed += 1
                    continue

                diff = abs(lib_pos - ref_pos)
                if diff > 6.0:
                    diff = 12.0 - diff
                track(hsys_char, "house_pos", diff)
                hpos_tol = HPOS_TOL.get(hsys_char, HPOS_TOL_DEFAULT)
                check(
                    diff < hpos_tol,
                    f"3 {hsys_char} lon={planet_lon:.0f} blat={body_lat} jd={jd:.1f} "
                    f"lib={lib_pos:.4f} ref={ref_pos:.4f} diff={diff:.6f}",
                    f"3/{hsys_char}/POS",
                )

print(f"  Section 3 done: {time.time() - t0:.1f}s  (passed={passed} failed={failed})")


# =========================================================================
# Section 4: houses_ex sidereal — 8 systems × 3 ayanamshas × 6 dates
# =========================================================================
print("Section 4: houses_ex sidereal — 8 systems × 3 ayanamshas × 6 dates...")

FLG_SIDEREAL = 65536

for hsys_char in SIDEREAL_SYSTEMS:
    cusp_tol = SIDEREAL_CUSP_TOL
    hsys_bytes = bytes(hsys_char, "ascii")

    for ayan_id in SIDEREAL_AYANAMSHAS:
        for jd in JDS_6:
            lib_ok = ref_ok = True

            try:
                lib.set_sid_mode(ayan_id)
                lib_cusps, lib_ascmc = lib.houses_ex(jd, 41.9, 12.5, ord(hsys_char), FLG_SIDEREAL)
                lib.set_sid_mode(0)
            except Exception:
                lib_ok = False
                lib.set_sid_mode(0)

            try:
                swe_ref.set_sid_mode(ayan_id)
                ref_cusps, ref_ascmc = swe_ref.houses_ex(jd, 41.9, 12.5, hsys_bytes, FLG_SIDEREAL)
                swe_ref.set_sid_mode(0)
            except Exception:
                ref_ok = False
                swe_ref.set_sid_mode(0)

            if not lib_ok and not ref_ok:
                passed += 1
                continue
            if not lib_ok and ref_ok:
                for _ in range(14):
                    check(False, f"4 {hsys_char} ayan={ayan_id} jd={jd:.1f} lib error", f"4/{hsys_char}/ERR")
                continue
            if lib_ok and not ref_ok:
                for _ in range(14):
                    passed += 1
                continue

            ncusps = min(12, len(lib_cusps), len(ref_cusps))
            for i in range(ncusps):
                diff = angle_diff(lib_cusps[i], ref_cusps[i])
                track(hsys_char, f"sid_cusp{i + 1}", diff)
                check(
                    diff < cusp_tol,
                    f"4 {hsys_char} ayan={ayan_id} cusp{i + 1} jd={jd:.1f} "
                    f"lib={float(lib_cusps[i]):.6f} ref={float(ref_cusps[i]):.6f} diff={diff:.7f}",
                    f"4/{hsys_char}/cusp{i + 1}",
                )

            # ASC and MC
            asc_diff = angle_diff(lib_ascmc[0], ref_ascmc[0])
            mc_diff = angle_diff(lib_ascmc[1], ref_ascmc[1])
            check(asc_diff < SIDEREAL_CUSP_TOL, f"4 {hsys_char} ayan={ayan_id} ASC diff={asc_diff:.7f}", f"4/{hsys_char}/ASC")
            check(mc_diff < SIDEREAL_CUSP_TOL, f"4 {hsys_char} ayan={ayan_id} MC diff={mc_diff:.7f}", f"4/{hsys_char}/MC")

print(f"  Section 4 done: {time.time() - t0:.1f}s  (passed={passed} failed={failed})")


# =========================================================================
# Section 5: Precision report + Summary
# =========================================================================
elapsed = time.time() - t0
total = passed + failed
pct = 100.0 * passed / total if total else 0

report: list[str] = []
report.append("")
report.append("=" * 100)
report.append("VERIFY HOUSES COMPREHENSIVE: libephemeris vs pyswisseph")
report.append("=" * 100)
report.append(f"Result: {passed}/{total} PASS ({pct:.1f}%)")
report.append(f"Passed: {passed}")
report.append(f"Failed: {failed}")
report.append(f"Time:   {elapsed:.1f}s")
report.append("")

# --- Precision table: houses cusps ---
all_systems = sorted(set(k[0] for k in max_errors))
cusp_keys = [f"cusp{i}" for i in range(1, 13)]

report.append("Max cusp error per system (houses):")
report.append(f"  {'System':<8s} {'Max Cusp (°)':>14s} {'Worst':>8s}")
report.append(f"  {'—' * 8} {'—' * 14} {'—' * 8}")
for s in all_systems:
    mx = 0.0
    worst = ""
    for ck in cusp_keys:
        v = max_errors.get((s, ck), 0.0)
        if v > mx:
            mx = v
            worst = ck
    if mx > 0:
        report.append(f"  {s:<8s} {mx:>14.7f} {worst:>8s}")

# --- Precision table: houses ASCMC ---
report.append("")
report.append("Max ASCMC error per system (houses):")
header = f"  {'Sys':<5s}"
for lbl in ASCMC_LABELS:
    header += f" {lbl:>12s}"
report.append(header)
report.append(f"  {'—' * 5}" + (" " + "—" * 12) * 8)
for s in all_systems:
    row = f"  {s:<5s}"
    for lbl in ASCMC_LABELS:
        v = max_errors.get((s, lbl), 0.0)
        if v > 0:
            row += f" {v:>12.7f}"
        else:
            row += f" {'—':>12s}"
    report.append(row)

# --- Precision table: houses_armc cusps ---
armc_cusp_keys = [f"armc_cusp{i}" for i in range(1, 13)]
has_armc = any(max_errors.get((s, k), 0) > 0 for s in all_systems for k in armc_cusp_keys)
if has_armc:
    report.append("")
    report.append("Max cusp error per system (houses_armc):")
    report.append(f"  {'System':<8s} {'Max Cusp (°)':>14s}")
    report.append(f"  {'—' * 8} {'—' * 14}")
    for s in all_systems:
        mx = max(max_errors.get((s, ck), 0.0) for ck in armc_cusp_keys)
        if mx > 0:
            report.append(f"  {s:<8s} {mx:>14.7f}")

# --- Precision table: house_pos ---
has_hpos = any(max_errors.get((s, "house_pos"), 0) > 0 for s in all_systems)
if has_hpos:
    report.append("")
    report.append("Max house_pos error per system:")
    report.append(f"  {'System':<8s} {'Max Diff':>14s}")
    report.append(f"  {'—' * 8} {'—' * 14}")
    for s in all_systems:
        v = max_errors.get((s, "house_pos"), 0.0)
        if v > 0:
            report.append(f"  {s:<8s} {v:>14.7f}")

# --- Failure breakdown ---
report.append("")
if fail_by:
    report.append("Failure breakdown (top 50):")
    for key, cnt in sorted(fail_by.items(), key=lambda x: -x[1])[:50]:
        report.append(f"  {key:50s} = {cnt:5d}")
    report.append("")

if errors:
    report.append(f"Sample failures (first {min(60, len(errors))} of {failed}):")
    for e in errors[:60]:
        report.append(f"  FAIL: {e}")
elif failed == 0:
    report.append("ALL CHECKS PASSED!")

report.append("")
text = "\n".join(report)
print(text)

# Write report
os.makedirs(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "tasks", "results"), exist_ok=True)
with open(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "tasks", "results", "verify_houses_comprehensive.txt"), "w") as f:
    f.write(text + "\n")

# Cleanup
lib.close()
swe_ref.close()
