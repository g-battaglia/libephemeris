#!/usr/bin/env python3
"""
Round 7: Heliacal Events Audit — FAST PARTS
=============================================
Parts that don't require slow day-by-day scanning:
  P4: Flag constants
  P5: API shape
  P3: vis_limit_mag (instant calculation)
  P2: heliacal_pheno_ut (instant calculation)
"""

from __future__ import annotations

import os
import sys
import time
import traceback

import swisseph as swe

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
import libephemeris as ephem

EPHE_PATH = os.path.join(os.path.dirname(__file__), "..", "swisseph", "ephe")
if os.path.exists(EPHE_PATH):
    swe.set_ephe_path(EPHE_PATH)

STANDARD_ATMO = (1013.25, 15.0, 50.0, 0.25)
STANDARD_OBSERVER = (25.0, 1.0, 1, 1, 0, 0)

total = 0
passed = 0
failed = 0
failures = []


def record(name, ok, detail=""):
    global total, passed, failed, failures
    total += 1
    if ok:
        passed += 1
        print(f"  [PASS] {name}: {detail}")
    else:
        failed += 1
        failures.append((name, detail))
        print(f"  [FAIL] {name}: {detail}")


# ============================================================================
# PART 4: Flag Constants
# ============================================================================
print("=" * 70)
print("PART 4: Flag Constants Verification")
print("=" * 70)

checks = [
    ("HELIACAL_RISING", ephem.HELIACAL_RISING, swe.HELIACAL_RISING),
    ("HELIACAL_SETTING", ephem.HELIACAL_SETTING, swe.HELIACAL_SETTING),
    ("EVENING_FIRST", ephem.EVENING_FIRST, swe.EVENING_FIRST),
    ("MORNING_LAST", ephem.MORNING_LAST, swe.MORNING_LAST),
    (
        "HELFLAG_OPTICAL_PARAMS",
        ephem.HELFLAG_OPTICAL_PARAMS,
        swe.HELFLAG_OPTICAL_PARAMS,
    ),
    ("HELFLAG_NO_DETAILS", ephem.HELFLAG_NO_DETAILS, swe.HELFLAG_NO_DETAILS),
    ("HELFLAG_VISLIM_DARK", ephem.HELFLAG_VISLIM_DARK, swe.HELFLAG_VISLIM_DARK),
    (
        "HELFLAG_VISLIM_NOMOON",
        ephem.HELFLAG_VISLIM_NOMOON,
        swe.HELFLAG_VISLIM_NOMOON,
    ),
]
for name, le_val, se_val in checks:
    record(
        f"const/{name}",
        le_val == se_val,
        f"LE={le_val} SE={se_val}" + ("" if le_val == se_val else " *** MISMATCH ***"),
    )

# Missing constants
missing = [
    ("HELFLAG_HIGH_PRECISION", swe.HELFLAG_HIGH_PRECISION),
    ("HELFLAG_LONG_SEARCH", swe.HELFLAG_LONG_SEARCH),
    ("HELFLAG_SEARCH_1_PERIOD", swe.HELFLAG_SEARCH_1_PERIOD),
    ("HELFLAG_VISLIM_PHOTOPIC", swe.HELFLAG_VISLIM_PHOTOPIC),
    ("HELFLAG_AV", swe.HELFLAG_AV),
    ("HELFLAG_AVKIND", swe.HELFLAG_AVKIND),
]
for name, se_val in missing:
    has = hasattr(ephem, f"SE_{name}") or hasattr(ephem, name)
    if has:
        le_val = getattr(ephem, f"SE_{name}", getattr(ephem, name, None))
        record(f"missing/{name}", le_val == se_val, f"LE={le_val} SE={se_val}")
    else:
        record(f"missing/{name}", False, f"Missing (SE={se_val})")

# ============================================================================
# PART 5: API Shape
# ============================================================================
print("\n" + "=" * 70)
print("PART 5: API Shape / Return Type Verification")
print("=" * 70)

jd_start = swe.julday(2024, 1, 1, 0.0)
geopos = (0.0, 30.0, 0)

# --- heliacal_pheno_ut ---
print("\n  --- heliacal_pheno_ut ---")
jd_pheno = swe.julday(2024, 6, 15, 12.0)
try:
    ret_se = swe.heliacal_pheno_ut(
        jd_pheno,
        geopos,
        STANDARD_ATMO,
        STANDARD_OBSERVER,
        "Venus",
        1,
        0,
    )
    ret_le = ephem.heliacal_pheno_ut(
        jd_pheno,
        geopos,
        STANDARD_ATMO,
        STANDARD_OBSERVER,
        "Venus",
        1,
    )

    se_flat = isinstance(ret_se, tuple) and len(ret_se) == 50
    if isinstance(ret_le, tuple) and len(ret_le) == 2 and isinstance(ret_le[0], tuple):
        record(
            "pheno/shape",
            False,
            f"SE returns flat 50-tuple, LE returns (data, retflag) — API MISMATCH",
        )
        data_le = ret_le[0]
    elif isinstance(ret_le, tuple) and len(ret_le) == 50:
        record("pheno/shape", True, "Both flat 50-tuple")
        data_le = ret_le
    else:
        record("pheno/shape", False, f"LE shape unexpected: len={len(ret_le)}")
        data_le = ret_le[0] if isinstance(ret_le[0], tuple) else ret_le

    # Compare key fields
    labels = {
        0: ("AltO", 1.0),
        1: ("AppAltO", 1.0),
        2: ("GeoAltO", 1.0),
        3: ("AziO", 2.0),
        4: ("AltS", 1.0),
        5: ("AziS", 2.0),
        6: ("TAVact", 1.0),
        7: ("ARCVact", 1.0),
        8: ("DAZact", 2.0),
        9: ("ARCLact", 1.0),
        10: ("kact", 0.15),
        19: ("ParO", 0.01),
        20: ("Magn", 1.0),
        26: ("CVAact", 2.0),
    }
    for idx, (label, tol) in labels.items():
        if idx >= len(ret_se) or idx >= len(data_le):
            continue
        vs = ret_se[idx]
        vl = data_le[idx]
        if vs == 0 and vl == 0:
            continue
        if vs == 99999999.0 or vl == 99999999.0:
            if vs == 99999999.0 != (vl == 99999999.0):
                record(
                    f"pheno/Venus/{label}", False, f"SE={vs} LE={vl} sentinel mismatch"
                )
            continue
        d = abs(vs - vl)
        record(f"pheno/Venus/{label}", d < tol, f"SE={vs:.4f} LE={vl:.4f} diff={d:.4f}")

except Exception as e:
    record(f"pheno/error", False, str(e))

# --- vis_limit_mag ---
print("\n  --- vis_limit_mag ---")
jd_night = swe.julday(2024, 8, 15, 22.0)
geopos_rome = (12.5, 41.9, 0)

for obj in ["Jupiter", "Venus", "Sirius", "Mars"]:
    try:
        ret_se = swe.vis_limit_mag(
            jd_night, geopos_rome, STANDARD_ATMO, STANDARD_OBSERVER, obj, 0
        )
        ret_le = ephem.vis_limit_mag(
            jd_night, geopos_rome, STANDARD_ATMO, STANDARD_OBSERVER, obj, 0
        )

        flag_se = int(ret_se[0])
        flag_le = ret_le[0]
        dret_se = ret_se[1]
        dret_le = ret_le[1]

        record(f"vlm/{obj}/flag", flag_se == flag_le, f"SE={flag_se} LE={flag_le}")
        record(
            f"vlm/{obj}/dret_len",
            len(dret_se) == len(dret_le),
            f"SE={len(dret_se)} LE={len(dret_le)}",
        )

        if flag_se == -2:
            continue

        labels = [
            "lim_mag",
            "obj_alt",
            "obj_az",
            "sun_alt",
            "sun_az",
            "moon_alt",
            "moon_az",
            "obj_mag",
        ]
        tols = [1.5, 0.5, 1.0, 0.5, 1.0, 0.5, 1.0, 1.0]
        for i, (lbl, tol) in enumerate(zip(labels, tols)):
            if i >= min(len(dret_se), len(dret_le)):
                break
            vs = dret_se[i]
            vl = dret_le[i]
            if vs == -100.0 and vl == 0.0:
                record(f"vlm/{obj}/{lbl}", False, f"SE=-100(sentinel) LE=0")
                continue
            d = abs(vs - vl)
            record(f"vlm/{obj}/{lbl}", d < tol, f"SE={vs:.4f} LE={vl:.4f} diff={d:.4f}")

    except Exception as e:
        record(f"vlm/{obj}/error", False, str(e))

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "=" * 70)
print("SUMMARY (fast parts only)")
print("=" * 70)
print(f"Total:  {total}")
print(f"Passed: {passed}")
print(f"Failed: {failed}")
if failures:
    print(f"\n--- {len(failures)} FAILURES ---")
    for n, d in failures:
        print(f"  {n}: {d}")
print(f"\nPass rate: {passed}/{total} = {100 * passed / max(total, 1):.1f}%")
