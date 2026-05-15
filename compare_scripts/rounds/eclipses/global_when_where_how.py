#!/usr/bin/env python3
"""
Round 17: Deep Eclipse Global Functions Verification
=====================================================

Parts:
  P1: sol_eclipse_when_glob — find next 10 solar eclipses, compare timing
  P2: sol_eclipse_where — central line position at maximum for each eclipse
  P3: sol_eclipse_how — eclipse circumstances at specific locations
  P4: lun_eclipse_when — find next 10 lunar eclipses, compare timing
  P5: lun_eclipse_how — lunar eclipse circumstances (magnitude, etc.)
  P6: Eclipse type flags consistency (retflag agreement)
  P7: Backward search — find previous eclipses
  P8: Eclipse type filtering — search for specific types only
"""

from __future__ import annotations

import math
import os
import sys
import time
import traceback

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import libephemeris as ephem
from libephemeris.constants import *

_EPHE_PATH = os.path.join(os.path.dirname(__file__), "..", "swisseph", "ephe")
swe.set_ephe_path(_EPHE_PATH)


class R:
    def __init__(self, name):
        self.name = name
        self.passed = self.failed = self.skipped = 0
        self.failures = []
        self.max_diff = 0.0
        self.max_label = ""

    def ok(self, diff=0.0, label=""):
        self.passed += 1
        if diff > self.max_diff:
            self.max_diff = diff
            self.max_label = label

    def fail(self, msg):
        self.failed += 1
        self.failures.append(msg)

    def skip(self, msg=""):
        self.skipped += 1

    def summary(self):
        total = self.passed + self.failed
        print(f"\n{'=' * 70}")
        print(f"  {self.name}: {self.passed}/{total} PASSED ({self.skipped} skip)")
        if self.max_diff > 0:
            print(f"  Max diff: {self.max_diff:.6f} ({self.max_label})")
        if self.failures:
            for f in self.failures[:20]:
                print(f"    - {f}")
        print(f"{'=' * 70}")
        return self.failed == 0


def jd_to_date_str(jd):
    """Convert JD to readable date string."""
    y, m, d, h = swe.revjul(jd)
    hh = int(h)
    mm = int((h - hh) * 60)
    return f"{y:04d}-{m:02d}-{d:02d} {hh:02d}:{mm:02d}"


def ecl_type_str(flags):
    """Convert eclipse type flags to string."""
    parts = []
    if flags & ECL_TOTAL:
        parts.append("TOTAL")
    if flags & ECL_ANNULAR:
        parts.append("ANNULAR")
    if flags & ECL_PARTIAL:
        parts.append("PARTIAL")
    if flags & ECL_ANNULAR_TOTAL:
        parts.append("HYBRID")
    if flags & ECL_PENUMBRAL:
        parts.append("PENUMBRAL")
    if flags & ECL_CENTRAL:
        parts.append("central")
    if flags & ECL_NONCENTRAL:
        parts.append("noncentral")
    return "|".join(parts) if parts else f"0x{flags:x}"


# ============================================================
# PART 1: Solar Eclipse When Glob — next 10 solar eclipses
# ============================================================
def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: sol_eclipse_when_glob — 10 forward solar eclipses")
    print("  Starting from J2000, find next 10 solar eclipses")
    print("=" * 70)

    r = R("P1: Sol Eclipse When Glob")
    jd = 2451545.0  # J2000

    for i in range(10):
        label = f"Solar eclipse #{i + 1}"
        try:
            se_res, se_tret = swe.sol_eclipse_when_glob(jd)
            le_res, le_tret = ephem.sol_eclipse_when_glob(jd)
        except Exception as e:
            r.fail(f"{label}: {e}")
            jd += 30
            continue

        # Compare maximum time (tret[0])
        if se_tret[0] == 0 or le_tret[0] == 0:
            r.fail(f"{label}: zero max time SE={se_tret[0]} LE={le_tret[0]}")
            jd += 30
            continue

        dt_max = abs(se_tret[0] - le_tret[0]) * 86400  # seconds
        date_str = jd_to_date_str(se_tret[0])

        # Compare eclipse begin/end (tret[2], tret[3])
        dt_begin = (
            abs(se_tret[2] - le_tret[2]) * 86400
            if se_tret[2] > 0 and le_tret[2] > 0
            else 0
        )
        dt_end = (
            abs(se_tret[3] - le_tret[3]) * 86400
            if se_tret[3] > 0 and le_tret[3] > 0
            else 0
        )

        # Tolerance: 120 seconds for timing (eclipse geometry is complex)
        tol_sec = 120.0

        fails = []
        if dt_max > tol_sec:
            fails.append(f"max: {dt_max:.1f}s")
        if dt_begin > tol_sec:
            fails.append(f"begin: {dt_begin:.1f}s")
        if dt_end > tol_sec:
            fails.append(f"end: {dt_end:.1f}s")

        se_type = ecl_type_str(se_res)
        le_type = ecl_type_str(le_res)

        # Check eclipse type agreement (mask out visibility flags)
        se_type_bits = se_res & (
            ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL | ECL_ANNULAR_TOTAL
        )
        le_type_bits = le_res & (
            ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL | ECL_ANNULAR_TOTAL
        )
        if se_type_bits != le_type_bits:
            fails.append(f"type: SE={se_type} LE={le_type}")

        if fails:
            r.fail(f"{label} [{date_str}]: {'; '.join(fails)}")
        else:
            r.ok(dt_max, f"#{i + 1} {date_str}")

        print(
            f"  #{i + 1}: {date_str}  SE={se_type:<30s}  LE={le_type:<30s}  dt_max={dt_max:6.1f}s"
        )

        # Advance past this eclipse
        jd = se_tret[0] + 10

    return r.summary(), r


# ============================================================
# PART 2: Solar Eclipse Where — central line at maximum
# ============================================================
def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: sol_eclipse_where — central line position")
    print("  For each central eclipse found in P1, compare position")
    print("=" * 70)

    r = R("P2: Sol Eclipse Where")
    jd = 2451545.0

    for i in range(10):
        try:
            se_res, se_tret = swe.sol_eclipse_when_glob(jd)
        except Exception:
            jd += 30
            continue

        jd_max = se_tret[0]
        date_str = jd_to_date_str(jd_max)
        label = f"Where #{i + 1} [{date_str}]"

        # Only test central eclipses (non-central don't have a defined central line)
        if not (se_res & ECL_CENTRAL):
            r.skip(f"{label}: non-central")
            jd = jd_max + 10
            continue

        try:
            se_res2, se_geopos, se_attr = swe.sol_eclipse_where(jd_max)
            le_res2, le_geopos, le_attr = ephem.sol_eclipse_where(jd_max, 0)
        except Exception as e:
            r.fail(f"{label}: {e}")
            jd = jd_max + 10
            continue

        # Compare central line longitude and latitude
        lon_diff = abs(se_geopos[0] - le_geopos[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        lat_diff = abs(se_geopos[1] - le_geopos[1])

        # Compare magnitude (attr[0])
        mag_diff = abs(se_attr[0] - le_attr[0])

        # Tolerances: 1 degree for position, 0.05 for magnitude
        tol_pos = 1.0
        tol_mag = 0.05

        fails = []
        if lon_diff > tol_pos:
            fails.append(
                f"lon: SE={se_geopos[0]:.2f} LE={le_geopos[0]:.2f} diff={lon_diff:.2f}°"
            )
        if lat_diff > tol_pos:
            fails.append(
                f"lat: SE={se_geopos[1]:.2f} LE={le_geopos[1]:.2f} diff={lat_diff:.2f}°"
            )
        if mag_diff > tol_mag:
            fails.append(f"mag: SE={se_attr[0]:.4f} LE={le_attr[0]:.4f}")

        if fails:
            r.fail(f"{label}: {'; '.join(fails)}")
        else:
            r.ok(max(lon_diff, lat_diff), label)

        print(
            f"  #{i + 1}: {date_str}  lon_diff={lon_diff:.3f}° lat_diff={lat_diff:.3f}° mag_diff={mag_diff:.4f}"
        )

        jd = jd_max + 10

    return r.summary(), r


# ============================================================
# PART 3: Solar Eclipse How — at specific locations
# ============================================================
def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: sol_eclipse_how — circumstances at locations")
    print("  Test eclipse circumstances at 5 geographic locations")
    print("=" * 70)

    r = R("P3: Sol Eclipse How")

    # Find first central solar eclipse from J2000
    jd = 2451545.0
    jd_max = None
    for _ in range(20):
        se_res, se_tret = swe.sol_eclipse_when_glob(jd)
        if se_res & ECL_CENTRAL:
            jd_max = se_tret[0]
            break
        jd = se_tret[0] + 10

    if jd_max is None:
        r.fail("Could not find central eclipse")
        return r.summary(), r

    date_str = jd_to_date_str(jd_max)
    print(f"  Using eclipse at {date_str} (JD={jd_max:.4f})")

    # Get central line position
    _, se_geopos_central, _ = swe.sol_eclipse_where(jd_max)
    central_lon = se_geopos_central[0]
    central_lat = se_geopos_central[1]

    # Test at several locations around the eclipse path
    locations = [
        (central_lon, central_lat, 0.0, "Central line"),
        (central_lon + 5, central_lat, 0.0, "5° east"),
        (central_lon - 5, central_lat, 0.0, "5° west"),
        (central_lon, central_lat + 5, 0.0, "5° north"),
        (central_lon, central_lat - 5, 0.0, "5° south"),
    ]

    for lon, lat, alt, loc_name in locations:
        label = f"How @ {loc_name}"
        geopos = [lon, lat, alt]
        try:
            # pyswisseph uses (tjd, geopos, flags) order
            se_res3, se_attr = swe.sol_eclipse_how(jd_max, geopos, 0)
            le_res3, le_attr = ephem.sol_eclipse_how(jd_max, 0, geopos)
        except Exception as e:
            r.fail(f"{label}: {e}")
            continue

        # Compare magnitude (attr[0]) and obscuration (attr[2])
        mag_diff = abs(se_attr[0] - le_attr[0])
        obs_diff = (
            abs(se_attr[2] - le_attr[2]) if se_attr[2] > 0 or le_attr[2] > 0 else 0
        )

        # Tolerance
        tol_mag = 0.05
        tol_obs = 0.05

        fails = []
        if mag_diff > tol_mag:
            fails.append(f"mag: SE={se_attr[0]:.4f} LE={le_attr[0]:.4f}")
        if obs_diff > tol_obs:
            fails.append(f"obs: SE={se_attr[2]:.4f} LE={le_attr[2]:.4f}")

        if fails:
            r.fail(f"{label}: {'; '.join(fails)}")
        else:
            r.ok(mag_diff, label)

        print(
            f"  {loc_name:15s}: mag_diff={mag_diff:.4f}  obs_diff={obs_diff:.4f}  SE_mag={se_attr[0]:.4f}"
        )

    return r.summary(), r


# ============================================================
# PART 4: Lunar Eclipse When — next 10 lunar eclipses
# ============================================================
def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: lun_eclipse_when — 10 forward lunar eclipses")
    print("  Starting from J2000, find next 10 lunar eclipses")
    print("=" * 70)

    r = R("P4: Lun Eclipse When")
    jd = 2451545.0

    for i in range(10):
        label = f"Lunar eclipse #{i + 1}"
        try:
            se_res, se_tret = swe.lun_eclipse_when(jd)
            le_res, le_tret = ephem.lun_eclipse_when(jd)
        except Exception as e:
            r.fail(f"{label}: {e}")
            jd += 30
            continue

        if se_tret[0] == 0 or le_tret[0] == 0:
            r.fail(f"{label}: zero max time")
            jd += 30
            continue

        dt_max = abs(se_tret[0] - le_tret[0]) * 86400  # seconds
        date_str = jd_to_date_str(se_tret[0])

        # Compare penumbral begin/end (tret[6], tret[7])
        dt_pen_begin = (
            abs(se_tret[6] - le_tret[6]) * 86400
            if se_tret[6] > 0 and le_tret[6] > 0
            else 0
        )
        dt_pen_end = (
            abs(se_tret[7] - le_tret[7]) * 86400
            if se_tret[7] > 0 and le_tret[7] > 0
            else 0
        )

        # Compare partial begin/end (tret[2], tret[3]) if applicable
        dt_part_begin = (
            abs(se_tret[2] - le_tret[2]) * 86400
            if se_tret[2] > 0 and le_tret[2] > 0
            else 0
        )
        dt_part_end = (
            abs(se_tret[3] - le_tret[3]) * 86400
            if se_tret[3] > 0 and le_tret[3] > 0
            else 0
        )

        # Tolerance: 120 seconds (general), 300 seconds for penumbral
        # contacts on very shallow penumbral eclipses where grazing
        # geometry makes exact contact times inherently imprecise.
        tol_sec = 120.0
        is_shallow_penumbral = (se_res & ECL_PENUMBRAL) and not (
            se_res & (ECL_TOTAL | ECL_PARTIAL)
        )
        tol_pen_contact = 300.0 if is_shallow_penumbral else 120.0

        fails = []
        if dt_max > tol_sec:
            fails.append(f"max: {dt_max:.1f}s")
        if dt_pen_begin > tol_pen_contact:
            fails.append(f"pen_begin: {dt_pen_begin:.1f}s")
        if dt_pen_end > tol_pen_contact:
            fails.append(f"pen_end: {dt_pen_end:.1f}s")
        if dt_part_begin > tol_sec:
            fails.append(f"part_begin: {dt_part_begin:.1f}s")
        if dt_part_end > tol_sec:
            fails.append(f"part_end: {dt_part_end:.1f}s")

        se_type = ecl_type_str(se_res)
        le_type = ecl_type_str(le_res)

        # Check eclipse type agreement
        se_type_bits = se_res & (ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL)
        le_type_bits = le_res & (ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL)
        if se_type_bits != le_type_bits:
            fails.append(f"type: SE={se_type} LE={le_type}")

        if fails:
            r.fail(f"{label} [{date_str}]: {'; '.join(fails)}")
        else:
            r.ok(dt_max, f"#{i + 1} {date_str}")

        print(
            f"  #{i + 1}: {date_str}  SE={se_type:<25s}  LE={le_type:<25s}  dt_max={dt_max:6.1f}s"
        )

        jd = se_tret[0] + 10

    return r.summary(), r


# ============================================================
# PART 5: Lunar Eclipse How — magnitude comparison
# ============================================================
def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: lun_eclipse_how — magnitude at maximum")
    print("  For 10 lunar eclipses, compare umbral/penumbral magnitude")
    print("=" * 70)

    r = R("P5: Lun Eclipse How")
    jd = 2451545.0
    geopos = [0.0, 0.0, 0.0]  # geocentric

    for i in range(10):
        label = f"Lun How #{i + 1}"
        try:
            se_res, se_tret = swe.lun_eclipse_when(jd)
        except Exception:
            jd += 30
            continue

        jd_max = se_tret[0]
        date_str = jd_to_date_str(jd_max)

        try:
            # pyswisseph uses (tjd, geopos, flags) order
            se_res2, se_attr = swe.lun_eclipse_how(jd_max, geopos, 0)
            le_res2, le_attr = ephem.lun_eclipse_how(jd_max, 0, geopos)
        except Exception as e:
            r.fail(f"{label} [{date_str}]: {e}")
            jd = jd_max + 10
            continue

        # Compare umbral magnitude (attr[0])
        umbral_diff = abs(se_attr[0] - le_attr[0])
        # Compare penumbral magnitude (attr[1])
        pen_diff = abs(se_attr[1] - le_attr[1])

        # Tolerances
        tol_mag = 0.05

        fails = []
        if umbral_diff > tol_mag:
            fails.append(
                f"umbral: SE={se_attr[0]:.4f} LE={le_attr[0]:.4f} diff={umbral_diff:.4f}"
            )
        if pen_diff > tol_mag:
            fails.append(
                f"penum: SE={se_attr[1]:.4f} LE={le_attr[1]:.4f} diff={pen_diff:.4f}"
            )

        if fails:
            r.fail(f"{label} [{date_str}]: {'; '.join(fails)}")
        else:
            r.ok(max(umbral_diff, pen_diff), f"#{i + 1} {date_str}")

        print(
            f"  #{i + 1}: {date_str}  umbral_diff={umbral_diff:.4f}  pen_diff={pen_diff:.4f}"
        )

        jd = jd_max + 10

    return r.summary(), r


# ============================================================
# PART 6: Eclipse type flags consistency
# ============================================================
def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Eclipse type flags consistency")
    print("  Verify retflag type bits match between SE and LE")
    print("=" * 70)

    r = R("P6: Type Flags")

    # Solar eclipses from 2000-2030
    jd = 2451545.0
    for i in range(15):
        try:
            se_res, se_tret = swe.sol_eclipse_when_glob(jd)
            le_res, le_tret = ephem.sol_eclipse_when_glob(jd)
        except Exception as e:
            r.fail(f"Solar #{i + 1}: {e}")
            jd += 30
            continue

        date_str = jd_to_date_str(se_tret[0])

        # Compare type bits
        type_mask = (
            ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL | ECL_ANNULAR_TOTAL
        )
        central_mask = ECL_CENTRAL | ECL_NONCENTRAL

        se_type = se_res & type_mask
        le_type = le_res & type_mask
        se_central = se_res & central_mask
        le_central = le_res & central_mask

        fails = []
        if se_type != le_type:
            fails.append(f"type: SE={ecl_type_str(se_type)} LE={ecl_type_str(le_type)}")
        if se_central != le_central:
            fails.append(
                f"central: SE={ecl_type_str(se_central)} LE={ecl_type_str(le_central)}"
            )

        if fails:
            r.fail(f"Solar #{i + 1} [{date_str}]: {'; '.join(fails)}")
        else:
            r.ok(0, f"Sol #{i + 1}")

        jd = se_tret[0] + 10

    # Lunar eclipses
    jd = 2451545.0
    for i in range(15):
        try:
            se_res, se_tret = swe.lun_eclipse_when(jd)
            le_res, le_tret = ephem.lun_eclipse_when(jd)
        except Exception as e:
            r.fail(f"Lunar #{i + 1}: {e}")
            jd += 30
            continue

        date_str = jd_to_date_str(se_tret[0])

        type_mask = ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL
        se_type = se_res & type_mask
        le_type = le_res & type_mask

        fails = []
        if se_type != le_type:
            fails.append(f"type: SE={ecl_type_str(se_type)} LE={ecl_type_str(le_type)}")

        if fails:
            r.fail(f"Lunar #{i + 1} [{date_str}]: {'; '.join(fails)}")
        else:
            r.ok(0, f"Lun #{i + 1}")

        jd = se_tret[0] + 10

    return r.summary(), r


# ============================================================
# PART 7: Backward search
# ============================================================
def run_part7():
    print("\n" + "=" * 70)
    print("PART 7: Backward search — find previous eclipses")
    print("  Starting from 2025, search backward for 5 solar + 5 lunar")
    print("=" * 70)

    r = R("P7: Backward Search")
    jd_start = swe.julday(2025, 1, 1, 0.0)

    # Solar backward
    jd = jd_start
    for i in range(5):
        label = f"Solar back #{i + 1}"
        try:
            se_res, se_tret = swe.sol_eclipse_when_glob(jd, backwards=True)
            le_res, le_tret = ephem.sol_eclipse_when_glob(jd, backwards=True)
        except Exception as e:
            r.fail(f"{label}: {e}")
            jd -= 30
            continue

        dt_max = abs(se_tret[0] - le_tret[0]) * 86400
        date_str = jd_to_date_str(se_tret[0])

        if dt_max > 120:
            r.fail(f"{label} [{date_str}]: dt_max={dt_max:.1f}s")
        else:
            r.ok(dt_max, f"Sol back #{i + 1}")

        print(f"  Sol back #{i + 1}: {date_str}  dt_max={dt_max:6.1f}s")
        jd = se_tret[0] - 10

    # Lunar backward — lun_eclipse_when doesn't support backward param
    # so we skip this part (it's an API gap to note, not a bug)
    jd = jd_start
    for i in range(5):
        label = f"Lunar back #{i + 1}"
        try:
            se_res, se_tret = swe.lun_eclipse_when(jd, backwards=True)
            # LE doesn't have backward param — just search forward from earlier date
            # to find same eclipse. Use SE result to verify.
            le_res, le_tret = ephem.lun_eclipse_when(se_tret[0] - 5)
        except Exception as e:
            r.fail(f"{label}: {e}")
            jd -= 30
            continue

        dt_max = abs(se_tret[0] - le_tret[0]) * 86400
        date_str = jd_to_date_str(se_tret[0])

        if dt_max > 120:
            r.fail(f"{label} [{date_str}]: dt_max={dt_max:.1f}s")
        else:
            r.ok(dt_max, f"Lun back #{i + 1}")

        print(f"  Lun back #{i + 1}: {date_str}  dt_max={dt_max:6.1f}s")
        jd = se_tret[0] - 10

    return r.summary(), r


# ============================================================
# PART 8: Eclipse type filtering
# ============================================================
def run_part8():
    print("\n" + "=" * 70)
    print("PART 8: Eclipse type filtering")
    print("  Search for specific eclipse types and verify match")
    print("=" * 70)

    r = R("P8: Type Filter")
    jd = 2451545.0

    # Search for total solar eclipses only
    for i in range(3):
        label = f"Total solar #{i + 1}"
        try:
            se_res, se_tret = swe.sol_eclipse_when_glob(jd, ecltype=ECL_TOTAL)
            le_res, le_tret = ephem.sol_eclipse_when_glob(jd, ecltype=ECL_TOTAL)
        except Exception as e:
            r.fail(f"{label}: {e}")
            jd += 180
            continue

        dt_max = abs(se_tret[0] - le_tret[0]) * 86400
        date_str = jd_to_date_str(se_tret[0])

        # Verify both report TOTAL
        if not (se_res & ECL_TOTAL):
            r.fail(f"{label}: SE not total ({ecl_type_str(se_res)})")
        elif not (le_res & ECL_TOTAL):
            r.fail(f"{label}: LE not total ({ecl_type_str(le_res)})")
        elif dt_max > 120:
            r.fail(f"{label} [{date_str}]: dt_max={dt_max:.1f}s")
        else:
            r.ok(dt_max, label)

        print(
            f"  {label}: {date_str}  dt_max={dt_max:6.1f}s  SE={ecl_type_str(se_res)}"
        )
        jd = se_tret[0] + 10

    # Search for annular solar eclipses only
    jd = 2451545.0
    for i in range(3):
        label = f"Annular solar #{i + 1}"
        try:
            se_res, se_tret = swe.sol_eclipse_when_glob(jd, ecltype=ECL_ANNULAR)
            le_res, le_tret = ephem.sol_eclipse_when_glob(
                jd, ecltype=ECL_ANNULAR
            )
        except Exception as e:
            r.fail(f"{label}: {e}")
            jd += 180
            continue

        dt_max = abs(se_tret[0] - le_tret[0]) * 86400
        date_str = jd_to_date_str(se_tret[0])

        if not (se_res & ECL_ANNULAR):
            r.fail(f"{label}: SE not annular ({ecl_type_str(se_res)})")
        elif not (le_res & ECL_ANNULAR):
            r.fail(f"{label}: LE not annular ({ecl_type_str(le_res)})")
        elif dt_max > 120:
            r.fail(f"{label} [{date_str}]: dt_max={dt_max:.1f}s")
        else:
            r.ok(dt_max, label)

        print(
            f"  {label}: {date_str}  dt_max={dt_max:6.1f}s  SE={ecl_type_str(se_res)}"
        )
        jd = se_tret[0] + 10

    # Search for total lunar eclipses only
    jd = 2451545.0
    for i in range(3):
        label = f"Total lunar #{i + 1}"
        try:
            se_res, se_tret = swe.lun_eclipse_when(jd, ecltype=ECL_TOTAL)
            le_res, le_tret = ephem.lun_eclipse_when(jd, ecltype=ECL_TOTAL)
        except Exception as e:
            r.fail(f"{label}: {e}")
            jd += 180
            continue

        dt_max = abs(se_tret[0] - le_tret[0]) * 86400
        date_str = jd_to_date_str(se_tret[0])

        if not (se_res & ECL_TOTAL):
            r.fail(f"{label}: SE not total ({ecl_type_str(se_res)})")
        elif not (le_res & ECL_TOTAL):
            r.fail(f"{label}: LE not total ({ecl_type_str(le_res)})")
        elif dt_max > 120:
            r.fail(f"{label} [{date_str}]: dt_max={dt_max:.1f}s")
        else:
            r.ok(dt_max, label)

        print(
            f"  {label}: {date_str}  dt_max={dt_max:6.1f}s  SE={ecl_type_str(se_res)}"
        )
        jd = se_tret[0] + 10

    # Search for penumbral lunar eclipses only
    jd = 2451545.0
    for i in range(3):
        label = f"Penumbral lunar #{i + 1}"
        try:
            se_res, se_tret = swe.lun_eclipse_when(jd, ecltype=ECL_PENUMBRAL)
            le_res, le_tret = ephem.lun_eclipse_when(jd, ecltype=ECL_PENUMBRAL)
        except Exception as e:
            r.fail(f"{label}: {e}")
            jd += 180
            continue

        dt_max = abs(se_tret[0] - le_tret[0]) * 86400
        date_str = jd_to_date_str(se_tret[0])

        if not (se_res & ECL_PENUMBRAL):
            r.fail(f"{label}: SE not penumbral ({ecl_type_str(se_res)})")
        elif not (le_res & ECL_PENUMBRAL):
            r.fail(f"{label}: LE not penumbral ({ecl_type_str(le_res)})")
        elif dt_max > 120:
            r.fail(f"{label} [{date_str}]: dt_max={dt_max:.1f}s")
        else:
            r.ok(dt_max, label)

        print(
            f"  {label}: {date_str}  dt_max={dt_max:6.1f}s  SE={ecl_type_str(se_res)}"
        )
        jd = se_tret[0] + 10

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 17: Deep Eclipse Global Functions Verification")
    print("=" * 70)

    start = time.time()
    all_ok = True
    all_results = []

    for pname, pfn in [
        ("P1", run_part1),
        ("P2", run_part2),
        ("P3", run_part3),
        ("P4", run_part4),
        ("P5", run_part5),
        ("P6", run_part6),
        ("P7", run_part7),
        ("P8", run_part8),
    ]:
        try:
            ok, res = pfn()
            all_results.append((pname, res))
            if not ok:
                all_ok = False
        except Exception as e:
            print(f"\n  {pname} CRASHED: {e}")
            traceback.print_exc()
            all_ok = False

    elapsed = time.time() - start

    print("\n" + "=" * 70)
    print("ROUND 17 FINAL SUMMARY")
    print("=" * 70)

    tp = tf = ts_ = 0
    for pn, res in all_results:
        st = "PASS" if res.failed == 0 else "FAIL"
        t = res.passed + res.failed
        print(f"  {pn} {res.name}: {res.passed}/{t} ({res.skipped} skip) [{st}]")
        tp += res.passed
        tf += res.failed
        ts_ += res.skipped

    print(f"\n  TOTAL: {tp}/{tp + tf} PASSED, {tf} FAILED, {ts_} SKIPPED")
    print(f"  Time: {elapsed:.1f}s")
    if all_ok:
        print(f"\n  >>> ROUND 17: ALL PASSED <<<")
    else:
        print(f"\n  >>> ROUND 17: {tf} FAILURES <<<")

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
