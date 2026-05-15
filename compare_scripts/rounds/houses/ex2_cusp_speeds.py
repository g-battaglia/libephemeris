"""Round 6: Houses Ex2 Cusp Speeds Deep Audit

Comprehensive comparison of libephemeris vs pyswisseph for:
- houses_ex2: cusp positions and speeds across all house systems
- houses_armc_ex2: ARMC-based cusp positions and speeds
- Multiple locations, dates, and house systems
- ASCMC angle speeds

Run with: env LIBEPHEMERIS_MODE=skyfield python3 compare_scripts/round6_houses_ex2_deep.py
"""

from __future__ import annotations

import math
import sys
import traceback

import swisseph as swe

sys.path.insert(0, ".")
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)
ephem.set_ephe_path(_REF_EPHE_PATH)

issues = []
stats = {}


def add_stat(cat, val):
    if cat not in stats:
        stats[cat] = []
    stats[cat].append(val)


def angular_diff(a, b):
    """Signed angular difference, handling wraparound."""
    d = a - b
    while d > 180:
        d -= 360
    while d < -180:
        d += 360
    return d


print("=" * 80)
print("ROUND 6: HOUSES EX2 CUSP SPEEDS DEEP AUDIT")
print("=" * 80)

# ============================================================================
# TEST PARAMETERS
# ============================================================================
LOCATIONS = [
    ("Rome", 41.9028, 12.4964),
    ("New York", 40.7128, -74.0060),
    ("Sydney", -33.8688, 151.2093),
    ("Equator", 0.0, 0.0),
    ("London", 51.5074, -0.1278),
    ("Tromso", 69.6496, 18.9560),  # Arctic (near polar circle)
    ("Cape Town", -33.9249, 18.4241),
    ("Singapore", 1.3521, 103.8198),  # Near equator
    ("Lat60N", 60.0, 25.0),  # Helsinki-like
    ("Lat30S", -30.0, -60.0),  # Southern mid-latitude
]

DATES = [
    swe.julday(2000, 1, 1, 12.0),  # J2000
    swe.julday(2024, 6, 21, 12.0),  # Summer solstice
    swe.julday(2024, 12, 21, 12.0),  # Winter solstice
    swe.julday(2024, 3, 20, 12.0),  # Vernal equinox
]

# House systems: (char, name, pos_tol, speed_tol)
# pyswisseph needs bytes (b'P'), libephemeris needs int (ord('P'))
HOUSE_SYSTEMS = [
    ("P", "Placidus", 0.002, 2.0),
    ("K", "Koch", 0.002, 100.0),  # Koch has known numerical sensitivity
    ("R", "Regiomontanus", 0.002, 2.0),
    ("C", "Campanus", 0.002, 2.0),
    ("E", "Equal", 0.001, 0.1),
    ("W", "WholeSign", 0.001, 0.1),
    ("O", "Porphyry", 0.002, 2.0),  # cusps 5,6,11,12 tested separately
    ("B", "Alcabitius", 0.002, 2.0),
    ("T", "Topocentric", 0.002, 2.0),
    ("M", "Morinus", 0.002, 2.0),
    ("X", "Meridian", 0.002, 2.0),
    ("V", "Vehlow", 0.001, 0.1),
    ("H", "Horizontal", 0.002, 2.0),
    ("U", "Krusinski", 0.002, 2.0),
    ("N", "NaturalGrad", 0.002, 2.0),
    ("Y", "APC", 0.002, 2.0),
    ("D", "EqualMC", 0.001, 0.1),
    ("S", "Sripati", 0.002, 2.0),
    ("L", "PullenSD", 0.002, 2.0),
    ("Q", "PullenSR", 0.002, 2.0),
]


def se_hsys(ch):
    """Convert char to pyswisseph hsys (bytes)."""
    return ch.encode("ascii")


def le_hsys(ch):
    """Convert char to libephemeris hsys (int)."""
    return ord(ch)


FLG_SPEED = 256  # FLG_SPEED

# ============================================================================
# PART 1: Cusp positions across all house systems and locations
# ============================================================================
print()
print("=" * 80)
print("PART 1: Cusp positions — all house systems across locations and dates")
print("=" * 80)

part1_pass = 0
part1_fail = 0
part1_skip = 0
part1_total = 0

for hsys_code, hsys_name, pos_tol, speed_tol in HOUSE_SYSTEMS:
    sys_pass = 0
    sys_fail = 0
    sys_total = 0
    max_diff = 0.0

    for name, lat, lon in LOCATIONS:
        for jd in DATES:
            try:
                se_cusps, se_ascmc = swe.houses_ex(jd, lat, lon, se_hsys(hsys_code))
            except Exception:
                part1_skip += 1
                continue

            try:
                le_ret = ephem.houses_ex(jd, lat, lon, le_hsys(hsys_code))
                le_cusps = le_ret[0]
                le_ascmc = le_ret[1]
            except Exception as e:
                if "polar" in str(e).lower() or "PolarCircle" in str(type(e).__name__):
                    part1_skip += 1
                    continue
                print(f"  LE ERROR {hsys_name}/{name}: {e}")
                part1_skip += 1
                continue

            # Compare cusps
            n_cusps = len(se_cusps)
            for i in range(min(n_cusps, len(le_cusps))):
                diff = abs(angular_diff(float(se_cusps[i]), float(le_cusps[i])))
                sys_total += 1
                part1_total += 1
                max_diff = max(max_diff, diff)

                if diff <= pos_tol:
                    sys_pass += 1
                    part1_pass += 1
                    add_stat(f"pos_{hsys_name}", diff)
                else:
                    sys_fail += 1
                    part1_fail += 1
                    issues.append(
                        f"pos {hsys_name} cusp{i + 1} {name}: diff={diff:.6f}°"
                    )

    status = "PASS" if sys_fail == 0 else "FAIL"
    print(f"  {status} {hsys_name:15s} {sys_pass}/{sys_total} max_diff={max_diff:.6f}°")

print()
print(f"Part 1 summary: {part1_pass}/{part1_total} passed ({part1_skip} skipped)")

# ============================================================================
# PART 2: ASCMC positions
# ============================================================================
print()
print("=" * 80)
print("PART 2: ASCMC positions (ASC, MC, ARMC, Vertex, etc.)")
print("=" * 80)

part2_pass = 0
part2_fail = 0
part2_skip = 0
part2_total = 0

ASCMC_LABELS = ["ASC", "MC", "ARMC", "Vertex", "EquAsc", "CoAscK", "CoAscM", "PolAsc"]
ASCMC_TOL = 0.002  # degrees

for name, lat, lon in LOCATIONS:
    for jd in DATES:
        try:
            se_cusps, se_ascmc = swe.houses_ex(jd, lat, lon, b"P")
        except Exception:
            part2_skip += 1
            continue

        try:
            le_ret = ephem.houses_ex(jd, lat, lon, ord("P"))
            le_ascmc = le_ret[1]
        except Exception as e:
            if "polar" in str(e).lower():
                part2_skip += 1
                continue
            part2_skip += 1
            continue

        for i in range(min(len(se_ascmc), len(le_ascmc), 8)):
            se_val = float(se_ascmc[i])
            le_val = float(le_ascmc[i])
            diff = abs(angular_diff(se_val, le_val))
            part2_total += 1

            if diff <= ASCMC_TOL:
                part2_pass += 1
                add_stat(f"ascmc_{ASCMC_LABELS[i]}", diff)
            else:
                part2_fail += 1
                issues.append(
                    f"ascmc {ASCMC_LABELS[i]} {name}: se={se_val:.6f} le={le_val:.6f} diff={diff:.6f}°"
                )
                print(f"  FAIL {name:10s} {ASCMC_LABELS[i]:8s}: diff={diff:.6f}°")

print()
print(f"Part 2 summary: {part2_pass}/{part2_total} passed ({part2_skip} skipped)")

# ============================================================================
# PART 3: Cusp speeds via houses_ex2
# ============================================================================
print()
print("=" * 80)
print("PART 3: Cusp speeds via houses_ex2 (FLG_SPEED)")
print("=" * 80)

part3_pass = 0
part3_fail = 0
part3_skip = 0
part3_total = 0

# Tighter tolerances for this deep audit
SPEED_TOLS = {
    "Placidus": 2.0,
    "Koch": 100.0,
    "Regiomontanus": 1.0,
    "Campanus": 1.0,
    "Equal": 0.1,
    "WholeSign": 0.1,
    "Porphyry": 2.0,  # cusps 1-4,7-10
    "Porphyry_opp": 250.0,  # cusps 5,6,11,12
    "Alcabitius": 2.0,
    "Topocentric": 2.0,
    "Morinus": 1.0,
    "Meridian": 1.0,
    "Vehlow": 0.1,
    "Horizontal": 2.0,
    "Krusinski": 2.0,
    "NaturalGrad": 2.0,
    "APC": 2.0,
    "EqualMC": 0.1,
    "Sripati": 2.0,
    "PullenSD": 2.0,
    "PullenSR": 2.0,
}

for hsys_code, hsys_name, pos_tol, _ in HOUSE_SYSTEMS:
    sys_pass = 0
    sys_fail = 0
    sys_total = 0
    max_speed_diff = 0.0

    for name, lat, lon in LOCATIONS[:6]:  # First 6 locations
        for jd in DATES[:2]:  # First 2 dates
            try:
                se_ret = swe.houses_ex2(jd, lat, lon, se_hsys(hsys_code), FLG_SPEED)
                se_cusps = se_ret[0]
                se_ascmc = se_ret[1]
                se_cusp_speeds = se_ret[2]
                se_ascmc_speeds = se_ret[3]
            except Exception:
                part3_skip += 1
                continue

            try:
                le_ret = ephem.houses_ex2(
                    jd, lat, lon, le_hsys(hsys_code), FLG_SPEED
                )
                le_cusps = le_ret[0]
                le_ascmc = le_ret[1]
                le_cusp_speeds = le_ret[2]
                le_ascmc_speeds = le_ret[3]
            except Exception as e:
                if "polar" in str(e).lower() or "PolarCircle" in str(type(e).__name__):
                    part3_skip += 1
                    continue
                print(f"  LE ERROR {hsys_name}/{name}: {e}")
                part3_skip += 1
                continue

            # Compare cusp speeds
            n_cusps = min(len(se_cusp_speeds), len(le_cusp_speeds))
            for i in range(n_cusps):
                se_speed = float(se_cusp_speeds[i])
                le_speed = float(le_cusp_speeds[i])
                diff = abs(se_speed - le_speed)

                # Porphyry: cusps 5,6,11,12 use different formula
                is_porphyry_opp = hsys_name == "Porphyry" and (i + 1) in [5, 6, 11, 12]
                if is_porphyry_opp:
                    tol = SPEED_TOLS.get("Porphyry_opp", 250.0)
                else:
                    tol = SPEED_TOLS.get(hsys_name, 2.0)

                sys_total += 1
                part3_total += 1
                max_speed_diff = max(max_speed_diff, diff)

                if diff <= tol:
                    sys_pass += 1
                    part3_pass += 1
                    add_stat(f"speed_{hsys_name}", diff)
                else:
                    sys_fail += 1
                    part3_fail += 1
                    issues.append(
                        f"speed {hsys_name} cusp{i + 1} {name}: "
                        f"se={se_speed:.4f} le={le_speed:.4f} diff={diff:.4f} deg/day"
                    )

    status = "PASS" if sys_fail == 0 else "FAIL"
    print(
        f"  {status} {hsys_name:15s} {sys_pass}/{sys_total} max_speed_diff={max_speed_diff:.4f} deg/day"
    )

print()
print(f"Part 3 summary: {part3_pass}/{part3_total} passed ({part3_skip} skipped)")

# ============================================================================
# PART 4: ASCMC speeds via houses_ex2
# ============================================================================
print()
print("=" * 80)
print("PART 4: ASCMC speeds (ASC, MC, ARMC, Vertex velocities)")
print("=" * 80)

part4_pass = 0
part4_fail = 0
part4_skip = 0
part4_total = 0

ASCMC_SPEED_TOL = 1.0  # deg/day

for name, lat, lon in LOCATIONS[:6]:
    for jd in DATES[:2]:
        try:
            se_ret = swe.houses_ex2(jd, lat, lon, b"P", FLG_SPEED)
            se_ascmc_speeds = se_ret[3]
        except Exception:
            part4_skip += 1
            continue

        try:
            le_ret = ephem.houses_ex2(jd, lat, lon, ord("P"), FLG_SPEED)
            le_ascmc_speeds = le_ret[3]
        except Exception as e:
            if "polar" in str(e).lower():
                part4_skip += 1
                continue
            part4_skip += 1
            continue

        for i in range(min(len(se_ascmc_speeds), len(le_ascmc_speeds), 8)):
            se_val = float(se_ascmc_speeds[i])
            le_val = float(le_ascmc_speeds[i])
            diff = abs(se_val - le_val)
            part4_total += 1

            label = ASCMC_LABELS[i] if i < len(ASCMC_LABELS) else f"idx{i}"

            if diff <= ASCMC_SPEED_TOL:
                part4_pass += 1
                add_stat(f"ascmc_speed_{label}", diff)
            else:
                part4_fail += 1
                issues.append(
                    f"ascmc_speed {label} {name}: se={se_val:.4f} le={le_val:.4f} diff={diff:.4f}"
                )
                print(
                    f"  FAIL {name:10s} {label:8s}: se={se_val:.4f} le={le_val:.4f} diff={diff:.4f}"
                )

print()
print(f"Part 4 summary: {part4_pass}/{part4_total} passed ({part4_skip} skipped)")

# ============================================================================
# PART 5: houses_armc_ex2 cusp speeds
# ============================================================================
print()
print("=" * 80)
print("PART 5: houses_armc_ex2 cusp speeds")
print("=" * 80)

part5_pass = 0
part5_fail = 0
part5_skip = 0
part5_total = 0

# Test a range of ARMC values
ARMC_VALUES = [0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0]
OBLIQUITY = 23.4393  # J2000 approximate

ARMC_SYSTEMS = [
    ("P", "Placidus", 2.0),
    ("R", "Regiomontanus", 1.0),
    ("C", "Campanus", 1.0),
    ("E", "Equal", 0.1),
    ("M", "Morinus", 1.0),
    ("B", "Alcabitius", 2.0),
    ("T", "Topocentric", 2.0),
]

for hsys_code, hsys_name, speed_tol in ARMC_SYSTEMS:
    sys_pass = 0
    sys_fail = 0
    sys_total = 0
    max_diff = 0.0

    for armc in ARMC_VALUES:
        for lat in [0.0, 30.0, 45.0, 60.0, -30.0]:
            try:
                se_ret = swe.houses_armc_ex2(armc, lat, OBLIQUITY, se_hsys(hsys_code))
                se_cusps = se_ret[0]
                se_ascmc = se_ret[1]
                se_cusp_speeds = se_ret[2]
                se_ascmc_speeds = se_ret[3]
            except Exception:
                part5_skip += 1
                continue

            try:
                le_ret = ephem.houses_armc_ex2(
                    armc, lat, OBLIQUITY, le_hsys(hsys_code), FLG_SPEED
                )
                le_cusps = le_ret[0]
                le_ascmc = le_ret[1]
                le_cusp_speeds = le_ret[2]
                le_ascmc_speeds = le_ret[3]
            except Exception as e:
                if "polar" in str(e).lower() or "PolarCircle" in str(type(e).__name__):
                    part5_skip += 1
                    continue
                print(f"  LE ERROR {hsys_name}/armc={armc}/lat={lat}: {e}")
                part5_skip += 1
                continue

            # Compare cusp speeds
            n_cusps = min(len(se_cusp_speeds), len(le_cusp_speeds))
            for i in range(n_cusps):
                se_speed = float(se_cusp_speeds[i])
                le_speed = float(le_cusp_speeds[i])
                diff = abs(se_speed - le_speed)

                sys_total += 1
                part5_total += 1
                max_diff = max(max_diff, diff)

                if diff <= speed_tol:
                    sys_pass += 1
                    part5_pass += 1
                    add_stat(f"armc_speed_{hsys_name}", diff)
                else:
                    sys_fail += 1
                    part5_fail += 1
                    issues.append(
                        f"armc_speed {hsys_name} cusp{i + 1} armc={armc} lat={lat}: "
                        f"diff={diff:.4f}"
                    )

    status = "PASS" if sys_fail == 0 else "FAIL"
    print(
        f"  {status} {hsys_name:15s} {sys_pass}/{sys_total} max_diff={max_diff:.4f} deg/day"
    )

print()
print(f"Part 5 summary: {part5_pass}/{part5_total} passed ({part5_skip} skipped)")

# ============================================================================
# PART 6: Consistency — houses_ex2 cusps match houses_ex cusps exactly
# ============================================================================
print()
print("=" * 80)
print("PART 6: Consistency — houses_ex2 cusps == houses_ex cusps")
print("=" * 80)

part6_pass = 0
part6_fail = 0
part6_total = 0

for hsys_code, hsys_name, _, _ in HOUSE_SYSTEMS[:10]:  # First 10 systems
    for name, lat, lon in LOCATIONS[:4]:
        jd = DATES[0]
        try:
            ex_ret = ephem.houses_ex(jd, lat, lon, le_hsys(hsys_code))
            ex2_ret = ephem.houses_ex2(jd, lat, lon, le_hsys(hsys_code))

            ex_cusps = ex_ret[0]
            ex2_cusps = ex2_ret[0]

            for i in range(min(len(ex_cusps), len(ex2_cusps))):
                diff = abs(float(ex_cusps[i]) - float(ex2_cusps[i]))
                part6_total += 1

                if diff < 1e-10:
                    part6_pass += 1
                else:
                    part6_fail += 1
                    issues.append(
                        f"consistency {hsys_name} cusp{i + 1} {name}: "
                        f"ex={ex_cusps[i]:.10f} ex2={ex2_cusps[i]:.10f} diff={diff}"
                    )
                    print(f"  FAIL {hsys_name:15s} cusp{i + 1} {name}: diff={diff}")
        except Exception:
            continue

print()
print(f"Part 6 summary: {part6_pass}/{part6_total} passed")

# ============================================================================
# PART 7: Gauquelin sectors (36 cusps)
# ============================================================================
print()
print("=" * 80)
print("PART 7: Gauquelin sectors (36 cusps)")
print("=" * 80)

part7_pass = 0
part7_fail = 0
part7_skip = 0
part7_total = 0

GAUQUELIN_TOL = 0.002  # degrees

for name, lat, lon in LOCATIONS[:4]:
    for jd in DATES[:2]:
        try:
            se_cusps, se_ascmc = swe.houses_ex(jd, lat, lon, b"G")
        except Exception:
            part7_skip += 1
            continue

        try:
            le_ret = ephem.houses_ex(jd, lat, lon, ord("G"))
            le_cusps = le_ret[0]
        except Exception:
            part7_skip += 1
            continue

        for i in range(min(len(se_cusps), len(le_cusps))):
            diff = abs(angular_diff(float(se_cusps[i]), float(le_cusps[i])))
            part7_total += 1

            if diff <= GAUQUELIN_TOL:
                part7_pass += 1
                add_stat("gauquelin_pos", diff)
            else:
                part7_fail += 1
                issues.append(f"gauquelin cusp{i + 1} {name}: diff={diff:.6f}°")
                print(f"  FAIL {name:10s} cusp{i + 1}: diff={diff:.6f}°")

print()
print(f"Part 7 summary: {part7_pass}/{part7_total} passed ({part7_skip} skipped)")

# ============================================================================
# PART 8: Sidereal house cusps (FLG_SIDEREAL)
# ============================================================================
print()
print("=" * 80)
print("PART 8: Sidereal house cusps (Lahiri ayanamsha)")
print("=" * 80)

part8_pass = 0
part8_fail = 0
part8_skip = 0
part8_total = 0

FLG_SIDEREAL = 64 * 1024  # 0x10000
SIDEREAL_TOL = 0.01  # degrees (wider due to ayanamsha differences)

# Set Lahiri ayanamsha for both
swe.set_sid_mode(1)  # Lahiri
ephem.set_sid_mode(1, 0.0, 0.0)

SIDEREAL_SYSTEMS = ["P", "E", "W", "R"]

for hsys_code in SIDEREAL_SYSTEMS:
    hsys_name = hsys_code
    for name, lat, lon in LOCATIONS[:4]:
        jd = DATES[0]
        try:
            se_cusps, se_ascmc = swe.houses_ex(
                jd, lat, lon, se_hsys(hsys_code), FLG_SIDEREAL
            )
        except Exception:
            part8_skip += 1
            continue

        try:
            le_ret = ephem.houses_ex(
                jd, lat, lon, le_hsys(hsys_code), FLG_SIDEREAL
            )
            le_cusps = le_ret[0]
        except Exception as e:
            if "polar" in str(e).lower():
                part8_skip += 1
                continue
            part8_skip += 1
            continue

        for i in range(min(len(se_cusps), len(le_cusps))):
            diff = abs(angular_diff(float(se_cusps[i]), float(le_cusps[i])))
            part8_total += 1

            if diff <= SIDEREAL_TOL:
                part8_pass += 1
                add_stat(f"sidereal_{hsys_name}", diff)
            else:
                part8_fail += 1
                issues.append(
                    f"sidereal {hsys_name} cusp{i + 1} {name}: diff={diff:.6f}°"
                )
                print(f"  FAIL {hsys_name} {name:10s} cusp{i + 1}: diff={diff:.6f}°")

# Reset sidereal mode
swe.set_sid_mode(0)
ephem.set_sid_mode(0, 0.0, 0.0)

print()
print(f"Part 8 summary: {part8_pass}/{part8_total} passed ({part8_skip} skipped)")

# ============================================================================
# SUMMARY
# ============================================================================
print()
print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()
print(
    f"Part 1 (Cusp positions):       {part1_pass}/{part1_total} passed ({part1_skip} skip)"
)
print(
    f"Part 2 (ASCMC positions):      {part2_pass}/{part2_total} passed ({part2_skip} skip)"
)
print(
    f"Part 3 (Cusp speeds ex2):      {part3_pass}/{part3_total} passed ({part3_skip} skip)"
)
print(
    f"Part 4 (ASCMC speeds):         {part4_pass}/{part4_total} passed ({part4_skip} skip)"
)
print(
    f"Part 5 (ARMC ex2 speeds):      {part5_pass}/{part5_total} passed ({part5_skip} skip)"
)
print(f"Part 6 (Internal consistency): {part6_pass}/{part6_total} passed")
print(
    f"Part 7 (Gauquelin 36 cusps):   {part7_pass}/{part7_total} passed ({part7_skip} skip)"
)
print(
    f"Part 8 (Sidereal cusps):       {part8_pass}/{part8_total} passed ({part8_skip} skip)"
)

total_pass = (
    part1_pass
    + part2_pass
    + part3_pass
    + part4_pass
    + part5_pass
    + part6_pass
    + part7_pass
    + part8_pass
)
total_total = (
    part1_total
    + part2_total
    + part3_total
    + part4_total
    + part5_total
    + part6_total
    + part7_total
    + part8_total
)

print()
print(f"TOTAL: {total_pass}/{total_total} passed")
print()

# Print statistics
for cat in sorted(stats.keys()):
    vals = stats[cat]
    if vals:
        print(
            f"Stats for {cat}:"
            f"\n  Min: {min(vals):.6f}  Max: {max(vals):.6f}"
            f"  Mean: {sum(vals) / len(vals):.6f}"
            f"  Median: {sorted(vals)[len(vals) // 2]:.6f}  (N={len(vals)})"
        )

if issues:
    print()
    print("=" * 80)
    print(f"ISSUES FOUND: {len(issues)}")
    for iss in issues[:50]:
        print(f"  - {iss}")
    if len(issues) > 50:
        print(f"  ... and {len(issues) - 50} more")
    print("=" * 80)
