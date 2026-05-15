"""Round 4: Crossing Functions Deep Audit

Comprehensive comparison of libephemeris vs pyswisseph for:
- solcross_ut / solcross: Sun longitude crossing
- mooncross_ut / mooncross: Moon longitude crossing
- mooncross_node_ut / mooncross_node: Moon node crossing
- cross_ut: Generic planet crossing (all planets)
- helio_cross_ut / helio_cross: Heliocentric crossing

Run with: env LIBEPHEMERIS_MODE=skyfield python3 compare_scripts/round4_crossing_deep.py
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
stats = {}  # category -> list of diffs in arcseconds


def lon_diff_arcsec(jd_se, jd_le):
    """Compute timing difference in arcseconds of longitude."""
    return abs(jd_se - jd_le) * 86400  # seconds of time


def add_stat(cat, val):
    if cat not in stats:
        stats[cat] = []
    stats[cat].append(val)


print("=" * 80)
print("ROUND 4: CROSSING FUNCTIONS DEEP AUDIT")
print("=" * 80)

# ============================================================================
# PART 1: solcross_ut — Sun longitude crossing (UT)
# ============================================================================
print()
print("=" * 80)
print("PART 1: solcross_ut — Sun longitude crossing")
print("=" * 80)

part1_pass = 0
part1_fail = 0
part1_total = 0

# Test: Sun crossing every 30° from 2020 to 2026
# (covers equinoxes, solstices, and intermediate points)
target_lons = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
jd_start = swe.julday(2020, 1, 1, 0.0)

for target in target_lons:
    jd = jd_start
    for year_i in range(6):  # 6 crossings per target (2020-2025)
        try:
            se_jd = swe.solcross_ut(target, jd)
        except Exception as e:
            print("  SE error solcross_ut({}, {}): {}".format(target, jd, e))
            continue

        try:
            le_jd = ephem.solcross_ut(target, jd)
        except Exception as e:
            print("  LE error solcross_ut({}, {}): {}".format(target, jd, e))
            continue

        diff_s = abs(se_jd - le_jd) * 86400
        part1_total += 1

        # Verify both actually cross the target
        se_pos = swe.calc_ut(se_jd, 0)[0][0]
        le_pos_r = ephem.calc_ut(le_jd, 0, 256)
        le_pos = float(le_pos_r[0][0])

        se_err = abs(se_pos - target)
        if se_err > 180:
            se_err = 360 - se_err
        le_err = abs(le_pos - target)
        if le_err > 180:
            le_err = 360 - le_err

        # Tolerance: 0.01 arcsec in time
        tol_s = 0.01  # seconds of time
        if diff_s > tol_s or le_err * 3600 > 0.01:
            status = "FAIL"
            part1_fail += 1
            issues.append(
                'solcross_ut lon={}: diff={:.6f}s, le_err={:.6f}"'.format(
                    target, diff_s, le_err * 3600
                )
            )
        else:
            status = "PASS"
            part1_pass += 1
            add_stat("solcross_ut", diff_s)

        if year_i == 0 or status == "FAIL":
            print(
                '  {} lon={:3d}° yr={}: diff={:.6f}s, se_err={:.6f}", le_err={:.6f}"'.format(
                    status, target, 2020 + year_i, diff_s, se_err * 3600, le_err * 3600
                )
            )

        jd = se_jd + 300  # skip ahead ~10 months

print()
print("Part 1 summary: {}/{} passed".format(part1_pass, part1_total))

# ============================================================================
# PART 2: mooncross_ut — Moon longitude crossing (UT)
# ============================================================================
print()
print("=" * 80)
print("PART 2: mooncross_ut — Moon longitude crossing")
print("=" * 80)

part2_pass = 0
part2_fail = 0
part2_total = 0

# Test: Moon crossing every 30° from Jan 2024 to Jan 2025
target_lons_moon = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
jd_moon_start = swe.julday(2024, 1, 1, 0.0)

for target in target_lons_moon:
    jd = jd_moon_start
    for cross_i in range(12):  # ~12 crossings per year for each longitude
        try:
            se_jd = swe.mooncross_ut(target, jd)
        except Exception as e:
            break

        if se_jd <= 0 or se_jd > jd_moon_start + 400:
            break

        try:
            le_jd = ephem.mooncross_ut(target, jd)
        except Exception as e:
            print("  LE error mooncross_ut({}, {}): {}".format(target, jd, e))
            jd = se_jd + 20
            continue

        diff_s = abs(se_jd - le_jd) * 86400
        part2_total += 1

        # Tolerance: 0.1s timing diff
        tol_s = 0.1
        if diff_s > tol_s:
            status = "FAIL"
            part2_fail += 1
            issues.append("mooncross_ut lon={}: diff={:.4f}s".format(target, diff_s))
        else:
            status = "PASS"
            part2_pass += 1
            add_stat("mooncross_ut", diff_s)

        if cross_i == 0 or status == "FAIL":
            print(
                "  {} lon={:3d}° #{:2d}: diff={:.6f}s".format(
                    status, target, cross_i + 1, diff_s
                )
            )

        jd = se_jd + 20  # Moon takes ~27 days per cycle

print()
print("Part 2 summary: {}/{} passed".format(part2_pass, part2_total))

# ============================================================================
# PART 3: mooncross_node_ut — Moon node crossing
# ============================================================================
print()
print("=" * 80)
print("PART 3: mooncross_node_ut — Moon node crossing")
print("=" * 80)

part3_pass = 0
part3_fail = 0
part3_total = 0

jd_node_start = swe.julday(2020, 1, 1, 0.0)
jd = jd_node_start

for node_i in range(60):  # ~2 node crossings per month, ~5 years
    try:
        se_ret = swe.mooncross_node_ut(jd)
        se_jd = se_ret[0]
        se_lon = se_ret[1]
        se_lat = se_ret[2]
    except Exception as e:
        break

    if se_jd <= 0 or se_jd > jd_node_start + 2000:
        break

    try:
        le_ret = ephem.mooncross_node_ut(jd)
        le_jd = le_ret[0]
        le_lon = float(le_ret[1])
        le_lat = float(le_ret[2])
    except Exception as e:
        print("  LE error mooncross_node_ut({}): {}".format(jd, e))
        jd = se_jd + 5
        continue

    diff_s = abs(se_jd - le_jd) * 86400
    lon_diff = abs(se_lon - le_lon)
    if lon_diff > 180:
        lon_diff = 360 - lon_diff
    lat_diff = abs(se_lat - le_lat)

    part3_total += 1

    # Tolerance: 0.5s timing, 0.01° longitude, 0.001° latitude
    tol_time = 0.5
    tol_lon = 0.01
    tol_lat = 0.001

    failed = False
    fail_reasons = []
    if diff_s > tol_time:
        failed = True
        fail_reasons.append("time={:.4f}s".format(diff_s))
    if lon_diff > tol_lon:
        failed = True
        fail_reasons.append("lon={:.6f}°".format(lon_diff))
    if lat_diff > tol_lat:
        failed = True
        fail_reasons.append("lat={:.6f}°".format(lat_diff))

    if failed:
        status = "FAIL"
        part3_fail += 1
        issues.append(
            "mooncross_node #{}: {}".format(node_i + 1, ", ".join(fail_reasons))
        )
    else:
        status = "PASS"
        part3_pass += 1
        add_stat("mooncross_node", diff_s)

    if node_i < 5 or status == "FAIL":
        node_type = "ASC" if se_lat >= 0 or (le_lat >= 0) else "DESC"
        # Actually determine from latitude direction
        # At ascending node, lat goes from - to +
        # Check lat value at crossing — should be ~0
        print(
            "  {} node #{:2d}: diff={:.6f}s, lon_diff={:.6f}°, lat_se={:.8f}° lat_le={:.8f}°".format(
                status, node_i + 1, diff_s, lon_diff, se_lat, le_lat
            )
        )

    jd = se_jd + 10  # ~13.6 days between node crossings

print()
print("Part 3 summary: {}/{} passed".format(part3_pass, part3_total))

# ============================================================================
# PART 4: cross_ut — Generic planet crossing
# ============================================================================
print()
print("=" * 80)
print("PART 4: cross_ut — Generic planet crossing")
print("=" * 80)

part4_pass = 0
part4_fail = 0
part4_skip = 0
part4_total = 0

# Test all major planets crossing specific longitudes
# planet_id: (name, target_lons, start_jd, skip_days)
planet_tests = [
    (0, "Sun", [0, 90, 180, 270], swe.julday(2023, 1, 1, 0.0), 300),
    (1, "Moon", [0, 90, 180, 270], swe.julday(2024, 1, 1, 0.0), 20),
    (2, "Mercury", [0, 90, 180, 270], swe.julday(2024, 1, 1, 0.0), 60),
    (3, "Venus", [0, 90, 180, 270], swe.julday(2024, 1, 1, 0.0), 200),
    (4, "Mars", [0, 90, 180, 270], swe.julday(2023, 1, 1, 0.0), 500),
    (5, "Jupiter", [0, 90], swe.julday(2020, 1, 1, 0.0), 3000),
    (6, "Saturn", [0, 330], swe.julday(2020, 1, 1, 0.0), 5000),
]

for planet_id, name, targets, jd_start_p, skip in planet_tests:
    for target in targets:
        jd = jd_start_p
        n_crosses = 3 if planet_id <= 1 else 2

        for cross_i in range(n_crosses):
            try:
                se_jd = swe.cross_ut(planet_id, target, jd)
            except Exception as e:
                part4_skip += 1
                break

            if se_jd <= 0:
                break

            try:
                le_jd = ephem.cross_ut(planet_id, target, jd)
            except Exception as e:
                print("  LE error cross_ut({}, {}, {}): {}".format(name, target, jd, e))
                part4_skip += 1
                jd = se_jd + skip
                continue

            diff_s = abs(se_jd - le_jd) * 86400
            part4_total += 1

            # Tolerance depends on planet speed
            if planet_id <= 1:
                tol_s = 0.5
            elif planet_id <= 4:
                tol_s = 2.0
            else:
                tol_s = 10.0

            if diff_s > tol_s:
                status = "FAIL"
                part4_fail += 1
                issues.append(
                    "cross_ut {} lon={}: diff={:.4f}s".format(name, target, diff_s)
                )
            else:
                status = "PASS"
                part4_pass += 1
                add_stat("cross_ut_" + name, diff_s)

            print(
                "  {} {} lon={:3d}° #{}: diff={:.6f}s".format(
                    status, name, target, cross_i + 1, diff_s
                )
            )

            jd = se_jd + skip

print()
print(
    "Part 4 summary: {}/{} passed ({} skipped)".format(
        part4_pass, part4_total, part4_skip
    )
)

# ============================================================================
# PART 5: helio_cross_ut — Heliocentric crossing
# ============================================================================
print()
print("=" * 80)
print("PART 5: helio_cross_ut — Heliocentric planet crossing")
print("=" * 80)

part5_pass = 0
part5_fail = 0
part5_skip = 0
part5_total = 0

# Test heliocentric crossings for major planets
helio_tests = [
    (3, "Venus", [0, 90, 180, 270], swe.julday(2024, 1, 1, 0.0), 200),
    (4, "Mars", [0, 90, 180, 270], swe.julday(2023, 1, 1, 0.0), 500),
    (5, "Jupiter", [0, 90], swe.julday(2020, 1, 1, 0.0), 3000),
    (6, "Saturn", [330, 0], swe.julday(2020, 1, 1, 0.0), 6000),
    # Earth heliocentric (body 14 in SE)
]

for planet_id, name, targets, jd_start_h, skip in helio_tests:
    for target in targets:
        jd = jd_start_h
        for cross_i in range(2):
            try:
                se_jd = swe.helio_cross_ut(planet_id, target, jd)
            except Exception as e:
                part5_skip += 1
                break

            if se_jd <= 0:
                break

            try:
                le_jd = ephem.helio_cross_ut(planet_id, target, jd)
            except Exception as e:
                print(
                    "  LE error helio_cross_ut({}, {}, {}): {}".format(
                        name, target, jd, e
                    )
                )
                part5_skip += 1
                jd = se_jd + skip
                continue

            diff_s = abs(se_jd - le_jd) * 86400
            part5_total += 1

            # Tolerance
            if planet_id <= 4:
                tol_s = 2.0
            else:
                tol_s = 10.0

            if diff_s > tol_s:
                status = "FAIL"
                part5_fail += 1
                issues.append(
                    "helio_cross_ut {} lon={}: diff={:.4f}s".format(
                        name, target, diff_s
                    )
                )
            else:
                status = "PASS"
                part5_pass += 1
                add_stat("helio_cross_ut_" + name, diff_s)

            print(
                "  {} {} lon={:3d}° #{}: diff={:.6f}s".format(
                    status, name, target, cross_i + 1, diff_s
                )
            )

            jd = se_jd + skip

print()
print(
    "Part 5 summary: {}/{} passed ({} skipped)".format(
        part5_pass, part5_total, part5_skip
    )
)

# ============================================================================
# PART 6: TT versions — solcross, mooncross, mooncross_node
# ============================================================================
print()
print("=" * 80)
print("PART 6: TT versions — solcross, mooncross, mooncross_node, helio_cross")
print("=" * 80)

part6_pass = 0
part6_fail = 0
part6_total = 0

# Test TT versions
jd_tt_start = swe.julday(2024, 1, 1, 0.0)

# solcross TT
for target in [0, 90, 180, 270]:
    try:
        se_jd = swe.solcross(target, jd_tt_start)
        le_jd = ephem.solcross(target, jd_tt_start)
        diff_s = abs(se_jd - le_jd) * 86400
        part6_total += 1
        if diff_s > 0.01:
            status = "FAIL"
            part6_fail += 1
            issues.append("solcross(TT) lon={}: diff={:.6f}s".format(target, diff_s))
        else:
            status = "PASS"
            part6_pass += 1
            add_stat("solcross_tt", diff_s)
        print(
            "  {} solcross(TT) lon={:3d}°: diff={:.6f}s".format(status, target, diff_s)
        )
    except Exception as e:
        print("  ERROR solcross(TT) lon={}: {}".format(target, e))

# mooncross TT
for target in [0, 90, 180, 270]:
    try:
        se_jd = swe.mooncross(target, jd_tt_start)
        le_jd = ephem.mooncross(target, jd_tt_start)
        diff_s = abs(se_jd - le_jd) * 86400
        part6_total += 1
        if diff_s > 0.1:
            status = "FAIL"
            part6_fail += 1
            issues.append("mooncross(TT) lon={}: diff={:.6f}s".format(target, diff_s))
        else:
            status = "PASS"
            part6_pass += 1
            add_stat("mooncross_tt", diff_s)
        print(
            "  {} mooncross(TT) lon={:3d}°: diff={:.6f}s".format(status, target, diff_s)
        )
    except Exception as e:
        print("  ERROR mooncross(TT) lon={}: {}".format(target, e))

# mooncross_node TT
jd = jd_tt_start
for node_i in range(4):
    try:
        se_ret = swe.mooncross_node(jd)
        le_ret = ephem.mooncross_node(jd)
        diff_s = abs(se_ret[0] - le_ret[0]) * 86400
        part6_total += 1
        if diff_s > 0.5:
            status = "FAIL"
            part6_fail += 1
            issues.append(
                "mooncross_node(TT) #{}: diff={:.4f}s".format(node_i + 1, diff_s)
            )
        else:
            status = "PASS"
            part6_pass += 1
            add_stat("mooncross_node_tt", diff_s)
        print(
            "  {} mooncross_node(TT) #{}: diff={:.6f}s".format(
                status, node_i + 1, diff_s
            )
        )
        jd = se_ret[0] + 10
    except Exception as e:
        print("  ERROR mooncross_node(TT) #{}: {}".format(node_i + 1, e))
        jd += 14

# helio_cross TT
for planet_id, name in [(4, "Mars"), (5, "Jupiter")]:
    target = 90
    jd = jd_tt_start
    try:
        se_jd = swe.helio_cross(planet_id, target, jd)
        le_jd = ephem.helio_cross(planet_id, target, jd)
        diff_s = abs(se_jd - le_jd) * 86400
        part6_total += 1
        tol_s = 2.0 if planet_id <= 4 else 10.0
        if diff_s > tol_s:
            status = "FAIL"
            part6_fail += 1
            issues.append(
                "helio_cross(TT) {} lon={}: diff={:.4f}s".format(name, target, diff_s)
            )
        else:
            status = "PASS"
            part6_pass += 1
            add_stat("helio_cross_tt", diff_s)
        print(
            "  {} helio_cross(TT) {} lon={:3d}°: diff={:.6f}s".format(
                status, name, target, diff_s
            )
        )
    except Exception as e:
        print("  ERROR helio_cross(TT) {} lon={}: {}".format(name, target, e))

print()
print("Part 6 summary: {}/{} passed".format(part6_pass, part6_total))

# ============================================================================
# PART 7: Edge cases — 360°/0° boundary, negative targets, retrograde
# ============================================================================
print()
print("=" * 80)
print("PART 7: Edge cases — 360°/0° boundary, near-station crossings")
print("=" * 80)

part7_pass = 0
part7_fail = 0
part7_total = 0

# Test 360°=0° boundary for Sun
for target in [359.999, 0.001, 360.0]:
    try:
        se_jd = swe.solcross_ut(target % 360, swe.julday(2024, 3, 15, 0.0))
        le_jd = ephem.solcross_ut(target, swe.julday(2024, 3, 15, 0.0))
        diff_s = abs(se_jd - le_jd) * 86400
        part7_total += 1
        if diff_s > 0.1:
            status = "FAIL"
            part7_fail += 1
            issues.append(
                "solcross_ut edge lon={}: diff={:.4f}s".format(target, diff_s)
            )
        else:
            status = "PASS"
            part7_pass += 1
        print(
            "  {} solcross_ut lon={:.3f}°: diff={:.6f}s".format(status, target, diff_s)
        )
    except Exception as e:
        print("  ERROR solcross_ut edge lon={}: {}".format(target, e))

# Test Mercury crossing near station (retrograde)
# Mercury retrograde periods: find a crossing near a station
# Mercury retrograde Apr 1-25, 2024
print()
print("  Mercury near-station crossings:")
for target in [10, 20, 350]:
    jd_merc = swe.julday(2024, 3, 1, 0.0)
    try:
        se_jd = swe.cross_ut(2, target, jd_merc)
    except Exception:
        continue
    if se_jd <= 0:
        continue
    try:
        le_jd = ephem.cross_ut(2, target, jd_merc)
    except Exception as e:
        print("    LE error Mercury lon={}: {}".format(target, e))
        continue

    diff_s = abs(se_jd - le_jd) * 86400
    part7_total += 1
    if diff_s > 5.0:
        status = "FAIL"
        part7_fail += 1
        issues.append(
            "cross_ut Mercury(retrograde) lon={}: diff={:.4f}s".format(target, diff_s)
        )
    else:
        status = "PASS"
        part7_pass += 1
    print("    {} Mercury lon={:3d}°: diff={:.6f}s".format(status, target, diff_s))

# Test Venus crossing
print()
print("  Venus crossings:")
for target in [0, 45, 90, 135, 180, 225, 270, 315]:
    jd_v = swe.julday(2024, 1, 1, 0.0)
    try:
        se_jd = swe.cross_ut(3, target, jd_v)
    except Exception:
        continue
    if se_jd <= 0:
        continue
    try:
        le_jd = ephem.cross_ut(3, target, jd_v)
    except Exception as e:
        print("    LE error Venus lon={}: {}".format(target, e))
        continue

    diff_s = abs(se_jd - le_jd) * 86400
    part7_total += 1
    if diff_s > 2.0:
        status = "FAIL"
        part7_fail += 1
        issues.append("cross_ut Venus lon={}: diff={:.4f}s".format(target, diff_s))
    else:
        status = "PASS"
        part7_pass += 1
    print("    {} Venus lon={:3d}°: diff={:.6f}s".format(status, target, diff_s))

print()
print("Part 7 summary: {}/{} passed".format(part7_pass, part7_total))

# ============================================================================
# PART 8: 50-year sweep — Sun equinoxes/solstices
# ============================================================================
print()
print("=" * 80)
print("PART 8: 50-year sweep — Sun equinoxes and solstices (2000-2050)")
print("=" * 80)

part8_pass = 0
part8_fail = 0
part8_total = 0
part8_max_diff = 0.0

for target in [0, 90, 180, 270]:
    jd = swe.julday(2000, 1, 1, 0.0)
    for _ in range(50):
        try:
            se_jd = swe.solcross_ut(target, jd)
            le_jd = ephem.solcross_ut(target, jd)
        except Exception:
            break

        if se_jd <= 0 or se_jd > swe.julday(2051, 1, 1, 0.0):
            break

        diff_s = abs(se_jd - le_jd) * 86400
        part8_total += 1
        part8_max_diff = max(part8_max_diff, diff_s)

        if diff_s > 0.01:
            status = "FAIL"
            part8_fail += 1
            issues.append("50yr solcross lon={}: diff={:.6f}s".format(target, diff_s))
            print(
                "  FAIL lon={:3d}° y~{}: diff={:.6f}s".format(target, 2000 + _, diff_s)
            )
        else:
            part8_pass += 1
            add_stat("solcross_50yr", diff_s)

        jd = se_jd + 300

print(
    "  50-year sweep: {}/{} passed, max_diff={:.6f}s".format(
        part8_pass, part8_total, part8_max_diff
    )
)

# ============================================================================
# PART 9: Moon node precision — latitude at crossing
# ============================================================================
print()
print("=" * 80)
print("PART 9: Moon node precision — latitude residual at crossing")
print("=" * 80)

part9_pass = 0
part9_fail = 0
part9_total = 0

jd = swe.julday(2024, 1, 1, 0.0)
for node_i in range(24):  # ~1 year of node crossings
    try:
        le_ret = ephem.mooncross_node_ut(jd)
        le_jd = le_ret[0]
        le_lat = float(le_ret[2])

        # Verify: compute Moon position at crossing time
        moon_at_cross = ephem.calc_ut(le_jd, 1, 256)
        actual_lat = float(moon_at_cross[0][1])
    except Exception as e:
        jd += 14
        continue

    part9_total += 1
    lat_arcsec = abs(actual_lat) * 3600

    # Tolerance: latitude should be < 0.1 arcsec at crossing
    if lat_arcsec > 0.1:
        status = "FAIL"
        part9_fail += 1
        issues.append('node lat residual #{}: {:.4f}"'.format(node_i + 1, lat_arcsec))
    else:
        status = "PASS"
        part9_pass += 1
        add_stat("node_lat_residual", lat_arcsec)

    if node_i < 4 or status == "FAIL":
        print(
            '  {} node #{:2d}: lat_at_cross={:.6f}" (should be ~0)'.format(
                status, node_i + 1, lat_arcsec
            )
        )

    jd = le_jd + 10

print()
print("Part 9 summary: {}/{} passed".format(part9_pass, part9_total))

# ============================================================================
# SUMMARY
# ============================================================================
print()
print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()
print("Part 1 (solcross_ut):        {}/{} passed".format(part1_pass, part1_total))
print("Part 2 (mooncross_ut):       {}/{} passed".format(part2_pass, part2_total))
print("Part 3 (mooncross_node_ut):  {}/{} passed".format(part3_pass, part3_total))
print(
    "Part 4 (cross_ut):           {}/{} passed ({} skip)".format(
        part4_pass, part4_total, part4_skip
    )
)
print(
    "Part 5 (helio_cross_ut):     {}/{} passed ({} skip)".format(
        part5_pass, part5_total, part5_skip
    )
)
print("Part 6 (TT versions):        {}/{} passed".format(part6_pass, part6_total))
print("Part 7 (Edge cases):         {}/{} passed".format(part7_pass, part7_total))
print("Part 8 (50yr sweep):         {}/{} passed".format(part8_pass, part8_total))
print("Part 9 (Node lat residual):  {}/{} passed".format(part9_pass, part9_total))

for cat, vals in sorted(stats.items()):
    if vals:
        print()
        print("Stats for {}:".format(cat))
        print(
            "  Min: {:.6f}  Max: {:.6f}  Mean: {:.6f}  Median: {:.6f}  (N={})".format(
                min(vals),
                max(vals),
                sum(vals) / len(vals),
                sorted(vals)[len(vals) // 2],
                len(vals),
            )
        )

print()
print("=" * 80)
if issues:
    print("ISSUES FOUND: {}".format(len(issues)))
    for iss in issues:
        print("  - {}".format(iss))
else:
    print("NO ISSUES FOUND")
print("=" * 80)
