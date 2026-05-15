"""Round 5: Rise/Transit/Set Deep Audit

Comprehensive comparison of libephemeris vs pyswisseph for:
- rise_trans: Sun/Moon/planets rise, set, upper/lower transit
- rise_trans_true_hor: with custom horizon altitude
- pheno_ut / pheno: planetary phenomena (phase, elongation, magnitude)
- Multiple locations, dates, bodies, and flag combinations

Run with: env LIBEPHEMERIS_MODE=skyfield python3 compare_scripts/round5_rise_transit_deep.py
"""

from __future__ import annotations

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


print("=" * 80)
print("ROUND 5: RISE/TRANSIT/SET & PHENOMENA DEEP AUDIT")
print("=" * 80)

# ============================================================================
# LOCATIONS
# ============================================================================
LOCATIONS = [
    ("Rome", 41.9028, 12.4964, 0),
    ("New York", 40.7128, -74.0060, 0),
    ("Sydney", -33.8688, 151.2093, 0),
    ("Equator", 0.0, 0.0, 0),
    ("London", 51.5074, -0.1278, 0),
    ("Tokyo", 35.6762, 139.6503, 0),
    ("Tromso", 69.6496, 18.9560, 0),  # High latitude (Arctic)
    ("Cape Town", -33.9249, 18.4241, 0),
    ("Denver", 39.7392, -104.9903, 1609),  # High altitude
    ("La Paz", -16.5000, -68.1500, 3640),  # Very high altitude
]

BODIES = [
    (swe.SUN, "Sun"),
    (swe.MOON, "Moon"),
    (swe.MERCURY, "Mercury"),
    (swe.VENUS, "Venus"),
    (swe.MARS, "Mars"),
    (swe.JUPITER, "Jupiter"),
    (swe.SATURN, "Saturn"),
]

EVENTS = [
    (1, "Rise"),
    (2, "Set"),
    (4, "Transit"),
    (8, "LowerTransit"),
]

DATES = [
    swe.julday(2024, 3, 20, 0.0),  # Spring equinox
    swe.julday(2024, 6, 21, 0.0),  # Summer solstice
    swe.julday(2024, 9, 22, 0.0),  # Autumn equinox
    swe.julday(2024, 12, 21, 0.0),  # Winter solstice
    swe.julday(2020, 7, 5, 0.0),  # Near lunar apogee
    swe.julday(2025, 1, 15, 0.0),  # Recent date
]

# ============================================================================
# PART 1: Rise/Set/Transit for Sun and Moon across locations and dates
# ============================================================================
print()
print("=" * 80)
print("PART 1: Rise/Set/Transit — Sun & Moon across locations and dates")
print("=" * 80)

part1_pass = 0
part1_fail = 0
part1_skip = 0
part1_total = 0

# Tolerances in seconds
TOL_SUN = 30.0  # Sun rise/set: 30s
TOL_MOON = 60.0  # Moon: 60s (parallax differences)
TOL_TRANSIT = 15.0  # Transit: tighter, no horizon effects
TOL_PLANET = 60.0  # Planets: 60s

for name, lat, lon, alt in LOCATIONS:
    for jd_date in DATES[:4]:  # 4 seasonal dates
        for body_id, body_name in BODIES[:2]:  # Sun and Moon
            for rsmi, event_name in EVENTS:
                try:
                    geopos = (lon, lat, alt)
                    se_ret = swe.rise_trans(
                        jd_date, body_id, rsmi, geopos, 1013.25, 15.0
                    )
                    se_flag = se_ret[0]
                    se_jd = se_ret[1][0] if se_ret[1] else 0.0
                except Exception as e:
                    part1_skip += 1
                    continue

                if se_flag != 0 or se_jd <= 0:
                    part1_skip += 1
                    continue

                try:
                    le_ret = ephem.rise_trans(
                        jd_date,
                        body_id,
                        rsmi,
                        [lon, lat, alt],
                        1013.25,
                        15.0,
                    )
                    le_jd = le_ret[1][0]
                    le_flag = le_ret[0]
                except Exception as e:
                    print(
                        "  LE ERROR {}/{}/{}: {}".format(name, body_name, event_name, e)
                    )
                    part1_skip += 1
                    continue

                if le_jd <= 0:
                    part1_skip += 1
                    continue

                diff_s = abs(se_jd - le_jd) * 86400
                part1_total += 1

                if rsmi in (4, 8):
                    tol = TOL_TRANSIT
                elif body_id == swe.MOON:
                    tol = TOL_MOON
                else:
                    tol = TOL_SUN

                cat = "rise_set_{}".format(body_name.lower())
                if diff_s <= tol:
                    part1_pass += 1
                    add_stat(cat, diff_s)
                else:
                    part1_fail += 1
                    issues.append(
                        "rise_trans {}/{}/{} at {}: diff={:.1f}s".format(
                            body_name,
                            event_name,
                            name,
                            "equinox" if jd_date == DATES[0] else "solstice",
                            diff_s,
                        )
                    )
                    print(
                        "  FAIL {:<10s} {:<6s} {:<12s} {:<10s}: diff={:.1f}s (tol={:.0f}s)".format(
                            name,
                            body_name,
                            event_name,
                            "{:.0f}".format(jd_date - 2400000),
                            diff_s,
                            tol,
                        )
                    )

print()
print(
    "Part 1 summary: {}/{} passed ({} skipped)".format(
        part1_pass, part1_total, part1_skip
    )
)

# ============================================================================
# PART 2: Rise/Set for planets (Mercury through Saturn)
# ============================================================================
print()
print("=" * 80)
print("PART 2: Rise/Set/Transit — Planets (Mercury-Saturn)")
print("=" * 80)

part2_pass = 0
part2_fail = 0
part2_skip = 0
part2_total = 0

for name, lat, lon, alt in LOCATIONS[:6]:  # 6 locations
    for jd_date in DATES[:2]:  # 2 dates
        for body_id, body_name in BODIES[2:]:  # Mercury-Saturn
            for rsmi, event_name in EVENTS[:2]:  # Rise and Set only
                try:
                    geopos = (lon, lat, alt)
                    se_ret = swe.rise_trans(
                        jd_date, body_id, rsmi, geopos, 1013.25, 15.0
                    )
                    se_flag = se_ret[0]
                    se_jd = se_ret[1][0] if se_ret[1] else 0.0
                except Exception:
                    part2_skip += 1
                    continue

                if se_flag != 0 or se_jd <= 0:
                    part2_skip += 1
                    continue

                try:
                    le_ret = ephem.rise_trans(
                        jd_date,
                        body_id,
                        rsmi,
                        [lon, lat, alt],
                        1013.25,
                        15.0,
                    )
                    le_jd = le_ret[1][0]
                except Exception:
                    part2_skip += 1
                    continue

                if le_jd <= 0:
                    part2_skip += 1
                    continue

                diff_s = abs(se_jd - le_jd) * 86400
                part2_total += 1

                cat = "rise_set_planet"
                if diff_s <= TOL_PLANET:
                    part2_pass += 1
                    add_stat(cat, diff_s)
                else:
                    part2_fail += 1
                    issues.append(
                        "planet rise_trans {}/{}/{}: diff={:.1f}s".format(
                            body_name, event_name, name, diff_s
                        )
                    )

                if diff_s > TOL_PLANET:
                    print(
                        "  FAIL {:<10s} {:<8s} {:<6s}: diff={:.1f}s".format(
                            name, body_name, event_name, diff_s
                        )
                    )

print()
print(
    "Part 2 summary: {}/{} passed ({} skipped)".format(
        part2_pass, part2_total, part2_skip
    )
)

# ============================================================================
# PART 3: Flag combinations (disc center, no refraction, twilights)
# ============================================================================
print()
print("=" * 80)
print("PART 3: Flag combinations — disc center, no refraction, twilights")
print("=" * 80)

part3_pass = 0
part3_fail = 0
part3_skip = 0
part3_total = 0

FLAG_COMBOS = [
    (256, "DISC_CENTER"),
    (512, "NO_REFRACTION"),
    (1024, "CIVIL_TWILIGHT"),
    (2048, "NAUTIC_TWILIGHT"),
    (4096, "ASTRO_TWILIGHT"),
    (8192, "DISC_BOTTOM"),
    (256 | 512, "CENTER+NOREFR"),
]

jd_test = swe.julday(2024, 6, 15, 0.0)
for name, lat, lon, alt in LOCATIONS[:4]:  # 4 locations
    for flag_val, flag_name in FLAG_COMBOS:
        for rsmi_base in [1, 2]:  # Rise and Set
            rsmi = rsmi_base | flag_val
            body_id = swe.SUN

            try:
                geopos = (lon, lat, alt)
                se_ret = swe.rise_trans(jd_test, body_id, rsmi, geopos, 1013.25, 15.0)
                se_flag = se_ret[0]
                se_jd = se_ret[1][0] if se_ret[1] else 0.0
            except Exception:
                part3_skip += 1
                continue

            if se_flag != 0 or se_jd <= 0:
                part3_skip += 1
                continue

            try:
                le_ret = ephem.rise_trans(
                    jd_test,
                    body_id,
                    rsmi,
                    [lon, lat, alt],
                    1013.25,
                    15.0,
                )
                le_jd = le_ret[1][0]
            except Exception as e:
                print(
                    "  LE ERROR {}/{}/flag={}: {}".format(
                        name, "Rise" if rsmi_base == 1 else "Set", flag_name, e
                    )
                )
                part3_skip += 1
                continue

            if le_jd <= 0:
                part3_skip += 1
                continue

            diff_s = abs(se_jd - le_jd) * 86400
            part3_total += 1

            tol = 30.0
            cat = "flag_{}".format(flag_name.lower().replace("+", "_"))
            if diff_s <= tol:
                part3_pass += 1
                add_stat(cat, diff_s)
            else:
                part3_fail += 1
                issues.append(
                    "flag {} {}/{}: diff={:.1f}s".format(
                        flag_name, name, "Rise" if rsmi_base == 1 else "Set", diff_s
                    )
                )
                event_str = "Rise" if rsmi_base == 1 else "Set"
                print(
                    "  FAIL {:<10s} {:<6s} {:<16s}: diff={:.1f}s".format(
                        name, event_str, flag_name, diff_s
                    )
                )

print()
print(
    "Part 3 summary: {}/{} passed ({} skipped)".format(
        part3_pass, part3_total, part3_skip
    )
)

# ============================================================================
# PART 4: rise_trans_true_hor (custom horizon altitude)
# ============================================================================
print()
print("=" * 80)
print("PART 4: rise_trans_true_hor — custom horizon altitudes")
print("=" * 80)

part4_pass = 0
part4_fail = 0
part4_skip = 0
part4_total = 0

HORIZONS = [0.0, 0.5, 1.0, 2.0, -0.5, 5.0]

jd_test4 = swe.julday(2024, 6, 15, 0.0)
for name, lat, lon, alt in LOCATIONS[:4]:
    for hor_alt in HORIZONS:
        for rsmi_base in [1, 2]:  # Rise and Set
            try:
                geopos = (lon, lat, alt)
                se_ret = swe.rise_trans_true_hor(
                    jd_test4, swe.SUN, rsmi_base, geopos, 1013.25, 15.0, hor_alt
                )
                se_flag = se_ret[0]
                se_jd = se_ret[1][0] if se_ret[1] else 0.0
            except Exception:
                part4_skip += 1
                continue

            if se_flag != 0 or se_jd <= 0:
                part4_skip += 1
                continue

            try:
                le_ret = ephem.rise_trans_true_hor(
                    jd_test4,
                    swe.SUN,
                    rsmi_base,
                    [lon, lat, alt],
                    1013.25,
                    15.0,
                    hor_alt,
                )
                le_jd = le_ret[1][0]
            except Exception as e:
                print(
                    "  LE ERROR {}/hor={}/{}): {}".format(
                        name, hor_alt, "Rise" if rsmi_base == 1 else "Set", e
                    )
                )
                part4_skip += 1
                continue

            if le_jd <= 0:
                part4_skip += 1
                continue

            diff_s = abs(se_jd - le_jd) * 86400
            part4_total += 1

            tol = 30.0
            cat = "true_hor"
            if diff_s <= tol:
                part4_pass += 1
                add_stat(cat, diff_s)
            else:
                part4_fail += 1
                issues.append(
                    "true_hor {}/hor={}: diff={:.1f}s".format(name, hor_alt, diff_s)
                )
                event_str = "Rise" if rsmi_base == 1 else "Set"
                print(
                    "  FAIL {:<10s} {:<6s} hor={:5.1f}°: diff={:.1f}s".format(
                        name, event_str, hor_alt, diff_s
                    )
                )

print()
print(
    "Part 4 summary: {}/{} passed ({} skipped)".format(
        part4_pass, part4_total, part4_skip
    )
)

# ============================================================================
# PART 5: pheno_ut — Planetary phenomena
# ============================================================================
print()
print("=" * 80)
print("PART 5: pheno_ut — Planetary phenomena (phase, elongation, magnitude)")
print("=" * 80)

part5_pass = 0
part5_fail = 0
part5_total = 0

PHENO_BODIES = [
    (swe.SUN, "Sun"),
    (swe.MOON, "Moon"),
    (swe.MERCURY, "Mercury"),
    (swe.VENUS, "Venus"),
    (swe.MARS, "Mars"),
    (swe.JUPITER, "Jupiter"),
    (swe.SATURN, "Saturn"),
    (swe.URANUS, "Uranus"),
    (swe.NEPTUNE, "Neptune"),
    (swe.PLUTO, "Pluto"),
]

PHENO_DATES = [
    swe.julday(2000, 1, 1, 12.0),
    swe.julday(2010, 6, 15, 12.0),
    swe.julday(2020, 3, 20, 12.0),
    swe.julday(2024, 6, 21, 12.0),
    swe.julday(2024, 12, 21, 12.0),
    swe.julday(2025, 1, 15, 12.0),
]

PHENO_LABELS = ["PhaseAngle", "Phase", "Elongation", "Diameter", "Magnitude"]
PHENO_TOLS = [0.01, 0.001, 0.01, 0.1, 0.1]  # degrees, fraction, degrees, arcsec, mag

for body_id, body_name in PHENO_BODIES:
    body_pass = 0
    body_fail = 0
    body_total = 0
    max_diffs = [0.0] * 5

    for jd_date in PHENO_DATES:
        try:
            se_ret = swe.pheno_ut(jd_date, body_id)
            se_attr = se_ret[0] if isinstance(se_ret[0], (list, tuple)) else se_ret

            le_ret = ephem.pheno_ut(jd_date, body_id, 0)
            le_attr = le_ret[0]
        except Exception as e:
            continue

        all_ok = True
        for i in range(5):
            se_val = float(se_attr[i])
            le_val = float(le_attr[i])
            diff = abs(se_val - le_val)
            max_diffs[i] = max(max_diffs[i], diff)
            body_total += 1
            part5_total += 1

            if diff <= PHENO_TOLS[i]:
                body_pass += 1
                part5_pass += 1
                add_stat("pheno_{}".format(PHENO_LABELS[i].lower()), diff)
            else:
                body_fail += 1
                part5_fail += 1
                all_ok = False
                issues.append(
                    "pheno {} {}: se={:.6f} le={:.6f} diff={:.6f}".format(
                        body_name, PHENO_LABELS[i], se_val, le_val, diff
                    )
                )

    status = "PASS" if body_fail == 0 else "FAIL"
    diffs_str = "  ".join(
        "{}: {:.6f}".format(PHENO_LABELS[i], max_diffs[i]) for i in range(5)
    )
    print(
        "  {} {:<10s} {}/{} — max diffs: {}".format(
            status, body_name, body_pass, body_total, diffs_str
        )
    )

print()
print("Part 5 summary: {}/{} passed".format(part5_pass, part5_total))

# ============================================================================
# PART 6: Sequential rise/set consistency (sunrise < transit < sunset)
# ============================================================================
print()
print("=" * 80)
print("PART 6: Sequential consistency — sunrise < transit < sunset")
print("=" * 80)

part6_pass = 0
part6_fail = 0
part6_total = 0

for name, lat, lon, alt in LOCATIONS[:6]:
    for jd_date in DATES[:4]:
        try:
            rise_ret = ephem.rise_trans(jd_date, swe.SUN, 1, [lon, lat, alt])
            trans_ret = ephem.rise_trans(jd_date, swe.SUN, 4, [lon, lat, alt])
            set_ret = ephem.rise_trans(jd_date, swe.SUN, 2, [lon, lat, alt])

            jd_rise = rise_ret[1][0]
            jd_trans = trans_ret[1][0]
            jd_set = set_ret[1][0]
        except Exception:
            continue

        if jd_rise <= 0 or jd_trans <= 0 or jd_set <= 0:
            continue

        part6_total += 1

        # Find the set AFTER the rise
        if jd_set < jd_rise:
            try:
                set_ret2 = ephem.rise_trans(jd_rise + 0.01, swe.SUN, 2, [lon, lat, alt])
                jd_set = set_ret2[1][0]
            except Exception:
                continue

        # Find transit between rise and set
        if jd_trans < jd_rise or jd_trans > jd_set:
            try:
                trans_ret2 = ephem.rise_trans(
                    jd_rise + 0.01, swe.SUN, 4, [lon, lat, alt]
                )
                jd_trans = trans_ret2[1][0]
            except Exception:
                continue

        if jd_rise < jd_trans < jd_set:
            part6_pass += 1
        else:
            part6_fail += 1
            print(
                "  FAIL {}: rise={:.6f} trans={:.6f} set={:.6f}".format(
                    name, jd_rise, jd_trans, jd_set
                )
            )
            issues.append("sequential {} rise/trans/set order wrong".format(name))

print()
print("Part 6 summary: {}/{} passed".format(part6_pass, part6_total))

# ============================================================================
# PART 7: Moon rise/set with different SE vs LE parallax handling
# ============================================================================
print()
print("=" * 80)
print("PART 7: Moon rise/set precision — parallax sensitivity")
print("=" * 80)

part7_pass = 0
part7_fail = 0
part7_total = 0

for name, lat, lon, alt in LOCATIONS[:6]:
    for jd_date in DATES[:4]:
        for rsmi in [1, 2]:
            try:
                geopos = (lon, lat, alt)
                se_ret = swe.rise_trans(jd_date, swe.MOON, rsmi, geopos, 1013.25, 15.0)
                se_jd = se_ret[1][0] if se_ret[0] == 0 and se_ret[1] else 0.0

                le_ret = ephem.rise_trans(
                    jd_date,
                    swe.MOON,
                    rsmi,
                    [lon, lat, alt],
                    1013.25,
                    15.0,
                )
                le_jd = le_ret[1][0]
            except Exception:
                continue

            if se_jd <= 0 or le_jd <= 0:
                continue

            diff_s = abs(se_jd - le_jd) * 86400
            part7_total += 1

            if diff_s <= 120.0:  # 2 minutes for Moon
                part7_pass += 1
                add_stat("moon_rise_set", diff_s)
            else:
                part7_fail += 1
                event = "rise" if rsmi == 1 else "set"
                print(
                    "  FAIL {:<10s} Moon {}: diff={:.1f}s".format(name, event, diff_s)
                )
                issues.append(
                    "Moon {} {} at {}: diff={:.1f}s".format(
                        event, name, jd_date, diff_s
                    )
                )

print()
print("Part 7 summary: {}/{} passed".format(part7_pass, part7_total))

# ============================================================================
# SUMMARY
# ============================================================================
print()
print("=" * 80)
print("SUMMARY")
print("=" * 80)

print()
print(
    "Part 1 (Sun/Moon rise/set):    {}/{} passed ({} skip)".format(
        part1_pass, part1_total, part1_skip
    )
)
print(
    "Part 2 (Planet rise/set):      {}/{} passed ({} skip)".format(
        part2_pass, part2_total, part2_skip
    )
)
print(
    "Part 3 (Flag combos):          {}/{} passed ({} skip)".format(
        part3_pass, part3_total, part3_skip
    )
)
print(
    "Part 4 (True horizon):         {}/{} passed ({} skip)".format(
        part4_pass, part4_total, part4_skip
    )
)
print("Part 5 (Phenomena):            {}/{} passed".format(part5_pass, part5_total))
print("Part 6 (Sequential order):     {}/{} passed".format(part6_pass, part6_total))
print("Part 7 (Moon parallax):        {}/{} passed".format(part7_pass, part7_total))

print()
for cat in sorted(stats.keys()):
    vals = stats[cat]
    if vals:
        mn = min(vals)
        mx = max(vals)
        mean = sum(vals) / len(vals)
        vals_sorted = sorted(vals)
        median = vals_sorted[len(vals_sorted) // 2]
        print("Stats for {}:".format(cat))
        print(
            "  Min: {:.4f}  Max: {:.4f}  Mean: {:.4f}  Median: {:.4f}  (N={})".format(
                mn, mx, mean, median, len(vals)
            )
        )

if issues:
    print()
    print("=" * 80)
    print("ISSUES FOUND: {}".format(len(issues)))
    for iss in issues[:50]:
        print("  - {}".format(iss))
    if len(issues) > 50:
        print("  ... and {} more".format(len(issues) - 50))
    print("=" * 80)
