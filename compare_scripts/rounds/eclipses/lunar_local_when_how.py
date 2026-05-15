"""Round 3: Lunar Eclipse Local Functions Deep Audit

Comprehensive comparison of libephemeris vs pyswisseph for:
- lun_eclipse_when: global lunar eclipse search
- lun_eclipse_when_loc: local lunar eclipse search with visibility
- lun_eclipse_how / lun_eclipse_how: lunar eclipse circumstances at a given time
- Attribute comparison: umbral magnitude, penumbral magnitude, Moon position
- Contact time comparison: all 7 contact times
- Sequential eclipse search: consecutive eclipses at multiple locations

Run with: env LIBEPHEMERIS_MODE=skyfield python3 compare_scripts/round3_lunar_eclipse_local_deep.py
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

# Test locations: (name, lon, lat, alt)
LOCATIONS = [
    ("Rome, Italy", 12.4964, 41.9028, 0.0),
    ("New York, NY", -74.006, 40.7128, 0.0),
    ("Tokyo, Japan", 139.6917, 35.6895, 0.0),
    ("Sydney, Australia", 151.2093, -33.8688, 0.0),
    ("London, UK", -0.1278, 51.5074, 0.0),
    ("Cairo, Egypt", 31.2357, 30.0444, 0.0),
    ("Sao Paulo, Brazil", -46.6333, -23.5505, 0.0),
    ("Mumbai, India", 72.8777, 19.0760, 0.0),
    ("Reykjavik, Iceland", -21.9426, 64.1466, 0.0),
    ("Santiago, Chile", -70.6693, -33.4489, 0.0),
]

# Eclipse epochs to test
ECLIPSE_YEARS = [2018, 2019, 2021, 2022, 2023, 2024, 2025]

issues = []
stats_timing = []
stats_mag = []
stats_penmag = []
stats_moonalt = []
stats_moonaz = []
stats_contact = []


def jd_to_utstr(jd):
    """Convert JD to readable UT string."""
    if jd <= 0:
        return "0"
    z = int(jd + 0.5)
    f = jd + 0.5 - z
    if z < 2299161:
        a = z
    else:
        alpha = int((z - 1867216.25) / 36524.25)
        a = z + 1 + alpha - int(alpha / 4)
    b = a + 1524
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)
    day = b - d - int(30.6001 * e) + f
    month = e - 1 if e < 14 else e - 13
    year = c - 4716 if month > 2 else c - 4715
    hour = (day % 1) * 24
    minute = (hour % 1) * 60
    second = (minute % 1) * 60
    return "{:04d}-{:02d}-{:02d} {:02d}:{:02d}:{:02d} UT".format(
        int(year), int(month), int(day), int(hour), int(minute), int(second)
    )


print("=" * 80)
print("ROUND 3: LUNAR ECLIPSE LOCAL FUNCTIONS DEEP AUDIT")
print("=" * 80)

# ============================================================================
# PART 1: lun_eclipse_when — global eclipse search comparison
# ============================================================================
print()
print("=" * 80)
print("PART 1: lun_eclipse_when — global eclipse search")
print("=" * 80)

part1_pass = 0
part1_fail = 0
part1_total = 0

for year in ECLIPSE_YEARS:
    jd_start = swe.julday(year, 1, 1, 0.0)

    # Find all eclipses in this year with SE
    se_eclipses = []
    jd = jd_start
    for _ in range(6):  # Max 6 eclipses per year
        try:
            ret = swe.lun_eclipse_when(jd)
            se_flags = ret[0]
            se_tret = ret[1]
            if se_tret[0] <= 0 or se_tret[0] > jd_start + 366:
                break
            se_eclipses.append((se_flags, se_tret))
            jd = se_tret[0] + 25
        except Exception:
            break

    # Find all eclipses in this year with LE
    le_eclipses = []
    jd = jd_start
    for _ in range(6):
        try:
            ret = ephem.lun_eclipse_when(jd, 0)
            le_flags = ret[0]
            le_tret = ret[1]
            if le_tret[0] <= 0 or le_tret[0] > jd_start + 366:
                break
            le_eclipses.append((le_flags, le_tret))
            jd = le_tret[0] + 25
        except Exception:
            break

    # Compare
    n_se = len(se_eclipses)
    n_le = len(le_eclipses)

    if n_se != n_le:
        print("  FAIL {}: SE found {} eclipses, LE found {}".format(year, n_se, n_le))
        issues.append("{}: eclipse count mismatch SE={} LE={}".format(year, n_se, n_le))
        part1_fail += 1
        part1_total += 1
        continue

    year_ok = True
    for i in range(n_se):
        se_f, se_t = se_eclipses[i]
        le_f, le_t = le_eclipses[i]

        max_diff = abs(se_t[0] - le_t[0]) * 86400

        # Compare all contact times
        contact_names = [
            "Maximum",
            "Partial begin",
            "Partial end",
            "Total begin",
            "Total end",
            "Penumbral begin",
            "Penumbral end",
        ]
        contact_diffs = []
        for ci in range(7):
            if ci < len(se_t) and ci < len(le_t):
                if se_t[ci] > 0 and le_t[ci] > 0:
                    cd = abs(se_t[ci] - le_t[ci]) * 86400
                    contact_diffs.append(cd)
                elif se_t[ci] > 0:
                    contact_diffs.append(-1)  # SE has, LE missing
                elif le_t[ci] > 0:
                    contact_diffs.append(-2)  # LE has, SE missing

        max_contact = max(contact_diffs) if contact_diffs else 0
        missing = [d for d in contact_diffs if d < 0]

        part1_total += 1
        # Tolerance: 120s for timing
        if max_diff > 120 or max_contact > 120 or missing:
            status = "FAIL"
            part1_fail += 1
            year_ok = False
            reason = ""
            if max_diff > 120:
                reason += " max_diff={:.1f}s".format(max_diff)
            if max_contact > 120:
                reason += " contact_diff={:.1f}s".format(max_contact)
            if missing:
                reason += " missing_contacts={}".format(len(missing))
            issues.append(
                "{} eclipse #{} @ {}: {}".format(
                    year, i + 1, jd_to_utstr(se_t[0]), reason.strip()
                )
            )
        else:
            status = "PASS"
            part1_pass += 1
            stats_timing.append(max_diff)
            for cd in contact_diffs:
                if cd >= 0:
                    stats_contact.append(cd)

        print(
            "  {} {} eclipse #{}: max_diff={:.1f}s, max_contact={:.1f}s @ {}".format(
                status,
                year,
                i + 1,
                max_diff,
                max_contact,
                jd_to_utstr(se_t[0]),
            )
        )

print()
print("Part 1 summary: {}/{} passed".format(part1_pass, part1_total))

# ============================================================================
# PART 2: lun_eclipse_how — attribute comparison at known eclipse times
# ============================================================================
print()
print("=" * 80)
print("PART 2: lun_eclipse_how — attribute comparison")
print("=" * 80)

part2_pass = 0
part2_fail = 0
part2_total = 0

# Known lunar eclipses with approximate maximum times
known_eclipses = [
    ("2024 Mar 25 Penumbral", swe.julday(2024, 3, 25, 7.0)),
    ("2024 Sep 18 Partial", swe.julday(2024, 9, 18, 2.5)),
    ("2025 Mar 14 Total", swe.julday(2025, 3, 14, 6.6)),
    ("2025 Sep 7 Total", swe.julday(2025, 9, 7, 18.3)),
    ("2022 May 16 Total", swe.julday(2022, 5, 16, 4.1)),
    ("2022 Nov 8 Total", swe.julday(2022, 11, 8, 11.0)),
    ("2021 May 26 Total", swe.julday(2021, 5, 26, 11.2)),
    ("2021 Nov 19 Partial", swe.julday(2021, 11, 19, 9.0)),
    ("2019 Jan 21 Total", swe.julday(2019, 1, 21, 5.1)),
    ("2019 Jul 16 Partial", swe.julday(2019, 7, 16, 21.3)),
    ("2018 Jan 31 Total", swe.julday(2018, 1, 31, 13.3)),
    ("2018 Jul 27 Total", swe.julday(2018, 7, 27, 20.2)),
    ("2023 May 5 Penumbral", swe.julday(2023, 5, 5, 17.2)),
    ("2023 Oct 28 Partial", swe.julday(2023, 10, 28, 20.1)),
]

for ecl_name, jd_approx in known_eclipses:
    print()
    print("  Eclipse: {}".format(ecl_name))

    # Find exact maximum from SE global search
    se_ret = swe.lun_eclipse_when(jd_approx - 5)
    se_max = se_ret[1][0]

    # Compare how at multiple locations
    for loc_name, lon, lat, alt in LOCATIONS[:6]:
        geopos = (lon, lat, alt)

        try:
            se_how = swe.lun_eclipse_how(se_max, geopos=geopos)
            se_attr = se_how[1]
        except Exception as e:
            print("    SE error at {}: {}".format(loc_name, e))
            continue

        try:
            le_how = ephem.lun_eclipse_how(se_max, 0, geopos=geopos)
            le_attr = le_how[1]
        except Exception as e:
            print("    LE error at {}: {}".format(loc_name, e))
            continue

        # Compare attributes
        mag_diff = abs(se_attr[0] - float(le_attr[0]))
        penmag_diff = abs(se_attr[1] - float(le_attr[1]))
        az_diff = abs(se_attr[4] - float(le_attr[4]))
        if az_diff > 180:
            az_diff = 360 - az_diff
        alt_diff = abs(se_attr[5] - float(le_attr[5]))
        appalt_diff = abs(se_attr[6] - float(le_attr[6]))
        dist_diff = abs(se_attr[7] - float(le_attr[7]))

        part2_total += 1

        # Tolerances
        tol_mag = 0.01
        tol_penmag = 0.01
        tol_az = 0.1
        tol_alt = 0.1

        failed = False
        fail_reasons = []
        if mag_diff > tol_mag:
            failed = True
            fail_reasons.append("mag={:.4f}".format(mag_diff))
        if penmag_diff > tol_penmag:
            failed = True
            fail_reasons.append("penmag={:.4f}".format(penmag_diff))
        if abs(se_attr[5]) < 80 and az_diff > tol_az:  # Skip azimuth near zenith/nadir
            failed = True
            fail_reasons.append("az={:.4f}".format(az_diff))
        if alt_diff > tol_alt:
            failed = True
            fail_reasons.append("alt={:.4f}".format(alt_diff))

        if failed:
            status = "FAIL"
            part2_fail += 1
            issues.append(
                "{} @ {}: how attrs {}".format(
                    ecl_name, loc_name, ", ".join(fail_reasons)
                )
            )
        else:
            status = "PASS"
            part2_pass += 1
            stats_mag.append(mag_diff)
            stats_penmag.append(penmag_diff)
            stats_moonalt.append(alt_diff)
            if abs(se_attr[5]) < 80:
                stats_moonaz.append(az_diff)

        print(
            "    {} {}: umag={:.4f}(d={:.4f}), pmag={:.4f}(d={:.4f}), "
            "alt={:.1f}°(d={:.2f}°), az_d={:.2f}°".format(
                status,
                loc_name,
                se_attr[0],
                mag_diff,
                se_attr[1],
                penmag_diff,
                se_attr[5],
                alt_diff,
                az_diff,
            )
        )

print()
print("Part 2 summary: {}/{} passed".format(part2_pass, part2_total))

# ============================================================================
# PART 3: lun_eclipse_when_loc — local eclipse search comparison
# ============================================================================
print()
print("=" * 80)
print("PART 3: lun_eclipse_when_loc — local visibility comparison")
print("=" * 80)

part3_pass = 0
part3_fail = 0
part3_skip = 0
part3_total = 0

for year in [2018, 2021, 2024, 2025]:
    jd_start = swe.julday(year, 1, 1, 0.0)

    for loc_name, lon, lat, alt in LOCATIONS:
        geopos_se = (lon, lat, alt)

        try:
            se_ret = swe.lun_eclipse_when_loc(jd_start, geopos=geopos_se)
            se_flags = se_ret[0]
            se_tret = se_ret[1]
            se_attr = se_ret[2]
        except Exception as e:
            print("  SKIP {} @ {}: SE error: {}".format(year, loc_name, e))
            part3_skip += 1
            continue

        if se_tret[0] <= 0 or se_tret[0] > jd_start + 366:
            # No eclipse found this year
            continue

        try:
            le_ret = ephem.lun_eclipse_when_loc(jd_start, (lon, lat, alt), 0)
            le_flags = le_ret[0]
            le_tret = le_ret[1]
            le_attr = le_ret[2]
        except Exception as e:
            print("  SKIP {} @ {}: LE error: {}".format(year, loc_name, e))
            part3_skip += 1
            continue

        if le_tret[0] <= 0 or le_tret[0] > jd_start + 366:
            # LE didn't find eclipse that SE found
            part3_total += 1
            part3_fail += 1
            max_diff = 0
            mag_diff = 0
            print(
                "  FAIL {} @ {}: LE did not find eclipse SE found at {}".format(
                    year, loc_name, jd_to_utstr(se_tret[0])
                )
            )
            issues.append(
                "when_loc {} @ {}: LE missed eclipse at {}".format(
                    year, loc_name, jd_to_utstr(se_tret[0])
                )
            )
            continue

        # Check if they found the same eclipse (within 2 days)
        if abs(se_tret[0] - le_tret[0]) > 2:
            # Different eclipses found
            part3_total += 1
            part3_fail += 1
            print(
                "  FAIL {} @ {}: different eclipses SE={} LE={}".format(
                    year,
                    loc_name,
                    jd_to_utstr(se_tret[0]),
                    jd_to_utstr(le_tret[0]),
                )
            )
            issues.append("when_loc {} @ {}: different eclipses".format(year, loc_name))
            continue

        part3_total += 1

        # Compare timing
        max_diff = abs(se_tret[0] - le_tret[0]) * 86400

        # Compare attributes
        # attr[0] = umbral magnitude
        # attr[1] = penumbral magnitude
        mag_diff = abs(float(se_attr[0]) - float(le_attr[0]))
        penmag_diff = abs(float(se_attr[1]) - float(le_attr[1]))

        # Contact time differences
        contact_names = [
            "Maximum",
            "Partial begin",
            "Partial end",
            "Total begin",
            "Total end",
            "Penumbral begin",
            "Penumbral end",
            "Moonrise",
            "Moonset",
        ]
        contact_issues_list = []
        contact_diffs_local = []
        for ci in range(min(9, len(se_tret), len(le_tret))):
            se_t = se_tret[ci]
            le_t = le_tret[ci]
            if se_t > 0 and le_t > 0:
                cd = abs(se_t - le_t) * 86400
                contact_diffs_local.append(cd)
            elif se_t > 0 and le_t <= 0:
                cn = (
                    contact_names[ci] if ci < len(contact_names) else "t[{}]".format(ci)
                )
                contact_issues_list.append(
                    "WARN {} — SE={}, LE=0".format(cn, jd_to_utstr(se_t))
                )
            elif le_t > 0 and se_t <= 0:
                cn = (
                    contact_names[ci] if ci < len(contact_names) else "t[{}]".format(ci)
                )
                contact_issues_list.append(
                    "WARN {} — SE=0, LE={}".format(cn, jd_to_utstr(le_t))
                )

        max_contact_diff = max(contact_diffs_local) if contact_diffs_local else 0

        # Tolerances
        tol_time = 120  # seconds
        tol_mag = 0.02
        tol_penmag = 0.02

        failed = False
        fail_reasons = []

        if max_diff > tol_time:
            failed = True
            fail_reasons.append("max_diff={:.1f}s".format(max_diff))
        if mag_diff > tol_mag:
            failed = True
            fail_reasons.append("mag_diff={:.4f}".format(mag_diff))
        if penmag_diff > tol_penmag:
            failed = True
            fail_reasons.append("penmag_diff={:.4f}".format(penmag_diff))

        if failed:
            status = "FAIL"
            part3_fail += 1
            issues.append(
                "when_loc {} @ {}: {}".format(year, loc_name, ", ".join(fail_reasons))
            )
        else:
            status = "PASS"
            part3_pass += 1
            stats_timing.append(max_diff)
            stats_mag.append(mag_diff)
            stats_penmag.append(penmag_diff)
            for cd in contact_diffs_local:
                stats_contact.append(cd)

        print(
            "  {} {} @ {}: max_diff={:.1f}s, mag_diff={:.4f}, "
            "penmag_diff={:.4f}".format(
                status, year, loc_name, max_diff, mag_diff, penmag_diff
            )
        )
        for ci_msg in contact_issues_list:
            print("    {}".format(ci_msg))

print()
print(
    "Part 3 summary: {}/{} passed ({} skipped)".format(
        part3_pass, part3_total, part3_skip
    )
)

# ============================================================================
# PART 4: Detailed contact time comparison — specific eclipses
# ============================================================================
print()
print("=" * 80)
print("PART 4: Detailed contact time comparison")
print("=" * 80)

detail_eclipses = [
    ("2025 Mar 14 Total", 2025, 3, 14, 6.6),
    ("2024 Sep 18 Partial", 2024, 9, 18, 2.5),
    ("2022 Nov 8 Total", 2022, 11, 8, 11.0),
    ("2019 Jan 21 Total", 2019, 1, 21, 5.1),
    ("2018 Jul 27 Total", 2018, 7, 27, 20.2),
]

detail_locs = [
    ("Rome, Italy", 12.4964, 41.9028, 0.0),
    ("New York, NY", -74.006, 40.7128, 0.0),
    ("Tokyo, Japan", 139.6917, 35.6895, 0.0),
    ("Sydney, Australia", 151.2093, -33.8688, 0.0),
    ("Cairo, Egypt", 31.2357, 30.0444, 0.0),
    ("Santiago, Chile", -70.6693, -33.4489, 0.0),
]

for ecl_name, year, month, day, hour in detail_eclipses:
    print()
    print("  {} — checking at 6 locations".format(ecl_name))

    # Global max from SE
    jd_approx = swe.julday(year, month, day, hour)
    se_global = swe.lun_eclipse_when(jd_approx - 5)
    se_max_global = se_global[1][0]

    for loc_name, lon, lat, alt in detail_locs:
        geopos_se = (lon, lat, alt)

        # SE when_loc
        try:
            se_ret = swe.lun_eclipse_when_loc(
                swe.julday(year, month, day - 3, 0.0), geopos=geopos_se
            )
            se_tret = se_ret[1]
            se_attr = se_ret[2]
            se_flags = se_ret[0]
        except Exception as e:
            print("    {} — SE error: {}".format(loc_name, e))
            continue

        # Check eclipse is in the right date range
        if se_tret[0] <= 0 or abs(se_tret[0] - se_max_global) > 2:
            print("    {} — not visible per SE".format(loc_name))
            continue

        # LE when_loc
        try:
            le_ret = ephem.lun_eclipse_when_loc(
                swe.julday(year, month, day - 3, 0.0), (lon, lat, alt), 0
            )
            le_tret = le_ret[1]
            le_attr = le_ret[2]
            le_flags = le_ret[0]
        except Exception as e:
            print("    {} — LE error: {}".format(loc_name, e))
            continue

        if le_tret[0] <= 0 or abs(le_tret[0] - se_max_global) > 2:
            print(
                "    {} — LE missed this eclipse (found {} instead)".format(
                    loc_name, jd_to_utstr(le_tret[0]) if le_tret[0] > 0 else "nothing"
                )
            )
            continue

        print("    {}:".format(loc_name))
        print("      Eclipse type: SE=0x{:04x}, LE=0x{:04x}".format(se_flags, le_flags))

        # Contact times
        contact_names = [
            "Maximum       ",
            "Partial begin ",
            "Partial end   ",
            "Total begin   ",
            "Total end     ",
            "Penumbral beg ",
            "Penumbral end ",
            "Moonrise      ",
            "Moonset       ",
        ]
        for ci in range(min(9, len(se_tret), len(le_tret))):
            se_t = se_tret[ci]
            le_t = le_tret[ci]
            cn = (
                contact_names[ci]
                if ci < len(contact_names)
                else "t[{}]           ".format(ci)
            )
            if se_t > 0 and le_t > 0:
                diff = abs(se_t - le_t) * 86400
                marker = "[OK]" if diff < 120 else "[!!]"
                print(
                    "      {}: diff={:8.1f}s  SE={}  LE={}  {}".format(
                        cn, diff, jd_to_utstr(se_t), jd_to_utstr(le_t), marker
                    )
                )
            elif se_t > 0:
                print("      {}: SE={}, LE=0 (MISSING)".format(cn, jd_to_utstr(se_t)))
            elif le_t > 0:
                print("      {}: SE=0, LE={} (EXTRA)".format(cn, jd_to_utstr(le_t)))

        # Attributes
        print(
            "      {:20s} {:>12s} {:>12s} {:>10s}".format(
                "Attribute", "SE", "LE", "Diff"
            )
        )
        attr_names = [
            "Umbral mag",
            "Penumbral mag",
            "(unused)",
            "(unused)",
            "Moon az",
            "Moon alt(true)",
            "Moon alt(app)",
            "Distance from opp",
        ]
        for ai in range(min(8, len(se_attr), len(le_attr))):
            an = attr_names[ai]
            sv = se_attr[ai]
            lv = float(le_attr[ai])
            diff = abs(sv - lv)
            print("      {:20s} {:12.4f} {:12.4f} {:10.4f}".format(an, sv, lv, diff))


# ============================================================================
# PART 5: Sequential local eclipse search — consecutive eclipses
# ============================================================================
print()
print("=" * 80)
print("PART 5: Sequential eclipse search — 5 consecutive eclipses per location")
print("=" * 80)

seq_locs = [
    ("Rome, Italy", 12.4964, 41.9028, 0.0),
    ("Tokyo, Japan", 139.6917, 35.6895, 0.0),
    ("New York, NY", -74.006, 40.7128, 0.0),
    ("Santiago, Chile", -70.6693, -33.4489, 0.0),
]

for loc_name, lon, lat, alt in seq_locs:
    geopos_se = (lon, lat, alt)
    print()
    print("  {}: searching 5 consecutive local lunar eclipses...".format(loc_name))

    jd_se = swe.julday(2020, 1, 1, 0.0)
    jd_le = jd_se

    for ecl_i in range(5):
        # SE search
        try:
            se_ret = swe.lun_eclipse_when_loc(jd_se, geopos=geopos_se)
            se_tret = se_ret[1]
            se_attr = se_ret[2]
            se_flags = se_ret[0]
        except Exception as e:
            print("    Eclipse #{}: SE error: {}".format(ecl_i + 1, e))
            break

        if se_tret[0] <= 0:
            print("    Eclipse #{}: SE found no more eclipses".format(ecl_i + 1))
            break

        # LE search
        try:
            le_ret = ephem.lun_eclipse_when_loc(jd_le, (lon, lat, alt), 0)
            le_tret = le_ret[1]
            le_attr = le_ret[2]
            le_flags = le_ret[0]
        except Exception as e:
            print("    Eclipse #{}: LE error: {}".format(ecl_i + 1, e))
            break

        if le_tret[0] <= 0:
            print(
                "    Eclipse #{}: LE found nothing, SE found {}".format(
                    ecl_i + 1, jd_to_utstr(se_tret[0])
                )
            )
            break

        # Check if same eclipse
        time_diff = abs(se_tret[0] - le_tret[0]) * 86400
        mag_diff = abs(float(se_attr[0]) - float(le_attr[0]))

        if abs(se_tret[0] - le_tret[0]) > 2:
            # Different eclipses — report and try to sync
            print(
                "    Eclipse #{}: DIFFERENT — SE={}, LE={}".format(
                    ecl_i + 1,
                    jd_to_utstr(se_tret[0]),
                    jd_to_utstr(le_tret[0]),
                )
            )
            issues.append(
                "seq {} #{}: different eclipses SE={} LE={}".format(
                    loc_name,
                    ecl_i + 1,
                    jd_to_utstr(se_tret[0]),
                    jd_to_utstr(le_tret[0]),
                )
            )
            # Sync forward from whichever is later
            jd_se = max(se_tret[0], le_tret[0]) + 25
            jd_le = jd_se
            continue

        status = "PASS" if time_diff < 120 and mag_diff < 0.02 else "FAIL"
        if status == "FAIL":
            issues.append(
                "seq {} #{} @ {}: time_diff={:.1f}s, mag_diff={:.4f}".format(
                    loc_name, ecl_i + 1, jd_to_utstr(se_tret[0]), time_diff, mag_diff
                )
            )

        print(
            "    Eclipse #{}: {}  time_diff={:6.1f}s  mag_diff={:.4f}  [{}]".format(
                ecl_i + 1,
                jd_to_utstr(se_tret[0]),
                time_diff,
                mag_diff,
                status,
            )
        )

        jd_se = se_tret[0] + 25
        jd_le = le_tret[0] + 25


# ============================================================================
# SUMMARY
# ============================================================================
print()
print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()
print("Part 1 (when — global): {}/{} passed".format(part1_pass, part1_total))
print("Part 2 (how — attrs):   {}/{} passed".format(part2_pass, part2_total))
print(
    "Part 3 (when_loc):      {}/{} passed ({} skipped)".format(
        part3_pass, part3_total, part3_skip
    )
)

if stats_timing:
    print()
    print("Timing diffs (maximum time):")
    print(
        "  Min: {:.1f}s  Max: {:.1f}s  Mean: {:.1f}s  Median: {:.1f}s".format(
            min(stats_timing),
            max(stats_timing),
            sum(stats_timing) / len(stats_timing),
            sorted(stats_timing)[len(stats_timing) // 2],
        )
    )

if stats_contact:
    print()
    print("Contact time diffs:")
    print(
        "  Min: {:.1f}s  Max: {:.1f}s  Mean: {:.1f}s  Median: {:.1f}s".format(
            min(stats_contact),
            max(stats_contact),
            sum(stats_contact) / len(stats_contact),
            sorted(stats_contact)[len(stats_contact) // 2],
        )
    )

if stats_mag:
    print()
    print("Umbral magnitude diffs:")
    print(
        "  Min: {:.4f}  Max: {:.4f}  Mean: {:.4f}  Median: {:.4f}".format(
            min(stats_mag),
            max(stats_mag),
            sum(stats_mag) / len(stats_mag),
            sorted(stats_mag)[len(stats_mag) // 2],
        )
    )

if stats_penmag:
    print()
    print("Penumbral magnitude diffs:")
    print(
        "  Min: {:.4f}  Max: {:.4f}  Mean: {:.4f}  Median: {:.4f}".format(
            min(stats_penmag),
            max(stats_penmag),
            sum(stats_penmag) / len(stats_penmag),
            sorted(stats_penmag)[len(stats_penmag) // 2],
        )
    )

if stats_moonalt:
    print()
    print("Moon altitude diffs:")
    print(
        "  Min: {:.2f}°  Max: {:.2f}°  Mean: {:.2f}°  Median: {:.2f}°".format(
            min(stats_moonalt),
            max(stats_moonalt),
            sum(stats_moonalt) / len(stats_moonalt),
            sorted(stats_moonalt)[len(stats_moonalt) // 2],
        )
    )

if stats_moonaz:
    print()
    print("Moon azimuth diffs:")
    print(
        "  Min: {:.2f}°  Max: {:.2f}°  Mean: {:.2f}°  Median: {:.2f}°".format(
            min(stats_moonaz),
            max(stats_moonaz),
            sum(stats_moonaz) / len(stats_moonaz),
            sorted(stats_moonaz)[len(stats_moonaz) // 2],
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
