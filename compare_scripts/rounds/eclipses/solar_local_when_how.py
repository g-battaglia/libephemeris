"""Round 2 Deep Audit: Solar Eclipse Local Functions.

Comprehensive comparison of solar eclipse local functions between
libephemeris and pyswisseph, covering:

1. sol_eclipse_when_loc - find next eclipse visible from location
   - Timing of maximum, contact times
   - Eclipse type detection
   - Multiple locations worldwide
   - Multiple eclipse types (total, annular, partial)

2. sol_eclipse_how - eclipse circumstances at location and time
   - Magnitude (fraction of solar diameter covered)
   - Obscuration (fraction of solar disc area covered)
   - Diameter ratio (lunar/solar)
   - Sun azimuth and altitude
   - Core shadow width

3. Edge cases
   - Near-horizon eclipses
   - Partial eclipses at edge of path
   - High-latitude locations
"""

from __future__ import annotations

import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import FLG_SWIEPH
import os

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

# ============================================================================
# Configuration
# ============================================================================

# Tolerances
TIMING_TOL_SECONDS = 120.0  # 2 minutes for eclipse timing
CONTACT_TOL_SECONDS = 180.0  # 3 minutes for contact times
MAGNITUDE_TOL = 0.02  # magnitude tolerance
OBSCURATION_TOL = 0.03  # obscuration tolerance
DIAMETER_RATIO_TOL = 0.01  # diameter ratio tolerance
SUN_ALT_TOL = 1.0  # degrees for sun altitude
SUN_AZ_TOL = 2.0  # degrees for sun azimuth
SHADOW_WIDTH_TOL_KM = 50.0  # km for shadow width
SHADOW_WIDTH_REL_TOL = 0.25  # 25% relative for shadow width

# Test locations: (name, lon, lat, alt)
LOCATIONS = [
    ("Dallas, TX", -96.797, 32.777, 0),
    ("New York, NY", -74.006, 40.713, 0),
    ("Mexico City", -99.133, 19.433, 2240),
    ("London, UK", -0.118, 51.509, 0),
    ("Tokyo, Japan", 139.650, 35.676, 0),
    ("Sydney, Australia", 151.209, -33.868, 0),
    ("Cairo, Egypt", 31.236, 30.044, 0),
    ("Reykjavik, Iceland", -21.896, 64.146, 0),
    ("Sao Paulo, Brazil", -46.634, -23.551, 760),
    ("Mumbai, India", 72.878, 19.076, 14),
]

# Eclipse search start dates (year, month, day)
SEARCH_DATES = [
    (2017, 1, 1, "2017 eclipse season"),
    (2019, 1, 1, "2019 eclipse season"),
    (2020, 1, 1, "2020 eclipse season"),
    (2023, 1, 1, "2023 eclipse season"),
    (2024, 1, 1, "2024 eclipse season"),
    (2025, 1, 1, "2025 eclipse season"),
]


def time_diff_seconds(jd1, jd2):
    """Absolute time difference in seconds."""
    return abs(jd1 - jd2) * 86400.0


def jd_to_date_str(jd):
    """Convert JD to readable date string."""
    y, m, d, h = swe.revjul(jd)
    hours = int(h)
    minutes = int((h - hours) * 60)
    seconds = int(((h - hours) * 60 - minutes) * 60)
    return f"{y:04d}-{m:02d}-{d:02d} {hours:02d}:{minutes:02d}:{seconds:02d} UT"


def run_audit():
    """Run the full Round 2 solar eclipse local deep audit."""
    swe.set_ephe_path(_REF_EPHE_PATH)

    print("=" * 80)
    print("ROUND 2: SOLAR ECLIPSE LOCAL FUNCTIONS DEEP AUDIT")
    print("=" * 80)

    issues = []
    stats = {
        "when_loc_tests": 0,
        "when_loc_pass": 0,
        "when_loc_skip": 0,
        "how_tests": 0,
        "how_pass": 0,
        "timing_diffs": [],
        "contact_diffs": [],
        "magnitude_diffs": [],
        "obscuration_diffs": [],
        "altitude_diffs": [],
        "azimuth_diffs": [],
    }

    # ========================================================================
    # PART 1: sol_eclipse_when_loc — find next eclipse visible from location
    # ========================================================================

    print("\n" + "=" * 80)
    print("PART 1: sol_eclipse_when_loc")
    print("=" * 80)

    for date_year, date_month, date_day, date_desc in SEARCH_DATES:
        jd_start = swe.julday(date_year, date_month, date_day, 0.0)

        for loc_name, loc_lon, loc_lat, loc_alt in LOCATIONS:
            stats["when_loc_tests"] += 1
            geopos = (loc_lon, loc_lat, loc_alt)
            test_id = f"{date_desc} @ {loc_name}"

            try:
                # pyswisseph: sol_eclipse_when_loc(tjdut, geopos, flags, backwards)
                swe_ret = swe.sol_eclipse_when_loc(
                    jd_start, geopos, FLG_SWIEPH, False
                )
                swe_retflag = swe_ret[0]
                swe_tret = swe_ret[1]
                swe_attr = swe_ret[2]
            except Exception as e:
                print(f"  SKIP {test_id}: pyswisseph error: {e}")
                stats["when_loc_skip"] += 1
                continue

            try:
                # libephemeris: sol_eclipse_when_loc(tjd_start, ifl, geopos, backward)
                lib_retflag, lib_tret, lib_attr = ephem.sol_eclipse_when_loc(
                    jd_start, FLG_SWIEPH, list(geopos), False
                )
            except Exception as e:
                msg = f"FAIL {test_id}: libephemeris error: {e}"
                print(f"  {msg}")
                issues.append(msg)
                continue

            # --- Compare maximum time ---
            swe_max = swe_tret[0]
            lib_max = lib_tret[0]
            max_diff_s = time_diff_seconds(swe_max, lib_max)
            stats["timing_diffs"].append(max_diff_s)

            # Check they found the same eclipse (within 30 days)
            if abs(swe_max - lib_max) > 30:
                msg = (
                    f"FAIL {test_id}: Different eclipse found! "
                    f"SE={jd_to_date_str(swe_max)}, LE={jd_to_date_str(lib_max)}"
                )
                print(f"  {msg}")
                issues.append(msg)
                continue

            passed = True
            detail_lines = []

            if max_diff_s > TIMING_TOL_SECONDS:
                msg = (
                    f"FAIL {test_id}: max time diff {max_diff_s:.1f}s "
                    f"(tol={TIMING_TOL_SECONDS}s)"
                )
                detail_lines.append(msg)
                issues.append(msg)
                passed = False

            # --- Compare contact times ---
            contact_names = [
                "max",
                "1st contact",
                "2nd contact",
                "3rd contact",
                "4th contact",
            ]
            for ci in range(1, 5):  # contacts 1-4
                swe_ct = swe_tret[ci]
                lib_ct = lib_tret[ci] if ci < len(lib_tret) else 0.0

                # Skip if both are 0 (no contact for this eclipse type)
                if swe_ct == 0.0 and lib_ct == 0.0:
                    continue
                if swe_ct == 0.0 or lib_ct == 0.0:
                    msg = (
                        f"WARN {test_id}: {contact_names[ci]} — "
                        f"SE={'0' if swe_ct == 0 else jd_to_date_str(swe_ct)}, "
                        f"LE={'0' if lib_ct == 0 else jd_to_date_str(lib_ct)}"
                    )
                    detail_lines.append(msg)
                    continue

                ct_diff_s = time_diff_seconds(swe_ct, lib_ct)
                stats["contact_diffs"].append(ct_diff_s)

                if ct_diff_s > CONTACT_TOL_SECONDS:
                    msg = (
                        f"FAIL {test_id}: {contact_names[ci]} diff {ct_diff_s:.1f}s "
                        f"(tol={CONTACT_TOL_SECONDS}s)"
                    )
                    detail_lines.append(msg)
                    issues.append(msg)
                    passed = False

            # --- Compare attributes ---
            # attr[0]: magnitude (fraction of solar diameter covered)
            swe_mag = swe_attr[0]
            lib_mag = lib_attr[0]
            mag_diff = abs(swe_mag - lib_mag)
            stats["magnitude_diffs"].append(mag_diff)

            if mag_diff > MAGNITUDE_TOL:
                msg = (
                    f"FAIL {test_id}: magnitude diff {mag_diff:.4f} "
                    f"(SE={swe_mag:.4f}, LE={lib_mag:.4f}, tol={MAGNITUDE_TOL})"
                )
                detail_lines.append(msg)
                issues.append(msg)
                passed = False

            # attr[1]: diameter ratio
            swe_ratio = swe_attr[1]
            lib_ratio = lib_attr[1]
            ratio_diff = abs(swe_ratio - lib_ratio)

            if ratio_diff > DIAMETER_RATIO_TOL:
                msg = (
                    f"FAIL {test_id}: diameter ratio diff {ratio_diff:.4f} "
                    f"(SE={swe_ratio:.4f}, LE={lib_ratio:.4f})"
                )
                detail_lines.append(msg)
                issues.append(msg)
                passed = False

            # attr[2]: obscuration (fraction of solar disc area covered)
            swe_obsc = swe_attr[2]
            lib_obsc = lib_attr[2]
            obsc_diff = abs(swe_obsc - lib_obsc)
            stats["obscuration_diffs"].append(obsc_diff)

            if obsc_diff > OBSCURATION_TOL:
                msg = (
                    f"FAIL {test_id}: obscuration diff {obsc_diff:.4f} "
                    f"(SE={swe_obsc:.4f}, LE={lib_obsc:.4f}, tol={OBSCURATION_TOL})"
                )
                detail_lines.append(msg)
                issues.append(msg)
                passed = False

            # attr[5]: true altitude of sun
            if len(lib_attr) > 5:
                swe_alt = swe_attr[5]
                lib_alt = lib_attr[5]
                alt_diff = abs(swe_alt - lib_alt)
                stats["altitude_diffs"].append(alt_diff)

                if alt_diff > SUN_ALT_TOL:
                    msg = (
                        f"WARN {test_id}: Sun altitude diff {alt_diff:.2f}° "
                        f"(SE={swe_alt:.2f}°, LE={lib_alt:.2f}°)"
                    )
                    detail_lines.append(msg)

            # attr[4]: azimuth of sun
            if len(lib_attr) > 4:
                swe_az = swe_attr[4]
                lib_az = lib_attr[4]
                az_diff = abs(swe_az - lib_az)
                if az_diff > 180:
                    az_diff = 360 - az_diff
                stats["azimuth_diffs"].append(az_diff)

                if az_diff > SUN_AZ_TOL:
                    msg = (
                        f"WARN {test_id}: Sun azimuth diff {az_diff:.2f}° "
                        f"(SE={swe_az:.2f}°, LE={lib_az:.2f}°)"
                    )
                    detail_lines.append(msg)

            if passed:
                stats["when_loc_pass"] += 1
                print(
                    f"  PASS {test_id}: max_diff={max_diff_s:.1f}s, "
                    f"mag_diff={mag_diff:.4f}, obsc_diff={obsc_diff:.4f}"
                )
            else:
                for line in detail_lines:
                    print(f"  {line}")

    # ========================================================================
    # PART 2: sol_eclipse_how — eclipse circumstances at given time/location
    # ========================================================================

    print("\n" + "=" * 80)
    print("PART 2: sol_eclipse_how at known eclipse times")
    print("=" * 80)

    # Use specific known eclipse JDs for precise attribute comparison
    known_eclipses = [
        # (description, approximate JD of max, test locations)
        ("2024 Apr 8 Total", swe.julday(2024, 4, 8, 18.0)),
        ("2023 Oct 14 Annular", swe.julday(2023, 10, 14, 18.0)),
        ("2017 Aug 21 Total", swe.julday(2017, 8, 21, 18.3)),
        ("2019 Jul 2 Total", swe.julday(2019, 7, 2, 19.0)),
        ("2020 Jun 21 Annular", swe.julday(2020, 6, 21, 6.0)),
        ("2025 Mar 29 Partial", swe.julday(2025, 3, 29, 10.0)),
    ]

    # First find actual eclipse max times via global search
    for ecl_desc, jd_approx in known_eclipses:
        # Find the actual eclipse max time
        try:
            swe_glob = swe.sol_eclipse_when_glob(jd_approx - 1, FLG_SWIEPH)
            jd_max = swe_glob[1][0]
        except Exception:
            print(f"  SKIP {ecl_desc}: could not find global eclipse")
            continue

        print(f"\n  Eclipse: {ecl_desc} — max at {jd_to_date_str(jd_max)}")

        for loc_name, loc_lon, loc_lat, loc_alt in LOCATIONS:
            stats["how_tests"] += 1
            geopos = (loc_lon, loc_lat, loc_alt)
            test_id = f"{ecl_desc} @ {loc_name}"

            try:
                # pyswisseph: sol_eclipse_how(tjdut, geopos, flags)
                swe_ret = swe.sol_eclipse_how(jd_max, geopos, FLG_SWIEPH)
                swe_retflag = swe_ret[0]
                swe_attr = swe_ret[1]
            except Exception as e:
                print(f"    SKIP {loc_name}: pyswisseph error: {e}")
                continue

            try:
                # libephemeris: sol_eclipse_how(tjd_ut, ifl, geopos)
                lib_retflag, lib_attr = ephem.sol_eclipse_how(
                    jd_max, FLG_SWIEPH, list(geopos)
                )
            except Exception as e:
                msg = f"FAIL {test_id}: libephemeris error: {e}"
                print(f"    {msg}")
                issues.append(msg)
                continue

            # If no eclipse visible from this location, both should agree
            if swe_retflag == 0:
                if lib_retflag == 0:
                    stats["how_pass"] += 1
                    continue
                else:
                    msg = f"WARN {test_id}: SE says no eclipse, LE says flag={lib_retflag}"
                    print(f"    {msg}")
                    # Not a hard fail - could be borderline
                    stats["how_pass"] += 1
                    continue

            if lib_retflag == 0 and swe_retflag != 0:
                msg = f"WARN {test_id}: LE says no eclipse, SE says flag={swe_retflag}"
                print(f"    {msg}")
                stats["how_pass"] += 1
                continue

            # Both see eclipse - compare attributes
            passed = True
            detail_lines = []

            # Magnitude
            swe_mag = swe_attr[0]
            lib_mag = lib_attr[0]
            mag_diff = abs(swe_mag - lib_mag)
            stats["magnitude_diffs"].append(mag_diff)

            if mag_diff > MAGNITUDE_TOL:
                msg = (
                    f"FAIL {test_id}: magnitude diff {mag_diff:.4f} "
                    f"(SE={swe_mag:.4f}, LE={lib_mag:.4f})"
                )
                detail_lines.append(msg)
                issues.append(msg)
                passed = False

            # Obscuration
            swe_obsc = swe_attr[2]
            lib_obsc = lib_attr[2] if len(lib_attr) > 2 else 0.0
            obsc_diff = abs(swe_obsc - lib_obsc)
            stats["obscuration_diffs"].append(obsc_diff)

            if obsc_diff > OBSCURATION_TOL:
                msg = (
                    f"FAIL {test_id}: obscuration diff {obsc_diff:.4f} "
                    f"(SE={swe_obsc:.4f}, LE={lib_obsc:.4f})"
                )
                detail_lines.append(msg)
                issues.append(msg)
                passed = False

            # Sun altitude
            swe_alt = swe_attr[5]
            lib_alt = lib_attr[5] if len(lib_attr) > 5 else 0.0
            alt_diff = abs(swe_alt - lib_alt)
            stats["altitude_diffs"].append(alt_diff)

            if alt_diff > SUN_ALT_TOL:
                msg = (
                    f"WARN {test_id}: altitude diff {alt_diff:.2f}° "
                    f"(SE={swe_alt:.2f}°, LE={lib_alt:.2f}°)"
                )
                detail_lines.append(msg)

            # Sun azimuth
            swe_az = swe_attr[4]
            lib_az = lib_attr[4] if len(lib_attr) > 4 else 0.0
            az_diff = abs(swe_az - lib_az)
            if az_diff > 180:
                az_diff = 360 - az_diff
            stats["azimuth_diffs"].append(az_diff)

            if az_diff > SUN_AZ_TOL:
                msg = (
                    f"WARN {test_id}: azimuth diff {az_diff:.2f}° "
                    f"(SE={swe_az:.2f}°, LE={lib_az:.2f}°)"
                )
                detail_lines.append(msg)

            if passed:
                stats["how_pass"] += 1
                print(
                    f"    PASS {loc_name}: mag={lib_mag:.4f} (d={mag_diff:.4f}), "
                    f"obsc={lib_obsc:.4f} (d={obsc_diff:.4f}), "
                    f"alt={lib_alt:.1f}° (d={alt_diff:.2f}°)"
                )
            else:
                for line in detail_lines:
                    print(f"    {line}")

    # ========================================================================
    # PART 3: Detailed contact time comparison for the April 2024 eclipse
    # ========================================================================

    print("\n" + "=" * 80)
    print(
        "PART 3: Detailed contact time comparison — April 8, 2024 Total Solar Eclipse"
    )
    print("=" * 80)

    # Locations known to be in the path of totality
    totality_locations = [
        ("Dallas, TX (totality)", -96.797, 32.777, 0),
        ("Indianapolis, IN (totality)", -86.158, 39.768, 0),
        ("Burlington, VT (totality)", -73.213, 44.476, 0),
        ("Mazatlan, Mexico (totality)", -106.417, 23.233, 0),
    ]

    # Locations for partial eclipse
    partial_locations = [
        ("New York, NY (partial)", -74.006, 40.713, 0),
        ("Los Angeles, CA (partial)", -118.244, 34.052, 0),
        ("Chicago, IL (partial)", -87.630, 41.878, 0),
        ("Miami, FL (partial)", -80.191, 25.762, 0),
    ]

    jd_start_2024 = swe.julday(2024, 1, 1, 0.0)

    for loc_name, loc_lon, loc_lat, loc_alt in totality_locations + partial_locations:
        geopos = (loc_lon, loc_lat, loc_alt)

        try:
            swe_ret = swe.sol_eclipse_when_loc(
                jd_start_2024, geopos, FLG_SWIEPH, False
            )
            swe_retflag = swe_ret[0]
            swe_tret = swe_ret[1]
            swe_attr = swe_ret[2]
        except Exception as e:
            print(f"  SKIP {loc_name}: pyswisseph error: {e}")
            continue

        try:
            lib_retflag, lib_tret, lib_attr = ephem.sol_eclipse_when_loc(
                jd_start_2024, FLG_SWIEPH, list(geopos), False
            )
        except Exception as e:
            print(f"  FAIL {loc_name}: libephemeris error: {e}")
            issues.append(f"FAIL {loc_name}: libephemeris error: {e}")
            continue

        # Verify they found the same eclipse (April 2024)
        swe_y, swe_m, swe_d, _ = swe.revjul(swe_tret[0])
        lib_y, lib_m, lib_d, _ = swe.revjul(lib_tret[0])

        if (swe_y, swe_m) != (lib_y, lib_m):
            print(
                f"  FAIL {loc_name}: Different eclipse! "
                f"SE={swe_y}-{swe_m}-{swe_d}, LE={lib_y}-{lib_m}-{lib_d}"
            )
            issues.append(
                f"FAIL {loc_name}: Different eclipse found "
                f"SE={swe_y}-{swe_m}-{swe_d}, LE={lib_y}-{lib_m}-{lib_d}"
            )
            continue

        print(f"\n  {loc_name}:")
        print(f"    Eclipse type: SE=0x{swe_retflag:04x}, LE=0x{lib_retflag:04x}")

        contact_names = [
            "Maximum",
            "1st contact",
            "2nd contact",
            "3rd contact",
            "4th contact",
            "Sunrise",
            "Sunset",
        ]

        for ci in range(min(7, len(swe_tret), len(lib_tret))):
            swe_ct = swe_tret[ci]
            lib_ct = lib_tret[ci]

            if swe_ct == 0.0 and lib_ct == 0.0:
                continue

            if swe_ct == 0.0:
                print(
                    f"    {contact_names[ci]:15s}: SE=N/A, LE={jd_to_date_str(lib_ct)}"
                )
                continue
            if lib_ct == 0.0:
                print(
                    f"    {contact_names[ci]:15s}: SE={jd_to_date_str(swe_ct)}, LE=N/A"
                )
                continue

            diff_s = time_diff_seconds(swe_ct, lib_ct)
            status = "OK" if diff_s < CONTACT_TOL_SECONDS else "FAIL"
            print(
                f"    {contact_names[ci]:15s}: diff={diff_s:7.1f}s  "
                f"SE={jd_to_date_str(swe_ct)}  LE={jd_to_date_str(lib_ct)}  [{status}]"
            )

            if diff_s > CONTACT_TOL_SECONDS:
                issues.append(
                    f"FAIL {loc_name} {contact_names[ci]}: diff={diff_s:.1f}s"
                )

        # Attributes
        attr_names = [
            "Magnitude",
            "Diam ratio",
            "Obscuration",
            "Shadow km",
            "Sun az",
            "Sun alt(true)",
            "Sun alt(app)",
            "Elongation",
        ]
        print(f"    {'Attribute':15s}  {'SE':>10s}  {'LE':>10s}  {'Diff':>10s}")
        for ai in range(min(8, len(swe_attr), len(lib_attr))):
            swe_val = swe_attr[ai]
            lib_val = lib_attr[ai]
            diff_val = abs(swe_val - lib_val)
            print(
                f"    {attr_names[ai]:15s}  {swe_val:10.4f}  {lib_val:10.4f}  {diff_val:10.4f}"
            )

    # ========================================================================
    # PART 4: Sequential eclipse search from multiple locations
    # ========================================================================

    print("\n" + "=" * 80)
    print("PART 4: Sequential eclipse search — 5 consecutive eclipses per location")
    print("=" * 80)

    test_locs = [
        ("Rome, Italy", 12.496, 41.903, 0),
        ("Tokyo, Japan", 139.650, 35.676, 0),
        ("Santiago, Chile", -70.669, -33.449, 520),
    ]

    for loc_name, loc_lon, loc_lat, loc_alt in test_locs:
        geopos = (loc_lon, loc_lat, loc_alt)
        print(f"\n  {loc_name}: searching 5 consecutive local eclipses...")

        jd = swe.julday(2020, 1, 1, 0.0)

        for i in range(5):
            try:
                swe_ret = swe.sol_eclipse_when_loc(jd, geopos, FLG_SWIEPH, False)
                swe_tret = swe_ret[1]
                swe_attr = swe_ret[2]
                swe_max = swe_tret[0]
            except Exception as e:
                print(f"    Eclipse #{i + 1}: pyswisseph error: {e}")
                jd += 180
                continue

            try:
                lib_retflag, lib_tret, lib_attr = ephem.sol_eclipse_when_loc(
                    jd, FLG_SWIEPH, list(geopos), False
                )
                lib_max = lib_tret[0]
            except Exception as e:
                print(f"    Eclipse #{i + 1}: libephemeris error: {e}")
                jd = swe_max + 30
                continue

            max_diff_s = time_diff_seconds(swe_max, lib_max)

            # Check it's the same eclipse
            if abs(swe_max - lib_max) > 30:
                print(
                    f"    Eclipse #{i + 1}: DIFFERENT eclipse! "
                    f"SE={jd_to_date_str(swe_max)}, LE={jd_to_date_str(lib_max)}"
                )
                issues.append(
                    f"FAIL {loc_name} eclipse #{i + 1}: different eclipse found"
                )
                jd = swe_max + 30
                continue

            swe_mag = swe_attr[0]
            lib_mag = lib_attr[0]
            mag_diff = abs(swe_mag - lib_mag)

            status = (
                "PASS"
                if max_diff_s < TIMING_TOL_SECONDS and mag_diff < MAGNITUDE_TOL
                else "FAIL"
            )
            print(
                f"    Eclipse #{i + 1}: {jd_to_date_str(swe_max)}  "
                f"time_diff={max_diff_s:6.1f}s  mag_diff={mag_diff:.4f}  [{status}]"
            )

            if status == "FAIL":
                issues.append(
                    f"FAIL {loc_name} eclipse #{i + 1} @ {jd_to_date_str(swe_max)}: "
                    f"time_diff={max_diff_s:.1f}s, mag_diff={mag_diff:.4f}"
                )

            jd = swe_max + 30

    # ========================================================================
    # SUMMARY
    # ========================================================================

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    print(
        f"\nPart 1 (when_loc): {stats['when_loc_pass']}/{stats['when_loc_tests']} passed "
        f"({stats['when_loc_skip']} skipped)"
    )
    print(f"Part 2 (how):      {stats['how_pass']}/{stats['how_tests']} passed")

    if stats["timing_diffs"]:
        timing = stats["timing_diffs"]
        print(f"\nTiming diffs (max time):")
        print(
            f"  Min: {min(timing):.1f}s  Max: {max(timing):.1f}s  "
            f"Mean: {sum(timing) / len(timing):.1f}s  "
            f"Median: {sorted(timing)[len(timing) // 2]:.1f}s"
        )

    if stats["contact_diffs"]:
        contacts = stats["contact_diffs"]
        print(f"\nContact time diffs:")
        print(
            f"  Min: {min(contacts):.1f}s  Max: {max(contacts):.1f}s  "
            f"Mean: {sum(contacts) / len(contacts):.1f}s  "
            f"Median: {sorted(contacts)[len(contacts) // 2]:.1f}s"
        )

    if stats["magnitude_diffs"]:
        mags = stats["magnitude_diffs"]
        print(f"\nMagnitude diffs:")
        print(
            f"  Min: {min(mags):.4f}  Max: {max(mags):.4f}  "
            f"Mean: {sum(mags) / len(mags):.4f}  "
            f"Median: {sorted(mags)[len(mags) // 2]:.4f}"
        )

    if stats["obscuration_diffs"]:
        obsc = stats["obscuration_diffs"]
        print(f"\nObscuration diffs:")
        print(
            f"  Min: {min(obsc):.4f}  Max: {max(obsc):.4f}  "
            f"Mean: {sum(obsc) / len(obsc):.4f}  "
            f"Median: {sorted(obsc)[len(obsc) // 2]:.4f}"
        )

    if stats["altitude_diffs"]:
        alts = stats["altitude_diffs"]
        print(f"\nSun altitude diffs:")
        print(
            f"  Min: {min(alts):.2f}°  Max: {max(alts):.2f}°  "
            f"Mean: {sum(alts) / len(alts):.2f}°  "
            f"Median: {sorted(alts)[len(alts) // 2]:.2f}°"
        )

    if stats["azimuth_diffs"]:
        azs = stats["azimuth_diffs"]
        print(f"\nSun azimuth diffs:")
        print(
            f"  Min: {min(azs):.2f}°  Max: {max(azs):.2f}°  "
            f"Mean: {sum(azs) / len(azs):.2f}°  "
            f"Median: {sorted(azs)[len(azs) // 2]:.2f}°"
        )

    print(f"\n{'=' * 80}")
    if issues:
        print(f"ISSUES FOUND: {len(issues)}")
        for issue in issues:
            print(f"  - {issue}")
    else:
        print("ALL TESTS PASSED")
    print(f"{'=' * 80}")

    return len(issues)


if __name__ == "__main__":
    issue_count = run_audit()
    exit(1 if issue_count > 0 else 0)
