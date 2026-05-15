#!/usr/bin/env python3
"""
Round 7: Deep Heliacal Events Audit
====================================
Compares libephemeris heliacal functions against pyswisseph:
  - heliacal_ut: timing of heliacal events (rising/setting/evening first/morning last)
  - heliacal_pheno_ut: detailed heliacal phenomena (50-element array)
  - vis_limit_mag: visual limiting magnitude

Tests multiple planets, stars, locations, and event types.
"""

from __future__ import annotations

import os
import sys
import time
import traceback

import swisseph as swe

# Force Skyfield mode
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import libephemeris as ephem

# Set SE ephemeris path
EPHE_PATH = os.path.join(os.path.dirname(__file__), "..", "swisseph", "ephe")
if os.path.exists(EPHE_PATH):
    swe.set_ephe_path(EPHE_PATH)

# ============================================================================
# CONSTANTS
# ============================================================================

HELIACAL_RISING = 1
HELIACAL_SETTING = 2
EVENING_FIRST = 3
MORNING_LAST = 4

EVENT_NAMES = {
    HELIACAL_RISING: "heliacal_rising",
    HELIACAL_SETTING: "heliacal_setting",
    EVENING_FIRST: "evening_first",
    MORNING_LAST: "morning_last",
}

# Planets for testing
PLANETS = [
    ("Mercury", swe.MERCURY),
    ("Venus", swe.VENUS),
    ("Mars", swe.MARS),
    ("Jupiter", swe.JUPITER),
    ("Saturn", swe.SATURN),
]

INNER_PLANETS = ["Mercury", "Venus"]

# Stars for testing
STARS = ["Sirius", "Arcturus", "Vega", "Aldebaran", "Regulus"]

# Locations (lon, lat, alt)
LOCATIONS = [
    ("Cairo", 31.24, 30.04, 75),
    ("Rome", 12.50, 41.90, 0),
    ("London", -0.12, 51.51, 0),
    ("Equator", 0.0, 0.0, 0),
]

# Standard atmospheric and observer conditions
STANDARD_ATMO = (1013.25, 15.0, 50.0, 0.25)
STANDARD_OBSERVER = (25.0, 1.0, 1, 1, 0, 0)

# Start dates
START_DATES = [
    ("2024-Jan", swe.julday(2024, 1, 1, 0.0)),
    ("2024-Jul", swe.julday(2024, 7, 1, 0.0)),
]

# ============================================================================
# COUNTERS
# ============================================================================

total = 0
passed = 0
failed = 0
skipped = 0
errors = 0
failures = []


def record(test_name, pass_fail, detail=""):
    global total, passed, failed, failures
    total += 1
    if pass_fail:
        passed += 1
        status = "PASS"
    else:
        failed += 1
        failures.append((test_name, detail))
        status = "FAIL"
    print(f"  [{status}] {test_name}: {detail}")


def record_skip(test_name, reason=""):
    global skipped
    skipped += 1
    print(f"  [SKIP] {test_name}: {reason}")


def record_error(test_name, error_msg):
    global errors
    errors += 1
    print(f"  [ERROR] {test_name}: {error_msg}")


# ============================================================================
# PART 1: heliacal_ut - Event Timing
# ============================================================================


def test_part1_heliacal_ut():
    """Test heliacal_ut event timing for planets and stars."""
    print("\n" + "=" * 70)
    print("PART 1: heliacal_ut — Event Timing")
    print("=" * 70)

    # Tolerances
    PLANET_TOL_DAYS = 3.0  # Start generous to characterize differences
    STAR_TOL_DAYS = 5.0

    for date_name, jd_start in START_DATES:
        for loc_name, lon, lat, alt in LOCATIONS:
            geopos = (lon, lat, alt)

            # Test planets
            for planet_name, planet_id in PLANETS:
                # Determine which event types to test
                event_types = [HELIACAL_RISING, HELIACAL_SETTING]
                if planet_name in INNER_PLANETS:
                    event_types.extend([EVENING_FIRST, MORNING_LAST])

                for event_type in event_types:
                    test_name = (
                        f"P1/{planet_name}/{EVENT_NAMES[event_type]}"
                        f"/{loc_name}/{date_name}"
                    )
                    try:
                        # pyswisseph
                        ret_se = swe.heliacal_ut(
                            jd_start,
                            geopos,
                            STANDARD_ATMO,
                            STANDARD_OBSERVER,
                            planet_name,
                            event_type,
                            0,
                        )
                        jd_se = ret_se[0]

                        if jd_se <= 0:
                            record_skip(test_name, "SE returned 0")
                            continue

                        # libephemeris
                        ret_le = ephem.heliacal_ut(
                            jd_start,
                            geopos,
                            STANDARD_ATMO,
                            STANDARD_OBSERVER,
                            planet_name,
                            event_type,
                        )
                        jd_le = ret_le[0]

                        if jd_le <= 0:
                            record(test_name, False, "LE returned 0, SE found event")
                            continue

                        diff_days = abs(jd_se - jd_le)
                        diff_hours = diff_days * 24

                        ok = diff_days < PLANET_TOL_DAYS
                        record(
                            test_name,
                            ok,
                            f"SE={jd_se:.5f} LE={jd_le:.5f} "
                            f"diff={diff_days:.3f}d ({diff_hours:.1f}h)",
                        )

                        # Also check jd2 and jd3 (optimum/end) vs SE
                        if len(ret_se) >= 3 and ret_se[1] > 0:
                            jd2_se = ret_se[1]
                            jd3_se = ret_se[2]
                            jd2_le = ret_le[1]
                            jd3_le = ret_le[2]

                            # Check if LE returns distinct jd2/jd3
                            if jd2_le == jd_le and jd3_le == jd_le:
                                record(
                                    f"{test_name}/jd2_jd3_stub",
                                    False,
                                    f"LE returns jd2=jd3=jd1 (stub). "
                                    f"SE jd2-jd1={abs(jd2_se - ret_se[0]) * 86400:.1f}s "
                                    f"jd3-jd1={abs(jd3_se - ret_se[0]) * 86400:.1f}s",
                                )

                    except Exception as e:
                        record_error(test_name, str(e))

            # Test stars (only rising and setting)
            for star_name in STARS:
                for event_type in [HELIACAL_RISING, HELIACAL_SETTING]:
                    test_name = (
                        f"P1/{star_name}/{EVENT_NAMES[event_type]}"
                        f"/{loc_name}/{date_name}"
                    )
                    try:
                        # pyswisseph
                        ret_se = swe.heliacal_ut(
                            jd_start,
                            geopos,
                            STANDARD_ATMO,
                            STANDARD_OBSERVER,
                            star_name,
                            event_type,
                            0,
                        )
                        jd_se = ret_se[0]

                        if jd_se <= 0:
                            record_skip(test_name, "SE returned 0")
                            continue

                        # libephemeris
                        ret_le = ephem.heliacal_ut(
                            jd_start,
                            geopos,
                            STANDARD_ATMO,
                            STANDARD_OBSERVER,
                            star_name,
                            event_type,
                        )
                        jd_le = ret_le[0]

                        if jd_le <= 0:
                            record(test_name, False, "LE returned 0, SE found event")
                            continue

                        diff_days = abs(jd_se - jd_le)
                        diff_hours = diff_days * 24

                        ok = diff_days < STAR_TOL_DAYS
                        record(
                            test_name,
                            ok,
                            f"SE={jd_se:.5f} LE={jd_le:.5f} "
                            f"diff={diff_days:.3f}d ({diff_hours:.1f}h)",
                        )

                    except Exception as e:
                        record_error(test_name, str(e))


# ============================================================================
# PART 2: heliacal_pheno_ut — Phenomena Details
# ============================================================================


def test_part2_heliacal_pheno():
    """Test heliacal_pheno_ut return values."""
    print("\n" + "=" * 70)
    print("PART 2: heliacal_pheno_ut — Phenomena Details")
    print("=" * 70)

    # Use a specific date and location for detailed comparison
    jd = swe.julday(2024, 6, 15, 12.0)
    loc_name = "Cairo"
    geopos = (31.24, 30.04, 75)

    PHENO_LABELS = {
        0: "AltO (topo alt)",
        1: "AppAltO (app alt)",
        2: "GeoAltO (geo alt)",
        3: "AziO (azimuth)",
        4: "AltS (Sun alt)",
        5: "AziS (Sun az)",
        6: "TAVact (topo AV)",
        7: "ARCVact (geo AV)",
        8: "DAZact (az diff)",
        9: "ARCLact (elong)",
        10: "kact (extinction)",
        11: "minTAV",
        12: "TfirstVR",
        13: "TbVR",
        14: "TlastVR",
        15: "TbYallop",
        16: "WMoon",
        17: "qYal",
        18: "qCrit",
        19: "ParO (parallax)",
        20: "Magn",
        21: "RiseO",
        22: "RiseS",
        23: "Lag",
        24: "TvisVR",
        25: "LMoon",
        26: "CVAact (phase angle)",
        27: "Illum",
    }

    # Tolerances per field
    PHENO_TOL = {
        0: 1.0,  # altitude: 1 deg
        1: 1.0,  # apparent alt
        2: 1.0,  # geo alt
        3: 2.0,  # azimuth: 2 deg
        4: 1.0,  # sun alt
        5: 2.0,  # sun az
        6: 1.0,  # topo AV
        7: 1.0,  # geo AV
        8: 2.0,  # az diff
        9: 1.0,  # elongation
        10: 0.1,  # extinction coeff
        19: 0.01,  # parallax
        20: 1.0,  # magnitude
        26: 2.0,  # phase angle
        27: 5.0,  # illumination %
    }

    for planet_name in ["Venus", "Mercury", "Mars", "Jupiter", "Saturn"]:
        for event_type in [HELIACAL_RISING]:
            test_base = f"P2/{planet_name}/{EVENT_NAMES[event_type]}/{loc_name}"
            try:
                # pyswisseph - returns flat 50-tuple
                ret_se = swe.heliacal_pheno_ut(
                    jd,
                    geopos,
                    STANDARD_ATMO,
                    STANDARD_OBSERVER,
                    planet_name,
                    event_type,
                    0,
                )
                data_se = ret_se  # flat tuple

                # libephemeris - returns (50-tuple, retflag)
                ret_le = ephem.heliacal_pheno_ut(
                    jd,
                    geopos,
                    STANDARD_ATMO,
                    STANDARD_OBSERVER,
                    planet_name,
                    event_type,
                )
                # Check API shape
                if isinstance(ret_le, tuple) and len(ret_le) == 2:
                    data_le, retflag_le = ret_le
                    record(
                        f"{test_base}/api_shape",
                        False,
                        "LE returns (data, retflag), SE returns flat tuple — API mismatch",
                    )
                else:
                    data_le = ret_le
                    record(f"{test_base}/api_shape", True, "API shape matches")

                # Compare field by field
                for idx in range(min(28, len(data_se), len(data_le))):
                    val_se = data_se[idx]
                    val_le = data_le[idx]
                    label = PHENO_LABELS.get(idx, f"field[{idx}]")

                    # Skip reserved fields (both zero)
                    if val_se == 0.0 and val_le == 0.0:
                        continue

                    # Skip 99999999.0 sentinel values
                    if val_se == 99999999.0 or val_le == 99999999.0:
                        if val_se == 99999999.0 and val_le != 99999999.0:
                            record(
                                f"{test_base}/{label}",
                                False,
                                f"SE=99999999 (sentinel), LE={val_le:.6f}",
                            )
                        elif val_le == 99999999.0 and val_se != 99999999.0:
                            record(
                                f"{test_base}/{label}",
                                False,
                                f"SE={val_se:.6f}, LE=99999999 (sentinel)",
                            )
                        continue

                    tol = PHENO_TOL.get(idx)
                    if tol is not None:
                        diff = abs(val_se - val_le)
                        ok = diff < tol
                        record(
                            f"{test_base}/{label}",
                            ok,
                            f"SE={val_se:.6f} LE={val_le:.6f} diff={diff:.6f} tol={tol}",
                        )

            except Exception as e:
                record_error(test_base, f"{e}\n{traceback.format_exc()}")


# ============================================================================
# PART 3: vis_limit_mag — Visual Limiting Magnitude
# ============================================================================


def test_part3_vis_limit_mag():
    """Test vis_limit_mag return values."""
    print("\n" + "=" * 70)
    print("PART 3: vis_limit_mag — Visual Limiting Magnitude")
    print("=" * 70)

    # Test at different times of night
    test_times = [
        ("midnight", swe.julday(2024, 6, 15, 22.0)),
        ("dusk_30N", swe.julday(2024, 6, 15, 19.5)),  # ~sunset at 30N
        ("dawn_30N", swe.julday(2024, 6, 16, 3.0)),  # ~pre-dawn at 30N
    ]

    objects = ["Venus", "Jupiter", "Mars", "Saturn", "Sirius", "Arcturus"]

    for time_name, jd in test_times:
        for loc_name, lon, lat, alt in LOCATIONS[:2]:
            geopos = (lon, lat, alt)

            for obj_name in objects:
                test_name = f"P3/{obj_name}/{loc_name}/{time_name}"
                try:
                    # pyswisseph
                    ret_se = swe.vis_limit_mag(
                        jd,
                        geopos,
                        STANDARD_ATMO,
                        STANDARD_OBSERVER,
                        obj_name,
                        0,
                    )
                    flag_se = int(ret_se[0])
                    dret_se = ret_se[1]

                    # libephemeris
                    ret_le = ephem.vis_limit_mag(
                        jd,
                        geopos,
                        STANDARD_ATMO,
                        STANDARD_OBSERVER,
                        obj_name,
                        0,
                    )
                    flag_le = ret_le[0]
                    dret_le = ret_le[1]

                    # Check dret length
                    len_se = len(dret_se)
                    len_le = len(dret_le)
                    if len_se != len_le:
                        record(
                            f"{test_name}/dret_len",
                            False,
                            f"SE returns {len_se} elements, LE returns {len_le}",
                        )

                    # Compare return flag
                    flag_ok = flag_se == flag_le
                    record(f"{test_name}/flag", flag_ok, f"SE={flag_se} LE={flag_le}")

                    # If below horizon, SE returns all zeros — skip data comparison
                    if flag_se == -2:
                        if flag_le != -2:
                            record(
                                f"{test_name}/below_horizon",
                                False,
                                "SE says below horizon, LE disagrees",
                            )
                        continue

                    # Compare key fields
                    field_labels = [
                        "lim_mag",
                        "obj_alt",
                        "obj_az",
                        "sun_alt",
                        "sun_az",
                        "moon_alt",
                        "moon_az",
                        "obj_mag",
                    ]

                    field_tols = [1.5, 0.5, 1.0, 0.5, 1.0, 0.5, 1.0, 1.0]

                    for i, (label, tol) in enumerate(zip(field_labels, field_tols)):
                        if i >= min(len_se, len_le):
                            break
                        val_se = dret_se[i]
                        val_le = dret_le[i]

                        # Skip if both are sentinel/zero
                        if val_se == -100.0 and val_le == 0.0:
                            record(
                                f"{test_name}/{label}",
                                False,
                                f"SE=-100.0 (sentinel), LE=0.0",
                            )
                            continue

                        diff = abs(val_se - val_le)
                        ok = diff < tol
                        record(
                            f"{test_name}/{label}",
                            ok,
                            f"SE={val_se:.4f} LE={val_le:.4f} diff={diff:.4f}",
                        )

                except Exception as e:
                    record_error(test_name, str(e))


# ============================================================================
# PART 4: Flag Constants Verification
# ============================================================================


def test_part4_flag_constants():
    """Verify heliacal flag constants match pyswisseph."""
    print("\n" + "=" * 70)
    print("PART 4: Flag Constants Verification")
    print("=" * 70)

    # Mapping of our constants to pyswisseph constants
    flag_checks = [
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

    for name, le_val, se_val in flag_checks:
        test_name = f"P4/const/{name}"
        ok = le_val == se_val
        record(
            test_name,
            ok,
            f"LE={le_val} SE={se_val}" + ("" if ok else " *** MISMATCH ***"),
        )

    # Check for missing constants
    missing_flags = [
        ("HELFLAG_HIGH_PRECISION", swe.HELFLAG_HIGH_PRECISION),
        ("HELFLAG_LONG_SEARCH", swe.HELFLAG_LONG_SEARCH),
        ("HELFLAG_SEARCH_1_PERIOD", swe.HELFLAG_SEARCH_1_PERIOD),
        ("HELFLAG_VISLIM_PHOTOPIC", swe.HELFLAG_VISLIM_PHOTOPIC),
        ("HELFLAG_AV", swe.HELFLAG_AV),
        ("HELFLAG_AVKIND", swe.HELFLAG_AVKIND),
        ("HELFLAG_AVKIND_MIN7", swe.HELFLAG_AVKIND_MIN7),
        ("HELFLAG_AVKIND_MIN9", swe.HELFLAG_AVKIND_MIN9),
        ("HELFLAG_AVKIND_PTO", swe.HELFLAG_AVKIND_PTO),
        ("HELFLAG_AVKIND_VR", swe.HELFLAG_AVKIND_VR),
    ]

    for name, se_val in missing_flags:
        test_name = f"P4/missing/{name}"
        has_it = hasattr(ephem, f"SE_{name}") or hasattr(ephem, name)
        if has_it:
            le_val = getattr(ephem, f"SE_{name}", getattr(ephem, name, None))
            ok = le_val == se_val
            record(test_name, ok, f"LE={le_val} SE={se_val}")
        else:
            record(test_name, False, f"Missing in libephemeris (SE value={se_val})")


# ============================================================================
# PART 5: API Shape / Return Type Verification
# ============================================================================


def test_part5_api_shape():
    """Verify API return types and shapes match pyswisseph."""
    print("\n" + "=" * 70)
    print("PART 5: API Shape / Return Type Verification")
    print("=" * 70)

    jd_start = swe.julday(2024, 1, 1, 0.0)
    geopos = (0.0, 30.0, 0)

    # --- heliacal_ut ---
    print("\n  --- heliacal_ut ---")
    try:
        ret_se = swe.heliacal_ut(
            jd_start,
            geopos,
            STANDARD_ATMO,
            STANDARD_OBSERVER,
            "Venus",
            HELIACAL_RISING,
            0,
        )
        ret_le = ephem.heliacal_ut(
            jd_start,
            geopos,
            STANDARD_ATMO,
            STANDARD_OBSERVER,
            "Venus",
            HELIACAL_RISING,
        )

        # Type check
        record(
            "P5/heliacal_ut/type",
            isinstance(ret_le, tuple),
            f"SE type={type(ret_se).__name__} LE type={type(ret_le).__name__}",
        )

        # Length check
        record(
            "P5/heliacal_ut/length",
            len(ret_le) == len(ret_se),
            f"SE len={len(ret_se)} LE len={len(ret_le)}",
        )

        # jd2 != jd1 check (optimum should differ from start)
        if ret_se[1] > 0 and ret_se[1] != ret_se[0]:
            jd2_distinct = ret_le[1] != ret_le[0]
            record(
                "P5/heliacal_ut/jd2_distinct",
                jd2_distinct,
                f"SE jd2-jd1={abs(ret_se[1] - ret_se[0]) * 86400:.1f}s, "
                f"LE jd2-jd1={abs(ret_le[1] - ret_le[0]) * 86400:.1f}s",
            )

        if ret_se[2] > 0 and ret_se[2] != ret_se[0]:
            jd3_distinct = ret_le[2] != ret_le[0]
            record(
                "P5/heliacal_ut/jd3_distinct",
                jd3_distinct,
                f"SE jd3-jd1={abs(ret_se[2] - ret_se[0]) * 86400:.1f}s, "
                f"LE jd3-jd1={abs(ret_le[2] - ret_le[0]) * 86400:.1f}s",
            )

    except Exception as e:
        record_error("P5/heliacal_ut", str(e))

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
            HELIACAL_RISING,
            0,
        )
        ret_le = ephem.heliacal_pheno_ut(
            jd_pheno,
            geopos,
            STANDARD_ATMO,
            STANDARD_OBSERVER,
            "Venus",
            HELIACAL_RISING,
        )

        # SE returns flat 50-tuple; check what LE returns
        se_is_flat = isinstance(ret_se, tuple) and len(ret_se) == 50
        if isinstance(ret_le, tuple) and len(ret_le) == 2:
            # LE returns (data, retflag) — API mismatch
            le_is_nested = isinstance(ret_le[0], tuple)
            record(
                "P5/heliacal_pheno_ut/shape",
                False,
                f"SE returns flat 50-tuple, LE returns (data_tuple, retflag) — "
                f"API mismatch. SE flat={se_is_flat}, LE nested={le_is_nested}",
            )
            data_le = ret_le[0]
        elif isinstance(ret_le, tuple) and len(ret_le) == 50:
            record("P5/heliacal_pheno_ut/shape", True, "Both return flat 50-tuple")
            data_le = ret_le
        else:
            record(
                "P5/heliacal_pheno_ut/shape",
                False,
                f"Unexpected LE shape: type={type(ret_le).__name__} len={len(ret_le)}",
            )

    except Exception as e:
        record_error("P5/heliacal_pheno_ut", str(e))

    # --- vis_limit_mag ---
    print("\n  --- vis_limit_mag ---")
    jd_night = swe.julday(2024, 8, 15, 22.0)
    geopos_vlm = (12.5, 41.9, 0)  # Rome
    try:
        ret_se = swe.vis_limit_mag(
            jd_night,
            geopos_vlm,
            STANDARD_ATMO,
            STANDARD_OBSERVER,
            "Jupiter",
            0,
        )
        ret_le = ephem.vis_limit_mag(
            jd_night,
            geopos_vlm,
            STANDARD_ATMO,
            STANDARD_OBSERVER,
            "Jupiter",
            0,
        )

        # SE returns (flag, 10-tuple)
        record(
            "P5/vis_limit_mag/type",
            isinstance(ret_le, tuple) and len(ret_le) == 2,
            f"SE: ({type(ret_se[0]).__name__}, {len(ret_se[1])}-tuple), "
            f"LE: ({type(ret_le[0]).__name__}, {len(ret_le[1])}-tuple)",
        )

        # Check dret length
        se_dret_len = len(ret_se[1])
        le_dret_len = len(ret_le[1])
        record(
            "P5/vis_limit_mag/dret_len",
            se_dret_len == le_dret_len,
            f"SE dret has {se_dret_len} elements, LE has {le_dret_len}",
        )

    except Exception as e:
        record_error("P5/vis_limit_mag", str(e))


# ============================================================================
# PART 6: Specific Edge Cases
# ============================================================================


def test_part6_edge_cases():
    """Test edge cases and known problematic scenarios."""
    print("\n" + "=" * 70)
    print("PART 6: Edge Cases")
    print("=" * 70)

    jd_start = swe.julday(2024, 1, 1, 0.0)

    # 6a: Test MORNING_LAST constant value (existing test has bug: MORNING_LAST=6)
    print("\n  --- 6a: MORNING_LAST constant ---")
    record(
        "P6/MORNING_LAST_value",
        swe.MORNING_LAST == 4,
        f"pyswisseph MORNING_LAST = {swe.MORNING_LAST} (should be 4)",
    )
    record(
        "P6/LE_MORNING_LAST_value",
        ephem.MORNING_LAST == 4,
        f"libephemeris MORNING_LAST = {ephem.MORNING_LAST}",
    )

    # 6b: Test high latitude (Tromsø 69.6°N) — may fail for polar regions
    print("\n  --- 6b: High latitude ---")
    geopos_tromso = (19.0, 69.6, 0)
    try:
        ret_se = swe.heliacal_ut(
            jd_start,
            geopos_tromso,
            STANDARD_ATMO,
            STANDARD_OBSERVER,
            "Venus",
            HELIACAL_RISING,
            0,
        )
        ret_le = ephem.heliacal_ut(
            jd_start,
            geopos_tromso,
            STANDARD_ATMO,
            STANDARD_OBSERVER,
            "Venus",
            HELIACAL_RISING,
        )
        if ret_se[0] > 0 and ret_le[0] > 0:
            diff = abs(ret_se[0] - ret_le[0])
            record("P6/high_lat_Venus_rising", diff < 5.0, f"diff={diff:.3f}d")
        elif ret_se[0] <= 0 and ret_le[0] <= 0:
            record("P6/high_lat_Venus_rising", True, "Both returned 0 (no event)")
        else:
            record(
                "P6/high_lat_Venus_rising",
                False,
                f"SE={ret_se[0]:.5f} LE={ret_le[0]:.5f}",
            )
    except Exception as e:
        record_error("P6/high_lat_Venus_rising", str(e))

    # 6c: Test Sirius heliacal rising at Cairo (classic Egyptian observation)
    print("\n  --- 6c: Sirius at Cairo ---")
    geopos_cairo = (31.24, 30.04, 75)
    try:
        ret_se = swe.heliacal_ut(
            jd_start,
            geopos_cairo,
            STANDARD_ATMO,
            STANDARD_OBSERVER,
            "Sirius",
            HELIACAL_RISING,
            0,
        )
        ret_le = ephem.heliacal_ut(
            jd_start,
            geopos_cairo,
            STANDARD_ATMO,
            STANDARD_OBSERVER,
            "Sirius",
            HELIACAL_RISING,
        )
        if ret_se[0] > 0 and ret_le[0] > 0:
            diff = abs(ret_se[0] - ret_le[0])
            record(
                "P6/Sirius_Cairo_rising",
                diff < 3.0,
                f"SE={ret_se[0]:.5f} LE={ret_le[0]:.5f} diff={diff:.3f}d",
            )
        else:
            record("P6/Sirius_Cairo_rising", False, f"SE={ret_se[0]} LE={ret_le[0]}")
    except Exception as e:
        record_error("P6/Sirius_Cairo_rising", str(e))

    # 6d: Test NO_DETAILS flag
    print("\n  --- 6d: NO_DETAILS flag ---")
    geopos_rome = (12.5, 41.9, 0)
    try:
        # Use correct flag value from pyswisseph
        no_details_flag = swe.HELFLAG_NO_DETAILS

        ret_se = swe.heliacal_ut(
            jd_start,
            geopos_rome,
            STANDARD_ATMO,
            STANDARD_OBSERVER,
            "Venus",
            HELIACAL_RISING,
            no_details_flag,
        )
        ret_le = ephem.heliacal_ut(
            jd_start,
            geopos_rome,
            STANDARD_ATMO,
            STANDARD_OBSERVER,
            "Venus",
            HELIACAL_RISING,
            no_details_flag,
        )

        # With NO_DETAILS, jd2 and jd3 should be 0
        if ret_se[0] > 0:
            # SE behavior with NO_DETAILS
            se_jd2_zero = ret_se[1] == 0.0
            se_jd3_zero = ret_se[2] == 0.0
            le_jd2_zero = ret_le[1] == 0.0
            le_jd3_zero = ret_le[2] == 0.0

            record(
                "P6/NO_DETAILS/jd2_zero",
                le_jd2_zero == se_jd2_zero,
                f"SE jd2=0: {se_jd2_zero}, LE jd2=0: {le_jd2_zero}",
            )
            record(
                "P6/NO_DETAILS/jd3_zero",
                le_jd3_zero == se_jd3_zero,
                f"SE jd3=0: {se_jd3_zero}, LE jd3=0: {le_jd3_zero}",
            )
    except Exception as e:
        record_error("P6/NO_DETAILS", str(e))


# ============================================================================
# MAIN
# ============================================================================


def main():
    global total, passed, failed, skipped, errors, failures

    print("=" * 70)
    print("ROUND 7: Deep Heliacal Events Audit")
    print("=" * 70)
    t0 = time.time()

    test_part4_flag_constants()
    test_part5_api_shape()
    test_part1_heliacal_ut()
    test_part2_heliacal_pheno()
    test_part3_vis_limit_mag()
    test_part6_edge_cases()

    elapsed = time.time() - t0

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total:   {total}")
    print(f"Passed:  {passed}")
    print(f"Failed:  {failed}")
    print(f"Skipped: {skipped}")
    print(f"Errors:  {errors}")
    print(f"Time:    {elapsed:.1f}s")

    if failures:
        print(f"\n--- {len(failures)} FAILURES ---")
        for name, detail in failures:
            print(f"  {name}: {detail}")

    print(f"\nPass rate: {passed}/{total} = {100 * passed / max(total, 1):.1f}%")


if __name__ == "__main__":
    main()
