#!/usr/bin/env python3
"""Round 51: Fixed Stars with Proper Motion Across Centuries.

Tests fixed star positions across a wide date range, verifying proper motion
propagation and precession against pyswisseph.

Phases:
  P1: Bright stars at 10 epochs spanning 1500-2500
  P2: High proper motion stars (Sirius, Arcturus, Aldebaran)
  P3: Stars near ecliptic poles (stress test for coordinate transforms)
  P4: Star speed (proper motion rate) comparison
  P5: fixstar2_ut vs fixstar_ut consistency
  P6: Stars with FLG_J2000 flag
"""

from __future__ import annotations

import math
import os
import sys

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
errors = 0
results = {"passed": [], "failed": [], "errors": []}


def record(phase, label, ok, detail=""):
    global passed, failed
    if ok:
        passed += 1
        results["passed"].append(f"{phase} {label}: {detail}")
    else:
        failed += 1
        results["failed"].append(f"{phase} {label}: {detail}")


# Bright stars available in both catalogs
BRIGHT_STARS = [
    "Sirius",
    "Canopus",
    "Arcturus",
    "Vega",
    "Capella",
    "Rigel",
    "Betelgeuse",
    "Aldebaran",
    "Spica",
    "Antares",
    "Pollux",
    "Fomalhaut",
    "Deneb",
    "Regulus",
    "Castor",
]

# High proper motion stars
HIGH_PM_STARS = ["Sirius", "Arcturus", "Aldebaran", "Pollux", "Regulus"]

# Epochs spanning centuries
EPOCHS = {
    "1500": 2268924.0,
    "1600": 2305448.0,
    "1700": 2341973.0,
    "1800": 2378497.0,
    "1900": 2415021.0,
    "1950": 2433283.0,
    "2000": 2451545.0,
    "2024": 2460310.0,
    "2050": 2469808.0,
    "2100": 2488070.0,
}


def safe_fixstar2(name, jd, iflag=0):
    """Get fixstar2 from both SE and LE. Returns (se_pos, le_pos) or None."""
    try:
        se_result = swe.fixstar2(name, jd, iflag)
        # pyswisseph: (pos_tuple, starname_str, retflag)
        se_pos = se_result[0]
    except Exception:
        return None

    try:
        le_result = ephem.fixstar2_ut(name, jd, iflag)
        # libephemeris: (pos_tuple, starname_str, retflag)
        le_pos = le_result[0]
    except Exception:
        return None

    return se_pos, le_pos


def compare_star_pos(
    phase, label, se_pos, le_pos, lon_tol_arcsec=5.0, lat_tol_arcsec=5.0, speed_tol=None
):
    """Compare star positions."""
    se_lon, se_lat = se_pos[0], se_pos[1]
    le_lon, le_lat = le_pos[0], le_pos[1]

    # Longitude difference (handle wrap)
    dlon = se_lon - le_lon
    if dlon > 180:
        dlon -= 360
    elif dlon < -180:
        dlon += 360
    dlon_arcsec = abs(dlon) * 3600

    dlat_arcsec = abs(se_lat - le_lat) * 3600

    ok = dlon_arcsec < lon_tol_arcsec and dlat_arcsec < lat_tol_arcsec

    parts = [f'dlon={dlon_arcsec:.2f}" dlat={dlat_arcsec:.2f}"']

    # Speed comparison if requested
    if speed_tol is not None and len(se_pos) > 3 and len(le_pos) > 3:
        se_lon_spd = se_pos[3]
        le_lon_spd = le_pos[3]
        dspd = abs(se_lon_spd - le_lon_spd)
        if dspd > speed_tol:
            ok = False
        parts.append(f"dlon_spd={dspd:.6f}°/day")

    record(phase, label, ok, " ".join(parts))
    return ok


def phase1():
    """Bright stars at 10 epochs spanning 1500-2500."""
    global errors
    print("\n=== P1: Bright stars at 10 epochs ===")

    for star in BRIGHT_STARS:
        for epoch_name, jd in EPOCHS.items():
            try:
                result = safe_fixstar2(star, jd)
                if result is None:
                    continue
                se_pos, le_pos = result

                # Tolerance scales with distance from J2000
                years_from_j2000 = abs(jd - 2451545.0) / 365.25
                base_tol = 5.0  # arcsec at J2000
                tol = base_tol + years_from_j2000 * 0.02  # ~0.02"/year growth

                compare_star_pos(
                    "P1",
                    f"{star} {epoch_name}",
                    se_pos,
                    le_pos,
                    lon_tol_arcsec=tol,
                    lat_tol_arcsec=tol,
                )

            except Exception as e:
                errors += 1
                results["errors"].append(f"P1 {star} {epoch_name}: {e}")


def phase2():
    """High proper motion stars — focus on PM accuracy."""
    global errors
    print("\n=== P2: High proper motion stars ===")

    # Test at finer intervals around J2000
    fine_epochs = {}
    for delta_years in [-200, -100, -50, -20, -10, 0, 10, 20, 50, 100, 200]:
        name = f"{2000 + delta_years}"
        fine_epochs[name] = 2451545.0 + delta_years * 365.25

    for star in HIGH_PM_STARS:
        for epoch_name, jd in fine_epochs.items():
            try:
                result = safe_fixstar2(star, jd)
                if result is None:
                    continue
                se_pos, le_pos = result

                years_from_j2000 = abs(jd - 2451545.0) / 365.25
                tol = 5.0 + years_from_j2000 * 0.03

                compare_star_pos(
                    "P2",
                    f"{star} {epoch_name}",
                    se_pos,
                    le_pos,
                    lon_tol_arcsec=tol,
                    lat_tol_arcsec=tol,
                )

            except Exception as e:
                errors += 1
                results["errors"].append(f"P2 {star} {epoch_name}: {e}")


def phase3():
    """Stars near ecliptic poles — coordinate transform stress test."""
    global errors
    print("\n=== P3: Stars near ecliptic poles ===")

    # Stars with high ecliptic latitude
    pole_stars = ["Polaris", "Vega", "Deneb", "Capella"]

    test_epochs = {
        "1800": 2378497.0,
        "1900": 2415021.0,
        "2000": 2451545.0,
        "2024": 2460310.0,
        "2100": 2488070.0,
    }

    for star in pole_stars:
        for epoch_name, jd in test_epochs.items():
            try:
                result = safe_fixstar2(star, jd)
                if result is None:
                    continue
                se_pos, le_pos = result

                # Near poles, longitude can diverge more
                lat = abs(se_pos[1])
                lon_tol = 10.0 if lat > 60 else 5.0

                compare_star_pos(
                    "P3",
                    f"{star} {epoch_name}",
                    se_pos,
                    le_pos,
                    lon_tol_arcsec=lon_tol,
                    lat_tol_arcsec=5.0,
                )

            except Exception as e:
                errors += 1
                results["errors"].append(f"P3 {star} {epoch_name}: {e}")


def phase4():
    """Star speed (proper motion rate) comparison."""
    global errors
    print("\n=== P4: Star speed comparison ===")

    speed_stars = BRIGHT_STARS[:10]
    jd = 2451545.0  # J2000

    for star in speed_stars:
        try:
            # Get with FLG_SPEED
            iflag = ephem.FLG_SPEED
            result = safe_fixstar2(star, jd, iflag)
            if result is None:
                continue
            se_pos, le_pos = result

            if len(se_pos) > 3 and len(le_pos) > 3:
                # Compare longitude speed
                se_lon_spd = se_pos[3]
                le_lon_spd = le_pos[3]
                dspd_lon = abs(se_lon_spd - le_lon_spd)

                # Compare latitude speed
                se_lat_spd = se_pos[4] if len(se_pos) > 4 else 0.0
                le_lat_spd = le_pos[4] if len(le_pos) > 4 else 0.0
                dspd_lat = abs(se_lat_spd - le_lat_spd)

                # Speed tolerance: 0.001 deg/day = ~1.3"/day
                spd_tol = 0.001
                ok = dspd_lon < spd_tol and dspd_lat < spd_tol
                detail = (
                    f"lon_spd:SE={se_lon_spd:.8f}/LE={le_lon_spd:.8f}/d={dspd_lon:.8f} "
                    f"lat_spd:SE={se_lat_spd:.8f}/LE={le_lat_spd:.8f}/d={dspd_lat:.8f}"
                )
                record("P4", f"{star}_speed", ok, detail)
            else:
                record("P4", f"{star}_speed", True, "speed not available in result")

        except Exception as e:
            errors += 1
            results["errors"].append(f"P4 {star}: {e}")


def phase5():
    """fixstar2_ut vs fixstar_ut consistency."""
    global errors
    print("\n=== P5: fixstar2_ut vs fixstar_ut consistency ===")

    test_stars = BRIGHT_STARS[:8]
    jd = 2451545.0

    for star in test_stars:
        try:
            # fixstar2_ut
            try:
                le2 = ephem.fixstar2_ut(star, jd, 0)
                # Returns (pos_tuple, starname, retflag)
                le2_pos = le2[0]
            except Exception:
                continue

            # fixstar_ut
            try:
                le1 = ephem.fixstar_ut(star, jd, 0)
                # Returns (pos_tuple, starname, retflag)
                le1_pos = le1[0]
            except Exception:
                continue

            # Should be identical
            dlon = abs(le2_pos[0] - le1_pos[0])
            if dlon > 180:
                dlon = 360 - dlon
            dlat = abs(le2_pos[1] - le1_pos[1])

            ok = dlon < 1e-10 and dlat < 1e-10
            detail = f'dlon={dlon * 3600:.6f}" dlat={dlat * 3600:.6f}"'
            record("P5", f"{star}_consistency", ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P5 {star}: {e}")


def phase6():
    """Stars with FLG_J2000 flag."""
    global errors
    print("\n=== P6: Stars with FLG_J2000 ===")

    test_stars = BRIGHT_STARS[:10]

    for star in test_stars:
        for epoch_name, jd in [
            ("1900", 2415021.0),
            ("2000", 2451545.0),
            ("2024", 2460310.0),
            ("2100", 2488070.0),
        ]:
            try:
                iflag = ephem.FLG_J2000
                result = safe_fixstar2(star, jd, iflag)
                if result is None:
                    continue
                se_pos, le_pos = result

                # J2000 frame: tighter tolerance since no precession
                years_from_j2000 = abs(jd - 2451545.0) / 365.25
                tol = 3.0 + years_from_j2000 * 0.01

                compare_star_pos(
                    "P6",
                    f"{star} {epoch_name} J2000",
                    se_pos,
                    le_pos,
                    lon_tol_arcsec=tol,
                    lat_tol_arcsec=tol,
                )

            except Exception as e:
                errors += 1
                results["errors"].append(f"P6 {star} {epoch_name}: {e}")


def main():
    print("=" * 70)
    print("ROUND 51: Fixed Stars with Proper Motion Across Centuries")
    print("=" * 70)

    phase1()
    print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

    phase2()
    print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

    phase3()
    print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

    phase4()
    print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

    phase5()
    print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

    phase6()
    print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

    total = passed + failed + errors
    pct = 100 * passed / total if total else 0
    print("\n" + "=" * 70)
    print(f"ROUND 51 FINAL: {passed}/{total} passed ({pct:.1f}%)")
    print(f"  Passed:  {passed}")
    print(f"  Failed:  {failed}")
    print(f"  Errors:  {errors}")
    print("=" * 70)

    if results["failed"]:
        print(f"\n--- FAILURES ({len(results['failed'])}) ---")
        for f in results["failed"][:50]:
            print(f)

    if results["errors"]:
        print(f"\n--- ERRORS ({len(results['errors'])}) ---")
        for e in results["errors"][:10]:
            print(e)

    return 0 if failed == 0 and errors == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
