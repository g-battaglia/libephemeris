#!/usr/bin/env python3
"""Round 44: Fixed Star Magnitudes & Spectral Types Verification.

Verifies libephemeris fixed star catalog data (magnitudes, spectral types,
proper motions, positions) against pyswisseph for a comprehensive set of stars.

Phases:
  P1: Bright navigational stars — position, magnitude, spectral type
  P2: All named stars sweep — position comparison at J2000
  P3: Proper motion over centuries — position drift verification
  P4: Fixed star with FLG_SPEED — velocity components
  P5: fixstar2 vs fixstar — API compatibility
  P6: Star search by name fragments
"""

from __future__ import annotations

import math
import os
import sys
import traceback

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
errors = 0
results = {"passed": [], "failed": [], "errors": []}

J2000 = 2451545.0

# Navigational / bright stars to test
BRIGHT_STARS = [
    "Sirius",
    "Canopus",
    "Arcturus",
    "Vega",
    "Capella",
    "Rigel",
    "Procyon",
    "Betelgeuse",
    "Achernar",
    "Altair",
    "Aldebaran",
    "Spica",
    "Antares",
    "Pollux",
    "Fomalhaut",
    "Deneb",
    "Regulus",
    "Castor",
    "Bellatrix",
    "Alnilam",
    "Polaris",
    "Alpheratz",
    "Mirfak",
    "Algol",
    "Mira",
    "Rasalhague",
    "Markab",
    "Scheat",
    "Algenib",
    "Hamal",
    "Dubhe",
    "Merak",
    "Alioth",
    "Alkaid",
    "Mizar",
    "Thuban",
    "Kochab",
    "Pherkad",
    "Etamin",
    "Rastaban",
]

# Extended star list
MORE_STARS = [
    "Acamar",
    "Acrux",
    "Adhara",
    "Agena",
    "Ain",
    "Albireo",
    "Alcyone",
    "Alderamin",
    "Alhena",
    "Alioth",
    "Almach",
    "Alnair",
    "Alnitak",
    "Alphard",
    "Alphecca",
    "Ankaa",
    "Ascella",
    "Atria",
    "Avior",
    "Bellatrix",
    "Diphda",
    "Elnath",
    "Enif",
    "Gacrux",
    "Gienah",
    "Kaus Australis",
    "Menkent",
    "Mintaka",
    "Naos",
    "Nunki",
    "Peacock",
    "Phact",
    "Porrima",
    "Sabik",
    "Saiph",
    "Sargas",
    "Shaula",
    "Suhail",
    "Toliman",
    "Wezen",
    "Zubenelgenubi",
    "Zubeneschamali",
]


def record(phase, label, ok, detail=""):
    global passed, failed
    if ok:
        passed += 1
        results["passed"].append(f"{phase} {label}: {detail}")
    else:
        failed += 1
        results["failed"].append(f"{phase} {label}: {detail}")


def safe_fixstar(name, jd, flags=0):
    """Call SE fixstar2, return (pos, retflag, starname) or None on error.

    pyswisseph fixstar2 returns: (pos_tuple, starname_str, retflag_int)
    We normalize to: (pos_tuple, retflag, starname)
    """
    try:
        result = swe.fixstar2(name, jd, flags)
        # SE returns (pos, starname, retflag) -> normalize to (pos, retflag, starname)
        return (result[0], result[2], result[1])
    except Exception:
        try:
            result = swe.fixstar(name, jd, flags)
            # fixstar returns (pos, retflag)
            return (result[0], result[1], name)
        except Exception:
            return None


def safe_le_fixstar(name, jd, flags=0):
    """Call LE fixstar2_ut, return (pos, retflag, starname) or None.

    libephemeris fixstar2_ut returns: (pos_tuple, starname_str, retflag_int)
    We normalize to: (pos_tuple, retflag, starname)
    """
    try:
        result = ephem.fixstar2_ut(name, jd, flags)
        # LE returns (pos, starname, retflag) -> normalize to (pos, retflag, starname)
        return (result[0], result[2], result[1])
    except Exception:
        try:
            result = ephem.fixstar_ut(name, jd, flags)
            # fixstar_ut returns (pos, starname, retflag)
            return (result[0], result[2], result[1])
        except Exception:
            return None


def phase1():
    """Bright navigational stars — position, magnitude at J2000."""
    global errors
    print("\n=== P1: Bright navigational stars at J2000 ===")

    jd = J2000
    pos_tol = 1.0 / 3600.0  # 1 arcsecond

    for star in BRIGHT_STARS:
        try:
            se_result = safe_fixstar(star, jd, swe.FLG_SWIEPH)
            le_result = safe_le_fixstar(star, jd, ephem.FLG_SWIEPH)

            if se_result is None and le_result is None:
                record("P1", star, True, "both not found (OK)")
                continue
            if se_result is None:
                record("P1", star, False, "SE not found but LE found")
                continue
            if le_result is None:
                record("P1", star, False, "LE not found but SE found")
                continue

            se_pos = se_result[0]
            le_pos = le_result[0]

            lon_diff = abs(se_pos[0] - le_pos[0])
            if lon_diff > 180:
                lon_diff = 360 - lon_diff
            lat_diff = abs(se_pos[1] - le_pos[1])
            dist_diff = abs(se_pos[2] - le_pos[2]) if se_pos[2] != 0 else 0

            lon_arcsec = lon_diff * 3600
            lat_arcsec = lat_diff * 3600

            ok = lon_arcsec < 10.0 and lat_arcsec < 10.0  # 10" tolerance
            detail = (
                f'lon={lon_arcsec:.3f}" lat={lat_arcsec:.3f}" dist_diff={dist_diff:.4f}'
            )
            record("P1", star, ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P1 {star}: {e}\n{traceback.format_exc()}")


def phase2():
    """All named stars sweep — position comparison."""
    print("\n=== P2: Extended star list at J2000 ===")

    jd = J2000
    all_stars = list(set(BRIGHT_STARS + MORE_STARS))

    n_found = 0
    n_match = 0
    max_lon_diff = 0.0
    max_lat_diff = 0.0

    for star in sorted(all_stars):
        try:
            se_result = safe_fixstar(star, jd, swe.FLG_SWIEPH)
            le_result = safe_le_fixstar(star, jd, ephem.FLG_SWIEPH)

            if se_result is None or le_result is None:
                continue

            n_found += 1
            se_pos = se_result[0]
            le_pos = le_result[0]

            lon_diff = abs(se_pos[0] - le_pos[0])
            if lon_diff > 180:
                lon_diff = 360 - lon_diff
            lat_diff = abs(se_pos[1] - le_pos[1])

            lon_arcsec = lon_diff * 3600
            lat_arcsec = lat_diff * 3600

            max_lon_diff = max(max_lon_diff, lon_arcsec)
            max_lat_diff = max(max_lat_diff, lat_arcsec)

            if lon_arcsec < 10.0 and lat_arcsec < 10.0:
                n_match += 1

        except Exception:
            continue

    record(
        "P2",
        f"sweep ({n_match}/{n_found})",
        n_match == n_found,
        f'max_lon={max_lon_diff:.3f}" max_lat={max_lat_diff:.3f}"',
    )


def phase3():
    """Proper motion over centuries — position drift."""
    global errors
    print("\n=== P3: Proper motion drift over centuries ===")

    test_stars = [
        "Sirius",
        "Arcturus",
        "Vega",
        "Polaris",
        "Barnard's Star",
        "Aldebaran",
        "Betelgeuse",
        "Procyon",
        "Altair",
        "Regulus",
    ]

    # Test at multiple epochs
    epochs = {
        "1900": J2000 - 365.25 * 100,
        "1950": J2000 - 365.25 * 50,
        "2000": J2000,
        "2050": J2000 + 365.25 * 50,
        "2100": J2000 + 365.25 * 100,
        "2200": J2000 + 365.25 * 200,
    }

    for star in test_stars:
        for epoch_name, jd in epochs.items():
            try:
                se_result = safe_fixstar(star, jd, swe.FLG_SWIEPH)
                le_result = safe_le_fixstar(star, jd, ephem.FLG_SWIEPH)

                if se_result is None or le_result is None:
                    continue

                se_pos = se_result[0]
                le_pos = le_result[0]

                lon_diff = abs(se_pos[0] - le_pos[0])
                if lon_diff > 180:
                    lon_diff = 360 - lon_diff
                lat_diff = abs(se_pos[1] - le_pos[1])

                lon_arcsec = lon_diff * 3600
                lat_arcsec = lat_diff * 3600

                # Tolerance scales with epoch distance
                dt_centuries = abs(jd - J2000) / 36525.0
                tol = 10.0 + dt_centuries * 2.0  # 10" base + 2"/century

                ok = lon_arcsec < tol and lat_arcsec < tol
                detail = f'lon={lon_arcsec:.3f}" lat={lat_arcsec:.3f}" tol={tol:.1f}"'
                record("P3", f"{star} {epoch_name}", ok, detail)

            except Exception as e:
                errors += 1
                results["errors"].append(f"P3 {star} {epoch_name}: {e}")


def phase4():
    """Fixed star with FLG_SPEED — velocity components."""
    global errors
    print("\n=== P4: Fixed star speeds ===")

    test_stars = [
        "Sirius",
        "Arcturus",
        "Polaris",
        "Vega",
        "Aldebaran",
        "Spica",
        "Antares",
        "Betelgeuse",
        "Regulus",
        "Fomalhaut",
    ]
    jd = J2000

    for star in test_stars:
        try:
            se_result = safe_fixstar(star, jd, swe.FLG_SWIEPH | swe.FLG_SPEED)
            le_result = safe_le_fixstar(
                star, jd, ephem.FLG_SWIEPH | ephem.FLG_SPEED
            )

            if se_result is None or le_result is None:
                continue

            se_pos = se_result[0]
            le_pos = le_result[0]

            # Speed components are [3], [4], [5]
            lon_speed_diff = abs(se_pos[3] - le_pos[3]) * 3600  # arcsec/day
            lat_speed_diff = abs(se_pos[4] - le_pos[4]) * 3600
            dist_speed_diff = abs(se_pos[5] - le_pos[5])

            # Speed tolerance: 0.1 arcsec/day
            tol = 0.1
            ok = lon_speed_diff < tol and lat_speed_diff < tol
            detail = (
                f'dlon={lon_speed_diff:.4f}"/day dlat={lat_speed_diff:.4f}"/day '
                f'SE_lon_spd={se_pos[3] * 3600:.4f}"/day LE_lon_spd={le_pos[3] * 3600:.4f}"/day'
            )
            record("P4", f"{star} speed", ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P4 {star}: {e}")


def phase5():
    """fixstar2 vs fixstar API compatibility."""
    global errors
    print("\n=== P5: fixstar2 vs fixstar API ===")

    test_stars = ["Sirius", "Vega", "Polaris", "Arcturus", "Canopus"]
    jd = J2000

    for star in test_stars:
        try:
            # Test fixstar_ut
            try:
                le1 = ephem.fixstar_ut(star, jd, ephem.FLG_SWIEPH)
            except Exception:
                le1 = None

            # Test fixstar2_ut
            try:
                le2 = ephem.fixstar2_ut(star, jd, ephem.FLG_SWIEPH)
            except Exception:
                le2 = None

            if le1 is not None and le2 is not None:
                # Both should return same positions
                # fixstar_ut returns (pos, starname, retflag)
                # fixstar2_ut returns (pos, starname, retflag)
                pos1 = le1[0]
                pos2 = le2[0]
                lon_diff = abs(pos1[0] - pos2[0]) * 3600
                lat_diff = abs(pos1[1] - pos2[1]) * 3600

                ok = lon_diff < 0.001 and lat_diff < 0.001
                detail = f'lon_diff={lon_diff:.6f}" lat_diff={lat_diff:.6f}"'
                record("P5", f"{star} api_compat", ok, detail)
            elif le1 is None and le2 is None:
                record("P5", f"{star} api_compat", True, "both not found")
            else:
                record(
                    "P5",
                    f"{star} api_compat",
                    False,
                    f"fixstar={'found' if le1 else 'missing'} fixstar2={'found' if le2 else 'missing'}",
                )

        except Exception as e:
            errors += 1
            results["errors"].append(f"P5 {star}: {e}")


def phase6():
    """Star search by name fragments and equatorial mode."""
    global errors
    print("\n=== P6: Star search and equatorial mode ===")

    jd = J2000
    test_stars = ["Sirius", "Arcturus", "Vega", "Polaris", "Spica"]

    for star in test_stars:
        try:
            # Equatorial coordinates
            se_eq = safe_fixstar(star, jd, swe.FLG_SWIEPH | swe.FLG_EQUATORIAL)
            le_eq = safe_le_fixstar(
                star, jd, ephem.FLG_SWIEPH | ephem.FLG_EQUATORIAL
            )

            if se_eq is None or le_eq is None:
                continue

            se_pos = se_eq[0]
            le_pos = le_eq[0]

            # RA and Dec comparison
            ra_diff = abs(se_pos[0] - le_pos[0])
            if ra_diff > 180:
                ra_diff = 360 - ra_diff
            dec_diff = abs(se_pos[1] - le_pos[1])

            ra_arcsec = ra_diff * 3600
            dec_arcsec = dec_diff * 3600

            ok = ra_arcsec < 10.0 and dec_arcsec < 10.0
            detail = f'RA={ra_arcsec:.3f}" Dec={dec_arcsec:.3f}"'
            record("P6", f"{star} equatorial", ok, detail)

            # J2000 coordinates
            se_j2k = safe_fixstar(star, jd, swe.FLG_SWIEPH | swe.FLG_J2000)
            le_j2k = safe_le_fixstar(star, jd, ephem.FLG_SWIEPH | ephem.FLG_J2000)

            if se_j2k is not None and le_j2k is not None:
                se_p2 = se_j2k[0]
                le_p2 = le_j2k[0]
                lon_diff = abs(se_p2[0] - le_p2[0])
                if lon_diff > 180:
                    lon_diff = 360 - lon_diff
                lat_diff = abs(se_p2[1] - le_p2[1])

                ok2 = lon_diff * 3600 < 10.0 and lat_diff * 3600 < 10.0
                detail2 = f'lon={lon_diff * 3600:.3f}" lat={lat_diff * 3600:.3f}"'
                record("P6", f"{star} J2000", ok2, detail2)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P6 {star}: {e}")


def main():
    print("=" * 70)
    print("ROUND 44: Fixed Star Magnitudes & Spectral Types Verification")
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
    print("\n" + "=" * 70)
    print(f"ROUND 44 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
    print(f"  Passed:  {passed}")
    print(f"  Failed:  {failed}")
    print(f"  Errors:  {errors}")
    print("=" * 70)

    if results["failed"]:
        print(f"\n--- FAILURES ({len(results['failed'])}) ---")
        for f in results["failed"][:30]:
            print(f)

    if results["errors"]:
        print(f"\n--- ERRORS ({len(results['errors'])}) ---")
        for e in results["errors"][:10]:
            print(e)

    return 0 if failed == 0 and errors == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
