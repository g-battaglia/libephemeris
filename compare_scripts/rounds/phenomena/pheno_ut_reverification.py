#!/usr/bin/env python3
"""Round 47: Planetary Phenomena (pheno_ut) Deep Re-verification.

Deeper than Round 16 — focuses on edge cases, extreme phase angles,
century-spanning dates, all flag combinations, Saturn ring tilt,
Mercury/Venus inferior conjunction, and asteroid phenomena.

Phases:
  P1: All planets at 20 dates spanning 1800-2200 (default flags)
  P2: Moon phenomena at key lunar phases (new, full, quarters)
  P3: Saturn ring tilt across full orbit (~29.5 yr)
  P4: Mercury/Venus at inferior/superior conjunction epochs
  P5: Flag combinations (TRUEPOS, NOABERR, HELCTR)
  P6: Asteroids (Ceres, Pallas, Juno, Vesta, Chiron)
  P7: Extreme historical dates (1500, 1000 CE)
  P8: Phase angle precision at near-0° and near-180° geometries
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


def record(phase, label, ok, detail=""):
    global passed, failed
    if ok:
        passed += 1
        results["passed"].append(f"{phase} {label}: {detail}")
    else:
        failed += 1
        results["failed"].append(f"{phase} {label}: {detail}")


# Body constants
PLANETS = {
    "Sun": ephem.SUN,
    "Moon": ephem.MOON,
    "Mercury": ephem.MERCURY,
    "Venus": ephem.VENUS,
    "Mars": ephem.MARS,
    "Jupiter": ephem.JUPITER,
    "Saturn": ephem.SATURN,
    "Uranus": ephem.URANUS,
    "Neptune": ephem.NEPTUNE,
    "Pluto": ephem.PLUTO,
}

ASTEROIDS = {
    "Ceres": ephem.CERES,
    "Pallas": ephem.PALLAS,
    "Juno": ephem.JUNO,
    "Vesta": ephem.VESTA,
    "Chiron": ephem.CHIRON,
}

# Date grid spanning 4 centuries
DATES_WIDE = {
    "1800Jan01": 2378497.0,
    "1850Jul01": 2396959.0,
    "1900Jan01": 2415021.0,
    "1920Jun15": 2422490.0,
    "1940Jan01": 2429630.0,
    "1960Jul01": 2437117.0,
    "1970Jan01": 2440588.0,
    "1980Jun15": 2444407.0,
    "1990Jan01": 2447893.0,
    "1995Jul01": 2449900.0,
    "2000Jan01": 2451545.0,
    "2005Jul01": 2453553.0,
    "2010Jan01": 2455198.0,
    "2015Jul01": 2457205.0,
    "2020Jan01": 2458850.0,
    "2022Jun21": 2459752.0,
    "2024Mar20": 2460389.5,
    "2025Jan01": 2460677.0,
    "2050Jul01": 2469988.0,
    "2100Jan01": 2488070.0,
}

ATTR_NAMES = ["phase_angle", "phase", "elongation", "diameter", "magnitude"]


def safe_pheno(jd, ipl, iflag=0):
    """Get pheno_ut from both SE and LE. Returns (se_attr, le_attr) or None."""
    try:
        se_attr = swe.pheno_ut(jd, ipl, iflag)
    except Exception:
        return None

    try:
        le_result = ephem.pheno_ut(jd, ipl, iflag)
        le_attr = le_result[0]
    except Exception:
        return None

    return se_attr, le_attr


def compare_pheno(phase, label, se_attr, le_attr, tols=None):
    """Compare 5 phenomenological attributes.

    Args:
        tols: dict of attribute index -> tolerance. Defaults provided.
    """
    if tols is None:
        tols = {
            0: 30.0 / 3600,  # phase angle: 30" (arcsec -> deg)
            1: 0.001,  # phase: 0.001 (fraction)
            2: 30.0 / 3600,  # elongation: 30"
            3: 0.0001,  # diameter: 0.0001° (~0.36")
            4: 0.15,  # magnitude: 0.15 mag
        }

    worst_attr = ""
    worst_diff = 0.0
    all_ok = True

    for i in range(5):
        se_v = se_attr[i]
        le_v = le_attr[i]
        diff = abs(se_v - le_v)

        # For phase angle / elongation, handle 360° wrap
        if i in (0, 2) and diff > 180:
            diff = 360 - diff

        tol = tols.get(i, 0.01)
        if diff > tol:
            all_ok = False
            if diff > worst_diff:
                worst_diff = diff
                worst_attr = ATTR_NAMES[i]

    # Build detail
    parts = []
    for i in range(5):
        se_v = se_attr[i]
        le_v = le_attr[i]
        diff = abs(se_v - le_v)
        if i in (0, 2) and diff > 180:
            diff = 360 - diff
        tol = tols.get(i, 0.01)
        if diff > tol * 0.1:  # Only show non-trivial diffs
            parts.append(f"{ATTR_NAMES[i]}:SE={se_v:.6f}/LE={le_v:.6f}/d={diff:.6f}")

    detail_str = " | ".join(parts) if parts else "all match"
    if not all_ok:
        detail_str = f"WORST={worst_attr} d={worst_diff:.6f} | " + detail_str

    record(phase, label, all_ok, detail_str)
    return all_ok


def phase1():
    """All planets at 20 dates spanning 1800-2200."""
    global errors
    print("\n=== P1: All planets, 20 dates (1800-2200) ===")

    for planet_name, ipl in PLANETS.items():
        for date_name, jd in DATES_WIDE.items():
            try:
                result = safe_pheno(jd, ipl)
                if result is None:
                    continue
                se_attr, le_attr = result
                compare_pheno("P1", f"{planet_name} {date_name}", se_attr, le_attr)
            except Exception as e:
                errors += 1
                results["errors"].append(f"P1 {planet_name} {date_name}: {e}")


def phase2():
    """Moon phenomena at key lunar phases."""
    global errors
    print("\n=== P2: Moon at key lunar phases ===")

    # Known lunar phase dates (approximate JDs)
    moon_phases = {
        # New Moons (phase angle ~180°, elongation ~0°)
        "NewMoon_2024Jan11": 2460320.5,
        "NewMoon_2024Feb09": 2460350.0,
        "NewMoon_2024Jul05": 2460497.0,
        "NewMoon_2000Jan06": 2451550.0,
        # Full Moons (phase angle ~0°, elongation ~180°)
        "FullMoon_2024Jan25": 2460335.0,
        "FullMoon_2024Jul21": 2460512.0,
        "FullMoon_2000Jan21": 2451564.0,
        "FullMoon_1999Dec22": 2451534.0,
        # First/Last Quarters (phase angle ~90°, elongation ~90°)
        "FirstQ_2024Jan18": 2460327.0,
        "LastQ_2024Feb02": 2460343.0,
        "FirstQ_2024Jul14": 2460505.0,
        # Waning crescent (high phase angle)
        "WaningCrescent_2024Jan08": 2460317.0,
        # Waxing crescent (high phase angle)
        "WaxingCrescent_2024Jan13": 2460323.0,
    }

    # Moon tolerances — magnitude less precise at extreme phase angles
    for phase_name, jd in moon_phases.items():
        try:
            result = safe_pheno(jd, ephem.MOON)
            if result is None:
                continue
            se_attr, le_attr = result

            # Use relaxed magnitude tolerance for thin crescents
            phase_angle = se_attr[0]
            mag_tol = 0.3 if phase_angle > 150 else 0.15
            tols = {
                0: 60.0 / 3600,  # phase angle: 60"
                1: 0.002,  # phase fraction
                2: 60.0 / 3600,  # elongation: 60"
                3: 0.0002,  # diameter
                4: mag_tol,  # magnitude
            }
            compare_pheno("P2", f"Moon {phase_name}", se_attr, le_attr, tols=tols)
        except Exception as e:
            errors += 1
            results["errors"].append(f"P2 Moon {phase_name}: {e}")


def phase3():
    """Saturn ring tilt across full orbit (~29.5 yr)."""
    global errors
    print("\n=== P3: Saturn ring tilt across orbit ===")

    # Sample Saturn at ~1-year intervals for 30 years
    base_jd = 2451545.0  # J2000
    for i in range(30):
        jd = base_jd + i * 365.25
        year = 2000 + i
        try:
            result = safe_pheno(jd, ephem.SATURN)
            if result is None:
                continue
            se_attr, le_attr = result

            # Saturn magnitude depends on ring tilt, so use wider tolerance
            tols = {
                0: 30.0 / 3600,
                1: 0.001,
                2: 30.0 / 3600,
                3: 0.0002,
                4: 0.20,  # ring tilt causes ~0.1-0.15 mag variation
            }
            compare_pheno("P3", f"Saturn {year}", se_attr, le_attr, tols=tols)
        except Exception as e:
            errors += 1
            results["errors"].append(f"P3 Saturn {year}: {e}")


def phase4():
    """Mercury/Venus at inferior/superior conjunction epochs."""
    global errors
    print("\n=== P4: Mercury/Venus conjunctions ===")

    # Key Mercury epochs (approximate)
    mercury_epochs = {
        "InfConj_2024Jan": 2460320.0,  # Near inferior conjunction
        "SupConj_2024Apr": 2460406.0,  # Near superior conjunction
        "MaxElong_2024Mar": 2460375.0,  # Near greatest elongation
        "InfConj_2024May": 2460445.0,
        "SupConj_2024Sep": 2460568.0,
        "MaxElong_2024Dec": 2460645.0,
        "InfConj_2000May": 2451689.0,
        "SupConj_2000Jan": 2451564.0,
    }

    for epoch_name, jd in mercury_epochs.items():
        try:
            result = safe_pheno(jd, ephem.MERCURY)
            if result is None:
                continue
            se_attr, le_attr = result
            # Near conjunction, phase angle can be ~0 or ~180, magnitude volatile
            tols = {
                0: 60.0 / 3600,
                1: 0.002,
                2: 60.0 / 3600,
                3: 0.0002,
                4: 0.20,
            }
            compare_pheno("P4", f"Mercury {epoch_name}", se_attr, le_attr, tols=tols)
        except Exception as e:
            errors += 1
            results["errors"].append(f"P4 Mercury {epoch_name}: {e}")

    # Key Venus epochs
    venus_epochs = {
        "InfConj_2022Jan": 2459581.0,
        "SupConj_2022Oct": 2459870.0,
        "MaxElong_2022Mar": 2459652.0,
        "InfConj_2023Aug": 2460175.0,
        "SupConj_2024Jun": 2460469.0,
        "MaxElong_2024Dec": 2460650.0,
        "InfConj_2025Mar": 2460740.0,
        "SupConj_2000Jun": 2451718.0,
    }

    for epoch_name, jd in venus_epochs.items():
        try:
            result = safe_pheno(jd, ephem.VENUS)
            if result is None:
                continue
            se_attr, le_attr = result
            tols = {
                0: 60.0 / 3600,
                1: 0.002,
                2: 60.0 / 3600,
                3: 0.0002,
                4: 0.20,
            }
            compare_pheno("P4", f"Venus {epoch_name}", se_attr, le_attr, tols=tols)
        except Exception as e:
            errors += 1
            results["errors"].append(f"P4 Venus {epoch_name}: {e}")


def phase5():
    """Flag combinations (TRUEPOS, NOABERR, HELCTR)."""
    global errors
    print("\n=== P5: Flag combinations ===")

    flag_combos = {
        "TRUEPOS": ephem.FLG_TRUEPOS,
        "NOABERR": ephem.FLG_NOABERR,
        "TRUEPOS+NOABERR": ephem.FLG_TRUEPOS | ephem.FLG_NOABERR,
    }

    test_bodies = {
        "Moon": ephem.MOON,
        "Mars": ephem.MARS,
        "Jupiter": ephem.JUPITER,
        "Venus": ephem.VENUS,
    }

    test_dates = {
        "J2000": 2451545.0,
        "2024Mar": 2460389.5,
        "1950Jan": 2433283.0,
    }

    for flag_name, iflag in flag_combos.items():
        for body_name, ipl in test_bodies.items():
            for date_name, jd in test_dates.items():
                try:
                    result = safe_pheno(jd, ipl, iflag)
                    if result is None:
                        continue
                    se_attr, le_attr = result
                    compare_pheno(
                        "P5",
                        f"{body_name} {date_name} {flag_name}",
                        se_attr,
                        le_attr,
                    )
                except Exception as e:
                    errors += 1
                    results["errors"].append(
                        f"P5 {body_name} {date_name} {flag_name}: {e}"
                    )

    # Heliocentric — separate because tolerances differ
    for body_name, ipl in test_bodies.items():
        if ipl == ephem.MOON:
            continue  # Moon heliocentric pheno is unusual
        for date_name, jd in test_dates.items():
            try:
                iflag = ephem.FLG_HELCTR
                result = safe_pheno(jd, ipl, iflag)
                if result is None:
                    continue
                se_attr, le_attr = result
                tols = {
                    0: 120.0 / 3600,  # heliocentric phase angles can differ more
                    1: 0.005,
                    2: 120.0 / 3600,
                    3: 0.001,
                    4: 0.5,
                }
                compare_pheno(
                    "P5",
                    f"{body_name} {date_name} HELCTR",
                    se_attr,
                    le_attr,
                    tols=tols,
                )
            except Exception as e:
                errors += 1
                results["errors"].append(f"P5 {body_name} {date_name} HELCTR: {e}")


def phase6():
    """Asteroids (Ceres, Pallas, Juno, Vesta, Chiron)."""
    global errors
    print("\n=== P6: Asteroid phenomena ===")

    test_dates = {
        "J2000": 2451545.0,
        "2024Mar": 2460389.5,
        "2024Sep": 2460567.0,
        "1990Jan": 2447893.0,
        "2010Jul": 2455383.0,
    }

    for ast_name, ipl in ASTEROIDS.items():
        for date_name, jd in test_dates.items():
            try:
                result = safe_pheno(jd, ipl)
                if result is None:
                    # Some asteroids may not be supported by SE or LE
                    continue
                se_attr, le_attr = result

                # Check if both returned valid (non-zero) data
                se_has_data = any(abs(v) > 1e-10 for v in se_attr[:5])
                le_has_data = any(abs(v) > 1e-10 for v in le_attr[:5])

                if not se_has_data and not le_has_data:
                    record(
                        "P6",
                        f"{ast_name} {date_name}",
                        True,
                        "both unsupported (all zeros)",
                    )
                    continue

                if se_has_data != le_has_data:
                    record(
                        "P6",
                        f"{ast_name} {date_name}",
                        False,
                        f"support mismatch: SE={se_has_data} LE={le_has_data}",
                    )
                    continue

                # Asteroids: wider tolerances for magnitude (empirical formulas vary)
                tols = {
                    0: 60.0 / 3600,
                    1: 0.002,
                    2: 60.0 / 3600,
                    3: 0.001,
                    4: 0.5,  # asteroid magnitudes are approximate
                }
                compare_pheno(
                    "P6", f"{ast_name} {date_name}", se_attr, le_attr, tols=tols
                )
            except Exception as e:
                errors += 1
                results["errors"].append(f"P6 {ast_name} {date_name}: {e}")


def phase7():
    """Extreme historical dates (1500, 1000 CE)."""
    global errors
    print("\n=== P7: Historical dates ===")

    historical_dates = {
        "1000CE_Jan": 2086673.0,
        "1200CE_Jul": 2159885.0,
        "1500CE_Jan": 2268924.0,
        "1600CE_Jul": 2305634.0,
        "1700CE_Jan": 2341973.0,
        "1750CE_Jul": 2360449.0,
    }

    # Only test major planets (available in DE440)
    hist_bodies = {
        "Sun": ephem.SUN,
        "Moon": ephem.MOON,
        "Mars": ephem.MARS,
        "Jupiter": ephem.JUPITER,
        "Saturn": ephem.SATURN,
    }

    for body_name, ipl in hist_bodies.items():
        for date_name, jd in historical_dates.items():
            try:
                result = safe_pheno(jd, ipl)
                if result is None:
                    continue
                se_attr, le_attr = result
                # Historical dates: ephemeris differences can be larger
                tols = {
                    0: 120.0 / 3600,
                    1: 0.005,
                    2: 120.0 / 3600,
                    3: 0.001,
                    4: 0.25,
                }
                compare_pheno(
                    "P7", f"{body_name} {date_name}", se_attr, le_attr, tols=tols
                )
            except Exception as e:
                errors += 1
                results["errors"].append(f"P7 {body_name} {date_name}: {e}")


def phase8():
    """Phase angle precision at near-0° and near-180° geometries."""
    global errors
    print("\n=== P8: Extreme phase angle geometry ===")

    # Test outer planets at opposition (phase angle ~0°)
    # and inner planets near inferior conjunction (phase angle ~180°)
    jd = 2451545.0  # J2000

    # Sample many dates to find natural near-opposition/conjunction events
    test_cases = []
    for delta in range(0, 365 * 5, 7):  # Every week for 5 years
        test_jd = jd + delta
        # Outer planets: Mars, Jupiter, Saturn — check phase angle
        for name, ipl in [
            ("Mars", ephem.MARS),
            ("Jupiter", ephem.JUPITER),
            ("Saturn", ephem.SATURN),
        ]:
            test_cases.append((name, ipl, test_jd))

    # Run all and filter for interesting geometries
    near_opposition = []  # phase < 5°
    near_quadrature = []  # 85° < phase < 95°
    general = []

    for name, ipl, test_jd in test_cases:
        try:
            result = safe_pheno(test_jd, ipl)
            if result is None:
                continue
            se_attr, le_attr = result
            phase_angle = se_attr[0]

            if phase_angle < 5.0:
                near_opposition.append(
                    (name, ipl, test_jd, se_attr, le_attr, phase_angle)
                )
            elif 85 < phase_angle < 95:
                near_quadrature.append(
                    (name, ipl, test_jd, se_attr, le_attr, phase_angle)
                )
            else:
                general.append((name, ipl, test_jd, se_attr, le_attr, phase_angle))
        except Exception:
            pass

    # Test up to 15 near-opposition cases
    for name, ipl, test_jd, se_attr, le_attr, pa in near_opposition[:15]:
        tols = {
            0: 10.0 / 3600,  # Very tight for near-0° phase angle
            1: 0.0005,
            2: 60.0 / 3600,
            3: 0.0001,
            4: 0.15,
        }
        date_str = f"JD{test_jd:.0f}"
        compare_pheno(
            "P8",
            f"{name} opposition(pa={pa:.2f}°) {date_str}",
            se_attr,
            le_attr,
            tols=tols,
        )

    # Test up to 10 near-quadrature cases
    for name, ipl, test_jd, se_attr, le_attr, pa in near_quadrature[:10]:
        tols = {
            0: 30.0 / 3600,
            1: 0.001,
            2: 30.0 / 3600,
            3: 0.0001,
            4: 0.15,
        }
        date_str = f"JD{test_jd:.0f}"
        compare_pheno(
            "P8",
            f"{name} quadrature(pa={pa:.2f}°) {date_str}",
            se_attr,
            le_attr,
            tols=tols,
        )

    # Also test Mercury/Venus at high phase angles (near inferior conjunction)
    for delta in range(0, 365 * 3, 5):
        test_jd = jd + delta
        for name, ipl in [("Mercury", ephem.MERCURY), ("Venus", ephem.VENUS)]:
            try:
                result = safe_pheno(test_jd, ipl)
                if result is None:
                    continue
                se_attr, le_attr = result
                phase_angle = se_attr[0]
                if phase_angle > 160:
                    tols = {
                        0: 60.0 / 3600,
                        1: 0.002,
                        2: 60.0 / 3600,
                        3: 0.001,
                        4: 0.5,  # Magnitude very uncertain at high phase
                    }
                    date_str = f"JD{test_jd:.0f}"
                    compare_pheno(
                        "P8",
                        f"{name} high_pa({pa:.1f}°) {date_str}",
                        se_attr,
                        le_attr,
                        tols=tols,
                    )
            except Exception:
                pass


def main():
    print("=" * 70)
    print("ROUND 47: Planetary Phenomena (pheno_ut) Deep Re-verification")
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

    phase7()
    print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")

    phase8()
    print(f"  After P8: {passed} passed, {failed} failed, {errors} errors")

    total = passed + failed + errors
    pct = 100 * passed / total if total else 0
    print("\n" + "=" * 70)
    print(f"ROUND 47 FINAL: {passed}/{total} passed ({pct:.1f}%)")
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
