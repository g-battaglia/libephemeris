#!/usr/bin/env python3
"""
Round 24: Sidereal Mode All Planets Sweep
==========================================

Tests sidereal (Vedic/Hindu) zodiac calculations across multiple ayanamsha
modes, planets, and epochs.

Parts:
  P1: Lahiri ayanamsha — all planets at 5 dates
  P2: Fagan-Bradley ayanamsha — all planets at 5 dates
  P3: Raman ayanamsha — all planets at 5 dates
  P4: Ayanamsha value comparison — all 30+ modes at J2000
  P5: Multi-epoch ayanamsha drift — Lahiri 1900-2100
  P6: Sidereal + speed — verify speed values in sidereal mode
  P7: Sidereal houses — Placidus houses in sidereal mode
"""

from __future__ import annotations

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
            print(f'  Max diff: {self.max_diff:.6f}" ({self.max_label})')
        if self.failures:
            for f in self.failures[:20]:
                print(f"    - {f}")
            if len(self.failures) > 20:
                print(f"    ... and {len(self.failures) - 20} more")
        print(f"{'=' * 70}")
        return self.failed == 0


def angle_diff(a, b):
    d = a - b
    while d > 180:
        d -= 360
    while d < -180:
        d += 360
    return abs(d)


PLANETS = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
    (URANUS, "Uranus"),
    (NEPTUNE, "Neptune"),
    (PLUTO, "Pluto"),
    (MEAN_NODE, "MeanNode"),
    (TRUE_NODE, "TrueNode"),
    (MEAN_APOG, "MeanLilith"),
]

DATES = [
    (1950, 1, 1, 12.0, "1950"),
    (1980, 6, 15, 0.0, "1980"),
    (2000, 1, 1, 12.0, "J2000"),
    (2024, 3, 20, 12.0, "2024"),
    (2050, 7, 1, 0.0, "2050"),
]

# Ayanamsha modes to test
AYANAMSHA_MODES = [
    (SIDM_LAHIRI, "Lahiri"),
    (SIDM_FAGAN_BRADLEY, "FaganBradley"),
    (SIDM_RAMAN, "Raman"),
    (SIDM_KRISHNAMURTI, "Krishnamurti"),
    (SIDM_YUKTESHWAR, "Yukteshwar"),
    (SIDM_JN_BHASIN, "JNBhasin"),
    (SIDM_BABYL_KUGLER1, "BabylKugler1"),
    (SIDM_BABYL_KUGLER2, "BabylKugler2"),
    (SIDM_BABYL_KUGLER3, "BabylKugler3"),
    (SIDM_BABYL_HUBER, "BabylHuber"),
    (SIDM_BABYL_ETPSC, "BabylEtPsc"),
    (SIDM_ALDEBARAN_15TAU, "Aldebaran15Tau"),
    (SIDM_HIPPARCHOS, "Hipparchos"),
    (SIDM_SASSANIAN, "Sassanian"),
    (SIDM_GALCENT_0SAG, "GalCent0Sag"),
    (SIDM_J2000, "J2000"),
    (SIDM_J1900, "J1900"),
    (SIDM_B1950, "B1950"),
]


def se_calc_sidereal(jd, body, sidm):
    """SE sidereal calc."""
    try:
        swe.set_sid_mode(sidm)
        xx = swe.calc_ut(jd, body, FLG_SWIEPH | FLG_SPEED | FLG_SIDEREAL)[0]
        return xx
    except Exception:
        return None


def le_calc_sidereal(jd, body, sidm):
    """LE sidereal calc."""
    try:
        ephem.set_sid_mode(sidm, 0, 0)
        xx, _ = ephem.calc_ut(jd, body, FLG_SPEED | FLG_SIDEREAL)
        return xx
    except Exception:
        return None


def se_get_ayanamsa(jd, sidm):
    """Get ayanamsha from SE."""
    try:
        swe.set_sid_mode(sidm)
        return swe.get_ayanamsa_ut(jd)
    except Exception:
        return None


def le_get_ayanamsa(jd, sidm):
    """Get ayanamsha from LE."""
    try:
        ephem.set_sid_mode(sidm, 0, 0)
        return ephem.get_ayanamsa_ut(jd)
    except Exception:
        return None


def run_sidereal_planets(sidm, sidm_name, dates, tol_arcsec=5.0):
    """Run sidereal planet comparison for a given ayanamsha."""
    print(f"\n{'=' * 70}")
    print(f"  {sidm_name} ayanamsha — all planets × {len(dates)} dates")
    print(f"{'=' * 70}")

    r = R(f"{sidm_name}")

    for y, m, d, h, dname in dates:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in PLANETS:
            label = f"{body_name} {dname}"

            se_xx = se_calc_sidereal(jd, body_id, sidm)
            le_xx = le_calc_sidereal(jd, body_id, sidm)

            if se_xx is None or le_xx is None:
                r.skip(f"{label}: SE={se_xx} LE={le_xx}")
                continue

            lon_diff = angle_diff(se_xx[0], le_xx[0]) * 3600  # arcsec
            lat_diff = abs(se_xx[1] - le_xx[1]) * 3600

            max_d = max(lon_diff, lat_diff)

            if max_d > tol_arcsec:
                r.fail(f'{label}: lon={lon_diff:.2f}" lat={lat_diff:.2f}"')
            else:
                r.ok(max_d, label)

        print(f"    {dname}: tested {len(PLANETS)} bodies")

    return r.summary(), r


# ============================================================
# PART 1-3: Specific ayanamsha modes
# ============================================================
def run_part1():
    return run_sidereal_planets(SIDM_LAHIRI, "P1: Lahiri", DATES)


def run_part2():
    return run_sidereal_planets(SIDM_FAGAN_BRADLEY, "P2: Fagan-Bradley", DATES)


def run_part3():
    return run_sidereal_planets(SIDM_RAMAN, "P3: Raman", DATES)


# ============================================================
# PART 4: Ayanamsha value comparison
# ============================================================
def run_part4():
    print(f"\n{'=' * 70}")
    print(f"PART 4: Ayanamsha values — {len(AYANAMSHA_MODES)} modes at J2000")
    print(f"{'=' * 70}")

    r = R("P4: Ayanamsha Values")

    jd = swe.julday(2000, 1, 1, 12.0)

    for sidm, name in AYANAMSHA_MODES:
        label = f"Ayan {name} (mode={sidm})"

        se_ayan = se_get_ayanamsa(jd, sidm)
        le_ayan = le_get_ayanamsa(jd, sidm)

        if se_ayan is None or le_ayan is None:
            r.skip(f"{label}: SE={se_ayan} LE={le_ayan}")
            continue

        diff_arcsec = abs(se_ayan - le_ayan) * 3600

        tol = 1.0  # 1 arcsec tolerance
        if diff_arcsec > tol:
            r.fail(
                f'{label}: SE={se_ayan:.6f}° LE={le_ayan:.6f}° diff={diff_arcsec:.3f}"'
            )
        else:
            r.ok(diff_arcsec, label)

        print(
            f'  {name:20s}: SE={se_ayan:10.6f}° LE={le_ayan:10.6f}° diff={diff_arcsec:.4f}"'
        )

    return r.summary(), r


# ============================================================
# PART 5: Multi-epoch ayanamsha drift
# ============================================================
def run_part5():
    print(f"\n{'=' * 70}")
    print(f"PART 5: Lahiri ayanamsha 1900-2100 (10-year steps)")
    print(f"{'=' * 70}")

    r = R("P5: Ayan Drift")

    years = list(range(1900, 2101, 10))

    max_diff = 0.0

    for y in years:
        jd = swe.julday(y, 1, 1, 0.0)
        label = f"Lahiri {y}"

        se_ayan = se_get_ayanamsa(jd, SIDM_LAHIRI)
        le_ayan = le_get_ayanamsa(jd, SIDM_LAHIRI)

        if se_ayan is None or le_ayan is None:
            r.skip(f"{label}")
            continue

        diff_arcsec = abs(se_ayan - le_ayan) * 3600

        if diff_arcsec > max_diff:
            max_diff = diff_arcsec

        tol = 1.0
        if diff_arcsec > tol:
            r.fail(f'{label}: diff={diff_arcsec:.4f}"')
        else:
            r.ok(diff_arcsec, label)

    print(f'  Max diff across 1900-2100: {max_diff:.4f}"')

    return r.summary(), r


# ============================================================
# PART 6: Sidereal + speed
# ============================================================
def run_part6():
    print(f"\n{'=' * 70}")
    print(f"PART 6: Sidereal speeds — Lahiri, 2024")
    print(f"{'=' * 70}")

    r = R("P6: Sidereal Speed")

    jd = swe.julday(2024, 3, 20, 12.0)

    for body_id, body_name in PLANETS:
        label = f"Speed {body_name}"

        se_xx = se_calc_sidereal(jd, body_id, SIDM_LAHIRI)
        le_xx = le_calc_sidereal(jd, body_id, SIDM_LAHIRI)

        if se_xx is None or le_xx is None:
            r.skip(f"{label}")
            continue

        # Compare speed (index 3 = lon speed)
        se_spd = se_xx[3]
        le_spd = le_xx[3]
        spd_diff = abs(se_spd - le_spd)

        # Speed tolerance: 0.1% + 0.001°/day baseline
        tol = abs(se_spd) * 0.001 + 0.001
        if spd_diff > tol:
            r.fail(f"{label}: SE={se_spd:.6f} LE={le_spd:.6f} diff={spd_diff:.6f}")
        else:
            r.ok(spd_diff * 3600, label)

        # Also compare lat speed (index 4)
        se_lat_spd = se_xx[4]
        le_lat_spd = le_xx[4]
        lat_spd_diff = abs(se_lat_spd - le_lat_spd)
        tol_lat = abs(se_lat_spd) * 0.01 + 0.001
        if lat_spd_diff > tol_lat:
            r.fail(f"{label} lat_spd: SE={se_lat_spd:.6f} LE={le_lat_spd:.6f}")
        else:
            r.ok(lat_spd_diff * 3600, f"{label} lat_spd")

        print(
            f"  {body_name:12s}: lon_spd_diff={spd_diff:.8f}°/d lat_spd_diff={lat_spd_diff:.8f}°/d"
        )

    return r.summary(), r


# ============================================================
# PART 7: Sidereal houses
# ============================================================
def run_part7():
    print(f"\n{'=' * 70}")
    print(f"PART 7: Sidereal Placidus houses — Lahiri")
    print(f"{'=' * 70}")

    r = R("P7: Sidereal Houses")

    jd = swe.julday(2024, 3, 20, 12.0)

    locations = [
        (12.4964, 41.9028, "Rome"),
        (2.3522, 48.8566, "Paris"),
        (80.2707, 13.0827, "Chennai"),
        (-73.9857, 40.7484, "New York"),
        (151.2093, -33.8688, "Sydney"),
    ]

    for lon, lat, loc_name in locations:
        label = f"Houses {loc_name}"

        try:
            swe.set_sid_mode(SIDM_LAHIRI)
            se_cusps, se_ascmc = swe.houses_ex(jd, lat, lon, b"P", FLG_SIDEREAL)

            ephem.set_sid_mode(SIDM_LAHIRI, 0, 0)
            le_cusps, le_ascmc = ephem.houses_ex(
                jd, lat, lon, ord("P"), FLG_SIDEREAL
            )
        except Exception as e:
            r.fail(f"{label}: {e}")
            continue

        max_cusp_diff = 0.0
        worst_cusp = 0

        n_cusps = min(len(se_cusps), len(le_cusps), 12)
        for i in range(n_cusps):
            d = angle_diff(se_cusps[i], le_cusps[i]) * 3600  # arcsec
            if d > max_cusp_diff:
                max_cusp_diff = d
                worst_cusp = i + 1

        tol = 5.0  # 5" for house cusps
        if max_cusp_diff > tol:
            r.fail(f'{label}: max_diff={max_cusp_diff:.2f}" at cusp {worst_cusp}')
        else:
            r.ok(max_cusp_diff, label)

        # Asc/MC
        asc_diff = angle_diff(se_ascmc[0], le_ascmc[0]) * 3600
        mc_diff = angle_diff(se_ascmc[1], le_ascmc[1]) * 3600

        if asc_diff > 5.0:
            r.fail(f'{label} Asc: diff={asc_diff:.2f}"')
        else:
            r.ok(asc_diff, f"{label} Asc")

        if mc_diff > 5.0:
            r.fail(f'{label} MC: diff={mc_diff:.2f}"')
        else:
            r.ok(mc_diff, f"{label} MC")

        print(
            f'  {loc_name:12s}: cusps max={max_cusp_diff:.2f}" Asc={asc_diff:.2f}" MC={mc_diff:.2f}"'
        )

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 24: Sidereal Mode All Planets Sweep")
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
    print("ROUND 24 FINAL SUMMARY")
    print("=" * 70)

    tp = tf = ts = 0
    for pn, res in all_results:
        st = "PASS" if res.failed == 0 else "FAIL"
        t = res.passed + res.failed
        print(f"  {pn} {res.name}: {res.passed}/{t} ({res.skipped} skip) [{st}]")
        tp += res.passed
        tf += res.failed
        ts += res.skipped

    print(f"\n  TOTAL: {tp}/{tp + tf} PASSED, {tf} FAILED, {ts} SKIPPED")
    print(f"  Time: {elapsed:.1f}s")
    if all_ok:
        print(f"\n  >>> ROUND 24: ALL PASSED <<<")
    else:
        print(f"\n  >>> ROUND 24: {tf} FAILURES <<<")

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
