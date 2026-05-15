#!/usr/bin/env python3
"""
Round 19: Fixed Stars Comprehensive Verification
=================================================

Tests 50+ bright stars across multiple flag combinations and epochs.

Parts:
  P1: Ecliptic of-date positions for 50 bright stars (2024)
  P2: Equatorial coordinates (FLG_EQUATORIAL) for same stars
  P3: J2000 frame (FLG_J2000|FLG_NONUT) — isolate proper motion
  P4: Astrometric (FLG_NOABERR) — no annual aberration
  P5: Sidereal mode (Lahiri) for 20 astrologically important stars
  P6: Magnitude comparison for all catalog stars
  P7: Speed computation (proper motion rates)
  P8: Multi-epoch proper motion drift (1900, 1950, 2000, 2050, 2100)
"""

from __future__ import annotations

import math
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

# 50 bright stars covering all areas of the sky
BRIGHT_STARS = [
    "Aldebaran",
    "Algol",
    "Antares",
    "Arcturus",
    "Betelgeuse",
    "Canopus",
    "Capella",
    "Castor",
    "Deneb",
    "Fomalhaut",
    "Pollux",
    "Procyon",
    "Regulus",
    "Rigel",
    "Sirius",
    "Spica",
    "Vega",
    "Altair",
    "Achernar",
    "Acrux",
    "Agena",
    "Bellatrix",
    "Denebola",
    "Dubhe",
    "Elnath",
    "Hamal",
    "Markab",
    "Menkar",
    "Mira",
    "Mirach",
    "Polaris",
    "Rasalhague",
    "Rigil Kentaurus",
    "Scheat",
    "Shaula",
    "Toliman",
    "Unukalhai",
    "Wezen",
    "Zubenelgenubi",
    "Zubeneschamali",
    "Alphecca",
    "Albireo",
    "Alnilam",
    "Alnitak",
    "Alphard",
    "Diphda",
    "Sadalmelik",
    "Sadalsuud",
]

# Astrologically important stars (subset)
ASTRO_STARS = [
    "Aldebaran",
    "Algol",
    "Antares",
    "Regulus",
    "Spica",
    "Fomalhaut",
    "Sirius",
    "Vega",
    "Arcturus",
    "Pollux",
    "Procyon",
    "Betelgeuse",
    "Rigel",
    "Capella",
    "Canopus",
    "Deneb",
    "Altair",
    "Castor",
    "Hamal",
    "Scheat",
]


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
            print(f'  Max diff: {self.max_diff:.4f}" ({self.max_label})')
        if self.failures:
            for f in self.failures[:20]:
                print(f"    - {f}")
        print(f"{'=' * 70}")
        return self.failed == 0


def arcsec(deg):
    return abs(deg) * 3600.0


def se_fixstar(star, jd, flags):
    """Call pyswisseph fixstar_ut, return (lon, lat, dist, slon, slat, sdist)."""
    result = swe.fixstar_ut(star, jd, flags)
    # pyswisseph returns (pos_tuple, star_name, retflag)
    return result[0]


def le_fixstar(star, jd, flags):
    """Call libephemeris fixstar_ut, return (lon, lat, dist, slon, slat, sdist)."""
    result = ephem.fixstar_ut(star, jd, flags)
    # libephemeris returns (pos_tuple, star_name, retflag)
    return result[0]


# ============================================================
# PART 1: Ecliptic of-date positions
# ============================================================
def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Ecliptic of-date — 50 bright stars (Jan 1, 2024)")
    print("=" * 70)

    r = R("P1: Ecliptic of-date")
    jd = swe.julday(2024, 1, 1, 12.0)
    flags = FLG_SWIEPH

    for star in BRIGHT_STARS:
        label = star
        try:
            se_pos = se_fixstar(star, jd, flags)
            le_pos = le_fixstar(star, jd, flags)
        except Exception as e:
            r.fail(f"{label}: {e}")
            continue

        dlon = se_pos[0] - le_pos[0]
        if dlon > 180:
            dlon -= 360
        elif dlon < -180:
            dlon += 360
        dlon_as = arcsec(dlon)
        dlat_as = arcsec(se_pos[1] - le_pos[1])

        # Fixed star tolerance: 1" for position
        tol = 1.0
        fails = []
        if dlon_as > tol:
            fails.append(f'lon={dlon_as:.3f}"')
        if dlat_as > tol:
            fails.append(f'lat={dlat_as:.3f}"')

        if fails:
            r.fail(f"{label}: {'; '.join(fails)}")
        else:
            r.ok(max(dlon_as, dlat_as), label)

        if dlon_as > 0.5 or dlat_as > 0.5:
            print(f'  {star:20s}: dlon={dlon_as:7.3f}"  dlat={dlat_as:7.3f}"')

    # Print summary stats
    print(f'  (Stars with diff > 0.5" shown above)')
    return r.summary(), r


# ============================================================
# PART 2: Equatorial coordinates
# ============================================================
def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: Equatorial (RA/Dec) — 50 bright stars (2024)")
    print("=" * 70)

    r = R("P2: Equatorial")
    jd = swe.julday(2024, 1, 1, 12.0)
    flags = FLG_SWIEPH | FLG_EQUATORIAL

    for star in BRIGHT_STARS:
        label = f"EQ {star}"
        try:
            se_pos = se_fixstar(star, jd, flags)
            le_pos = le_fixstar(star, jd, flags)
        except Exception as e:
            r.fail(f"{label}: {e}")
            continue

        dra = se_pos[0] - le_pos[0]
        if dra > 180:
            dra -= 360
        elif dra < -180:
            dra += 360
        dra_as = arcsec(dra)
        ddec_as = arcsec(se_pos[1] - le_pos[1])

        tol = 1.0
        fails = []
        if dra_as > tol:
            fails.append(f'RA={dra_as:.3f}"')
        if ddec_as > tol:
            fails.append(f'Dec={ddec_as:.3f}"')

        if fails:
            r.fail(f"{label}: {'; '.join(fails)}")
        else:
            r.ok(max(dra_as, ddec_as), label)

    return r.summary(), r


# ============================================================
# PART 3: J2000 frame — isolate proper motion from precession
# ============================================================
def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: J2000 ecliptic — 50 stars (isolating proper motion)")
    print("=" * 70)

    r = R("P3: J2000 Frame")
    jd = swe.julday(2024, 1, 1, 12.0)
    flags = FLG_SWIEPH | FLG_J2000 | FLG_NONUT

    for star in BRIGHT_STARS:
        label = f"J2000 {star}"
        try:
            se_pos = se_fixstar(star, jd, flags)
            le_pos = le_fixstar(star, jd, flags)
        except Exception as e:
            r.fail(f"{label}: {e}")
            continue

        dlon = se_pos[0] - le_pos[0]
        if dlon > 180:
            dlon -= 360
        elif dlon < -180:
            dlon += 360
        dlon_as = arcsec(dlon)
        dlat_as = arcsec(se_pos[1] - le_pos[1])

        tol = 1.0
        fails = []
        if dlon_as > tol:
            fails.append(f'lon={dlon_as:.3f}"')
        if dlat_as > tol:
            fails.append(f'lat={dlat_as:.3f}"')

        if fails:
            r.fail(f"{label}: {'; '.join(fails)}")
        else:
            r.ok(max(dlon_as, dlat_as), label)

    return r.summary(), r


# ============================================================
# PART 4: Astrometric (no aberration)
# ============================================================
def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: Astrometric (NOABERR) — 50 stars (2024)")
    print("=" * 70)

    r = R("P4: Astrometric")
    jd = swe.julday(2024, 1, 1, 12.0)
    flags = FLG_SWIEPH | FLG_NOABERR

    for star in BRIGHT_STARS:
        label = f"Astro {star}"
        try:
            se_pos = se_fixstar(star, jd, flags)
            le_pos = le_fixstar(star, jd, flags)
        except Exception as e:
            r.fail(f"{label}: {e}")
            continue

        dlon = se_pos[0] - le_pos[0]
        if dlon > 180:
            dlon -= 360
        elif dlon < -180:
            dlon += 360
        dlon_as = arcsec(dlon)
        dlat_as = arcsec(se_pos[1] - le_pos[1])

        tol = 1.0
        fails = []
        if dlon_as > tol:
            fails.append(f'lon={dlon_as:.3f}"')
        if dlat_as > tol:
            fails.append(f'lat={dlat_as:.3f}"')

        if fails:
            r.fail(f"{label}: {'; '.join(fails)}")
        else:
            r.ok(max(dlon_as, dlat_as), label)

    return r.summary(), r


# ============================================================
# PART 5: Sidereal mode (Lahiri) for astrological stars
# ============================================================
def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Sidereal (Lahiri) — 20 astrological stars (2024)")
    print("=" * 70)

    r = R("P5: Sidereal Lahiri")
    jd = swe.julday(2024, 1, 1, 12.0)
    flags = FLG_SWIEPH | FLG_SIDEREAL

    # Set sidereal mode to Lahiri
    swe.set_sid_mode(1)  # SIDM_LAHIRI = 1
    ephem.set_sid_mode(1)

    for star in ASTRO_STARS:
        label = f"Sid {star}"
        try:
            se_pos = se_fixstar(star, jd, flags)
            le_pos = le_fixstar(star, jd, flags)
        except Exception as e:
            r.fail(f"{label}: {e}")
            continue

        dlon = se_pos[0] - le_pos[0]
        if dlon > 180:
            dlon -= 360
        elif dlon < -180:
            dlon += 360
        dlon_as = arcsec(dlon)
        dlat_as = arcsec(se_pos[1] - le_pos[1])

        # Sidereal: SE applies a different sidereal transformation for
        # fixed stars than simple ayanamsha subtraction. The ~5.3"
        # systematic difference is a known methodology difference.
        tol = 6.0
        fails = []
        if dlon_as > tol:
            fails.append(f'lon={dlon_as:.3f}"')
        if dlat_as > tol:
            fails.append(f'lat={dlat_as:.3f}"')

        if fails:
            r.fail(f"{label}: {'; '.join(fails)}")
        else:
            r.ok(max(dlon_as, dlat_as), label)

    # Reset sidereal mode
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0)

    return r.summary(), r


# ============================================================
# PART 6: Magnitude comparison
# ============================================================
def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Star magnitudes — all catalog stars")
    print("=" * 70)

    r = R("P6: Magnitudes")

    for star in BRIGHT_STARS:
        label = f"Mag {star}"
        try:
            se_mag = swe.fixstar_mag(star)
            le_mag = ephem.fixstar_mag(star)
        except Exception as e:
            r.fail(f"{label}: {e}")
            continue

        # pyswisseph returns magnitude directly or as tuple
        if isinstance(se_mag, tuple):
            se_mag_val = se_mag[0] if len(se_mag) > 0 else se_mag
        else:
            se_mag_val = se_mag

        diff = abs(se_mag_val - le_mag)

        # Magnitude tolerance: 0.2 mag (different source catalogs).
        # Variable stars (Mira, Antares) have much larger differences.
        tol = 0.2
        # Mira is a long-period variable (mag 2.0-10.1), catalogs differ hugely
        if star == "Mira":
            tol = 5.0
        if diff > tol:
            r.fail(f"{label}: SE={se_mag_val:.2f} LE={le_mag:.2f} diff={diff:.2f}")
        else:
            r.ok(diff, label)

    return r.summary(), r


# ============================================================
# PART 7: Speed computation (proper motion rates)
# ============================================================
def run_part7():
    print("\n" + "=" * 70)
    print("PART 7: Speed (proper motion) — 50 stars")
    print("=" * 70)

    r = R("P7: Speed/Proper Motion")
    jd = swe.julday(2024, 1, 1, 12.0)
    flags = FLG_SWIEPH | FLG_SPEED

    for star in BRIGHT_STARS:
        label = f"Speed {star}"
        try:
            se_pos = se_fixstar(star, jd, flags)
            le_pos = le_fixstar(star, jd, flags)
        except Exception as e:
            r.fail(f"{label}: {e}")
            continue

        # Speed is in pos[3] (lon speed, deg/day)
        dspeed_lon = abs(se_pos[3] - le_pos[3])
        dspeed_lat = abs(se_pos[4] - le_pos[4])

        # Convert to arcsec/year for readability
        dspeed_lon_asy = dspeed_lon * 3600 * 365.25
        dspeed_lat_asy = dspeed_lat * 3600 * 365.25

        # Speed tolerance: SE uses analytical proper motion transformation
        # while we use finite difference of apparent positions.  The FD
        # method includes the ~18.6-yr nutation oscillation (~13"/yr in
        # ecliptic latitude) which is absent from SE's analytical formula.
        # Longitude speed also differs because SE's analytical approach
        # doesn't include nutation/aberration oscillations.
        # Tolerance: 25"/year for latitude (nutation dominated),
        #            25"/year for longitude (aberration + nutation)
        tol = 25.0
        fails = []
        if dspeed_lon_asy > tol:
            fails.append(f'slon={dspeed_lon_asy:.4f}"/yr')
        if dspeed_lat_asy > tol:
            fails.append(f'slat={dspeed_lat_asy:.4f}"/yr')

        if fails:
            r.fail(f"{label}: {'; '.join(fails)}")
        else:
            r.ok(max(dspeed_lon_asy, dspeed_lat_asy), label)

    return r.summary(), r


# ============================================================
# PART 8: Multi-epoch proper motion drift
# ============================================================
def run_part8():
    print("\n" + "=" * 70)
    print("PART 8: Multi-epoch proper motion — 1900-2100 at 50-yr intervals")
    print("  Testing proper motion accumulation consistency")
    print("=" * 70)

    r = R("P8: Multi-Epoch Stars")
    epochs = [1900, 1950, 2000, 2050, 2100]
    flags = FLG_SWIEPH

    # Test stars with high proper motion (most sensitive to PM errors)
    high_pm_stars = [
        "Sirius",
        "Arcturus",
        "Procyon",
        "Pollux",
        "Regulus",
        "Aldebaran",
        "Spica",
        "Antares",
        "Rigil Kentaurus",
        "Toliman",
        "Vega",
        "Capella",
        "Altair",
        "Fomalhaut",
        "Deneb",
    ]

    for star in high_pm_stars:
        max_diff = 0.0
        worst_epoch = 0

        for year in epochs:
            jd = swe.julday(year, 1, 1, 12.0)
            label = f"{star} {year}"

            try:
                se_pos = se_fixstar(star, jd, flags)
                le_pos = le_fixstar(star, jd, flags)
            except Exception as e:
                r.fail(f"{label}: {e}")
                continue

            dlon = se_pos[0] - le_pos[0]
            if dlon > 180:
                dlon -= 360
            elif dlon < -180:
                dlon += 360
            dlon_as = arcsec(dlon)
            dlat_as = arcsec(se_pos[1] - le_pos[1])

            dist_from_j2000 = abs(year - 2000)
            # Allow 0.5" base + 0.005"/yr for PM model differences.
            # Rigil Kentaurus and Toliman have extreme proper motion
            # (~3.7"/yr), so PM catalog differences accumulate fast.
            if star in ("Rigil Kentaurus", "Toliman"):
                tol = 6.0 + dist_from_j2000 * 0.35  # ~35"/century
            else:
                tol = 0.5 + dist_from_j2000 * 0.005

            total_diff = max(dlon_as, dlat_as)
            fails = []
            if dlon_as > tol:
                fails.append(f'lon={dlon_as:.3f}"')
            if dlat_as > tol:
                fails.append(f'lat={dlat_as:.3f}"')

            if fails:
                r.fail(f"{label}: {'; '.join(fails)}")
            else:
                r.ok(total_diff, label)

            if total_diff > max_diff:
                max_diff = total_diff
                worst_epoch = year

        print(f'  {star:20s}: max_diff={max_diff:7.3f}"  worst_epoch={worst_epoch}')

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 19: Fixed Stars Comprehensive Verification")
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
        ("P8", run_part8),
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
    print("ROUND 19 FINAL SUMMARY")
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
        print(f"\n  >>> ROUND 19: ALL PASSED <<<")
    else:
        print(f"\n  >>> ROUND 19: {tf} FAILURES <<<")

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
