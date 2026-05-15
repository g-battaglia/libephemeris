#!/usr/bin/env python3
"""
Round 16: Deep Planetary Phenomena (pheno_ut) Verification
===============================================================

Parts:
  P1: All planets phase angle/elongation/phase — 5 epochs × 10 bodies
  P2: Apparent diameter comparison — 5 epochs × 10 bodies
  P3: Apparent magnitude comparison — 5 epochs × 10 bodies
  P4: Return structure (20 elements, reserved zeros)
  P5: Phase formula consistency: phase = (1+cos(phase_angle))/2
  P6: TT vs UT consistency
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


BODIES = [
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
]

EPOCHS = [
    (1950, 6, 15, 0.0, "1950"),
    (1980, 3, 21, 12.0, "1980"),
    (2000, 1, 1, 12.0, "J2000"),
    (2020, 6, 21, 0.0, "2020"),
    (2025, 3, 14, 0.0, "2025"),
]


def jd_for(e):
    return swe.julday(e[0], e[1], e[2], e[3])


class R:
    def __init__(self, name):
        self.name = name
        self.passed = self.failed = self.skipped = 0
        self.failures = []
        self.max_diff = 0.0
        self.max_body = ""

    def ok(self, diff=0.0, body=""):
        self.passed += 1
        if diff > self.max_diff:
            self.max_diff = diff
            self.max_body = body

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
            print(f"  Max diff: {self.max_diff:.6f} ({self.max_body})")
        if self.failures:
            for f in self.failures[:15]:
                print(f"    - {f}")
        print(f"{'=' * 70}")
        return self.failed == 0


def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Phase Angle / Phase / Elongation")
    print("  10 bodies × 5 epochs = 50 tests")
    print("=" * 70)

    r = R("P1: Phase/Elongation")

    for epoch in EPOCHS:
        jd = jd_for(epoch)
        for ipl, name in BODIES:
            label = f"{name} @ {epoch[4]}"
            try:
                se_attr = swe.pheno_ut(jd, ipl, 0)
                le_attr, _ = ephem.pheno_ut(jd, ipl, 0)
            except Exception as e:
                r.fail(f"{label}: {e}")
                continue

            # Phase angle (degrees)
            pa_diff = abs(se_attr[0] - le_attr[0])
            # Phase (illuminated fraction)
            ph_diff = abs(se_attr[1] - le_attr[1])
            # Elongation (degrees)
            el_diff = abs(se_attr[2] - le_attr[2])
            if el_diff > 180:
                el_diff = 360 - el_diff

            tol_pa = 0.01  # degrees
            tol_ph = 0.001
            tol_el = 0.01  # degrees

            fails = []
            if pa_diff >= tol_pa:
                fails.append(
                    f"phase_angle: SE={se_attr[0]:.4f} LE={le_attr[0]:.4f} diff={pa_diff:.4f}"
                )
            if ph_diff >= tol_ph:
                fails.append(
                    f"phase: SE={se_attr[1]:.6f} LE={le_attr[1]:.6f} diff={ph_diff:.6f}"
                )
            if el_diff >= tol_el:
                fails.append(
                    f"elongation: SE={se_attr[2]:.4f} LE={le_attr[2]:.4f} diff={el_diff:.4f}"
                )

            if fails:
                r.fail(f"{label}: {'; '.join(fails)}")
            else:
                r.ok(max(pa_diff, el_diff), name)

    return r.summary(), r


def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: Apparent Diameter")
    print("  10 bodies × 5 epochs = 50 tests")
    print("=" * 70)

    r = R("P2: Diameter")

    for epoch in EPOCHS:
        jd = jd_for(epoch)
        for ipl, name in BODIES:
            label = f"{name} @ {epoch[4]}"
            try:
                se_attr = swe.pheno_ut(jd, ipl, 0)
                le_attr, _ = ephem.pheno_ut(jd, ipl, 0)
            except Exception as e:
                r.fail(f"{label}: {e}")
                continue

            se_diam = se_attr[3]
            le_diam = le_attr[3]

            if se_diam > 0 and le_diam > 0:
                rel_diff = abs(se_diam - le_diam) / se_diam
                tol = 0.02  # 2% relative
                if rel_diff >= tol:
                    r.fail(
                        f"{label}: diam SE={se_diam:.8f} LE={le_diam:.8f} "
                        f"rel={rel_diff:.4f}"
                    )
                else:
                    r.ok(rel_diff, name)
            elif se_diam == 0 and le_diam == 0:
                r.ok(0, name)
            else:
                r.fail(f"{label}: diam SE={se_diam} LE={le_diam}")

    return r.summary(), r


def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: Apparent Magnitude")
    print("  10 bodies × 5 epochs = 50 tests")
    print("=" * 70)

    r = R("P3: Magnitude")

    # Magnitude tolerances per body
    mag_tol = {
        "Sun": 0.05,
        "Moon": 0.5,
        "Mercury": 0.5,
        "Venus": 0.3,
        "Mars": 0.3,
        "Jupiter": 0.2,
        "Saturn": 0.5,
        "Uranus": 0.3,
        "Neptune": 0.3,
        "Pluto": 0.5,
    }

    for epoch in EPOCHS:
        jd = jd_for(epoch)
        for ipl, name in BODIES:
            label = f"{name} @ {epoch[4]}"
            try:
                se_attr = swe.pheno_ut(jd, ipl, 0)
                le_attr, _ = ephem.pheno_ut(jd, ipl, 0)
            except Exception as e:
                r.fail(f"{label}: {e}")
                continue

            se_mag = se_attr[4]
            le_mag = le_attr[4]
            mag_diff = abs(se_mag - le_mag)
            tol = mag_tol.get(name, 0.5)

            if mag_diff >= tol:
                r.fail(
                    f"{label}: mag SE={se_mag:.2f} LE={le_mag:.2f} diff={mag_diff:.2f} (tol={tol})"
                )
            else:
                r.ok(mag_diff, name)

    return r.summary(), r


def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: Return Structure (20 elements, reserved zeros)")
    print("=" * 70)

    r = R("P4: Structure")

    jd = jd_for(EPOCHS[2])  # J2000

    for ipl, name in BODIES:
        label = f"{name} structure"
        try:
            se_attr = swe.pheno_ut(jd, ipl, 0)
            le_attr, le_flag = ephem.pheno_ut(jd, ipl, 0)
        except Exception as e:
            r.fail(f"{label}: {e}")
            continue

        # Check length
        if len(le_attr) != 20:
            r.fail(f"{label}: len={len(le_attr)} (expected 20)")
            continue

        # Check reserved zeros (indices 5-19)
        nonzero = [i for i in range(5, 20) if le_attr[i] != 0.0]
        if nonzero:
            r.fail(f"{label}: non-zero reserved at indices {nonzero}")
        else:
            r.ok(0, name)

    return r.summary(), r


def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Phase Formula Consistency")
    print("  phase = (1 + cos(phase_angle)) / 2")
    print("=" * 70)

    r = R("P5: Phase Formula")

    jd = jd_for(EPOCHS[2])

    for ipl, name in BODIES:
        if ipl == SUN:
            r.skip("Sun")
            continue
        label = f"{name} formula"
        try:
            le_attr, _ = ephem.pheno_ut(jd, ipl, 0)
        except Exception as e:
            r.fail(f"{label}: {e}")
            continue

        phase_angle_rad = math.radians(le_attr[0])
        expected_phase = (1.0 + math.cos(phase_angle_rad)) / 2.0
        actual_phase = le_attr[1]
        diff = abs(expected_phase - actual_phase)

        if diff >= 1e-6:
            r.fail(
                f"{label}: expected={expected_phase:.8f} actual={actual_phase:.8f} diff={diff:.8f}"
            )
        else:
            r.ok(diff, name)

    return r.summary(), r


def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: TT vs UT Consistency")
    print("=" * 70)

    r = R("P6: TT/UT")

    for epoch in EPOCHS:
        jd_ut = jd_for(epoch)
        dt = ephem.deltat(jd_ut)
        jd_tt = jd_ut + dt

        for ipl, name in [
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (MOON, "Moon"),
        ]:
            label = f"{name} @ {epoch[4]}"
            try:
                ut_attr, _ = ephem.pheno_ut(jd_ut, ipl, 0)
                tt_attr, _ = ephem.pheno(jd_tt, ipl, 0)
            except Exception as e:
                r.fail(f"{label}: {e}")
                continue

            diff = abs(ut_attr[0] - tt_attr[0])
            if diff >= 1e-4:
                r.fail(
                    f"{label}: UT pa={ut_attr[0]:.6f} TT pa={tt_attr[0]:.6f} diff={diff:.6f}"
                )
            else:
                r.ok(diff, name)

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 16: Deep Planetary Phenomena (pheno_ut) Verification")
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
    print("ROUND 16 FINAL SUMMARY")
    print("=" * 70)

    tp = tf = ts_ = 0
    for pn, res in all_results:
        st = "PASS" if res.failed == 0 else "FAIL"
        t = res.passed + res.failed
        print(f"  {pn} {res.name}: {res.passed}/{t} ({res.skipped} skip) [{st}]")
        tp += res.passed
        tf += res.failed
        ts_ += res.skipped

    print(f"\n  TOTAL: {tp}/{tp + tf} PASSED, {tf} FAILED, {ts_} SKIPPED")
    print(f"  Time: {elapsed:.1f}s")
    if all_ok:
        print("\n  >>> ROUND 16: ALL PASSED <<<")
    else:
        print(f"\n  >>> ROUND 16: {tf} FAILURES <<<")

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
