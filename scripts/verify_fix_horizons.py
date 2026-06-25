"""Verify the Horizons HTTP backend against the local LEB/Skyfield pipeline.

Review finding C2 (REVIEW-2026-06-10.md): every fetch-based Horizons
calculation crashed with TypeError; behind that sat a deflection constant
wrong by 3e6, a nonexistent `ayanamsha` module import, and a missing
REF_PLANE in the VECTORS query (23.4 deg frame error).

This script needs INTERNET access (queries the real JPL Horizons API).
It compares horizons_calc_ut directly against the default calc_ut backend
(LEB/Skyfield, both DE440/DE441-based, expected agreement ~milliarcsec).

Usage: .venv/bin/python scripts/verify_fix_horizons.py
"""

from __future__ import annotations

import sys

import libephemeris as eph
from libephemeris.constants import (
    FLG_EQUATORIAL,
    FLG_HELCTR,
    FLG_J2000,
    FLG_SPEED,
    FLG_SWIEPH,
)
from libephemeris.horizons_backend import HorizonsClient, horizons_calc_ut

JD_UT = 2460500.5  # 2024-07-09
ARCSEC = 1.0 / 3600.0

CASES = [
    ("Mars geocentric", 4, FLG_SWIEPH, 0.05),
    ("Mars geocentric + speed", 4, FLG_SWIEPH | FLG_SPEED, 0.05),
    ("Jupiter heliocentric", 5, FLG_SWIEPH | FLG_HELCTR, 0.05),
    ("Venus equatorial", 3, FLG_SWIEPH | FLG_EQUATORIAL, 0.05),
    ("Mercury J2000", 2, FLG_SWIEPH | FLG_J2000, 0.05),
    ("Saturn geocentric", 6, FLG_SWIEPH, 0.05),
    ("Moon geocentric", 1, FLG_SWIEPH, 0.5),  # Moon: topocentric-scale effects
]


def main() -> int:
    client = HorizonsClient()
    failures = 0
    for name, body, flags, tol_arcsec in CASES:
        ref, _ = eph.calc_ut(JD_UT, body, flags)
        try:
            got, _ = horizons_calc_ut(client, JD_UT, body, flags)
        except Exception as exc:  # noqa: BLE001
            print(f"FAIL {name}: raised {type(exc).__name__}: {exc}")
            failures += 1
            continue
        dlon = abs(got[0] - ref[0])
        dlon = min(dlon, 360.0 - dlon) * 3600.0
        dlat = abs(got[1] - ref[1]) * 3600.0
        ddist_rel = abs(got[2] - ref[2]) / max(abs(ref[2]), 1e-12)
        speed_note = ""
        speed_ok = True
        if flags & FLG_SPEED:
            dspd = abs(got[3] - ref[3]) * 3600.0
            speed_note = f" dspeed={dspd:.4f}\"/d"
            # 0.5"/day allows the known LEB-analytic vs apparent-derivative
            # difference (~0.1-0.35"/day) while catching real regressions.
            speed_ok = dspd <= 0.5
        ok = dlon <= tol_arcsec and dlat <= tol_arcsec and ddist_rel < 1e-6 and speed_ok
        status = "ok  " if ok else "FAIL"
        if not ok:
            failures += 1
        print(
            f"{status} {name}: dlon={dlon:.4f}\" dlat={dlat:.4f}\" "
            f"ddist_rel={ddist_rel:.2e}{speed_note}"
        )
    if failures:
        print(f"{failures} case(s) failed")
        return 1
    print("PASS — Horizons pipeline matches the local backend")
    return 0


if __name__ == "__main__":
    sys.exit(main())
