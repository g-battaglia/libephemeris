"""Verify the Lieske two-epoch precession in astrometry._precess_ecliptic.

Compares the function against a reference built from erfa.pmat76 (the IAU
1976/Lieske precession matrix), using the module's own mean obliquity for
the ecliptic<->equatorial rotations so that only the zeta/z/theta matrix is
under test.

Review finding C1 (REVIEW-2026-06-10.md): the quadratic/cubic coefficients
of zeta and z were transposed, giving ~1" error at +-5 cy and ~906" at
+-50 cy from J2000. After the fix the function must match erfa.pmat76 to
<= 0.005" across +-50 centuries.

Usage: .venv/bin/python scripts/verify_fix_precession.py

Provenance:
    Project-authored conformance test against ERFA's BSD-licensed implementation
    of the published IAU 1976/Lieske matrix. ERFA supplies an independent
    standard realization; sampled deltas merely accept or reject local algebra
    and never become fitted precession coefficients.
"""

from __future__ import annotations

import math
import sys

import erfa
import numpy as np

from libephemeris.astrometry import J2000, _mean_obliquity, _precess_ecliptic

DEG = math.pi / 180.0


def reference_precess(lon: float, lat: float, from_jd: float, to_jd: float):
    """Ecliptic precession via erfa.pmat76, mirroring _precess_ecliptic."""
    eps1 = _mean_obliquity(from_jd) * DEG
    eps2 = _mean_obliquity(to_jd) * DEG

    lon_r, lat_r = lon * DEG, lat * DEG
    vec_ecl = np.array(
        [
            math.cos(lat_r) * math.cos(lon_r),
            math.cos(lat_r) * math.sin(lon_r),
            math.sin(lat_r),
        ]
    )
    # ecliptic -> equatorial (rotate about x by -eps1)
    rot1 = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, math.cos(eps1), -math.sin(eps1)],
            [0.0, math.sin(eps1), math.cos(eps1)],
        ]
    )
    vec_eq = rot1 @ vec_ecl
    # precess (erfa matrices multiply column vectors: r2 = P @ r1)
    pmat = np.array(erfa.pmat76(from_jd, 0.0)).T  # J2000 <- from
    pmat2 = np.array(erfa.pmat76(to_jd, 0.0)).T  # J2000 <- to
    vec_eq2 = pmat2.T @ (pmat @ vec_eq)
    # equatorial -> ecliptic of target epoch
    rot2 = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, math.cos(eps2), math.sin(eps2)],
            [0.0, -math.sin(eps2), math.cos(eps2)],
        ]
    )
    out = rot2 @ vec_eq2
    new_lon = math.degrees(math.atan2(out[1], out[0])) % 360.0
    new_lat = math.degrees(math.asin(max(-1.0, min(1.0, out[2]))))
    return new_lon, new_lat


def main() -> int:
    """Conformance-check local two-epoch precession against ERFA ``pmat76``.

    The fixed century/coordinate matrix spans both ordinary and extreme
    conditions.  ERFA is the independently licensed realization of the cited
    IAU 1976/Lieske standard; measured deltas only decide this test verdict.

    Returns:
        Zero when the worst sampled error is at most 0.005 arcseconds.
    """
    worst = 0.0
    worst_case = ""
    centuries = [-50, -20, -5, -1, 0.5, 1, 5, 20, 50]
    for cy in centuries:
        to_jd = J2000 + cy * 36525.0
        for lon in (0.0, 45.0, 123.456, 200.0, 300.0):
            for lat in (-60.0, -5.0, 0.0, 5.0, 60.0):
                got_lon, got_lat = _precess_ecliptic(lon, lat, J2000, to_jd)
                ref_lon, ref_lat = reference_precess(lon, lat, J2000, to_jd)
                dlon = abs(got_lon - ref_lon)
                dlon = min(dlon, 360.0 - dlon) * math.cos(math.radians(ref_lat))
                dlat = abs(got_lat - ref_lat)
                err = math.hypot(dlon, dlat) * 3600.0
                if err > worst:
                    worst = err
                    worst_case = f"cy={cy} lon={lon} lat={lat}"
    print(f"worst |error| vs erfa.pmat76: {worst:.5f} arcsec @ {worst_case}")
    tol = 0.005
    if worst > tol:
        print(f"FAIL (tolerance {tol} arcsec)")
        return 1
    print(f"PASS (tolerance {tol} arcsec)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
