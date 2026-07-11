# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Fixed-epoch sidereal modes (SIDM_J2000 / J1900 / B1950 / GALALIGN_MARDYKS).

The reference API implements these ayanamsha modes as *frame*
transformations, not as a scalar longitude offset: a sidereal position is
the apparent position expressed on the mean ecliptic and equinox of the
mode's reference epoch t0 (measured black-box, exact to sub-milliarcsecond):

- ``SIDM_J2000``: bit-identical to a ``FLG_J2000 | FLG_NONUT`` request —
  positions, speeds, equatorial and XYZ representations alike.
- ``SIDM_J1900`` / ``SIDM_B1950``: the ``FLG_J2000 | FLG_NONUT`` position
  additionally precessed (Vondrák 2011 chain) to the mean frame of t0 and,
  for ecliptic output, rotated onto the mean ecliptic of t0 (IAU 2006 mean
  obliquity at t0).
- ``SIDM_GALALIGN_MARDYKS``: same frame projection with t0 fitted black-box
  as JD 2451079.771 (the September-1998 "galactic alignment" equinox;
  0.0000" residual over 8 stars, planets, the lunar node and epochs
  1900-2050), plus a constant 30-degree longitude offset on the ecliptic
  and XYZ channels (the mode's ayanamsha is defined as exactly 30 degrees
  at t0). The equatorial channel carries no offset.

The returned flags echo the caller's flags with ``FLG_NONUT`` added (the
reference does the same; the internal ``FLG_J2000`` bit is not echoed).

Callers rewrite the request via :func:`fixed_epoch_request_flags`, run the
normal machinery, then map the result back with
:func:`transform_fixed_epoch_result`.
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import Tuple

import erfa
import numpy as np

from .constants import (
    FLG_EQUATORIAL,
    FLG_J2000,
    FLG_NONUT,
    FLG_RADIANS,
    FLG_SIDEREAL,
    FLG_XYZ,
    SIDM_B1950,
    SIDM_GALALIGN_MARDYKS,
    SIDM_J1900,
    SIDM_J2000,
)
from .precession_vondrak import vondrak_precession_matrix

_J2000 = 2451545.0

# Fitted reference epoch of the Mardyks galactic-alignment mode (see module
# docstring): the September-1998 equinox.
GALALIGN_MARDYKS_T0 = 2451079.771

# Reference epochs of the fixed-epoch ayanamsha modes.
FIXED_EPOCH_T0 = {
    SIDM_J2000: _J2000,
    SIDM_J1900: 2415020.0,
    SIDM_B1950: 2433282.42345905,
    SIDM_GALALIGN_MARDYKS: GALALIGN_MARDYKS_T0,
}

# Constant longitude offset subtracted from the ecliptic/XYZ channels after
# the t0-frame projection (never from the equatorial channel).
FIXED_EPOCH_LON_OFFSET = {
    SIDM_GALALIGN_MARDYKS: 30.0,
}


def is_fixed_epoch_request(iflag: int, sid_mode: int) -> bool:
    """True when this is a FLG_SIDEREAL request under a fixed-epoch mode."""
    return bool(iflag & FLG_SIDEREAL) and sid_mode in FIXED_EPOCH_T0


def fixed_epoch_request_flags(iflag: int) -> int:
    """Rewrite a fixed-epoch sidereal request to its equivalent frame request."""
    return (iflag & ~(FLG_SIDEREAL | FLG_RADIANS)) | FLG_J2000 | FLG_NONUT


def fixed_epoch_retflag(retflag: int, iflag: int) -> int:
    """Echo flags for a fixed-epoch sidereal call.

    The reference echoes the caller's flags plus FLG_NONUT; the internal
    FLG_J2000 rewrite is not exposed. Bits the backend added on top of the
    rewritten request (ephemeris echo bits, implied bits) are kept, minus
    the internal FLG_J2000.
    """
    return (retflag & ~FLG_J2000) | (iflag & (FLG_SIDEREAL | FLG_RADIANS)) | FLG_NONUT


@lru_cache(maxsize=8)
def _epoch_matrices(sid_mode: int) -> Tuple[np.ndarray, np.ndarray]:
    """(M_eq, M_ecl) rotating J2000 frames to the t0 frames of ``sid_mode``.

    M_eq: mean equatorial J2000 -> mean equatorial of t0.
    M_ecl: mean ecliptic J2000 -> mean ecliptic of t0.
    Both are identity for SIDM_J2000.
    """
    t0 = FIXED_EPOCH_T0[sid_mode]
    if t0 == _J2000:
        eye = np.eye(3)
        return eye, eye
    v_t0 = np.array(vondrak_precession_matrix(t0))
    v_j2k = np.array(vondrak_precession_matrix(_J2000))
    m_eq = v_t0 @ v_j2k.T  # ICRS frame bias cancels in the chain
    eps_j2k = erfa.obl06(_J2000, 0.0)
    eps_t0 = erfa.obl06(t0, 0.0)
    m_ecl = _rot_x(eps_t0) @ m_eq @ _rot_x(-eps_j2k)
    return m_eq, m_ecl


def _rot_x(angle_rad: float) -> np.ndarray:
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]])


def _rot_z(angle_rad: float) -> np.ndarray:
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])


def _rotate_spherical(
    xx: Tuple[float, float, float, float, float, float], m: np.ndarray
) -> Tuple[float, float, float, float, float, float]:
    """Rotate a spherical (lon, lat, r, dlon, dlat, dr) state by matrix ``m``.

    The rotation is time-independent, so the Cartesian velocity rotates with
    the same matrix as the position.
    """
    lon, lat, r, dlon, dlat, dr = xx
    lon_r, lat_r = math.radians(lon), math.radians(lat)
    dlon_r, dlat_r = math.radians(dlon), math.radians(dlat)
    cl, sl = math.cos(lon_r), math.sin(lon_r)
    cb, sb = math.cos(lat_r), math.sin(lat_r)
    if r == 0.0:
        r_eff = 1.0
        dr_eff = 0.0
    else:
        r_eff, dr_eff = r, dr
    pos = np.array([r_eff * cb * cl, r_eff * cb * sl, r_eff * sb])
    vel = np.array(
        [
            dr_eff * cb * cl - r_eff * sb * cl * dlat_r - r_eff * cb * sl * dlon_r,
            dr_eff * cb * sl - r_eff * sb * sl * dlat_r + r_eff * cb * cl * dlon_r,
            dr_eff * sb + r_eff * cb * dlat_r,
        ]
    )
    p = m @ pos
    v = m @ vel
    rho2 = p[0] * p[0] + p[1] * p[1]
    rho = math.sqrt(rho2)
    rn = math.sqrt(rho2 + p[2] * p[2])
    lon2 = math.atan2(p[1], p[0]) % (2.0 * math.pi)
    lat2 = math.asin(max(-1.0, min(1.0, p[2] / rn)))
    dlon2 = (p[0] * v[1] - p[1] * v[0]) / rho2 if rho2 > 0.0 else 0.0
    dr2 = (p @ v) / rn
    dlat2 = (v[2] - dr2 * p[2] / rn) / rho if rho > 0.0 else 0.0
    if r == 0.0:
        rn, dr2 = 0.0, 0.0
        if dlon == 0.0 and dlat == 0.0 and dr == 0.0:
            dlon2 = dlat2 = 0.0
    return (
        math.degrees(lon2),
        math.degrees(lat2),
        float(rn),
        math.degrees(dlon2),
        math.degrees(dlat2),
        float(dr2),
    )


def transform_fixed_epoch_result(
    xx: Tuple[float, ...], iflag: int, sid_mode: int
) -> Tuple[float, ...]:
    """Map a FLG_J2000|FLG_NONUT result onto the t0 frame of ``sid_mode``.

    ``xx`` is the raw 6-tuple returned by the rewritten request (degrees /
    AU; FLG_RADIANS was stripped from the rewritten request). Handles the
    ecliptic-spherical, equatorial-spherical and XYZ representations, and
    re-applies FLG_RADIANS at the end when the caller asked for it.
    """
    m_eq, m_ecl = _epoch_matrices(sid_mode)
    # Constant longitude offset (ecliptic/XYZ channels only): folded into the
    # rotation as an extra spin about the t0 ecliptic pole, so positions and
    # velocities pick it up uniformly in every representation.
    off = FIXED_EPOCH_LON_OFFSET.get(sid_mode, 0.0)
    if off and not (iflag & FLG_EQUATORIAL):
        m_ecl = _rot_z(math.radians(off)) @ m_ecl
    if sid_mode != SIDM_J2000 or off:
        if iflag & FLG_XYZ:
            m = m_eq if iflag & FLG_EQUATORIAL else m_ecl
            pos = m @ np.array(xx[:3])
            vel = m @ np.array(xx[3:6])
            xx = tuple(float(c) for c in pos) + tuple(float(c) for c in vel)
        else:
            m = m_eq if iflag & FLG_EQUATORIAL else m_ecl
            xx = _rotate_spherical(tuple(xx), m)  # type: ignore[arg-type]
    if iflag & FLG_RADIANS and not (iflag & FLG_XYZ):
        xx = (
            math.radians(xx[0]),
            math.radians(xx[1]),
            xx[2],
            math.radians(xx[3]),
            math.radians(xx[4]),
            xx[5],
        )
    return xx
