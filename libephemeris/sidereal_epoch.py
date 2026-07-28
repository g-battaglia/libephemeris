# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Fixed-epoch sidereal modes (SIDM_J2000 / J1900 / B1950 / GALALIGN_MARDYKS).

These modes express positions on the mean ecliptic and equinox of a reference
epoch rather than applying only a scalar longitude offset:

- ``SIDM_J2000``: bit-identical to a ``FLG_J2000 | FLG_NONUT`` request —
  positions, speeds, equatorial and XYZ representations alike.
- ``SIDM_J1900`` / ``SIDM_B1950``: the ``FLG_J2000 | FLG_NONUT`` position
  additionally precessed (Vondrák 2011 chain) to the mean frame of t0 and,
  for ecliptic output, rotated onto the mean ecliptic of t0 (IAU 2006 mean
  obliquity at t0).
- ``SIDM_GALALIGN_MARDYKS``: the same t0-frame construction anchored at the
  mode's defining epoch (September equinox 1998), with the defining 30
  degree ayanamsha subtracted from the ecliptic channels (spherical
  longitude and ecliptic XYZ alike). The equatorial channels are the plain
  mean equator/equinox of t0 without the longitude offset.
The returned flags echo the caller's flags with ``FLG_NONUT`` added (the
reference does the same; the internal ``FLG_J2000`` bit is not echoed).

Callers rewrite the request via :func:`fixed_epoch_request_flags`, run the
normal machinery, then map the result back with
:func:`transform_fixed_epoch_result`.

Provenance:
    Frame transformations use ERFA/IAU and the registered Vondrak 2011
    precession chain. J2000, J1900, and B1950 epochs are public astronomical
    conventions; request rewriting and flag echoing are project API plumbing.
    No sidereal zero point, polynomial, or fitted compatibility correction is
    introduced here.
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import Tuple

import erfa
import numpy as np

from .ayanamsha_definitions import MARDYKS_DEFINING
from .constants import (
    FLG_EQUATORIAL,
    FLG_J2000,
    FLG_NONUT,
    FLG_RADIANS,
    FLG_SIDEREAL,
    FLG_XYZ,
    SIDBIT_ECL_T0,
    SIDBIT_SSY_PLANE,
    SIDM_B1950,
    SIDM_GALALIGN_MARDYKS,
    SIDM_J1900,
    SIDM_J2000,
)
from .precession_vondrak import vondrak_precession_matrix

# Solar-system invariable plane relative to the J2000 mean ecliptic and equinox:
# inclination and ascending-node longitude from Souami & Souchay (2012),
# "New determination of the invariable plane of the Solar System", A&A 543,
# A133 (i = 1.578690 deg, node = 107.582222 deg). Used by SIDBIT_SSY_PLANE.
INVARIABLE_PLANE_INCL_DEG = 1.578690
INVARIABLE_PLANE_NODE_DEG = 107.582222


def _julian_epoch_jd(epoch: float) -> float:
    """Return a Julian-epoch date using the ERFA/SOFA definition."""
    date1, date2 = erfa.epj2jd(epoch)
    return float(date1 + date2)


def _besselian_epoch_jd(epoch: float) -> float:
    """Return a Besselian-epoch date using the ERFA/SOFA definition."""
    date1, date2 = erfa.epb2jd(epoch)
    return float(date1 + date2)


_J2000 = _julian_epoch_jd(2000.0)

# Reference epochs of the fixed-epoch ayanamsha modes.
FIXED_EPOCH_T0 = {
    SIDM_J2000: _J2000,
    SIDM_J1900: _julian_epoch_jd(1900.0),
    SIDM_B1950: _besselian_epoch_jd(1950.0),
    SIDM_GALALIGN_MARDYKS: MARDYKS_DEFINING[1],
}

# Constant ecliptic-longitude offset of the t0 frame in degrees (the mode's
# defining ayanamsha at t0). Subtracted from the ecliptic channels only; the
# equatorial output stays the plain mean equator/equinox of t0.
FIXED_EPOCH_LON_OFFSET = {
    SIDM_GALALIGN_MARDYKS: MARDYKS_DEFINING[0],
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
    lon_offset = FIXED_EPOCH_LON_OFFSET.get(sid_mode, 0.0)
    if t0 == _J2000 and lon_offset == 0.0:
        eye = np.eye(3)
        return eye, eye
    v_t0 = np.array(vondrak_precession_matrix(t0))
    v_j2k = np.array(vondrak_precession_matrix(_J2000))
    m_eq = v_t0 @ v_j2k.T  # ICRS frame bias cancels in the chain
    eps_j2k = erfa.obl06(_J2000, 0.0)
    eps_t0 = erfa.obl06(t0, 0.0)
    m_ecl = _rot_x(eps_t0) @ m_eq @ _rot_x(-eps_j2k)
    if lon_offset:
        m_ecl = _rot_z(math.radians(lon_offset)) @ m_ecl
    return m_eq, m_ecl


def _rot_x(angle_rad: float) -> np.ndarray:
    """Return the passive x-axis rotation matrix used for ecliptic frames."""
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]])


def _rot_z(angle_rad: float) -> np.ndarray:
    """Return the passive z-axis rotation (subtracts the angle from longitudes)."""
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


def _invariable_plane_matrix() -> np.ndarray:
    """J2000 mean ecliptic -> solar-system invariable-plane frame.

    Passive rotation: ascending node onto the x-axis, then tilt by the
    inclination so the invariable plane becomes the new x-y plane.
    """
    return _rot_x(math.radians(INVARIABLE_PLANE_INCL_DEG)) @ _rot_z(
        math.radians(INVARIABLE_PLANE_NODE_DEG)
    )


def _ecliptic_of_t0_matrix(t0_jd: float) -> np.ndarray:
    """J2000 mean ecliptic -> mean ecliptic and equinox of t0 (IAU 2006)."""
    m_eq = np.array(erfa.pmat06(_J2000, t0_jd - _J2000))
    eps0 = erfa.obl06(_J2000, 0.0)
    epst = erfa.obl06(t0_jd, 0.0)
    return _rot_x(epst) @ m_eq @ _rot_x(-eps0)


def sidbit_ecliptic_matrix(
    sid_bits: int, t0_jd: float, zero_point_deg: float
) -> np.ndarray | None:
    """Ecliptic-frame rotation for a SIDBIT projection, or None.

    The returned matrix maps a J2000 mean-ecliptic Cartesian state to the
    projection frame with the sidereal zero point rotated onto the x-axis, so
    that applying it and reading spherical longitude yields the sidereal
    longitude directly. ``zero_point_deg`` is the longitude of the sidereal
    zero point *in the projection frame* (for SIDBIT_ECL_T0 it is the mode's
    ayanamsha at t0; for SIDBIT_SSY_PLANE it is the projection of the J2000
    ayanamsha direction onto the invariable plane, computed by the caller).

    Args:
        sid_bits: The active SIDBIT projection flags.
        t0_jd: Reference epoch (JD TT) for SIDBIT_ECL_T0.
        zero_point_deg: Sidereal zero-point longitude in the projection frame.

    Returns:
        A 3x3 numpy rotation matrix, or None when no projection bit is set.
    """
    # Measured reference behavior: with BOTH projection bits set,
    # SIDBIT_ECL_T0 takes precedence over SIDBIT_SSY_PLANE (must mirror the
    # callers' precedence, which choose t0_jd/zero_point the same way).
    if sid_bits & SIDBIT_ECL_T0:
        base = _ecliptic_of_t0_matrix(t0_jd)
    elif sid_bits & SIDBIT_SSY_PLANE:
        base = _invariable_plane_matrix()
    else:
        return None
    return _rot_z(math.radians(zero_point_deg)) @ base


def ssy_plane_zero_point_deg(ayan_j2000_deg: float) -> float:
    """Invariable-plane longitude of the J2000 ayanamsha zero-point direction.

    The sidereal zero point sits at J2000 mean-ecliptic longitude
    ``ayan_j2000_deg`` (latitude 0); projecting that unit direction onto the
    invariable plane gives the longitude from which SIDBIT_SSY_PLANE sidereal
    longitudes are measured.
    """
    m = _invariable_plane_matrix()
    a = math.radians(ayan_j2000_deg)
    p = m @ np.array([math.cos(a), math.sin(a), 0.0])
    return math.degrees(math.atan2(p[1], p[0])) % 360.0


def transform_sidbit_result(
    xx: Tuple[float, ...], iflag: int, m_ecl: np.ndarray
) -> Tuple[float, ...]:
    """Project a J2000|NONUT *ecliptic* result onto a SIDBIT frame.

    Only ecliptic output is handled (spherical or FLG_XYZ); equatorial requests
    are not routed here. ``m_ecl`` already folds in the sidereal zero-point
    rotation, so the rotated longitude is the sidereal longitude. FLG_RADIANS
    is re-applied at the end for spherical output.

    Args:
        xx: The 6-tuple from the rewritten J2000|NONUT request (degrees/AU;
            FLG_RADIANS was stripped from the rewrite).
        iflag: The original caller flags (for FLG_XYZ / FLG_RADIANS).
        m_ecl: Ecliptic-frame rotation from :func:`sidbit_ecliptic_matrix`.

    Returns:
        The projected 6-tuple.
    """
    if iflag & FLG_XYZ:
        pos = m_ecl @ np.array(xx[:3])
        vel = m_ecl @ np.array(xx[3:6])
        return tuple(float(c) for c in pos) + tuple(float(c) for c in vel)
    out = _rotate_spherical(tuple(xx), m_ecl)  # type: ignore[arg-type]
    if iflag & FLG_RADIANS:
        out = (
            math.radians(out[0]),
            math.radians(out[1]),
            out[2],
            math.radians(out[3]),
            math.radians(out[4]),
            out[5],
        )
    return out


def equatorial_epoch_matrix(t0_jd: float) -> np.ndarray:
    """J2000 mean equatorial -> mean equatorial and equinox of ``t0_jd``.

    Same construction as the fixed-epoch chain in :func:`_epoch_matrices`
    (Vondrak et al. 2011 long-term precession; the ICRS frame bias cancels
    in the product), for an arbitrary epoch instead of a mode-table one.
    """
    if t0_jd == _J2000:
        return np.eye(3)
    v_t0 = np.array(vondrak_precession_matrix(t0_jd))
    v_j2k = np.array(vondrak_precession_matrix(_J2000))
    return v_t0 @ v_j2k.T


def transform_equatorial_epoch_result(
    xx: Tuple[float, ...], iflag: int, t0_jd: float
) -> Tuple[float, ...]:
    """Map a FLG_J2000|FLG_NONUT *equatorial* result onto the mean equator
    and equinox of ``t0_jd``.

    ``xx`` is the raw 6-tuple from the rewritten request (degrees / AU;
    FLG_RADIANS was stripped from the rewrite). Handles the spherical and
    FLG_XYZ representations and re-applies FLG_RADIANS at the end. The
    rotation is time-independent, so position and velocity rotate with the
    same matrix.
    """
    m_eq = equatorial_epoch_matrix(t0_jd)
    if iflag & FLG_XYZ:
        pos = m_eq @ np.array(xx[:3])
        vel = m_eq @ np.array(xx[3:6])
        return tuple(float(c) for c in pos) + tuple(float(c) for c in vel)
    out = _rotate_spherical(tuple(xx), m_eq)  # type: ignore[arg-type]
    if iflag & FLG_RADIANS:
        out = (
            math.radians(out[0]),
            math.radians(out[1]),
            out[2],
            math.radians(out[3]),
            math.radians(out[4]),
            out[5],
        )
    return out


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
    if sid_mode != SIDM_J2000:
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
