# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Lunar node and apsis calculations.

Mean node and mean apogee delegate to :mod:`libephemeris.mean_lunar_apse`,
which derives their geometry from the IERS 2003 Delaunay arguments implemented
by BSD-licensed ERFA.

Interpolated apogee and perigee are independent, non-antipodal curves anchored
to every DE440 lunar-distance extremum over the medium kernel interval. Their
runtime series and residual tables are regenerated solely by
``scripts/generate_lunar_apse_model.py`` from DE440 states and the same IERS
arguments. No compatibility-ephemeris samples enter that model.

True node, osculating apogee, and osculating perigee come directly from the
active JPL geocentric lunar position and velocity through angular-momentum and
eccentricity-vector geometry.

Inputs are TT Julian days. Mean points begin on the mean ecliptic of date;
interpolated and osculating points begin on the true ecliptic of date. The
shared calculation pipeline applies the requested output-frame and flag
transformations.

Primary sources are IERS Conventions (2010), Chapter 5, Eq. 5.43, and Park et
al. (2021), "The JPL Planetary and Lunar Ephemerides DE440 and DE441",
Astronomical Journal 161:105.

Provenance:
    Source/derivation boundaries are intentionally separate: ERFA evaluates the
    IERS mean arguments; JPL states feed angular-momentum and eccentricity-vector
    geometry; and the checked-in interpolated curves are produced only by
    ``scripts/generate_lunar_apse_model.py``. Grid spacing, residual
    interpolation, edge taper, cache size, and finite-difference speed steps are
    project choices with tests and documented bounds. Compatibility observations
    never enter either generated model or runtime correction.
"""

from __future__ import annotations

import math
import warnings
from functools import lru_cache
from typing import Tuple
from .mean_lunar_apse import (
    EARTH_ECCENTRICITY_J2000,
    EARTH_ECCENTRICITY_T,
    EARTH_ECCENTRICITY_T2,
    _mean_lunar_apogee_position_unchecked,
    _mean_lunar_node_position_unchecked,
    lunar_delaunay_arguments,
)
from .state import _get_computation_ephemeris, get_planets, get_timescale
from .exceptions import EphemerisRangeError

try:
    from .lunar_corrections import (
        PERIGEE_PERTURBATION_CORRECTIONS,
        PERIGEE_CORRECTION_START_YEAR,
        PERIGEE_CORRECTION_END_YEAR,
        PERIGEE_CORRECTION_STEP_YEARS,
    )

    _PERIGEE_CORRECTIONS_AVAILABLE = True
except ImportError:
    _PERIGEE_CORRECTIONS_AVAILABLE = False
    PERIGEE_PERTURBATION_CORRECTIONS = ()
    PERIGEE_CORRECTION_START_YEAR = 0
    PERIGEE_CORRECTION_END_YEAR = 0
    PERIGEE_CORRECTION_STEP_YEARS = 2

try:
    from .lunar_apse_model import (
        APOGEE_CORRECTIONS,
        APOGEE_CT_COUNT,
        APOGEE_CT_JD_START,
        APOGEE_CT_STEP_DAYS,
        APOGEE_LAT_COEFFS,
        APOGEE_MEAN_DIST_AU,
        APOGEE_TRIG_TERMS,
        PERIGEE_CORRECTIONS,
        PERIGEE_CT_COUNT,
        PERIGEE_CT_JD_START,
        PERIGEE_CT_STEP_DAYS,
        PERIGEE_LAT_COEFFS,
        PERIGEE_MEAN_DIST_AU,
        PERIGEE_TRIG_TERMS,
    )

    _APSE_CORRECTIONS_AVAILABLE = True
except ImportError:
    _APSE_CORRECTIONS_AVAILABLE = False
    APOGEE_CT_JD_START = 0.0
    APOGEE_CT_STEP_DAYS = 10.0
    APOGEE_CT_COUNT = 0
    APOGEE_CORRECTIONS = ()
    PERIGEE_CT_JD_START = 0.0
    PERIGEE_CT_STEP_DAYS = 2.0
    PERIGEE_CT_COUNT = 0
    PERIGEE_CORRECTIONS = ()
    APOGEE_TRIG_TERMS = ()
    PERIGEE_TRIG_TERMS = ()
    # Published mean lunar-orbit inclination (Meeus ch. 47) as the
    # single-harmonic fallback; the shipped model carries the DE440 fit.
    APOGEE_LAT_COEFFS = (5.145396, 0.0)
    PERIGEE_LAT_COEFFS = (5.145396, 0.0)
    # DE440 mean apsis distances (see scripts/generate_lunar_apse_model.py).
    APOGEE_MEAN_DIST_AU = 0.0027099776307836668
    PERIGEE_MEAN_DIST_AU = 0.0024236498902102843


# Earth-Moon system gravitational parameter in AU^3/day^2. The two inputs are
# the JPL DE430 ephemeris constants (Folkner et al. 2014, IPN Progress Report
# 42-196): GM_Earth = 398600.435436 km^3/s^2 and Earth/Moon mass ratio
# 81.3005691 (DE430 EMRAT = 81.30056907...). DE440 revised these marginally
# (GM_Earth = 398600.435507, EMRAT = 81.3005682...), a ~2e-10 relative change
# that is negligible for the osculating-element math below, which is the single
# consumer of this constant.
GM_EARTH_MOON_AU3_DAY2: float = (
    (398600.435436 * (1.0 + 1.0 / 81.3005691)) / (149597870.7**3) * (86400.0**2)
)

# Validity range constants for Meeus polynomial approximations
# The polynomials are optimized for dates near J2000.0 (year 2000)
# Precision degrades for dates far from J2000 due to truncated polynomial terms
MEEUS_OPTIMAL_CENTURIES = 2.0  # ±200 years: <0.001° error
MEEUS_VALID_CENTURIES = 10.0  # ±1000 years: <0.01° error
MEEUS_MAX_CENTURIES = 20.0  # ±2000 years: error grows significantly beyond


def _jd_to_year(jd_tt: float) -> float:
    """Convert Julian Day (TT) to year (floating point)."""
    return (jd_tt - 2451545.0) / 365.25 + 2000.0


# Linear taper window applied to correction-table edge values outside the
# table range.  Stepping straight to 0.0 at the boundary used to produce
# discontinuities up to the edge-correction magnitude (hundreds of arcsec
# around 1549/2651 CE) and speed spikes from central differences that
# straddled the boundary; holding the edge value forever would instead
# extrapolate a fit beyond its data.  The taper keeps positions continuous
# and decays the correction to zero within a year of the table edge.
_APSE_EDGE_TAPER_YEARS: float = 1.0
_APSE_EDGE_TAPER_DAYS: float = 365.25


def _interpolated_apse_model_window(body_id: int) -> tuple[float, float] | None:
    """Return the exact inclusive TT window of an interpolated-apse model.

    The bounds come directly from the generated residual-grid metadata in
    ``lunar_apse_model``: the first grid Julian Day and the Julian Day of its
    last sample (``start + step * (count - 1)``).  They describe the model,
    not a binary file, so they remain available when no ``apogee`` group is
    installed.

    Args:
        body_id: ``INTP_APOG`` (21) or ``INTP_PERG`` (22).

    Returns:
        ``(jd_start, jd_end)`` in TT, or ``None`` for another body or when the
        generated residual grid is unavailable.  A caller holding a UT day
        number may compare it directly: the delta-T offset is minutes against
        a grid sampled every few days, and the interpolation tapers rather
        than steps past the last sample.
    """
    from .constants import INTP_APOG, INTP_PERG

    if not _APSE_CORRECTIONS_AVAILABLE:
        return None
    if body_id == INTP_APOG:
        start, step, count = (
            APOGEE_CT_JD_START,
            APOGEE_CT_STEP_DAYS,
            APOGEE_CT_COUNT,
        )
    elif body_id == INTP_PERG:
        start, step, count = (
            PERIGEE_CT_JD_START,
            PERIGEE_CT_STEP_DAYS,
            PERIGEE_CT_COUNT,
        )
    else:
        return None
    if count < 1:
        return None
    return float(start), float(start) + float(step) * float(count - 1)


def _interpolate_perigee_correction(jd_tt: float) -> float:
    """
    Interpolate perigee perturbation correction from precomputed table.

    Uses linear interpolation between table entries. The perigee correction
    table has its own start/end years and step size, independent from the
    mean element correction tables.

    Args:
        jd_tt: Julian Day in TT

    Returns:
        Interpolated correction in degrees.  Outside the table range the
        edge value is linearly tapered to 0.0 over one year (no step).
    """
    if not _PERIGEE_CORRECTIONS_AVAILABLE or not PERIGEE_PERTURBATION_CORRECTIONS:
        return 0.0

    year = _jd_to_year(jd_tt)

    if year < PERIGEE_CORRECTION_START_YEAR:
        weight = 1.0 - (PERIGEE_CORRECTION_START_YEAR - year) / _APSE_EDGE_TAPER_YEARS
        if weight <= 0.0:
            return 0.0
        return float(PERIGEE_PERTURBATION_CORRECTIONS[0]) * weight
    if year > PERIGEE_CORRECTION_END_YEAR:
        weight = 1.0 - (year - PERIGEE_CORRECTION_END_YEAR) / _APSE_EDGE_TAPER_YEARS
        if weight <= 0.0:
            return 0.0
        return float(PERIGEE_PERTURBATION_CORRECTIONS[-1]) * weight

    idx_float = (year - PERIGEE_CORRECTION_START_YEAR) / PERIGEE_CORRECTION_STEP_YEARS
    idx_low = int(idx_float)

    if idx_low < 0:
        return float(PERIGEE_PERTURBATION_CORRECTIONS[0])
    if idx_low >= len(PERIGEE_PERTURBATION_CORRECTIONS) - 1:
        return float(PERIGEE_PERTURBATION_CORRECTIONS[-1])

    frac = idx_float - idx_low
    return float(PERIGEE_PERTURBATION_CORRECTIONS[idx_low]) + frac * (
        float(PERIGEE_PERTURBATION_CORRECTIONS[idx_low + 1])
        - float(PERIGEE_PERTURBATION_CORRECTIONS[idx_low])
    )


def _interpolate_apse_correction(
    jd_tt: float, jd_start: float, step_days: float, corrections: tuple, n: int
) -> float:
    """
    Interpolate apse correction from a JD-indexed correction table.

    Uses linear interpolation between table entries. Returns the correction
    in arcseconds.

    Args:
        jd_tt: Julian Day in TT
        jd_start: JD of the first table entry
        step_days: Spacing between table entries in days
        corrections: Tuple of correction values in arcseconds
        n: Number of entries in the table

    Returns:
        Interpolated correction in arcseconds.  Outside the table range
        the edge value is linearly tapered to 0.0 over one year — a hard
        0.0 step here used to jump positions by the edge correction
        (up to a few hundred arcsec at the 1549/2651 CE table limits).
    """
    if not corrections or n == 0:
        return 0.0

    idx_float = (jd_tt - jd_start) / step_days

    if idx_float < 0:
        dist_days = -idx_float * step_days
        weight = 1.0 - dist_days / _APSE_EDGE_TAPER_DAYS
        if weight <= 0.0:
            return 0.0
        return float(corrections[0]) * weight
    if idx_float >= n - 1:
        dist_days = (idx_float - (n - 1)) * step_days
        weight = 1.0 - dist_days / _APSE_EDGE_TAPER_DAYS
        if weight <= 0.0:
            return 0.0
        return float(corrections[n - 1]) * weight

    idx_low = int(idx_float)
    frac = idx_float - idx_low

    return float(corrections[idx_low]) + frac * (
        float(corrections[idx_low + 1]) - float(corrections[idx_low])
    )


def _calc_lunar_fundamental_arguments(
    jd_tt: float,
) -> Tuple[float, float, float, float]:
    """
    Calculate the fundamental arguments for lunar perturbation theory.

    These are the core angular arguments used in lunar perturbation series
    (Meeus "Astronomical Algorithms", Chapter 47).

    Extended with T^4 and T^5 terms from Simon et al. (1994) and
    Chapront et al. (2002) for improved accuracy at historical dates.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (D, M, M_prime, F) in radians, where:
            - D: Mean elongation of Moon from Sun
            - M: Mean anomaly of the Sun (solar perturbation)
            - M_prime: Mean anomaly of the Moon
            - F: Mean argument of latitude of the Moon

    References:
        - Meeus, J. "Astronomical Algorithms" (2nd ed., 1998), Chapter 47
        - Chapront-Touzé, M. & Chapront, J. "ELP 2000-85"
        - Simon, J.L. et al. (1994) "Numerical expressions for precession
          formulae and mean elements for the Moon and planets", A&A 282
        - Chapront, J. et al. (2002) "A new determination of lunar orbital
          parameters, precession constant and tidal acceleration", A&A 387
    """
    T = (jd_tt - 2451545.0) / 36525.0  # Julian centuries from J2000.0
    T2 = T * T
    T3 = T2 * T
    T4 = T3 * T
    T5 = T4 * T

    # Mean elongation of Moon from Sun (D)
    # Extended with T^5 term from Chapront et al. (2002)
    D = (
        297.8501921
        + 445267.1114034 * T
        - 0.0018819 * T2
        + T3 / 545868.0
        - T4 / 113065000.0
        + T5 / 18999000000.0  # T^5 correction for historical dates
    )

    # Mean anomaly of Sun (M) - solar perturbation argument
    # Extended with T^4 and T^5 terms from Simon et al. (1994)
    M = (
        357.5291092
        + 35999.0502909 * T
        - 0.0001536 * T2
        + T3 / 24490000.0
        - T4 / 992300000.0  # T^4 term for improved historical accuracy
        + T5 / 189900000000.0  # T^5 correction
    )

    # Mean anomaly of Moon (M')
    # Extended with T^5 term from Chapront et al. (2002)
    M_prime = (
        134.9633964
        + 477198.8675055 * T
        + 0.0087414 * T2
        + T3 / 69699.0
        - T4 / 14712000.0
        + T5 / 2520410000.0  # T^5 correction for historical dates
    )

    # Mean argument of latitude of Moon (F)
    # Extended with T^5 term from Chapront et al. (2002)
    F = (
        93.2720950
        + 483202.0175233 * T
        - 0.0036539 * T2
        - T3 / 3526000.0
        + T4 / 863310000.0
        - T5 / 142650000000.0  # T^5 correction for historical dates
    )

    # Convert to radians and normalize to [0, 2π)
    D = math.radians(D % 360.0)
    M = math.radians(M % 360.0)
    M_prime = math.radians(M_prime % 360.0)
    F = math.radians(F % 360.0)

    return D, M, M_prime, F


def _calc_de440_apogee_passage_terms(jd_tt: float) -> float:
    """Evaluate the IERS-basis terms fitted to DE440 apogee passages.

    The project-designed harmonic basis and every coefficient are generated by
    ``scripts/generate_lunar_apse_model.py``. The input angles are the ERFA
    implementations of the IERS 2003 Delaunay arguments. Trig-only residuals
    against the 14,581 DE440 passages have RMS 163 arcsec; the separately
    generated residual table supplies the remaining low-frequency correction.
    """
    T = (jd_tt - 2451545.0) / 36525.0  # Julian centuries from J2000.0
    D, M, M_prime, F = lunar_delaunay_arguments(jd_tt)
    eccentricity = (
        EARTH_ECCENTRICITY_J2000
        + EARTH_ECCENTRICITY_T * T
        + EARTH_ECCENTRICITY_T2 * T**2
    )
    E = eccentricity / EARTH_ECCENTRICITY_J2000

    perturbation = 0.0
    for coeff, c_d, c_m, c_mp, c_f, e_pow in APOGEE_TRIG_TERMS:
        term = coeff * math.sin(c_d * D + c_m * M + c_mp * M_prime + c_f * F)
        if e_pow == 1:
            term *= E
        elif e_pow == 2:
            term *= E * E
        perturbation += term
    return perturbation


def _calc_de440_perigee_passage_terms(jd_tt: float) -> float:
    """Evaluate the IERS-basis terms fitted to DE440 perigee passages.

    The perigee basis is generated independently of the apogee basis and keeps
    harmonics through order 30. Every coefficient comes from the 14,581 DE440
    passages. Trig-only residuals have RMS 403 arcsec; the generated residual
    table supplies the remaining low-frequency correction.
    """
    T = (jd_tt - 2451545.0) / 36525.0  # Julian centuries from J2000.0
    D, M, M_prime, F = lunar_delaunay_arguments(jd_tt)
    eccentricity = (
        EARTH_ECCENTRICITY_J2000
        + EARTH_ECCENTRICITY_T * T
        + EARTH_ECCENTRICITY_T2 * T**2
    )
    E = eccentricity / EARTH_ECCENTRICITY_J2000

    perturbation = 0.0
    for coeff, c_d, c_m, c_mp, c_f, e_pow in PERIGEE_TRIG_TERMS:
        term = coeff * math.sin(c_d * D + c_m * M + c_mp * M_prime + c_f * F)
        if e_pow == 1:
            term *= E
        elif e_pow == 2:
            term *= E * E
        perturbation += term
    return perturbation


def _mean_obliquity_radians(jd_tt: float) -> float:
    """
    Calculate mean obliquity of the ecliptic (IAU 2006).

    The obliquity of the ecliptic changes over time due to precession.
    Using a fixed J2000 value introduces an error that grows with distance
    from J2000.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        float: Mean obliquity in radians

    Formula (IAU 2006):
        ε = 84381.406" - 46.836769"×T - 0.0001831"×T² + 0.00200340"×T³
            - 0.000000576"×T⁴ - 0.0000000434"×T⁵
        where T = Julian centuries since J2000.0

    References:
        - Capitaine, N., Wallace, P. T. & Chapront, J. (2003), "Expressions for
          IAU 2000 precession quantities", Astronomy & Astrophysics 412, 567-586,
          eq. 37 (the P03 obliquity polynomial adopted as IAU 2006).
        - CALCS.md, section "Conversione coordinate ICRS → Eclittiche"
    """
    T = (jd_tt - 2451545.0) / 36525.0  # Julian centuries from J2000.0

    # Mean obliquity in arcseconds (IAU 2006)
    eps_arcsec = (
        84381.406
        - 46.836769 * T
        - 0.0001831 * T**2
        + 0.00200340 * T**3
        - 0.000000576 * T**4
        - 0.0000000434 * T**5
    )

    # Convert arcseconds to radians: arcsec -> degrees -> radians
    return math.radians(eps_arcsec / 3600.0)


class MeeusPolynomialWarning(UserWarning):
    """Warning issued when Meeus polynomial is used outside its optimal range.

    The warning thresholds are based on distance from J2000.0. They never
    raise a range error: beyond ±1000 years the message reports degraded
    precision, and beyond ±2000 years it warns that the error may exceed one
    degree. These are project warning bands around the retained polynomial,
    not universal accuracy guarantees for every lunar quantity.
    """

    pass


class MeeusRangeError(ValueError):
    """Deprecated: this exception is defined but never raised.

    Kept for backward compatibility with code that may import it.
    Use ``EphemerisRangeError`` for date-range errors instead.
    """

    pass


def _check_meeus_range(T: float) -> None:
    """Check if T (centuries from J2000) is within valid/optimal range and emit warnings.

    Args:
        T: Centuries from J2000.0

    Emits:
        MeeusPolynomialWarning if outside optimal (±10) or valid (±20) range
    """
    abs_T = abs(T)
    if abs_T > 20.0:
        warnings.warn(
            f"Date is {abs_T:.1f} centuries from J2000.0, outside the valid range "
            f"(±20 centuries). Error may exceed 1 degree. "
            f"Results should be used with caution.",
            MeeusPolynomialWarning,
            stacklevel=3,
        )
    elif abs_T > 10.0:
        warnings.warn(
            f"Date is {abs_T:.1f} centuries from J2000.0, outside the optimal range "
            f"(±10 centuries). Precision is degraded but still usable.",
            MeeusPolynomialWarning,
            stacklevel=3,
        )


def calc_true_lunar_node(jd_tt: float) -> Tuple[float, float, float]:
    """
        Calculate True (osculating) Lunar Node using orbital mechanics approach.

        The True Lunar Node represents the instantaneous ascending node of the Moon's
        osculating orbit - the point where the Moon crosses the ecliptic plane from
        south to north at the given moment. Unlike the Mean Node (which moves smoothly
        at ~19.3°/year retrograde), the True Node oscillates around the mean position
        with amplitudes up to ±1.5° on timescales of days to weeks.

        **Calculation method**

        This function uses a rigorous orbital mechanics approach, computing the
        angular momentum vector directly in the true ecliptic frame of date:

        **Step 1: Obtain Moon State Vectors in Ecliptic Frame**
            - Query JPL DE ephemeris (DE440/DE441) via Skyfield
            - Get geocentric position r and velocity v in the true ecliptic
              frame of date (Skyfield's ``ecliptic_frame``)
            - This frame automatically includes IAU 2006 precession and
              IAU 2000A nutation via Skyfield's internal rotation matrices

        **Step 2: Compute Angular Momentum Vector**
            - h = r × v (cross product) in ecliptic coordinates
            - h is perpendicular to the instantaneous orbital plane
            - Since r and v are already in the ecliptic frame, no further
              coordinate transformation is needed

        **Step 3: Compute Ascending Node Longitude**
            - The ascending node direction n = k × h (where k is ecliptic pole)
            - Simplifies to: n = (-h_y, h_x, 0)
            - Longitude = atan2(h_x, -h_y), normalized to [0°, 360°)
            - Result is directly in the true ecliptic of date

        **Mathematical foundation**

        The osculating orbital plane is defined by the angular momentum vector:

            h = r × v = norm(r) norm(v) sin(θ) n̂

        where θ is the angle between r and v, and n̂ is the unit normal to the
        orbital plane. The ascending node is where the orbital plane intersects
        the ecliptic plane, moving from south to north.

        For the ascending node direction:
            n_node = k̂_ecliptic × ĥ

        The longitude of the ascending node Ω is:
            Ω = atan2(n_x, n_y) = atan2(h_x, -h_y)

        This geometric approach captures the instantaneous orbital geometry,
        including all perturbations affecting the Moon's position and velocity
        at the given moment.

        Args:
            jd_tt: Julian Day in Terrestrial Time (TT).
                   TT is the uniform time scale used for ephemeris calculations,
                   approximately TT = UTC + 32.184 seconds + leap seconds.

        Returns:
            Tuple[float, float, float]: (longitude, latitude, distance) where:
                - longitude: Ecliptic longitude of ascending node in degrees [0, 360),
                            referenced to true ecliptic of date (includes nutation)
                - latitude: Always 0.0 (the node lies on the ecliptic by definition)
                - distance: Geocentric distance of the osculating orbit at the
                           ascending node in AU, from the conic r = p/(1 + e cos nu)
                           evaluated at the node (p, e, omega from r, v)

    **Precision and accuracy**

    **Compared to JPL DE ephemeris geometric method (1000 random dates, 1950-2050):**
        - Mean error: ~8.9 arcsec (~0.0025 degrees)
        - RMS error: ~11.8 arcsec (~0.0033 degrees)
        - Maximum error: ~52 arcsec (~0.014 degrees)
        - 100% of dates within 60 arcsec

    **Across the full DE440 range (1550-2650):**
        - Typical error: 2-13 arcsec
        - Maximum observed: ~23 arcsec

    **Temporal Behavior:**
        - The true node oscillates +/-1.5 degrees around the mean node
        - Primary oscillation period: ~27.2 days (draconic month)
        - Secondary oscillations: fortnightly (~14.8 days), monthly (~29.5 days)
        - Long-term motion: retrograde ~19.3 degrees/year (18.6 year period)

        **Physical interpretation**

        The True Node represents the actual intersection of the Moon's instantaneous
        orbit with the ecliptic. Key points:

        - **Eclipses**: Solar and lunar eclipses occur when the Sun is near a node
          during New/Full Moon. The True Node gives the instantaneous position.

        - **Oscillation**: The ±1.5° oscillation is primarily caused by:
          1. The Moon's orbital eccentricity (e ≈ 0.0549)
          2. Solar gravitational perturbations (evection, variation)
          3. The tilt of the Moon's orbit (~5.145° to ecliptic)

        - **Astrological Use**: Many systems prefer the True Node for its
          astronomical accuracy; others use the Mean Node for smoother motion.

        Note:
            The true node can move rapidly (several arcminutes per hour) and
            occasionally appears to reverse direction briefly due to the
            complex perturbation interplay. This is physically correct behavior,
            not a calculation artifact.

        See Also:
            - calc_mean_lunar_node: Smoothed average node position

    References:
        Primary:
            - Park, R.S. et al. (2021) "The JPL Planetary and Lunar Ephemerides DE440 and DE441"
            - Vallado, D. "Fundamentals of Astrodynamics and Applications"
              (4th ed., 2013), Chapter 2: Orbit Determination

        Orbital Mechanics:
            - Bate, Mueller, White "Fundamentals of Astrodynamics" (1971)
            - Roy, A.E. "Orbital Motion" (4th ed., 2005)
    """
    from skyfield.framelib import ecliptic_frame

    from .cache import get_cached_time_tt

    planets = _get_computation_ephemeris()
    t = get_cached_time_tt(jd_tt)

    earth = planets["earth"]
    moon = planets["moon"]

    # Get geocentric Moon position and velocity directly in the true ecliptic
    # frame of date. Skyfield's ecliptic_frame applies IAU 2006 precession
    # and IAU 2000A nutation internally, so no manual corrections are needed.
    moon_pos = (moon - earth).at(t)
    r = moon_pos.frame_xyz(ecliptic_frame).au
    v = moon_pos.frame_xyz_and_velocity(ecliptic_frame)[1].au_per_d

    # Angular momentum vector h = r × v (perpendicular to orbital plane)
    # Since r and v are already in the ecliptic frame, h components give
    # us the orbital plane orientation relative to the ecliptic directly.
    h_x = float(r[1] * v[2] - r[2] * v[1])
    h_y = float(r[2] * v[0] - r[0] * v[2])
    h_z = float(r[0] * v[1] - r[1] * v[0])

    # The ascending node longitude in the true ecliptic of date
    # n = k × h = (-h_y, h_x, 0), longitude = atan2(h_x, -h_y)
    node_lon = math.degrees(math.atan2(h_x, -h_y)) % 360.0

    # No analytical perturbation series is applied here: the geometric h = r x v
    # approach already captures perturbation effects through the JPL DE
    # ephemeris state vectors.

    # Distance at the ascending node from osculating orbit elements.
    # Uses the vis-viva relation to derive semi-latus rectum p = h²/GM,
    # eccentricity from e_vec = (v×h)/GM - r̂, and argument of perigee ω.
    # Distance at ascending node: r = p / (1 + e·cos(2π - ω))
    _GM_EARTH_MOON = GM_EARTH_MOON_AU3_DAY2

    r_x, r_y, r_z = float(r[0]), float(r[1]), float(r[2])
    v_x, v_y, v_z = float(v[0]), float(v[1]), float(v[2])
    r_mag = math.sqrt(r_x * r_x + r_y * r_y + r_z * r_z)
    h_mag = math.sqrt(h_x * h_x + h_y * h_y + h_z * h_z)

    # Semi-latus rectum
    p_orbit = h_mag * h_mag / _GM_EARTH_MOON

    # Eccentricity vector: e = (v × h) / GM - r̂
    # v × h components
    vxh_x = v_y * h_z - v_z * h_y
    vxh_y = v_z * h_x - v_x * h_z
    vxh_z = v_x * h_y - v_y * h_x
    e_x = vxh_x / _GM_EARTH_MOON - r_x / r_mag
    e_y = vxh_y / _GM_EARTH_MOON - r_y / r_mag
    e_z = vxh_z / _GM_EARTH_MOON - r_z / r_mag
    ecc = math.sqrt(e_x * e_x + e_y * e_y + e_z * e_z)

    # Argument of perigee: angle from ascending node direction to eccentricity vector
    # Node direction: n = (-h_y, h_x, 0)
    n_mag = math.sqrt(h_y * h_y + h_x * h_x)
    if n_mag > 0 and ecc > 1e-15:
        cos_omega = (-h_y * e_x + h_x * e_y) / (n_mag * ecc)
        cos_omega = max(-1.0, min(1.0, cos_omega))
        omega = math.acos(cos_omega)
        if e_z < 0:
            omega = 2.0 * math.pi - omega

        # True anomaly at ascending node: f = 2π - ω
        f_asc = 2.0 * math.pi - omega
        dist = p_orbit / (1.0 + ecc * math.cos(f_asc))
    else:
        # Degenerate case: use mean distance
        dist = 384400.0 / 149597870.7

    return node_lon, 0.0, dist


def calc_true_lilith(jd_tt: float) -> Tuple[float, float, float]:
    """
    Calculate True Lilith (osculating lunar apogee).

    Computes the osculating lunar apogee using the eccentricity vector method
    in the true ecliptic frame of date. The eccentricity vector, derived from
    the Moon's instantaneous position and velocity from JPL DE ephemeris,
    points toward perigee. The apogee direction is 180° from perigee.

    **Algorithm**

    **Step 1: Obtain Moon State Vectors in Ecliptic Frame**
        - Query JPL DE ephemeris via Skyfield
        - Get geocentric position r and velocity v in the true ecliptic
          frame of date (Skyfield's ``ecliptic_frame``)
        - This frame automatically includes IAU 2006 precession and
          IAU 2000A nutation

    **Step 2: Compute Eccentricity Vector**
        - h = r × v (angular momentum)
        - e = (v × h)/μ - r/norm(r) (points toward perigee)
        - μ = G(M_Earth + M_Moon) for the two-body problem
        - Apogee direction = -e (opposite to perigee)

    **Step 3: Convert to Ecliptic Coordinates**
        - Since r and v are already in the ecliptic frame, the apogee
          vector is directly in ecliptic coordinates of date
        - Convert from Cartesian to spherical (longitude, latitude)

    **Physical background**

    The osculating lunar apogee is the apogee direction of the instantaneous
    Keplerian orbit that passes through the Moon's current position with its
    current velocity. This differs from the mean apogee due to solar and
    planetary perturbations that continuously modify the Moon's orbital shape.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT).

    Returns:
        Tuple[float, float, float]: (longitude, latitude, distance) where:
            - longitude: Ecliptic longitude of apogee in degrees [0, 360)
            - latitude: Ecliptic latitude in degrees (small, typically < 5°)
            - distance: Apogee distance from Earth in AU

    **Precision**

    **Computed from JPL DE ephemeris state vectors (500 random dates, 1950-2050):**
        - Mean internal consistency: sub-arcsecond
        - Maximum numerical error: ~1 milliarcsecond

    References:
        - Vallado, D. "Fundamentals of Astrodynamics and Applications"
        - Park, R.S. et al. (2021) "The JPL Planetary and Lunar Ephemerides DE440 and DE441"
    """
    from skyfield.framelib import ecliptic_frame

    from .cache import get_cached_time_tt

    planets = _get_computation_ephemeris()
    t = get_cached_time_tt(jd_tt)

    earth = planets["earth"]
    moon = planets["moon"]

    # Get geocentric Moon position and velocity directly in the true ecliptic
    # frame of date. Skyfield's ecliptic_frame applies IAU 2006 precession
    # and IAU 2000A nutation internally, so no manual corrections are needed.
    moon_pos = (moon - earth).at(t)
    r = moon_pos.frame_xyz(ecliptic_frame).au
    v = moon_pos.frame_xyz_and_velocity(ecliptic_frame)[1].au_per_d

    r_mag = math.sqrt(float(r[0] ** 2 + r[1] ** 2 + r[2] ** 2))

    # Angular momentum h = r × v
    h_x = r[1] * v[2] - r[2] * v[1]
    h_y = r[2] * v[0] - r[0] * v[2]
    h_z = r[0] * v[1] - r[1] * v[0]

    # Earth-Moon system gravitational parameter (shared module constant)
    mu = GM_EARTH_MOON_AU3_DAY2

    # Eccentricity vector e = (v × h)/μ - r/|r| (points toward perigee)
    vxh_x = v[1] * h_z - v[2] * h_y
    vxh_y = v[2] * h_x - v[0] * h_z
    vxh_z = v[0] * h_y - v[1] * h_x

    e_x = float(vxh_x / mu - r[0] / r_mag)
    e_y = float(vxh_y / mu - r[1] / r_mag)
    e_z = float(vxh_z / mu - r[2] / r_mag)

    e_mag = math.sqrt(e_x**2 + e_y**2 + e_z**2)

    # Apogee is opposite to perigee (180° from eccentricity vector)
    apogee_x, apogee_y, apogee_z = -e_x, -e_y, -e_z
    apogee_mag = math.sqrt(apogee_x**2 + apogee_y**2 + apogee_z**2)

    # Convert to spherical ecliptic coordinates (already in ecliptic of date)
    longitude = math.degrees(math.atan2(apogee_y, apogee_x)) % 360.0
    lat = math.degrees(math.asin(apogee_z / apogee_mag))

    # Compute apogee distance in AU from orbital elements.
    # OSCU_APOG returns Earth-to-apogee distance, not eccentricity magnitude.
    # Semi-latus rectum: p = h²/μ
    # Apogee distance: r_apogee = p / (1 - e)  [= a(1+e)]
    h_mag = math.sqrt(float(h_x**2 + h_y**2 + h_z**2))
    p = h_mag**2 / mu
    r_apogee = p / (1.0 - e_mag)

    return longitude, lat, r_apogee


def calc_true_lilith_orbital_elements(jd_tt: float) -> Tuple[float, float, float]:
    """
    Calculate True Lilith using the classical orbital elements method.

    This is an alias for calc_true_lilith() maintained for backward
    compatibility. Both functions now use the same ecliptic-frame
    eccentricity vector approach.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT).

    Returns:
        Tuple[float, float, float]: (longitude, latitude, distance) where:
            - longitude: Ecliptic longitude of apogee in degrees [0, 360)
            - latitude: Ecliptic latitude in degrees (small, typically < 5°)
            - distance: Apogee distance from Earth in AU

    See Also:
        calc_true_lilith: Primary implementation.
    """
    return calc_true_lilith(jd_tt)


def calc_osculating_perigee(jd_tt: float) -> Tuple[float, float, float]:
    """
    Calculate the osculating lunar perigee directly from orbital mechanics.

    Computes the osculating lunar perigee using the eccentricity vector method
    in the true ecliptic frame of date. The eccentricity vector, derived from
    the Moon's instantaneous position and velocity from JPL DE ephemeris,
    naturally points toward perigee.

    Physical Background
    ===================

    The osculating perigee is calculated from instantaneous orbital elements
    that change rapidly due to solar perturbations. The perigee can deviate
    significantly from being exactly opposite to apogee depending on Sun-Moon
    geometry.

    This function computes perigee independently from apogee by using the
    eccentricity vector directly (which points toward perigee by definition)
    rather than negating the apogee direction.

    Algorithm
    =========

    **Step 1: Obtain Moon State Vectors in Ecliptic Frame**
        - Query JPL DE ephemeris via Skyfield
        - Get geocentric position r and velocity v in the true ecliptic
          frame of date (Skyfield's ``ecliptic_frame``)
        - This frame automatically includes IAU 2006 precession and
          IAU 2000A nutation

    **Step 2: Compute Eccentricity Vector**
        - h = r x v (angular momentum)
        - e = (v x h)/mu - r/|r| (points toward perigee)
        - mu = G(M_Earth + M_Moon) for the two-body problem
        - Perigee direction = +e (the eccentricity vector itself)

    **Step 3: Convert to Ecliptic Coordinates**
        - Since r and v are already in the ecliptic frame, the perigee
          vector is directly in ecliptic coordinates of date
        - Convert from Cartesian to spherical (longitude, latitude)

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT).

    Returns:
        Tuple[float, float, float]: (longitude, latitude, distance) where:
            - longitude: Ecliptic longitude of perigee in degrees [0, 360)
            - latitude: Ecliptic latitude in degrees (small, typically < 5 degrees)
            - distance: Perigee distance from Earth in AU

    References:
        - Vallado, D. "Fundamentals of Astrodynamics and Applications"
        - Park, R.S. et al. (2021) "The JPL Planetary and Lunar Ephemerides DE440 and DE441"
    """
    from skyfield.framelib import ecliptic_frame

    from .cache import get_cached_time_tt

    planets = _get_computation_ephemeris()
    t = get_cached_time_tt(jd_tt)

    earth = planets["earth"]
    moon = planets["moon"]

    # Get geocentric Moon position and velocity directly in the true ecliptic
    # frame of date. Skyfield's ecliptic_frame applies IAU 2006 precession
    # and IAU 2000A nutation internally, so no manual corrections are needed.
    moon_pos = (moon - earth).at(t)
    r = moon_pos.frame_xyz(ecliptic_frame).au
    v = moon_pos.frame_xyz_and_velocity(ecliptic_frame)[1].au_per_d

    r_mag = math.sqrt(float(r[0] ** 2 + r[1] ** 2 + r[2] ** 2))

    # Angular momentum h = r × v
    h_x = r[1] * v[2] - r[2] * v[1]
    h_y = r[2] * v[0] - r[0] * v[2]
    h_z = r[0] * v[1] - r[1] * v[0]

    # Earth-Moon system gravitational parameter (shared module constant)
    mu = GM_EARTH_MOON_AU3_DAY2

    # Eccentricity vector e = (v × h)/μ - r/|r| (points toward perigee)
    vxh_x = v[1] * h_z - v[2] * h_y
    vxh_y = v[2] * h_x - v[0] * h_z
    vxh_z = v[0] * h_y - v[1] * h_x

    e_x = float(vxh_x / mu - r[0] / r_mag)
    e_y = float(vxh_y / mu - r[1] / r_mag)
    e_z = float(vxh_z / mu - r[2] / r_mag)

    e_mag = math.sqrt(e_x**2 + e_y**2 + e_z**2)

    # Perigee is in the direction of the eccentricity vector (no negation)
    # This is different from calc_true_lilith which negates to get apogee
    perigee_x, perigee_y, perigee_z = e_x, e_y, e_z

    # Convert to spherical ecliptic coordinates (already in ecliptic of date)
    longitude = math.degrees(math.atan2(perigee_y, perigee_x)) % 360.0
    lat = math.degrees(math.asin(perigee_z / e_mag))

    # Compute perigee distance in AU from orbital elements.
    # Semi-latus rectum: p = h²/μ
    # Perigee distance: r_perigee = p / (1 + e)  [= a(1-e)]
    h_mag = math.sqrt(float(h_x**2 + h_y**2 + h_z**2))
    p = h_mag**2 / mu
    r_perigee = p / (1.0 + e_mag)

    return longitude, lat, r_perigee


def _get_ephemeris_range() -> Tuple[float, float]:
    """
    Get the valid Julian Date range for the current ephemeris.

    Returns:
        Tuple[float, float]: (min_jd, max_jd) Julian Date range in TT.

    Note:
        The range depends on the loaded ephemeris file:
        - de440s.bsp: 1849 to 2150
        - de440.bsp (default): 1550 to 2650
        - de441.bsp: -13200 to 17191 (split into two segments)

        In sealed calculation mode ``leb`` the range comes from the active
        LEB reader (Moon/Earth channel intersection) instead of opening the
        JPL kernel, mirroring ``_get_computation_ephemeris``.
    """
    from .state import _JPL_SOURCE_ACCESS, get_calc_mode, get_leb_reader

    # _JPL_SOURCE_ACCESS is the sanctioned provisioning window of
    # state._allow_jpl_source(): only offline LEB GENERATORS enter it (their
    # purpose is reading the registered source kernels), never a runtime
    # calculation - so this gate cannot open JPL under sealed runtime use.
    # With no reader/coverage the generic fallback below is a conservative
    # SUBSET of every tier's real coverage (it under-promises; an actual
    # out-of-coverage computation still raises its typed error later).
    if get_calc_mode() == "leb" and not _JPL_SOURCE_ACCESS.get():
        try:
            reader = get_leb_reader()
        except RuntimeError:
            reader = None
        if reader is not None:
            from .inventory import get_reader_body_coverage

            moon_cov = get_reader_body_coverage(reader, 1)
            earth_cov = get_reader_body_coverage(reader, 14)
            if moon_cov is not None and earth_cov is not None:
                return (
                    max(moon_cov.jd_start, earth_cov.jd_start),
                    min(moon_cov.jd_end, earth_cov.jd_end),
                )
        return (2415020.0, 2471184.0)

    planets = get_planets()
    ts = get_timescale()

    try:
        min_jd = float("inf")
        max_jd = float("-inf")

        for segment in planets.segments:
            try:
                start_time, end_time = segment.time_range(ts)
                min_jd = min(min_jd, float(start_time.tt))
                max_jd = max(max_jd, float(end_time.tt))
            except (AttributeError, TypeError, ValueError):
                continue

        if min_jd == float("inf"):
            return (2415020.0, 2471184.0)

        return (min_jd, max_jd)
    except (AttributeError, TypeError, ValueError, KeyError):
        return (2415020.0, 2471184.0)


def _sample_osculating_apogee_with_fallback(
    jd_tt: float,
    half_window_days: float,
    num_samples: int,
) -> Tuple[list, list, list, list, int]:
    """
    Sample osculating apogee positions with fallback for ephemeris boundaries.

    This function handles edge cases where the sampling window extends beyond
    the ephemeris range by:
    1. First trying symmetric sampling around the target date
    2. If that fails, shifting the window to stay within the ephemeris range
    3. If asymmetric sampling still fails, reducing the number of samples
    4. As a last resort, returning just the central sample

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT) - the target date.
        half_window_days: Half of the total sampling window in days.
        num_samples: Number of samples to take.

    Returns:
        Tuple containing:
            - sample_times: List of Julian dates for successful samples
            - sample_lons: List of longitudes in degrees
            - sample_lats: List of latitudes in degrees
            - sample_eccs: List of eccentricities
            - target_idx: Index of the sample closest to jd_tt
    """
    # Get ephemeris range
    min_jd, max_jd = _get_ephemeris_range()

    # Check if symmetric window is within range
    if jd_tt - half_window_days >= min_jd and jd_tt + half_window_days <= max_jd:
        # Standard symmetric sampling
        sample_times = []
        for i in range(num_samples):
            offset = -half_window_days + (2 * half_window_days * i / (num_samples - 1))
            sample_times.append(jd_tt + offset)
        target_idx = num_samples // 2
    else:
        # Need to adjust window for boundary
        if jd_tt - half_window_days < min_jd:
            # Near start of ephemeris - shift window forward
            window_start = max(min_jd, jd_tt - half_window_days)
            window_end = min(max_jd, window_start + 2 * half_window_days)
            # Ensure we don't exceed the end
            if window_end > max_jd:
                window_end = max_jd
                window_start = max(min_jd, window_end - 2 * half_window_days)
        else:
            # Near end of ephemeris - shift window backward
            window_end = min(max_jd, jd_tt + half_window_days)
            window_start = max(min_jd, window_end - 2 * half_window_days)
            # Ensure we don't go before the start
            if window_start < min_jd:
                window_start = min_jd
                window_end = min(max_jd, window_start + 2 * half_window_days)

        # Calculate actual window size
        actual_window = window_end - window_start

        # If window is very small, reduce samples proportionally
        if actual_window < half_window_days:
            # Very constrained window - use fewer samples
            num_samples = max(
                3, int(num_samples * actual_window / (2 * half_window_days))
            )

        # Generate sample times within the adjusted window
        sample_times = []
        if num_samples > 1:
            for i in range(num_samples):
                t = window_start + (actual_window * i / (num_samples - 1))
                sample_times.append(t)
        else:
            sample_times = [jd_tt]

        # Find which sample is closest to the target date
        target_idx = 0
        min_dist = abs(sample_times[0] - jd_tt)
        for i in range(1, len(sample_times)):
            dist = abs(sample_times[i] - jd_tt)
            if dist < min_dist:
                min_dist = dist
                target_idx = i

    # Sample the osculating apogee at each time, with fallback for failures
    sample_lons = []
    sample_lats = []
    sample_eccs = []
    valid_times = []

    for sample_jd in sample_times:
        try:
            lon, lat, ecc = calc_true_lilith(sample_jd)
            sample_lons.append(lon)
            sample_lats.append(lat)
            sample_eccs.append(ecc)
            valid_times.append(sample_jd)
        except (EphemerisRangeError, ValueError, ZeroDivisionError):
            # Skip samples that fail (outside ephemeris range)
            continue

    # If we lost samples, recalculate target_idx
    if len(valid_times) < len(sample_times):
        # Find which valid sample is closest to the target date
        if valid_times:
            target_idx = 0
            min_dist = abs(valid_times[0] - jd_tt)
            for i in range(1, len(valid_times)):
                dist = abs(valid_times[i] - jd_tt)
                if dist < min_dist:
                    min_dist = dist
                    target_idx = i
        else:
            # No valid samples - try just the target date
            try:
                lon, lat, ecc = calc_true_lilith(jd_tt)
                return [jd_tt], [lon], [lat], [ecc], 0
            except (EphemerisRangeError, ValueError, ZeroDivisionError):
                # Even target date fails - re-raise the original error
                raise

    # Need at least 2 samples for linear regression
    if len(valid_times) < 2:
        # Fall back to osculating apogee for the target date
        lon, lat, ecc = calc_true_lilith(jd_tt)
        return [jd_tt], [lon], [lat], [ecc], 0

    return valid_times, sample_lons, sample_lats, sample_eccs, target_idx


def _unwrap_longitudes(longitudes: list) -> list:
    """
    Unwrap a sequence of longitudes to handle 0°/360° discontinuity.

    This ensures polynomial fitting works correctly when the apogee crosses
    the 0°/360° boundary. The unwrapped values may exceed [0, 360) but
    maintain continuous change between consecutive values.

    Args:
        longitudes: List of longitude values in degrees [0, 360).

    Returns:
        List of unwrapped longitude values (may be outside [0, 360)).
    """
    if not longitudes:
        return []

    unwrapped = [longitudes[0]]
    for i in range(1, len(longitudes)):
        diff = longitudes[i] - longitudes[i - 1]
        # Detect wraparound: if diff > 180, we crossed from ~360 to ~0
        # if diff < -180, we crossed from ~0 to ~360
        if diff > 180:
            diff -= 360
        elif diff < -180:
            diff += 360
        unwrapped.append(unwrapped[-1] + diff)

    return unwrapped


def calc_interpolated_apogee(jd_tt: float) -> Tuple[float, float, float]:
    """Return the independently generated interpolated lunar apogee.

    Longitude combines the ERFA/IERS mean-apogee baseline, an IERS Delaunay
    series fitted to all DE440 apogee passages over the medium-kernel interval,
    and a fixed-grid interpolation of the remaining passage residual. Latitude
    is a two-harmonic inclined-plane fit to the same passages; distance is their
    mean geocentric distance.

    The result begins on the true ecliptic of date. Beyond the fitted
    approximately 1550--2650 interval, the residual table tapers to zero over
    one year and the trig-only longitude has about 0.045 degree RMS error
    against the DE440 passage track. This abstract coordinate is not a
    lunar-distance event finder.

    Args:
        jd_tt: Julian day in Terrestrial Time.

    Returns:
        Native Python floats ``(longitude_deg, latitude_deg, distance_au)``.

    References:
        IERS Conventions (2010), Eq. 5.43; Park et al. (2021), AJ 161:105.
    """
    # Standards-derived IERS mean apogee used by the DE440 model generator.
    mean_apogee = _mean_lunar_apogee_position_unchecked(jd_tt)[0]

    # Add the passage-fitted terms that model the apsidal oscillation.
    perturbation = _calc_de440_apogee_passage_terms(jd_tt)

    # Combine mean position and perturbations
    interp_lon = mean_apogee + perturbation

    # Add correction table residual (arcseconds -> degrees) for high precision
    if _APSE_CORRECTIONS_AVAILABLE and APOGEE_CORRECTIONS:
        correction_arcsec = _interpolate_apse_correction(
            jd_tt,
            APOGEE_CT_JD_START,
            APOGEE_CT_STEP_DAYS,
            APOGEE_CORRECTIONS,
            APOGEE_CT_COUNT,
        )
        interp_lon += correction_arcsec / 3600.0

    interp_lon = interp_lon % 360.0

    # Latitude: inclination * sin(longitude - ascending node)
    # The interpolated apogee lies near the lunar orbital plane, which is
    # inclined ~5.145° to the ecliptic. Latitude varies sinusoidally as the
    # apogee moves relative to the ascending node.
    node_lon = _mean_lunar_node_position_unchecked(jd_tt)[0]
    _a1, _a3 = APOGEE_LAT_COEFFS
    _omega = math.radians(interp_lon - node_lon)
    interp_lat = _a1 * math.sin(_omega) + _a3 * math.sin(3.0 * _omega)

    # Distance: mean DE440 apogee-passage distance (~405,500 km).
    interp_dist = APOGEE_MEAN_DIST_AU

    return interp_lon, interp_lat, interp_dist


def _cubic_spline_interpolate(x: list, y: list, x_eval: float) -> float:
    """
    Natural cubic spline interpolation.

    Computes the interpolated value at x_eval using natural cubic spline
    through the given (x, y) data points. Natural spline has zero second
    derivative at the endpoints.

    The cubic spline ensures:
    1. The interpolant passes through all data points
    2. First and second derivatives are continuous at interior points
    3. Second derivative is zero at endpoints (natural boundary condition)

    Algorithm
    =========

    For n+1 data points (x_0, y_0), ..., (x_n, y_n), the cubic spline S(x)
    is defined piecewise on each interval [x_i, x_{i+1}]:

        S_i(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3

    where the coefficients are determined by:
    1. S_i(x_i) = y_i (interpolation condition)
    2. S_i(x_{i+1}) = S_{i+1}(x_{i+1}) (continuity)
    3. S_i'(x_{i+1}) = S_{i+1}'(x_{i+1}) (continuous first derivative)
    4. S_i''(x_{i+1}) = S_{i+1}''(x_{i+1}) (continuous second derivative)
    5. S_0''(x_0) = 0 and S_{n-1}''(x_n) = 0 (natural boundary conditions)

    Args:
        x: List of x-coordinates (must be sorted in ascending order)
        y: List of y-coordinates (same length as x)
        x_eval: The x-coordinate at which to evaluate the spline

    Returns:
        float: The interpolated y-value at x_eval

    References:
        - Press et al., "Numerical Recipes", Chapter 3.3
        - Burden & Faires, "Numerical Analysis", Chapter 3.5
    """
    n = len(x) - 1  # Number of intervals

    if n < 1:
        return y[0] if len(y) > 0 else 0.0

    # Compute the differences
    h = [x[i + 1] - x[i] for i in range(n)]

    # Set up the tridiagonal system for the second derivatives (c values)
    # For natural spline: c[0] = 0, c[n] = 0

    # Right-hand side of the tridiagonal system
    alpha = [0.0] * n
    for i in range(1, n):
        alpha[i] = (3.0 / h[i]) * (y[i + 1] - y[i]) - (3.0 / h[i - 1]) * (
            y[i] - y[i - 1]
        )

    # Solve the tridiagonal system using Thomas algorithm
    # Diagonal elements: 2*(h[i-1] + h[i]) for i = 1, ..., n-1
    # Sub/super diagonal: h[i]

    diag = [0.0] * (n + 1)
    mu = [0.0] * (n + 1)
    z = [0.0] * (n + 1)

    diag[0] = 1.0
    mu[0] = 0.0
    z[0] = 0.0

    for i in range(1, n):
        diag[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1]
        if abs(diag[i]) < 1e-15:
            diag[i] = 1e-15  # Prevent division by zero
        mu[i] = h[i] / diag[i]
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / diag[i]

    diag[n] = 1.0
    z[n] = 0.0

    # Back substitution to get c values (second derivatives / 2)
    c = [0.0] * (n + 1)
    c[n] = 0.0

    for j in range(n - 1, -1, -1):
        c[j] = z[j] - mu[j] * c[j + 1]

    # Compute b and d coefficients
    b = [0.0] * n
    d = [0.0] * n

    for i in range(n):
        b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (c[i + 1] + 2.0 * c[i]) / 3.0
        d[i] = (c[i + 1] - c[i]) / (3.0 * h[i])

    # Find the interval containing x_eval
    # Clamp to valid range if outside
    if x_eval <= x[0]:
        idx = 0
        dx = x_eval - x[0]
    elif x_eval >= x[n]:
        idx = n - 1
        dx = x_eval - x[n - 1]
    else:
        # Binary search for the interval
        idx = 0
        for i in range(n):
            if x[i] <= x_eval < x[i + 1]:
                idx = i
                break
        dx = x_eval - x[idx]

    # Evaluate the spline at x_eval
    # S_i(x) = y[i] + b[i]*dx + c[i]*dx^2 + d[i]*dx^3
    result = y[idx] + b[idx] * dx + c[idx] * dx * dx + d[idx] * dx * dx * dx

    return result


def _sample_osculating_perigee_with_fallback(
    jd_tt: float,
    half_window_days: float,
    num_samples: int,
) -> Tuple[list, list, list, list, int]:
    """
    Sample osculating perigee positions with fallback for ephemeris boundaries.

    This function handles edge cases where the sampling window extends beyond
    the ephemeris range by:
    1. First trying symmetric sampling around the target date
    2. If that fails, shifting the window to stay within the ephemeris range
    3. If asymmetric sampling still fails, reducing the number of samples
    4. As a last resort, returning just the central sample

    The perigee is computed directly from the eccentricity vector using
    calc_osculating_perigee(), NOT by adding 180 degrees to apogee. This is
    important because osculating apogee and perigee are NOT exactly opposite -
    they can deviate by up to 28 degrees depending on Sun-Moon geometry, and
    are only roughly opposite when the Sun is in conjunction with one of them
    or at a 90 degree angle.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT) - the target date.
        half_window_days: Half of the total sampling window in days.
        num_samples: Number of samples to take.

    Returns:
        Tuple containing:
            - sample_times: List of Julian dates for successful samples
            - sample_lons: List of perigee longitudes in degrees
            - sample_lats: List of perigee latitudes in degrees
            - sample_eccs: List of eccentricities
            - target_idx: Index of the sample closest to jd_tt
    """
    # Get ephemeris range
    min_jd, max_jd = _get_ephemeris_range()

    # Check if symmetric window is within range
    if jd_tt - half_window_days >= min_jd and jd_tt + half_window_days <= max_jd:
        # Standard symmetric sampling
        sample_times = []
        for i in range(num_samples):
            offset = -half_window_days + (2 * half_window_days * i / (num_samples - 1))
            sample_times.append(jd_tt + offset)
        target_idx = num_samples // 2
    else:
        # Need to adjust window for boundary
        if jd_tt - half_window_days < min_jd:
            # Near start of ephemeris - shift window forward
            window_start = max(min_jd, jd_tt - half_window_days)
            window_end = min(max_jd, window_start + 2 * half_window_days)
            # Ensure we don't exceed the end
            if window_end > max_jd:
                window_end = max_jd
                window_start = max(min_jd, window_end - 2 * half_window_days)
        else:
            # Near end of ephemeris - shift window backward
            window_end = min(max_jd, jd_tt + half_window_days)
            window_start = max(min_jd, window_end - 2 * half_window_days)
            # Ensure we don't go before the start
            if window_start < min_jd:
                window_start = min_jd
                window_end = min(max_jd, window_start + 2 * half_window_days)

        # Calculate actual window size
        actual_window = window_end - window_start

        # If window is very small, reduce samples proportionally
        if actual_window < half_window_days:
            # Very constrained window - use fewer samples
            num_samples = max(
                3, int(num_samples * actual_window / (2 * half_window_days))
            )

        # Generate sample times within the adjusted window
        sample_times = []
        if num_samples > 1:
            for i in range(num_samples):
                t = window_start + (actual_window * i / (num_samples - 1))
                sample_times.append(t)
        else:
            sample_times = [jd_tt]

        # Find which sample is closest to the target date
        target_idx = 0
        min_dist = abs(sample_times[0] - jd_tt)
        for i in range(1, len(sample_times)):
            dist = abs(sample_times[i] - jd_tt)
            if dist < min_dist:
                min_dist = dist
                target_idx = i

    # Sample the osculating perigee at each time, with fallback for failures
    # Use calc_osculating_perigee to compute perigee directly from the
    # eccentricity vector, rather than deriving it from apogee + 180 degrees.
    # This is important because apogee and perigee are NOT exactly 180 degrees
    # apart - they can deviate by up to 28 degrees depending on Sun-Moon geometry.
    sample_lons = []
    sample_lats = []
    sample_eccs = []
    valid_times = []

    for sample_jd in sample_times:
        try:
            perigee_lon, perigee_lat, ecc = calc_osculating_perigee(sample_jd)
            sample_lons.append(perigee_lon)
            sample_lats.append(perigee_lat)
            sample_eccs.append(ecc)
            valid_times.append(sample_jd)
        except (EphemerisRangeError, ValueError, ZeroDivisionError):
            # Skip samples that fail (outside ephemeris range)
            continue

    # If we lost samples, recalculate target_idx
    if len(valid_times) < len(sample_times):
        # Find which valid sample is closest to the target date
        if valid_times:
            target_idx = 0
            min_dist = abs(valid_times[0] - jd_tt)
            for i in range(1, len(valid_times)):
                dist = abs(valid_times[i] - jd_tt)
                if dist < min_dist:
                    min_dist = dist
                    target_idx = i
        else:
            # No valid samples - try just the target date
            try:
                perigee_lon, perigee_lat, ecc = calc_osculating_perigee(jd_tt)
                return [jd_tt], [perigee_lon], [perigee_lat], [ecc], 0
            except (EphemerisRangeError, ValueError, ZeroDivisionError):
                # Even target date fails - re-raise the original error
                raise

    # Need at least 2 samples for linear regression
    if len(valid_times) < 2:
        # Fall back to osculating perigee for the target date
        perigee_lon, perigee_lat, ecc = calc_osculating_perigee(jd_tt)
        return [jd_tt], [perigee_lon], [perigee_lat], [ecc], 0

    return valid_times, sample_lons, sample_lats, sample_eccs, target_idx


def calc_interpolated_perigee(jd_tt: float) -> Tuple[float, float, float]:
    """Return the independently generated interpolated lunar perigee.

    This curve is fitted independently from apogee and is not forced to be
    antipodal. Longitude combines the ERFA/IERS mean-perigee baseline, a
    Delaunay series fitted to all DE440 perigee passages over the medium-kernel
    interval, and a fixed-grid interpolation of the remaining passage residual.
    Latitude is a two-harmonic fit to the same passages; distance is their mean
    geocentric distance.

    The result begins on the true ecliptic of date. Beyond the fitted
    approximately 1550--2650 interval, the residual table tapers to zero over
    one year and the trig-only longitude has about 0.112 degree RMS error
    against the DE440 passage track. This abstract coordinate is not a
    lunar-distance or tidal event finder.

    Args:
        jd_tt: Julian day in Terrestrial Time.

    Returns:
        Native Python floats ``(longitude_deg, latitude_deg, distance_au)``.

    References:
        IERS Conventions (2010), Eq. 5.43; Park et al. (2021), AJ 161:105.
    """
    mean_perigee = (_mean_lunar_apogee_position_unchecked(jd_tt)[0] + 180.0) % 360.0

    perturbation = _calc_de440_perigee_passage_terms(jd_tt)

    interp_lon = mean_perigee + perturbation

    # Add correction table residual (arcseconds -> degrees) for high precision
    if _APSE_CORRECTIONS_AVAILABLE and PERIGEE_CORRECTIONS:
        correction_arcsec = _interpolate_apse_correction(
            jd_tt,
            PERIGEE_CT_JD_START,
            PERIGEE_CT_STEP_DAYS,
            PERIGEE_CORRECTIONS,
            PERIGEE_CT_COUNT,
        )
        interp_lon += correction_arcsec / 3600.0

    interp_lon = interp_lon % 360.0

    # Latitude: inclination * sin(longitude - ascending node)
    # The interpolated perigee lies near the lunar orbital plane, which is
    # inclined ~5.145° to the ecliptic.
    node_lon = _mean_lunar_node_position_unchecked(jd_tt)[0]
    _p1, _p3 = PERIGEE_LAT_COEFFS
    _omega = math.radians(interp_lon - node_lon)
    interp_lat = _p1 * math.sin(_omega) + _p3 * math.sin(3.0 * _omega)

    # Distance: mean DE440 perigee-passage distance (~362,570 km).
    interp_dist = PERIGEE_MEAN_DIST_AU

    return interp_lon, interp_lat, interp_dist


def _solve_linear_system(A: list, b: list) -> list:
    """
    Solve a linear system Ax = b using Gaussian elimination with partial pivoting.

    Args:
        A: Coefficient matrix (list of lists)
        b: Right-hand side vector (list)

    Returns:
        Solution vector x (list)
    """
    n = len(b)

    # Create augmented matrix
    aug = [row[:] + [b[i]] for i, row in enumerate(A)]

    # Forward elimination with partial pivoting
    for col in range(n):
        # Find pivot
        max_row = col
        for row in range(col + 1, n):
            if abs(aug[row][col]) > abs(aug[max_row][col]):
                max_row = row

        # Swap rows
        aug[col], aug[max_row] = aug[max_row], aug[col]

        # Check for singular matrix
        if abs(aug[col][col]) < 1e-12:
            # Matrix is singular or nearly singular, return zeros
            return [0.0] * n

        # Eliminate column
        for row in range(col + 1, n):
            factor = aug[row][col] / aug[col][col]
            for j in range(col, n + 1):
                aug[row][j] -= factor * aug[col][j]

    # Back substitution
    x = [0.0] * n
    for row in range(n - 1, -1, -1):
        x[row] = aug[row][n]
        for col in range(row + 1, n):
            x[row] -= aug[row][col] * x[col]
        x[row] /= aug[row][row]

    return x


def compare_true_lilith_methods(
    jd_tt: float,
) -> dict:
    """
    Compare the two True Lilith calculation methods.

    This function computes True Lilith using both the eccentricity vector
    method and the orbital elements method, returning detailed results for
    comparison and analysis.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT).

    Returns:
        dict: Dictionary containing:
            - 'eccentricity_vector': (lon, lat, e_mag) from calc_true_lilith
            - 'orbital_elements': (lon, lat, e_mag) from calc_true_lilith_orbital_elements
            - 'lon_diff': Difference in longitude (degrees)
            - 'lat_diff': Difference in latitude (degrees)
            - 'e_diff': Difference in eccentricity

    Example:
        >>> result = compare_true_lilith_methods(2451545.0)  # J2000.0
        >>> print(f"Longitude difference: {result['lon_diff']:.4f}°")
    """
    # Compute using eccentricity vector method
    ev_lon, ev_lat, ev_e = calc_true_lilith(jd_tt)

    # Compute using orbital elements method
    oe_lon, oe_lat, oe_e = calc_true_lilith_orbital_elements(jd_tt)

    # Calculate differences
    lon_diff = ev_lon - oe_lon
    if lon_diff > 180:
        lon_diff -= 360
    if lon_diff < -180:
        lon_diff += 360

    lat_diff = ev_lat - oe_lat
    e_diff = ev_e - oe_e

    return {
        "eccentricity_vector": (ev_lon, ev_lat, ev_e),
        "orbital_elements": (oe_lon, oe_lat, oe_e),
        "lon_diff": lon_diff,
        "lat_diff": lat_diff,
        "e_diff": e_diff,
    }


def calc_mean_lunar_node(jd_tt: float) -> float:
    """Return the independent ERFA/IERS mean lunar-node longitude."""
    from .mean_lunar_apse import mean_lunar_node_position

    longitude, _, _ = mean_lunar_node_position(float(jd_tt))
    return float(longitude)


def calc_mean_lilith(jd_tt: float) -> float:
    """Return the independent ERFA/IERS mean lunar-apogee longitude."""
    from .mean_lunar_apse import mean_lunar_apogee_position

    longitude, _, _ = mean_lunar_apogee_position(float(jd_tt))
    return float(longitude)


def calc_mean_lilith_with_latitude(jd_tt: float) -> Tuple[float, float]:
    """Return the independent mean lunar-apogee longitude and latitude."""
    from .mean_lunar_apse import mean_lunar_apogee_position

    longitude, latitude, _ = mean_lunar_apogee_position(float(jd_tt))
    return float(longitude), float(latitude)


def _finite_difference_polar_state(
    position_function,
    jd_tt: float,
    *,
    step_days: float = 0.002,
) -> Tuple[float, float, float, float, float, float]:
    """Return a native-float polar state and centered derivative.

    The half-step matches the osculating-apogee stencil (~3 minutes): the
    interpolated-apsis longitude carries short-period Delaunay/residual
    structure that a half-day chord misrepresents by up to ~10"/day at
    fast-swing phases (the derivative is step-stable from 0.05 down to
    1e-4 days, so 0.002 sits safely inside the plateau while float64
    roundoff stays negligible).
    """
    jd = float(jd_tt)
    longitude, latitude, distance = position_function(jd)
    before = position_function(jd - step_days)
    after = position_function(jd + step_days)
    longitude_delta = (after[0] - before[0] + 180.0) % 360.0 - 180.0
    span = 2.0 * step_days
    return (
        float(longitude),
        float(latitude),
        float(distance),
        float(longitude_delta / span),
        float((after[1] - before[1]) / span),
        float((after[2] - before[2]) / span),
    )


def calc_mean_lunar_node_state(
    jd_tt: float, *, speed3: bool = False
) -> Tuple[float, float, float, float, float, float]:
    """Return the independent mean-node polar state on the mean ecliptic."""
    from .mean_lunar_apse import mean_lunar_node_state

    return mean_lunar_node_state(float(jd_tt), speed3=speed3)


def calc_mean_lilith_state(
    jd_tt: float, *, speed3: bool = False
) -> Tuple[float, float, float, float, float, float]:
    """Return the independent mean-apogee polar state on the mean ecliptic."""
    from .mean_lunar_apse import mean_lunar_apogee_state

    return mean_lunar_apogee_state(float(jd_tt), speed3=speed3)


def calc_interpolated_apse_state(
    jd_tt: float, body_id: int, *, speed3: bool = False
) -> Tuple[float, float, float, float, float, float]:
    """Return the independently generated interpolated lunar-apsis state."""
    del speed3
    from .constants import INTP_APOG, INTP_PERG

    if body_id == INTP_APOG:
        function = calc_interpolated_apogee
    elif body_id == INTP_PERG:
        function = calc_interpolated_perigee
    else:
        raise ValueError(f"Unsupported interpolated lunar apsis body: {body_id}")
    return _finite_difference_polar_state(function, float(jd_tt))


@lru_cache(maxsize=1)
def _get_interpolated_apsides_reader():
    """Retain the asset-free evaluator lifecycle API used by session state."""
    from .interpolated_lunar_apsides import open_interpolated_lunar_apsides

    return open_interpolated_lunar_apsides()


def release_interpolated_apsides_reader() -> None:
    """Close the auxiliary evaluator retained for lifecycle compatibility."""
    if _get_interpolated_apsides_reader.cache_info().currsize:
        _get_interpolated_apsides_reader().close()
    _get_interpolated_apsides_reader.cache_clear()
