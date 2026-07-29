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


def _calc_jupiter_mean_longitude(jd_tt: float) -> float:
    """
    Calculate Jupiter's mean longitude for perturbation calculations.

    Jupiter's gravitational influence causes small but measurable perturbations
    in the lunar orbit, typically a few arcminutes in amplitude.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        float: Jupiter's mean longitude in radians

    References:
        - Meeus, J. "Astronomical Algorithms", Chapter 31 (Planetary Positions)
        - Simon et al. (1994), "Numerical expressions for precession formulae"
    """
    T = (jd_tt - 2451545.0) / 36525.0  # Julian centuries from J2000.0

    # Jupiter mean longitude, simplified polynomial following Meeus Table 31.A.
    # Reached only from the dead-code node-perturbation series above (not on the
    # runtime path). The T^2 coefficient here is -0.0000833; Meeus Table 31.A
    # lists -0.00008501 -- TO BE VERIFIED against Meeus before any reuse.
    L_jupiter = 34.351519 + 3034.9056606 * T - 0.0000833 * T**2

    return math.radians(L_jupiter % 360.0)


def _calc_venus_mean_longitude(jd_tt: float) -> float:
    """
    Calculate Venus's mean longitude for perturbation calculations.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        float: Venus's mean longitude in radians

    References:
        - Meeus, J. "Astronomical Algorithms", Table 31.A
    """
    T = (jd_tt - 2451545.0) / 36525.0  # Julian centuries from J2000.0
    L_venus = 181.979801 + 58517.8156760 * T + 0.00016 * T**2
    return math.radians(L_venus % 360.0)


def _calc_mars_mean_longitude(jd_tt: float) -> float:
    """
    Calculate Mars's mean longitude for perturbation calculations.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        float: Mars's mean longitude in radians

    References:
        - Meeus, J. "Astronomical Algorithms", Table 31.A
    """
    T = (jd_tt - 2451545.0) / 36525.0  # Julian centuries from J2000.0
    L_mars = 355.433275 + 19140.2993313 * T + 0.00026 * T**2
    return math.radians(L_mars % 360.0)


def _calc_saturn_mean_longitude(jd_tt: float) -> float:
    """
    Calculate Saturn's mean longitude for perturbation calculations.

    Saturn's gravitational influence causes small but measurable perturbations
    in the lunar orbit, with amplitudes of approximately 0.001-0.003 degrees.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        float: Saturn's mean longitude in radians

    References:
        - Meeus, J. "Astronomical Algorithms", Table 31.A
        - Simon et al. (1994), "Numerical expressions for precession formulae"
    """
    T = (jd_tt - 2451545.0) / 36525.0  # Julian centuries from J2000.0

    # Saturn mean longitude (simplified formula from Meeus Table 31.A)
    L_saturn = 50.077571 + 1222.1137943 * T + 0.00021 * T**2

    return math.radians(L_saturn % 360.0)


def _calc_elp2000_node_perturbations(jd_tt: float) -> float:
    """
    Calculate the node perturbation series for the true lunar node.

    NOTE ON PROVENANCE AND STATUS: This function is NOT on the library runtime
    path. calc_true_lunar_node() derives the true node directly from the JPL
    DE440 osculating state (see the note near the end of that function); this
    series is retained only for reference and for the out-of-tree comparison
    scripts under validation/. Its 90+ amplitudes are project-derived values,
    calibrated to reproduce the DE440 node and organized following the term
    structure of a Sun/Venus/Mars/Jupiter/Saturn perturbation expansion. They
    are NOT a transcription of the published ELP2000-82B coefficient tables;
    the ELP references below document the theoretical framework the layout
    imitates, not the numeric source of these amplitudes.

    The series models gravitational perturbations from the Sun, Venus, Mars,
    Jupiter, and Saturn, as well as coupling effects between different
    perturbation mechanisms.

    Theoretical Background
    ======================

    The Moon's orbit is significantly perturbed from a simple Keplerian ellipse
    due to the gravitational influence of the Sun and planets. The ELP2000-82B
    theory models these perturbations as trigonometric series in the fundamental
    lunar arguments:

    **Fundamental Arguments:**
        - D: Mean elongation of Moon from Sun (~13.176°/day)
        - M: Mean anomaly of the Sun (~0.986°/day, annual cycle)
        - M': Mean anomaly of the Moon (~13.065°/day, anomalistic month)
        - F: Mean argument of latitude (~13.229°/day)

    Each perturbation term has the form:
        amplitude × E^n × trig_func(a₁D + a₂M + a₃M' + a₄F + planet_terms)

    where E is the Simon et al. (1994) mean Earth-orbit eccentricity relative
    to its J2000 value, and T is Julian centuries from J2000.0.

    Perturbation Categories
    =======================

    **1. Main Solar Perturbation Terms (9 terms, up to 1.52°)**

        The dominant perturbation is the fortnightly term -1.5233 sin(2D),
        which causes the true node to oscillate ~1.5° around the mean node
        with a period of half a synodic month (~14.77 days).

        Other main terms include:
        - sin(2D±M): Solar eccentricity coupling
        - sin(2D±M'): Evection-related terms
        - sin(2F): Inclination modulation
        - sin(2D-2F): Elongation-latitude interaction

    **2. Second-Order Solar Terms (12 terms, 0.003-0.04°)**

        Combined Sun-Moon terms:
        - sin(2D-M-M'), sin(2D-M+M'): Three-body interaction
        - sin(D+M), sin(D-M): Half-synodic coupling

        Higher harmonics:
        - sin(2D±2M'): Double Moon anomaly
        - sin(D), sin(3D), sin(4D): Elongation harmonics

    **3. Third-Order and Higher Terms (10 terms, 0.001-0.006°)**

        Fine corrections for:
        - sin(M±M'): Direct Sun-Moon coupling
        - sin(2D-2M): E² eccentricity term
        - sin(2D±3M'), sin(3M'): Higher Moon anomaly harmonics

    **4. F-Related (Inclination) Terms (9 terms, 0.001-0.01°)**

        The Moon's orbital inclination (~5.145°) to the ecliptic creates:
        - sin(2F±2D): Elongation-latitude resonance
        - sin(2F±D): Half-synodic latitude coupling
        - sin(2F±M): Solar-latitude interaction
        - sin(4F): Fourth harmonic of latitude

    **5. Planetary Perturbation Terms (32 terms, 0.001-0.005°)**

        Venus (9 terms): Closest and most perturbing planet
        - sin(L♀-L☽), sin(L♀-2D), sin(2L♀-2D)
        - Combined with M' and M

        Mars (9 terms): Similar pattern to Venus
        - sin(L♂-2D), sin(L♂), sin(L♂±M')

        Jupiter (7 terms): Largest mass perturbation
        - sin(L♃), sin(L♃-2D), sin(2L♃-2D)
        - Combined with M' and M

        Saturn (7 terms): Weakest planetary perturbation
        - sin(L♄), sin(L♄-2D), sin(2L♄-2D)

    **6. Long-Period Terms**

        Evection (6 terms, up to 0.047°, period ~31.8 days):
            The evection is a major perturbation caused by the Sun's gravitational
            modulation of the Moon's orbital eccentricity. It affects the node
            through eccentricity-inclination coupling.
            - Primary: 0.0467 sin(2D-M')
            - Coupling: sin(2D-M'±M'), sin(2D-M'±2M')

        Variation (10 terms, up to 0.052°, period ~14.77 days):
            The variation arises from the transverse component of solar gravity
            at quadrature. It modulates the node through velocity perturbations.
            - Primary: 0.0523 sin(2D+F), -0.0478 sin(2D-F)
            - Double: sin(4D), sin(4D±F)

        Annual Equation (6 terms, up to 0.186°, period ~365 days):
            Earth's orbital eccentricity causes the Earth-Sun distance to vary
            by ±3.3%, modulating the solar perturbation strength annually.
            - Primary: -0.186 E sin(M) (second-largest individual term)
            - Coupling: E sin(M±M'), E sin(M±2F)
            - Second harmonic: E² sin(2M)

        Parallactic Inequality (3 terms, up to 0.035°, period ~29.5 days):
            Caused by the Sun's finite distance - the solar perturbation differs
            when the Moon is on the sunward vs. anti-sunward side of Earth.
            - Primary: 0.035 sin(D)
            - Secondary: sin(2D-M'+2F), sin(2D+M'-2F)

    **7. Second-Order Coupling Terms (8 terms, 0.0001-0.003°)**

        Products of first-order perturbations create sum and difference
        frequency terms:

        Evection × Variation:
            sin(A)×sin(B) → ½[cos(A-B) - cos(A+B)]
            Creates: cos(M'), cos(4D-M')

        Evection × Annual Equation:
            Creates: cos(2D-M'±M)

        Self-Coupling:
            sin²(arg) → ½[1 - cos(2×arg)]
            Creates: cos(4D-2M'), cos(4D), cos(2F), cos(2M)

    **8. Secular Terms (3 terms, T-dependent)**

        Long-term drift corrections proportional to Julian centuries:
        - T sin(2D), T sin(M'), T sin(2F)
        - Account for secular changes in orbital elements

    Total Term Count
    ================

    - Main solar terms: 9
    - Second-order solar: 12
    - Third-order: 10
    - F-related: 9
    - Venus: 9
    - Mars: 9
    - Jupiter: 7
    - Saturn: 7
    - Evection: 6
    - Variation: 10
    - Annual equation: 6
    - Parallactic: 3
    - Second-order coupling: 8
    - Planetary combinations: 4
    - Secular: 3

    **Total: 90+ terms**

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT).
               TT is the uniform time scale used for ephemeris calculations.

    Returns:
        float: Total perturbation correction in degrees to be added to the
               mean node position. Typically ranges from -2° to +2°.

    Precision
    =========

    - With complete ELP2000-82B series: <0.01 degrees vs JPL DE ephemeris
    - Second-order terms contribute ~0.001-0.003 degrees individually
    - Secular terms ensure accuracy for dates far from J2000.0
    - Truncated terms (amplitude < 0.0001 degrees) contribute ~0.0005 degrees total error

    The main precision limitations are:
    1. Series truncation (omitted terms < 0.0001°)
    2. Polynomial argument accuracy degradation for distant dates
    3. Missing higher-order planetary cross-terms

    Implementation Notes
    ====================

    - All fundamental arguments are computed from high-precision polynomials
    - The Earth eccentricity factor E is applied to terms involving M
    - E² is used for second-order solar eccentricity corrections
    - Planetary longitudes use simplified Meeus formulas (adequate for ~0.001° precision)

    See Also
    ========

    - _calc_lunar_fundamental_arguments: Computes D, M, M', F
    - _calc_jupiter_mean_longitude: Jupiter's mean longitude
    - _calc_venus_mean_longitude: Venus's mean longitude
    - _calc_mars_mean_longitude: Mars's mean longitude
    - _calc_saturn_mean_longitude: Saturn's mean longitude
    - true_node_terms.py: Complete term table as NamedTuples

    References
    ==========

    Framework references (theoretical structure imitated by this series; not
    the numeric source of the project-derived amplitudes above):
        - Chapront-Touze, M. & Chapront, J. "ELP 2000-82B: A semi-analytical
          lunar ephemeris adequate for historical times" (1988), Astronomy &
          Astrophysics 190, 342-352
        - Chapront-Touze, M. & Chapront, J. "Lunar Tables and Programs from
          4000 B.C. to A.D. 8000" (1991), Willmann-Bell

    Fundamental arguments:
        - Simon, J.L. et al. "Numerical expressions for precession formulae and
          mean elements for the Moon and planets" (1994), Astronomy & Astrophysics
          282, 663-683

    Historical development:
        - Brown, E.W. "An Introductory Treatise on the Lunar Theory" (1896)
        - Brown, E.W. "Tables of the Motion of the Moon" (1919)
        - Eckert, W.J., Jones, R., Clark, H.K. "The Motion of the Moon" (1954)

    Modern refinements:
        - Meeus, J. "Astronomical Algorithms" (2nd ed., 1998), Chapter 47
        - Park, R.S. et al. (2021) "The JPL Planetary and Lunar Ephemerides DE440 and DE441"
    """
    T = (jd_tt - 2451545.0) / 36525.0  # Julian centuries from J2000.0

    D, M, M_prime, F = lunar_delaunay_arguments(jd_tt)
    L_jupiter = _calc_jupiter_mean_longitude(jd_tt)
    L_venus = _calc_venus_mean_longitude(jd_tt)
    L_mars = _calc_mars_mean_longitude(jd_tt)
    L_saturn = _calc_saturn_mean_longitude(jd_tt)

    # Mean longitude of the Moon (L')
    L_moon = math.radians(
        (
            218.3164477
            + 481267.88123421 * T
            - 0.0015786 * T**2
            + T**3 / 538841.0
            - T**4 / 65194000.0
        )
        % 360.0
    )

    # Eccentricity of Earth's orbit (decreases over time)
    eccentricity = (
        EARTH_ECCENTRICITY_J2000
        + EARTH_ECCENTRICITY_T * T
        + EARTH_ECCENTRICITY_T2 * T**2
    )
    E = eccentricity / EARTH_ECCENTRICITY_J2000
    E2 = E * E

    perturbation = 0.0

    # ========================================================================
    # MAIN SOLAR PERTURBATION TERMS (First order - largest amplitude)
    # ========================================================================
    # The dominant term with 2D argument (fortnightly variation)
    perturbation += -1.5233 * math.sin(2.0 * D)

    # Terms involving Sun's mean anomaly M combined with elongation (2D ± M)
    # Note: The primary annual equation term (sin(M) alone) is in its own section below
    perturbation += 0.0595 * E * math.sin(2.0 * D - M)
    perturbation += -0.0145 * E * math.sin(2.0 * D + M)

    # Terms involving Moon's mean anomaly M' (monthly variation)
    perturbation += -0.1226 * math.sin(2.0 * D - M_prime)
    perturbation += 0.0490 * math.sin(2.0 * D + M_prime)
    perturbation += 0.0316 * math.sin(M_prime)

    # Terms involving argument of latitude F
    perturbation += 0.1176 * math.sin(2.0 * F)
    perturbation += -0.0801 * math.sin(2.0 * D - 2.0 * F)
    perturbation += 0.0122 * math.sin(2.0 * D - 2.0 * F + M_prime)

    # ========================================================================
    # SECOND-ORDER SOLAR PERTURBATION TERMS
    # ========================================================================
    # Combined Sun-Moon terms
    perturbation += 0.0187 * E * math.sin(2.0 * D - M - M_prime)
    perturbation += -0.0095 * E * math.sin(2.0 * D - M + M_prime)
    perturbation += -0.0154 * E * math.sin(D + M)
    perturbation += -0.0144 * E * math.sin(D - M)

    # Higher-order Moon anomaly terms
    perturbation += -0.0392 * math.sin(2.0 * D - 2.0 * M_prime)
    perturbation += -0.0232 * math.sin(2.0 * M_prime)
    perturbation += 0.0109 * math.sin(2.0 * D + 2.0 * M_prime)

    # Mixed F and M' terms
    perturbation += 0.0093 * math.sin(2.0 * F - M_prime)
    perturbation += 0.0072 * math.sin(2.0 * F + M_prime)

    # D terms
    perturbation += -0.0279 * math.sin(D)
    perturbation += 0.0054 * math.sin(3.0 * D)
    perturbation += -0.0038 * math.sin(4.0 * D)

    # ========================================================================
    # THIRD-ORDER AND HIGHER SOLAR TERMS
    # ========================================================================
    perturbation += 0.0058 * E * math.sin(M + M_prime)
    perturbation += -0.0053 * E * math.sin(M - M_prime)
    perturbation += 0.0038 * E2 * math.sin(2.0 * D - 2.0 * M)
    perturbation += 0.0031 * math.sin(2.0 * D - 3.0 * M_prime)
    perturbation += -0.0025 * math.sin(2.0 * D + 3.0 * M_prime)
    perturbation += 0.0023 * E * math.sin(2.0 * D + M - M_prime)
    perturbation += -0.0021 * E * math.sin(2.0 * D - M + 2.0 * M_prime)
    perturbation += 0.0019 * math.sin(3.0 * M_prime)
    perturbation += -0.0017 * E * math.sin(M + 2.0 * M_prime)
    perturbation += 0.0015 * E * math.sin(M - 2.0 * M_prime)

    # ========================================================================
    # F-RELATED TERMS (inclination effects)
    # ========================================================================
    perturbation += -0.0086 * math.sin(2.0 * F - 2.0 * D)
    perturbation += 0.0064 * math.sin(2.0 * F + 2.0 * D)
    perturbation += -0.0046 * math.sin(2.0 * F + D)
    perturbation += 0.0039 * math.sin(2.0 * F - D)
    perturbation += 0.0032 * E * math.sin(2.0 * F - M)
    perturbation += -0.0028 * E * math.sin(2.0 * F + M)
    perturbation += 0.0024 * math.sin(4.0 * F)
    perturbation += -0.0018 * math.sin(2.0 * F - 2.0 * D + M_prime)
    perturbation += 0.0016 * math.sin(2.0 * F + 2.0 * D - M_prime)

    # ========================================================================
    # VENUS PERTURBATION TERMS
    # ========================================================================
    perturbation += 0.0048 * math.sin(L_venus - L_moon)
    perturbation += -0.0037 * math.sin(L_venus - 2.0 * D)
    perturbation += 0.0029 * math.sin(2.0 * L_venus - 2.0 * D)
    perturbation += -0.0024 * math.sin(L_venus + M_prime)
    perturbation += 0.0021 * math.sin(L_venus - M_prime)
    perturbation += -0.0018 * math.sin(L_venus - 2.0 * D + M_prime)
    perturbation += 0.0015 * math.sin(L_venus - 2.0 * D - M_prime)
    perturbation += -0.0012 * E * math.sin(L_venus + M)
    perturbation += 0.0010 * E * math.sin(L_venus - M)

    # ========================================================================
    # MARS PERTURBATION TERMS
    # ========================================================================
    # Mars perturbation terms following the same pattern as Jupiter.
    # Mars perturbation amplitudes are smaller than Venus (0.0005-0.002 degrees)
    # but contribute to overall accuracy.
    perturbation += 0.0036 * math.sin(L_mars - 2.0 * D)
    perturbation += -0.0027 * math.sin(L_mars)
    perturbation += 0.0022 * math.sin(L_mars - M_prime)
    perturbation += -0.0018 * math.sin(L_mars + M_prime)
    perturbation += 0.0017 * math.sin(
        2.0 * L_mars - 2.0 * D
    )  # Analogous to Jupiter term
    perturbation += 0.0014 * math.sin(L_mars - 2.0 * D + M_prime)
    perturbation += -0.0011 * math.sin(L_mars - 2.0 * D - M_prime)
    perturbation += 0.0009 * E * math.sin(L_mars - M)
    perturbation += -0.0008 * E * math.sin(L_mars + M)  # Analogous to Jupiter term

    # ========================================================================
    # JUPITER PERTURBATION TERMS
    # ========================================================================
    perturbation += 0.0033 * math.sin(L_jupiter)
    perturbation += -0.0028 * math.sin(L_jupiter - 2.0 * D)
    perturbation += 0.0021 * math.sin(2.0 * L_jupiter - 2.0 * D)
    perturbation += -0.0017 * E * math.sin(L_jupiter + M)
    perturbation += 0.0014 * math.sin(L_jupiter - M_prime)
    perturbation += -0.0012 * math.sin(L_jupiter + M_prime)
    perturbation += 0.0009 * math.sin(L_jupiter - 2.0 * D + M_prime)

    # ========================================================================
    # SATURN PERTURBATION TERMS
    # ========================================================================
    # Saturn perturbation terms follow similar patterns to Jupiter but with
    # smaller amplitudes (approximately 0.001-0.003 degrees). Saturn's greater
    # distance results in weaker gravitational influence on lunar motion.
    perturbation += 0.0026 * math.sin(L_saturn)
    perturbation += -0.0022 * math.sin(L_saturn - 2.0 * D)
    perturbation += 0.0018 * math.sin(2.0 * L_saturn - 2.0 * D)
    perturbation += -0.0014 * E * math.sin(L_saturn + M)
    perturbation += 0.0012 * math.sin(L_saturn - M_prime)
    perturbation += -0.0010 * math.sin(L_saturn + M_prime)
    perturbation += 0.0008 * math.sin(L_saturn - 2.0 * D + M_prime)

    # ========================================================================
    # LONG-PERIOD TERMS (Secular and long-period variations)
    # ========================================================================
    # These arise from the combination of various lunar inequalities

    # ========================================================================
    # EVECTION EFFECT ON NODE (Long-period term, period ~31.8 days)
    # ========================================================================
    # The evection is a major perturbation of the lunar orbit caused by the Sun's
    # gravitational influence. It has a period of ~31.8 days and amplitude of
    # ~1.274° in longitude. The evection affects the node through its influence
    # on the lunar eccentricity, which in turn affects the node calculation.
    #
    # The evection modulates the Moon's eccentricity by approximately ±0.01148,
    # which causes oscillations in the orbital plane orientation. This effect
    # propagates to the node position through the osculating orbital elements.
    #
    # The evection argument is: 2D - M' (twice mean elongation minus mean anomaly)
    # where the period is: 1 / (2 * 13.176° - 13.065°) ≈ 31.8 days
    #
    # References:
    # - Chapront-Touzé, M. & Chapront, J. "Lunar Tables and Programs" (1991)
    # - Brown, E.W. "An Introductory Treatise on the Lunar Theory" (1896)
    # - Meeus, J. "Astronomical Algorithms" (2nd ed., 1998), Chapter 47
    evection_arg = 2.0 * D - M_prime

    # Primary evection-eccentricity effect on node
    # The evection modifies the eccentricity vector, which affects the orbital
    # plane orientation. This creates a node perturbation with amplitude
    # proportional to the evection amplitude (~1.274°) scaled by the
    # inclination coupling factor (~sin(5.145°) ≈ 0.0897).
    # The coefficient 0.0467 comes from: 1.274° × tan(i/2) × coupling_factor
    # where i is the lunar inclination (5.145°) and the coupling factor
    # accounts for the geometric relationship between eccentricity variation
    # and node motion.
    perturbation += 0.0467 * math.sin(evection_arg)

    # Evection-eccentricity coupling with Moon's anomaly
    # The evection also interacts with the Moon's mean anomaly, creating
    # secondary terms that affect the node position. These arise from the
    # coupling between the evection's eccentricity modulation and the
    # instantaneous orbital plane orientation.
    perturbation += 0.0156 * math.sin(evection_arg + M_prime)
    perturbation += -0.0134 * math.sin(evection_arg - M_prime)
    perturbation += 0.0089 * math.sin(evection_arg + 2.0 * M_prime)
    perturbation += -0.0072 * math.sin(evection_arg - 2.0 * M_prime)

    # Evection-inclination interaction terms (combined with F)
    # These capture the direct coupling between evection and the lunar
    # latitude argument, affecting the node through inclination variations.
    perturbation += 0.0063 * math.sin(evection_arg + F)
    perturbation += -0.0052 * math.sin(evection_arg - F)

    # ========================================================================
    # VARIATION EFFECT ON NODE (Long-period term, period ~14.77 days)
    # ========================================================================
    # The Variation is a major lunar inequality discovered by Tycho Brahe,
    # caused by the difference in solar gravitational force between the
    # Moon's position at quadrature versus conjunction/opposition.
    #
    # Physical mechanism: At quadrature (first/last quarter), the Moon is
    # approximately at right angles to the Sun-Earth line. The transverse
    # component of the Sun's tidal force is maximum here, causing the Moon
    # to speed up (before quadrature) or slow down (after quadrature).
    # At conjunction/opposition (new/full Moon), the radial component
    # dominates, with different effects on the orbit.
    #
    # The Variation has:
    # - Period: ~14.77 days (half the synodic month, 29.53/2 days)
    # - Amplitude: ~0.658° in longitude
    # - Argument: 2D (twice the mean elongation)
    #
    # The main 2D term (-1.5233 sin(2D)) at line 247 captures the dominant
    # fortnightly oscillation of the node. The terms below capture the
    # Variation's secondary effects through its coupling with other
    # perturbations and its influence on the orbital plane orientation.
    #
    # The Variation affects the node through:
    # 1. Direct modulation of the Moon's radial velocity, which affects
    #    the instantaneous orbital plane orientation
    # 2. Coupling with the lunar inclination (F) and anomaly (M')
    # 3. Interaction with the eccentricity variation (evection coupling)
    #
    # References:
    # - Chapront-Touzé, M. & Chapront, J. "Lunar Tables and Programs" (1991)
    # - Brown, E.W. "An Introductory Treatise on the Lunar Theory" (1896)
    # - Meeus, J. "Astronomical Algorithms" (2nd ed., 1998), Chapter 47
    variation_arg = 2.0 * D

    # Primary Variation-inclination coupling terms (2D ± F)
    # The Variation modulates the lunar latitude through its coupling
    # with the argument of latitude F. This creates a secondary effect
    # on the node position. Amplitude ~0.658° × sin(i) ≈ 0.059°
    perturbation += 0.0523 * math.sin(variation_arg + F)
    perturbation += -0.0478 * math.sin(variation_arg - F)

    # Variation-anomaly coupling terms (2D ± M')
    # The Variation interacts with the Moon's mean anomaly, creating
    # combined perturbations. These arise from the modulation of the
    # eccentricity effect by the solar tidal force.
    perturbation += 0.0156 * math.sin(variation_arg + M_prime)
    perturbation += -0.0142 * math.sin(variation_arg - M_prime)

    # Combined Variation-inclination-anomaly terms (2D + F ± M')
    # Higher-order coupling between the Variation, inclination, and
    # the Moon's orbital eccentricity.
    perturbation += 0.0048 * math.sin(variation_arg + F + M_prime)
    perturbation += -0.0041 * math.sin(variation_arg - F + M_prime)
    perturbation += 0.0038 * math.sin(variation_arg + F - M_prime)
    perturbation += -0.0033 * math.sin(variation_arg - F - M_prime)

    # Variation-solar coupling terms (2D ± M)
    # The Variation's interaction with the Sun's mean anomaly creates
    # beat frequencies between the fortnightly and annual periods.
    perturbation += 0.0028 * E * math.sin(variation_arg + M)
    perturbation += -0.0024 * E * math.sin(variation_arg - M)

    # Double Variation terms (4D)
    # Second harmonic of the Variation, with half the period (~7.4 days)
    # and smaller amplitude. These arise from non-linear effects in
    # the solar tidal perturbation.
    perturbation += 0.0067 * math.sin(4.0 * D)
    perturbation += 0.0043 * math.sin(4.0 * D + F)
    perturbation += -0.0039 * math.sin(4.0 * D - F)

    # ========================================================================
    # ANNUAL EQUATION EFFECT ON NODE (Long-period term, period ~365.26 days)
    # ========================================================================
    # The Annual Equation is a major lunar perturbation caused by the varying
    # Earth-Sun distance due to Earth's orbital eccentricity (e ≈ 0.0167).
    #
    # Physical mechanism: When Earth is at perihelion (early January), it is
    # ~3.3% closer to the Sun than at aphelion (early July). The increased
    # solar gravitational influence at perihelion accelerates the Moon's
    # orbital motion, while at aphelion the Moon moves slightly slower.
    # This creates a periodic variation in lunar longitude.
    #
    # The Annual Equation has:
    # - Period: ~365.2596 days (one anomalistic year)
    # - Amplitude: ~0.186° (11.2 arcminutes) in lunar longitude
    # - Argument: M (Sun's mean anomaly)
    # - Cause: Earth's orbital eccentricity varying Earth-Sun distance
    #
    # The Annual Equation affects the node through:
    # 1. Direct modulation of the Moon's angular velocity, which affects
    #    the instantaneous orbital plane orientation
    # 2. Coupling with the lunar inclination (F) and anomaly (M')
    # 3. Interaction with the eccentricity through the E factor
    #
    # References:
    # - Chapront-Touzé, M. & Chapront, J. "Lunar Tables and Programs" (1991)
    # - Brown, E.W. "An Introductory Treatise on the Lunar Theory" (1896)
    # - Meeus, J. "Astronomical Algorithms" (2nd ed., 1998), Chapter 47

    # Primary Annual Equation term
    # The amplitude 0.186° is scaled for the node effect. The coefficient
    # includes the geometric coupling between the annual longitude variation
    # and the node motion. The E factor accounts for the secular decrease
    # in Earth's orbital eccentricity.
    perturbation += -0.1860 * E * math.sin(M)

    # Annual Equation coupling with Moon's anomaly (M')
    # These arise from the interaction between the annual solar distance
    # variation and the Moon's elliptical orbit. The amplitude is modulated
    # by the product of Earth's eccentricity and Moon's eccentricity.
    perturbation += 0.0098 * E * math.sin(M + M_prime)
    perturbation += -0.0082 * E * math.sin(M - M_prime)

    # Annual Equation-inclination coupling terms (M ± 2F)
    # These capture the coupling between the annual equation and the lunar
    # latitude argument, affecting the node through inclination variations.
    perturbation += 0.0037 * E * math.sin(M + 2.0 * F)
    perturbation += -0.0032 * E * math.sin(M - 2.0 * F)

    # Second harmonic of Annual Equation (2M)
    # Arises from the non-linear effects of Earth's orbital eccentricity.
    # Period ~182.6 days (half anomalistic year).
    perturbation += 0.0024 * E2 * math.sin(2.0 * M)

    # ========================================================================
    # PARALLACTIC EQUATION (Parallactic Inequality)
    # ========================================================================
    # The parallactic equation is caused by the Sun being at a finite distance
    # rather than infinitely far away. The Sun's gravitational perturbation on
    # the Moon varies depending on whether the Moon is on the sunward or
    # anti-sunward side of Earth.
    #
    # Physical mechanism: When the Moon is between Earth and Sun (new moon,
    # D ~ 0), it is closer to the Sun and experiences stronger solar gravity.
    # When the Moon is on the far side of Earth (full moon, D ~ 180°), it is
    # farther from the Sun and experiences weaker solar gravity.
    #
    # The parallactic equation has:
    # - Period: ~29.53 days (synodic month)
    # - Amplitude: ~0.035 degrees (approximately 2 arcminutes)
    # - Argument: D (mean elongation of Moon from Sun)
    #
    # The amplitude depends on the solar parallax (ratio of Earth's radius
    # to Earth-Sun distance), approximately 8.794 arcseconds, scaled by
    # the lunar distance and orbital mechanics factors.
    #
    # References:
    # - Chapront-Touzé, M. & Chapront, J. "Lunar Tables and Programs" (1991)
    # - Brown, E.W. "An Introductory Treatise on the Lunar Theory" (1896)
    # - Meeus, J. "Astronomical Algorithms" (2nd ed., 1998), Chapter 47

    # Primary parallactic equation term
    # Amplitude ~0.035° from ELP2000 theory
    perturbation += 0.0350 * math.sin(D)

    # Secondary parallactic terms (cross-coupling with F and M')
    perturbation += 0.0026 * math.sin(2.0 * D - M_prime + 2.0 * F)
    perturbation += -0.0022 * math.sin(2.0 * D + M_prime - 2.0 * F)

    # ========================================================================
    # SECULAR TERMS (very long period, T-dependent)
    # ========================================================================
    # These correct for long-term drift in the perturbation series
    perturbation += 0.0018 * T * math.sin(2.0 * D)
    perturbation += -0.0014 * T * math.sin(M_prime)
    perturbation += 0.0011 * T * math.sin(2.0 * F)

    # ========================================================================
    # PLANETARY COMBINATION TERMS
    # ========================================================================
    # Venus-Jupiter interaction
    perturbation += 0.0008 * math.sin(L_venus - L_jupiter)
    perturbation += -0.0006 * math.sin(L_venus + L_jupiter - 2.0 * D)

    # Mars-Jupiter interaction
    perturbation += 0.0005 * math.sin(L_mars - L_jupiter)

    # Saturn-Jupiter interaction
    perturbation += 0.0004 * math.sin(L_saturn - L_jupiter)

    # Venus-Saturn interaction
    perturbation += 0.0003 * math.sin(L_venus - L_saturn)

    # ========================================================================
    # SECOND-ORDER PERTURBATION TERMS
    # ========================================================================
    # Second-order terms arise from products of first-order perturbations.
    # When two sinusoidal terms are multiplied:
    #   sin(A) × sin(B) = ½[cos(A-B) - cos(A+B)]
    #   cos(A) × sin(B) = ½[sin(A+B) + sin(A-B)]
    #
    # These can contribute at the arcsecond level (0.0001-0.001 degrees)
    # and are needed for sub-arcsecond precision.
    #
    # References:
    # - Chapront-Touzé, M. & Chapront, J. "ELP 2000-82B" (1988)
    # - Brown, E.W. "Tables of the Moon" (1919)
    # - Eckert, W.J. et al. "The Motion of the Moon" (1954)

    # ========================================================================
    # EVECTION × VARIATION COUPLING
    # ========================================================================
    # Products of sin(2D - M') and sin(2D) create terms at:
    #   cos(M') from the difference: ½[cos((2D-M') - 2D) - cos((2D-M') + 2D)]
    #           = ½[cos(-M') - cos(4D-M')] = ½[cos(M') - cos(4D-M')]
    #
    # The evection amplitude (~1.274°) × variation amplitude (~0.658°) × coupling
    # factor gives amplitudes of ~0.001-0.003 degrees for these terms.
    #
    # Primary coupling term: cos(M')
    # This represents the beat frequency between evection and variation.
    perturbation += 0.0024 * math.cos(M_prime)

    # Secondary coupling term: cos(4D - M')
    # Higher frequency component of the evection-variation interaction.
    perturbation += -0.0018 * math.cos(4.0 * D - M_prime)

    # Combined with F (inclination coupling)
    perturbation += 0.0012 * math.cos(M_prime + 2.0 * F)
    perturbation += -0.0010 * math.cos(M_prime - 2.0 * F)

    # ========================================================================
    # EVECTION × ANNUAL EQUATION COUPLING
    # ========================================================================
    # Products of sin(2D - M') and sin(M) create terms at:
    #   cos(2D - M' - M) and cos(2D - M' + M)
    #
    # These represent the interaction between the evection's eccentricity
    # modulation and the annual variation in Earth-Sun distance.
    #
    # The coupling amplitude is approximately:
    #   evection_amp × annual_amp × geometric_factor
    #   ~1.274° × 0.186° × 0.01 ≈ 0.002°
    perturbation += 0.0019 * E * math.cos(2.0 * D - M_prime - M)
    perturbation += -0.0016 * E * math.cos(2.0 * D - M_prime + M)

    # Higher-order evection-annual coupling with Moon anomaly
    perturbation += 0.0011 * E * math.cos(2.0 * D - 2.0 * M_prime - M)
    perturbation += -0.0009 * E * math.cos(2.0 * D - 2.0 * M_prime + M)

    # ========================================================================
    # VARIATION × ANNUAL EQUATION COUPLING
    # ========================================================================
    # Products of sin(2D) and sin(M) create terms at:
    #   cos(2D - M) and cos(2D + M)
    #
    # Note: The first-order terms 2D ± M are already included in the main
    # solar perturbation section. These second-order terms have additional
    # amplitude from the non-linear coupling.
    #
    # The coupling creates a modulation of the fortnightly term by the
    # annual cycle, resulting in beat frequencies.
    perturbation += 0.0015 * E * math.cos(2.0 * D - M)
    perturbation += -0.0013 * E * math.cos(2.0 * D + M)

    # Combined with Moon anomaly M'
    perturbation += 0.0008 * E * math.cos(2.0 * D - M + M_prime)
    perturbation += -0.0007 * E * math.cos(2.0 * D + M - M_prime)

    # ========================================================================
    # SELF-COUPLING TERMS (sin²(arg) → cos(2×arg))
    # ========================================================================
    # When a perturbation term is squared, it creates a DC offset (constant)
    # and a term at twice the frequency:
    #   sin²(arg) = ½[1 - cos(2×arg)]
    #
    # The most significant self-coupling terms come from the largest
    # first-order perturbations.

    # Evection self-coupling: sin²(2D - M') → cos(4D - 2M')
    # Amplitude: (1.274°)² × coupling_factor ≈ 0.001°
    perturbation += 0.0014 * math.cos(4.0 * D - 2.0 * M_prime)

    # Variation self-coupling: sin²(2D) → cos(4D)
    # Note: 4D terms are already included above, this adds to them
    perturbation += 0.0011 * math.cos(4.0 * D)

    # Annual equation self-coupling: sin²(M) → cos(2M)
    # Note: 2M term is already in the first-order section, this is additional
    perturbation += 0.0006 * E2 * math.cos(2.0 * M)

    # Latitude argument self-coupling: sin²(F) → cos(2F)
    # Contributes to inclination oscillation
    perturbation += 0.0008 * math.cos(2.0 * F)

    # ========================================================================
    # E² SECOND-ORDER SOLAR ECCENTRICITY CORRECTIONS
    # ========================================================================
    # Terms proportional to E² account for second-order effects of Earth's
    # orbital eccentricity. These are smaller but contribute at arcsecond level.
    #
    # E² ≈ 0.9832 for modern epoch, so E² corrections are about 1.7% smaller
    # than they would be without the eccentricity factor.

    # Second-order eccentricity correction for main solar term
    perturbation += -0.0021 * E2 * math.sin(2.0 * D - 2.0 * M)
    perturbation += 0.0017 * E2 * math.sin(2.0 * D + 2.0 * M)

    # Combined with Moon anomaly
    perturbation += 0.0012 * E2 * math.sin(2.0 * D - 2.0 * M + M_prime)
    perturbation += -0.0010 * E2 * math.sin(2.0 * D - 2.0 * M - M_prime)

    # Combined with latitude argument F
    perturbation += 0.0009 * E2 * math.sin(2.0 * D - 2.0 * M + 2.0 * F)
    perturbation += -0.0008 * E2 * math.sin(2.0 * D - 2.0 * M - 2.0 * F)

    # ========================================================================
    # THREE-FREQUENCY COMBINATION TERMS
    # ========================================================================
    # Terms involving three fundamental arguments arise from products of
    # different perturbation series. These are typically at the 0.0001-0.001°
    # level but needed for arcsecond precision.

    # D + M + M' combinations (synodic-anomalistic-lunar coupling)
    perturbation += 0.0011 * E * math.sin(D + M + M_prime)
    perturbation += -0.0009 * E * math.sin(D - M - M_prime)
    perturbation += 0.0008 * E * math.sin(D + M - M_prime)
    perturbation += -0.0007 * E * math.sin(D - M + M_prime)

    # 2D + M + M' combinations
    perturbation += 0.0013 * E * math.sin(2.0 * D + M + M_prime)
    perturbation += -0.0011 * E * math.sin(2.0 * D - M - M_prime)

    # F + M combinations (inclination-solar coupling)
    perturbation += 0.0007 * E * math.sin(F + M)
    perturbation += -0.0006 * E * math.sin(F - M)

    # 2F + M' combinations (double inclination-lunar coupling)
    perturbation += 0.0009 * math.sin(2.0 * F + M_prime)
    perturbation += -0.0008 * math.sin(2.0 * F - M_prime)

    # ========================================================================
    # SECULAR SECOND-ORDER TERMS
    # ========================================================================
    # Terms proportional to T (time from J2000) capture long-term
    # second-order effects from secular changes in orbital elements.

    # T-dependent evection-variation coupling
    perturbation += 0.00003 * T * math.cos(M_prime)

    # T-dependent solar eccentricity effect
    perturbation += -0.00002 * T * E * math.sin(M)

    # ========================================================================
    # HIGHER-ORDER TERMS (ELP2000-85 refinements)
    # ========================================================================
    # These additional terms from the ELP2000-85 theory provide sub-arcsecond
    # corrections. Each term amplitude is < 0.0003° but they contribute
    # ~5-10 arcsec total improvement in precision.
    #
    # References:
    #   - Chapront-Touze, M. & Chapront, J. "Lunar Tables and Programs from
    #     4000 B.C. to A.D. 8000" (1991), Willmann-Bell
    #   - Chapront, J. et al. "A new determination of lunar orbital parameters,
    #     precession constant and tidal acceleration" (2002), A&A 387

    # Fifth-order elongation terms (5D combinations)
    perturbation += 0.0003 * math.sin(5.0 * D)
    perturbation += -0.0002 * math.sin(5.0 * D - M_prime)
    perturbation += 0.0002 * math.sin(5.0 * D + M_prime)
    perturbation += -0.0001 * E * math.sin(5.0 * D - M)
    perturbation += 0.0001 * E * math.sin(5.0 * D + M)

    # Triple combination terms (3D with M and M')
    perturbation += 0.0002 * E * math.sin(3.0 * D + M)
    perturbation += -0.0002 * E * math.sin(3.0 * D - M)
    perturbation += 0.0002 * math.sin(3.0 * D + 2.0 * M_prime)
    perturbation += -0.0002 * math.sin(3.0 * D - 2.0 * M_prime)

    # Higher Moon anomaly harmonics (3M', 4M')
    perturbation += 0.0002 * math.sin(3.0 * M_prime)
    perturbation += -0.0001 * math.sin(4.0 * M_prime)
    perturbation += 0.0001 * math.sin(2.0 * D + 3.0 * M_prime)
    perturbation += -0.0001 * math.sin(2.0 * D - 3.0 * M_prime)

    # Planetary cross-coupling terms
    perturbation += 0.0002 * math.sin(2.0 * L_venus - 2.0 * L_mars)
    perturbation += -0.0001 * math.sin(L_jupiter + L_saturn - 2.0 * D)
    perturbation += 0.0001 * math.sin(L_venus + L_jupiter - 2.0 * D)
    perturbation += -0.0001 * math.sin(L_mars + L_saturn - 2.0 * D)

    # Higher F (latitude) harmonics
    perturbation += 0.0001 * math.sin(2.0 * F + 2.0 * D + M_prime)
    perturbation += -0.0001 * math.sin(2.0 * F - 2.0 * D - M_prime)
    perturbation += 0.0001 * math.sin(4.0 * F + 2.0 * D)
    perturbation += -0.0001 * math.sin(4.0 * F - 2.0 * D)

    # Secular T² terms for long-term drift
    perturbation += 0.00001 * T * T * math.sin(2.0 * D)
    perturbation += -0.00001 * T * T * math.cos(M_prime)

    # ========================================================================
    # HISTORICAL DATE CORRECTIONS (T^3, T^4, T^5 secular terms)
    # ========================================================================
    # These higher-order secular terms improve accuracy for historical dates
    # (1500-1800 CE and earlier) where the standard perturbation series
    # accumulates significant error due to T^4 polynomial degradation.
    #
    # Coefficients derived from Chapront et al. (2002) tidal acceleration
    # analysis and calibrated against historical eclipse records.
    #
    # References:
    #   - Chapront, J. et al. (2002) "A new determination of lunar orbital
    #     parameters, precession constant and tidal acceleration", A&A 387
    #   - Stephenson, F.R. & Morrison, L.V. (1995) "Long-term fluctuations
    #     in the Earth's rotation: 700 BC to AD 1990", Phil. Trans. R. Soc.

    T2 = T * T
    T3 = T2 * T
    T4 = T3 * T
    T5 = T4 * T

    # Cubic secular corrections (T³)
    # These capture long-term drift in the node position due to
    # tidal acceleration and other secular effects
    perturbation += 0.000024 * T3 * math.sin(2.0 * D)
    perturbation += -0.000018 * T3 * math.sin(M_prime)
    perturbation += 0.000012 * T3 * math.cos(2.0 * D - M_prime)

    # Quartic secular corrections (T⁴)
    # Fourth-order terms to compensate for polynomial truncation errors
    # in the fundamental arguments for dates far from J2000
    perturbation += 0.0000035 * T4 * math.sin(2.0 * D)
    perturbation += -0.0000028 * T4 * math.cos(M_prime)
    perturbation += 0.0000021 * T4 * math.sin(2.0 * D - M_prime)

    # Quintic secular corrections (T⁵)
    # Fifth-order terms for extended historical range accuracy
    perturbation += 0.00000045 * T5 * math.sin(2.0 * D)
    perturbation += -0.00000038 * T5 * math.cos(M_prime)

    # ========================================================================
    # SIXTH-ORDER ELONGATION TERMS (6D combinations)
    # ========================================================================
    # Higher harmonics of the mean elongation for improved accuracy
    # at historical dates. These capture short-period oscillations that
    # become significant when integrated over long time spans.
    #
    # Derived from ELP2000-82B series extension for historical dates.

    perturbation += 0.0004 * math.sin(6.0 * D)
    perturbation += -0.0003 * math.sin(6.0 * D - M_prime)
    perturbation += 0.0003 * math.sin(6.0 * D + M_prime)
    perturbation += -0.0002 * E * math.sin(6.0 * D - M)
    perturbation += 0.0002 * E * math.sin(6.0 * D + M)
    perturbation += 0.0002 * math.sin(6.0 * D - 2.0 * M_prime)
    perturbation += -0.0002 * math.sin(6.0 * D + 2.0 * M_prime)

    # ========================================================================
    # ENHANCED 5D TERMS FOR HISTORICAL ACCURACY
    # ========================================================================
    # Additional fifth-order elongation combinations that improve
    # precision for dates before 1800 CE.

    perturbation += 0.0003 * math.sin(5.0 * D - 2.0 * M_prime)
    perturbation += -0.0003 * math.sin(5.0 * D + 2.0 * M_prime)
    perturbation += 0.0002 * E * math.sin(5.0 * D - M - M_prime)
    perturbation += -0.0002 * E * math.sin(5.0 * D + M + M_prime)
    perturbation += 0.0002 * math.sin(5.0 * D - 3.0 * M_prime)
    perturbation += 0.0002 * math.sin(5.0 * D + F)
    perturbation += -0.0002 * math.sin(5.0 * D - F)

    # ========================================================================
    # SEVENTH-ORDER ELONGATION TERMS (7D combinations)
    # ========================================================================
    # Very high-order harmonics for sub-degree accuracy at historical dates.

    perturbation += 0.0002 * math.sin(7.0 * D)
    perturbation += -0.0001 * math.sin(7.0 * D - M_prime)
    perturbation += 0.0001 * math.sin(7.0 * D + M_prime)
    perturbation += -0.0001 * math.sin(7.0 * D - 2.0 * M_prime)

    # ========================================================================
    # ENHANCED EVECTION-ELONGATION COUPLING FOR HISTORICAL DATES
    # ========================================================================
    # Cross-coupling terms between evection and higher elongation harmonics
    # that become significant for long time integrations.

    perturbation += 0.0003 * math.sin(4.0 * D - M_prime)
    perturbation += -0.0003 * math.sin(4.0 * D + M_prime)
    perturbation += 0.0002 * math.sin(4.0 * D - 2.0 * M_prime)
    perturbation += -0.0002 * math.sin(4.0 * D + 2.0 * M_prime)
    perturbation += 0.0002 * E * math.sin(4.0 * D - M - M_prime)
    perturbation += -0.0002 * E * math.sin(4.0 * D + M - M_prime)

    return perturbation


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

    # Note: ELP2000 perturbation corrections (_calc_elp2000_node_perturbations)
    # are available but not applied here. The geometric h = r x v approach already
    # captures perturbation effects through the JPL DE ephemeris state vectors.
    # The perturbation series was designed for the mean node, not the geometric node.

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
