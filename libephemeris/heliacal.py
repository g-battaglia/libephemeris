# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Heliacal event calculations for libephemeris.

Calculates heliacal risings and settings, visual limiting magnitude,
and related heliacal phenomena for celestial bodies.

Functions:
- _heliacal_ut_pythonic: Find heliacal rising/setting time for a body
- heliacal_ut: Alias for _heliacal_ut_pythonic
- _heliacal_pheno_ut_pythonic: Calculate detailed heliacal phenomena
- heliacal_pheno_ut: Alias for _heliacal_pheno_ut_pythonic
- vis_limit_mag: Calculate visual limiting magnitude
- vis_limit_mag: Alias for vis_limit_mag

Historical Note:
    Heliacal risings were crucial for ancient calendars. The heliacal
    rising of Sirius marked the Egyptian new year and predicted the
    Nile flood. Babylonians used heliacal events to track planetary
    positions without modern instruments.

References:
    - Schaefer, B.E. (1990), "Telescopic Limiting Magnitudes", PASP 102,
      212-229, DOI 10.1086/132629
    - Schaefer, B.E. (1993), "Astronomy and the Limits of Vision", Vistas in
      Astronomy 36, 311-361; and Schaefer, B.E. (1998), "To the Visual
      Limits", Sky & Telescope 95(5), 57-60. Together these define the
      VISLIMIT V-band limiting-magnitude algorithm and its per-component
      relative-airmass relations used by this module's active path (the
      low-altitude airmass form follows Rozenberg, G.V. (1966), "Twilight: A
      Study in Atmospheric Optics"; Kasten & Young (1989) is NOT used here).
    - Krisciunas, K. & Schaefer, B.E. (1991), "A model of the brightness of
      moonlight", PASP 103, 1033-1039, DOI 10.1086/132921 (scattered
      moonlight and the daylight sky-scattering terms)
    - Yallop, B.D. (1997), NAO Technical Note No. 69 (lunar-crescent q
      criterion); Bruin, F. (1977), "The first visibility of the lunar
      crescent", Vistas in Astronomy 21, 331-358 (crescent width)

Provenance:
    Body positions and photometry come from the registered ephemeris pipeline;
    atmospheric extinction, sky brightness, contrast, observer optics, and the
    lunar-crescent classes come from the cited public models. This module labels
    all project decisions separately: event scan cadence, root tolerance,
    neutral visibility margin, boundary clamps, fallback order, and fixed result
    slots. No prose, constant table, or hidden visibility curve is imported from
    another ephemeris implementation.
"""

from __future__ import annotations

import math
from typing import Callable, Tuple, NamedTuple, Optional, Sequence, Union

import numpy as np

from .constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    FIXSTAR_OFFSET,
    FLG_SPEED,
    FLG_SWIEPH,
    HELFLAG_HIGH_PRECISION,
    HELFLAG_VISLIM_PHOTOPIC,
    HELFLAG_VISLIM_SCOTOPIC,
)

from .exceptions import Error, LEBCorruptionError

# Inner planets (orbit inside Earth's orbit)
# These have both inferior conjunction (between Earth and Sun)
# and superior conjunction (behind the Sun)
INNER_PLANETS = {MERCURY, VENUS}


def _yallop_visibility_code(q: float) -> float:
    """Map a Yallop q-test value to its new-crescent visibility class code.

    Yallop (1997), NAO Technical Note No. 69, classifies first-visibility of
    the new lunar crescent into six bands from the q-test value
    ``q = (ARCV - criterion) / 10``:

        A: ``q > +0.216``            -> 1  (easily visible)
        B: ``+0.216 >= q > -0.014``  -> 2  (visible under perfect conditions)
        C: ``-0.014 >= q > -0.160``  -> 3  (may need optical aid to find)
        D: ``-0.160 >= q > -0.232``  -> 4  (will need optical aid to find)
        E: ``-0.232 >= q > -0.293``  -> 5  (not visible with a telescope)
        F: ``q <= -0.293``           -> 6  (not visible, below Danjon limit)

    The public result contract places this discrete band code (not the raw
    arcus-visionis criterion polynomial) in ``heliacal_pheno_ut`` element 18.

    Args:
        q: Yallop q-test value (``heliacal_pheno_ut`` element 17).

    Returns:
        Visibility class code 1..6, returned as a float like every result slot.
    """
    if q > 0.216:
        return 1.0
    if q > -0.014:
        return 2.0
    if q > -0.160:
        return 3.0
    if q > -0.232:
        return 4.0
    if q > -0.293:
        return 5.0
    return 6.0


def _arc_of_light(arcv_deg: float, daz_deg: float) -> float:
    """Arc of light (ARCL) from the arcus visionis and azimuth separation.

    ARCV (vertical), DAZ (horizontal), and ARCL (the hypotenuse arc between the
    body and the Sun) form a right-angled spherical triangle, so
    ``cos ARCL = cos ARCV * cos DAZ`` (see, e.g., Yallop, B.D. (1997), NAO
    Technical Note No. 69, and Bruin, F. (1977), Vistas in Astronomy 21, 331).
    The result is a non-negative arc in [0, 180] degrees.

    This is the exact spherical relation, not the ``sqrt(ARCV^2 + DAZ^2)``
    small-angle (planar) approximation, which underestimates the arc as the
    body moves away from the Sun. ``cos`` is even and 360-periodic, so the
    result is invariant to how ``daz_deg`` is normalised.

    Args:
        arcv_deg: Geocentric arcus visionis (object-minus-Sun altitude), deg.
        daz_deg: Object-minus-Sun azimuth separation, deg.

    Returns:
        Arc of light between the body and the Sun, in degrees.
    """
    c = math.cos(math.radians(arcv_deg)) * math.cos(math.radians(daz_deg))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


# Detection margin (magnitudes) for the twilight visibility test used by the
# heliacal-event search. With the published VISLIMIT limiting-magnitude model
# (dret[0] is a catalog-magnitude limit, extinction folded in) an object is at
# the visibility threshold when its catalog magnitude equals the limiting
# magnitude, so the physically neutral margin is zero.
_HELIACAL_VIS_MARGIN = 0.0


# =============================================================================
# SCHAEFER (1990) ATMOSPHERIC MODEL
# =============================================================================
#
# Project implementation of the cited Schaefer visibility-model family:
#   - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
#   - Schaefer, B.E. (1993) "Astronomy and the Limits of Vision"
#
# The model calculates:
#   1. Atmospheric extinction (Rayleigh + Aerosol + Ozone)
#   2. Sky brightness (twilight + moonlight + airglow)
#   3. Numerically solved arcus-visionis threshold over the limiting-magnitude
#      model (historical terminology does not supply a coefficient table)
#   4. Limiting visual magnitude for naked-eye observation
# =============================================================================


class SchaeferConstants:
    """
    Legacy descriptive constants for a Schaefer-style atmospheric model.

    NOTE: these values are documentation/compatibility defaults only. The
    active VISLIMIT limiting-magnitude path (``SchaeferModel`` and the
    ``_vl_*`` helpers above) does not read this class; it uses its own
    coefficients. Citations below are given for the published origin of each
    representative value, not because this class drives any result.
    """

    # Rayleigh (molecular) scattering coefficient at sea level, V band
    # (~550 nm), per unit airmass. The 0.1451 mag/airmass V-band value derives
    # from the atmospheric-extinction model of Hayes, D.S. & Latham, D.W.
    # (1975), "A rediscussion of the atmospheric extinction and the absolute
    # spectral-energy distribution of Vega", ApJ 197, 593-601.
    K_RAYLEIGH_SEA_LEVEL = 0.1451

    # Ozone absorption coefficient (per airmass). Small, roughly constant
    # contribution at visual wavelengths (Chappuis bands); representative
    # V-band value used by Schaefer (1990), PASP 102, 212.
    K_OZONE = 0.016

    # Aerosol (Mie) scattering base coefficient at V, clear-air baseline.
    # Aerosol extinction follows the wavelength power law of Angstrom, A.
    # (1929), "On the atmospheric transmission of sun radiation and on dust in
    # the air", Geografiska Annaler 11, 156-166.
    K_AEROSOL_BASE = 0.08

    # Scale height for pressure (km) - standard atmospheric density scale.
    SCALE_HEIGHT_PRESSURE = 8.5

    # Scale height for aerosols (km) - aerosols are concentrated in the lower
    # troposphere (order 1-2 km; e.g. Allen, "Astrophysical Quantities").
    SCALE_HEIGHT_AEROSOL = 1.5

    # Airglow brightness (nanoLamberts) - natural minimum sky glow. Standard
    # dark-night airglow level (Allen, "Astrophysical Quantities"; Roach, F.E.
    # & Gordon, J.L. (1973), "The Light of the Night Sky").
    AIRGLOW = 145.0

    # Zodiacal light brightness (nanoLamberts) - representative near-ecliptic
    # value (Allen, "Astrophysical Quantities"; Roach & Gordon 1973).
    ZODIACAL_LIGHT = 100.0

    # Dark sky brightness at zenith (mag/arcsec^2) - representative excellent
    # dark-site V-band value.
    DARK_SKY_MAG = 21.6

    # Full Moon brightness relative to Sun (magnitude difference). The Sun-Moon
    # magnitude offset near full phase is ~14 mag (Allen, "Astrophysical
    # Quantities").
    FULL_MOON_MAG_DIFF = 14.0

    # Visual limiting magnitude for perfect dark-sky conditions - classic
    # naked-eye limit for a young dark-adapted observer (Schaefer (1990),
    # PASP 102, 212).
    PERFECT_SKY_LIM_MAG = 6.5

    # Legacy descriptive defaults retained for callers inspecting this class.
    # The active arcus-visionis calculation below root-solves the Schaefer
    # limiting-magnitude model and does not read this mapping. These exact
    # object values are project defaults; they lie in the same range as the
    # classical Schoch/Ptolemaic arcus-visionis thresholds discussed in the
    # heliacal-visibility literature (e.g. Schaefer, B.E. (1987), "Heliacal
    # Rise Phenomena", Journal for the History of Astronomy 18,
    # Archaeoastronomy Suppl. 11, S19-S33) but are not a verbatim
    # transcription of any historical table.
    ARCUS_VISIONIS = {
        "venus": 5.0,
        "mercury": 10.0,
        "mars": 11.0,
        "jupiter": 9.0,
        "saturn": 11.0,
        "star_bright": 7.0,  # Stars brighter than mag 1
        "star_medium": 10.0,  # Stars mag 1-3
        "star_faint": 13.0,  # Stars fainter than mag 3
    }

    # Representative threshold-contrast constant for naked-eye point-source
    # detection from the Schaefer visibility framework (Schaefer (1990),
    # PASP 102, 212; Schaefer (1993), Vistas in Astronomy 36, 311). Legacy
    # descriptive value; not read by the active VISLIMIT path.
    THRESHOLD_CONTRAST = 0.0094

    # Eye pupil diameter in dark-adapted conditions (mm) - classical maximum
    # dark-adapted pupil; the active path instead uses the age-dependent
    # relation _vl_eye_pupil (Schaefer 1990/1998).
    DARK_PUPIL_DIAMETER = 7.0

    # Eye pupil diameter in bright (photopic) conditions (mm).
    BRIGHT_PUPIL_DIAMETER = 2.0


class AtmosphericConditions(NamedTuple):
    """Atmospheric conditions for visibility calculations."""

    pressure: float  # Atmospheric pressure (mbar)
    temperature: float  # Temperature (Celsius)
    humidity: float  # Relative humidity (0-100%)
    met_range: float  # Meteorological range (km), or ktot if < 1
    altitude: float  # Observer altitude (meters)


class ObserverParams(NamedTuple):
    """Observer parameters for visibility calculations."""

    age: float  # Observer age (years)
    snellen: float  # Snellen ratio (1.0 = normal)
    binocular: bool  # True if binocular, False if monocular
    telescope_mag: float  # Telescope magnification (0 = naked eye)
    aperture: float  # Telescope aperture (mm)
    transmission: float  # Optical transmission coefficient


# =============================================================================
# Schaefer VISLIMIT (1998) visual-limiting-magnitude model, V band
# =============================================================================
#
# Independent implementation of Bradley E. Schaefer's VISLIMIT algorithm
# ("To the Visual Limits", Sky & Telescope, May 1998; "Astronomy and the
# limits of vision", Vistas in Astronomy 36, 1993), restricted to the V
# photometric band (0.55 um). The extinction, airmass, twilight, moonlight
# and dark-night-sky components and the contrast-threshold (Hecht)
# magnitude conversion follow Schaefer's published algorithm. The three
# atmospheric scale heights are parameters of this library's visibility
# model (``_VL_SH_*`` below), and the dark-night-sky normalisation uses a
# published natural-sky luminance (Falchi et al. 2016, see below).
#
# V-band VISLIMIT constants (I = 3 in Schaefer's DATA statements):
_VL_MOI = -11.05  # extra-atmospheric reference magnitude constant
_VL_MS = -26.74  # Sun V magnitude constant
_VL_OZ = 0.031  # ozone coefficient
_VL_WT = 0.031  # water-vapour coefficient
# Falchi et al. (2016, Science Advances 2:e1600377) use a natural V-band
# background of 22.0 mag arcsec^-2 (174 micro-cd m^-2).  The photometric
# zero point below is the standard visual relation
# L[cd m^-2] = 108000 * 10**(-0.4 * mag_arcsec2).  One nanolambert is
# 1e-5/pi cd m^-2.  Keeping the source quantity in magnitudes makes the
# provenance and unit conversion explicit.
_VL_DARK_SKY_MAG_ARCSEC2 = 22.0
_VL_DARK_SKY_NL = (
    108000.0 * 10.0 ** (-0.4 * _VL_DARK_SKY_MAG_ARCSEC2) * math.pi / 1.0e-5
)
# Unit-scale factor converting nanoLamberts to the VISLIMIT model's internal
# radiance scale, on which the additive daylight/twilight/moonlight scattering
# terms below are expressed (their absolute coefficients such as 43.27 and the
# 6.2e7/10**6.15/10**5.36 scattering constants come from Schaefer's VISLIMIT
# algorithm, Sky & Telescope 95(5), 57). It cancels for the dark-sky term
# (which is multiplied and then divided by it) and fixes the relative weight of
# the scattered-light terms; the sum is divided back out to return nanoLamberts.
_VL_MODEL_RADIANCE_PER_NL = 1.11e-15
# Exponential scale heights (m) of the three extinction components that
# depend on the observer's elevation h: molecular (Rayleigh) scattering,
# aerosol scattering and water-vapour absorption. Each coefficient falls off
# as exp(-h/H) in _vl_extinction_components(). They are parameters of this
# library's visibility model, conserved for result stability; adopted
# values, source not verified.
_VL_SH_RAY = 8515.0
_VL_SH_AER = 3745.0
_VL_SH_WAT = 3000.0
_RD = math.pi / 180.0


def _vl_month_cont(jd: float) -> float:
    """Continuous month (1..13) of a Julian Day, for the seasonal terms."""
    Z = math.floor(jd + 0.5)
    F = jd + 0.5 - Z
    A = Z
    if Z >= 2299161:
        alpha = math.floor((Z - 1867216.25) / 36524.25)
        A = Z + 1 + alpha - math.floor(alpha / 4)
    B = A + 1524
    C = math.floor((B - 122.1) / 365.25)
    D = math.floor(365.25 * C)
    E = math.floor((B - D) / 30.6001)
    day = B - D - math.floor(30.6001 * E) + F
    m = E - 1 if E < 14 else E - 13
    return m + (day - 1.0) / 30.4375


def _vl_extinction_components(
    temperature: float,
    humidity_pct: float,
    altitude_m: float,
    latitude: float,
    jd: float,
    pressure_mbar: float = 0.0,
) -> Tuple[float, float, float, float]:
    """Return the four V-band extinction coefficients (KR, KA, KO, KW).

    Schaefer VISLIMIT extinction with this library's adopted scale heights
    (``_VL_SH_*``) and the seasonal ozone term. Coefficients are per unit
    airmass.  Aerosol extinction varies with humidity and altitude; no
    empirically recovered seasonal multiplier is applied.

    Rayleigh scattering is proportional to the molecular column mass, i.e.
    to the site pressure (Schaefer 1993): an explicit barometric pressure
    (``pressure_mbar > 0``) therefore scales the molecular coefficient as
    P/1013.25; with the sentinel 0 the standard-atmosphere altitude form
    ``exp(-h/H)`` is used, which is the same quantity for a standard site.
    """
    Mc = _vl_month_cont(jd)
    RA = (Mc - 3.0) * 30.0 * _RD
    LT = latitude * _RD
    RH = min(max(humidity_pct, 1.0), 99.5)
    if pressure_mbar > 0.0:
        kr = 0.1066 * (pressure_mbar / 1013.25)
    else:
        kr = 0.1066 * math.exp(-altitude_m / _VL_SH_RAY)
    ka = (
        0.1
        * math.exp(-altitude_m / _VL_SH_AER)
        * ((1.0 - 0.32 / math.log(RH / 100.0)) ** 1.33)
    )
    ko = _VL_OZ * (3.0 + 0.4 * (LT * math.cos(RA) - math.cos(3.0 * LT))) / 3.0
    kw = (
        _VL_WT
        * 0.94
        * (RH / 100.0)
        * math.exp(temperature / 15.0)
        * math.exp(-altitude_m / _VL_SH_WAT)
    )
    return kr, ka, ko, kw


def _vl_airmass(alt_deg: float) -> float:
    """Sky airmass vs altitude, Rozenberg (1966) low-altitude form.

    ``X = 1 / (cos z + 0.025 * exp(-11 cos z))``, valid to the horizon, from
    Rozenberg, G.V. (1966), "Twilight: A Study in Atmospheric Optics". This is
    the general sky airmass used by Schaefer's VISLIMIT model (NOT the Kasten &
    Young relation, which appears in the standalone extinction helper module).
    """
    if alt_deg <= 0.0:
        return 40.0
    c = math.cos((90.0 - alt_deg) * _RD)
    return 1.0 / (c + 0.025 * math.exp(-11.0 * c))


def _vl_airmass_gas(alt_deg: float) -> float:
    """Gas/Rayleigh relative airmass, Schaefer VISLIMIT (S&T 95(5), 57)."""
    c = math.cos((90.0 - alt_deg) * _RD)
    return 1.0 / (c + 0.0286 * math.exp(-10.5 * c))


def _vl_airmass_aer(alt_deg: float) -> float:
    """Aerosol relative airmass, Schaefer VISLIMIT (S&T 95(5), 57)."""
    c = math.cos((90.0 - alt_deg) * _RD)
    return 1.0 / (c + 0.0123 * math.exp(-24.5 * c))


def _vl_airmass_ozone(alt_deg: float) -> float:
    """Spherical-shell ozone airmass, Schaefer VISLIMIT (S&T 95(5), 57)."""
    s = math.sin((90.0 - alt_deg) * _RD)
    val = 1.0 - (s / (1.0 + 20.0 / 6378.0)) ** 2
    return 1.0 / math.sqrt(val) if val > 0 else 40.0


def _vl_hecht_threshold(bl: float) -> float:
    """Hecht/Knoll-Schaefer contrast threshold (foot-candles) for sky ``bl``.

    Two photoreceptor regimes (photopic above / scotopic below 1500 nL), exactly
    as in Schaefer's VISLIMIT algorithm (Schaefer (1998), Sky & Telescope
    95(5), 57): ``TH = C1 * (1 + sqrt(C2 * BL))**2``, with the C1/C2 pairs the
    published photopic and scotopic Hecht-model constants.
    """
    if bl > 1500.0:
        c1 = 10.0**-8.350001
        c2 = 10.0**-5.9
    else:
        c1 = 10.0**-9.8
        c2 = 10.0**-1.9
    return c1 * ((1.0 + math.sqrt(c2 * bl)) ** 2)


# =============================================================================
# Observer-dependent corrections (eye age + optical instrument)
# =============================================================================
#
# Observer corrections are derived from pupil area, optical throughput, exit
# pupil, and the Hecht contrast threshold.  They are deltas relative to the
# naked-eye default (age 36, Snellen 1, no instrument).
#
# Eye dark-adapted pupil diameter vs age (mm), from Schaefer's published
# visual/telescopic limiting-magnitude work (PASP 102, 212; Sky & Telescope,
# May 1998): DE = 7 * exp(-age**2 / 20000). The pupil is capped at its
# age-23 value (younger observers do not gain further dark-adapted diameter),
# following Schaefer's published observer model.
_VL_PUPIL_D0 = 7.0
_VL_PUPIL_K = 20000.0
_VL_PUPIL_AGE_MIN = 23.0
_VL_PUPIL_REF_AGE = 36.0


def _vl_eye_pupil(age: float) -> float:
    """Dark-adapted eye pupil diameter (mm) for an observer age (Schaefer)."""
    a = age if age > _VL_PUPIL_AGE_MIN else _VL_PUPIL_AGE_MIN
    return _VL_PUPIL_D0 * math.exp(-(a * a) / _VL_PUPIL_K)


# Reference pupil at age 36 (~6.56079 mm); the corrections are relative to it.
_VL_PUPIL_REF = _vl_eye_pupil(_VL_PUPIL_REF_AGE)

_NAKED_OBSERVER: Tuple[float, ...] = (36.0, 1.0, 0.0, 0.0, 0.0, 0.0)


def _parse_observer_optics(
    observer: Sequence[float], flags: int
) -> Tuple[float, float, bool, float, float, float]:
    """Parse an observer tuple into Schaefer-model parameters.

    Age (``observer[0]``) and Snellen ratio (``observer[1]``) are always read.
    The optical-instrument parameters (``observer[2..5]``: binocular flag,
    magnification, aperture in mm, transmission) are read ONLY when
    ``HELFLAG_OPTICAL_PARAMS`` is set in ``flags``; otherwise they are left at
    their naked-eye no-op defaults.

    Returns:
        (age, snellen, binocular, telescope_mag, aperture, transmission)
    """
    from .constants import HELFLAG_OPTICAL_PARAMS

    age = float(observer[0]) if len(observer) > 0 else 36.0
    snellen = float(observer[1]) if len(observer) > 1 else 1.0
    binocular = False
    telescope_mag = 0.0
    aperture = 0.0
    transmission = 1.0
    if flags & HELFLAG_OPTICAL_PARAMS:
        binocular = bool(observer[2]) if len(observer) > 2 else False
        telescope_mag = float(observer[3]) if len(observer) > 3 else 0.0
        aperture = float(observer[4]) if len(observer) > 4 else 0.0
        transmission = float(observer[5]) if len(observer) > 5 else 1.0
    return age, snellen, binocular, telescope_mag, aperture, transmission


class SchaeferModel:
    """
    Complete Schaefer atmospheric model for heliacal visibility.

    This class implements the Schaefer (1990) model for calculating:
    - Atmospheric extinction
    - Sky brightness
    - Limiting visual magnitude
    - Visibility thresholds

    The implementation follows the cited Schaefer visibility model and public
    result layout.
    """

    def __init__(
        self,
        atmo: Optional[AtmosphericConditions] = None,
        observer: Optional[ObserverParams] = None,
        latitude: float = 0.0,
        jd: float = 2451545.0,
        pressure_scales_extinction: bool = False,
    ):
        """
        Initialize the Schaefer model with atmospheric and observer conditions.

        Args:
            atmo: Atmospheric conditions (defaults to standard atmosphere)
            observer: Observer parameters (defaults to standard observer)
            latitude: Observer latitude in degrees (for the seasonal/latitude
                extinction terms of Schaefer's VISLIMIT model).
            jd: Julian Day (UT) used to derive the season (month) for the
                seasonal aerosol/ozone extinction terms.
        """
        if atmo is None:
            atmo = AtmosphericConditions(
                pressure=1013.25,
                temperature=15.0,
                humidity=50.0,
                met_range=0.0,
                altitude=0.0,
            )
        if observer is None:
            observer = ObserverParams(
                age=36.0,
                snellen=1.0,
                binocular=False,
                telescope_mag=0.0,
                aperture=0.0,
                transmission=1.0,
            )

        self.atmo = atmo

        # The published VISLIMIT Rayleigh form is the standard-atmosphere
        # altitude law 0.1066*exp(-h/H) (Schaefer, S&T 95(5), 57; Schaefer
        # 1993) — pressure-independent by construction. Scaling with
        # P/1013.25 is the equivalent parametrization when barometric
        # pressure is an explicit observation, used by the
        # limiting-magnitude surface. The constructor toggle selects the
        # parametrization for this model instance; the sentinel 0.0 selects
        # the published altitude form, which the heliacal model uses.
        self._extinction_pressure = atmo.pressure if pressure_scales_extinction else 0.0
        self.observer = observer
        self.latitude = latitude
        self.jd = jd

        # Precompute extinction coefficients
        self._compute_extinction_coefficients()

    def _compute_extinction_coefficients(self) -> None:
        """Compute Schaefer VISLIMIT (V-band) extinction coefficients.

        Uses the observer temperature/humidity/altitude/latitude and the
        season (from ``self.jd``) exactly as in Schaefer's VISLIMIT model
        with Reijs's refined scale heights. ``met_range`` (datm[3]) is
        honoured as an override: 0 < v < 1 sets the total coefficient
        directly (ktot); v >= 1 is a meteorological visual range in km
        whose aerosol contribution follows Koschmieder's relation.
        """
        kr, ka, ko, kw = _vl_extinction_components(
            self.atmo.temperature,
            self.atmo.humidity,
            self.atmo.altitude,
            self.latitude,
            self.jd,
            pressure_mbar=self._extinction_pressure,
        )
        mr = self.atmo.met_range
        if 0.0 < mr < 1.0:
            # Direct total-extinction (ktot) override. Distribute the
            # requested total across the components in the model's own
            # proportions so the per-component airmass weighting still
            # applies; clamp to the molecular+ozone+water floor.
            floor = kr + ko + kw
            k_total = max(mr, floor)
            ka = max(0.0, k_total - floor)
        elif mr >= 1.0:
            # Meteorological visual range V (km) -> aerosol coefficient via
            # the Koschmieder (1924) relation, k = 3.912 / V, with the
            # molecular part kr removed. Koschmieder, H. (1924), "Theorie der
            # horizontalen Sichtweite", Beitr. Phys. freien Atmos. 12, 33-53.
            ka = max(0.0, 3.912 / mr - kr)
        self.k_rayleigh = kr
        self.k_aerosol = ka
        self.k_ozone = ko
        self.k_water = kw
        self.k_total = kr + ka + ko + kw

    def update_season(self, jd: float) -> None:
        """Refresh the extinction for a new date (season) if it changed.

        The aerosol/ozone extinction terms depend on the season (month).
        A heliacal search spans many days, so the visibility check must
        use the extinction of the evaluated day, not of the search start.
        """
        if abs(jd - self.jd) < 3.0:
            return
        self.jd = jd
        self._compute_extinction_coefficients()

    def airmass(self, altitude_deg: float) -> float:
        """Schaefer VISLIMIT sky airmass (gas component) at an altitude."""
        return _vl_airmass(altitude_deg)

    def extinction(self, altitude_deg: float) -> float:
        """
        Total atmospheric extinction in magnitudes at an altitude.

        Uses Schaefer's per-component airmasses (gas, aerosol, ozone,
        water) so the near-horizon behaviour matches VISLIMIT:
        DM = KR*XG + KA*XA + KO*XO + KW*XG.
        """
        xg = _vl_airmass_gas(altitude_deg)
        xa = _vl_airmass_aer(altitude_deg)
        xo = _vl_airmass_ozone(altitude_deg)
        return (
            self.k_rayleigh * xg
            + self.k_aerosol * xa
            + self.k_ozone * xo
            + self.k_water * xg
        )

    def _sky_brightness_bl(
        self,
        sun_alt: float,
        moon_alt: float,
        moon_phase: float,
        obj_alt: float,
        sun_obj_angle: float,
        moon_obj_angle: float,
    ) -> float:
        """
        Total sky brightness at the object, in nanoLamberts (Schaefer VISLIMIT).

        Sums the dark night sky (BN), the twilight OR daylight component
        (BT/BD, whichever is brighter, with the daylight term gated to a
        Sun above the horizon), and the moonlight component (BM).

        Args:
            sun_alt: Sun altitude in degrees.
            moon_alt: Moon altitude in degrees.
            moon_phase: Moon illuminated fraction (0 = new, 1 = full).
            obj_alt: Object altitude in degrees.
            sun_obj_angle: Sky angular separation object-Sun in degrees.
            moon_obj_angle: Sky angular separation object-Moon in degrees.

        Returns:
            Total sky brightness in nanoLamberts.
        """
        K = self.k_total
        Z = 90.0 - obj_alt
        X = _vl_airmass(obj_alt)
        ZZ = Z * _RD
        year = 2000.0 + (self.jd - 2451545.0) / 365.25

        # Dark night sky (Schaefer VISLIMIT dark-sky term, Schaefer (1998),
        # Sky & Telescope 95(5), 57): a 0.3-amplitude 11-year solar-cycle
        # modulation of airglow, a van Rhijn zenith-distance brightening
        # (0.4 + 0.6/sqrt(1 - 0.96 sin^2 Z)), and extinction 10**(-0.4 K X)
        # to the object altitude. The absolute level is set by the published
        # natural dark-sky luminance in _VL_DARK_SKY_NL.
        bn = _VL_DARK_SKY_NL * _VL_MODEL_RADIANCE_PER_NL
        bn *= 1.0 + 0.3 * math.cos(6.283 * (year - 1992.0) / 11.0)
        bn *= 0.4 + 0.6 / math.sqrt(1.0 - 0.96 * (math.sin(ZZ)) ** 2)
        bn *= 10.0 ** (-0.4 * K * X)
        b = bn

        one_minus = 1.0 - 10.0 ** (-0.4 * K * X)
        rs = max(sun_obj_angle, 1.0)
        if sun_alt > 0.0:
            # Daylight scattering (Sun above the horizon). The twilight term
            # BT is only valid for a Sun below the horizon (its -HS factor
            # diverges for HS > 0), so daylight uses BD alone.
            # fs is the angular scattering function f(rho) of Krisciunas, K. &
            # Schaefer, B.E. (1991), PASP 103, 1033 (their scattered-light
            # relations, eqs. 15-21): an aureole/forward term (6.2e7*rho^-2),
            # an intermediate-angle term (10**(6.15 - rho/40)), and a
            # Rayleigh+Mie term (10**5.36*(1.06 + cos^2 rho)), with rho the
            # Sun-object sky angle in degrees.
            xs = _vl_airmass(sun_alt)
            c4 = 10.0 ** (-0.4 * K * xs)
            fs = (
                6.2e7 * (rs**-2)
                + 10.0 ** (6.15 - rs / 40.0)
                + (10.0**5.36) * (1.06 + (math.cos(rs * _RD)) ** 2)
            )
            b += (
                10.0 ** (-0.4 * (_VL_MS - _VL_MOI + 43.27))
                * one_minus
                * (fs * c4 + 440000.0 * (1.0 - c4))
            )
        else:
            # Twilight brightness (Sun below the horizon), from Schaefer's
            # VISLIMIT twilight term (Schaefer (1998), Sky & Telescope 95(5),
            # 57): the 32.5 offset and the -HS / -(Z/(360 K)) dependence on
            # Sun depression and object zenith distance are the published
            # BT relation.
            bt = 10.0 ** (
                -0.4 * (_VL_MS - _VL_MOI + 32.5 - sun_alt - (Z / (360.0 * K)))
            )
            b += bt * (100.0 / rs) * one_minus

        # Moonlight, following Krisciunas, K. & Schaefer, B.E. (1991),
        # "A model of the brightness of moonlight", PASP 103, 1033-1039,
        # DOI 10.1086/132921.
        if moon_alt > 0.0:
            rm = max(moon_obj_angle, 1.0)
            am = math.degrees(math.acos(min(1.0, max(-1.0, 2.0 * moon_phase - 1.0))))
            # Moon V magnitude as a function of phase angle am (deg): their
            # eq. 9, m = -12.73 + 0.026*|am| + 4e-9*am**4.
            mm = -12.73 + 0.026 * abs(am) + 4e-9 * (am**4)
            xm = _vl_airmass(moon_alt)
            c3 = 10.0 ** (-0.4 * K * xm)
            # fm is the same K&S (1991) scattering function f(rho) as fs above
            # (eqs. 15-21), evaluated at the Moon-object sky angle rm.
            fm = (
                6.2e7 * (rm**-2)
                + 10.0 ** (6.15 - rm / 40.0)
                + (10.0**5.36) * (1.06 + (math.cos(rm * _RD)) ** 2)
            )
            bm = (
                10.0 ** (-0.4 * (mm - _VL_MOI + 43.27))
                * one_minus
                * (fm * c3 + 440000.0 * (1.0 - c3))
            )
            b += bm

        return b / _VL_MODEL_RADIANCE_PER_NL

    def sky_brightness_bl(
        self,
        sun_alt: float,
        moon_alt: float,
        moon_phase: float,
        obj_alt: float,
        sun_obj_angle: float,
        moon_obj_angle: float,
    ) -> float:
        """Total sky brightness in nanoLamberts (Schaefer VISLIMIT)."""
        return self._sky_brightness_bl(
            sun_alt, moon_alt, moon_phase, obj_alt, sun_obj_angle, moon_obj_angle
        )

    def limiting_magnitude(
        self,
        sun_alt: float,
        moon_alt: float = -90.0,
        moon_phase: float = 0.0,
        obj_alt: float = 45.0,
        sun_obj_angle: float = 180.0,
        moon_obj_angle: float = 180.0,
    ) -> float:
        """
        Calculate the limiting visual magnitude (Schaefer VISLIMIT).

        Returns the faintest *catalog* magnitude visible at the object's
        position and the given sky conditions. Atmospheric extinction to
        the object's altitude is already folded in (via the -DM term), so
        this may be compared directly against a body's catalog magnitude
        as required by the public ``dret[0]`` contract.

        Args:
            sun_alt: Sun altitude in degrees
            moon_alt: Moon altitude in degrees
            moon_phase: Moon illuminated fraction (0 = new, 1 = full)
            obj_alt: Object altitude in degrees
            sun_obj_angle: Sky angular separation object-Sun in degrees
            moon_obj_angle: Sky angular separation object-Moon in degrees

        Returns:
            Limiting visual magnitude (fainter = larger number)
        """
        bl = self._sky_brightness_bl(
            sun_alt, moon_alt, moon_phase, obj_alt, sun_obj_angle, moon_obj_angle
        )
        # Contrast threshold (Hecht/Knoll-Schaefer): photopic vs scotopic.
        th = _vl_hecht_threshold(bl)  # foot-candles

        # Extinction to the object (magnitudes) is included in the limit.
        # The -16.57 magnitude zero point converts the Hecht/Knoll-Schaefer
        # threshold illuminance (foot-candles) into a V-band limiting
        # magnitude, per Schaefer's VISLIMIT algorithm (Schaefer (1998),
        # Sky & Telescope 95(5), 57).
        dm = self.extinction(obj_alt)
        m_lim = -16.57 - 2.5 * math.log10(th) - dm

        # Snellen ratio (better eyes see fainter); VISLIMIT uses 5 log10(SN).
        # Always applied and exactly zero at SN = 1 (naked-eye default).
        if self.observer.snellen > 0:
            m_lim += 5.0 * math.log10(self.observer.snellen)

        # Observer age / optical-aid correction. Both are deltas that are
        # exactly zero for the naked-eye default (age 36, no instrument).
        m_lim += self._observer_correction(bl)

        return m_lim

    def _observer_correction(self, bl: float) -> float:
        """Observer age OR optical-instrument delta (magnitudes) at sky ``bl``.

        When an optical instrument is used, its effective aperture includes the
        observer's age-dependent pupil.  Otherwise the pupil-area age term is
        applied directly.  The Snellen term is handled separately.
        """
        obs = self.observer
        if obs.telescope_mag > 1.0 and obs.aperture > 0.0:
            return self._optics_delta(bl)
        return self._age_delta(bl)

    def _age_delta(self, _bl: float) -> float:
        """Pupil-area magnitude delta (zero at the age-36 reference pupil)."""
        r = _vl_eye_pupil(self.observer.age) / _VL_PUPIL_REF
        return 5.0 * math.log10(r)

    def _optics_delta(self, bl: float) -> float:
        """Limiting-magnitude delta from an optical instrument (zero for none).

        Uses the physical point-source light grasp and exit-pupil relations from
        Schaefer (1990, PASP 102, 212).  The entrance aperture admitted by the
        eye is ``min(D, M * eye_pupil)``.  Throughput scales source photons and
        retinal sky luminance; binocular viewing receives the usual ``sqrt(2)``
        independent-eye signal-to-noise gain.
        """
        obs = self.observer
        d = obs.aperture
        mg = obs.telescope_mag
        t = max(obs.transmission, 1.0e-12)
        dp = d / mg  # exit pupil (mm)
        eye_pupil = _vl_eye_pupil(obs.age)
        effective_aperture = min(d, mg * eye_pupil)
        source_ratio = t * (effective_aperture / _VL_PUPIL_REF) ** 2
        if obs.binocular:
            source_ratio *= math.sqrt(2.0)
        retinal_pupil = min(dp, eye_pupil)
        background_ratio = t * (retinal_pupil / _VL_PUPIL_REF) ** 2
        th0 = _vl_hecht_threshold(bl)
        thr = _vl_hecht_threshold(bl * background_ratio)
        return 2.5 * math.log10(source_ratio) - 2.5 * math.log10(thr / th0)

    def arcus_visionis_required(
        self, body_mag: float, body_type: str = "planet"
    ) -> float:
        """
        Calculate the required arcus visionis for visibility.

        The arcus visionis is the minimum altitude difference between
        a celestial body and the Sun for the body to be visible.

        The historical term *arcus visionis* describes the requested angular
        quantity. The numeric result is independently root-solved over the
        cited Schaefer limiting-magnitude model; no Ptolemaic coefficient
        table is used.

        Args:
            body_mag: Visual magnitude of the body
            body_type: Type of body ("planet" or "star")

        Returns:
            Required arcus visionis in degrees
        """

        # Root-find on the VISLIMIT limiting-magnitude model. For a body of
        # magnitude ``body_mag`` placed directly above the Sun (azimuth
        # aligned, so the Sun-object sky angle equals the arcus visionis),
        # find, at each Sun depression, the object altitude at which the
        # object is exactly at the visibility limit (limiting mag == body
        # mag), then take the minimum required altitude difference. This
        # is the topocentric arcus-visionis construction.
        def _obj_alt_at_limit(sun_alt: float) -> Optional[float]:
            lo, hi = 0.2, 60.0
            f_lo = (
                self.limiting_magnitude(sun_alt, -90.0, 0.0, lo, lo - sun_alt, 180.0)
                - body_mag
            )
            f_hi = (
                self.limiting_magnitude(sun_alt, -90.0, 0.0, hi, hi - sun_alt, 180.0)
                - body_mag
            )
            if f_lo > 0:  # already visible at the horizon
                return lo
            if f_hi < 0:  # never visible up to 60 deg
                return None
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                fm = (
                    self.limiting_magnitude(
                        sun_alt, -90.0, 0.0, mid, mid - sun_alt, 180.0
                    )
                    - body_mag
                )
                if fm >= 0:
                    hi = mid
                else:
                    lo = mid
                if hi - lo < 1e-3:
                    break
            return 0.5 * (lo + hi)

        best_av = None
        sa = -1.0
        while sa >= -22.0:
            h = _obj_alt_at_limit(sa)
            if h is not None:
                av = h - sa
                if best_av is None or av < best_av:
                    best_av = av
            sa -= 0.5
        if best_av is None:
            return 0.0
        return best_av

    def is_visible(
        self,
        body_alt: float,
        body_mag: float,
        sun_alt: float,
        elongation: float,
        moon_alt: float = -90.0,
        moon_phase: float = 0.0,
        moon_obj_angle: float = 180.0,
        margin: float = 0.5,
    ) -> bool:
        """
        Determine if a celestial body is visible.

        Uses a Schaefer-style limiting-magnitude model. The body is
        visible when its apparent magnitude (catalog + extinction)
        is brighter than the sky's limiting apparent magnitude.

        A minimum elongation check prevents false positives when
        the body is geometrically above the horizon but lost in
        the Sun's glare.

        Args:
            body_alt: Body altitude in degrees
            body_mag: Body visual magnitude (catalog)
            sun_alt: Sun altitude in degrees
            elongation: Elongation from Sun in degrees
            moon_alt: Moon altitude in degrees
            moon_phase: Moon phase (0 = new, 1 = full)
            moon_obj_angle: Angular separation from Moon in degrees
            margin: Detection margin in magnitudes. Higher values
                require the body to be brighter relative to the
                limiting magnitude. Default 0.5 corresponds to
                ~90% detection probability.

        Returns:
            True if body is visible, False otherwise
        """
        # Body must be above horizon
        if body_alt < 0:
            return False

        # Sun must be below horizon for heliacal visibility
        if sun_alt > 0:
            return False

        # Minimum elongation: below ~7° even the brightest objects
        # are lost in the Sun's glare
        if elongation < 7.0:
            return False

        # Limiting (catalog) magnitude at the body position; this already
        # folds in atmospheric extinction to the body's altitude, so it is
        # compared directly against the body's catalog magnitude.
        lim_mag = self.limiting_magnitude(
            sun_alt=sun_alt,
            moon_alt=moon_alt,
            moon_phase=moon_phase,
            obj_alt=body_alt,
            sun_obj_angle=elongation,
            moon_obj_angle=moon_obj_angle,
        )

        # Body is visible if its catalog magnitude is brighter than the
        # limiting magnitude minus the detection margin.
        return body_mag <= lim_mag - margin

    def heliacal_visibility_angle(
        self,
        body_mag: float,
        sun_alt: float = -8.0,
    ) -> float:
        """
        Calculate the required body altitude for heliacal visibility.

        This is the altitude above horizon at which a body of given
        magnitude becomes visible, given the Sun's altitude.

        Args:
            body_mag: Body visual magnitude
            sun_alt: Sun altitude in degrees

        Returns:
            Required body altitude in degrees
        """
        # Required arcus visionis
        av_required = self.arcus_visionis_required(body_mag)

        # Body altitude = Sun altitude + arcus visionis
        body_alt = sun_alt + av_required

        # Minimum altitude above horizon
        return max(body_alt, 0.0)


def create_schaefer_model(
    pressure: float = 1013.25,
    temperature: float = 15.0,
    humidity: float = 50.0,
    met_range: float = 0.0,
    altitude: float = 0.0,
    observer_age: float = 36.0,
    snellen: float = 1.0,
    binocular: bool = False,
    telescope_mag: float = 0.0,
    aperture: float = 0.0,
    transmission: float = 1.0,
    latitude: float = 0.0,
    jd: float = 2451545.0,
    pressure_scales_extinction: bool = False,
) -> SchaeferModel:
    """
    Create a SchaeferModel instance with given parameters.

    Convenience function for creating the atmospheric model.

    Args:
        pressure: Atmospheric pressure in mbar (default 1013.25)
        temperature: Temperature in Celsius (default 15.0)
        humidity: Relative humidity in percent (default 50.0)
        met_range: Meteorological range in km, or ktot if < 1 (default 0.0)
        altitude: Observer altitude in meters (default 0.0)
        observer_age: Observer age in years (default 36.0)
        snellen: Snellen ratio, 1.0 = normal vision (default 1.0)
        binocular: True for binocular (two-eyed) viewing through the instrument
            (default False). No-op unless an instrument is given.
        telescope_mag: Optical magnification (0 = naked eye, default 0.0).
        aperture: Optical aperture / objective diameter in mm (default 0.0).
        transmission: Optical transmission coefficient, 0-1 (default 1.0).
        latitude: Observer latitude in degrees, for the seasonal/latitude
            extinction terms (default 0.0)
        jd: Julian Day (UT) for the season of the extinction terms
            (default J2000.0)

    Returns:
        SchaeferModel instance configured with given parameters
    """
    atmo = AtmosphericConditions(
        pressure=pressure,
        temperature=temperature,
        humidity=humidity,
        met_range=met_range,
        altitude=altitude,
    )
    observer = ObserverParams(
        age=observer_age,
        snellen=snellen,
        binocular=binocular,
        telescope_mag=telescope_mag,
        aperture=aperture,
        transmission=transmission,
    )
    return SchaeferModel(
        atmo,
        observer,
        latitude=latitude,
        jd=jd,
        pressure_scales_extinction=pressure_scales_extinction,
    )


# Outer planets (orbit outside Earth's orbit)
# These only have one type of conjunction (behind the Sun)
# EVENING_FIRST and MORNING_LAST are not applicable to these


def is_fixed_star(body: int) -> bool:
    """
    Check if a body ID corresponds to a fixed star.

    Args:
        body: Body ID constant

    Returns:
        True if the body is a fixed star (ID >= FIXSTAR_OFFSET), False otherwise
    """
    return body >= FIXSTAR_OFFSET


def _get_star_magnitude(star_id: int) -> float:
    """
    Get the visual magnitude of a fixed star from the catalog.

    Args:
        star_id: Star ID (body constant, e.g., SIRIUS)

    Returns:
        Visual magnitude of the star. Brighter stars have lower/negative values.
        Returns 6.0 (faint) if star not found.
    """
    from .fixed_stars import STAR_CATALOG

    for entry in STAR_CATALOG:
        if entry.id == star_id:
            return entry.magnitude
    return 6.0  # Default to faint star if not found


def is_inner_planet(body: int) -> bool:
    """
    Check if a body is an inner planet (Mercury or Venus).

    Inner planets orbit inside Earth's orbit and have both inferior
    and superior conjunctions with the Sun. They can appear as both
    morning and evening stars.

    Args:
        body: Planet ID constant (MERCURY, VENUS, etc.)

    Returns:
        True if the body is Mercury or Venus, False otherwise
    """
    return body in INNER_PLANETS


# =============================================================================
# LEB (Binary Ephemeris) paths for heliacal functions
# =============================================================================
#
# These private functions replicate the Skyfield-based heliacal logic
# using the LEB fast_calc pipeline.  They are called from the three public
# entry-points (_heliacal_ut_pythonic, _heliacal_pheno_ut_pythonic, vis_limit_mag) when a
# LEBReader is available, avoiding the 3 GB DE-kernel load entirely.
# =============================================================================


def _star_name_from_id(star_id: int) -> str:
    """Return the canonical star name for a star body ID."""
    from .fixed_stars import STAR_CATALOG

    for entry in STAR_CATALOG:
        if entry.id == star_id:
            return entry.name
    raise Error(f"unknown fixed star id {star_id}")


def _altaz_from_radec(
    ra_hours: float,
    dec_deg: float,
    lon_deg: float,
    lat_deg: float,
    gast_hours: float,
) -> tuple[float, float]:
    """Altitude/azimuth from of-date equatorial coordinates via the hour angle.

    Standard spherical-astronomy reduction (Explanatory Supplement to the
    Astronomical Almanac; Meeus, Astronomical Algorithms, ch. 13). The azimuth
    is returned in the north-zero, eastward-positive convention to match
    Skyfield's ``altaz()`` and the public heliacal azimuth slots.

    Args:
        ra_hours: Right ascension of date, in hours.
        dec_deg: Declination of date, in degrees.
        lon_deg: Observer geographic longitude (east positive), degrees.
        lat_deg: Observer geographic latitude, degrees.
        gast_hours: Greenwich apparent sidereal time, in hours.

    Returns:
        ``(altitude_deg, azimuth_deg)``; altitude is geometric (unrefracted).
    """
    lst_hours = gast_hours + lon_deg / 15.0
    ha_rad = math.radians((lst_hours - ra_hours) * 15.0)
    dec_rad = math.radians(dec_deg)
    lat_rad = math.radians(lat_deg)
    sin_alt = math.sin(lat_rad) * math.sin(dec_rad) + math.cos(lat_rad) * math.cos(
        dec_rad
    ) * math.cos(ha_rad)
    alt_deg = math.degrees(math.asin(max(-1.0, min(1.0, sin_alt))))
    # Azimuth from South (classic hour-angle form), rotated by 180 deg to the
    # north-zero convention.
    y = math.sin(ha_rad) * math.cos(dec_rad)
    x = math.cos(ha_rad) * math.sin(lat_rad) * math.cos(dec_rad) - math.sin(
        dec_rad
    ) * math.cos(lat_rad)
    az_deg = (math.degrees(math.atan2(y, x)) + 180.0) % 360.0
    return alt_deg, az_deg


def _leb_body_altaz(
    reader,
    jd_ut: float,
    body_id: int,
    geopos_lonlat: tuple,
    atpress: float,
    attemp: float,
    is_star: bool = False,
    star_name: str | None = None,
):
    """Compute horizontal coordinates (az, true_alt, app_alt) for a body via LEB.

    Args:
        reader: LEBReader instance.
        jd_ut: Julian Day UT1.
        body_id: Body constant.
        geopos_lonlat: (lon_deg, lat_deg, alt_m) for the observer.
        atpress: Atmospheric pressure in mbar.
        attemp: Atmospheric temperature in Celsius.
        is_star: True when computing a fixed star position.
        star_name: Star name string (required when is_star is True).

    Returns:
        (azimuth, true_altitude, apparent_altitude) in degrees.
    """
    from .utils import azalt, ECL2HOR

    if is_star:
        from .fixed_stars import fixstar_ut

        assert star_name is not None
        pos, _name, _flags = fixstar_ut(star_name, jd_ut, FLG_SPEED)
        ecl_lon, ecl_lat, ecl_dist = pos[0], pos[1], pos[2]
    else:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat

        jd_tt = jd_ut + deltat(jd_ut)
        topo = _topo_ecliptic(reader, jd_tt, jd_ut, body_id, geopos_lonlat)
        ecl_lon, ecl_lat, ecl_dist = topo[0], topo[1], topo[2]

    az, alt_true, alt_app = azalt(
        jd_ut,
        ECL2HOR,
        geopos_lonlat,
        atpress,
        attemp,
        (ecl_lon, ecl_lat, ecl_dist),
    )
    # Convert the public south-zero, westward azimuth convention to
    # Skyfield's north-zero, eastward convention.
    az_north = (az + 180.0) % 360.0
    return az_north, alt_true, alt_app


def _leb_ecliptic_pos(
    reader,
    jd_ut: float,
    body_id: int,
    geopos_lonlat: tuple,
    is_star: bool = False,
    star_name: str | None = None,
):
    """Return ecliptic (lon, lat, dist) for a body via LEB.

    For planets this is topocentric ecliptic of date; for stars it
    uses the already-LEB-ported fixstar_ut.
    """
    if is_star:
        from .fixed_stars import fixstar_ut

        assert star_name is not None
        pos, _name, _flags = fixstar_ut(star_name, jd_ut, FLG_SPEED)
        return pos[0], pos[1], pos[2]
    else:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat

        jd_tt = jd_ut + deltat(jd_ut)
        topo = _topo_ecliptic(reader, jd_tt, jd_ut, body_id, geopos_lonlat)
        return topo[0], topo[1], topo[2]


def _heliacal_ut_leb(
    reader,
    jd_start: float,
    lat: float,
    lon: float,
    altitude: float,
    pressure: float,
    temperature: float,
    humidity: float,
    body: int,
    event_type: int,
    flags: int,
    met_range: float = 0.0,
    observer: Sequence[float] = _NAKED_OBSERVER,
) -> tuple:
    """LEB-backed implementation of _heliacal_ut_pythonic().

    Mirrors the Skyfield version but obtains all positions via
    _topo_ecliptic / fixstar_ut / azalt / angular_separation.
    """
    from .constants import (
        HELIACAL_RISING,
        HELIACAL_SETTING,
        EVENING_FIRST,
        MORNING_LAST,
    )
    from .planets import _PLANET_MAP, pheno_ut
    from .utils import angular_separation

    # --- validation (same as Skyfield path) ---
    if event_type not in (
        HELIACAL_RISING,
        HELIACAL_SETTING,
        EVENING_FIRST,
        MORNING_LAST,
    ):
        raise Error(
            f"invalid event type {event_type}; use HELIACAL_RISING, "
            "HELIACAL_SETTING, EVENING_FIRST or MORNING_LAST"
        )
    if event_type in (EVENING_FIRST, MORNING_LAST):
        if is_fixed_star(body):
            raise Error(
                "EVENING_FIRST and MORNING_LAST do not apply to fixed stars; "
                "use HELIACAL_RISING or HELIACAL_SETTING"
            )
        if not is_inner_planet(body) and body != MOON:
            raise Error(
                "EVENING_FIRST and MORNING_LAST apply only to the inner planets "
                "(Mercury, Venus) and the Moon; for the outer planets use "
                "HELIACAL_RISING or HELIACAL_SETTING"
            )
    if body == SUN:
        raise Error("the Sun is not a valid object for heliacal calculations")
    if body == MOON and event_type in (HELIACAL_RISING, HELIACAL_SETTING):
        raise Error(
            "heliacal rising and setting (event types 1 and 2) are not defined "
            "for the Moon"
        )

    is_star = is_fixed_star(body)
    if not is_star and body not in _PLANET_MAP:
        raise Error(f"body id {body} is neither a planet nor a fixed star")

    star_name: str | None = None
    star_magnitude = 0.0
    if is_star:
        star_name = _star_name_from_id(body)
        star_magnitude = _get_star_magnitude(body)

    # Observer geopos in (lon, lat, alt) order for azalt / _topo_ecliptic
    geopos = (lon, lat, altitude)

    # --- inner helpers (closures over reader / geopos / etc.) ---

    def _get_altitudes(jd: float):
        """Return (sun_alt, body_alt, body_az)."""
        _, sun_alt, _ = _leb_body_altaz(
            reader,
            jd,
            SUN,
            geopos,
            pressure,
            temperature,
        )
        az, body_alt, _ = _leb_body_altaz(
            reader,
            jd,
            body,
            geopos,
            pressure,
            temperature,
            is_star=is_star,
            star_name=star_name,
        )
        return sun_alt, body_alt, az

    def _get_sun_alt(jd: float) -> float:
        """Sun altitude only (avoids the body position for twilight scans)."""
        _, sun_alt, _ = _leb_body_altaz(
            reader,
            jd,
            SUN,
            geopos,
            pressure,
            temperature,
        )
        return sun_alt

    def _get_elongation(jd: float) -> float:
        if not is_star:
            try:
                pheno = pheno_ut(jd, body, _heliacal_eph_flags(flags))
                return pheno[2]
            except (ValueError, TypeError, ArithmeticError):
                pass
        # Manual elongation from ecliptic coords
        sun_ecl = _leb_ecliptic_pos(reader, jd, SUN, geopos)
        body_ecl = _leb_ecliptic_pos(
            reader,
            jd,
            body,
            geopos,
            is_star=is_star,
            star_name=star_name,
        )
        return angular_separation(sun_ecl[0], sun_ecl[1], body_ecl[0], body_ecl[1])

    def _sun_body_lon_diff(jd: float) -> float:
        """Signed ecliptic longitude of the body minus the Sun, in (-180, 180].

        Negative means the body is WEST of the Sun (rises before it -- a
        morning apparition); positive means EAST (sets after it -- an evening
        apparition). Mirrors the Skyfield twin so the two backends reject the
        same wrong-apparition twilight detections (see that helper).
        """
        sun_ecl = _leb_ecliptic_pos(reader, jd, SUN, geopos)
        body_ecl = _leb_ecliptic_pos(
            reader,
            jd,
            body,
            geopos,
            is_star=is_star,
            star_name=star_name,
        )
        return (body_ecl[0] - sun_ecl[0] + 180.0) % 360.0 - 180.0

    def _event_side_ok(jd: float, morning: bool) -> bool:
        """Whether a visibility at ``jd`` is on the side the event requires.

        Morning events require the body WEST of the Sun; evening events require
        it EAST. Structural mirror of the Skyfield twin.
        """
        lon_diff = _sun_body_lon_diff(jd)
        return (lon_diff < 0.0) if morning else (lon_diff > 0.0)

    def _get_body_magnitude(jd: float) -> float:
        if is_star:
            return star_magnitude
        try:
            pheno = pheno_ut(jd, body, _heliacal_eph_flags(flags))
            return pheno[4]
        except (ValueError, TypeError, ArithmeticError):
            return 0.0

    (
        _obs_age,
        _obs_snellen,
        _obs_bino,
        _obs_mag,
        _obs_aper,
        _obs_trans,
    ) = _parse_observer_optics(observer, flags)
    schaefer = create_schaefer_model(
        pressure=pressure,
        temperature=temperature,
        humidity=humidity * 100.0 if humidity <= 1.0 else humidity,
        altitude=altitude,
        met_range=met_range,
        observer_age=_obs_age,
        snellen=_obs_snellen,
        binocular=_obs_bino,
        telescope_mag=_obs_mag,
        aperture=_obs_aper,
        transmission=_obs_trans,
        latitude=lat,
        jd=jd_start,
    )

    def _get_moon_data(jd: float):
        """Return (moon_alt, phase, moon_body_sep)."""
        _, moon_alt, _ = _leb_body_altaz(
            reader,
            jd,
            MOON,
            geopos,
            pressure,
            temperature,
        )
        # Moon phase from geocentric elongation (matching Skyfield path)
        from .planets import calc_ut as _scu_md

        _sun_geo, _ = _scu_md(jd, SUN, FLG_SPEED)
        _moon_geo, _ = _scu_md(jd, MOON, FLG_SPEED)
        elong_moon = angular_separation(
            _sun_geo[0],
            _sun_geo[1],
            _moon_geo[0],
            _moon_geo[1],
        )
        phase = (1.0 - math.cos(math.radians(elong_moon))) / 2.0
        moon_ecl = _leb_ecliptic_pos(reader, jd, MOON, geopos)

        body_ecl = _leb_ecliptic_pos(
            reader,
            jd,
            body,
            geopos,
            is_star=is_star,
            star_name=star_name,
        )
        moon_body_sep = angular_separation(
            moon_ecl[0],
            moon_ecl[1],
            body_ecl[0],
            body_ecl[1],
        )
        return moon_alt, phase, moon_body_sep

    def _is_body_visible(jd: float):
        sun_alt, body_alt, _body_az = _get_altitudes(jd)
        elongation = _get_elongation(jd)
        if body_alt < 0:
            return False, sun_alt, body_alt, elongation
        if sun_alt > 0:
            return False, sun_alt, body_alt, elongation
        body_mag = _get_body_magnitude(jd)
        moon_alt, moon_phase, moon_body_sep = _get_moon_data(jd)
        schaefer.update_season(jd)
        visible = schaefer.is_visible(
            body_alt=body_alt,
            body_mag=body_mag,
            sun_alt=sun_alt,
            elongation=elongation,
            moon_alt=moon_alt,
            moon_phase=moon_phase,
            moon_obj_angle=moon_body_sep,
            margin=_HELIACAL_VIS_MARGIN,
        )
        return visible, sun_alt, body_alt, elongation

    def _is_body_visible_no_moon(jd: float, margin: float = 0.0) -> bool:
        sun_alt, body_alt, _ = _get_altitudes(jd)
        elongation = _get_elongation(jd)
        if body_alt < 0 or sun_alt > 0:
            return False
        body_mag = _get_body_magnitude(jd)
        schaefer.update_season(jd)
        return schaefer.is_visible(
            body_alt=body_alt,
            body_mag=body_mag,
            sun_alt=sun_alt,
            elongation=elongation,
            moon_alt=-90.0,
            moon_phase=0.0,
            moon_obj_angle=180.0,
            margin=margin,
        )

    def _find_twilight_center(jd_day: float, morning: bool) -> float:
        best_ut = -1.0
        best_score = -999.0
        # Precompute the 25 hourly Sun altitudes once (Sun only: the body
        # position is not needed to locate the twilight centre).
        sun_hourly = [_get_sun_alt(jd_day + h / 24.0) for h in range(25)]
        for h in range(24):
            sun_alt = sun_hourly[h]
            if -22.0 < sun_alt < 2.0:
                sun_next = sun_hourly[h + 1]
                if morning and sun_next > sun_alt:
                    score = -abs(sun_alt + 8.0)
                    if best_ut < 0 or score > best_score:
                        best_score = score
                        best_ut = float(h)
                elif not morning and sun_next < sun_alt:
                    score = -abs(sun_alt + 8.0)
                    if best_ut < 0 or score > best_score:
                        best_score = score
                        best_ut = float(h)
        return best_ut

    def _check_twilight_visibility(jd_day: float, morning: bool):
        # Anchor the scan grid to 0h UT of jd_day's civil day. The search
        # loops pass jd_start + k, carrying jd_start's arbitrary time-of-day
        # into every scanned day, so the 15-min sampling grid shifted with
        # the requested start time and a few minutes' difference in jd_start
        # could flip a marginal detection (the first days of an apparition
        # are visible only for minutes), flipping in turn the apparition
        # classification of the whole search. A fixed per-day anchor makes
        # each day's verdict a function of the day alone.
        jd_day = math.floor(jd_day + 0.5) - 0.5
        center_ut = _find_twilight_center(jd_day, morning)
        if center_ut < 0:
            return False, 0.0
        # STRUCTURAL MIRROR of the Skyfield batched detector: same scan grid
        # (centre ±3 h at 15-min steps), same gates, same accept rule (FIRST
        # visible sample), same per-day magnitude epoch (0h UT of the civil
        # day) and the same alt/az-derived Sun-body elongation per scan
        # point. Any structural difference here yields backend-dependent
        # heliacal dates: a former ±45-min/3-min max-margin variant with a
        # centre-evaluated elongation systematically reported days the
        # Skyfield grid stepped over (~1-day drift on a third of cases and a
        # whole spurious Mercury apparition). Do not "improve" one twin
        # without changing the other identically.
        body_mag = _get_body_magnitude(jd_day)
        # Same morning/evening upper bound as the Skyfield twins: gate at
        # -5 deg in the morning to prevent false detections during civil
        # twilight where the sky-brightness model is unreliable, and at
        # -2 deg in the evening where setting bodies are only briefly
        # visible at shallow Sun depressions.
        sun_upper = -5.0 if morning else -2.0
        for dt_min in range(-180, 181, 15):
            ut_hour = center_ut + dt_min / 60.0
            jd_check = jd_day + ut_hour / 24.0
            sun_az, sun_alt, _sun_app = _leb_body_altaz(
                reader, jd_check, SUN, geopos, pressure, temperature
            )
            body_az, body_alt, _body_app = _leb_body_altaz(
                reader,
                jd_check,
                body,
                geopos,
                pressure,
                temperature,
                is_star=is_star,
                star_name=star_name,
            )
            if not (-18.0 < sun_alt < sun_upper and body_alt > 0.5):
                continue
            if body_alt < 0 or sun_alt > 0:
                continue
            # Elongation from topocentric alt/az, exactly like the Skyfield
            # batched detector (which derives it from the same coordinates
            # instead of a geocentric pheno call).
            _sa_r = math.radians(sun_alt)
            _ba_r = math.radians(body_alt)
            _daz_r = math.radians(sun_az - body_az)
            _cos_elong = math.sin(_sa_r) * math.sin(_ba_r) + math.cos(_sa_r) * math.cos(
                _ba_r
            ) * math.cos(_daz_r)
            elong = math.degrees(math.acos(max(-1.0, min(1.0, _cos_elong))))
            schaefer.update_season(jd_check)
            visible = schaefer.is_visible(
                body_alt=body_alt,
                body_mag=body_mag,
                sun_alt=sun_alt,
                elongation=elong,
                moon_alt=-90.0,
                moon_phase=0.0,
                moon_obj_angle=180.0,
                margin=_HELIACAL_VIS_MARGIN,
            )
            if visible and _event_side_ok(jd_check, morning):
                return True, jd_check
        return False, 0.0

    # --- Vectorized batch is NOT used in the LEB path ---
    # The LEB path is already ~14x faster per call than Skyfield,
    # so the simpler per-day loop is adequate.

    _HELIACAL_BATCH = 100

    def _batch_check_twilight_visibility(jd_days_list: list, morning: bool) -> list:
        """Per-day twilight visibility check (non-vectorized LEB version)."""
        return [_check_twilight_visibility(jd_day, morning) for jd_day in jd_days_list]

    def _refine_heliacal_time(jd_approx: float, is_morning: bool) -> float:
        jd_low = jd_approx - 0.1
        jd_high = jd_approx + 0.1
        for _ in range(30):
            jd_mid = (jd_low + jd_high) / 2
            visible, _sun_alt, _body_alt, _elong = _is_body_visible(jd_mid)
            if is_morning:
                if visible:
                    jd_high = jd_mid
                else:
                    jd_low = jd_mid
            else:
                if visible:
                    jd_low = jd_mid
                else:
                    jd_high = jd_mid
            if jd_high - jd_low < 1e-6:
                break
        return (jd_low + jd_high) / 2

    def _search_heliacal_rising(jd_start_inner: float) -> float:
        # A heliacal rising is the FIRST morning visibility after the body was
        # invisible because it was too close to the Sun (i.e. it emerged from a
        # conjunction gap). Requiring only N consecutive invisible days is not
        # enough: mid-apparition a bright planet can flicker invisible for a few
        # days (Moon glare / marginal twilight) and that would be mis-reported
        # as a rising. So we also require the invisibility streak to have passed
        # through low solar elongation (a real conjunction gap), which a
        # mid-apparition dip never does. There is one more legitimate case the
        # elongation test alone misses: when the conjunction happened *before*
        # the search window (the body starts the search already past conjunction
        # but still invisible, e.g. Mars at 17 deg and climbing). Then the streak
        # never dips below ELONG_GAP inside the window, yet the first visibility
        # is still a true heliacal rising. We accept it whenever the body has not
        # been visible at all since the start of the search (seen_visible False),
        # which a mid-apparition flicker -- visible before the dip -- never is.
        ELONG_GAP = 10.0  # deg: streak min elongation below this == near conjunction
        ELONG_START_GAP = 30.0
        max_days = 800
        lookback_jds = [jd_start_inner - i for i in range(1, 7)]
        lookback_vis = _batch_check_twilight_visibility(lookback_jds, morning=True)
        consecutive_invisible = 0
        min_elong = 999.0
        seen_visible = False
        for jd_lb, (vis, _) in zip(lookback_jds, lookback_vis):
            if not vis:
                consecutive_invisible += 1
                min_elong = min(min_elong, _get_elongation(jd_lb))
            else:
                seen_visible = True
                break
        if not seen_visible and min_elong > ELONG_GAP:
            extended_jds = [jd_start_inner - i for i in range(7, 31)]
            extended_vis = _batch_check_twilight_visibility(extended_jds, morning=True)
            for jd_lb, (vis, _) in zip(extended_jds, extended_vis):
                if vis:
                    seen_visible = True
                    break
                min_elong = min(min_elong, _get_elongation(jd_lb))
        for batch_start in range(0, max_days, _HELIACAL_BATCH):
            batch_end = min(batch_start + _HELIACAL_BATCH, max_days)
            jd_days = [jd_start_inner + d for d in range(batch_start, batch_end)]
            vis_results = _batch_check_twilight_visibility(jd_days, morning=True)
            for (vis, jd_vis), jd_day in zip(vis_results, jd_days):
                if not vis:
                    consecutive_invisible += 1
                    min_elong = min(min_elong, _get_elongation(jd_day))
                else:
                    current_elong = _get_elongation(jd_day)
                    pre_window_gap = (
                        not seen_visible
                        and not (flags & HELFLAG_HIGH_PRECISION)
                        and min_elong <= ELONG_START_GAP
                        and current_elong >= min_elong + 2.0
                    )
                    if consecutive_invisible >= 5 and (
                        is_star or min_elong <= ELONG_GAP or pre_window_gap
                    ):
                        return _refine_heliacal_time(jd_vis, is_morning=True)
                    consecutive_invisible = 0
                    min_elong = 999.0
                    seen_visible = True
        return 0.0

    def _search_heliacal_setting(jd_start_inner: float) -> float:
        max_days = 800
        last_visible_jd = 0.0
        consecutive_invisible = 0
        for batch_start in range(0, max_days, _HELIACAL_BATCH):
            batch_end = min(batch_start + _HELIACAL_BATCH, max_days)
            jd_days = [jd_start_inner + d for d in range(batch_start, batch_end)]
            vis_results = _batch_check_twilight_visibility(jd_days, morning=False)
            for vis, jd_vis in vis_results:
                if vis:
                    last_visible_jd = jd_vis
                    consecutive_invisible = 0
                else:
                    consecutive_invisible += 1
                    if consecutive_invisible >= 5 and last_visible_jd > 0:
                        return _refine_heliacal_time(last_visible_jd, is_morning=False)
        if last_visible_jd > 0:
            return _refine_heliacal_time(last_visible_jd, is_morning=False)
        return 0.0

    def _search_evening_first(jd_start_inner: float) -> float:
        ELONG_GAP = 10.0
        ELONG_START_GAP = 30.0
        max_days = 800

        # Fast-path only when the body is near conjunction (elongation below the
        # gap) AND invisible both this evening and on the 6 preceding evenings:
        # that is a genuine conjunction the search starts inside, so the next
        # evening visibility IS the evening-first event. If the body is visible
        # tonight (ongoing apparition) or was visible on a recent evening (a
        # marginal-conditions flicker inside the current apparition), this is
        # not a conjunction gap and we must fall through to the normal search
        # for the NEXT event after the coming superior conjunction.
        if not is_star and _get_elongation(jd_start_inner) <= ELONG_GAP:
            gate_jds = [jd_start_inner - i for i in range(0, 7)]
            gate_vis = _batch_check_twilight_visibility(gate_jds, morning=False)
            if not any(vis for vis, _ in gate_vis):
                for batch_start in range(1, max_days, _HELIACAL_BATCH):
                    batch_end = min(batch_start + _HELIACAL_BATCH, max_days)
                    jd_days = [
                        jd_start_inner + d for d in range(batch_start, batch_end)
                    ]
                    vis_results = _batch_check_twilight_visibility(
                        jd_days, morning=False
                    )
                    for vis, jd_vis in vis_results:
                        if vis:
                            return _refine_heliacal_time(jd_vis, is_morning=False)
                return 0.0

        lookback_jds = [jd_start_inner - i for i in range(1, 7)]
        lookback_vis = _batch_check_twilight_visibility(lookback_jds, morning=False)
        consecutive_invisible = 0
        min_elong = 999.0
        seen_visible = False
        for jd_lb, (vis, _) in zip(lookback_jds, lookback_vis):
            if not vis:
                consecutive_invisible += 1
                min_elong = min(min_elong, _get_elongation(jd_lb))
            else:
                seen_visible = True
                break
        for batch_start in range(0, max_days, _HELIACAL_BATCH):
            batch_end = min(batch_start + _HELIACAL_BATCH, max_days)
            jd_days = [jd_start_inner + d for d in range(batch_start, batch_end)]
            vis_results = _batch_check_twilight_visibility(jd_days, morning=False)
            for (vis, jd_vis), jd_day in zip(vis_results, jd_days):
                if not vis:
                    consecutive_invisible += 1
                    min_elong = min(min_elong, _get_elongation(jd_day))
                else:
                    current_elong = _get_elongation(jd_day)
                    pre_window_gap = (
                        not seen_visible
                        and not (flags & HELFLAG_HIGH_PRECISION)
                        and min_elong <= ELONG_START_GAP
                        and current_elong >= min_elong + 2.0
                    )
                    if consecutive_invisible >= 5 and (
                        is_star or min_elong <= ELONG_GAP or pre_window_gap
                    ):
                        return _refine_heliacal_time(jd_vis, is_morning=False)
                    consecutive_invisible = 0
                    min_elong = 999.0
                    seen_visible = True
        return 0.0

    def _search_morning_last(jd_start_inner: float) -> float:
        ELONG_GAP = 10.0
        max_days = 800
        last_visible_jd = 0.0
        found_visible = False
        consecutive_invisible = 0
        streak_min_elong = 999.0
        for batch_start in range(0, max_days, _HELIACAL_BATCH):
            batch_end = min(batch_start + _HELIACAL_BATCH, max_days)
            jd_days = [jd_start_inner + d for d in range(batch_start, batch_end)]
            vis_results = _batch_check_twilight_visibility(jd_days, morning=True)
            for (vis, jd_vis), jd_day in zip(vis_results, jd_days):
                if vis:
                    last_visible_jd = jd_vis
                    found_visible = True
                    consecutive_invisible = 0
                    streak_min_elong = 999.0
                elif found_visible:
                    consecutive_invisible += 1
                    streak_min_elong = min(streak_min_elong, _get_elongation(jd_day))
                    if (
                        consecutive_invisible >= 3
                        and streak_min_elong <= ELONG_GAP
                        and last_visible_jd > 0
                    ):
                        return _refine_heliacal_time(last_visible_jd, is_morning=True)
        return last_visible_jd if last_visible_jd > 0 else 0.0

    # --- main dispatch ---
    if event_type == HELIACAL_RISING:
        _search = _search_heliacal_rising
    elif event_type == HELIACAL_SETTING:
        _search = _search_heliacal_setting
    elif event_type == EVENING_FIRST:
        _search = _search_evening_first
    elif event_type == MORNING_LAST:
        _search = _search_morning_last
    else:
        _search = None

    jd_event = 0.0
    if _search is not None:
        jd_event = _search(jd_start)
        # Contract gate: the returned event must be >= jd_start. Day 0 of
        # the search is the civil day containing jd_start and its twilight
        # scan covers that whole day, so when jd_start falls later in the
        # day than the visibility optimum, the streak gate can classify the
        # already-past optimum as the event and the refinement converges
        # before jd_start. Restarting one day later moves that day into the
        # lookback, where its visibility marks the apparition as ongoing --
        # the same branch a search started the following day takes -- and
        # the search skips to the next apparition. The retry cannot violate
        # the contract again: every instant it scans lies on civil days
        # strictly after jd_start's.
        if 0.0 < jd_event < jd_start - 1e-6:
            jd_event = _search(jd_start + 1.0)

    if jd_event > 0:
        return jd_event, event_type
    else:
        return 0.0, -1


def _heliacal_pheno_ut_leb(
    reader,
    jd: float,
    lat: float,
    lon: float,
    altitude: float,
    pressure: float,
    temperature: float,
    humidity: float,
    body: int,
    event_type: int,
    flags: int,
    met_range: float = 0.0,
    observer: Sequence[float] = _NAKED_OBSERVER,
) -> tuple:
    """LEB-backed twin of ``_heliacal_pheno_ut_pythonic``.

    Same inputs (preceded by the LEB ``reader``) and the same 50-slot result
    layout, documented there; positions come from the LEB reader instead of
    the Skyfield backend.
    """
    from .constants import (
        HELIACAL_RISING,
        HELIACAL_SETTING,
        EVENING_FIRST,
        MORNING_LAST,
    )
    from .planets import _PLANET_MAP, pheno_ut
    from .utils import angular_separation

    if event_type not in (
        HELIACAL_RISING,
        HELIACAL_SETTING,
        EVENING_FIRST,
        MORNING_LAST,
    ):
        raise Error(
            f"invalid event type {event_type}; use HELIACAL_RISING, "
            "HELIACAL_SETTING, EVENING_FIRST or MORNING_LAST"
        )

    is_star = is_fixed_star(body)
    if not is_star and body not in _PLANET_MAP:
        raise Error(f"body id {body} is neither a planet nor a fixed star")

    star_name: str | None = None
    star_magnitude = 0.0
    if is_star:
        star_name = _star_name_from_id(body)
        star_magnitude = _get_star_magnitude(body)

    geopos = (lon, lat, altitude)

    # --- positions via LEB ---
    # The geometric slots (0, 2, 3, 4, 5 and the arcs derived from them) use
    # geometric (TRUEPOS-like) topocentric places, i.e. the astrometric
    # direction without aberration or light-time, exactly as the Skyfield twin
    # does. Only slot 1 starts from the apparent altitude, to which the
    # refraction model below is applied; the visibility detector (which calls
    # the shared _leb_body_altaz) is unaffected.
    from .utils import azalt as _azalt, ECL2HOR as _ECL2HOR
    from .constants import FLG_TRUEPOS as _FLG_TRUEPOS
    from .fast_calc import _topo_ecliptic as _topo_ecliptic
    from .time_utils import deltat as _deltat

    jd_tt = jd + _deltat(jd)

    def _geometric_topo_altaz(body_id: int, star_flag: bool) -> tuple[float, float]:
        """Geometric topocentric (north-zero azimuth, unrefracted altitude)."""
        if star_flag:
            from .fixed_stars import fixstar_ut as _fixstar_ut

            assert star_name is not None
            star_pos, _, _ = _fixstar_ut(star_name, jd, _FLG_TRUEPOS)
            ecl_lon, ecl_lat, ecl_dist = star_pos[0], star_pos[1], star_pos[2]
        else:
            topo_ecl = _topo_ecliptic(
                reader, jd_tt, jd, body_id, geopos, iflag=_FLG_TRUEPOS
            )
            ecl_lon, ecl_lat, ecl_dist = topo_ecl[0], topo_ecl[1], topo_ecl[2]
        az_south_zero, alt_unrefracted, _ = _azalt(
            jd, _ECL2HOR, geopos, 0.0, 0.0, (ecl_lon, ecl_lat, ecl_dist)
        )
        return (az_south_zero + 180.0) % 360.0, alt_unrefracted

    # Sun: geometric topocentric altitude/azimuth (slots 4 and 5).
    sun_az_deg, sun_alt_deg = _geometric_topo_altaz(SUN, False)

    # Object apparent topocentric altitude, used only as the base of the
    # refracted altitude in slot 1.
    _body_app_az_deg, body_app_alt_deg, _ = _leb_body_altaz(
        reader,
        jd,
        body,
        geopos,
        pressure,
        temperature,
        is_star=is_star,
        star_name=star_name,
    )

    # Object: geometric topocentric altitude/azimuth (slots 0 and 3).
    body_az_deg, body_alt_deg = _geometric_topo_altaz(body, is_star)

    # Atmospheric refraction from the APPARENT true altitude, in arcminutes:
    # R = 1.02 / tan(h + 10.3/(h + 5.11)).  Sæmundsson, Þ. (1986),
    # "Astronomical Refraction", Sky & Telescope 72, 70; as given in
    # Meeus (1998), Astronomical Algorithms 2nd ed., ch. 16, eq. 16.4.
    # Below h = -1° the 0.5° fallback approximates horizon refraction; slot 1
    # is this refracted apparent altitude.
    if body_app_alt_deg > -1:
        refraction = 1.02 / math.tan(
            math.radians(body_app_alt_deg + 10.3 / (body_app_alt_deg + 5.11))
        )
        refraction /= 60.0
    else:
        refraction = 0.5
    app_alt_deg = body_app_alt_deg + refraction

    # Get ecliptic positions for parallax and elongation calculations
    body_ecliptic = _leb_ecliptic_pos(
        reader,
        jd,
        body,
        geopos,
        is_star=is_star,
        star_name=star_name,
    )
    sun_ecliptic = _leb_ecliptic_pos(reader, jd, SUN, geopos)

    # Geocentric altitude from equatorial coords (matching Skyfield path).
    # GAST is referenced to the true equinox of date, so the RA/Dec must be of
    # date (FLG_EQUATORIAL alone, NOT FLG_J2000) to keep the hour-angle frame
    # consistent. Mixing J2000 RA with of-date GAST injected the full
    # J2000→date precession+nutation into the hour angle. FLG_TRUEPOS selects
    # the astrometric (aberration- and light-time-free) direction used for the
    # geometric slots.
    from .constants import FLG_EQUATORIAL

    if is_star:
        from .fixed_stars import fixstar_ut

        assert star_name is not None
        equatorial_pos, _, _ = fixstar_ut(
            star_name, jd, FLG_EQUATORIAL | _FLG_TRUEPOS | FLG_SPEED
        )
        body_ra_deg, body_dec_deg = equatorial_pos[0], equatorial_pos[1]
    else:
        from .planets import calc_ut as _calc_ut

        equatorial_pos, _ = _calc_ut(
            jd, body, FLG_EQUATORIAL | _FLG_TRUEPOS | FLG_SPEED
        )
        body_ra_deg, body_dec_deg = equatorial_pos[0], equatorial_pos[1]

    # Use GAST to match the Skyfield path's of-date RA + t.gast
    from .state import get_timescale as _get_timescale

    timescale = _get_timescale()
    time_ut1 = timescale.ut1_jd(jd)
    gast_hours = time_ut1.gast
    lst_hours = gast_hours + geopos[0] / 15.0
    hour_angle_deg = (lst_hours * 15.0 - body_ra_deg) % 360.0
    hour_angle_rad = math.radians(hour_angle_deg)
    dec_rad = math.radians(body_dec_deg)
    lat_rad = math.radians(geopos[1])
    sin_geo_alt = math.sin(lat_rad) * math.sin(dec_rad) + math.cos(lat_rad) * math.cos(
        dec_rad
    ) * math.cos(hour_angle_rad)
    geo_alt_deg = math.degrees(math.asin(max(-1.0, min(1.0, sin_geo_alt))))

    # Arcus visionis
    topo_arcus_visionis = body_alt_deg - sun_alt_deg
    geo_arcus_visionis = geo_alt_deg - sun_alt_deg

    # Public signed-azimuth convention: Sun minus object.
    azimuth_diff_deg = sun_az_deg - body_az_deg
    while azimuth_diff_deg > 180:
        azimuth_diff_deg -= 360
    while azimuth_diff_deg < -180:
        azimuth_diff_deg += 360

    # Elongation and magnitude
    if is_star:
        elongation = angular_separation(
            sun_ecliptic[0],
            sun_ecliptic[1],
            body_ecliptic[0],
            body_ecliptic[1],
        )
        magnitude = star_magnitude
        phase_angle = 0.0
    else:
        try:
            body_pheno = pheno_ut(jd, body, _heliacal_eph_flags(flags))
            elongation = body_pheno[2]
            magnitude = body_pheno[4]
            phase_angle = body_pheno[0]
        except (ValueError, TypeError, ArithmeticError):
            elongation = angular_separation(
                sun_ecliptic[0],
                sun_ecliptic[1],
                body_ecliptic[0],
                body_ecliptic[1],
            )
            magnitude = 0.0
            phase_angle = 0.0

    # Arc of light (slot 9) between the body and the Sun via the exact
    # right-spherical-triangle relation cos ARCL = cos ARCV * cos DAZ.
    arc_of_light_deg = _arc_of_light(geo_arcus_visionis, azimuth_diff_deg)

    (
        observer_age,
        snellen_ratio,
        binocular_flag,
        telescope_magnification,
        aperture_mm,
        optics_transmission,
    ) = _parse_observer_optics(observer, flags)
    visibility_model = create_schaefer_model(
        pressure=pressure,
        temperature=temperature,
        humidity=humidity * 100.0 if humidity <= 1.0 else humidity,
        altitude=altitude,
        met_range=met_range,
        observer_age=observer_age,
        snellen=snellen_ratio,
        binocular=binocular_flag,
        telescope_mag=telescope_magnification,
        aperture=aperture_mm,
        transmission=optics_transmission,
        latitude=lat,
        jd=jd,
    )
    extinction_coeff = visibility_model.k_total
    required_arcus_visionis = visibility_model.arcus_visionis_required(magnitude)

    # Parallax
    if is_star:
        altitude_parallax_deg = 0.0
    else:
        # Altitude parallax is the geocentric-minus-topocentric altitude, not
        # the horizontal parallax angle.
        altitude_parallax_deg = geo_alt_deg - body_alt_deg

    # Rise/set of the disc centers and the visibility window via the
    # shared event-time helper. Missing instants use the public sentinel.
    catalog_star_name = ""
    if is_star:
        from .fixed_stars import STAR_CATALOG as _catalog

        for catalog_entry in _catalog:
            if catalog_entry.id == body:
                catalog_star_name = catalog_entry.name
                break
    (
        object_rise_set_jd,
        sun_rise_set_jd,
        lag_days,
        visibility_start_jd,
        visibility_best_jd,
        visibility_end_jd,
        visibility_duration_days,
        yallop_best_jd,
    ) = _pheno_rise_window(
        jd,
        body,
        catalog_star_name,
        event_type,
        (lon, lat, altitude),
        (
            pressure,
            temperature,
            humidity * 100.0 if humidity <= 1.0 else humidity,
            met_range,
        ),
        # Caller's observer: the window slots depend on it (see the
        # Skyfield-path twin).
        tuple(observer[:6])
        if observer is not None
        else (36.0, 1.0, 0.0, 0.0, 0.0, 0.0),
        flags,
    )

    if is_star:
        # Fixed stars are unresolved point sources; the public illumination
        # field uses the full-illumination convention for them.
        illuminated_pct = 100.0
    elif phase_angle > 0:
        illuminated_pct = (1.0 + math.cos(math.radians(phase_angle))) / 2.0 * 100.0
    else:
        illuminated_pct = 0.0

    crescent_width_deg = 0.0
    crescent_length_deg = 0.0
    if body == MOON:
        try:
            moon_pheno = pheno_ut(jd, MOON, _heliacal_eph_flags(flags))
            illuminated_pct = moon_pheno[1] * 100.0  # [1] is illuminated fraction 0-1
            moon_diameter_deg = (
                moon_pheno[3] if len(moon_pheno) > 3 and moon_pheno[3] > 0 else 0.5
            )
            # Crescent width W = SD * (1 - cos ARCL) from Bruin, F. (1977),
            # "The first visibility of the lunar crescent", Vistas in Astronomy
            # 21, 331-358 (adopted by Yallop's q-test below), where
            # ARCL is the arc of light (Sun-Moon elongation, moon_pheno[2]) and
            # SD is the Moon's apparent semidiameter in arcmin (moon_pheno[3]
            # is the apparent diameter in degrees). Near new moon ARCL -> 0 so
            # W -> 0; near full moon ARCL -> 180 deg so W -> the disc diameter.
            moon_arc_of_light_deg = moon_pheno[2]
            semidiameter_arcmin = moon_diameter_deg / 2.0 * 60.0
            crescent_width_arcmin = semidiameter_arcmin * (
                1.0 - math.cos(math.radians(moon_arc_of_light_deg))
            )
            crescent_width_deg = crescent_width_arcmin / 60.0  # slot 16, degrees
            crescent_length_deg = math.pi * moon_diameter_deg / 2
        except (ValueError, TypeError, ArithmeticError):
            pass

    if body == MOON:
        # Yallop q-test (Yallop, B.D. (1997), NAO Technical Note No. 69):
        # the cubic criterion in the crescent width (arcmin) and the
        # q = (ARCV - criterion) / 10 statistic.
        width_arcmin = crescent_width_deg * 60.0
        yallop_threshold = (
            11.8371
            - 6.3226 * width_arcmin
            + 0.7319 * width_arcmin**2
            - 0.1018 * width_arcmin**3
        )
        yallop_q = (geo_arcus_visionis - yallop_threshold) / 10.0
        # Slot 18 carries the Yallop visibility class code (1..6), not the raw
        # arcus-visionis criterion polynomial (which stays internal above).
        yallop_class_code = _yallop_visibility_code(yallop_q)
    else:
        yallop_q = 0.0
        yallop_class_code = 0.0

    values = [0.0] * 50
    values[0] = body_alt_deg
    values[1] = app_alt_deg
    values[2] = geo_alt_deg
    values[3] = body_az_deg
    values[4] = sun_alt_deg
    values[5] = sun_az_deg
    values[6] = topo_arcus_visionis
    values[7] = geo_arcus_visionis
    values[8] = azimuth_diff_deg
    values[9] = arc_of_light_deg
    values[10] = extinction_coeff
    values[11] = required_arcus_visionis
    values[12] = visibility_start_jd
    values[13] = visibility_best_jd
    values[14] = visibility_end_jd
    values[15] = yallop_best_jd
    values[16] = crescent_width_deg
    values[17] = yallop_q
    values[18] = yallop_class_code
    values[19] = altitude_parallax_deg
    values[20] = magnitude
    values[21] = object_rise_set_jd
    values[22] = sun_rise_set_jd
    values[23] = lag_days
    values[24] = visibility_duration_days
    values[25] = crescent_length_deg
    values[26] = elongation
    values[27] = illuminated_pct

    return tuple(values), flags


def _vislim_scotopic_flag(
    schaefer,
    sun_alt: float,
    moon_alt: float,
    moon_phase: float,
    obj_alt: float,
    sun_obj_angle: float,
    moon_obj_angle: float,
    flags: int,
) -> int:
    """Public scotopic/photopic return flag for ``vis_limit_mag``.

    Schaefer's VISLIMIT contrast threshold switches from the photopic to
    the scotopic (dark, night-vision) branch when the total sky
    background at the object falls below 1500 nanolamberts. This mirrors
    that switch directly on the model's total sky brightness (nL).
    HELFLAG_VISLIM_PHOTOPIC / HELFLAG_VISLIM_SCOTOPIC override it.
    """
    from .constants import HELFLAG_PHOTOPIC, HELFLAG_SCOTOPIC

    bsk_nl = schaefer.sky_brightness_bl(
        sun_alt, moon_alt, moon_phase, obj_alt, sun_obj_angle, moon_obj_angle
    )
    is_scotopic = bsk_nl < 1500.0
    if flags & HELFLAG_VISLIM_PHOTOPIC:
        is_scotopic = False
    if flags & HELFLAG_VISLIM_SCOTOPIC:
        is_scotopic = True
    # The public result contract uses a binary photopic/scotopic encoding; the
    # below-horizon case is handled by the caller.
    return HELFLAG_SCOTOPIC if is_scotopic else HELFLAG_PHOTOPIC


def _vis_limit_mag_leb(
    reader,
    tjdut: float,
    geopos: tuple,
    atmo: tuple,
    observer_tup: tuple,
    objname: str,
    flags: int,
) -> tuple:
    """LEB-backed implementation of vis_limit_mag()."""
    from .planets import _PLANET_MAP, pheno_ut
    from .fixed_stars import fixstar2_ut, fixstar2_mag
    from .utils import azalt, angular_separation, ECL2HOR
    from .constants import (
        SUN,
        MOON,
        HELFLAG_VISLIM_DARK,
        HELFLAG_VISLIM_NOMOON,
        HELFLAG_BELOW_HORIZON,
    )

    if not objname:
        raise Error("the object name is empty")

    lon = geopos[0] if len(geopos) > 0 else 0.0
    lat = geopos[1] if len(geopos) > 1 else 0.0
    alt_m = geopos[2] if len(geopos) > 2 else 0.0

    pressure = atmo[0] if len(atmo) > 0 else 1013.25
    temperature = atmo[1] if len(atmo) > 1 else 15.0
    humidity_pct = atmo[2] if len(atmo) > 2 else 50.0
    met_range = atmo[3] if len(atmo) > 3 else 0.0

    (
        observer_age,
        snellen_ratio,
        obs_binocular,
        obs_telescope_mag,
        obs_aperture,
        obs_transmission,
    ) = _parse_observer_optics(observer_tup, flags)

    geopos_ll = (lon, lat, alt_m)

    # Sun position
    sun_az, sun_alt, _ = _leb_body_altaz(
        reader,
        tjdut,
        SUN,
        geopos_ll,
        pressure,
        temperature,
    )

    # Moon position
    moon_az, moon_alt, _ = _leb_body_altaz(
        reader,
        tjdut,
        MOON,
        geopos_ll,
        pressure,
        temperature,
    )

    # Determine object
    body_id = None
    is_star_flag = False

    try:
        body_id = int(objname)
    except ValueError:
        name_upper = objname.upper().strip()
        planet_names = {
            "SUN": SUN,
            "MOON": MOON,
            "MERCURY": 2,
            "VENUS": 3,
            "MARS": 4,
            "JUPITER": 5,
            "SATURN": 6,
            "URANUS": 7,
            "NEPTUNE": 8,
            "PLUTO": 9,
        }
        if name_upper in planet_names:
            body_id = planet_names[name_upper]
        else:
            is_star_flag = True

    # The Sun has no meaningful limiting magnitude in this model.
    if not is_star_flag and body_id == SUN:
        raise Error("vis_limit_mag() is not defined for the Sun")

    obj_alt = 0.0
    obj_az = 0.0
    obj_mag = 0.0

    if is_star_flag:
        try:
            star_result, _star_name_out, _retflag = fixstar2_ut(
                objname, tjdut, _heliacal_eph_flags(flags)
            )
            star_lon = star_result[0]
            star_lat = star_result[1]

            try:
                star_mag_val, _star_name_mag = fixstar2_mag(objname)
                obj_mag = star_mag_val
            except (ValueError, TypeError, ArithmeticError):
                obj_mag = 2.0

            hor_result = azalt(
                tjdut,
                ECL2HOR,
                geopos_ll,
                pressure,
                temperature,
                (star_lon, star_lat, 1.0),
            )
            obj_az = (hor_result[0] + 180.0) % 360.0
            obj_alt = hor_result[1]
        except ValueError:
            raise
        except (KeyError, IndexError, OSError) as e:
            raise Error(
                f"no fixed star matches the name {objname.lower()!r}: {e}"
            ) from e
    else:
        if body_id is None:
            raise Error(f"unknown object name {objname!r}")

        if body_id in _PLANET_MAP:
            az, alt_true, _ = _leb_body_altaz(
                reader,
                tjdut,
                body_id,
                geopos_ll,
                pressure,
                temperature,
            )
            obj_alt = alt_true
            obj_az = az

            try:
                pheno_result = pheno_ut(tjdut, body_id, _heliacal_eph_flags(flags))
                obj_mag = pheno_result[4]
            except (ValueError, TypeError, ArithmeticError):
                obj_mag = 0.0
        else:
            raise Error(f"body id {body_id} is neither a planet nor a fixed star")

    if obj_alt < 0:
        # Public result contract: the below-horizon flag uses a sentinel in
        # the limiting-magnitude slot and zeros in the remaining slots.
        dret = (-100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        return float(HELFLAG_BELOW_HORIZON), dret

    use_dark_sky = bool(flags & HELFLAG_VISLIM_DARK)
    exclude_moon = bool(flags & HELFLAG_VISLIM_NOMOON)

    if use_dark_sky:
        sun_alt = -90.0
    if exclude_moon:
        moon_alt = -90.0

    if body_id == MOON:
        # Compatibility contract: when the observed object IS the Moon,
        # its own light is not sky glow hindering the observation, so the
        # moonlight term is removed exactly as with VISLIM_NOMOON (the sky
        # stays dark and the limiting magnitude reflects the object alone).
        moon_alt = -90.0

    schaefer = create_schaefer_model(
        pressure=pressure,
        temperature=temperature,
        humidity=humidity_pct,
        met_range=met_range,
        altitude=alt_m,
        observer_age=observer_age,
        snellen=snellen_ratio,
        binocular=obs_binocular,
        telescope_mag=obs_telescope_mag,
        aperture=obs_aperture,
        transmission=obs_transmission,
        latitude=lat,
        jd=tjdut,
        pressure_scales_extinction=True,
    )

    # Moon phase and angular separations via LEB ecliptic positions
    moon_phase = 0.0
    moon_obj_angle = 180.0
    sun_obj_angle = 180.0
    try:
        moon_pheno = pheno_ut(tjdut, MOON, _heliacal_eph_flags(flags))
        # pheno_ut[1] is the illuminated fraction (0 = new, 1 = full), exactly
        # what SchaeferModel.sky_brightness_moon expects. Deriving it from the
        # phase angle (pheno_ut[0], 0 deg = full) as (1 - cos)/2 inverts the
        # convention — full moon would be passed as new and vice versa.
        moon_phase = moon_pheno[1]

        if is_star_flag:
            from .fixed_stars import fixstar_ut as _sfut_vlm

            _star_ecl, _, _ = _sfut_vlm(objname, tjdut, FLG_SPEED)
            _star_pos = (_star_ecl[0], _star_ecl[1])
        else:
            assert body_id is not None
            _body_ecl = _leb_ecliptic_pos(reader, tjdut, body_id, geopos_ll)
            _star_pos = (_body_ecl[0], _body_ecl[1])
        moon_ecl = _leb_ecliptic_pos(reader, tjdut, MOON, geopos_ll)
        sun_ecl = _leb_ecliptic_pos(reader, tjdut, SUN, geopos_ll)
        moon_obj_angle = angular_separation(
            _star_pos[0],
            _star_pos[1],
            moon_ecl[0],
            moon_ecl[1],
        )
        sun_obj_angle = angular_separation(
            _star_pos[0],
            _star_pos[1],
            sun_ecl[0],
            sun_ecl[1],
        )
    except (ValueError, TypeError, AttributeError):
        pass

    limiting_mag = schaefer.limiting_magnitude(
        sun_alt=sun_alt,
        moon_alt=moon_alt if not exclude_moon else -90.0,
        moon_phase=moon_phase,
        obj_alt=obj_alt,
        sun_obj_angle=sun_obj_angle,
        moon_obj_angle=moon_obj_angle,
    )

    vision_type = _vislim_scotopic_flag(
        schaefer,
        sun_alt if not use_dark_sky else -90.0,
        moon_alt if not exclude_moon else -90.0,
        moon_phase,
        obj_alt,
        sun_obj_angle,
        moon_obj_angle,
        flags,
    )

    dret = (
        limiting_mag,
        obj_alt,
        obj_az,
        sun_alt,
        sun_az,
        moon_alt,
        moon_az,
        obj_mag,
        0.0,
        0.0,
    )
    return float(vision_type), dret


def _heliacal_ut_pythonic(
    jd_start: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    pressure: float = 1013.25,
    temperature: float = 15.0,
    humidity: float = 0.5,
    body: int = SUN,
    event_type: int = 1,
    flags: int = FLG_SWIEPH,
    met_range: float = 0.0,
    observer: Sequence[float] = _NAKED_OBSERVER,
) -> Tuple[float, int]:
    """
    Calculate heliacal rising or setting time for a celestial body.

    Heliacal events are the first/last visibility of a celestial body
    at dawn or dusk. These were fundamental for ancient calendars:
    - Heliacal rising: First morning visibility after a period of invisibility
    - Heliacal setting: Last evening visibility before becoming invisible

    Args:
        jd_start: Julian Day (UT) to start search from
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        pressure: Atmospheric pressure in mbar/hPa for refraction (default 1013.25)
        temperature: Temperature in Celsius for refraction (default 15)
        humidity: Relative humidity 0.0-1.0 for atmospheric extinction (default 0.5)
        body: Planet/body ID (MERCURY, VENUS, MARS, JUPITER, SATURN,
              or fixed stars). Note: SUN and MOON are not valid for heliacal events.
        event_type: Type of heliacal event:
            - HELIACAL_RISING (1): Morning first visibility (heliacal rising)
            - HELIACAL_SETTING (2): Evening last visibility (heliacal setting)
            - EVENING_FIRST (3): First evening visibility (after superior conjunction)
              Note: Only valid for inner planets (Mercury, Venus)
            - MORNING_LAST (4): Last morning visibility (before superior conjunction)
              Note: Only valid for inner planets (Mercury, Venus)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - jd_event: Julian Day (UT) of the heliacal event, or 0.0 if not found
            - retflag: Return flag (event_type on success, negative on error)

    Raises:
        ValueError: If invalid body ID, event_type, or if EVENING_FIRST/MORNING_LAST
            is used with an outer planet (Mars, Jupiter, Saturn, etc.)

    Algorithm:
        The algorithm searches for the moment when:
        1. The body is at a specific altitude above the horizon (arcus visionis)
        2. The Sun is at twilight position (typically -6 to -12 below horizon)
        3. The body's apparent magnitude is brighter than the sky's limiting magnitude

        For heliacal rising (morning first):
        - Search forward for when the body first becomes visible at dawn
        - Body must be above horizon while Sun is still below
        - Sky must be dark enough for the body to be seen

        For heliacal setting (evening last):
        - Search forward for when the body is last visible at dusk
        - Body must be above horizon while Sun is setting
        - Sky brightness must not overwhelm the body's light

    Historical Note:
        Heliacal risings were crucial for ancient calendars. The heliacal
        rising of Sirius marked the Egyptian new year and predicted the
        Nile flood. Babylonians used heliacal events to track planetary
        positions without modern instruments.

    Example:
        >>> from libephemeris import julday, _heliacal_ut_pythonic, VENUS, HELIACAL_RISING
        >>> jd = julday(2024, 1, 1, 0)
        >>> # Find next heliacal rising of Venus from Rome
        >>> jd_event, flag = _heliacal_ut_pythonic(jd, 41.9, 12.5, body=VENUS,
        ...                              event_type=HELIACAL_RISING)
        >>> print(f"Heliacal rising at JD {jd_event:.5f}")

    References:
        - Schaefer (1990), PASP 102, 212-229, DOI 10.1086/132629,
          limiting-magnitude/observer model.
        - Schaefer (1993), Vistas in Astronomy 36, 311-361; Schaefer (1998),
          Sky & Telescope 95(5), 57-60 (the VISLIMIT V-band model and its
          per-component relative airmasses actually used here).
        - Rozenberg (1966), "Twilight: A Study in Atmospheric Optics"
          (low-altitude sky airmass form used by the VISLIMIT path).

        The historical note above is context only; it does not provide the
        numerical event algorithm or its thresholds.
    """
    # --- LEB fast path ---
    from .state import get_leb_reader as _get_leb_reader

    _leb_rdr = _get_leb_reader()
    if _leb_rdr is not None:
        try:
            return _heliacal_ut_leb(
                _leb_rdr,
                jd_start,
                lat,
                lon,
                altitude,
                pressure,
                temperature,
                humidity,
                body,
                event_type,
                flags,
                met_range,
                observer,
            )
        except LEBCorruptionError:
            raise
        except (KeyError, ValueError) as _leb_err:
            from .eclipse import _raise_if_sealed_leb_miss

            # Sealed leb mode raises the documented typed error and never
            # continues past the LEB path.
            _raise_if_sealed_leb_miss(_leb_err)
            # Auto mode may continue through its normal backend chain.
            from .leb_reader import log_leb_fallback

            log_leb_fallback("heliacal", _leb_err)
    # --- END LEB fast path ---

    from .constants import (
        HELIACAL_RISING,
        HELIACAL_SETTING,
        EVENING_FIRST,
        MORNING_LAST,
    )
    from .planets import _PLANET_MAP, pheno_ut
    from .state import get_planets, get_timescale
    from skyfield.api import wgs84

    # Validate event type
    if event_type not in (
        HELIACAL_RISING,
        HELIACAL_SETTING,
        EVENING_FIRST,
        MORNING_LAST,
    ):
        raise Error(
            f"invalid event type {event_type}; use HELIACAL_RISING, "
            "HELIACAL_SETTING, EVENING_FIRST or MORNING_LAST"
        )

    # EVENING_FIRST and MORNING_LAST are only valid for inner planets
    # (Mercury and Venus) because they relate to superior conjunction visibility.
    # Outer planets only have one type of conjunction and these events don't apply.
    # Fixed stars only have heliacal rising and setting.
    if event_type in (EVENING_FIRST, MORNING_LAST):
        if is_fixed_star(body):
            raise Error(
                "EVENING_FIRST and MORNING_LAST do not apply to fixed stars; "
                "use HELIACAL_RISING or HELIACAL_SETTING"
            )
        if not is_inner_planet(body) and body != MOON:
            raise Error(
                "EVENING_FIRST and MORNING_LAST apply only to the inner planets "
                "(Mercury, Venus) and the Moon; for the outer planets use "
                "HELIACAL_RISING or HELIACAL_SETTING"
            )

    # Sun and Moon are not valid for heliacal events
    if body == SUN:
        raise Error("the Sun is not a valid object for heliacal calculations")
    if body == MOON and event_type in (HELIACAL_RISING, HELIACAL_SETTING):
        raise Error(
            "heliacal rising and setting (event types 1 and 2) are not defined "
            "for the Moon"
        )

    # Check if this is a fixed star
    is_star = is_fixed_star(body)

    # Validate body - must be either a known planet or a fixed star
    if not is_star and body not in _PLANET_MAP:
        raise Error(f"body id {body} is neither a planet nor a fixed star")

    # Get ephemeris and timescale
    eph = get_planets()
    ts = get_timescale()

    sun = eph["sun"]
    earth = eph["earth"]

    # Get target - either planet or star
    target = None
    star_object = None
    star_magnitude = 0.0  # Default, will be set for stars

    if is_star:
        # For fixed stars, create a Skyfield Star object
        from skyfield.api import Star
        from .fixed_stars import STAR_CATALOG

        star_magnitude = _get_star_magnitude(body)

        # Find star data from catalog
        star_data = None
        for entry in STAR_CATALOG:
            if entry.id == body:
                star_data = entry.data
                break

        if star_data is None:
            raise Error(f"unknown fixed star id {body}")

        # Create Skyfield Star object for position calculations
        star_object = Star(
            ra_hours=star_data.ra_j2000 / 15.0,
            dec_degrees=star_data.dec_j2000,
            ra_mas_per_year=star_data.pm_ra * 1000.0,
            dec_mas_per_year=star_data.pm_dec * 1000.0,
        )
    else:
        # For planets, get the target from ephemeris
        target_name = _PLANET_MAP[body]
        from .planets import _PLANET_FALLBACK

        try:
            target = eph[target_name]
        except KeyError:
            if target_name in _PLANET_FALLBACK:
                target = eph[_PLANET_FALLBACK[target_name]]
            else:
                raise

    # Parse the observer tuple BEFORE the name is rebound to the Skyfield
    # GeographicPosition below (the tuple parameter and the Skyfield
    # observer share the name in this function).
    (
        _obs_age,
        _obs_snellen,
        _obs_bino,
        _obs_mag,
        _obs_aper,
        _obs_trans,
    ) = _parse_observer_optics(observer, flags)

    # Create observer location
    observer = wgs84.latlon(lat, lon, altitude)

    # Note: Arcus visionis (minimum altitude for visibility) and
    # sun altitude thresholds are used in the visibility check functions.
    # Typical values: arcus visionis ~7 for bright objects, sun altitude ~-8.

    def _get_altitudes(jd: float) -> Tuple[float, float, float]:
        """Get Sun altitude, body altitude, and body azimuth at given JD."""
        t = ts.ut1_jd(jd)
        observer_at = earth + observer

        # Sun position
        sun_app = observer_at.at(t).observe(sun).apparent()
        sun_alt, _, _ = sun_app.altaz()

        # Body position (handle both planets and stars)
        if is_star and star_object is not None:
            body_app = observer_at.at(t).observe(star_object).apparent()
        else:
            body_app = observer_at.at(t).observe(target).apparent()
        body_alt, body_az, _ = body_app.altaz()

        return sun_alt.degrees, body_alt.degrees, body_az.degrees

    def _get_elongation(jd: float) -> float:
        """Get the elongation of body from Sun in degrees."""
        if not is_star:
            try:
                pheno = pheno_ut(jd, body, _heliacal_eph_flags(flags))
                return pheno[2]  # Elongation
            except (ValueError, TypeError, ArithmeticError):
                pass

        # For stars or fallback: calculate elongation manually
        t = ts.ut1_jd(jd)
        sun_app = earth.at(t).observe(sun).apparent()
        if is_star and star_object is not None:
            body_app = earth.at(t).observe(star_object).apparent()
        else:
            body_app = earth.at(t).observe(target).apparent()
        return body_app.separation_from(sun_app).degrees

    def _sun_body_lon_diff(jd: float) -> float:
        """Signed geocentric ecliptic longitude of the body minus the Sun.

        Normalised to (-180, 180]. Negative means the body is WEST of the Sun
        (rises before it -- a morning apparition); positive means EAST of the
        Sun (sets after it -- an evening apparition). Used to reject a twilight
        detection whose Sun-relative side contradicts the requested event: at
        high latitudes in summer the morning and evening twilights are only a
        few hours apart, so a scan window centred near local midnight can spill
        into the opposite twilight and latch the wrong apparition.
        """
        t = ts.ut1_jd(jd)
        _, sun_lon, _ = earth.at(t).observe(sun).apparent().ecliptic_latlon()
        if is_star and star_object is not None:
            _, body_lon, _ = (
                earth.at(t).observe(star_object).apparent().ecliptic_latlon()
            )
        else:
            _, body_lon, _ = earth.at(t).observe(target).apparent().ecliptic_latlon()
        return (body_lon.degrees - sun_lon.degrees + 180.0) % 360.0 - 180.0

    def _event_side_ok(jd: float, morning: bool) -> bool:
        """Whether a visibility at ``jd`` is on the side the event requires.

        Morning events (heliacal rising, morning last) require the body WEST of
        the Sun; evening events (heliacal setting, evening first) require it
        EAST. Rejects an opposite-apparition detection that a spilled twilight
        window would otherwise accept.
        """
        lon_diff = _sun_body_lon_diff(jd)
        return (lon_diff < 0.0) if morning else (lon_diff > 0.0)

    def _get_body_magnitude(jd: float) -> float:
        """Get the visual magnitude of the body."""
        if is_star:
            # For fixed stars, return the catalog magnitude
            return star_magnitude
        try:
            pheno = pheno_ut(jd, body, _heliacal_eph_flags(flags))
            return pheno[4]  # Visual magnitude
        except (ValueError, TypeError, ArithmeticError):
            return 0.0  # Default to bright magnitude

    # Create Schaefer model for visibility calculations (observer tuple
    # already parsed above, before the name was rebound).
    schaefer = create_schaefer_model(
        pressure=pressure,
        temperature=temperature,
        humidity=humidity * 100.0 if humidity <= 1.0 else humidity,  # Convert to %
        altitude=altitude,
        met_range=met_range,
        observer_age=_obs_age,
        snellen=_obs_snellen,
        binocular=_obs_bino,
        telescope_mag=_obs_mag,
        aperture=_obs_aper,
        transmission=_obs_trans,
        latitude=lat,
        jd=jd_start,
    )

    def _get_moon_data(jd: float) -> Tuple[float, float, float]:
        """Get Moon altitude, phase, and angular separation from body."""
        t = ts.ut1_jd(jd)
        observer_at = earth + observer

        # Moon position
        moon = eph["moon"]
        moon_app = observer_at.at(t).observe(moon).apparent()
        moon_alt, moon_az, _ = moon_app.altaz()

        # Moon phase (0 = new, 1 = full)
        sun_pos = earth.at(t).observe(sun).apparent()
        moon_geo = earth.at(t).observe(moon).apparent()
        elongation_moon = moon_geo.separation_from(sun_pos).degrees
        phase = (1.0 - math.cos(math.radians(elongation_moon))) / 2.0

        # Angular separation between body and Moon
        if is_star and star_object is not None:
            body_app = observer_at.at(t).observe(star_object).apparent()
        else:
            body_app = observer_at.at(t).observe(target).apparent()
        moon_body_sep = body_app.separation_from(moon_app).degrees

        return moon_alt.degrees, phase, moon_body_sep

    def _calculate_limiting_magnitude(
        sun_alt: float, body_alt: float, jd: float
    ) -> float:
        """
        Calculate the limiting magnitude using Schaefer (1990) model.

        Uses complete atmospheric model considering:
        - Twilight sky brightness
        - Moon contribution
        - Atmospheric extinction
        - Observer characteristics
        """
        # Get Moon data for more accurate calculation
        moon_alt, moon_phase, moon_body_sep = _get_moon_data(jd)
        elongation = _get_elongation(jd)

        return schaefer.limiting_magnitude(
            sun_alt=sun_alt,
            moon_alt=moon_alt,
            moon_phase=moon_phase,
            obj_alt=body_alt,
            sun_obj_angle=elongation,
            moon_obj_angle=moon_body_sep,
        )

    def _is_body_visible(jd: float) -> Tuple[bool, float, float, float]:
        """
        Check if body is visible at given time using Schaefer model.

        Returns: (is_visible, sun_alt, body_alt, elongation)
        """
        sun_alt, body_alt, body_az = _get_altitudes(jd)
        elongation = _get_elongation(jd)

        # Body must be above minimum altitude
        if body_alt < 0:
            return False, sun_alt, body_alt, elongation

        # Sun must be below horizon
        if sun_alt > 0:
            return False, sun_alt, body_alt, elongation

        # Get body magnitude
        body_mag = _get_body_magnitude(jd)

        # Get Moon data
        moon_alt, moon_phase, moon_body_sep = _get_moon_data(jd)

        # Use Schaefer model for visibility check (with the day's season)
        schaefer.update_season(jd)
        is_visible = schaefer.is_visible(
            body_alt=body_alt,
            body_mag=body_mag,
            sun_alt=sun_alt,
            elongation=elongation,
            moon_alt=moon_alt,
            moon_phase=moon_phase,
            moon_obj_angle=moon_body_sep,
            margin=_HELIACAL_VIS_MARGIN,
        )

        return is_visible, sun_alt, body_alt, elongation

    def _is_body_visible_no_moon(jd: float, margin: float = 0.0) -> bool:
        """
        Check if body would be visible ignoring moonlight.

        Used for heliacal transition detection: moonlight can cause
        temporary multi-day invisibility that should not be confused
        with a conjunction passage. By checking visibility without
        Moon contribution, we detect only Sun-body geometry transitions.

        Args:
            jd: Julian day to check
            margin: Detection margin in magnitudes (default 0.5)
        """
        sun_alt, body_alt, _ = _get_altitudes(jd)
        elongation = _get_elongation(jd)

        if body_alt < 0 or sun_alt > 0:
            return False

        body_mag = _get_body_magnitude(jd)

        # Check visibility with Moon forced below horizon (day's season)
        schaefer.update_season(jd)
        return schaefer.is_visible(
            body_alt=body_alt,
            body_mag=body_mag,
            sun_alt=sun_alt,
            elongation=elongation,
            moon_alt=-90.0,
            moon_phase=0.0,
            moon_obj_angle=180.0,
            margin=margin,
        )

    def _find_twilight_time(jd: float, sun_target_alt: float, rising: bool) -> float:
        """
        Find when Sun crosses target altitude (morning or evening).

        Args:
            jd: Starting JD
            sun_target_alt: Target Sun altitude (negative for below horizon)
            rising: True for morning (Sun rising), False for evening (Sun setting)
        """
        # Search within one day
        for _ in range(50):
            sun_alt, _, _ = _get_altitudes(jd)

            # Check if we're at the right phase of day
            if rising:
                # Morning: looking for Sun rising through target altitude
                if abs(sun_alt - sun_target_alt) < 0.01:
                    return jd
            else:
                # Evening: looking for Sun setting through target altitude
                if abs(sun_alt - sun_target_alt) < 0.01:
                    return jd

            # Adjust time based on Sun's position
            # Sun moves ~15/hour = 0.25/minute = 360/day
            sun_rate = 360.0 / 1.0  # degrees per day
            diff = sun_target_alt - sun_alt

            # Crude estimate of time adjustment
            dt = diff / sun_rate

            # Limit step size
            dt = max(-0.1, min(0.1, dt))

            if abs(dt) < 1e-6:
                return jd

            jd += dt

        return jd

    def _find_twilight_center(jd_day: float, morning: bool) -> float:
        """
        Find the approximate UT hour when the Sun is near a target
        depression angle for heliacal scanning.

        Instead of using fixed local-hour windows (which fail at high
        latitudes where dawn can be at 1 AM local in summer), this
        dynamically locates the twilight window by scanning the full
        24-hour period for when the Sun is between 0° and -20°.

        Args:
            jd_day: JD at 0h UT of the day
            morning: True for pre-sunrise twilight, False for post-sunset

        Returns:
            UT fractional hour of the best scan center, or -1 if no
            valid twilight found (e.g. polar summer, sun never below -3°).
        """
        best_ut = -1.0
        best_sun = -999.0

        # Scan 0-24h UT in 1-hour steps. Use a wide acceptance range
        # (-22° to +2°) so we never miss the twilight window even
        # when the Sun transitions rapidly (equatorial latitudes).
        for h in range(24):
            jd_check = jd_day + h / 24.0
            sun_alt, _, _ = _get_altitudes(jd_check)

            if -22.0 < sun_alt < 2.0:
                # Check if Sun is rising or setting
                jd_next = jd_day + (h + 1) / 24.0
                sun_next, _, _ = _get_altitudes(jd_next)

                if morning and sun_next > sun_alt:
                    # Sun is rising → morning twilight region
                    score = -abs(sun_alt + 8.0)
                    if best_ut < 0 or score > best_sun:
                        best_sun = score
                        best_ut = float(h)
                elif not morning and sun_next < sun_alt:
                    # Sun is setting → evening twilight region
                    score = -abs(sun_alt + 8.0)
                    if best_ut < 0 or score > best_sun:
                        best_sun = score
                        best_ut = float(h)

        return best_ut

    def _check_twilight_visibility(jd_day: float, morning: bool) -> Tuple[bool, float]:
        """
        Check if the body is visible during twilight on the given day,
        ignoring moonlight for transition-detection robustness.

        Uses dynamic twilight scanning: first finds the twilight center,
        then samples ±3 hours around it in 15-minute steps.

        Moon brightness is excluded because a bright Moon near full can
        create multi-day fake invisible streaks that would be mistaken
        for conjunction passages. Heliacal events are fundamentally
        about Sun-body geometry, not moonlight.

        The decision uses the neutral zero-magnitude margin produced by the
        independent Falchi/Schaefer physical visibility model; it has no
        fitted dawn/dusk or elongation threshold.

        Args:
            jd_day: JD at 0h UT of the day
            morning: True for dawn, False for dusk

        Returns:
            (visible, jd_best): whether body was found visible, and
            the JD of the first visibility moment found.
        """
        # Delegate to the batched detector so there is exactly ONE decision
        # function per backend (and the batch twin mirrors the LEB detector
        # sample-for-sample). A former hand-rolled loop here evaluated the
        # magnitude per scan point instead of per day and used the
        # geocentric pheno elongation instead of the alt/az-derived one,
        # which made this variant disagree with the batch twin on razor-thin
        # margins.
        return _batch_check_twilight_visibility([jd_day], morning)[0]

    # --- Vectorized batch visibility computation (~30-40x speedup) ---
    # Instead of calling Skyfield once per time-point in a Python loop,
    # pass arrays of JDs through Skyfield in a single vectorized call.

    _HELIACAL_BATCH = 100  # days per batch

    def _batch_check_twilight_visibility(jd_days_list: list, morning: bool) -> list:
        """Vectorized twilight visibility check for multiple days.

        Replaces per-day calls to ``_check_twilight_visibility`` with
        batched Skyfield computation, achieving ~30-40x speedup.

        Returns list of ``(visible, jd_best)`` tuples, one per input day.
        """
        n_days = len(jd_days_list)
        if n_days == 0:
            return []

        jd_days_arr = np.asarray(jd_days_list, dtype=np.float64)
        # Anchor the scan grid to 0h UT of each day's civil day. The search
        # loops pass jd_start + k, carrying jd_start's arbitrary time-of-day
        # into every scanned day, so the 15-min sampling grid shifted with
        # the requested start time and a few minutes' difference in jd_start
        # could flip a marginal detection (the first days of an apparition
        # are visible only for minutes), flipping in turn the apparition
        # classification of the whole search. A fixed per-day anchor makes
        # each day's verdict a function of the day alone.
        jd_days_arr = np.floor(jd_days_arr + 0.5) - 0.5
        sun_upper = -5.0 if morning else -2.0
        observer_at = earth + observer

        # -- Phase 1: find twilight centres (sun alt at 25 hourly points) --
        hours = np.arange(25, dtype=np.float64) / 24.0  # 0h..24h
        hourly_jds = (jd_days_arr[:, np.newaxis] + hours[np.newaxis, :]).ravel()

        sun_alts_grid = (
            observer_at.at(ts.ut1_jd(hourly_jds))
            .observe(sun)
            .apparent()
            .altaz()[0]
            .degrees
        ).reshape(n_days, 25)

        twilight_centers = np.full(n_days, -1.0)
        for i in range(n_days):
            best_ut = -1.0
            best_score = -999.0
            for h in range(24):
                sa = float(sun_alts_grid[i, h])
                if -22.0 < sa < 2.0:
                    sa_next = float(sun_alts_grid[i, h + 1])
                    if morning and sa_next > sa:
                        score = -abs(sa + 8.0)
                        if best_ut < 0 or score > best_score:
                            best_score = score
                            best_ut = float(h)
                    elif not morning and sa_next < sa:
                        score = -abs(sa + 8.0)
                        if best_ut < 0 or score > best_score:
                            best_score = score
                            best_ut = float(h)
            twilight_centers[i] = best_ut

        # -- Phase 2: build scan JDs around each valid twilight centre --
        scan_offsets_min = np.arange(-180, 181, 15, dtype=np.float64)
        n_scans = len(scan_offsets_min)

        valid_indices = np.where(twilight_centers >= 0)[0]
        if len(valid_indices) == 0:
            return [(False, 0.0)] * n_days

        scan_jds_list: list[float] = []
        scan_day_map: list[int] = []
        for idx in valid_indices:
            center = twilight_centers[idx]
            for s in range(n_scans):
                ut_h = center + scan_offsets_min[s] / 60.0
                scan_jds_list.append(float(jd_days_arr[idx]) + ut_h / 24.0)
                scan_day_map.append(int(idx))

        scan_jds = np.asarray(scan_jds_list, dtype=np.float64)
        scan_day_idx = np.asarray(scan_day_map, dtype=np.intp)

        # -- Phase 3: batch-compute positions at all scan points --
        # Only 2 vectorised Skyfield calls instead of 4: compute alt+az
        # for both Sun and body, then derive elongation from those coords.
        t_scan = ts.ut1_jd(scan_jds)

        sun_app_scan = observer_at.at(t_scan).observe(sun).apparent()
        sun_altaz = sun_app_scan.altaz()
        sun_alts_scan = sun_altaz[0].degrees
        sun_azs_scan = sun_altaz[1].degrees

        _body_target = star_object if (is_star and star_object is not None) else target
        body_app_scan = observer_at.at(t_scan).observe(_body_target).apparent()
        body_altaz = body_app_scan.altaz()
        body_alts_scan = body_altaz[0].degrees
        body_azs_scan = body_altaz[1].degrees

        # Elongation from topocentric alt/az (avoids 2 extra geocentric calls)
        _sa_r = np.radians(np.asarray(sun_alts_scan))
        _ba_r = np.radians(np.asarray(body_alts_scan))
        _daz_r = np.radians(np.asarray(sun_azs_scan) - np.asarray(body_azs_scan))
        _cos_elong = np.sin(_sa_r) * np.sin(_ba_r) + np.cos(_sa_r) * np.cos(
            _ba_r
        ) * np.cos(_daz_r)
        np.clip(_cos_elong, -1.0, 1.0, out=_cos_elong)
        elongations = np.degrees(np.arccos(_cos_elong))

        # -- Phase 4: body magnitude --
        # EXACT per-day magnitude at 0h UT of each civil day, matching the
        # LEB detector twin sample-for-sample. An earlier 5-day
        # interpolation saved ~1 ms/day of pheno_ut but introduced up to
        # ~0.01 mag of curvature error near Mercury's fast light-curve
        # turns — enough to flip a razor-thin (~0.03 mag) twilight
        # visibility margin by a day RELATIVE TO THE LEB BACKEND, i.e.
        # backend-dependent heliacal dates. Magnitudes are only computed
        # for days with a valid twilight centre.
        day_mags: dict[int, float] = {}
        if is_star:
            for idx in valid_indices:
                day_mags[int(idx)] = star_magnitude
        else:
            for idx in valid_indices:
                day_mags[int(idx)] = _get_body_magnitude(float(jd_days_arr[int(idx)]))

        # -- Phase 5: check Schaefer visibility per scan point --
        results: list[tuple[bool, float]] = [(False, 0.0)] * n_days
        found: set[int] = set()

        for k in range(len(scan_jds)):
            day_i = int(scan_day_idx[k])
            if day_i in found:
                continue

            sa = float(sun_alts_scan[k])
            ba = float(body_alts_scan[k])

            if not (-18.0 < sa < sun_upper and ba > 0.5):
                continue
            if ba < 0 or sa > 0:
                continue

            schaefer.update_season(float(scan_jds[k]))
            visible = schaefer.is_visible(
                body_alt=ba,
                body_mag=day_mags[day_i],
                sun_alt=sa,
                elongation=float(elongations[k]),
                moon_alt=-90.0,
                moon_phase=0.0,
                moon_obj_angle=180.0,
                margin=_HELIACAL_VIS_MARGIN,
            )
            if visible and _event_side_ok(float(scan_jds[k]), morning):
                results[day_i] = (True, float(scan_jds[k]))
                found.add(day_i)

        return results

    def _search_heliacal_rising(jd_start: float) -> float:
        """Search for heliacal rising using batched Skyfield computation.

        A heliacal rising is the first morning visibility after the body
        emerged from a conjunction gap (was too close to the Sun). Requiring
        only N consecutive invisible days is insufficient: mid-apparition a
        bright planet can flicker invisible for a few days (Moon glare /
        marginal twilight), which would be mis-reported as a rising. So the
        invisibility streak must also have passed through low solar elongation
        (a real conjunction gap), which a mid-apparition dip never does.

        One legitimate case the elongation test alone misses: when the
        conjunction happened *before* the search window, the body starts the
        search already past conjunction but still invisible (e.g. Mars at 17 deg
        and climbing), so the streak never dips below ELONG_GAP inside the
        window even though the first visibility is a true rising. Accept it when
        the body has not been visible at all since the search start
        (seen_visible False) -- a mid-apparition flicker, visible before the
        dip, never satisfies that.
        """
        ELONG_GAP = 10.0  # deg
        ELONG_START_GAP = 30.0
        max_days = 800

        # Look back to establish initial visibility state.
        lookback_jds = [jd_start - i for i in range(1, 7)]
        lookback_vis = _batch_check_twilight_visibility(lookback_jds, morning=True)
        consecutive_invisible = 0
        min_elong = 999.0
        seen_visible = False
        for jd_lb, (vis, _) in zip(lookback_jds, lookback_vis):
            if not vis:
                consecutive_invisible += 1
                min_elong = min(min_elong, _get_elongation(jd_lb))
            else:
                seen_visible = True
                break
        if not seen_visible and min_elong > ELONG_GAP:
            extended_jds = [jd_start - i for i in range(7, 31)]
            extended_vis = _batch_check_twilight_visibility(extended_jds, morning=True)
            for jd_lb, (vis, _) in zip(extended_jds, extended_vis):
                if vis:
                    seen_visible = True
                    break
                min_elong = min(min_elong, _get_elongation(jd_lb))

        for batch_start in range(0, max_days, _HELIACAL_BATCH):
            batch_end = min(batch_start + _HELIACAL_BATCH, max_days)
            jd_days = [jd_start + d for d in range(batch_start, batch_end)]
            vis_results = _batch_check_twilight_visibility(jd_days, morning=True)

            for (vis, jd_vis), jd_day in zip(vis_results, jd_days):
                if not vis:
                    consecutive_invisible += 1
                    # use the day JD (always valid); jd_vis is a sentinel when
                    # there is no twilight visibility window.
                    min_elong = min(min_elong, _get_elongation(jd_day))
                else:
                    current_elong = _get_elongation(jd_day)
                    pre_window_gap = (
                        not seen_visible
                        and not (flags & HELFLAG_HIGH_PRECISION)
                        and min_elong <= ELONG_START_GAP
                        and current_elong >= min_elong + 2.0
                    )
                    if consecutive_invisible >= 5 and (
                        is_star or min_elong <= ELONG_GAP or pre_window_gap
                    ):
                        return _refine_heliacal_time(jd_vis, is_morning=True)
                    consecutive_invisible = 0
                    min_elong = 999.0
                    seen_visible = True

        return 0.0

    def _search_heliacal_setting(jd_start: float) -> float:
        """Search for heliacal setting using batched Skyfield computation."""
        max_days = 800
        last_visible_jd = 0.0
        consecutive_invisible = 0

        for batch_start in range(0, max_days, _HELIACAL_BATCH):
            batch_end = min(batch_start + _HELIACAL_BATCH, max_days)
            jd_days = [jd_start + d for d in range(batch_start, batch_end)]
            vis_results = _batch_check_twilight_visibility(jd_days, morning=False)

            for vis, jd_vis in vis_results:
                if vis:
                    last_visible_jd = jd_vis
                    consecutive_invisible = 0
                else:
                    consecutive_invisible += 1
                    if consecutive_invisible >= 5 and last_visible_jd > 0:
                        return _refine_heliacal_time(last_visible_jd, is_morning=False)

        if last_visible_jd > 0:
            return _refine_heliacal_time(last_visible_jd, is_morning=False)
        return 0.0

    def _search_evening_first(jd_start: float) -> float:
        """Search for evening first visibility using batched computation."""
        ELONG_GAP = 10.0
        ELONG_START_GAP = 30.0
        max_days = 800

        # Fast-path only when the body is near conjunction (elongation below the
        # gap) AND invisible both this evening and on the 6 preceding evenings:
        # that is a genuine conjunction the search starts inside, so the next
        # evening visibility IS the evening-first event. If the body is visible
        # tonight (ongoing apparition) or was visible on a recent evening (a
        # marginal-conditions flicker inside the current apparition), this is
        # not a conjunction gap and we must fall through to the normal search
        # for the NEXT event after the coming superior conjunction.
        if not is_star and _get_elongation(jd_start) <= ELONG_GAP:
            gate_jds = [jd_start - i for i in range(0, 7)]
            gate_vis = _batch_check_twilight_visibility(gate_jds, morning=False)
            if not any(vis for vis, _ in gate_vis):
                for batch_start in range(1, max_days, _HELIACAL_BATCH):
                    batch_end = min(batch_start + _HELIACAL_BATCH, max_days)
                    jd_days = [jd_start + d for d in range(batch_start, batch_end)]
                    vis_results = _batch_check_twilight_visibility(
                        jd_days, morning=False
                    )
                    for vis, jd_vis in vis_results:
                        if vis:
                            return _refine_heliacal_time(jd_vis, is_morning=False)
                return 0.0

        lookback_jds = [jd_start - i for i in range(1, 7)]
        lookback_vis = _batch_check_twilight_visibility(lookback_jds, morning=False)
        consecutive_invisible = 0
        min_elong = 999.0
        seen_visible = False
        for jd_lb, (vis, _) in zip(lookback_jds, lookback_vis):
            if not vis:
                consecutive_invisible += 1
                min_elong = min(min_elong, _get_elongation(jd_lb))
            else:
                seen_visible = True
                break

        for batch_start in range(0, max_days, _HELIACAL_BATCH):
            batch_end = min(batch_start + _HELIACAL_BATCH, max_days)
            jd_days = [jd_start + d for d in range(batch_start, batch_end)]
            vis_results = _batch_check_twilight_visibility(jd_days, morning=False)

            for (vis, jd_vis), jd_day in zip(vis_results, jd_days):
                if not vis:
                    consecutive_invisible += 1
                    min_elong = min(min_elong, _get_elongation(jd_day))
                else:
                    current_elong = _get_elongation(jd_day)
                    pre_window_gap = (
                        not seen_visible
                        and not (flags & HELFLAG_HIGH_PRECISION)
                        and min_elong <= ELONG_START_GAP
                        and current_elong >= min_elong + 2.0
                    )
                    if consecutive_invisible >= 5 and (
                        is_star or min_elong <= ELONG_GAP or pre_window_gap
                    ):
                        return _refine_heliacal_time(jd_vis, is_morning=False)
                    consecutive_invisible = 0
                    min_elong = 999.0
                    seen_visible = True

        return 0.0

    def _search_morning_last(jd_start: float) -> float:
        """Search for morning last visibility using batched computation."""
        ELONG_GAP = 10.0
        max_days = 800
        last_visible_jd = 0.0
        found_visible = False
        consecutive_invisible = 0
        streak_min_elong = 999.0

        for batch_start in range(0, max_days, _HELIACAL_BATCH):
            batch_end = min(batch_start + _HELIACAL_BATCH, max_days)
            jd_days = [jd_start + d for d in range(batch_start, batch_end)]
            vis_results = _batch_check_twilight_visibility(jd_days, morning=True)

            for (vis, jd_vis), jd_day in zip(vis_results, jd_days):
                if vis:
                    last_visible_jd = jd_vis
                    found_visible = True
                    consecutive_invisible = 0
                    streak_min_elong = 999.0
                elif found_visible:
                    consecutive_invisible += 1
                    streak_min_elong = min(streak_min_elong, _get_elongation(jd_day))
                    if (
                        consecutive_invisible >= 3
                        and streak_min_elong <= ELONG_GAP
                        and last_visible_jd > 0
                    ):
                        return _refine_heliacal_time(last_visible_jd, is_morning=True)

        return last_visible_jd if last_visible_jd > 0 else 0.0

    def _refine_heliacal_time(jd_approx: float, is_morning: bool) -> float:
        """
        Refine the heliacal event time using binary search.

        Find the exact moment when the body becomes just visible/invisible.
        """
        # Use binary search to find the transition point
        jd_low = jd_approx - 0.1  # ~2.4 hours before
        jd_high = jd_approx + 0.1  # ~2.4 hours after

        for _ in range(30):  # ~30 iterations gives very high precision
            jd_mid = (jd_low + jd_high) / 2

            visible, _sun_alt, _body_alt, _elong = _is_body_visible(jd_mid)

            # For heliacal rising: looking for first visibility
            # For heliacal setting: looking for last visibility
            if is_morning:
                # Morning: visibility increases as time progresses (Sun rises)
                # But we want first visibility, so search backwards
                if visible:
                    jd_high = jd_mid  # Look earlier
                else:
                    jd_low = jd_mid  # Look later
            else:
                # Evening: visibility decreases as time progresses (sky darkens)
                # For last visibility, we want the last moment still visible
                if visible:
                    jd_low = jd_mid  # Look later for last visibility
                else:
                    jd_high = jd_mid  # Look earlier

            if jd_high - jd_low < 1e-6:  # ~0.1 second precision
                break

        return (jd_low + jd_high) / 2

    # Main search logic based on event type
    if event_type == HELIACAL_RISING:
        _search = _search_heliacal_rising
    elif event_type == HELIACAL_SETTING:
        _search = _search_heliacal_setting
    elif event_type == EVENING_FIRST:
        _search = _search_evening_first
    elif event_type == MORNING_LAST:
        _search = _search_morning_last
    else:
        _search = None

    jd_event = 0.0
    if _search is not None:
        jd_event = _search(jd_start)
        # Contract gate: the returned event must be >= jd_start. Day 0 of
        # the search is the civil day containing jd_start and its twilight
        # scan covers that whole day, so when jd_start falls later in the
        # day than the visibility optimum, the streak gate can classify the
        # already-past optimum as the event and the refinement converges
        # before jd_start. Restarting one day later moves that day into the
        # lookback, where its visibility marks the apparition as ongoing --
        # the same branch a search started the following day takes -- and
        # the search skips to the next apparition. The retry cannot violate
        # the contract again: every instant it scans lies on civil days
        # strictly after jd_start's.
        if 0.0 < jd_event < jd_start - 1e-6:
            jd_event = _search(jd_start + 1.0)

    if jd_event > 0:
        return jd_event, event_type
    else:
        return 0.0, -1  # Not found


def heliacal_ut(
    tjdut: float,
    geopos: tuple,
    atmo: tuple,
    observer: tuple,
    objname: str,
    eventtype: int,
    flags: int = FLG_SWIEPH,
) -> Tuple[float, float, float]:
    """Find the next requested heliacal event after a UT start date.

    Args:
        tjdut: Search start as a Julian Day in UT.
        geopos: Longitude, latitude, and observer altitude.
        atmo: Pressure, temperature, humidity, and either meteorological
            range or a total atmospheric coefficient.
        observer: Observer age and visual-acuity values, followed optionally
            by binocular/telescope parameters.
        objname: Planet or fixed-star name. The Sun and Moon are invalid for
            this API.
        eventtype: ``HELIACAL_RISING``, ``HELIACAL_SETTING``,
            ``EVENING_FIRST``, or ``MORNING_LAST``.
        flags: Ephemeris and heliacal calculation flags.

    Returns:
        Three Julian Days giving visibility start, optimum, and end. Detail
        slots may be zero when ``HELFLAG_NO_DETAILS`` is set.

    Raises:
        ValueError: If the object or event type is invalid.

    Notes:
        The visibility search combines body altitude, solar twilight, object
        magnitude, atmospheric extinction, sky brightness, and observer
        characteristics.
    """
    from .constants import (
        HELIACAL_RISING,
        HELIACAL_SETTING,
        EVENING_FIRST,
        MORNING_LAST,
        HELFLAG_NO_DETAILS,
    )

    # Validate event type
    if eventtype not in (
        HELIACAL_RISING,
        HELIACAL_SETTING,
        EVENING_FIRST,
        MORNING_LAST,
    ):
        raise Error(
            f"invalid event type {eventtype}; use HELIACAL_RISING, "
            "HELIACAL_SETTING, EVENING_FIRST or MORNING_LAST"
        )

    # Parse geopos
    lon = geopos[0] if len(geopos) > 0 else 0.0
    lat = geopos[1] if len(geopos) > 1 else 0.0
    altitude = geopos[2] if len(geopos) > 2 else 0.0

    # Parse the atmospheric tuple with the public input defaults.
    pressure = atmo[0] if len(atmo) > 0 and atmo[0] > 0 else 1013.25
    temperature = atmo[1] if len(atmo) > 1 else 15.0
    humidity_pct = atmo[2] if len(atmo) > 2 else 40.0
    # atmo[3] is the meteorological range in km (or ktot if 0 < value < 1).
    # It drives atmospheric extinction and materially shifts the event date,
    # so thread it into the search (0.0 selects the clear-sky default).
    met_range = atmo[3] if len(atmo) > 3 else 0.0

    # datm[2] is relative humidity in percent in the public input tuple;
    # convert to the 0-1 range used internally.
    humidity = humidity_pct / 100.0

    # Parse objname to get body ID. The Moon is a valid heliacal body for the
    # evening-first / morning-last events (crescent visibility) but not for the
    # rising/setting events under the public accepted-body contract.
    from .constants import MOON

    body_id = _parse_object_name(objname, allow_moon=True)
    if body_id == MOON and eventtype in (HELIACAL_RISING, HELIACAL_SETTING):
        # Event names follow the archaeoastronomy usage in which the heliacal
        # rising is the "morning first" and the heliacal setting the "evening
        # last" sighting of a body in its apparition cycle: Schaefer, B.E.
        # (1987), "Heliacal Rise Phenomena", Journal for the History of
        # Astronomy 18, Archaeoastronomy Suppl. 11, S19-S33; Purrington, R.D.
        # (1988), "Heliacal Rising and Setting: Quantitative Aspects", Journal
        # for the History of Astronomy 19, Archaeoastronomy Suppl. 12, S72-S85.
        ev_name = {
            HELIACAL_RISING: "morning first",
            HELIACAL_SETTING: "evening last",
        }[eventtype]
        raise Error(
            f"the {ev_name} event (type {eventtype}) is not defined for the Moon"
        )

    # Call the internal _heliacal_ut_pythonic function
    jd_event, retflag = _heliacal_ut_pythonic(
        jd_start=tjdut,
        lat=lat,
        lon=lon,
        altitude=altitude,
        pressure=pressure,
        temperature=temperature,
        humidity=humidity,
        body=body_id,
        event_type=eventtype,
        flags=flags,
        met_range=met_range,
        observer=observer,
    )

    # Build the public three-float result.
    jd1 = 0.0  # Start of visibility
    jd2 = 0.0  # Optimum visibility
    jd3 = 0.0  # End of visibility

    if jd_event > 0:
        jd1 = jd_event

        if not (flags & HELFLAG_NO_DETAILS):
            # Construct the optimum (maximum visibility margin) and the two
            # crossings where that margin vanishes.
            jd1, jd2, jd3 = _heliacal_visibility_window(
                jd_event, geopos, atmo, observer, objname, flags
            )

    return (jd1, jd2, jd3)


def _heliacal_visibility_window(
    jd_event: float,
    geopos: tuple,
    atmo: tuple,
    observer: tuple,
    objname: str,
    flags: int,
) -> Tuple[float, float, float]:
    """(start, optimum, end) of a heliacal visibility window.

    The visibility margin is the limiting visual magnitude minus the
    body's apparent magnitude with extinction; the body is visible
    where the margin is positive. The optimum maximizes the margin, the
    window limits are its zero crossings on either side. When a
    crossing cannot be bracketed (e.g. the object stays visible into
    darkness), the event time itself is reported for that limit.
    """

    def _margin(jd: float) -> float:
        # dret[0] is the limiting (catalog) magnitude at the object, with
        # atmospheric extinction already folded in, so the visibility
        # margin is simply limiting mag minus the body's catalog magnitude.
        # A heliacal event is a *twilight* phenomenon: reject a Sun above
        # the horizon so the optimum/window search cannot climb into the
        # daytime sky (a very bright object such as Venus can otherwise
        # register a secondary daylight visibility maximum).
        retval, dret = vis_limit_mag(jd, geopos, atmo, observer, objname, flags)
        if retval < 0 or dret[3] >= 0.0:
            return -99.0
        return dret[0] - dret[7]

    # Optimum: the maximum margin within ~3 hours of the event. The margin
    # is a narrow twilight peak flanked by flat -99 plateaus (the body below
    # the horizon before dawn / after dusk, and the daytime sky); a unimodal
    # golden-section search is not robust to that shape, so scan the window on
    # a ~2.5-min grid and take the maximum, then refine locally.
    win = 0.12
    n_steps = 144
    jd_opt = jd_event
    best = _margin(jd_event)
    for i in range(n_steps + 1):
        t = jd_event - win + (2.0 * win) * i / n_steps
        m = _margin(t)
        if m > best:
            best = m
            jd_opt = t
    # Golden-section refine within +/- one grid step of the scan maximum.
    step = (2.0 * win) / n_steps
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    lo, hi = jd_opt - step, jd_opt + step
    a = hi - (hi - lo) / phi
    b = lo + (hi - lo) / phi
    fa, fb = _margin(a), _margin(b)
    for _ in range(20):
        if fa > fb:
            hi, b, fb = b, a, fa
            a = hi - (hi - lo) / phi
            fa = _margin(a)
        else:
            lo, a, fa = a, b, fb
            b = lo + (hi - lo) / phi
            fb = _margin(b)
        if hi - lo < 2e-5:
            break
    if _margin(0.5 * (lo + hi)) > best:
        jd_opt = 0.5 * (lo + hi)

    def _crossing(t_in: float, t_out: float) -> float:
        # t_in is always the optimum (visible, margin > 0); t_out is the far
        # edge of the search bracket on one side of it. The return value is
        # bracketed strictly between t_out and t_in, so start <= opt <= end is
        # guaranteed by construction: on a failed bracket we clamp to the
        # optimum (degenerate window) or to the search edge, never to an
        # unrelated instant outside the interval.
        f_in = _margin(t_in)
        f_out = _margin(t_out)
        if f_in <= 0.0:
            # Optimum itself is not visible: report a degenerate window at it.
            return t_in
        if f_out > 0.0:
            # Still visible at the search edge: the crossing is beyond it, so
            # report the edge (keeps ordering; widths saturate at the bracket).
            return t_out
        for _ in range(25):
            mid = 0.5 * (t_in + t_out)
            if _margin(mid) > 0.0:
                t_in = mid
            else:
                t_out = mid
            if abs(t_out - t_in) < 2e-5:
                break
        return 0.5 * (t_in + t_out)

    jd_start_vis = _crossing(jd_opt, jd_opt - win)
    jd_end_vis = _crossing(jd_opt, jd_opt + win)
    # Guard against any residual ordering inversion from a degenerate search.
    lo_vis = min(jd_start_vis, jd_opt)
    hi_vis = max(jd_end_vis, jd_opt)
    return lo_vis, jd_opt, hi_vis


def _parse_object_name(
    object_name: str, allow_moon: bool = False, allow_sun: bool = False
) -> int:
    """
    Parse an object name string to get the corresponding body ID.

    Args:
        object_name: Name of planet or fixed star (e.g., "Venus", "Sirius")
        allow_moon: If True, allow Moon as a valid body (for _heliacal_pheno_ut_pythonic
            which supports Moon crescent calculations)
        allow_sun: If True, allow the Sun for ``heliacal_pheno_ut``;
            ``heliacal_ut`` and ``vis_limit_mag`` keep the default rejection.

    Returns:
        Body ID (body constant) for planets

    Raises:
        Error: If object name is not recognized or not valid for heliacal
    """
    # Import planet constants
    from .constants import (
        SUN,
        MOON,
        MERCURY,
        VENUS,
        MARS,
        JUPITER,
        SATURN,
        URANUS,
        NEPTUNE,
        PLUTO,
    )

    if not isinstance(object_name, str):
        parsed_body_id = int(object_name)
        if parsed_body_id == SUN and not allow_sun:
            raise Error("the Sun is not a valid object for heliacal calculations")
        if parsed_body_id == MOON and not allow_moon:
            raise Error("the Moon is not a valid object for this heliacal calculation")
        return parsed_body_id

    # Normalize name
    name_upper = object_name.upper().strip()

    # Planet name to ID mapping
    planet_map = {
        "SUN": SUN,
        "MOON": MOON,
        "MERCURY": MERCURY,
        "VENUS": VENUS,
        "MARS": MARS,
        "JUPITER": JUPITER,
        "SATURN": SATURN,
        "URANUS": URANUS,
        "NEPTUNE": NEPTUNE,
        "PLUTO": PLUTO,
    }

    # Check if it's a planet
    if name_upper in planet_map:
        body_id = planet_map[name_upper]

        # Sun is valid only for the public phenomena result when requested by
        # its wrapper.
        if body_id == SUN and not allow_sun:
            raise Error("the Sun is not a valid object for heliacal calculations")
        if body_id == MOON and not allow_moon:
            raise Error("the Moon is not a valid object for this heliacal calculation")

        return body_id

    # Try to resolve as a fixed star name
    from .fixed_stars import resolve_star_name

    star_id = resolve_star_name(object_name)
    if star_id is not None:
        return star_id

    # Neither a planet nor a catalogued star.
    raise Error(
        f"unknown object name {object_name!r}; use a planet name (Mercury, "
        "Venus, Mars, Jupiter, Saturn, ...), an integer planet id (2-9 for "
        "Mercury through Pluto) or a fixed star name (Sirius, Regulus, "
        "Aldebaran, ...)"
    )


# Placeholder written to the time-valued slots of the heliacal phenomena
# result (rise/set instants, the visibility window and its Yallop optimum)
# when the instant does not exist or could not be located. Same value as the
# public ``TJD_INVALID`` constant; it lies far beyond any Julian Day the
# library can evaluate, so a slot can be tested against it directly.
_TJD_INVALID = 99999999.0


def _heliacal_eph_flags(flags: int) -> int:
    """Reduce a heliacal flags word to its ephemeris-selection bits.

    HELFLAG_* values share bit positions with FLG_* calculation flags
    (e.g. HELFLAG_VISLIM_DARK == FLG_XYZ), so heliacal flags must never
    reach calc/pheno/fixstar calls unmasked.
    """
    from .constants import FLG_JPLEPH, FLG_MOSEPH

    return flags & (FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH)


_PHENO_BODY_NAMES = {
    1: "Moon",
    2: "Mercury",
    3: "Venus",
    4: "Mars",
    5: "Jupiter",
    6: "Saturn",
    7: "Uranus",
    8: "Neptune",
    9: "Pluto",
}


# ---------------------------------------------------------------------------
# The window a heliacal sighting happens in
# ---------------------------------------------------------------------------

#: How far before the instant a caller asks about a horizon crossing may still
#: be the one that date is about, in days. The reach backwards is not
#: cosmetic: the instant is routinely a civil midnight, and the evening event
#: that belongs to that date begins before it.
_PHENO_CROSSING_REACH = 4.0 / 24.0

#: Width, in days, of the twilight a sighting can happen in, measured from the
#: Sun's crossing: backwards for a morning event, forwards for an evening one.
_PHENO_TWILIGHT_WIDTH = 4.0 / 24.0

#: The event types framed by a rise (first and last morning visibility). The
#: other two, 2 and 3, are framed by a set.
_PHENO_MORNING_EVENTS = frozenset((1, 4))

#: The event types that describe an apparition on the far side of a
#: conjunction (first evening and last morning visibility).
_PHENO_FAR_SIDE_EVENTS = frozenset((3, 4))

#: The objects whose apparition on the far side of a conjunction still carries
#: a visibility interval. For anything else those two event types answer no
#: interval and an exactly zero duration.
_PHENO_FAR_SIDE_OBJECTS = frozenset((MOON, MERCURY, VENUS))

#: Fraction of the lag, after the Sun's crossing, at which the lunar crescent
#: is classically best looked for. Bruin, F. (1977), "The First Visibility of
#: the Lunar Crescent", Vistas in Astronomy 21, 331-358, derives the optimum
#: from the competition between the fading twilight and the sinking crescent;
#: Yallop, B.D. (1997), NAO Technical Note No. 69, states it as four ninths.
_PHENO_BEST_TIME_FRACTION = 4.0 / 9.0

#: Number of equal steps the twilight is first walked in. This library's own
#: choice, not a published quantity: half a minute is finer than the narrowest
#: interval the visibility model produces, so a window cannot fall whole
#: between two samples of the scan.
_PHENO_SCAN_STEPS = 480

#: Width in days below which a bracketed instant is converged, about nine
#: milliseconds. Also this library's own choice: it is far below the
#: minute-scale uncertainty of the visibility model that furnishes the margin.
_PHENO_SEARCH_TOL = 1.0e-7


def _pheno_event_is_morning(event_type: int) -> bool:
    """Say whether an event type is framed by a rise rather than by a set.

    Args:
        event_type: Heliacal event type, 1 to 4.

    Returns:
        True when the two horizon crossings of the event are rises.
    """
    return event_type in _PHENO_MORNING_EVENTS


def _pheno_horizon_crossing(
    jd_from: float,
    target: Union[int, str],
    morning: bool,
    geopos3: tuple,
    eph_flags: int,
) -> float:
    """Locate a disc-centre horizon crossing at or after an instant.

    Args:
        jd_from: Earliest Julian Day (UT) the crossing may fall on.
        target: Numeric body id, or the catalogue name of a fixed star.
        morning: True for a rise, False for a set.
        geopos3: Longitude and latitude in degrees, height above sea level in
            metres.
        eph_flags: Ephemeris-selection bits only.

    Returns:
        The Julian Day (UT) of the crossing, or the invalid-time sentinel when
        the object never crosses the horizon at this site - it is either
        circumpolar or permanently below it. The two cases are not
        distinguished.
    """
    from .constants import BIT_DISC_CENTER, CALC_RISE, CALC_SET
    from .eclipse import rise_trans

    rsmi = (CALC_RISE if morning else CALC_SET) | BIT_DISC_CENTER
    # The crossings are geometry: the caller's atmosphere is deliberately not
    # passed on, so the rise/set path applies its own default for the site.
    status, tret = rise_trans(jd_from, target, rsmi, geopos3, flags=eph_flags)
    if status < 0 or not tret or tret[0] <= 0.0:
        return _TJD_INVALID
    return float(tret[0])


def _pheno_visibility_margin(
    tjd: float,
    objname: str,
    geopos3: tuple,
    atmo4: tuple,
    obs6: tuple,
    flags: int,
) -> float:
    """Measure how far the object is from being hidden by the sky, in magnitudes.

    Args:
        tjd: Instant as a Julian Day in UT.
        objname: Planet or catalogue star name the limiting-magnitude entry
            point accepts.
        geopos3: Longitude and latitude in degrees, height above sea level in
            metres.
        atmo4: Pressure in hPa, temperature in Celsius, relative humidity in
            percent, meteorological range or total extinction.
        obs6: Age in years, Snellen ratio, binocular flag, magnification,
            aperture in millimetres, transmission.
        flags: The caller's whole flags word - the visibility options live in
            it and the margin is meant to feel them.

    Returns:
        The limiting magnitude of the sky at the object's place, extinction
        already folded in, minus the object's own magnitude, so positive
        exactly where the object can be seen. Negative infinity where the
        margin has no meaning: the limiting magnitude reports a failure, which
        is what it does for an object below the horizon, or the Sun is above
        the horizon - a heliacal sighting is a twilight phenomenon, and a
        bright object otherwise registers a second, daytime maximum that has
        nothing to do with the event.
    """
    status, detail = vis_limit_mag(tjd, geopos3, atmo4, obs6, objname, flags)
    if status < 0 or detail[3] > 0.0:
        return -math.inf
    return float(detail[0]) - float(detail[7])


def _pheno_brightest_instant(
    margin: Callable[[float], float],
    lo: float,
    hi: float,
    tol: float,
) -> float:
    """Locate the largest margin inside a bracket by golden-section search.

    Args:
        margin: The visibility margin as a function of a Julian Day in UT.
        lo: Lower end of the bracket, a Julian Day in UT.
        hi: Upper end of the bracket, a Julian Day in UT.
        tol: Width in days below which the bracket is small enough.

    Returns:
        A Julian Day (UT) strictly inside the bracket, at or next to the
        largest margin the bracket holds.
    """
    ratio = (math.sqrt(5.0) - 1.0) / 2.0
    left, right = lo, hi
    inner_lo = right - (right - left) * ratio
    inner_hi = left + (right - left) * ratio
    value_lo = margin(inner_lo)
    value_hi = margin(inner_hi)
    while right - left > tol:
        if value_lo >= value_hi:
            right, inner_hi, value_hi = inner_hi, inner_lo, value_lo
            inner_lo = right - (right - left) * ratio
            value_lo = margin(inner_lo)
        else:
            left, inner_lo, value_lo = inner_lo, inner_hi, value_hi
            inner_hi = left + (right - left) * ratio
            value_hi = margin(inner_hi)
    return inner_lo if value_lo >= value_hi else inner_hi


def _pheno_margin_edge(
    margin: Callable[[float], float],
    inside: float,
    outside: float,
    tol: float,
) -> float:
    """Bisect the instant at which the margin ceases to be positive.

    Args:
        margin: The visibility margin as a function of a Julian Day in UT.
        inside: An instant whose margin is positive.
        outside: An instant on the other side of the sign change, whose margin
            is not.
        tol: Width in days below which the bracket is small enough.

    Returns:
        The last Julian Day (UT) on the ``inside`` side that still carries a
        positive margin.
    """
    while abs(outside - inside) > tol:
        middle = 0.5 * (inside + outside)
        if margin(middle) > 0.0:
            inside = middle
        else:
            outside = middle
    return inside


def _pheno_visibility_interval(
    margin: Callable[[float], float],
    first: float,
    last: float,
) -> Tuple[float, float, float]:
    """Cut out the interval in which the sky can show the object.

    The twilight is walked on a uniform grid, the largest margin is refined
    inside the cell that holds it, and the two sign changes next to that
    optimum are bisected. Only those two bound the interval: a margin that
    goes positive again elsewhere in the twilight is another apparition, not
    this one.

    Args:
        margin: The visibility margin as a function of a Julian Day in UT.
        first: Earlier end of the twilight bracket, a Julian Day in UT.
        last: Later end of the twilight bracket, a Julian Day in UT.

    Returns:
        Beginning, optimum and end as Julian Days in UT, or three invalid-time
        sentinels when the margin is never positive anywhere in the bracket.
    """
    steps = _PHENO_SCAN_STEPS
    tol = _PHENO_SEARCH_TOL
    samples = [first + (last - first) * (index / steps) for index in range(steps)]
    samples.append(last)
    values = [margin(instant) for instant in samples]
    peak = max(range(len(samples)), key=values.__getitem__)

    cell_lo = samples[peak - 1] if peak > 0 else samples[peak]
    cell_hi = samples[peak + 1] if peak + 1 < len(samples) else samples[peak]
    optimum, best = samples[peak], values[peak]
    if cell_hi > cell_lo:
        refined = _pheno_brightest_instant(margin, cell_lo, cell_hi, tol)
        refined_margin = margin(refined)
        if refined_margin > best:
            optimum, best = refined, refined_margin
    if best <= 0.0:
        return _TJD_INVALID, _TJD_INVALID, _TJD_INVALID

    if values[peak] > 0.0:
        # The scan already stands inside the interval: walk out along it.
        low = peak
        while low > 0 and values[low - 1] > 0.0:
            low -= 1
        high = peak
        while high + 1 < len(samples) and values[high + 1] > 0.0:
            high += 1
        inner_lo, inner_hi = samples[low], samples[high]
        outer_lo = samples[low - 1] if low > 0 else None
        outer_hi = samples[high + 1] if high + 1 < len(samples) else None
    else:
        # The whole interval lies between two samples of the scan.
        inner_lo = inner_hi = optimum
        outer_lo, outer_hi = cell_lo, cell_hi

    # A limit that is not bounded by a sign change is bounded by the twilight
    # itself, and lands exactly on the edge of the bracket.
    begin = (
        inner_lo
        if outer_lo is None
        else _pheno_margin_edge(margin, inner_lo, outer_lo, tol)
    )
    end = (
        inner_hi
        if outer_hi is None
        else _pheno_margin_edge(margin, inner_hi, outer_hi, tol)
    )
    return begin, optimum, end


def _pheno_rise_window(
    jd: float,
    body: int,
    star_name: str,
    event_type: int,
    geopos3: tuple,
    atmo4: tuple,
    obs6: tuple,
    flags: int,
) -> Tuple[float, float, float, float, float, float, float, float]:
    """Place the heliacal event of one date, at one site, in time.

    A heliacal event is an interval, not an instant: it is bounded at one end
    by the object clearing the horizon or going down, and at the other by the
    Sun making the sky too bright to show it. This locates the two horizon
    crossings that frame the requested kind of event and the part of the
    twilight between them in which the object is actually above the detection
    threshold of the sky.

    Args:
        jd: The instant the caller asked about, a Julian Day in UT.
        body: Numeric body id; not the object when ``star_name`` is given.
        star_name: Catalogue name of a fixed star, empty for anything else.
        event_type: 1 to 4. Types 1 and 4 are framed by a rise, 2 and 3 by a
            set.
        geopos3: Longitude and latitude in degrees, height above sea level in
            metres.
        atmo4: Pressure in hPa, temperature in Celsius, relative humidity in
            percent, meteorological range or total extinction.
        obs6: Age in years, Snellen ratio, binocular flag, magnification,
            aperture in millimetres, transmission.
        flags: The caller's whole flags word.

    Returns:
        Eight floats: the object's horizon crossing, the Sun's, the lag
        between them in days, the beginning, optimum and end of the visibility
        interval, its duration in days, and the Moon's best time. Every
        instant that does not exist answers the invalid-time sentinel; the lag
        is exactly zero when either crossing is missing, and the duration is
        exactly zero when the event type carries no interval for this object.
    """
    morning = _pheno_event_is_morning(event_type)
    # HELFLAG_* options share bit positions with the calculation flags, so the
    # crossings see the ephemeris selection and nothing else: a visibility
    # option must never be able to move a rise time.
    eph_flags = _heliacal_eph_flags(flags)
    is_star = bool(star_name)
    target: Union[int, str] = star_name if is_star else body

    # The crossing this date is about is the first of the requested kind that
    # is not more than four hours earlier than the instant asked about.
    search_from = jd - _PHENO_CROSSING_REACH
    object_crossing = _pheno_horizon_crossing(
        search_from, target, morning, geopos3, eph_flags
    )
    sun_crossing = _pheno_horizon_crossing(
        search_from, SUN, morning, geopos3, eph_flags
    )

    # The lag is the object's crossing minus the Sun's, so it is negative when
    # the object crosses first. Where one of them does not exist there is no
    # difference to take and the lag is exactly zero rather than a difference
    # against a sentinel.
    both_crossings = _TJD_INVALID not in (object_crossing, sun_crossing)
    lag = object_crossing - sun_crossing if both_crossings else 0.0

    best_time = (
        sun_crossing + lag * _PHENO_BEST_TIME_FRACTION
        if both_crossings and not is_star and body == MOON
        else _TJD_INVALID
    )
    absent = _TJD_INVALID

    # An apparition on the far side of a conjunction is an observable event
    # only for the Moon and the two inner planets. For every other object the
    # interval does not exist and the duration says so with an exact zero,
    # which is what tells this case apart from the other empty answers. The
    # rule holds even where the Sun does not cross the horizon at all.
    if event_type in _PHENO_FAR_SIDE_EVENTS and (
        is_star or body not in _PHENO_FAR_SIDE_OBJECTS
    ):
        return (
            object_crossing,
            sun_crossing,
            lag,
            absent,
            absent,
            absent,
            0.0,
            best_time,
        )

    # The Sun is not an object this sky can be asked to show, and without the
    # Sun's crossing there is no twilight to place a sighting in.
    if (not is_star and body == SUN) or sun_crossing == _TJD_INVALID:
        return (
            object_crossing,
            sun_crossing,
            lag,
            absent,
            absent,
            absent,
            absent,
            best_time,
        )

    if morning:
        first, last = sun_crossing - _PHENO_TWILIGHT_WIDTH, sun_crossing
    else:
        first, last = sun_crossing, sun_crossing + _PHENO_TWILIGHT_WIDTH

    objname = star_name if is_star else _PHENO_BODY_NAMES.get(body, str(body))
    seen: dict = {}

    def margin(instant: float) -> float:
        """Return the visibility margin, computing each instant only once."""
        value = seen.get(instant)
        if value is None:
            value = _pheno_visibility_margin(
                instant, objname, geopos3, atmo4, obs6, flags
            )
            seen[instant] = value
        return value

    begin, optimum, end = _pheno_visibility_interval(margin, first, last)
    if begin == _TJD_INVALID:
        return (
            object_crossing,
            sun_crossing,
            lag,
            absent,
            absent,
            absent,
            absent,
            best_time,
        )

    # The interval may not run while the object is on the wrong side of the
    # horizon, so a morning beginning is pushed to the object's rise and an
    # evening end pulled back to its set. When that crossing lies outside the
    # twilight altogether the limits come out in the wrong order and the
    # duration negative: that is the answer, not a fault to repair.
    if object_crossing != _TJD_INVALID:
        if morning:
            begin = max(begin, object_crossing)
        else:
            end = min(end, object_crossing)

    return (
        object_crossing,
        sun_crossing,
        lag,
        begin,
        optimum,
        end,
        end - begin,
        best_time,
    )


def _heliacal_pheno_ut_pythonic(
    jd: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    pressure: float = 1013.25,
    temperature: float = 15.0,
    humidity: float = 0.5,
    body: int = SUN,
    event_type: int = 1,
    flags: int = FLG_SWIEPH,
    met_range: float = 0.0,
    observer: Sequence[float] = _NAKED_OBSERVER,
) -> Tuple[Tuple[float, ...], int]:
    """Evaluate the observing geometry and visibility metrics of a heliacal event.

    The 50-slot result preserves the public compatibility layout.  Each populated
    slot is documented below using independently worded definitions; unused slots
    are returned as zero.

    Args:
        jd: Julian Day (UT) for the calculation
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        pressure: Atmospheric pressure in mbar/hPa for refraction (default 1013.25)
        temperature: Temperature in Celsius for refraction (default 15)
        humidity: Relative humidity 0.0-1.0 for atmospheric extinction (default 0.5)
        body: Planet/body ID (MERCURY, VENUS, MARS, JUPITER, SATURN, etc.)
        event_type: Type of heliacal event:
            - HELIACAL_RISING (1): Morning first visibility (heliacal rising)
            - HELIACAL_SETTING (2): Evening last visibility (heliacal setting)
            - EVENING_FIRST (3): First evening visibility (after superior conjunction)
            - MORNING_LAST (4): Last morning visibility (before superior conjunction)
        flags: Calculation flags (FLG_SWIEPH, etc.)
        met_range: Meteorological range in km, or a total extinction
            coefficient when 0 < value < 1; 0 selects the clear-sky default.
        observer: Observer tuple (age, Snellen ratio, binocular flag,
            magnification, aperture in mm, transmission); the optical values
            are read only when ``HELFLAG_OPTICAL_PARAMS`` is set.

    Returns:
        A pair ``(values, status)``. ``values`` is a tuple of 50 floats whose
        layout is fixed by the public compatibility contract. The populated
        slots are, in order:

        - 0: topocentric altitude of the object before refraction (deg).
        - 1: apparent altitude of the object after refraction (deg).
        - 2: altitude of the object from its geocentric direction (deg).
        - 3: azimuth of the object (deg, north-zero, east-positive).
        - 4, 5: topocentric altitude and azimuth of the Sun (deg).
        - 6: vertical separation of the object from the Sun using the
          topocentric altitude, the topocentric arcus visionis (deg).
        - 7: the same separation using the geocentric altitude (deg).
        - 8: signed azimuth difference, Sun minus object, in -180..180 (deg).
        - 9: arc of light, the great-circle distance between the object and
          the Sun from cos ARCL = cos ARCV * cos DAZ (deg).
        - 10: total atmospheric extinction coefficient (mag per airmass).
        - 11: smallest vertical separation at which the object would be
          detectable for its magnitude under the given sky (deg).
        - 12, 13, 14: beginning, optimum and end of the visibility window
          (Julian days, UT).
        - 15: Yallop's best-time estimate for the lunar crescent (Julian day).
        - 16: width of the lunar crescent (deg).
        - 17: Yallop q statistic of the crescent.
        - 18: Yallop visibility class as a number, 1 (A) through 6 (F).
        - 19: altitude parallax of the object, geocentric minus topocentric
          altitude (deg).
        - 20: apparent visual magnitude of the object.
        - 21, 22: rise or set instant of the object and of the Sun for the
          requested event (Julian days, UT).
        - 23: lag, slot 21 minus slot 22 when both exist, else 0.0 (days).
        - 24: duration of the visibility window (days).
        - 25: length of the lunar crescent (deg).
        - 26: elongation of the object from the Sun (deg).
        - 27: illuminated fraction of the disc, in percent.
        - 28-49: reserved, always 0.0.

        Time-valued slots that could not be determined hold the sentinel
        ``TJD_INVALID`` (99999999.0). ``status`` echoes the calculation
        flags.

    Raises:
        Error: If the body id or the event type is not accepted.

    Example:
        >>> from libephemeris import julday, _heliacal_pheno_ut_pythonic, VENUS, HELIACAL_RISING
        >>> jd = julday(2024, 1, 1, 0)
        >>> # Get heliacal phenomena for Venus at Rome
        >>> dret, flag = _heliacal_pheno_ut_pythonic(jd, 41.9, 12.5, body=VENUS,
        ...                                 event_type=HELIACAL_RISING)
        >>> print(f"Object altitude: {dret[0]:.2f}, Sun altitude: {dret[4]:.2f}")

    References:
        - Schaefer (1990), PASP 102, 212-229, DOI 10.1086/132629,
          limiting-magnitude/observer model; Schaefer (1998), Sky & Telescope
          95(5), 57-60 (VISLIMIT V-band model and its relative airmasses).
        - Krisciunas & Schaefer (1991), PASP 103, 1033-1039,
          DOI 10.1086/132921 (scattered-moonlight sky brightness).
        - Yallop (1997), NAO Technical Note 69, only for the lunar-crescent
          q statistic and class boundaries placed in result slots 17-18;
          Bruin (1977), Vistas in Astronomy 21, 331-358 (crescent width).
    """
    # --- LEB fast path ---
    from .state import get_leb_reader as _get_leb_reader

    leb_reader = _get_leb_reader()
    if leb_reader is not None:
        try:
            return _heliacal_pheno_ut_leb(
                leb_reader,
                jd,
                lat,
                lon,
                altitude,
                pressure,
                temperature,
                humidity,
                body,
                event_type,
                flags,
                met_range,
                observer,
            )
        except LEBCorruptionError:
            raise
        except (KeyError, ValueError) as leb_miss:
            from .eclipse import _raise_if_sealed_leb_miss

            # Sealed leb mode raises the documented typed error and never
            # continues past the LEB path.
            _raise_if_sealed_leb_miss(leb_miss)
            # Auto mode may continue through its normal backend chain.
            from .leb_reader import log_leb_fallback

            log_leb_fallback("heliacal", leb_miss)
    # --- END LEB fast path ---

    from .constants import (
        HELIACAL_RISING,
        HELIACAL_SETTING,
        EVENING_FIRST,
        MORNING_LAST,
    )
    from .planets import _PLANET_MAP, pheno_ut
    from .state import get_planets, get_timescale
    from skyfield.api import wgs84

    # Validate event type
    if event_type not in (
        HELIACAL_RISING,
        HELIACAL_SETTING,
        EVENING_FIRST,
        MORNING_LAST,
    ):
        raise Error(
            f"invalid event type {event_type}; use HELIACAL_RISING, "
            "HELIACAL_SETTING, EVENING_FIRST or MORNING_LAST"
        )

    # Check if this is a fixed star
    is_star = is_fixed_star(body)

    # Validate body - must be either a known planet or a fixed star
    if not is_star and body not in _PLANET_MAP:
        raise Error(f"body id {body} is neither a planet nor a fixed star")

    # Get ephemeris and timescale
    ephemeris = get_planets()
    timescale = get_timescale()

    sun = ephemeris["sun"]
    earth = ephemeris["earth"]

    # Get target - either planet or star
    planet_target = None
    skyfield_star = None
    star_magnitude = 0.0

    if is_star:
        # For fixed stars, create a Skyfield Star object
        from skyfield.api import Star
        from .fixed_stars import STAR_CATALOG

        star_magnitude = _get_star_magnitude(body)

        # Find star data from catalog
        catalog_star = None
        for catalog_entry in STAR_CATALOG:
            if catalog_entry.id == body:
                catalog_star = catalog_entry.data
                break

        if catalog_star is None:
            raise Error(f"unknown fixed star id {body}")

        # Create Skyfield Star object for position calculations
        skyfield_star = Star(
            ra_hours=catalog_star.ra_j2000 / 15.0,
            dec_degrees=catalog_star.dec_j2000,
            ra_mas_per_year=catalog_star.pm_ra * 1000.0,
            dec_mas_per_year=catalog_star.pm_dec * 1000.0,
        )
    else:
        # For planets, get the target from ephemeris
        planet_key = _PLANET_MAP[body]
        from .planets import _PLANET_FALLBACK

        try:
            planet_target = ephemeris[planet_key]
        except KeyError:
            if planet_key in _PLANET_FALLBACK:
                planet_target = ephemeris[_PLANET_FALLBACK[planet_key]]
            else:
                raise

    # Parse the observer tuple into the visibility-model parameters.
    (
        observer_age,
        snellen_ratio,
        binocular_flag,
        telescope_magnification,
        aperture_mm,
        optics_transmission,
    ) = _parse_observer_optics(observer, flags)
    # Keep the first six raw values too: the rise-window search below needs
    # them.
    observer_values = tuple(float(v) for v in tuple(observer)[:6])
    observer_tuple6 = observer_values + (0.0,) * (6 - len(observer_values))

    # Observing site as a Skyfield geographic position
    site = wgs84.latlon(lat, lon, altitude)

    # Initialize result array with 50 zeros
    values = [0.0] * 50

    # Calculate positions at the given time
    time_ut1 = timescale.ut1_jd(jd)
    site_at = earth + site

    # --- Geometric (TRUEPOS-like) places for the geometric output slots -------
    # Slots 0, 2, 3, 4, 5 and the arcs derived from them use the astrometric
    # (aberration- and light-time-free) direction, whereas an apparent place
    # would carry the ~20" annual aberration (plus, for the fast planets, the
    # light-time shift of the direction). Only slot 1 starts from the apparent
    # altitude and applies the refraction model below; the visibility detector
    # is unaffected.
    gast_hours = time_ut1.gast

    # Sun: geometric topocentric altitude/azimuth (slots 4 and 5).
    sun_geometric = (sun - site_at).at(time_ut1)
    sun_alt, sun_az, _ = sun_geometric.altaz()
    sun_alt_deg = sun_alt.degrees
    sun_az_deg = sun_az.degrees

    # Object apparent topocentric altitude, used only as the base of the
    # refracted altitude in slot 1 (see below).
    if is_star and skyfield_star is not None:
        body_apparent = site_at.at(time_ut1).observe(skyfield_star).apparent()
    else:
        body_apparent = site_at.at(time_ut1).observe(planet_target).apparent()
    body_app_alt_deg = body_apparent.altaz()[0].degrees

    # Object geometric topocentric altitude/azimuth (slots 0 and 3) and the
    # of-date geocentric RA/Dec for slot 2. A Star cannot form a geometric
    # difference vector and altaz() rejects an astrometric position, so its
    # of-date astrometric place is reduced by hand; a fixed star's diurnal
    # parallax is nil, so its topocentric and geocentric altitudes coincide.
    if is_star and skyfield_star is not None:
        topo_ra, topo_dec, _ = (
            site_at.at(time_ut1).observe(skyfield_star).radec(epoch=time_ut1)
        )
        body_alt_deg, body_az_deg = _altaz_from_radec(
            topo_ra.hours, topo_dec.degrees, lon, lat, gast_hours
        )
        geo_ra, geo_dec, _ = (
            earth.at(time_ut1).observe(skyfield_star).radec(epoch=time_ut1)
        )
        geo_ra_hours = geo_ra.hours
        geo_dec_deg = geo_dec.degrees
    else:
        body_geometric = (planet_target - site_at).at(time_ut1)
        body_alt, body_az, _ = body_geometric.altaz()
        body_alt_deg = body_alt.degrees
        body_az_deg = body_az.degrees
        geo_ra, geo_dec, _ = (planet_target - earth).at(time_ut1).radec(epoch=time_ut1)
        geo_ra_hours = geo_ra.hours
        geo_dec_deg = geo_dec.degrees

    # Geocentric altitude via hour angle. GAST is referenced to the true equinox
    # of date, so the RA/Dec must be of date too (radec(epoch=t)); mixing a
    # J2000 RA with an of-date LST would inject the J2000->date precession+
    # nutation into the hour angle.
    lst_hours = gast_hours + lon / 15.0  # Local sidereal time in hours
    hour_angle_deg = (lst_hours - geo_ra_hours) * 15.0  # Hour angle in degrees
    dec_rad = math.radians(geo_dec_deg)
    lat_rad = math.radians(lat)
    hour_angle_rad = math.radians(hour_angle_deg)
    sin_geo_alt = math.sin(lat_rad) * math.sin(dec_rad) + math.cos(lat_rad) * math.cos(
        dec_rad
    ) * math.cos(hour_angle_rad)
    sin_geo_alt = max(-1.0, min(1.0, sin_geo_alt))
    geo_alt_deg = math.degrees(math.asin(sin_geo_alt))

    # Apparent geocentric place kept for the fixed-star / fallback elongation
    # separation below (angular separation is aberration-insensitive).
    if is_star and skyfield_star is not None:
        body_geocentric_apparent = earth.at(time_ut1).observe(skyfield_star).apparent()
    else:
        body_geocentric_apparent = earth.at(time_ut1).observe(planet_target).apparent()

    # Atmospheric refraction from the APPARENT true altitude, in arcminutes:
    # R = 1.02 / tan(h + 10.3/(h + 5.11)).  Sæmundsson, Þ. (1986),
    # "Astronomical Refraction", Sky & Telescope 72, 70; as given in
    # Meeus (1998), Astronomical Algorithms 2nd ed., ch. 16, eq. 16.4. Slot 1
    # is this refracted apparent altitude.
    if body_app_alt_deg > -1:
        refraction = 1.02 / math.tan(
            math.radians(body_app_alt_deg + 10.3 / (body_app_alt_deg + 5.11))
        )
        refraction /= 60.0  # Convert arcminutes to degrees
    else:
        refraction = 0.5  # Near-horizon fallback (below the formula's range)

    # Apparent altitude (with refraction)
    app_alt_deg = body_app_alt_deg + refraction

    # Calculate arcus visionis (altitude difference between body and Sun)
    # Topocentric arcus visionis
    topo_arcus_visionis = body_alt_deg - sun_alt_deg

    # Geocentric arcus visionis
    geo_arcus_visionis = geo_alt_deg - sun_alt_deg

    # Public signed-azimuth convention: Sun minus object.
    azimuth_diff_deg = sun_az_deg - body_az_deg
    # Normalize to -180 to +180
    while azimuth_diff_deg > 180:
        azimuth_diff_deg -= 360
    while azimuth_diff_deg < -180:
        azimuth_diff_deg += 360

    # Get elongation (longitude difference) from Sun
    # For fixed stars, always calculate manually since pheno_ut doesn't support them
    if is_star:
        # Calculate elongation manually for stars
        sun_geocentric_apparent = earth.at(time_ut1).observe(sun).apparent()
        elongation = body_geocentric_apparent.separation_from(
            sun_geocentric_apparent
        ).degrees
        magnitude = star_magnitude
        phase_angle = 0.0  # Stars don't have phase angle
    else:
        try:
            body_pheno = pheno_ut(jd, body, _heliacal_eph_flags(flags))
            elongation = body_pheno[2]  # Elongation
            magnitude = body_pheno[4]  # Visual magnitude
            phase_angle = body_pheno[0]  # Phase angle (index 0, not 1)
        except (ValueError, TypeError, ArithmeticError):
            # Calculate elongation manually
            sun_geocentric_apparent = earth.at(time_ut1).observe(sun).apparent()
            elongation = body_geocentric_apparent.separation_from(
                sun_geocentric_apparent
            ).degrees
            magnitude = 0.0
            phase_angle = 0.0

    # Arc of light (slot 9) between the body and the Sun via the exact
    # right-spherical-triangle relation cos ARCL = cos ARCV * cos DAZ (ARCV
    # vertical, DAZ horizontal), not the sqrt(ARCV^2 + DAZ^2) planar form.
    arc_of_light_deg = _arc_of_light(geo_arcus_visionis, azimuth_diff_deg)

    # Schaefer model for the extinction and the required arcus visionis
    # (observer tuple parsed above).
    visibility_model = create_schaefer_model(
        pressure=pressure,
        temperature=temperature,
        humidity=humidity * 100.0 if humidity <= 1.0 else humidity,
        altitude=altitude,
        met_range=met_range,
        observer_age=observer_age,
        snellen=snellen_ratio,
        binocular=binocular_flag,
        telescope_mag=telescope_magnification,
        aperture=aperture_mm,
        transmission=optics_transmission,
        latitude=lat,
        jd=jd,
    )
    extinction_coeff = visibility_model.k_total
    required_arcus_visionis = visibility_model.arcus_visionis_required(magnitude)

    # Parallax of object (in degrees)
    # Fixed stars have essentially zero parallax
    if is_star:
        altitude_parallax_deg = 0.0
    else:
        # Altitude parallax is the geocentric-minus-topocentric altitude, not
        # the horizontal parallax angle.
        altitude_parallax_deg = geo_alt_deg - body_alt_deg

    # Rise/set of the disc centers and the visibility window use the shared
    # event-time helper and its public invalid-time sentinel.
    catalog_star_name = ""
    if is_star:
        from .fixed_stars import STAR_CATALOG as _catalog

        for catalog_entry in _catalog:
            if catalog_entry.id == body:
                catalog_star_name = catalog_entry.name
                break
    site_geopos = (lon, lat, altitude)
    atmosphere = (
        pressure,
        temperature,
        humidity * 100.0 if humidity <= 1.0 else humidity,
        met_range,
    )
    # The visibility threshold is observer-dependent, so its window slots use
    # the caller's age, Snellen ratio, and optical parameters just like slot 11.
    window_observer = observer_tuple6
    (
        object_rise_set_jd,
        sun_rise_set_jd,
        lag_days,
        visibility_start_jd,
        visibility_best_jd,
        visibility_end_jd,
        visibility_duration_days,
        yallop_best_jd,
    ) = _pheno_rise_window(
        jd,
        body,
        catalog_star_name,
        event_type,
        site_geopos,
        atmosphere,
        window_observer,
        flags,
    )

    # Illumination percentage for all bodies
    # For planets: (1 + cos(phase_angle)) / 2 * 100
    if is_star:
        # Fixed stars are unresolved point sources; the public illumination
        # field uses the full-illumination convention for them.
        illuminated_pct = 100.0
    elif phase_angle > 0:
        illuminated_pct = (1.0 + math.cos(math.radians(phase_angle))) / 2.0 * 100.0
    else:
        illuminated_pct = 0.0

    # For Moon-specific calculations
    crescent_width_deg = 0.0  # Crescent width
    crescent_length_deg = 0.0  # Crescent length

    if body == MOON:
        # Calculate Moon phase and crescent geometry
        try:
            moon_pheno = pheno_ut(jd, MOON, _heliacal_eph_flags(flags))
            illuminated_pct = moon_pheno[1] * 100.0  # [1] = illuminated fraction 0-1

            # Crescent width W = SD * (1 - cos ARCL) from Bruin, F. (1977),
            # "The first visibility of the lunar crescent", Vistas in Astronomy
            # 21, 331-358 (adopted by Yallop's q-test below), where
            # ARCL is the arc of light (Sun-Moon elongation, moon_pheno[2]) and
            # SD is the Moon's apparent semidiameter in arcmin (moon_pheno[3]
            # is the apparent diameter in degrees). Near new moon ARCL -> 0 so
            # W -> 0; near full moon ARCL -> 180 deg so W -> the disc diameter.
            #
            # Crescent length uses the actual distance-based apparent diameter
            # from pheno_ut (varies 0.49° at apogee to 0.56° at perigee):
            # L = pi * D / 2 (semicircle approximation).
            moon_diameter_deg = (
                moon_pheno[3] if len(moon_pheno) > 3 and moon_pheno[3] > 0 else 0.5
            )
            moon_arc_of_light_deg = moon_pheno[2]
            semidiameter_arcmin = moon_diameter_deg / 2.0 * 60.0
            crescent_width_arcmin = semidiameter_arcmin * (
                1.0 - math.cos(math.radians(moon_arc_of_light_deg))
            )
            crescent_width_deg = crescent_width_arcmin / 60.0  # slot 16, degrees
            crescent_length_deg = math.pi * moon_diameter_deg / 2
        except (ValueError, TypeError, ArithmeticError):
            pass

    # Yallop q-test for lunar crescent visibility, Yallop, B.D. (1997), NAO
    # Technical Note No. 69:
    # q = (ARCV - (11.8371 - 6.3226*W + 0.7319*W^2 - 0.1018*W^3)) / 10
    if body == MOON:
        width_arcmin = crescent_width_deg * 60.0  # Convert to arcminutes for formula
        yallop_threshold = (
            11.8371
            - 6.3226 * width_arcmin
            + 0.7319 * width_arcmin**2
            - 0.1018 * width_arcmin**3
        )
        yallop_q = (geo_arcus_visionis - yallop_threshold) / 10.0
        # Slot 18 carries the Yallop visibility class code (1..6), not the raw
        # arcus-visionis criterion polynomial (which stays internal above).
        yallop_class_code = _yallop_visibility_code(yallop_q)
    else:
        yallop_q = 0.0
        yallop_class_code = 0.0

    # Fill in the result array
    values[0] = body_alt_deg  # topocentric altitude, unrefracted
    values[1] = app_alt_deg  # apparent altitude, refracted
    values[2] = geo_alt_deg  # altitude from the geocentric direction
    values[3] = body_az_deg  # azimuth of the object
    values[4] = sun_alt_deg  # topocentric altitude of the Sun
    values[5] = sun_az_deg  # azimuth of the Sun
    values[6] = topo_arcus_visionis  # object minus Sun, topocentric altitudes
    values[7] = geo_arcus_visionis  # object minus Sun, geocentric altitude
    values[8] = azimuth_diff_deg  # Sun minus object azimuth, signed
    values[9] = arc_of_light_deg  # great-circle distance object-Sun
    values[10] = extinction_coeff  # total extinction, mag per airmass
    values[11] = required_arcus_visionis  # detection threshold on slot 6
    values[12] = visibility_start_jd  # visibility window begins
    values[13] = visibility_best_jd  # visibility window optimum
    values[14] = visibility_end_jd  # visibility window ends
    values[15] = yallop_best_jd  # Yallop best time (Moon)
    values[16] = crescent_width_deg  # lunar crescent width
    values[17] = yallop_q  # Yallop q statistic (Moon)
    values[18] = yallop_class_code  # Yallop class 1..6 (Moon)
    values[19] = altitude_parallax_deg  # geocentric minus topocentric altitude
    values[20] = magnitude  # apparent visual magnitude
    values[21] = object_rise_set_jd  # object rise/set for this event
    values[22] = sun_rise_set_jd  # Sun rise/set for this event
    values[23] = lag_days  # object minus Sun rise/set lag
    values[24] = visibility_duration_days  # visibility window length
    values[25] = crescent_length_deg  # lunar crescent length
    values[26] = elongation  # elongation from the Sun
    values[27] = illuminated_pct  # illuminated fraction, percent
    # slots 28-49 are reserved and stay 0.0

    # Coerce to native Python floats: several entries (alt/az from the Skyfield
    # backend) are numpy.float64, and the public contract returns native floats.
    return tuple(float(v) for v in values), flags


def heliacal_pheno_ut(
    tjdut: float,
    geopos: tuple,
    atmo: tuple,
    observer: tuple,
    objname: str,
    eventtype: int,
    flags: int = FLG_SWIEPH,
) -> Tuple[float, ...]:
    """
    Provides data relevant for the calculation of heliacal risings and settings.

    This public wrapper accepts the compatibility parameter layout and returns
    a flat 50-element tuple.

    Args:
        tjdut: Julian Day (UT) for the calculation
        geopos: Geographic position (lon, lat, alt_m)
        atmo: Atmospheric conditions (pressure, temperature, humidity%, met_range)
        observer: Observer description (age, snellen, binocular, mag, aperture, transmission)
        objname: Name of planet or fixed star (e.g., "Venus", "Sirius")
        eventtype: Type of heliacal event (1-4)
        flags: Calculation flags

    Returns:
        Flat tuple of 50 floats with heliacal phenomena data; the slot layout
        is documented on ``_heliacal_pheno_ut_pythonic``.
    """
    # Parse geopos
    lon = geopos[0] if len(geopos) > 0 else 0.0
    lat = geopos[1] if len(geopos) > 1 else 0.0
    altitude = geopos[2] if len(geopos) > 2 else 0.0

    # Parse atmo with defaults
    pressure = atmo[0] if len(atmo) > 0 and atmo[0] > 0 else 1013.25
    temperature = atmo[1] if len(atmo) > 1 else 15.0
    humidity_pct = atmo[2] if len(atmo) > 2 else 40.0
    # atmo[3]: meteorological range (km) or ktot (0<v<1); drives extinction.
    met_range = atmo[3] if len(atmo) > 3 else 0.0

    # datm[2] is relative humidity in percent in the public input tuple.
    humidity = humidity_pct / 100.0

    # Parse objname to body ID. The public phenomena result accepts the Moon
    # and Sun; heliacal_ut/vis_limit_mag apply narrower physical domains.
    body_id = _parse_object_name(objname, allow_moon=True, allow_sun=True)

    # Call internal function
    dret, retflag = _heliacal_pheno_ut_pythonic(
        jd=tjdut,
        lat=lat,
        lon=lon,
        altitude=altitude,
        pressure=pressure,
        temperature=temperature,
        humidity=humidity,
        body=body_id,
        event_type=eventtype,
        flags=flags,
        met_range=met_range,
        observer=observer,
    )

    # Return the public flat 50-tuple.
    return dret


def vis_limit_mag(
    tjdut: float,
    geopos: tuple,
    atmo: tuple,
    observer: tuple,
    objname: str,
    flags: int = FLG_SWIEPH,
) -> Tuple[float, Tuple[float, ...]]:
    """Evaluate the limiting visual magnitude for an object.

    Args:
        tjdut: Observation time as a Julian Day in UT.
        geopos: Longitude, latitude, and observer altitude.
        atmo: Pressure, temperature, humidity, and meteorological-range or
            total-extinction input.
        observer: Observer age and acuity, with optional optical-instrument
            parameters.
        objname: Planet, fixed-star name, or supported numeric body string.
        flags: Ephemeris and heliacal visibility flags.

    Returns:
        A status code and an eight-value detail tuple containing limiting
        magnitude, object/Sun/Moon horizontal coordinates, and object
        magnitude.

    Raises:
        ValueError: If the object name is empty or invalid.

    Notes:
        The Schaefer visibility model includes solar and lunar sky brightness,
        atmospheric extinction, eye adaptation, and optional optics.
    """
    # --- LEB fast path ---
    from .state import get_leb_reader as _get_leb_reader

    _leb_rdr = _get_leb_reader()
    if _leb_rdr is not None:
        try:
            return _vis_limit_mag_leb(
                _leb_rdr,
                tjdut,
                geopos,
                atmo,
                observer,
                objname,
                flags,
            )
        except LEBCorruptionError:
            raise
        except (KeyError, ValueError) as _leb_err:
            from .eclipse import _raise_if_sealed_leb_miss

            # Sealed leb mode raises the documented typed error and never
            # continues past the LEB path.
            _raise_if_sealed_leb_miss(_leb_err)
            # Auto mode may continue through its normal backend chain.
            from .leb_reader import log_leb_fallback

            log_leb_fallback("heliacal", _leb_err)
    # --- END LEB fast path ---

    from .planets import _PLANET_MAP, pheno_ut
    from .fixed_stars import fixstar2_ut, fixstar2_mag
    from .state import get_planets, get_timescale
    from .constants import (
        SUN,
        MOON,
        HELFLAG_VISLIM_DARK,
        HELFLAG_VISLIM_NOMOON,
        HELFLAG_BELOW_HORIZON,
    )
    from skyfield.api import wgs84

    if not objname:
        raise Error("the object name is empty")

    # Parse geographic position
    lon = geopos[0] if len(geopos) > 0 else 0.0
    lat = geopos[1] if len(geopos) > 1 else 0.0
    alt_m = geopos[2] if len(geopos) > 2 else 0.0

    # Parse atmospheric conditions
    pressure = atmo[0] if len(atmo) > 0 else 1013.25
    temperature = atmo[1] if len(atmo) > 1 else 15.0
    humidity_pct = atmo[2] if len(atmo) > 2 else 50.0
    met_range = atmo[3] if len(atmo) > 3 else 0.0

    # Parse observer data (optical params only when HELFLAG_OPTICAL_PARAMS set).
    (
        observer_age,
        snellen_ratio,
        obs_binocular,
        obs_telescope_mag,
        obs_aperture,
        obs_transmission,
    ) = _parse_observer_optics(observer, flags)

    # Get ephemeris and timescale
    eph = get_planets()
    ts = get_timescale()

    # Create observer location
    obs_location = wgs84.latlon(lat, lon, alt_m)
    earth = eph["earth"]
    sun = eph["sun"]
    moon = eph["moon"]

    t = ts.ut1_jd(tjdut)
    observer_at = earth + obs_location

    # Calculate Sun position
    sun_app = observer_at.at(t).observe(sun).apparent()
    sun_alt_deg, sun_az_deg, _ = sun_app.altaz()
    sun_alt = sun_alt_deg.degrees
    sun_az = sun_az_deg.degrees

    # Calculate Moon position
    moon_app = observer_at.at(t).observe(moon).apparent()
    moon_alt_deg, moon_az_deg, _ = moon_app.altaz()
    moon_alt = moon_alt_deg.degrees
    moon_az = moon_az_deg.degrees

    # Determine if objname is a planet ID or name
    body_id = None
    is_fixed_star = False

    # Try parsing as integer (planet ID)
    try:
        body_id = int(objname)
    except ValueError:
        # Try to find planet by name
        name_upper = objname.upper().strip()
        planet_names = {
            "SUN": SUN,
            "MOON": MOON,
            "MERCURY": 2,
            "VENUS": 3,
            "MARS": 4,
            "JUPITER": 5,
            "SATURN": 6,
            "URANUS": 7,
            "NEPTUNE": 8,
            "PLUTO": 9,
        }
        if name_upper in planet_names:
            body_id = planet_names[name_upper]
        else:
            # Assume it's a fixed star
            is_fixed_star = True

    # The Sun has no meaningful limiting magnitude in this model.
    if not is_fixed_star and body_id == SUN:
        raise Error("vis_limit_mag() is not defined for the Sun")

    # Calculate object position and magnitude
    obj_alt = 0.0
    obj_az = 0.0
    obj_mag = 0.0

    if is_fixed_star:
        # Fixed star calculation
        try:
            star_result, _star_name_out, _retflag = fixstar2_ut(
                objname, tjdut, _heliacal_eph_flags(flags)
            )

            # star_result is (lon, lat, dist, lon_speed, lat_speed, dist_speed)
            # We need to convert ecliptic to horizontal
            star_lon = star_result[0]
            star_lat = star_result[1]

            # Get star magnitude
            try:
                star_mag_val, _star_name_mag = fixstar2_mag(objname)
                obj_mag = star_mag_val
            except (ValueError, TypeError, ArithmeticError):
                obj_mag = 2.0  # Default magnitude if not found

            # Convert ecliptic to equatorial then to horizontal
            # Simplified: use azalt function if available
            from .utils import azalt, ECL2HOR

            hor_result = azalt(
                tjdut,
                ECL2HOR,
                (lon, lat, alt_m),
                pressure,
                temperature,
                (star_lon, star_lat, 1.0),
            )
            obj_az = (hor_result[0] + 180.0) % 360.0
            obj_alt = hor_result[1]

        except ValueError:
            raise
        except (KeyError, IndexError, OSError) as e:
            # Star not found or other error
            raise Error(
                f"no fixed star matches the name {objname.lower()!r}: {e}"
            ) from e
    else:
        # Planet calculation
        if body_id is None:
            raise Error(f"unknown object name {objname!r}")

        # Get planet name from _PLANET_MAP
        if body_id in _PLANET_MAP:
            target_name = _PLANET_MAP[body_id]
            # Try planet center first, fall back to barycenter if not available
            from .planets import _PLANET_FALLBACK

            try:
                target = eph[target_name]
            except KeyError:
                if target_name in _PLANET_FALLBACK:
                    target = eph[_PLANET_FALLBACK[target_name]]
                else:
                    raise

            # Calculate position
            body_app = observer_at.at(t).observe(target).apparent()
            body_alt_deg, body_az_deg, _ = body_app.altaz()
            obj_alt = body_alt_deg.degrees
            obj_az = body_az_deg.degrees

            # Get magnitude from pheno
            try:
                pheno_result = pheno_ut(tjdut, body_id, _heliacal_eph_flags(flags))
                obj_mag = pheno_result[4]  # Visual magnitude
            except (ValueError, TypeError, ArithmeticError):
                obj_mag = 0.0  # Default bright
        else:
            raise Error(f"body id {body_id} is neither a planet nor a fixed star")

    # Check if object is below horizon
    if obj_alt < 0:
        # Public result contract: the below-horizon flag uses a sentinel in
        # the limiting-magnitude slot and zeros in the remaining slots.
        dret = (-100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        return float(HELFLAG_BELOW_HORIZON), dret

    # Apply HELFLAG options
    use_dark_sky = bool(flags & HELFLAG_VISLIM_DARK)
    exclude_moon = bool(flags & HELFLAG_VISLIM_NOMOON)

    if use_dark_sky:
        sun_alt = -90.0  # Assume Sun at nadir

    if exclude_moon:
        moon_alt = -90.0  # Assume Moon at nadir

    if body_id == MOON:
        # Compatibility contract: when the observed object IS the Moon,
        # its own light is not sky glow hindering the observation, so the
        # moonlight term is removed exactly as with VISLIM_NOMOON (the sky
        # stays dark and the limiting magnitude reflects the object alone).
        moon_alt = -90.0

    # Create Schaefer model for visibility calculations
    schaefer = create_schaefer_model(
        pressure=pressure,
        temperature=temperature,
        humidity=humidity_pct,
        met_range=met_range,
        altitude=alt_m,
        observer_age=observer_age,
        snellen=snellen_ratio,
        binocular=obs_binocular,
        telescope_mag=obs_telescope_mag,
        aperture=obs_aperture,
        transmission=obs_transmission,
        latitude=lat,
        jd=tjdut,
        pressure_scales_extinction=True,
    )

    # Calculate Moon phase for sky brightness
    moon_phase = 0.0
    moon_obj_angle = 180.0
    sun_obj_angle = 180.0
    try:
        moon_pheno = pheno_ut(tjdut, MOON, _heliacal_eph_flags(flags))
        # pheno_ut[1] is the illuminated fraction (0 = new, 1 = full), exactly
        # what SchaeferModel.sky_brightness_moon expects. Deriving it from the
        # phase angle (pheno_ut[0], 0 deg = full) as (1 - cos)/2 inverts the
        # convention — full moon would be passed as new and vice versa.
        moon_phase = moon_pheno[1]

        # Calculate angular separation between object and Moon
        if is_fixed_star:
            # Real separations from the star's ecliptic position vs the
            # topocentric Moon and Sun (the old fixed 90-degree
            # placeholders distorted the sky-brightness model).
            from skyfield.framelib import ecliptic_frame

            from .utils import angular_separation

            m_lat, m_lon, _md = moon_app.frame_latlon(ecliptic_frame)
            s_lat, s_lon, _sd = sun_app.frame_latlon(ecliptic_frame)
            moon_obj_angle = angular_separation(
                star_lon, star_lat, m_lon.degrees % 360.0, m_lat.degrees
            )
            sun_obj_angle = angular_separation(
                star_lon, star_lat, s_lon.degrees % 360.0, s_lat.degrees
            )
        else:
            body_app_geo = observer_at.at(t).observe(target).apparent()
            moon_obj_angle = body_app_geo.separation_from(moon_app).degrees
            sun_obj_angle = body_app_geo.separation_from(sun_app).degrees
    except (ValueError, TypeError, AttributeError):
        pass

    # Calculate limiting magnitude using Schaefer model
    limiting_mag = schaefer.limiting_magnitude(
        sun_alt=sun_alt,
        moon_alt=moon_alt if not exclude_moon else -90.0,
        moon_phase=moon_phase,
        obj_alt=obj_alt,
        sun_obj_angle=sun_obj_angle,
        moon_obj_angle=moon_obj_angle,
    )

    vision_type = _vislim_scotopic_flag(
        schaefer,
        sun_alt if not use_dark_sky else -90.0,
        moon_alt if not exclude_moon else -90.0,
        moon_phase,
        obj_alt,
        sun_obj_angle,
        moon_obj_angle,
        flags,
    )

    # Build the public 10-element result; dret[7] is the body's magnitude
    # without extinction. The alt/az values come from the Skyfield backend as
    # numpy.float64; coerce to native Python floats.
    dret = (
        float(limiting_mag),  # 0: Limiting visual magnitude
        float(obj_alt),  # 1: Altitude of object
        float(obj_az),  # 2: Azimuth of object
        float(sun_alt),  # 3: Altitude of Sun
        float(sun_az),  # 4: Azimuth of Sun
        float(moon_alt),  # 5: Altitude of Moon
        float(moon_az),  # 6: Azimuth of Moon
        float(obj_mag),  # 7: Magnitude of object
        0.0,  # 8: Reserved
        0.0,  # 9: Reserved
    )

    return float(vision_type), dret


# Public compatibility alias.
vis_limit_mag = vis_limit_mag
