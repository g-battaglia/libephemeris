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
      212-229, DOI 10.1086/132629; and Schaefer (1993), "Astronomy and the
      Limits of Vision", Vistas in Astronomy 36, 311-361
    - Yallop (1997), NAO Technical Note 69, lunar-crescent q criterion
    - Kasten & Young (1989), "Revised optical air mass tables and
      approximation formula", DOI 10.1364/AO.28.004735

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
from typing import Tuple, NamedTuple, Optional, Sequence

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
    arcus-visionis criterion polynomial) in ``heliacal_pheno_ut`` element 18
    (qCrit).

    Args:
        q: Yallop q-test value (``heliacal_pheno_ut`` element 17, qYal).

    Returns:
        Visibility class code 1..6, returned as a float to match the dret slot.
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
    Constants for the Schaefer atmospheric model.

    Based on Schaefer (1990).
    """

    # Rayleigh scattering coefficient at sea level (per airmass)
    # At 550nm (V-band), k_r = 0.1451 at sea level
    K_RAYLEIGH_SEA_LEVEL = 0.1451

    # Ozone absorption coefficient (per airmass)
    # Small contribution at visual wavelengths
    K_OZONE = 0.016

    # Aerosol scattering base coefficient
    # Typical value 0.05-0.15 depending on conditions
    K_AEROSOL_BASE = 0.08

    # Scale height for pressure (km)
    SCALE_HEIGHT_PRESSURE = 8.5

    # Scale height for aerosols (km)
    SCALE_HEIGHT_AEROSOL = 1.5

    # Airglow brightness (nanoLamberts) - natural sky glow
    AIRGLOW = 145.0

    # Zodiacal light brightness (nanoLamberts) - typical value
    ZODIACAL_LIGHT = 100.0

    # Dark sky brightness at zenith (mag/arcsec^2)
    DARK_SKY_MAG = 21.6

    # Full Moon brightness relative to Sun (magnitude difference)
    FULL_MOON_MAG_DIFF = 14.0

    # Visual limiting magnitude for perfect dark sky conditions
    PERFECT_SKY_LIM_MAG = 6.5

    # Legacy descriptive defaults retained for callers inspecting this class.
    # The active arcus-visionis calculation below root-solves the Schaefer
    # limiting-magnitude model and does not read this mapping. These exact
    # object values are project defaults, not a transcription from Ptolemy.
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

    # Threshold contrast for naked-eye visibility (log units)
    THRESHOLD_CONTRAST = 0.0094  # Schaefer's value

    # Eye pupil diameter in dark-adapted conditions (mm)
    DARK_PUPIL_DIAMETER = 7.0

    # Eye pupil diameter in bright conditions (mm)
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
# photometric band (0.55 um), with the atmospheric scale heights refined
# per V. Reijs (archaeocosmology.org / ARCHAEOCOSMO): Rayleigh 8515 m,
# aerosol 3745 m, water 3000 m. The extinction, airmass, twilight,
# moonlight and dark-night-sky components and the contrast-threshold
# (Hecht) magnitude conversion are exactly Schaefer's; the dark-night-sky
# normalisation uses a published natural-sky luminance instead of a value
# inferred from another ephemeris implementation.
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
_VL_MODEL_RADIANCE_PER_NL = 1.11e-15
# Reijs refined scale heights (m)
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
) -> Tuple[float, float, float, float]:
    """Return the four V-band extinction coefficients (KR, KA, KO, KW).

    Schaefer VISLIMIT extinction with Reijs scale heights and the seasonal
    ozone term. Coefficients are per unit airmass.  Aerosol extinction varies
    with humidity and altitude; no empirically recovered seasonal multiplier
    is applied.
    """
    Mc = _vl_month_cont(jd)
    RA = (Mc - 3.0) * 30.0 * _RD
    LT = latitude * _RD
    RH = min(max(humidity_pct, 1.0), 99.5)
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
    """Schaefer sky airmass (gas) as a function of altitude."""
    if alt_deg <= 0.0:
        return 40.0
    c = math.cos((90.0 - alt_deg) * _RD)
    return 1.0 / (c + 0.025 * math.exp(-11.0 * c))


def _vl_airmass_gas(alt_deg: float) -> float:
    """Return Schaefer VISLIMIT's gas/Rayleigh relative airmass."""
    c = math.cos((90.0 - alt_deg) * _RD)
    return 1.0 / (c + 0.0286 * math.exp(-10.5 * c))


def _vl_airmass_aer(alt_deg: float) -> float:
    """Return Schaefer VISLIMIT's aerosol relative airmass."""
    c = math.cos((90.0 - alt_deg) * _RD)
    return 1.0 / (c + 0.0123 * math.exp(-24.5 * c))


def _vl_airmass_ozone(alt_deg: float) -> float:
    """Return the spherical-shell ozone airmass used by VISLIMIT."""
    s = math.sin((90.0 - alt_deg) * _RD)
    val = 1.0 - (s / (1.0 + 20.0 / 6378.0)) ** 2
    return 1.0 / math.sqrt(val) if val > 0 else 40.0


def _vl_hecht_threshold(bl: float) -> float:
    """Hecht/Knoll-Schaefer contrast threshold (foot-candles) for sky ``bl``.

    Two photoreceptor regimes (photopic above / scotopic below 1500 nL), exactly
    as in Schaefer's VISLIMIT: ``TH = C1 * (1 + sqrt(C2 * BL))**2``.
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
            # Meteorological visual range (km): Koschmieder aerosol.
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

        # Dark night sky
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
            # Twilight brightness (Sun below the horizon).
            bt = 10.0 ** (
                -0.4 * (_VL_MS - _VL_MOI + 32.5 - sun_alt - (Z / (360.0 * K)))
            )
            b += bt * (100.0 / rs) * one_minus

        # Moonlight
        if moon_alt > 0.0:
            rm = max(moon_obj_angle, 1.0)
            am = math.degrees(math.acos(min(1.0, max(-1.0, 2.0 * moon_phase - 1.0))))
            mm = -12.73 + 0.026 * abs(am) + 4e-9 * (am**4)
            xm = _vl_airmass(moon_alt)
            c3 = 10.0 ** (-0.4 * K * xm)
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
    return SchaeferModel(atmo, observer, latitude=latitude, jd=jd)


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
        star_id: Star ID (SE_* constant, e.g., SIRIUS)

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
    raise Error(f"Star ID {star_id} not found in catalog")


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
        body_id: SE_* body constant.
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
            f"Invalid event_type: {event_type}. Use HELIACAL_RISING, "
            "HELIACAL_SETTING, EVENING_FIRST, or MORNING_LAST."
        )
    if event_type in (EVENING_FIRST, MORNING_LAST):
        if is_fixed_star(body):
            raise Error(
                "EVENING_FIRST and MORNING_LAST are not valid for fixed stars. "
                "For fixed stars, use HELIACAL_RISING or HELIACAL_SETTING."
            )
        if not is_inner_planet(body) and body != MOON:
            raise Error(
                "EVENING_FIRST and MORNING_LAST are only valid for inner planets "
                "(Mercury, Venus) and the Moon. For outer planets, use "
                "HELIACAL_RISING or HELIACAL_SETTING."
            )
    if body == SUN:
        raise Error("SUN is not valid for heliacal calculations")
    if body == MOON and event_type in (HELIACAL_RISING, HELIACAL_SETTING):
        raise Error("the Moon has no heliacal rising or setting (event types 1, 2)")

    is_star = is_fixed_star(body)
    if not is_star and body not in _PLANET_MAP:
        raise Error(f"illegal planet number {body}.")

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
            if visible:
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
    """LEB-backed implementation of _heliacal_pheno_ut_pythonic()."""
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
            f"Invalid event_type: {event_type}. Use HELIACAL_RISING, "
            "HELIACAL_SETTING, EVENING_FIRST, or MORNING_LAST."
        )

    is_star = is_fixed_star(body)
    if not is_star and body not in _PLANET_MAP:
        raise Error(f"illegal planet number {body}.")

    star_name: str | None = None
    star_magnitude = 0.0
    if is_star:
        star_name = _star_name_from_id(body)
        star_magnitude = _get_star_magnitude(body)

    geopos = (lon, lat, altitude)

    # --- positions via LEB ---
    # Sun
    sun_az, sun_alt_deg, _ = _leb_body_altaz(
        reader,
        jd,
        SUN,
        geopos,
        pressure,
        temperature,
    )

    # Body
    body_az_deg, body_alt_deg, _body_app_alt = _leb_body_altaz(
        reader,
        jd,
        body,
        geopos,
        pressure,
        temperature,
        is_star=is_star,
        star_name=star_name,
    )

    # Atmospheric refraction from true altitude, in arcminutes:
    # R = 1.02 / tan(h + 10.3/(h + 5.11)).  Sæmundsson, Þ. (1986),
    # "Astronomical Refraction", Sky & Telescope 72, 70; as given in
    # Meeus (1998), Astronomical Algorithms 2nd ed., ch. 16, eq. 16.4.
    # Below h = -1° the 0.5° fallback approximates horizon refraction.
    if body_alt_deg > -1:
        refraction = 1.02 / math.tan(
            math.radians(body_alt_deg + 10.3 / (body_alt_deg + 5.11))
        )
        refraction /= 60.0
    else:
        refraction = 0.5
    app_alt_deg = body_alt_deg + refraction

    # Get ecliptic positions for parallax and elongation calculations
    body_ecl = _leb_ecliptic_pos(
        reader,
        jd,
        body,
        geopos,
        is_star=is_star,
        star_name=star_name,
    )
    sun_ecl = _leb_ecliptic_pos(reader, jd, SUN, geopos)

    # Geocentric altitude from equatorial coords (matching Skyfield path).
    # GAST is referenced to the true equinox of date, so the RA/Dec must be of
    # date (FLG_EQUATORIAL alone, NOT FLG_J2000) to keep the hour-angle frame
    # consistent. Mixing J2000 RA with of-date GAST injected the full
    # J2000→date precession+nutation into the hour angle.
    from .constants import FLG_EQUATORIAL

    if is_star:
        from .fixed_stars import fixstar_ut

        assert star_name is not None
        _eq_pos, _, _ = fixstar_ut(star_name, jd, FLG_EQUATORIAL | FLG_SPEED)
        _body_ra_deg, _body_dec_deg = _eq_pos[0], _eq_pos[1]
    else:
        from .planets import calc_ut as _scu_hp

        _eq_pos, _ = _scu_hp(jd, body, FLG_EQUATORIAL | FLG_SPEED)
        _body_ra_deg, _body_dec_deg = _eq_pos[0], _eq_pos[1]

    # Use GAST to match the Skyfield path's of-date RA + t.gast
    from .state import get_timescale as _gts_hp

    _ts_hp = _gts_hp()
    _t_hp = _ts_hp.ut1_jd(jd)
    _gast_hours = _t_hp.gast
    _lst_hours = _gast_hours + geopos[0] / 15.0
    _ha_deg = (_lst_hours * 15.0 - _body_ra_deg) % 360.0
    _ha_rad = math.radians(_ha_deg)
    _dec_rad = math.radians(_body_dec_deg)
    _lat_rad = math.radians(geopos[1])
    _sin_alt = math.sin(_lat_rad) * math.sin(_dec_rad) + math.cos(_lat_rad) * math.cos(
        _dec_rad
    ) * math.cos(_ha_rad)
    geo_alt_deg = math.degrees(math.asin(max(-1.0, min(1.0, _sin_alt))))

    # Arcus visionis
    tav_act = body_alt_deg - sun_alt_deg
    arcv_act = geo_alt_deg - sun_alt_deg

    # Public signed-azimuth convention: Sun minus object.
    daz_act = sun_az - body_az_deg
    while daz_act > 180:
        daz_act -= 360
    while daz_act < -180:
        daz_act += 360

    # Elongation and magnitude
    if is_star:
        elongation = angular_separation(
            sun_ecl[0],
            sun_ecl[1],
            body_ecl[0],
            body_ecl[1],
        )
        magnitude = star_magnitude
        phase_angle = 0.0
    else:
        try:
            pheno = pheno_ut(jd, body, _heliacal_eph_flags(flags))
            elongation = pheno[2]
            magnitude = pheno[4]
            phase_angle = pheno[0]
        except (ValueError, TypeError, ArithmeticError):
            elongation = angular_separation(
                sun_ecl[0],
                sun_ecl[1],
                body_ecl[0],
                body_ecl[1],
            )
            magnitude = 0.0
            phase_angle = 0.0

    arcl_act = math.sqrt(arcv_act**2 + daz_act**2)

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
        jd=jd,
    )
    k_act = schaefer.k_total
    min_tav = schaefer.arcus_visionis_required(magnitude)

    # Parallax
    if is_star:
        parallax = 0.0
    else:
        # Altitude parallax is the geocentric-minus-topocentric altitude, not
        # the horizontal parallax angle.
        parallax = geo_alt_deg - body_alt_deg

    # Rise/set of the disc centers and the visibility window via the
    # shared event-time helper. Missing instants use the public sentinel.
    _star_name_p = ""
    if is_star:
        from .fixed_stars import STAR_CATALOG as _SC

        for _entry in _SC:
            if _entry.id == body:
                _star_name_p = _entry.name
                break
    (
        rise_o,
        rise_s,
        lag,
        t_first_vr,
        t_best_vr,
        t_last_vr,
        tvis_vr,
        t_b_yallop,
    ) = _pheno_rise_window(
        jd,
        body,
        _star_name_p,
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
        illumination = 100.0
    elif phase_angle > 0:
        illumination = (1.0 + math.cos(math.radians(phase_angle))) / 2.0 * 100.0
    else:
        illumination = 0.0

    w_moon = 0.0
    l_moon = 0.0
    if body == MOON:
        try:
            moon_pheno = pheno_ut(jd, MOON, _heliacal_eph_flags(flags))
            illumination = moon_pheno[1] * 100.0  # [1] is illuminated fraction 0-1
            moon_diameter = (
                moon_pheno[3] if len(moon_pheno) > 3 and moon_pheno[3] > 0 else 0.5
            )
            # Crescent width (Yallop/Bruin): W = SD * (1 - cos ARCL), where
            # ARCL is the arc of light (Sun-Moon elongation, moon_pheno[2]) and
            # SD is the Moon's apparent semidiameter in arcmin (moon_pheno[3]
            # is the apparent diameter in degrees). Near new moon ARCL -> 0 so
            # W -> 0; near full moon ARCL -> 180 deg so W -> the disc diameter.
            arcl_deg = moon_pheno[2]
            sd_arcmin = moon_diameter / 2.0 * 60.0
            w_arcmin = sd_arcmin * (1.0 - math.cos(math.radians(arcl_deg)))
            w_moon = w_arcmin / 60.0  # stored in degrees (dret[16])
            l_moon = math.pi * moon_diameter / 2
        except (ValueError, TypeError, ArithmeticError):
            pass

    if body == MOON:
        w = w_moon * 60.0
        q_criterion = 11.8371 - 6.3226 * w + 0.7319 * w**2 - 0.1018 * w**3
        q_yallop = (arcv_act - q_criterion) / 10.0
        # dret[18] (qCrit) is the Yallop visibility CLASS code (1..6), not the
        # raw arcus-visionis criterion polynomial (which stays internal above).
        q_crit = _yallop_visibility_code(q_yallop)
    else:
        q_yallop = 0.0
        q_crit = 0.0

    dret = [0.0] * 50
    dret[0] = body_alt_deg
    dret[1] = app_alt_deg
    dret[2] = geo_alt_deg
    dret[3] = body_az_deg
    dret[4] = sun_alt_deg
    dret[5] = sun_az
    dret[6] = tav_act
    dret[7] = arcv_act
    dret[8] = daz_act
    dret[9] = arcl_act
    dret[10] = k_act
    dret[11] = min_tav
    dret[12] = t_first_vr
    dret[13] = t_best_vr
    dret[14] = t_last_vr
    dret[15] = t_b_yallop
    dret[16] = w_moon
    dret[17] = q_yallop
    dret[18] = q_crit
    dret[19] = parallax
    dret[20] = magnitude
    dret[21] = rise_o
    dret[22] = rise_s
    dret[23] = lag
    dret[24] = tvis_vr
    dret[25] = l_moon
    dret[26] = elongation
    dret[27] = illumination

    return tuple(dret), flags


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
        raise Error("objname cannot be empty")

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
        raise Error("it makes no sense to call vis_limit_mag() for the Sun")

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
            raise Error(f"could not find star name {objname.lower()}: {e}") from e
    else:
        if body_id is None:
            raise Error(f"Unknown object: {objname}")

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
            raise Error(f"illegal planet number {body_id}.")

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
        - Schaefer (1993), Vistas in Astronomy 36, 311-361, visibility
          context and observer effects.
        - Kasten & Young (1989), DOI 10.1364/AO.28.004735, optical airmass.

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
            from .state import get_calc_mode

            if get_calc_mode() == "leb":
                raise
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
            f"Invalid event_type: {event_type}. Use HELIACAL_RISING, "
            "HELIACAL_SETTING, EVENING_FIRST, or MORNING_LAST."
        )

    # EVENING_FIRST and MORNING_LAST are only valid for inner planets
    # (Mercury and Venus) because they relate to superior conjunction visibility.
    # Outer planets only have one type of conjunction and these events don't apply.
    # Fixed stars only have heliacal rising and setting.
    if event_type in (EVENING_FIRST, MORNING_LAST):
        if is_fixed_star(body):
            raise Error(
                "EVENING_FIRST and MORNING_LAST are not valid for fixed stars. "
                "For fixed stars, use HELIACAL_RISING or HELIACAL_SETTING."
            )
        if not is_inner_planet(body) and body != MOON:
            raise Error(
                "EVENING_FIRST and MORNING_LAST are only valid for inner planets "
                "(Mercury, Venus) and the Moon. For outer planets, use "
                "HELIACAL_RISING or HELIACAL_SETTING."
            )

    # Sun and Moon are not valid for heliacal events
    if body == SUN:
        raise Error("SUN is not valid for heliacal calculations")
    if body == MOON and event_type in (HELIACAL_RISING, HELIACAL_SETTING):
        raise Error("the Moon has no heliacal rising or setting (event types 1, 2)")

    # Check if this is a fixed star
    is_star = is_fixed_star(body)

    # Validate body - must be either a known planet or a fixed star
    if not is_star and body not in _PLANET_MAP:
        raise Error(f"illegal planet number {body}.")

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
            raise Error(f"Star ID {body} not found in catalog")

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
            if visible:
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
            f"Invalid event_type: {eventtype}. Use HELIACAL_RISING, "
            "HELIACAL_SETTING, EVENING_FIRST, or MORNING_LAST."
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
        ev_name = {
            HELIACAL_RISING: "morning first",
            HELIACAL_SETTING: "evening last",
        }[eventtype]
        raise Error(f"{ev_name} (event type {eventtype}) does not exist for the Moon")

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
        Body ID (SE_* constant) for planets

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
            raise Error("Sun is not valid for heliacal calculations")
        if parsed_body_id == MOON and not allow_moon:
            raise Error("Moon is not valid for heliacal calculations")
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
            raise Error("Sun is not valid for heliacal calculations")
        if body_id == MOON and not allow_moon:
            raise Error("Moon is not valid for heliacal calculations")

        return body_id

    # Try to resolve as a fixed star name
    from .fixed_stars import resolve_star_name

    star_id = resolve_star_name(object_name)
    if star_id is not None:
        return star_id

    # Object not recognized
    raise Error(
        f"Object '{object_name}' not recognized. "
        "Use planet names (Mercury, Venus, Mars, Jupiter, Saturn, etc.), "
        "integer planet IDs (2-9 for Mercury-Pluto), or fixed star names "
        "(Sirius, Regulus, Aldebaran, etc.)."
    )


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
    """Rise/set and visibility-window instants for heliacal_pheno_ut.

    The Sun's and the object's disc-center rise (event types 1/4) or set
    (2/3) are searched from four hours before ``jd``;
    the first/optimum/last visibility instants come from the visibility
    margin scanned on the night side of the Sun's event, and the
    public invalid-time sentinel marks instants that do not exist. Returns
    (rise_o, rise_s, lag, t_first, t_best, t_last, tvis, tb_yallop).
    """
    from .constants import (
        BIT_DISC_CENTER,
        CALC_RISE,
        CALC_SET,
        FLG_JPLEPH,
        FLG_MOSEPH,
    )
    from .eclipse import rise_trans

    eph_flags = flags & (FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH)
    rs = CALC_RISE if event_type in (1, 4) else CALC_SET
    is_star = bool(star_name)

    rc_s, tr_s = rise_trans(
        jd - 4.0 / 24.0, SUN, rs | BIT_DISC_CENTER, geopos3, 0.0, 0.0, eph_flags
    )
    rise_s = tr_s[0] if rc_s == 0 else _TJD_INVALID
    target = star_name if is_star else body
    rc_o, tr_o = rise_trans(
        jd - 4.0 / 24.0, target, rs | BIT_DISC_CENTER, geopos3, 0.0, 0.0, eph_flags
    )
    norise = rc_o != 0
    rise_o = tr_o[0] if rc_o == 0 else _TJD_INVALID

    lag = 0.0 if (norise or rc_s != 0) else rise_o - rise_s
    tb_yallop = _TJD_INVALID
    if not is_star and body == MOON and not norise and rc_s == 0:
        tb_yallop = (rise_o * 4.0 + rise_s * 5.0) / 9.0

    t_first = _TJD_INVALID
    t_best = _TJD_INVALID
    t_last = _TJD_INVALID
    tvis = _TJD_INVALID

    # Acronychal event types carry a visibility window only for the Moon and
    # inner planets under the public result contract.
    if event_type in (3, 4) and (is_star or body not in (1, 2, 3)):
        return rise_o, rise_s, lag, t_first, t_best, t_last, 0.0, tb_yallop
    if rc_s != 0:
        return rise_o, rise_s, lag, t_first, t_best, t_last, tvis, tb_yallop

    objname = star_name if is_star else _PHENO_BODY_NAMES.get(body)
    if objname is None:
        return rise_o, rise_s, lag, t_first, t_best, t_last, tvis, tb_yallop

    def _margin(t: float) -> float:
        # dret[0] already includes extinction to the object, so the margin is
        # limiting magnitude minus catalog magnitude.
        # Reject a Sun above the horizon: the visibility window is a
        # twilight phenomenon and must not follow a daytime maximum.
        retval, dret_v = vis_limit_mag(t, geopos3, atmo4, obs6, objname, flags)
        if retval < 0 or dret_v[3] >= 0.0:
            return -99.0
        return dret_v[0] - dret_v[7]

    win = 4.0 / 24.0
    if rs == CALC_RISE:
        lo, hi = rise_s - win, rise_s
    else:
        lo, hi = rise_s, rise_s + win

    # The visibility margin is not unimodal over the 4-hour bracket: it is
    # a narrow positive spike anchored at the twilight edge sitting on a
    # flat -99 plateau (Sun above horizon / object invisible), so a bare
    # golden-section search whose first probes land on the plateau slides
    # to the wrong edge and misses a real window. Locate the spike with a
    # coarse grid first (5-minute steps across the bracket), then refine
    # the best sample with a golden-section pass confined to its
    # neighbouring samples, where the margin IS unimodal.
    n_grid = 48
    step = (hi - lo) / n_grid
    best_i, best_m = 0, -1e99
    for i in range(n_grid + 1):
        m_i = _margin(lo + i * step)
        if m_i > best_m:
            best_i, best_m = i, m_i
    g_lo = lo + max(best_i - 1, 0) * step
    g_hi = lo + min(best_i + 1, n_grid) * step

    phi = (1.0 + math.sqrt(5.0)) / 2.0
    a = g_hi - (g_hi - g_lo) / phi
    b = g_lo + (g_hi - g_lo) / phi
    fa, fb = _margin(a), _margin(b)
    for _ in range(22):
        if fa > fb:
            g_hi, b, fb = b, a, fa
            a = g_hi - (g_hi - g_lo) / phi
            fa = _margin(a)
        else:
            g_lo, a, fa = a, b, fb
            b = g_lo + (g_hi - g_lo) / phi
            fb = _margin(b)
        if g_hi - g_lo < 3e-5:
            break
    t_peak = 0.5 * (g_lo + g_hi)
    m_peak = _margin(t_peak)
    if m_peak <= 0.0:
        return rise_o, rise_s, lag, t_first, t_best, t_last, tvis, tb_yallop
    t_best = t_peak

    def _cross(t_in: float, t_out: float) -> float:
        # If the object remains visible at a scan boundary, that boundary
        # itself bounds the finite search window.
        if _margin(t_out) > 0.0:
            return t_out
        for _ in range(22):
            mid = 0.5 * (t_in + t_out)
            if _margin(mid) > 0.0:
                t_in = mid
            else:
                t_out = mid
            if abs(t_out - t_in) < 3e-5:
                break
        return 0.5 * (t_in + t_out)

    t_first = _cross(t_peak, lo)
    t_last = _cross(t_peak, hi)
    if not norise:
        if rs == CALC_RISE and t_first != _TJD_INVALID:
            t_first = max(t_first, rise_o)
        elif rs == CALC_SET and t_last != _TJD_INVALID:
            t_last = min(t_last, rise_o)
    if t_first != _TJD_INVALID and t_last != _TJD_INVALID:
        tvis = t_last - t_first
    return rise_o, rise_s, lag, t_first, t_best, t_last, tvis, tb_yallop


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

    Returns:
        A pair ``(values, status)``. ``values`` has the following fixed layout:

        - 0 ``AltO``: object altitude at the observing site before refraction (deg).
        - 1 ``AppAltO``: refracted object altitude seen by the observer (deg).
        - 2 ``GeoAltO``: altitude obtained from the geocentric direction (deg).
        - 3 ``AziO``: object azimuth (deg).
        - 4 ``AltS`` and 5 ``AziS``: topocentric solar altitude and azimuth (deg).
        - 6 ``TAVact``: topocentric vertical separation from the Sun (deg).
        - 7 ``ARCVact``: geocentric vertical separation from the Sun (deg).
        - 8 ``DAZact``: signed object-minus-Sun azimuth separation (deg).
        - 9 ``ARCLact``: signed object-minus-Sun longitude separation (deg).
        - 10 ``kact``: atmospheric extinction coefficient.
        - 11 ``minTAV``: limiting topocentric vertical separation (deg).
        - 12-14 ``TfirstVR``, ``TbVR``, ``TlastVR``: beginning, optimum, and
          end of the visibility window (Julian days).
        - 15 ``TbYallop``: Yallop optimum-visibility time (Julian day).
        - 16 ``WMoon``: angular width of the lunar crescent (deg).
        - 17 ``qYal``: Yallop q statistic.
        - 18 ``qCrit``: numeric Yallop class, 1 through 6 (A through F).
        - 19 ``ParO``: object parallax (deg).
        - 20 ``Magn``: apparent magnitude.
        - 21 ``RiseO`` and 22 ``RiseS``: relevant object and solar horizon
          crossings (Julian days).
        - 23 ``Lag``: difference ``RiseO - RiseS`` (days).
        - 24 ``TvisVR``: duration of the visibility interval (days).
        - 25 ``LMoon``: angular length of the lunar crescent (deg).
        - 26-49: compatibility padding reserved for later fields.

        ``status`` contains the accepted calculation flags, or a negative error
        indicator.

    Raises:
        ValueError: If invalid body ID or event_type

    Example:
        >>> from libephemeris import julday, _heliacal_pheno_ut_pythonic, VENUS, HELIACAL_RISING
        >>> jd = julday(2024, 1, 1, 0)
        >>> # Get heliacal phenomena for Venus at Rome
        >>> dret, flag = _heliacal_pheno_ut_pythonic(jd, 41.9, 12.5, body=VENUS,
        ...                                 event_type=HELIACAL_RISING)
        >>> print(f"Object altitude: {dret[0]:.2f}, Sun altitude: {dret[4]:.2f}")

    References:
        - Schaefer (1990), PASP 102, 212-229, DOI 10.1086/132629,
          limiting-magnitude/observer model.
        - Yallop (1997), NAO Technical Note 69, only for the lunar-crescent
          q statistic and class boundaries placed in result slots 17-18.
        - Kasten & Young (1989), DOI 10.1364/AO.28.004735, optical airmass.
    """
    # --- LEB fast path ---
    from .state import get_leb_reader as _get_leb_reader

    _leb_rdr = _get_leb_reader()
    if _leb_rdr is not None:
        try:
            return _heliacal_pheno_ut_leb(
                _leb_rdr,
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
        except (KeyError, ValueError) as _leb_err:
            from .state import get_calc_mode

            if get_calc_mode() == "leb":
                raise
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
            f"Invalid event_type: {event_type}. Use HELIACAL_RISING, "
            "HELIACAL_SETTING, EVENING_FIRST, or MORNING_LAST."
        )

    # Check if this is a fixed star
    is_star = is_fixed_star(body)

    # Validate body - must be either a known planet or a fixed star
    if not is_star and body not in _PLANET_MAP:
        raise Error(f"illegal planet number {body}.")

    # Get ephemeris and timescale
    eph = get_planets()
    ts = get_timescale()

    sun = eph["sun"]
    earth = eph["earth"]

    # Get target - either planet or star
    target = None
    star_object = None
    star_magnitude = 0.0

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
            raise Error(f"Star ID {body} not found in catalog")

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
    # Keep the raw tuple too: the rise-window computation below needs it
    # after ``observer`` is rebound.
    _obs_raw = tuple(float(v) for v in tuple(observer)[:6])
    _obs_tuple6 = _obs_raw + (0.0,) * (6 - len(_obs_raw))

    # Create observer location
    observer = wgs84.latlon(lat, lon, altitude)

    # Initialize result array with 50 zeros
    dret = [0.0] * 50

    # Calculate positions at the given time
    t = ts.ut1_jd(jd)
    observer_at = earth + observer

    # Calculate Sun position
    sun_app = observer_at.at(t).observe(sun).apparent()
    sun_alt_topo, sun_az, _ = sun_app.altaz()
    sun_alt_deg = sun_alt_topo.degrees
    sun_az_deg = sun_az.degrees

    # Calculate body position (handle both planets and stars)
    if is_star and star_object is not None:
        body_app = observer_at.at(t).observe(star_object).apparent()
    else:
        body_app = observer_at.at(t).observe(target).apparent()
    body_alt_topo, body_az, body_dist = body_app.altaz()
    body_alt_deg = body_alt_topo.degrees
    body_az_deg = body_az.degrees

    # Get geocentric altitude (without refraction or topocentric correction)
    if is_star and star_object is not None:
        body_geo = earth.at(t).observe(star_object).apparent()
    else:
        body_geo = earth.at(t).observe(target).apparent()
    body_geo_ra, body_geo_dec, body_geo_dist = body_geo.radec(epoch=t)

    # Calculate geocentric altitude using hour angle
    # First get the local sidereal time. GAST is referenced to the true equinox
    # of date, so the RA must be of date too: radec(epoch=t) reduces to the
    # equator of date (apparent). Using the default J2000 radec() here mixed a
    # J2000 RA with an of-date LST, injecting the full J2000→date
    # precession+nutation in RA into the hour angle (GeoAltO/ARCVact/ARCLact off
    # by ~0.2-0.3 deg in 2024, growing with distance from J2000).
    gast = t.gast
    lst = gast + lon / 15.0  # Local sidereal time in hours
    ra_hours = body_geo_ra.hours
    ha = (lst - ra_hours) * 15.0  # Hour angle in degrees

    # Geocentric altitude calculation
    dec_rad = math.radians(body_geo_dec.degrees)
    lat_rad = math.radians(lat)
    ha_rad = math.radians(ha)

    sin_alt = math.sin(lat_rad) * math.sin(dec_rad) + math.cos(lat_rad) * math.cos(
        dec_rad
    ) * math.cos(ha_rad)
    sin_alt = max(-1.0, min(1.0, sin_alt))
    geo_alt_deg = math.degrees(math.asin(sin_alt))

    # Atmospheric refraction from true altitude, in arcminutes:
    # R = 1.02 / tan(h + 10.3/(h + 5.11)).  Sæmundsson, Þ. (1986),
    # "Astronomical Refraction", Sky & Telescope 72, 70; as given in
    # Meeus (1998), Astronomical Algorithms 2nd ed., ch. 16, eq. 16.4.
    if body_alt_deg > -1:
        refraction = 1.02 / math.tan(
            math.radians(body_alt_deg + 10.3 / (body_alt_deg + 5.11))
        )
        refraction /= 60.0  # Convert arcminutes to degrees
    else:
        refraction = 0.5  # Near-horizon fallback (below the formula's range)

    # Apparent altitude (with refraction)
    app_alt_deg = body_alt_deg + refraction

    # Calculate arcus visionis (altitude difference between body and Sun)
    # Topocentric arcus visionis
    tav_act = body_alt_deg - sun_alt_deg

    # Geocentric arcus visionis
    arcv_act = geo_alt_deg - sun_alt_deg

    # Public signed-azimuth convention: Sun minus object.
    daz_act = sun_az_deg - body_az_deg
    # Normalize to -180 to +180
    while daz_act > 180:
        daz_act -= 360
    while daz_act < -180:
        daz_act += 360

    # Get elongation (longitude difference) from Sun
    # For fixed stars, always calculate manually since pheno_ut doesn't support them
    if is_star:
        # Calculate elongation manually for stars
        sun_geo = earth.at(t).observe(sun).apparent()
        elongation = body_geo.separation_from(sun_geo).degrees
        magnitude = star_magnitude
        phase_angle = 0.0  # Stars don't have phase angle
    else:
        try:
            pheno = pheno_ut(jd, body, _heliacal_eph_flags(flags))
            elongation = pheno[2]  # Elongation
            magnitude = pheno[4]  # Visual magnitude
            phase_angle = pheno[0]  # Phase angle (index 0, not 1)
        except (ValueError, TypeError, ArithmeticError):
            # Calculate elongation manually
            sun_geo = earth.at(t).observe(sun).apparent()
            elongation = body_geo.separation_from(sun_geo).degrees
            magnitude = 0.0
            phase_angle = 0.0

    # ARCLact: actual arc length between body and Sun in horizontal coords
    # Computed as great-circle distance: sqrt(ARCV² + DAZ²)
    arcl_act = math.sqrt(arcv_act**2 + daz_act**2)

    # Use Schaefer model for extinction and arcus visionis (observer tuple
    # already parsed above, before the name was rebound).
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
        jd=jd,
    )
    k_act = schaefer.k_total
    min_tav = schaefer.arcus_visionis_required(magnitude)

    # Parallax of object (in degrees)
    # Fixed stars have essentially zero parallax
    if is_star:
        parallax = 0.0
    else:
        # Altitude parallax is the geocentric-minus-topocentric altitude, not
        # the horizontal parallax angle.
        parallax = geo_alt_deg - body_alt_deg

    # Rise/set of the disc centers and the visibility window use the shared
    # event-time helper and its public invalid-time sentinel.
    _star_name = ""
    if is_star:
        from .fixed_stars import STAR_CATALOG as _SC

        for _entry in _SC:
            if _entry.id == body:
                _star_name = _entry.name
                break
    _geopos3 = (lon, lat, altitude)
    _atmo4 = (
        pressure,
        temperature,
        humidity * 100.0 if humidity <= 1.0 else humidity,
        met_range,
    )
    # The visibility threshold is observer-dependent, so its window slots use
    # the caller's age, Snellen ratio, and optical parameters just like minTAV.
    _obs6 = _obs_tuple6
    (
        rise_o,
        rise_s,
        lag,
        t_first_vr,
        t_best_vr,
        t_last_vr,
        tvis_vr,
        t_b_yallop,
    ) = _pheno_rise_window(
        jd, body, _star_name, event_type, _geopos3, _atmo4, _obs6, flags
    )

    # Illumination percentage for all bodies
    # For planets: (1 + cos(phase_angle)) / 2 * 100
    if is_star:
        # Fixed stars are unresolved point sources; the public illumination
        # field uses the full-illumination convention for them.
        illumination = 100.0
    elif phase_angle > 0:
        illumination = (1.0 + math.cos(math.radians(phase_angle))) / 2.0 * 100.0
    else:
        illumination = 0.0

    # For Moon-specific calculations
    w_moon = 0.0  # Crescent width
    l_moon = 0.0  # Crescent length

    if body == MOON:
        # Calculate Moon phase and crescent geometry
        try:
            moon_pheno = pheno_ut(jd, MOON, _heliacal_eph_flags(flags))
            illumination = moon_pheno[1] * 100.0  # [1] = illuminated fraction 0-1

            # Crescent width (Yallop/Bruin): W = SD * (1 - cos ARCL), where
            # ARCL is the arc of light (Sun-Moon elongation, moon_pheno[2]) and
            # SD is the Moon's apparent semidiameter in arcmin (moon_pheno[3]
            # is the apparent diameter in degrees). Near new moon ARCL -> 0 so
            # W -> 0; near full moon ARCL -> 180 deg so W -> the disc diameter.
            #
            # Crescent length uses the actual distance-based apparent diameter
            # from pheno_ut (varies 0.49° at apogee to 0.56° at perigee):
            # L = pi * D / 2 (semicircle approximation).
            moon_diameter = (
                moon_pheno[3] if len(moon_pheno) > 3 and moon_pheno[3] > 0 else 0.5
            )
            arcl_deg = moon_pheno[2]
            sd_arcmin = moon_diameter / 2.0 * 60.0
            w_arcmin = sd_arcmin * (1.0 - math.cos(math.radians(arcl_deg)))
            w_moon = w_arcmin / 60.0  # stored in degrees (dret[16])
            l_moon = math.pi * moon_diameter / 2
        except (ValueError, TypeError, ArithmeticError):
            pass

    # Yallop q-test (for lunar crescent visibility)
    # q = (ARCV - (11.8371 - 6.3226*W + 0.7319*W^2 - 0.1018*W^3)) / 10
    # Simplified version
    if body == MOON:
        w = w_moon * 60.0  # Convert to arcminutes for formula
        q_criterion = 11.8371 - 6.3226 * w + 0.7319 * w**2 - 0.1018 * w**3
        q_yallop = (arcv_act - q_criterion) / 10.0
        # dret[18] (qCrit) is the Yallop visibility CLASS code (1..6), not the
        # raw arcus-visionis criterion polynomial (which stays internal above).
        q_crit = _yallop_visibility_code(q_yallop)
    else:
        q_yallop = 0.0
        q_crit = 0.0

    # Fill in the result array
    dret[0] = body_alt_deg  # AltO - topocentric altitude (unrefracted)
    dret[1] = app_alt_deg  # AppAltO - apparent altitude (refracted)
    dret[2] = geo_alt_deg  # GeoAltO - geocentric altitude
    dret[3] = body_az_deg  # AziO - azimuth of object
    dret[4] = sun_alt_deg  # AltS - topocentric altitude of Sun
    dret[5] = sun_az_deg  # AziS - azimuth of Sun
    dret[6] = tav_act  # TAVact - topocentric arcus visionis
    dret[7] = arcv_act  # ARCVact - geocentric arcus visionis
    dret[8] = daz_act  # DAZact - azimuth difference
    dret[9] = arcl_act  # ARCLact - actual arc distance between body and Sun
    dret[10] = k_act  # kact - extinction coefficient
    dret[11] = min_tav  # minTAV - minimum topocentric arcus visionis
    dret[12] = t_first_vr  # TfirstVR - first visibility time
    dret[13] = t_best_vr  # TbVR - best visibility time
    dret[14] = t_last_vr  # TlastVR - last visibility time
    dret[15] = t_b_yallop  # TbYallop - best time according to Yallop
    dret[16] = w_moon  # WMoon - crescent width
    dret[17] = q_yallop  # qYal - Yallop q-test value
    dret[18] = q_crit  # qCrit - Yallop visibility class code (1..6)
    dret[19] = parallax  # ParO - parallax of object
    dret[20] = magnitude  # Magn - magnitude
    dret[21] = rise_o  # RiseO - rise/set time of object
    dret[22] = rise_s  # RiseS - rise/set time of Sun
    dret[23] = lag  # Lag - time difference
    dret[24] = tvis_vr  # TvisVR - visibility duration
    dret[25] = l_moon  # LMoon - crescent length
    dret[26] = elongation  # elong - elongation from the Sun
    dret[27] = illumination  # Illum - illumination percentage
    # dret[28] onwards are reserved, already 0.0

    # Coerce to native Python floats: several entries (alt/az from the Skyfield
    # backend) are numpy.float64, and the public contract returns native floats.
    return tuple(float(v) for v in dret), flags


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
        Flat tuple of 50 floats with heliacal phenomena data.
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
            from .state import get_calc_mode

            if get_calc_mode() == "leb":
                raise
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
        raise Error("objname cannot be empty")

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
        raise Error("it makes no sense to call vis_limit_mag() for the Sun")

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
            raise Error(f"could not find star name {objname.lower()}: {e}") from e
    else:
        # Planet calculation
        if body_id is None:
            raise Error(f"Unknown object: {objname}")

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
            raise Error(f"illegal planet number {body_id}.")

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
