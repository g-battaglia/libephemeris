"""
Heliacal event calculations for libephemeris.

Calculates heliacal risings and settings, visual limiting magnitude,
and related heliacal phenomena for celestial bodies.

Functions:
- heliacal_ut: Find heliacal rising/setting time for a body
- swe_heliacal_ut: Alias for heliacal_ut
- heliacal_pheno_ut: Calculate detailed heliacal phenomena
- swe_heliacal_pheno_ut: Alias for heliacal_pheno_ut
- vis_limit_mag: Calculate visual limiting magnitude
- swe_vis_limit_mag: Alias for vis_limit_mag

Historical Note:
    Heliacal risings were crucial for ancient calendars. The heliacal
    rising of Sirius marked the Egyptian new year and predicted the
    Nile flood. Babylonians used heliacal events to track planetary
    positions without modern instruments.

References:
    - Reference documentation
    - Schoch "Planets in Mesopotamian Astral Science"
    - Ptolemy's criteria for heliacal visibility
    - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
"""

from __future__ import annotations

import math
from typing import Tuple, NamedTuple, Optional

import numpy as np

from .constants import (
    SE_SUN,
    SE_MOON,
    SE_MERCURY,
    SE_VENUS,
    SE_FIXSTAR_OFFSET,
    SEFLG_SPEED,
    SEFLG_SWIEPH,
)

# Inner planets (orbit inside Earth's orbit)
# These have both inferior conjunction (between Earth and Sun)
# and superior conjunction (behind the Sun)
INNER_PLANETS = {SE_MERCURY, SE_VENUS}


# =============================================================================
# SCHAEFER (1990) ATMOSPHERIC MODEL
# =============================================================================
#
# Complete implementation of Schaefer's visibility model. Based on:
#   - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
#   - Schaefer, B.E. (1993) "Astronomy and the Limits of Vision"
#
# The model calculates:
#   1. Atmospheric extinction (Rayleigh + Aerosol + Ozone)
#   2. Sky brightness (twilight + moonlight + airglow)
#   3. Ptolemaic visibility thresholds (arcus visionis)
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

    # Ptolemaic arcus visionis thresholds (degrees)
    # Based on ancient observations and Schaefer (1990) visibility model
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


class SchaeferModel:
    """
    Complete Schaefer atmospheric model for heliacal visibility.

    This class implements the Schaefer (1990) model for calculating:
    - Atmospheric extinction
    - Sky brightness
    - Limiting visual magnitude
    - Visibility thresholds

    The implementation follows the reference approach closely
    to ensure compatibility for heliacal event calculations.
    """

    def __init__(
        self,
        atmo: Optional[AtmosphericConditions] = None,
        observer: Optional[ObserverParams] = None,
    ):
        """
        Initialize the Schaefer model with atmospheric and observer conditions.

        Args:
            atmo: Atmospheric conditions (defaults to standard atmosphere)
            observer: Observer parameters (defaults to standard observer)
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

        # Precompute extinction coefficients
        self._compute_extinction_coefficients()

    def _compute_extinction_coefficients(self) -> None:
        """Compute atmospheric extinction coefficients."""
        C = SchaeferConstants

        # Pressure factor (relative to sea level)
        self.pressure_factor = self.atmo.pressure / 1013.25

        # Altitude factor for Rayleigh scattering
        alt_km = self.atmo.altitude / 1000.0
        self.altitude_rayleigh = math.exp(-alt_km / C.SCALE_HEIGHT_PRESSURE)

        # Altitude factor for aerosols
        self.altitude_aerosol = math.exp(-alt_km / C.SCALE_HEIGHT_AEROSOL)

        # Rayleigh scattering coefficient (depends on pressure)
        self.k_rayleigh = (
            C.K_RAYLEIGH_SEA_LEVEL * self.pressure_factor * self.altitude_rayleigh
        )

        # Aerosol coefficient
        if self.atmo.met_range >= 1.0:
            # Meteorological range given (km)
            # k_aerosol = 3.912 / V - k_rayleigh (Koschmieder's formula)
            self.k_aerosol = max(0.0, 3.912 / self.atmo.met_range - 0.106)
        elif 0 < self.atmo.met_range < 1.0:
            # ktot given directly
            self.k_total_override = self.atmo.met_range
            self.k_aerosol = max(
                0.0, self.k_total_override - self.k_rayleigh - C.K_OZONE
            )
        else:
            # Estimate from humidity
            # Schaefer model: aerosol increases with humidity
            # For standard 50% humidity, target k_aerosol ~ 0.1
            humidity_frac = self.atmo.humidity / 100.0
            self.k_aerosol = (
                C.K_AEROSOL_BASE * (1.0 + humidity_frac) * self.altitude_aerosol
            )

        # Ozone coefficient (constant)
        self.k_ozone = C.K_OZONE

        # Total extinction coefficient per airmass
        self.k_total = self.k_rayleigh + self.k_aerosol + self.k_ozone

    def airmass(self, altitude_deg: float) -> float:
        """
        Calculate airmass using Rozenberg's formula.

        This formula is more accurate near the horizon than the
        simple sec(z) formula.

        Args:
            altitude_deg: Altitude above horizon in degrees

        Returns:
            Airmass (1.0 at zenith, ~38 at horizon)
        """
        if altitude_deg <= -10.0:
            return 100.0  # Below horizon / extreme

        # Rozenberg (1966) formula for better horizon accuracy
        alt_rad = math.radians(max(altitude_deg, -5.0))
        sin_alt = math.sin(alt_rad)
        cos_alt = math.cos(alt_rad)

        # Rozenberg formula
        if altitude_deg > 0:
            airmass = 1.0 / (sin_alt + 0.025 * math.exp(-11.0 * sin_alt))
        else:
            # Extended formula for negative altitudes (refraction)
            airmass = 40.0  # Maximum reasonable airmass

        return min(airmass, 100.0)

    def extinction(self, altitude_deg: float) -> float:
        """
        Calculate total atmospheric extinction in magnitudes.

        Args:
            altitude_deg: Altitude above horizon in degrees

        Returns:
            Extinction in magnitudes
        """
        X = self.airmass(altitude_deg)
        return self.k_total * X

    def _ks_scattering(self, rho_deg: float) -> float:
        """
        Krisciunas & Schaefer (1991) Eq. 16 scattering function.

        Combines Rayleigh scattering (1 + cos²ρ dipole pattern) and
        Mie forward scattering (aerosol peak at small angles).

        This function determines how the sky brightness varies with
        angular distance ρ from an illuminating source (Sun or Moon).

        Args:
            rho_deg: Angular distance from light source in degrees

        Returns:
            Scattering intensity (nanoLamberts per foot-candle)
        """
        rho_deg = max(rho_deg, 0.5)  # Avoid singularity at 0°
        rho_rad = math.radians(rho_deg)
        cos_rho = math.cos(rho_rad)

        # Rayleigh scattering: symmetric pattern with minima at 90°
        rayleigh = 10**5.36 * (1.06 + cos_rho**2)

        # Mie scattering: strong forward peak, exponential falloff
        mie = 10 ** (6.15 - rho_deg / 40.0)

        return rayleigh + mie

    def sky_brightness_twilight(
        self, sun_alt: float, obj_alt: float, elongation: float
    ) -> float:
        """
        Calculate sky brightness contribution from twilight.

        Based on Schaefer (1993) and Krisciunas & Schaefer (1991).
        Returns brightness reduction in magnitudes relative to dark sky.

        The model accounts for:
        - Sun depression angle (primary twilight gradient)
        - Angular distance from Sun (Mie forward-scattering peak)
        - Object altitude (airmass through illuminated atmosphere)

        Args:
            sun_alt: Sun altitude in degrees (negative for below horizon)
            obj_alt: Object altitude in degrees
            elongation: Angular separation from Sun in degrees

        Returns:
            Sky brightness factor (magnitudes of limiting-mag reduction)
        """
        # Daylight case - return very large value
        if sun_alt >= 0:
            return 10.0  # Very bright

        # Base twilight brightness from sun depression angle.
        # Empirical relationship calibrated against reference heliacal
        # event dates:
        #   sun_alt = -6°  → sky = 3.5 (civil twilight, lim_mag ~ 3.0)
        #   sun_alt = -10° → sky = 1.6 (mid-nautical, lim_mag ~ 4.9)
        #   sun_alt = -14° → sky = 0.5 (late nautical, lim_mag ~ 6.0)
        #   sun_alt = -18° → sky = 0.0 (astronomical, lim_mag ~ 6.5)
        if sun_alt >= -6:
            # Civil twilight: very bright, rapid change
            base = 5.5 + sun_alt * 0.333  # 5.5 at horizon, 3.5 at -6
        elif sun_alt >= -10:
            # Early nautical twilight: still brightening rapidly
            # 3.5 at -6°, 1.6 at -10° (steeper than before)
            base = 3.5 + (sun_alt + 6) * 0.475
        elif sun_alt >= -14:
            # Late nautical twilight: moderate to dim
            # 1.6 at -10°, 0.5 at -14°
            base = 1.6 + (sun_alt + 10) * 0.275
        elif sun_alt >= -18:
            # Astronomical twilight: approaching dark
            # 0.5 at -14°, 0.0 at -18°
            base = 0.5 + (sun_alt + 14) * 0.125
        else:
            # Full night: dark sky
            base = 0.0

        # Horizon brightness: objects near the horizon look through
        # more of the illuminated atmosphere. During twilight, the
        # sky near the horizon is significantly brighter than at
        # higher altitudes. This is because scattered sunlight
        # concentrates in the lower atmosphere layers.
        if obj_alt < 15.0 and sun_alt > -18.0:
            # Objects below 15° altitude see extra brightness
            # ~0.4 mag at alt=2°, ~0.1 at alt=10°, 0 at alt=15°
            horizon_factor = 0.4 * ((15.0 - max(obj_alt, 0.5)) / 15.0) ** 1.5
            # Scale by twilight intensity
            twilight_scale = min(1.0, max(0.0, (sun_alt + 18.0) / 12.0))
            base += horizon_factor * twilight_scale

        # Mie forward-scattering penalty: aerosol-scattered sunlight
        # creates a bright glow around the Sun's position that reduces
        # the limiting magnitude for objects near the Sun.
        # Calibrated against reference heliacal event dates.
        if elongation < 60.0 and sun_alt > -18.0:
            # Quadratic ramp: ~1.0 mag at elong=10°, ~0.25 at 30°, 0 at 60°
            elong_factor = 1.0 * ((60.0 - elongation) / 60.0) ** 2
            # Scale by twilight intensity (no effect in full darkness)
            twilight_scale = min(1.0, max(0.0, (sun_alt + 18.0) / 12.0))
            base += elong_factor * twilight_scale

        return base

    def sky_brightness_moon(
        self,
        moon_alt: float,
        moon_phase: float,
        obj_alt: float,
        moon_obj_angle: float,
    ) -> float:
        """
        Calculate sky brightness contribution from the Moon.

        Uses the Krisciunas & Schaefer (1991) model. The Moon's
        contribution depends on its phase (illuminance), altitude
        (airmass extinction), angular distance to the object
        (scattering function), and the object's zenith distance
        (atmospheric path length for scattered moonlight).

        The illuminance model uses Allen (1976) lunar magnitudes
        as a function of phase angle, converted to ground-level
        illuminance in foot-candles after atmospheric extinction.

        Args:
            moon_alt: Moon altitude in degrees
            moon_phase: Moon phase fraction (0 = new, 1 = full)
            obj_alt: Object altitude in degrees
            moon_obj_angle: Angular separation between Moon and object in degrees

        Returns:
            Sky brightness contribution in magnitudes of limiting-mag reduction
        """
        if moon_alt <= 0:
            return 0.0

        # Convert phase fraction to phase angle (degrees).
        # phase=0 → α=180° (new), phase=1 → α=0° (full)
        alpha = (1.0 - moon_phase) * 180.0

        # Moon visual magnitude as function of phase angle
        # (Allen 1976, Krisciunas & Schaefer 1991 Eq. 9).
        # V_moon = -12.73 + 0.026|α| + 4e-9 * α^4
        abs_alpha = abs(alpha)
        moon_mag = -12.73 + 0.026 * abs_alpha + 4.0e-9 * abs_alpha**4

        # Convert magnitude to illuminance in foot-candles
        # (K&S 1991 Eq. 8): I = 10^(-0.4 * (V + 16.57))
        I_star = 10.0 ** (-0.4 * (moon_mag + 16.57))

        # Extinction of moonlight through atmosphere
        X_moon = self.airmass(moon_alt)
        I_ground = I_star * 10.0 ** (-0.4 * self.k_total * X_moon)

        # Scattering function at angle rho from Moon
        f_rho = self._ks_scattering(moon_obj_angle)

        # Airmass at object position
        X_obj = self.airmass(obj_alt)

        # Sky brightness from Moon (K&S 1991 Eq. 15/20):
        # B_moon = f(ρ) * I * 10^(-0.4*k*X_obj) * (1 - 10^(-0.4*k*X_obj))
        #
        # The term (1 - 10^(-0.4*k*X_obj)) represents the fraction of
        # atmosphere along the line of sight available for scattering.
        extinction_factor = 10.0 ** (-0.4 * self.k_total * X_obj)
        scatter_depth = 1.0 - extinction_factor

        # B_moon in nanoLamberts
        B_moon = f_rho * I_ground * extinction_factor * scatter_depth

        # Convert B_moon to magnitude reduction.
        # Dark sky brightness B_dark ≈ 145 nL (airglow) + 100 nL (zodiacal)
        # = 245 nL total. A full Moon at zenith adds ~400 nL at 60° away,
        # giving ~2.0 mag reduction. This matches observations.
        B_dark = 245.0  # nanoLamberts (airglow + zodiacal)
        if B_moon <= 0:
            return 0.0

        reduction = 2.5 * math.log10(1.0 + B_moon / B_dark)
        return max(0.0, min(5.0, reduction))

    def sky_brightness_total(
        self,
        sun_alt: float,
        moon_alt: float,
        moon_phase: float,
        obj_alt: float,
        sun_obj_angle: float,
        moon_obj_angle: float,
    ) -> float:
        """
        Calculate total sky brightness at object position.

        Combines twilight, moonlight, and airglow contributions
        additively in linear brightness space (K&S 1991 Eq. 1),
        then converts back to magnitude reduction.

        Returns brightness reduction in magnitudes (higher = brighter
        sky = lower limiting mag).

        Args:
            sun_alt: Sun altitude in degrees
            moon_alt: Moon altitude in degrees
            moon_phase: Moon phase (0 = new, 1 = full)
            obj_alt: Object altitude in degrees
            sun_obj_angle: Elongation from Sun in degrees
            moon_obj_angle: Angular separation from Moon in degrees

        Returns:
            Total sky brightness reduction in magnitudes
        """
        # Each component is already in magnitude reduction.
        # Convert each to linear brightness ratio (B/B_dark),
        # add them, then convert back to magnitude reduction.
        b_twi_mag = self.sky_brightness_twilight(sun_alt, obj_alt, sun_obj_angle)
        b_moon_mag = self.sky_brightness_moon(
            moon_alt, moon_phase, obj_alt, moon_obj_angle
        )

        # Convert from mag reduction to linear brightness excess:
        # mag_reduction = 2.5 * log10(1 + B/B_dark)
        # → B/B_dark = 10^(mag/2.5) - 1
        B_twi = 10.0 ** (b_twi_mag / 2.5) - 1.0 if b_twi_mag > 0 else 0.0
        B_moon = 10.0 ** (b_moon_mag / 2.5) - 1.0 if b_moon_mag > 0 else 0.0

        # Total brightness = sum of individual contributions
        B_total = B_twi + B_moon

        if B_total <= 0:
            return 0.0

        return 2.5 * math.log10(1.0 + B_total)

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
        Calculate the limiting apparent visual magnitude.

        Returns the faintest apparent magnitude visible at the given
        sky conditions. Compare directly against apparent magnitude
        (catalog magnitude + extinction).

        Based on Schaefer (1993) model considering:
        - Sky brightness from twilight, moonlight, and airglow
        - Observer eye characteristics

        Note: This returns a limit on *apparent* magnitude. The caller
        must add atmospheric extinction to the body's catalog magnitude
        before comparing. Extinction is NOT included here to avoid
        double-counting.

        Args:
            sun_alt: Sun altitude in degrees
            moon_alt: Moon altitude in degrees
            moon_phase: Moon phase (0 = new, 1 = full)
            obj_alt: Object altitude in degrees
            sun_obj_angle: Elongation from Sun in degrees
            moon_obj_angle: Angular separation from Moon in degrees

        Returns:
            Limiting apparent visual magnitude (fainter = larger number)
        """
        C = SchaeferConstants

        # Base limiting magnitude for perfect dark sky conditions
        m_lim_base = C.PERFECT_SKY_LIM_MAG  # 6.5

        # Sky brightness reduction (in magnitudes)
        sky_reduction = self.sky_brightness_total(
            sun_alt, moon_alt, moon_phase, obj_alt, sun_obj_angle, moon_obj_angle
        )

        # Calculate limiting magnitude
        # Brighter sky = lower limiting magnitude (can see fewer faint objects)
        # NOTE: extinction is NOT subtracted here — it is applied to the body's
        # catalog magnitude by the caller (is_visible). This avoids the previous
        # double-counting bug where extinction was subtracted from the limit
        # AND added to the body magnitude.
        m_lim = m_lim_base - sky_reduction

        # Observer corrections
        # Age effect (eyes deteriorate with age)
        age_factor = max(0, (self.observer.age - 30) / 10.0) * 0.1
        m_lim -= age_factor

        # Snellen ratio (better eyes see fainter)
        if self.observer.snellen > 0:
            snellen_factor = 2.5 * math.log10(self.observer.snellen)
            m_lim += snellen_factor

        # Binocular gain (about 0.8 mag)
        if self.observer.binocular:
            m_lim += 0.8

        # Telescope gain
        if self.observer.telescope_mag > 1.0 and self.observer.aperture > 0:
            # Gain = 5 * log10(D/7) where D is aperture in mm
            aperture_gain = 5.0 * math.log10(self.observer.aperture / 7.0)
            m_lim += aperture_gain * self.observer.transmission

        return m_lim

    def arcus_visionis_required(
        self, body_mag: float, body_type: str = "planet"
    ) -> float:
        """
        Calculate the required arcus visionis for visibility.

        The arcus visionis is the minimum altitude difference between
        a celestial body and the Sun for the body to be visible.

        Based on Ptolemaic criteria and Schaefer's analysis.

        Args:
            body_mag: Visual magnitude of the body
            body_type: Type of body ("planet" or "star")

        Returns:
            Required arcus visionis in degrees
        """
        C = SchaeferConstants

        # Base arcus visionis depends on magnitude
        # Brighter objects need less arcus visionis
        if body_mag < -3.0:
            # Very bright (Venus at brightest)
            base_av = 5.0
        elif body_mag < -1.0:
            # Bright (Jupiter, Sirius)
            base_av = 7.0
        elif body_mag < 0.5:
            # Moderately bright (Saturn, Canopus)
            base_av = 9.0
        elif body_mag < 1.5:
            # Medium brightness
            base_av = 10.0
        elif body_mag < 3.0:
            # Faint
            base_av = 12.0
        else:
            # Very faint
            base_av = 14.0 + (body_mag - 3.0) * 0.5

        # Adjust for atmospheric conditions
        # Poor conditions require larger arcus visionis
        av = base_av * (1.0 + 0.5 * (self.k_total - 0.25))

        return av

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

        # Calculate limiting apparent magnitude at body position
        lim_mag = self.limiting_magnitude(
            sun_alt=sun_alt,
            moon_alt=moon_alt,
            moon_phase=moon_phase,
            obj_alt=body_alt,
            sun_obj_angle=elongation,
            moon_obj_angle=moon_obj_angle,
        )

        # Apply extinction to body magnitude to get apparent magnitude
        apparent_mag = body_mag + self.extinction(body_alt)

        # Body is visible if apparent magnitude is brighter than limiting
        # magnitude minus the detection margin.
        return apparent_mag <= lim_mag - margin

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
        binocular=False,
        telescope_mag=0.0,
        aperture=0.0,
        transmission=1.0,
    )
    return SchaeferModel(atmo, observer)


# Outer planets (orbit outside Earth's orbit)
# These only have one type of conjunction (behind the Sun)
# SE_EVENING_FIRST and SE_MORNING_LAST are not applicable to these


def is_fixed_star(body: int) -> bool:
    """
    Check if a body ID corresponds to a fixed star.

    Args:
        body: Body ID constant

    Returns:
        True if the body is a fixed star (ID >= SE_FIXSTAR_OFFSET), False otherwise
    """
    return body >= SE_FIXSTAR_OFFSET


def _get_star_magnitude(star_id: int) -> float:
    """
    Get the visual magnitude of a fixed star from the catalog.

    Args:
        star_id: Star ID (SE_* constant, e.g., SE_SIRIUS)

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
        body: Planet ID constant (SE_MERCURY, SE_VENUS, etc.)

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
# entry-points (heliacal_ut, heliacal_pheno_ut, vis_limit_mag) when a
# LEBReader is available, avoiding the 3 GB DE-kernel load entirely.
# =============================================================================


def _star_name_from_id(star_id: int) -> str:
    """Return the canonical star name for a star body ID."""
    from .fixed_stars import STAR_CATALOG

    for entry in STAR_CATALOG:
        if entry.id == star_id:
            return entry.name
    raise ValueError(f"Star ID {star_id} not found in catalog")


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
    from .utils import azalt, SE_ECL2HOR

    if is_star:
        from .fixed_stars import swe_fixstar_ut

        assert star_name is not None
        pos, _name, _flags = swe_fixstar_ut(star_name, jd_ut, SEFLG_SPEED)
        ecl_lon, ecl_lat, ecl_dist = pos[0], pos[1], pos[2]
    else:
        from .fast_calc import _topo_ecliptic
        from .time_utils import swe_deltat

        jd_tt = jd_ut + swe_deltat(jd_ut)
        topo = _topo_ecliptic(reader, jd_tt, jd_ut, body_id, geopos_lonlat)
        ecl_lon, ecl_lat, ecl_dist = topo[0], topo[1], topo[2]

    az, alt_true, alt_app = azalt(
        jd_ut, SE_ECL2HOR, geopos_lonlat, atpress, attemp,
        (ecl_lon, ecl_lat, ecl_dist),
    )
    return az, alt_true, alt_app


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
    uses the already-LEB-ported swe_fixstar_ut.
    """
    if is_star:
        from .fixed_stars import swe_fixstar_ut

        assert star_name is not None
        pos, _name, _flags = swe_fixstar_ut(star_name, jd_ut, SEFLG_SPEED)
        return pos[0], pos[1], pos[2]
    else:
        from .fast_calc import _topo_ecliptic
        from .time_utils import swe_deltat

        jd_tt = jd_ut + swe_deltat(jd_ut)
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
) -> tuple:
    """LEB-backed implementation of heliacal_ut().

    Mirrors the Skyfield version but obtains all positions via
    _topo_ecliptic / swe_fixstar_ut / azalt / angular_separation.
    """
    from .constants import (
        SE_HELIACAL_RISING,
        SE_HELIACAL_SETTING,
        SE_EVENING_FIRST,
        SE_MORNING_LAST,
    )
    from .planets import _PLANET_MAP, swe_pheno_ut
    from .utils import angular_separation

    # --- validation (same as Skyfield path) ---
    if event_type not in (
        SE_HELIACAL_RISING, SE_HELIACAL_SETTING,
        SE_EVENING_FIRST, SE_MORNING_LAST,
    ):
        raise ValueError(
            f"Invalid event_type: {event_type}. Use SE_HELIACAL_RISING, "
            "SE_HELIACAL_SETTING, SE_EVENING_FIRST, or SE_MORNING_LAST."
        )
    if event_type in (SE_EVENING_FIRST, SE_MORNING_LAST):
        if is_fixed_star(body):
            raise ValueError(
                "SE_EVENING_FIRST and SE_MORNING_LAST are not valid for fixed stars. "
                "For fixed stars, use SE_HELIACAL_RISING or SE_HELIACAL_SETTING."
            )
        if not is_inner_planet(body):
            raise ValueError(
                "SE_EVENING_FIRST and SE_MORNING_LAST are only valid for inner planets "
                "(Mercury, Venus). For outer planets, use SE_HELIACAL_RISING or "
                "SE_HELIACAL_SETTING."
            )
    if body == SE_SUN:
        raise ValueError("SE_SUN is not valid for heliacal calculations")
    if body == SE_MOON:
        raise ValueError("SE_MOON is not valid for heliacal calculations")

    is_star = is_fixed_star(body)
    if not is_star and body not in _PLANET_MAP:
        raise ValueError(f"illegal planet number {body}.")

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
            reader, jd, SE_SUN, geopos, pressure, temperature,
        )
        az, body_alt, _ = _leb_body_altaz(
            reader, jd, body, geopos, pressure, temperature,
            is_star=is_star, star_name=star_name,
        )
        return sun_alt, body_alt, az

    def _get_elongation(jd: float) -> float:
        if not is_star:
            try:
                pheno = swe_pheno_ut(jd, body, flags)
                return pheno[2]
            except (ValueError, TypeError, ArithmeticError):
                pass
        # Manual elongation from ecliptic coords
        sun_ecl = _leb_ecliptic_pos(reader, jd, SE_SUN, geopos)
        body_ecl = _leb_ecliptic_pos(
            reader, jd, body, geopos,
            is_star=is_star, star_name=star_name,
        )
        return angular_separation(sun_ecl[0], sun_ecl[1], body_ecl[0], body_ecl[1])

    def _get_body_magnitude(jd: float) -> float:
        if is_star:
            return star_magnitude
        try:
            pheno = swe_pheno_ut(jd, body, flags)
            return pheno[4]
        except (ValueError, TypeError, ArithmeticError):
            return 0.0

    schaefer = create_schaefer_model(
        pressure=pressure,
        temperature=temperature,
        humidity=humidity * 100.0 if humidity <= 1.0 else humidity,
        altitude=altitude,
    )

    def _get_moon_data(jd: float):
        """Return (moon_alt, phase, moon_body_sep)."""
        _, moon_alt, _ = _leb_body_altaz(
            reader, jd, SE_MOON, geopos, pressure, temperature,
        )
        # Moon phase from ecliptic separation Sun-Moon
        sun_ecl = _leb_ecliptic_pos(reader, jd, SE_SUN, geopos)
        moon_ecl = _leb_ecliptic_pos(reader, jd, SE_MOON, geopos)
        elong_moon = angular_separation(
            sun_ecl[0], sun_ecl[1], moon_ecl[0], moon_ecl[1],
        )
        phase = (1.0 - math.cos(math.radians(elong_moon))) / 2.0

        body_ecl = _leb_ecliptic_pos(
            reader, jd, body, geopos,
            is_star=is_star, star_name=star_name,
        )
        moon_body_sep = angular_separation(
            moon_ecl[0], moon_ecl[1], body_ecl[0], body_ecl[1],
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
        visible = schaefer.is_visible(
            body_alt=body_alt, body_mag=body_mag, sun_alt=sun_alt,
            elongation=elongation, moon_alt=moon_alt, moon_phase=moon_phase,
            moon_obj_angle=moon_body_sep,
        )
        return visible, sun_alt, body_alt, elongation

    def _is_body_visible_no_moon(jd: float, margin: float = 0.5) -> bool:
        sun_alt, body_alt, _ = _get_altitudes(jd)
        elongation = _get_elongation(jd)
        if body_alt < 0 or sun_alt > 0:
            return False
        body_mag = _get_body_magnitude(jd)
        return schaefer.is_visible(
            body_alt=body_alt, body_mag=body_mag, sun_alt=sun_alt,
            elongation=elongation, moon_alt=-90.0, moon_phase=0.0,
            moon_obj_angle=180.0, margin=margin,
        )

    def _find_twilight_center(jd_day: float, morning: bool) -> float:
        best_ut = -1.0
        best_score = -999.0
        for h in range(24):
            jd_check = jd_day + h / 24.0
            sun_alt, _, _ = _get_altitudes(jd_check)
            if -22.0 < sun_alt < 2.0:
                jd_next = jd_day + (h + 1) / 24.0
                sun_next, _, _ = _get_altitudes(jd_next)
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
        center_ut = _find_twilight_center(jd_day, morning)
        if center_ut < 0:
            return False, 0.0
        sun_upper = -5.0 if morning else -2.0
        for dt_min in range(-180, 181, 15):
            ut_hour = center_ut + dt_min / 60.0
            jd_check = jd_day + ut_hour / 24.0
            sun_alt, body_alt, _ = _get_altitudes(jd_check)
            if -18.0 < sun_alt < sun_upper and body_alt > 0.5:
                if morning:
                    vis_margin = 0.70 if sun_alt <= -10.0 else 0.50
                else:
                    elong = _get_elongation(jd_check)
                    vis_margin = min(0.63 + elong * 0.006, 0.85)
                visible = _is_body_visible_no_moon(jd_check, margin=vis_margin)
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
            visible, sun_alt, body_alt, elong = _is_body_visible(jd_mid)
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
        max_days = 800
        lookback_jds = [jd_start_inner - i for i in range(1, 7)]
        lookback_vis = _batch_check_twilight_visibility(lookback_jds, morning=True)
        consecutive_invisible = 0
        for vis, _ in lookback_vis:
            if not vis:
                consecutive_invisible += 1
            else:
                break
        for batch_start in range(0, max_days, _HELIACAL_BATCH):
            batch_end = min(batch_start + _HELIACAL_BATCH, max_days)
            jd_days = [jd_start_inner + d for d in range(batch_start, batch_end)]
            vis_results = _batch_check_twilight_visibility(jd_days, morning=True)
            for vis, jd_vis in vis_results:
                if not vis:
                    consecutive_invisible += 1
                else:
                    if consecutive_invisible >= 5:
                        return _refine_heliacal_time(jd_vis, is_morning=True)
                    consecutive_invisible = 0
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
        max_days = 800
        was_invisible = False
        for batch_start in range(0, max_days, _HELIACAL_BATCH):
            batch_end = min(batch_start + _HELIACAL_BATCH, max_days)
            jd_days = [jd_start_inner + d for d in range(batch_start, batch_end)]
            vis_results = _batch_check_twilight_visibility(jd_days, morning=False)
            for vis, jd_vis in vis_results:
                if not vis:
                    was_invisible = True
                elif was_invisible:
                    return _refine_heliacal_time(jd_vis, is_morning=False)
        return 0.0

    def _search_morning_last(jd_start_inner: float) -> float:
        max_days = 800
        last_visible_jd = 0.0
        found_visible = False
        consecutive_invisible = 0
        for batch_start in range(0, max_days, _HELIACAL_BATCH):
            batch_end = min(batch_start + _HELIACAL_BATCH, max_days)
            jd_days = [jd_start_inner + d for d in range(batch_start, batch_end)]
            vis_results = _batch_check_twilight_visibility(jd_days, morning=True)
            for vis, jd_vis in vis_results:
                if vis:
                    last_visible_jd = jd_vis
                    found_visible = True
                    consecutive_invisible = 0
                elif found_visible:
                    consecutive_invisible += 1
                    if consecutive_invisible >= 3 and last_visible_jd > 0:
                        return _refine_heliacal_time(last_visible_jd, is_morning=True)
        return last_visible_jd if last_visible_jd > 0 else 0.0

    # --- main dispatch ---
    if event_type == SE_HELIACAL_RISING:
        jd_event = _search_heliacal_rising(jd_start)
    elif event_type == SE_HELIACAL_SETTING:
        jd_event = _search_heliacal_setting(jd_start)
    elif event_type == SE_EVENING_FIRST:
        jd_event = _search_evening_first(jd_start)
    elif event_type == SE_MORNING_LAST:
        jd_event = _search_morning_last(jd_start)
    else:
        jd_event = 0.0

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
) -> tuple:
    """LEB-backed implementation of heliacal_pheno_ut()."""
    from .constants import (
        SE_HELIACAL_RISING,
        SE_HELIACAL_SETTING,
        SE_EVENING_FIRST,
        SE_MORNING_LAST,
    )
    from .planets import _PLANET_MAP, swe_pheno_ut
    from .utils import angular_separation

    if event_type not in (
        SE_HELIACAL_RISING, SE_HELIACAL_SETTING,
        SE_EVENING_FIRST, SE_MORNING_LAST,
    ):
        raise ValueError(
            f"Invalid event_type: {event_type}. Use SE_HELIACAL_RISING, "
            "SE_HELIACAL_SETTING, SE_EVENING_FIRST, or SE_MORNING_LAST."
        )

    is_star = is_fixed_star(body)
    if not is_star and body not in _PLANET_MAP:
        raise ValueError(f"illegal planet number {body}.")

    star_name: str | None = None
    star_magnitude = 0.0
    if is_star:
        star_name = _star_name_from_id(body)
        star_magnitude = _get_star_magnitude(body)

    geopos = (lon, lat, altitude)

    # --- positions via LEB ---
    # Sun
    sun_az, sun_alt_deg, _ = _leb_body_altaz(
        reader, jd, SE_SUN, geopos, pressure, temperature,
    )

    # Body
    body_az_deg, body_alt_deg, body_app_alt = _leb_body_altaz(
        reader, jd, body, geopos, pressure, temperature,
        is_star=is_star, star_name=star_name,
    )

    # Atmospheric refraction (same formula as Skyfield path)
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
        reader, jd, body, geopos,
        is_star=is_star, star_name=star_name,
    )
    sun_ecl = _leb_ecliptic_pos(reader, jd, SE_SUN, geopos)

    # Geocentric altitude from equatorial coords (matching Skyfield path)
    from .constants import SEFLG_EQUATORIAL

    if is_star:
        from .fixed_stars import swe_fixstar_ut
        _eq_pos, _, _ = swe_fixstar_ut(star_name, jd, SEFLG_EQUATORIAL | SEFLG_SPEED)
        _body_ra_deg, _body_dec_deg = _eq_pos[0], _eq_pos[1]
    else:
        from .planets import swe_calc_ut as _scu_hp
        _eq_pos, _ = _scu_hp(jd, body, SEFLG_EQUATORIAL | SEFLG_SPEED)
        _body_ra_deg, _body_dec_deg = _eq_pos[0], _eq_pos[1]

    from .time_utils import _sidtime_internal
    from .state import get_timescale as _gts_hp
    from .cache import get_true_obliquity as _gto_hp

    _ts_hp = _gts_hp()
    _t_hp = _ts_hp.ut1_jd(jd)
    _jd_tt_hp = _t_hp.tt
    import erfa as _erfa_hp
    _dpsi_hp, _deps_hp = _erfa_hp.nut06a(2451545.0, _jd_tt_hp - 2451545.0)
    _eps0_hp = _erfa_hp.obl06(2451545.0, _jd_tt_hp - 2451545.0)
    _eps_hp = math.degrees(_eps0_hp + _deps_hp)
    _dpsi_deg_hp = math.degrees(_dpsi_hp)
    _lst_hours = _sidtime_internal(jd, geopos[0], _eps_hp, _dpsi_deg_hp)
    _ha_deg = (_lst_hours * 15.0 - _body_ra_deg) % 360.0
    _ha_rad = math.radians(_ha_deg)
    _dec_rad = math.radians(_body_dec_deg)
    _lat_rad = math.radians(geopos[1])
    _sin_alt = (math.sin(_lat_rad) * math.sin(_dec_rad)
                + math.cos(_lat_rad) * math.cos(_dec_rad) * math.cos(_ha_rad))
    geo_alt_deg = math.degrees(math.asin(max(-1.0, min(1.0, _sin_alt))))

    # Arcus visionis
    tav_act = body_alt_deg - sun_alt_deg
    arcv_act = geo_alt_deg - sun_alt_deg

    # Azimuth difference
    daz_act = body_az_deg - sun_az
    while daz_act > 180:
        daz_act -= 360
    while daz_act < -180:
        daz_act += 360
    daz_act = abs(daz_act)

    # Elongation and magnitude
    if is_star:
        elongation = angular_separation(
            sun_ecl[0], sun_ecl[1], body_ecl[0], body_ecl[1],
        )
        magnitude = star_magnitude
        phase_angle = 0.0
    else:
        try:
            pheno = swe_pheno_ut(jd, body, flags)
            elongation = pheno[2]
            magnitude = pheno[4]
            phase_angle = pheno[0]
        except (ValueError, TypeError, ArithmeticError):
            elongation = angular_separation(
                sun_ecl[0], sun_ecl[1], body_ecl[0], body_ecl[1],
            )
            magnitude = 0.0
            phase_angle = 0.0

    arcl_act = math.sqrt(arcv_act**2 + daz_act**2)

    schaefer = create_schaefer_model(
        pressure=pressure,
        temperature=temperature,
        humidity=humidity * 100.0 if humidity <= 1.0 else humidity,
        altitude=altitude,
    )
    k_act = schaefer.k_total
    min_tav = schaefer.arcus_visionis_required(magnitude)

    # Parallax
    if is_star:
        parallax = 0.0
    else:
        earth_radius_au = 6371.0 / 149597870.7
        dist_au = body_ecl[2]
        if dist_au > 0:
            parallax = math.degrees(math.asin(earth_radius_au / dist_au))
        else:
            parallax = 0.0

    # Rise/set estimates (same simplified logic as Skyfield path)
    is_morning = event_type in (SE_HELIACAL_RISING, SE_MORNING_LAST)
    rise_set_correction = -0.833
    lat_rad = math.radians(lat)

    if body_alt_deg != 0:
        alt_rate = 15.0 * math.cos(lat_rad)
        if alt_rate > 0:
            if is_morning:
                time_to_horizon = (body_alt_deg - rise_set_correction) / alt_rate
                rise_o = jd - time_to_horizon / 24.0
            else:
                time_to_horizon = (body_alt_deg - rise_set_correction) / alt_rate
                rise_o = jd + time_to_horizon / 24.0
        else:
            rise_o = jd
    else:
        rise_o = jd

    if sun_alt_deg != 0:
        alt_rate = 15.0 * math.cos(lat_rad)
        if alt_rate > 0:
            if is_morning:
                time_to_horizon = (sun_alt_deg - rise_set_correction) / alt_rate
                rise_s = jd - time_to_horizon / 24.0
            else:
                time_to_horizon = (sun_alt_deg - rise_set_correction) / alt_rate
                rise_s = jd + time_to_horizon / 24.0
        else:
            rise_s = jd
    else:
        rise_s = jd

    lag = rise_o - rise_s

    if sun_alt_deg < -6 and body_alt_deg > 0:
        tvis_vr = abs(sun_alt_deg + 6) / 15.0 / 24.0
    else:
        tvis_vr = 0.0

    if not is_star and phase_angle > 0:
        illumination = (1.0 + math.cos(math.radians(phase_angle))) / 2.0 * 100.0
    else:
        illumination = 0.0

    w_moon = 0.0
    l_moon = 0.0
    if body == SE_MOON:
        try:
            moon_pheno = swe_pheno_ut(jd, SE_MOON, flags)
            phase = moon_pheno[0]
            illumination = phase * 100.0
            if len(moon_pheno) > 1:
                pa_rad = math.radians(moon_pheno[1])
                w_moon = 15.0 * (1 - math.cos(pa_rad / 2)) / 60.0
            moon_diameter = (
                moon_pheno[3] if len(moon_pheno) > 3 and moon_pheno[3] > 0 else 0.5
            )
            l_moon = math.pi * moon_diameter / 2
        except (ValueError, TypeError, ArithmeticError):
            pass

    if body == SE_MOON:
        w = w_moon * 60.0
        q_criterion = 11.8371 - 6.3226 * w + 0.7319 * w**2 - 0.1018 * w**3
        q_yallop = (arcv_act - q_criterion) / 10.0
        q_crit = q_criterion
    else:
        q_yallop = 0.0
        q_crit = 0.0

    sun_rate = 15.0 * math.cos(lat_rad)
    if sun_rate > 0:
        if is_morning:
            t_first_vr = jd + (sun_alt_deg + 10) / sun_rate / 24.0
            t_best_vr = jd + (sun_alt_deg + 8) / sun_rate / 24.0
            t_last_vr = jd + (sun_alt_deg + 6) / sun_rate / 24.0
        else:
            t_first_vr = jd + (-6 - sun_alt_deg) / sun_rate / 24.0
            t_best_vr = jd + (-8 - sun_alt_deg) / sun_rate / 24.0
            t_last_vr = jd + (-10 - sun_alt_deg) / sun_rate / 24.0
    else:
        t_first_vr = jd
        t_best_vr = jd
        t_last_vr = jd

    t_b_yallop = t_best_vr

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
    dret[26] = phase_angle
    dret[27] = illumination

    return tuple(dret), flags


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
    from .planets import _PLANET_MAP, swe_pheno_ut
    from .fixed_stars import swe_fixstar2_ut, swe_fixstar2_mag
    from .utils import azalt, angular_separation, SE_ECL2HOR
    from .constants import (
        SE_SUN,
        SE_MOON,
        SE_HELFLAG_VISLIM_DARK,
        SE_HELFLAG_VISLIM_NOMOON,
        SE_HELFLAG_BELOW_HORIZON,
        SE_HELFLAG_PHOTOPIC,
        SE_HELFLAG_SCOTOPIC,
        SE_HELFLAG_MIXED,
    )

    if not objname:
        raise ValueError("objname cannot be empty")

    lon = geopos[0] if len(geopos) > 0 else 0.0
    lat = geopos[1] if len(geopos) > 1 else 0.0
    alt_m = geopos[2] if len(geopos) > 2 else 0.0

    pressure = atmo[0] if len(atmo) > 0 else 1013.25
    temperature = atmo[1] if len(atmo) > 1 else 15.0
    humidity_pct = atmo[2] if len(atmo) > 2 else 50.0
    met_range = atmo[3] if len(atmo) > 3 else 0.0

    observer_age = observer_tup[0] if len(observer_tup) > 0 else 36.0
    snellen_ratio = observer_tup[1] if len(observer_tup) > 1 else 1.0

    geopos_ll = (lon, lat, alt_m)

    # Sun position
    sun_az, sun_alt, _ = _leb_body_altaz(
        reader, tjdut, SE_SUN, geopos_ll, pressure, temperature,
    )

    # Moon position
    moon_az, moon_alt, _ = _leb_body_altaz(
        reader, tjdut, SE_MOON, geopos_ll, pressure, temperature,
    )

    # Determine object
    body_id = None
    is_star_flag = False

    try:
        body_id = int(objname)
    except ValueError:
        name_upper = objname.upper().strip()
        planet_names = {
            "SUN": SE_SUN, "MOON": SE_MOON,
            "MERCURY": 2, "VENUS": 3, "MARS": 4,
            "JUPITER": 5, "SATURN": 6, "URANUS": 7,
            "NEPTUNE": 8, "PLUTO": 9,
        }
        if name_upper in planet_names:
            body_id = planet_names[name_upper]
        else:
            is_star_flag = True

    obj_alt = 0.0
    obj_az = 0.0
    obj_mag = 0.0

    if is_star_flag:
        try:
            star_result, star_name_out, retflag = swe_fixstar2_ut(
                objname, tjdut, flags & 0xFF
            )
            star_lon = star_result[0]
            star_lat = star_result[1]

            try:
                star_mag_val, _star_name_mag = swe_fixstar2_mag(objname)
                obj_mag = star_mag_val
            except (ValueError, TypeError, ArithmeticError):
                obj_mag = 2.0

            hor_result = azalt(
                tjdut, SE_ECL2HOR, geopos_ll, pressure, temperature,
                (star_lon, star_lat, 1.0),
            )
            obj_az = hor_result[0]
            obj_alt = hor_result[1]
        except ValueError:
            raise
        except (KeyError, IndexError, OSError) as e:
            raise ValueError(f"could not find star name {objname.lower()}: {e}") from e
    else:
        if body_id is None:
            raise ValueError(f"Unknown object: {objname}")

        if body_id in _PLANET_MAP:
            az, alt_true, _ = _leb_body_altaz(
                reader, tjdut, body_id, geopos_ll, pressure, temperature,
            )
            obj_alt = alt_true
            obj_az = az

            try:
                pheno_result = swe_pheno_ut(tjdut, body_id, flags)
                obj_mag = pheno_result[4]
            except (ValueError, TypeError, ArithmeticError):
                obj_mag = 0.0
        else:
            raise ValueError(f"illegal planet number {body_id}.")

    if obj_alt < 0:
        dret = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        return float(SE_HELFLAG_BELOW_HORIZON), dret

    use_dark_sky = bool(flags & SE_HELFLAG_VISLIM_DARK)
    exclude_moon = bool(flags & SE_HELFLAG_VISLIM_NOMOON)

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
    )

    # Moon phase and angular separations via LEB ecliptic positions
    moon_phase = 0.0
    moon_obj_angle = 180.0
    sun_obj_angle = 180.0
    try:
        moon_pheno = swe_pheno_ut(tjdut, SE_MOON, flags & 0xFF)
        phase_angle = moon_pheno[0]
        moon_phase = (1.0 - math.cos(math.radians(phase_angle))) / 2.0

        if is_star_flag:
            moon_obj_angle = 90.0
            sun_obj_angle = 90.0
        else:
            body_ecl = _leb_ecliptic_pos(reader, tjdut, body_id, geopos_ll)
            moon_ecl = _leb_ecliptic_pos(reader, tjdut, SE_MOON, geopos_ll)
            sun_ecl = _leb_ecliptic_pos(reader, tjdut, SE_SUN, geopos_ll)
            moon_obj_angle = angular_separation(
                body_ecl[0], body_ecl[1], moon_ecl[0], moon_ecl[1],
            )
            sun_obj_angle = angular_separation(
                body_ecl[0], body_ecl[1], sun_ecl[0], sun_ecl[1],
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

    apparent_obj_mag = obj_mag + schaefer.extinction(obj_alt)

    if sun_alt >= -6:
        vision_type = SE_HELFLAG_PHOTOPIC
    elif sun_alt >= -12:
        vision_type = SE_HELFLAG_MIXED
    else:
        vision_type = SE_HELFLAG_SCOTOPIC

    dret = (
        limiting_mag, obj_alt, obj_az,
        sun_alt, sun_az, moon_alt, moon_az,
        apparent_obj_mag, 0.0, 0.0,
    )
    return float(vision_type), dret


def heliacal_ut(
    jd_start: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    pressure: float = 1013.25,
    temperature: float = 15.0,
    humidity: float = 0.5,
    body: int = SE_SUN,
    event_type: int = 1,
    flags: int = SEFLG_SWIEPH,
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
        body: Planet/body ID (SE_MERCURY, SE_VENUS, SE_MARS, SE_JUPITER, SE_SATURN,
              or fixed stars). Note: SE_SUN and SE_MOON are not valid for heliacal events.
        event_type: Type of heliacal event:
            - SE_HELIACAL_RISING (1): Morning first visibility (heliacal rising)
            - SE_HELIACAL_SETTING (2): Evening last visibility (heliacal setting)
            - SE_EVENING_FIRST (3): First evening visibility (after superior conjunction)
              Note: Only valid for inner planets (Mercury, Venus)
            - SE_MORNING_LAST (4): Last morning visibility (before superior conjunction)
              Note: Only valid for inner planets (Mercury, Venus)
        flags: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - jd_event: Julian Day (UT) of the heliacal event, or 0.0 if not found
            - retflag: Return flag (event_type on success, negative on error)

    Raises:
        ValueError: If invalid body ID, event_type, or if SE_EVENING_FIRST/SE_MORNING_LAST
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
        >>> from libephemeris import julday, heliacal_ut, SE_VENUS, SE_HELIACAL_RISING
        >>> jd = julday(2024, 1, 1, 0)
        >>> # Find next heliacal rising of Venus from Rome
        >>> jd_event, flag = heliacal_ut(jd, 41.9, 12.5, body=SE_VENUS,
        ...                              event_type=SE_HELIACAL_RISING)
        >>> print(f"Heliacal rising at JD {jd_event:.5f}")

    References:
        - Reference API: swe_heliacal_ut()
        - Schoch "Planets in Mesopotamian Astral Science"
        - Ptolemy's criteria for heliacal visibility
    """
    # --- LEB fast path ---
    from .state import get_leb_reader as _get_leb_reader

    _leb_rdr = _get_leb_reader()
    if _leb_rdr is not None:
        try:
            return _heliacal_ut_leb(
                _leb_rdr, jd_start, lat, lon, altitude,
                pressure, temperature, humidity, body, event_type, flags,
            )
        except KeyError:
            pass  # Body not in LEB file
        except ValueError:
            pass  # Date outside LEB range
    # --- END LEB fast path ---

    from .constants import (
        SE_HELIACAL_RISING,
        SE_HELIACAL_SETTING,
        SE_EVENING_FIRST,
        SE_MORNING_LAST,
    )
    from .planets import _PLANET_MAP, swe_pheno_ut
    from .state import get_planets, get_timescale
    from skyfield.api import wgs84

    # Validate event type
    if event_type not in (
        SE_HELIACAL_RISING,
        SE_HELIACAL_SETTING,
        SE_EVENING_FIRST,
        SE_MORNING_LAST,
    ):
        raise ValueError(
            f"Invalid event_type: {event_type}. Use SE_HELIACAL_RISING, "
            "SE_HELIACAL_SETTING, SE_EVENING_FIRST, or SE_MORNING_LAST."
        )

    # SE_EVENING_FIRST and SE_MORNING_LAST are only valid for inner planets
    # (Mercury and Venus) because they relate to superior conjunction visibility.
    # Outer planets only have one type of conjunction and these events don't apply.
    # Fixed stars only have heliacal rising and setting.
    if event_type in (SE_EVENING_FIRST, SE_MORNING_LAST):
        if is_fixed_star(body):
            raise ValueError(
                "SE_EVENING_FIRST and SE_MORNING_LAST are not valid for fixed stars. "
                "For fixed stars, use SE_HELIACAL_RISING or SE_HELIACAL_SETTING."
            )
        if not is_inner_planet(body):
            raise ValueError(
                "SE_EVENING_FIRST and SE_MORNING_LAST are only valid for inner planets "
                "(Mercury, Venus). For outer planets, use SE_HELIACAL_RISING or "
                "SE_HELIACAL_SETTING."
            )

    # Sun and Moon are not valid for heliacal events
    if body == SE_SUN:
        raise ValueError("SE_SUN is not valid for heliacal calculations")
    if body == SE_MOON:
        raise ValueError("SE_MOON is not valid for heliacal calculations")

    # Check if this is a fixed star
    is_star = is_fixed_star(body)

    # Validate body - must be either a known planet or a fixed star
    if not is_star and body not in _PLANET_MAP:
        raise ValueError(f"illegal planet number {body}.")

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
            raise ValueError(f"Star ID {body} not found in catalog")

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
                pheno = swe_pheno_ut(jd, body, flags)
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
            pheno = swe_pheno_ut(jd, body, flags)
            return pheno[4]  # Visual magnitude
        except (ValueError, TypeError, ArithmeticError):
            return 0.0  # Default to bright magnitude

    # Create Schaefer model for visibility calculations
    schaefer = create_schaefer_model(
        pressure=pressure,
        temperature=temperature,
        humidity=humidity * 100.0 if humidity <= 1.0 else humidity,  # Convert to %
        altitude=altitude,
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

        # Use Schaefer model for visibility check
        is_visible = schaefer.is_visible(
            body_alt=body_alt,
            body_mag=body_mag,
            sun_alt=sun_alt,
            elongation=elongation,
            moon_alt=moon_alt,
            moon_phase=moon_phase,
            moon_obj_angle=moon_body_sep,
        )

        return is_visible, sun_alt, body_alt, elongation

    def _is_body_visible_no_moon(jd: float, margin: float = 0.5) -> bool:
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

        # Check visibility with Moon forced below horizon
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

        Uses asymmetric detection margins:

        Morning (rising): sun-altitude-dependent margin.  At deep
        twilight (sun <= -10 deg) the empirical sky brightness model
        underestimates actual brightness, inflating the computed
        margin; a higher threshold (0.70 mag) compensates.  At
        shallower twilight (sun > -10 deg) the standard Schaefer
        (1993) 0.50 mag threshold applies, corresponding to ~90%
        detection probability.

        Evening (setting): elongation-dependent margin. At small
        elongations (<20°) scattered sunlight near the Sun makes
        the sky brighter than the model predicts, inflating the
        computed visibility margin. A lower threshold compensates.
        At large elongations the model is more accurate and a
        higher threshold matches reference transition points.

        Args:
            jd_day: JD at 0h UT of the day
            morning: True for dawn, False for dusk

        Returns:
            (visible, jd_best): whether body was found visible, and
            the JD of the first visibility moment found.
        """
        center_ut = _find_twilight_center(jd_day, morning)
        if center_ut < 0:
            return False, 0.0

        # Scan ±3 hours around center in 15-minute steps.
        # Use different sun altitude upper bounds for morning vs evening:
        # - Morning (rising): gate at -5° to prevent false detections during
        #   civil twilight where the sky brightness model is unreliable
        # - Evening (setting): gate at -2° because setting bodies are only
        #   visible briefly after sunset at shallow sun depressions
        sun_upper = -5.0 if morning else -2.0

        for dt_min in range(-180, 181, 15):
            ut_hour = center_ut + dt_min / 60.0
            jd_check = jd_day + ut_hour / 24.0
            sun_alt, body_alt, _ = _get_altitudes(jd_check)

            if -18.0 < sun_alt < sun_upper and body_alt > 0.5:
                if morning:
                    # Sun-altitude-dependent morning threshold.
                    #
                    # At deep twilight (sun <= -10 deg) the empirical
                    # sky brightness model underestimates the actual
                    # sky brightness, producing inflated visibility
                    # margins that cause premature detection by 1-2
                    # days.  A higher threshold (0.70 mag) compensates
                    # for the model's bias in the late nautical
                    # twilight regime.
                    #
                    # At shallower twilight (sun > -10 deg) the
                    # standard Schaefer (1993) 0.50 mag threshold
                    # applies, corresponding to ~90% detection
                    # probability.
                    if sun_alt <= -10.0:
                        vis_margin = 0.70
                    else:
                        vis_margin = 0.50
                else:
                    # Elongation-dependent evening threshold.
                    # At small elongations (<20°) scattered sunlight
                    # near the Sun is brighter than the model predicts,
                    # producing inflated visibility margins. A lower
                    # threshold compensates. At large elongations the
                    # model is more accurate, so a higher threshold
                    # matches reference transition points.
                    elong = _get_elongation(jd_check)
                    vis_margin = min(0.63 + elong * 0.006, 0.85)

                visible = _is_body_visible_no_moon(jd_check, margin=vis_margin)
                if visible:
                    return True, jd_check

        return False, 0.0

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
        # For planets swe_pheno_ut is expensive (~5ms). Sample every 10 days
        # and interpolate; magnitude changes by <0.1 mag over 10 days for
        # most planets, well within the Schaefer detection margin.
        day_mags: dict[int, float] = {}
        if is_star:
            for idx in valid_indices:
                day_mags[int(idx)] = star_magnitude
        else:
            _mag_step = 10
            _sample_idx = list(range(0, n_days, _mag_step))
            if _sample_idx[-1] != n_days - 1:
                _sample_idx.append(n_days - 1)
            _sample_mags = np.array(
                [_get_body_magnitude(float(jd_days_arr[i])) for i in _sample_idx]
            )
            _all_mags = np.interp(np.arange(n_days), _sample_idx, _sample_mags)
            for idx in valid_indices:
                day_mags[int(idx)] = float(_all_mags[int(idx)])

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

            if morning:
                vis_margin = 0.70 if sa <= -10.0 else 0.50
            else:
                elong = float(elongations[k])
                vis_margin = min(0.63 + elong * 0.006, 0.85)

            visible = schaefer.is_visible(
                body_alt=ba,
                body_mag=day_mags[day_i],
                sun_alt=sa,
                elongation=float(elongations[k]),
                moon_alt=-90.0,
                moon_phase=0.0,
                moon_obj_angle=180.0,
                margin=vis_margin,
            )
            if visible:
                results[day_i] = (True, float(scan_jds[k]))
                found.add(day_i)

        return results

    def _search_heliacal_rising(jd_start: float) -> float:
        """Search for heliacal rising using batched Skyfield computation."""
        max_days = 800

        # Look back to establish initial visibility state.
        lookback_jds = [jd_start - i for i in range(1, 7)]
        lookback_vis = _batch_check_twilight_visibility(lookback_jds, morning=True)
        consecutive_invisible = 0
        for vis, _ in lookback_vis:
            if not vis:
                consecutive_invisible += 1
            else:
                break

        for batch_start in range(0, max_days, _HELIACAL_BATCH):
            batch_end = min(batch_start + _HELIACAL_BATCH, max_days)
            jd_days = [jd_start + d for d in range(batch_start, batch_end)]
            vis_results = _batch_check_twilight_visibility(jd_days, morning=True)

            for vis, jd_vis in vis_results:
                if not vis:
                    consecutive_invisible += 1
                else:
                    if consecutive_invisible >= 5:
                        return _refine_heliacal_time(jd_vis, is_morning=True)
                    consecutive_invisible = 0

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
        max_days = 800
        was_invisible = False

        for batch_start in range(0, max_days, _HELIACAL_BATCH):
            batch_end = min(batch_start + _HELIACAL_BATCH, max_days)
            jd_days = [jd_start + d for d in range(batch_start, batch_end)]
            vis_results = _batch_check_twilight_visibility(jd_days, morning=False)

            for vis, jd_vis in vis_results:
                if not vis:
                    was_invisible = True
                elif was_invisible:
                    return _refine_heliacal_time(jd_vis, is_morning=False)

        return 0.0

    def _search_morning_last(jd_start: float) -> float:
        """Search for morning last visibility using batched computation."""
        max_days = 800
        last_visible_jd = 0.0
        found_visible = False
        consecutive_invisible = 0

        for batch_start in range(0, max_days, _HELIACAL_BATCH):
            batch_end = min(batch_start + _HELIACAL_BATCH, max_days)
            jd_days = [jd_start + d for d in range(batch_start, batch_end)]
            vis_results = _batch_check_twilight_visibility(jd_days, morning=True)

            for vis, jd_vis in vis_results:
                if vis:
                    last_visible_jd = jd_vis
                    found_visible = True
                    consecutive_invisible = 0
                elif found_visible:
                    consecutive_invisible += 1
                    if consecutive_invisible >= 3 and last_visible_jd > 0:
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

            visible, sun_alt, body_alt, elong = _is_body_visible(jd_mid)

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
    if event_type == SE_HELIACAL_RISING:
        jd_event = _search_heliacal_rising(jd_start)
    elif event_type == SE_HELIACAL_SETTING:
        jd_event = _search_heliacal_setting(jd_start)
    elif event_type == SE_EVENING_FIRST:
        jd_event = _search_evening_first(jd_start)
    elif event_type == SE_MORNING_LAST:
        jd_event = _search_morning_last(jd_start)
    else:
        jd_event = 0.0

    if jd_event > 0:
        return jd_event, event_type
    else:
        return 0.0, -1  # Not found


def swe_heliacal_ut(
    tjdut: float,
    geopos: tuple,
    atmo: tuple,
    observer: tuple,
    objname: str,
    eventtype: int,
    flags: int = SEFLG_SWIEPH,
) -> Tuple[float, float, float]:
    """
    Calculate heliacal rising or setting time for a celestial body.

    This function implements the swe_heliacal_ut() API,
    finding the Julian day of the next heliacal phenomenon after a given
    start date. It works between geographic latitudes 60S - 60N.

    Args:
        tjdut: Julian Day (UT) to start search from for the heliacal event.
        geopos: Geographic position as a sequence of at least 3 values:
            - [0]: Geographic longitude in degrees (east positive)
            - [1]: Geographic latitude in degrees (north positive)
            - [2]: Altitude above sea level in meters (eye height)
        atmo: Atmospheric conditions as a sequence of at least 4 values:
            - [0]: Atmospheric pressure in mbar/hPa (default: 1013.25)
            - [1]: Atmospheric temperature in degrees Celsius (default: 15)
            - [2]: Relative humidity in percent (0-100) (default: 40)
            - [3]: If >= 1: Meteorological Range in km
                   If 0 < value < 1: Total atmospheric coefficient (ktot)
                   If 0: Calculate ktot from other parameters
        observer: Observer description as a sequence of at least 6 values:
            - [0]: Age of observer in years (default: 36)
            - [1]: Snellen ratio of observer's eyes (default: 1.0 = normal)
            - [2]: 0 = monocular, 1 = binocular (used if SE_HELFLAG_OPTICAL_PARAMS)
            - [3]: Telescope magnification (0 = naked eye)
            - [4]: Optical aperture (telescope diameter) in mm
            - [5]: Optical transmission coefficient
        objname: Name of the celestial body. Can be:
            - Planet name: "Sun", "Moon", "Mercury", "Venus", "Mars",
              "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"
            - Fixed star name: "Sirius", "Aldebaran", "Regulus", etc.
            Note: Sun and Moon are not valid for heliacal calculations
              and will raise a ValueError.
        eventtype: Type of heliacal event:
            - SE_HELIACAL_RISING (1): Morning first visibility (heliacal rising)
              Exists for all visible planets and stars.
            - SE_HELIACAL_SETTING (2): Evening last visibility (heliacal setting)
              Exists for all visible planets and stars.
            - SE_EVENING_FIRST (3): First evening visibility after superior
              conjunction. Only valid for inner planets (Mercury, Venus).
            - SE_MORNING_LAST (4): Last morning visibility before superior
              conjunction. Only valid for inner planets (Mercury, Venus).
        flags: Calculation flags (bitmap). Contains ephemeris flags like
            SEFLG_SWIEPH, plus heliacal-specific flags:
            - SE_HELFLAG_OPTICAL_PARAMS (512): Use optical instrument parameters
              from observer[2-5]. Without this flag, those values are ignored.
            - SE_HELFLAG_NO_DETAILS (1024): Provide date only, skip visibility
              start/optimum/end details. Makes calculation faster.
            - SE_HELFLAG_VISLIM_DARK (4096): Function behaves as if Sun at nadir.
            - SE_HELFLAG_VISLIM_NOMOON (8192): Exclude Moon's brightness
              contribution (useful for calculating epoch heliacal dates).

    Returns:
        Tuple of 3 floats:
            - [0]: Start of visibility (Julian day number)
            - [1]: Optimum visibility (Julian day number),
                   zero if SE_HELFLAG_NO_DETAILS is set
            - [2]: End of visibility (Julian day number),
                   zero if SE_HELFLAG_NO_DETAILS is set

    Raises:
        ValueError: If invalid object_name, body ID, or event_type.
            Also raised if object_name is "Sun" or "Moon" which are not
            valid for heliacal calculations.

    Algorithm:
        The function searches for the moment when:
        1. The body is at a specific altitude above the horizon (arcus visionis)
        2. The Sun is at twilight position (typically -6 to -12 degrees below)
        3. The body's apparent magnitude is brighter than sky's limiting magnitude

        For heliacal rising (morning first):
        - Search forward for when body first becomes visible at dawn
        - Body must be above horizon while Sun is still below
        - Sky must be dark enough for the body to be seen

        For heliacal setting (evening last):
        - Search forward for when body is last visible at dusk
        - Body must be above horizon while Sun is setting
        - Sky brightness must not overwhelm the body's light

    Historical Note:
        Heliacal risings were crucial for ancient calendars. The heliacal
        rising of Sirius marked the Egyptian new year and predicted the
        Nile flood. Babylonians used heliacal events to track planetary
        positions without modern instruments.

    Example:
        >>> from libephemeris import julday, swe_heliacal_ut
        >>> from libephemeris import SE_HELIACAL_RISING, SEFLG_SWIEPH
        >>> jd = julday(2024, 1, 1, 0)
        >>> # Geographic position: Rome (lon, lat, altitude)
        >>> geopos = (12.5, 41.9, 0)
        >>> # Atmospheric conditions: standard
        >>> datm = (1013.25, 15.0, 40.0, 0.0)
        >>> # Observer: age 36, normal vision
        >>> dobs = (36.0, 1.0, 0, 0, 0, 0)
        >>> # Find next heliacal rising of Venus
        >>> jd1, jd2, jd3 = swe_heliacal_ut(jd, geopos, datm, dobs,
        ...                                  "Venus", SE_HELIACAL_RISING)
        >>> if jd1 > 0:
        ...     print(f"Venus heliacal rising at JD {jd1:.5f}")

    See Also:
        - heliacal_ut: Internal function using planet ID instead of name
        - heliacal_pheno_ut: Detailed heliacal phenomena calculation
        - vis_limit_mag: Calculate limiting visual magnitude

    References:
        - Reference documentation
        - Schoch "Planets in Mesopotamian Astral Science"
        - Ptolemy's criteria for heliacal visibility
    """
    from .constants import (
        SE_HELIACAL_RISING,
        SE_HELIACAL_SETTING,
        SE_EVENING_FIRST,
        SE_MORNING_LAST,
        SE_HELFLAG_NO_DETAILS,
    )

    # Validate event type
    if eventtype not in (
        SE_HELIACAL_RISING,
        SE_HELIACAL_SETTING,
        SE_EVENING_FIRST,
        SE_MORNING_LAST,
    ):
        raise ValueError(
            f"Invalid event_type: {eventtype}. Use SE_HELIACAL_RISING, "
            "SE_HELIACAL_SETTING, SE_EVENING_FIRST, or SE_MORNING_LAST."
        )

    # Parse geopos
    lon = geopos[0] if len(geopos) > 0 else 0.0
    lat = geopos[1] if len(geopos) > 1 else 0.0
    altitude = geopos[2] if len(geopos) > 2 else 0.0

    # Parse atmo with defaults per reference documentation
    pressure = atmo[0] if len(atmo) > 0 and atmo[0] > 0 else 1013.25
    temperature = atmo[1] if len(atmo) > 1 else 15.0
    humidity_pct = atmo[2] if len(atmo) > 2 else 40.0
    # atmo[3] is meteorological range / ktot, handled internally

    # Convert humidity from percent to 0-1 range for internal use
    humidity = humidity_pct / 100.0 if humidity_pct > 1.0 else humidity_pct

    # Parse objname to get body ID
    body_id = _parse_object_name(objname)

    # Call the internal heliacal_ut function
    jd_event, retflag = heliacal_ut(
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
    )

    # Build the result as 3 floats matching pyswisseph API
    jd1 = 0.0  # Start of visibility
    jd2 = 0.0  # Optimum visibility
    jd3 = 0.0  # End of visibility

    if jd_event > 0:
        jd1 = jd_event  # Start visibility

        # Calculate optimum and end visibility if details requested
        if not (flags & SE_HELFLAG_NO_DETAILS):
            # For detailed calculation, we estimate optimum and end times
            # based on typical visibility window durations
            # This is an approximation; full implementation would require
            # more complex calculations
            jd2 = jd_event  # Optimum (same as start for now)
            jd3 = jd_event  # End (same as start for now)

    return (jd1, jd2, jd3)


def _parse_object_name(object_name: str, allow_moon: bool = False) -> int:
    """
    Parse an object name string to get the corresponding body ID.

    Args:
        object_name: Name of planet or fixed star (e.g., "Venus", "Sirius")
        allow_moon: If True, allow Moon as a valid body (for heliacal_pheno_ut
            which supports Moon crescent calculations)

    Returns:
        Body ID (SE_* constant) for planets

    Raises:
        ValueError: If object name is not recognized or not valid for heliacal
    """
    # Import planet constants
    from .constants import (
        SE_SUN,
        SE_MOON,
        SE_MERCURY,
        SE_VENUS,
        SE_MARS,
        SE_JUPITER,
        SE_SATURN,
        SE_URANUS,
        SE_NEPTUNE,
        SE_PLUTO,
    )

    # Normalize name
    name_upper = object_name.upper().strip()

    # Planet name to ID mapping
    planet_map = {
        "SUN": SE_SUN,
        "MOON": SE_MOON,
        "MERCURY": SE_MERCURY,
        "VENUS": SE_VENUS,
        "MARS": SE_MARS,
        "JUPITER": SE_JUPITER,
        "SATURN": SE_SATURN,
        "URANUS": SE_URANUS,
        "NEPTUNE": SE_NEPTUNE,
        "PLUTO": SE_PLUTO,
    }

    # Check if it's a planet
    if name_upper in planet_map:
        body_id = planet_map[name_upper]

        # Sun is not valid for heliacal calculations
        if body_id == SE_SUN:
            raise ValueError("Sun is not valid for heliacal calculations")
        if body_id == SE_MOON and not allow_moon:
            raise ValueError("Moon is not valid for heliacal calculations")

        return body_id

    # Try to parse as planet number
    try:
        body_id = int(object_name)
        if body_id == SE_SUN:
            raise ValueError("Sun is not valid for heliacal calculations")
        if body_id == SE_MOON and not allow_moon:
            raise ValueError("Moon is not valid for heliacal calculations")
        return body_id
    except ValueError:
        pass

    # Try to resolve as a fixed star name
    from .fixed_stars import resolve_star_name

    star_id = resolve_star_name(object_name)
    if star_id is not None:
        return star_id

    # Object not recognized
    raise ValueError(
        f"Object '{object_name}' not recognized. "
        "Use planet names (Mercury, Venus, Mars, Jupiter, Saturn, etc.), "
        "planet IDs (2-9 for Mercury-Pluto), or fixed star names "
        "(Sirius, Regulus, Aldebaran, etc.)."
    )


def heliacal_pheno_ut(
    jd: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    pressure: float = 1013.25,
    temperature: float = 15.0,
    humidity: float = 0.5,
    body: int = SE_SUN,
    event_type: int = 1,
    flags: int = SEFLG_SWIEPH,
) -> Tuple[Tuple[float, ...], int]:
    """
    Provides data relevant for the calculation of heliacal risings and settings.

    This function calculates detailed phenomena associated with heliacal events,
    including altitudes, azimuths, arcus visionis, magnitude, visibility times,
    and other parameters used in heliacal visibility calculations.

    Args:
        jd: Julian Day (UT) for the calculation
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        pressure: Atmospheric pressure in mbar/hPa for refraction (default 1013.25)
        temperature: Temperature in Celsius for refraction (default 15)
        humidity: Relative humidity 0.0-1.0 for atmospheric extinction (default 0.5)
        body: Planet/body ID (SE_MERCURY, SE_VENUS, SE_MARS, SE_JUPITER, SE_SATURN, etc.)
        event_type: Type of heliacal event:
            - SE_HELIACAL_RISING (1): Morning first visibility (heliacal rising)
            - SE_HELIACAL_SETTING (2): Evening last visibility (heliacal setting)
            - SE_EVENING_FIRST (3): First evening visibility (after superior conjunction)
            - SE_MORNING_LAST (4): Last morning visibility (before superior conjunction)
        flags: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - dret: Tuple of 50 floats with heliacal phenomena data:
                - 0: AltO [deg] topocentric altitude of object (unrefracted)
                - 1: AppAltO [deg] apparent altitude of object (refracted)
                - 2: GeoAltO [deg] geocentric altitude of object
                - 3: AziO [deg] azimuth of object
                - 4: AltS [deg] topocentric altitude of Sun
                - 5: AziS [deg] azimuth of Sun
                - 6: TAVact [deg] actual topocentric arcus visionis
                - 7: ARCVact [deg] actual (geocentric) arcus visionis
                - 8: DAZact [deg] actual difference between object's and sun's azimuth
                - 9: ARCLact [deg] actual longitude difference between object and sun
                - 10: kact [-] extinction coefficient
                - 11: minTAV [deg] smallest topocentric arcus visionis
                - 12: TfirstVR [JDN] first time object is visible, according to VR
                - 13: TbVR [JDN] optimum time the object is visible, according to VR
                - 14: TlastVR [JDN] last time object is visible, according to VR
                - 15: TbYallop [JDN] best time the object is visible, according to Yallop
                - 16: WMoon [deg] crescent width of Moon
                - 17: qYal [-] q-test value of Yallop
                - 18: qCrit [-] q-test criterion of Yallop
                - 19: ParO [deg] parallax of object
                - 20: Magn [-] magnitude of object
                - 21: RiseO [JDN] rise/set time of object
                - 22: RiseS [JDN] rise/set time of Sun
                - 23: Lag [JDN] rise/set time of object minus rise/set time of Sun
                - 24: TvisVR [JDN] visibility duration
                - 25: LMoon [deg] crescent length of Moon
                - 26-49: Reserved for future use
            - retflag: Return flag (flags on success, negative on error)

    Raises:
        ValueError: If invalid body ID or event_type

    Example:
        >>> from libephemeris import julday, heliacal_pheno_ut, SE_VENUS, SE_HELIACAL_RISING
        >>> jd = julday(2024, 1, 1, 0)
        >>> # Get heliacal phenomena for Venus at Rome
        >>> dret, flag = heliacal_pheno_ut(jd, 41.9, 12.5, body=SE_VENUS,
        ...                                 event_type=SE_HELIACAL_RISING)
        >>> print(f"Object altitude: {dret[0]:.2f}, Sun altitude: {dret[4]:.2f}")

    References:
        - Reference API: swe_heliacal_pheno_ut()
        - Schoch "Planets in Mesopotamian Astral Science"
    """
    # --- LEB fast path ---
    from .state import get_leb_reader as _get_leb_reader

    _leb_rdr = _get_leb_reader()
    if _leb_rdr is not None:
        try:
            return _heliacal_pheno_ut_leb(
                _leb_rdr, jd, lat, lon, altitude,
                pressure, temperature, humidity, body, event_type, flags,
            )
        except KeyError:
            pass  # Body not in LEB file
        except ValueError:
            pass  # Date outside LEB range
    # --- END LEB fast path ---

    from .constants import (
        SE_HELIACAL_RISING,
        SE_HELIACAL_SETTING,
        SE_EVENING_FIRST,
        SE_MORNING_LAST,
    )
    from .planets import _PLANET_MAP, swe_pheno_ut
    from .state import get_planets, get_timescale
    from skyfield.api import wgs84

    # Validate event type
    if event_type not in (
        SE_HELIACAL_RISING,
        SE_HELIACAL_SETTING,
        SE_EVENING_FIRST,
        SE_MORNING_LAST,
    ):
        raise ValueError(
            f"Invalid event_type: {event_type}. Use SE_HELIACAL_RISING, "
            "SE_HELIACAL_SETTING, SE_EVENING_FIRST, or SE_MORNING_LAST."
        )

    # Check if this is a fixed star
    is_star = is_fixed_star(body)

    # Validate body - must be either a known planet or a fixed star
    if not is_star and body not in _PLANET_MAP:
        raise ValueError(f"illegal planet number {body}.")

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
            raise ValueError(f"Star ID {body} not found in catalog")

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
    body_geo_ra, body_geo_dec, body_geo_dist = body_geo.radec()

    # Calculate geocentric altitude using hour angle
    # First get the local sidereal time
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

    # Calculate atmospheric refraction
    # Use simplified formula: R = 1.02 / tan(h + 10.3/(h + 5.11)) in arcminutes
    if body_alt_deg > -1:
        refraction = 1.02 / math.tan(
            math.radians(body_alt_deg + 10.3 / (body_alt_deg + 5.11))
        )
        refraction /= 60.0  # Convert arcminutes to degrees
    else:
        refraction = 0.5  # Near horizon approximation

    # Apparent altitude (with refraction)
    app_alt_deg = body_alt_deg + refraction

    # Calculate arcus visionis (altitude difference between body and Sun)
    # Topocentric arcus visionis
    tav_act = body_alt_deg - sun_alt_deg

    # Geocentric arcus visionis
    arcv_act = geo_alt_deg - sun_alt_deg

    # Azimuth difference (absolute value, matching reference API convention)
    daz_act = body_az_deg - sun_az_deg
    # Normalize to -180 to +180
    while daz_act > 180:
        daz_act -= 360
    while daz_act < -180:
        daz_act += 360
    daz_act = abs(daz_act)

    # Get elongation (longitude difference) from Sun
    # For fixed stars, always calculate manually since swe_pheno_ut doesn't support them
    if is_star:
        # Calculate elongation manually for stars
        sun_geo = earth.at(t).observe(sun).apparent()
        elongation = body_geo.separation_from(sun_geo).degrees
        magnitude = star_magnitude
        phase_angle = 0.0  # Stars don't have phase angle
    else:
        try:
            pheno = swe_pheno_ut(jd, body, flags)
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

    # Use Schaefer model for extinction and arcus visionis
    schaefer = create_schaefer_model(
        pressure=pressure,
        temperature=temperature,
        humidity=humidity * 100.0 if humidity <= 1.0 else humidity,
        altitude=altitude,
    )
    k_act = schaefer.k_total
    min_tav = schaefer.arcus_visionis_required(magnitude)

    # Parallax of object (in degrees)
    # Fixed stars have essentially zero parallax
    if is_star:
        parallax = 0.0
    else:
        # Parallax = arcsin(Earth_radius / distance)
        # Earth radius ~ 6371 km, distance in AU (1 AU ~ 149,597,870.7 km)
        earth_radius_au = 6371.0 / 149597870.7
        if body_geo_dist.au > 0:
            parallax = math.degrees(math.asin(earth_radius_au / body_geo_dist.au))
        else:
            parallax = 0.0

    # Calculate rise/set times for object and Sun
    # For morning events, we look for rise times; for evening, set times
    is_morning = event_type in (SE_HELIACAL_RISING, SE_MORNING_LAST)

    # Use a simple estimate for rise/set time based on altitude
    # Rise/set occurs when altitude crosses 0 (corrected for refraction)
    # Time rate: approximately 1 degree of altitude per 4 minutes (at mid-latitudes)
    rise_set_correction = -0.833  # Standard refraction + semidiameter for Sun

    # Estimate object rise/set time
    if body_alt_deg != 0:
        # Rough estimate: time to rise/set based on altitude rate
        # Altitude rate ~ 15 cos(lat) per hour at the horizon
        alt_rate = 15.0 * math.cos(lat_rad)  # degrees per hour
        if alt_rate > 0:
            if is_morning:
                # Object rising: how long until it crosses horizon
                time_to_horizon = (
                    body_alt_deg - rise_set_correction
                ) / alt_rate  # hours
                rise_o = jd - time_to_horizon / 24.0
            else:
                # Object setting
                time_to_horizon = (
                    body_alt_deg - rise_set_correction
                ) / alt_rate  # hours
                rise_o = jd + time_to_horizon / 24.0
        else:
            rise_o = jd
    else:
        rise_o = jd

    # Estimate Sun rise/set time
    if sun_alt_deg != 0:
        alt_rate = 15.0 * math.cos(lat_rad)
        if alt_rate > 0:
            if is_morning:
                time_to_horizon = (sun_alt_deg - rise_set_correction) / alt_rate
                rise_s = jd - time_to_horizon / 24.0
            else:
                time_to_horizon = (sun_alt_deg - rise_set_correction) / alt_rate
                rise_s = jd + time_to_horizon / 24.0
        else:
            rise_s = jd
    else:
        rise_s = jd

    # Lag time (object rise - Sun rise)
    lag = rise_o - rise_s

    # Visibility duration estimate (simplified)
    # Based on how long the object is above horizon while Sun is below
    if sun_alt_deg < -6 and body_alt_deg > 0:
        # Object visible during civil twilight or darker
        # Estimate visibility window
        tvis_vr = abs(sun_alt_deg + 6) / 15.0 / 24.0  # days
    else:
        tvis_vr = 0.0

    # Illumination percentage for all bodies
    # For planets: (1 + cos(phase_angle)) / 2 * 100
    if not is_star and phase_angle > 0:
        illumination = (1.0 + math.cos(math.radians(phase_angle))) / 2.0 * 100.0
    else:
        illumination = 0.0

    # For Moon-specific calculations
    w_moon = 0.0  # Crescent width
    l_moon = 0.0  # Crescent length

    if body == SE_MOON:
        # Calculate Moon phase and crescent geometry
        try:
            moon_pheno = swe_pheno_ut(jd, SE_MOON, flags)
            phase = moon_pheno[0]  # Phase 0-1
            illumination = phase * 100.0  # Percentage

            # Crescent width approximation (Danjon's formula)
            # W = 15 * (1 - cos(phase_angle/2)) arcminutes (approximate)
            if len(moon_pheno) > 1:
                pa_rad = math.radians(moon_pheno[1])  # Phase angle in radians
                w_moon = 15.0 * (1 - math.cos(pa_rad / 2)) / 60.0  # In degrees

            # Crescent length (semicircle approximation)
            # L = pi * D / 2 where D is diameter
            # Use actual distance-based apparent diameter from swe_pheno_ut
            # (varies 0.49° at apogee to 0.56° at perigee)
            moon_diameter = (
                moon_pheno[3] if len(moon_pheno) > 3 and moon_pheno[3] > 0 else 0.5
            )
            l_moon = math.pi * moon_diameter / 2
        except (ValueError, TypeError, ArithmeticError):
            pass

    # Yallop q-test (for lunar crescent visibility)
    # q = (ARCV - (11.8371 - 6.3226*W + 0.7319*W^2 - 0.1018*W^3)) / 10
    # Simplified version
    if body == SE_MOON:
        w = w_moon * 60.0  # Convert to arcminutes for formula
        q_criterion = 11.8371 - 6.3226 * w + 0.7319 * w**2 - 0.1018 * w**3
        q_yallop = (arcv_act - q_criterion) / 10.0
        q_crit = q_criterion
    else:
        q_yallop = 0.0
        q_crit = 0.0

    # First, best, and last visibility times (VR = visibility rule)
    # Simplified: based on Sun altitude thresholds
    # First visibility: Sun at about -10
    # Best visibility: Sun at about -8
    # Last visibility: Sun at about -6
    sun_rate = 15.0 * math.cos(lat_rad)  # degrees per hour
    if sun_rate > 0:
        if is_morning:
            t_first_vr = jd + (sun_alt_deg + 10) / sun_rate / 24.0
            t_best_vr = jd + (sun_alt_deg + 8) / sun_rate / 24.0
            t_last_vr = jd + (sun_alt_deg + 6) / sun_rate / 24.0
        else:
            t_first_vr = jd + (-6 - sun_alt_deg) / sun_rate / 24.0
            t_best_vr = jd + (-8 - sun_alt_deg) / sun_rate / 24.0
            t_last_vr = jd + (-10 - sun_alt_deg) / sun_rate / 24.0
    else:
        t_first_vr = jd
        t_best_vr = jd
        t_last_vr = jd

    # Best time according to Yallop (for Moon)
    t_b_yallop = t_best_vr  # Use same as best VR for simplicity

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
    dret[18] = q_crit  # qCrit - Yallop criterion
    dret[19] = parallax  # ParO - parallax of object
    dret[20] = magnitude  # Magn - magnitude
    dret[21] = rise_o  # RiseO - rise/set time of object
    dret[22] = rise_s  # RiseS - rise/set time of Sun
    dret[23] = lag  # Lag - time difference
    dret[24] = tvis_vr  # TvisVR - visibility duration
    dret[25] = l_moon  # LMoon - crescent length
    dret[26] = phase_angle  # CVAact (using phase angle)
    dret[27] = illumination  # Illum - illumination percentage
    # dret[28] onwards are reserved, already 0.0

    return tuple(dret), flags


def swe_heliacal_pheno_ut(
    tjdut: float,
    geopos: tuple,
    atmo: tuple,
    observer: tuple,
    objname: str,
    eventtype: int,
    flags: int = SEFLG_SWIEPH,
) -> Tuple[float, ...]:
    """
    Provides data relevant for the calculation of heliacal risings and settings.

    This is the reference-compatible wrapper around heliacal_pheno_ut(). It
    accepts the same parameter layout as pyswisseph's swe_heliacal_pheno_ut
    and returns a flat 50-element tuple.

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

    # Convert humidity from percent to 0-1 range
    humidity = humidity_pct / 100.0 if humidity_pct > 1.0 else humidity_pct

    # Parse objname to body ID (Moon is allowed for pheno calculations)
    body_id = _parse_object_name(objname, allow_moon=True)

    # Call internal function
    dret, retflag = heliacal_pheno_ut(
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
    )

    # Return flat 50-tuple (matching reference API)
    return dret


def vis_limit_mag(
    tjdut: float,
    geopos: tuple,
    atmo: tuple,
    observer: tuple,
    objname: str,
    flags: int = SEFLG_SWIEPH,
) -> Tuple[float, Tuple[float, ...]]:
    """
    Calculate the limiting visual magnitude for observing a celestial body.

    This function determines whether a celestial body (planet, star, etc.)
    is visible given the current sky brightness conditions, atmospheric
    parameters, and observer characteristics. It returns both the visibility
    status and detailed information about the observation conditions.

    Args:
        tjdut: Julian Day (UT) for the observation time
        geopos: Geographic position as a sequence:
            - [0]: Geographic longitude in degrees (east positive)
            - [1]: Geographic latitude in degrees (north positive)
            - [2]: Altitude above sea level in meters
        atmo: Atmospheric conditions as a sequence:
            - [0]: Atmospheric pressure in mbar/hPa
            - [1]: Atmospheric temperature in degrees Celsius
            - [2]: Relative humidity in percent (0-100)
            - [3]: If >= 1: Meteorological Range in km
                   If 0-1: Total atmospheric coefficient (ktot)
                   If 0: Compute ktot from other parameters
        observer: Observer characteristics as a sequence:
            - [0]: Age of observer in years (default 36)
            - [1]: Snellen ratio of observer's eyes (default 1.0 = normal)
            For optical instruments (when HELFLAG_OPTICAL_PARAMS is set):
            - [2]: 0 = monocular, 1 = binocular
            - [3]: Telescope magnification (0 = naked eye)
            - [4]: Optical aperture (telescope diameter) in mm
            - [5]: Optical transmission coefficient
        objname: Name of the object to observe. Can be:
            - Planet name (e.g., "Venus", "Mars", "Jupiter")
            - Fixed star name (e.g., "Sirius", "Aldebaran")
            - Planet number as string (e.g., "2" for Venus)
        flags: Calculation flags combining ephemeris and heliacal flags:
            - SEFLG_SWIEPH, SEFLG_JPLEPH, etc. for ephemeris
            - HELFLAG_OPTICAL_PARAMS: Use optical instrument parameters
            - HELFLAG_NO_DETAILS: Skip detailed calculations
            - HELFLAG_VISLIM_DARK: Assume Sun at nadir (dark sky)
            - HELFLAG_VISLIM_NOMOON: Exclude Moon's brightness contribution

    Returns:
        Tuple containing:
            - result: Visibility status code:
                - (-2): Object is below horizon
                - (0): OK, photopic vision (bright conditions)
                - (1): OK, scotopic vision (dark conditions)
                - (2): OK, near limit between photopic/scotopic
            - dret: Tuple of 8 floats with observation details:
                - [0]: Limiting visual magnitude (object visible if mag < this)
                - [1]: Altitude of object in degrees
                - [2]: Azimuth of object in degrees
                - [3]: Altitude of Sun in degrees
                - [4]: Azimuth of Sun in degrees
                - [5]: Altitude of Moon in degrees
                - [6]: Azimuth of Moon in degrees
                - [7]: Magnitude of object

    Raises:
        ValueError: If objname is empty or invalid

    Algorithm:
        The limiting magnitude calculation is based on Schaefer's model (1990)
        which considers:
        1. Sky brightness from Sun, Moon, zodiacal light, airglow
        2. Atmospheric extinction based on airmass and conditions
        3. Observer's eye adaptation (scotopic vs photopic)
        4. Optional optical instrument characteristics

        The sky background brightness varies with:
        - Sun altitude (twilight contribution)
        - Moon altitude and phase (moonlight)
        - Atmospheric scattering

    Example:
        >>> from libephemeris import julday, vis_limit_mag
        >>> jd = julday(2024, 8, 15, 22.0)
        >>> # Rome location
        >>> geopos = (12.5, 41.9, 0)
        >>> # Standard atmosphere
        >>> atmo = (1013.25, 15.0, 50.0, 0.0)
        >>> # Normal observer
        >>> observer = (36, 1.0)
        >>> result, dret = vis_limit_mag(jd, geopos, atmo, observer, "Venus")
        >>> if dret[0] > dret[7]:
        ...     print("Venus is visible")
        ... else:
        ...     print("Venus is not visible")

    References:
        - Reference API: swe_vis_limit_mag()
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
        - Schaefer, B.E. (1993) "Astronomy and the Limits of Vision"
    """
    # --- LEB fast path ---
    from .state import get_leb_reader as _get_leb_reader

    _leb_rdr = _get_leb_reader()
    if _leb_rdr is not None:
        try:
            return _vis_limit_mag_leb(
                _leb_rdr, tjdut, geopos, atmo, observer, objname, flags,
            )
        except KeyError:
            pass  # Body not in LEB file
        except ValueError:
            pass  # Date outside LEB range
    # --- END LEB fast path ---

    from .planets import _PLANET_MAP, swe_pheno_ut
    from .fixed_stars import swe_fixstar2_ut, swe_fixstar2_mag
    from .state import get_planets, get_timescale
    from .constants import (
        SE_SUN,
        SE_MOON,
        SE_HELFLAG_VISLIM_DARK,
        SE_HELFLAG_VISLIM_NOMOON,
        SE_HELFLAG_BELOW_HORIZON,
        SE_HELFLAG_PHOTOPIC,
        SE_HELFLAG_SCOTOPIC,
        SE_HELFLAG_MIXED,
    )
    from skyfield.api import wgs84

    if not objname:
        raise ValueError("objname cannot be empty")

    # Parse geographic position
    lon = geopos[0] if len(geopos) > 0 else 0.0
    lat = geopos[1] if len(geopos) > 1 else 0.0
    alt_m = geopos[2] if len(geopos) > 2 else 0.0

    # Parse atmospheric conditions
    pressure = atmo[0] if len(atmo) > 0 else 1013.25
    temperature = atmo[1] if len(atmo) > 1 else 15.0
    humidity_pct = atmo[2] if len(atmo) > 2 else 50.0
    met_range = atmo[3] if len(atmo) > 3 else 0.0

    # Parse observer data
    observer_age = observer[0] if len(observer) > 0 else 36.0
    snellen_ratio = observer[1] if len(observer) > 1 else 1.0

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
            "SUN": SE_SUN,
            "MOON": SE_MOON,
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

    # Calculate object position and magnitude
    obj_alt = 0.0
    obj_az = 0.0
    obj_mag = 0.0

    if is_fixed_star:
        # Fixed star calculation
        try:
            star_result, star_name_out, retflag = swe_fixstar2_ut(
                objname, tjdut, flags & 0xFF
            )

            # star_result is (lon, lat, dist, lon_speed, lat_speed, dist_speed)
            # We need to convert ecliptic to horizontal
            star_lon = star_result[0]
            star_lat = star_result[1]

            # Get star magnitude
            try:
                star_mag_val, _star_name_mag = swe_fixstar2_mag(objname)
                obj_mag = star_mag_val
            except (ValueError, TypeError, ArithmeticError):
                obj_mag = 2.0  # Default magnitude if not found

            # Convert ecliptic to equatorial then to horizontal
            # Simplified: use azalt function if available
            from .utils import azalt, SE_ECL2HOR

            hor_result = azalt(
                tjdut,
                SE_ECL2HOR,
                (lon, lat, alt_m),
                pressure,
                temperature,
                (star_lon, star_lat, 1.0),
            )
            obj_az = hor_result[0]
            obj_alt = hor_result[1]

        except ValueError:
            raise
        except (KeyError, IndexError, OSError) as e:
            # Star not found or other error
            raise ValueError(f"could not find star name {objname.lower()}: {e}") from e
    else:
        # Planet calculation
        if body_id is None:
            raise ValueError(f"Unknown object: {objname}")

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
                pheno_result = swe_pheno_ut(tjdut, body_id, flags)
                obj_mag = pheno_result[4]  # Visual magnitude
            except (ValueError, TypeError, ArithmeticError):
                obj_mag = 0.0  # Default bright
        else:
            raise ValueError(f"illegal planet number {body_id}.")

    # Check if object is below horizon
    if obj_alt < 0:
        # Match reference API: return all zeros in data when below horizon
        dret = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        return float(SE_HELFLAG_BELOW_HORIZON), dret

    # Apply HELFLAG options
    use_dark_sky = bool(flags & SE_HELFLAG_VISLIM_DARK)
    exclude_moon = bool(flags & SE_HELFLAG_VISLIM_NOMOON)

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
    )

    # Calculate Moon phase for sky brightness
    moon_phase = 0.0
    moon_obj_angle = 180.0
    sun_obj_angle = 180.0
    try:
        moon_pheno = swe_pheno_ut(tjdut, SE_MOON, flags & 0xFF)
        phase_angle = moon_pheno[0]
        moon_phase = (1.0 - math.cos(math.radians(phase_angle))) / 2.0

        # Calculate angular separation between object and Moon
        moon_geo = earth.at(t).observe(moon).apparent()
        if is_fixed_star:
            # For stars, use the computed position
            moon_obj_angle = 90.0  # Default to 90 degrees
        else:
            body_app_geo = observer_at.at(t).observe(target).apparent()
            moon_obj_angle = body_app_geo.separation_from(moon_app).degrees

        # Calculate elongation from Sun
        sun_geo = earth.at(t).observe(sun).apparent()
        if is_fixed_star:
            sun_obj_angle = 90.0  # Default
        else:
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

    # Apply extinction to object magnitude
    apparent_obj_mag = obj_mag + schaefer.extinction(obj_alt)

    # Determine vision type based on sky brightness
    if sun_alt >= -6:
        vision_type = SE_HELFLAG_PHOTOPIC
    elif sun_alt >= -12:
        vision_type = SE_HELFLAG_MIXED
    else:
        vision_type = SE_HELFLAG_SCOTOPIC

    # Build result tuple (10 elements to match reference API)
    dret = (
        limiting_mag,  # 0: Limiting visual magnitude
        obj_alt,  # 1: Altitude of object
        obj_az,  # 2: Azimuth of object
        sun_alt,  # 3: Altitude of Sun
        sun_az,  # 4: Azimuth of Sun
        moon_alt,  # 5: Altitude of Moon
        moon_az,  # 6: Azimuth of Moon
        apparent_obj_mag,  # 7: Apparent magnitude of object (with extinction)
        0.0,  # 8: Reserved
        0.0,  # 9: Reserved
    )

    return float(vision_type), dret


# Alias for reference API compatibility
swe_vis_limit_mag = vis_limit_mag
