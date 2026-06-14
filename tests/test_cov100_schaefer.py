# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for ``libephemeris.schaefer``.

The Schaefer (1990) atmospheric visibility model is pure-math with no
external dependencies, so these tests exercise every public function and
every conditional branch using deterministic numeric inputs.
"""

from __future__ import annotations

import math

import pytest

from libephemeris import schaefer


# ---------------------------------------------------------------------------
# Extinction
# ---------------------------------------------------------------------------


def test_calc_rayleigh_extinction() -> None:
    """Rayleigh extinction scales with pressure and altitude."""
    sea = schaefer.calc_rayleigh_extinction(1013.25, 0.0)
    assert sea == pytest.approx(schaefer.RAYLEIGH_COEFF_SEA_LEVEL)
    high = schaefer.calc_rayleigh_extinction(1013.25, 8500.0)
    assert high < sea


def test_calc_aerosol_extinction_met_range_ge_one() -> None:
    """met_range_km >= 1.0 uses the visual-range formula."""
    k = schaefer.calc_aerosol_extinction(0.5, met_range_km=25.0)
    assert k == pytest.approx(max(0.0, 3.912 / 25.0 - 0.106))


def test_calc_aerosol_extinction_met_range_ge_one_clamped() -> None:
    """Very large met range clamps the aerosol coefficient at zero."""
    k = schaefer.calc_aerosol_extinction(0.5, met_range_km=1000.0)
    assert k == 0.0


def test_calc_aerosol_extinction_direct_coeff() -> None:
    """0 < met_range_km < 1.0 is treated as a direct coefficient."""
    k = schaefer.calc_aerosol_extinction(0.5, met_range_km=0.3, altitude_m=0.0)
    assert k == pytest.approx(0.3)


def test_calc_aerosol_extinction_from_humidity() -> None:
    """met_range_km == 0 estimates the coefficient from humidity."""
    k = schaefer.calc_aerosol_extinction(1.0, met_range_km=0.0, altitude_m=0.0)
    assert k == pytest.approx(0.04 + 0.18 + 0.12)


def test_calc_aerosol_extinction_altitude_reduces() -> None:
    """Higher altitude reduces aerosol extinction."""
    low = schaefer.calc_aerosol_extinction(0.5, met_range_km=0.3, altitude_m=0.0)
    high = schaefer.calc_aerosol_extinction(0.5, met_range_km=0.3, altitude_m=1500.0)
    assert high < low


def test_calc_ozone_extinction() -> None:
    """Ozone coefficient is the module constant."""
    assert schaefer.calc_ozone_extinction() == schaefer.OZONE_COEFF


def test_calc_water_vapor_extinction() -> None:
    """Water vapor scales linearly with humidity."""
    assert schaefer.calc_water_vapor_extinction(0.5) == pytest.approx(
        schaefer.WATER_VAPOR_COEFF * 0.5
    )


def test_calc_total_extinction() -> None:
    """Total extinction is the sum of the four components."""
    k = schaefer.calc_total_extinction()
    expected = (
        schaefer.calc_rayleigh_extinction(1013.25, 0.0)
        + schaefer.calc_aerosol_extinction(0.5, 0.0, 0.0)
        + schaefer.calc_ozone_extinction()
        + schaefer.calc_water_vapor_extinction(0.5)
    )
    assert k == pytest.approx(expected)


# ---------------------------------------------------------------------------
# Airmass
# ---------------------------------------------------------------------------


def test_calc_airmass_below_horizon() -> None:
    """At/below horizon returns the maximum airmass."""
    assert schaefer.calc_airmass(0.0) == 40.0
    assert schaefer.calc_airmass(-5.0) == 40.0


def test_calc_airmass_zenith() -> None:
    """At/above the zenith returns 1.0."""
    assert schaefer.calc_airmass(90.0) == 1.0
    assert schaefer.calc_airmass(95.0) == 1.0


def test_calc_airmass_normal() -> None:
    """A mid-altitude airmass falls between 1 and 40."""
    am = schaefer.calc_airmass(45.0)
    assert 1.0 < am < 40.0


def test_calc_airmass_near_horizon() -> None:
    """Very low (positive) altitude yields a large but bounded airmass."""
    am = schaefer.calc_airmass(0.001)
    assert 1.0 < am <= 40.0


# ---------------------------------------------------------------------------
# Extinction magnitude
# ---------------------------------------------------------------------------


def test_calc_extinction_magnitude() -> None:
    """Extinction magnitude equals coefficient times airmass."""
    mag = schaefer.calc_extinction_magnitude(45.0, 0.3)
    assert mag == pytest.approx(0.3 * schaefer.calc_airmass(45.0))


# ---------------------------------------------------------------------------
# Twilight brightness
# ---------------------------------------------------------------------------


def test_calc_twilight_brightness_daytime() -> None:
    """Sun at or above the horizon is effectively infinite brightness."""
    assert schaefer.calc_twilight_brightness(0.0) == 1e10
    assert schaefer.calc_twilight_brightness(5.0) == 1e10


def test_calc_twilight_brightness_night() -> None:
    """Astronomical night has no twilight contribution."""
    assert schaefer.calc_twilight_brightness(-18.0) == 0.0
    assert schaefer.calc_twilight_brightness(-25.0) == 0.0


def test_calc_twilight_brightness_civil() -> None:
    """Civil twilight (h <= 6) uses the steep exponential."""
    b = schaefer.calc_twilight_brightness(-3.0)
    assert b == pytest.approx(1e8 * math.exp(-3.0 * 1.5))


def test_calc_twilight_brightness_nautical() -> None:
    """Nautical twilight (6 < h <= 12) uses the moderate exponential."""
    b = schaefer.calc_twilight_brightness(-9.0)
    assert b == pytest.approx(5e4 * math.exp(-3.0 * 0.6))


def test_calc_twilight_brightness_astronomical() -> None:
    """Astronomical twilight (12 < h < 18) uses the shallow exponential."""
    b = schaefer.calc_twilight_brightness(-15.0)
    assert b == pytest.approx(1e3 * math.exp(-3.0 * 0.4))


# ---------------------------------------------------------------------------
# Moon brightness
# ---------------------------------------------------------------------------


def test_calc_moon_brightness_below_horizon() -> None:
    """Moon below the horizon contributes nothing."""
    assert schaefer.calc_moon_brightness(-5.0, 0.0, 45.0) == 0.0


def test_calc_moon_brightness_new_moon() -> None:
    """Near new moon (illuminated fraction < 0.01) contributes nothing."""
    assert schaefer.calc_moon_brightness(45.0, 180.0, 45.0) == 0.0


def test_calc_moon_brightness_close_distance() -> None:
    """angular_distance < 1 is clamped (no division by zero) and < 10 path."""
    b = schaefer.calc_moon_brightness(45.0, 0.0, 0.0)
    assert b > 0.0


def test_calc_moon_brightness_mid_distance() -> None:
    """10 <= angular_distance < 45 uses the linear falloff branch."""
    b = schaefer.calc_moon_brightness(45.0, 0.0, 20.0)
    assert b > 0.0


def test_calc_moon_brightness_far_distance() -> None:
    """angular_distance >= 45 uses the flat distance factor."""
    b = schaefer.calc_moon_brightness(45.0, 0.0, 90.0)
    moon_airmass = schaefer.calc_airmass(45.0)
    alt_factor = math.exp(-0.4 * (moon_airmass - 1.0))
    expected = 3000.0 * 1.0 * alt_factor * 1.0
    assert b == pytest.approx(expected)


# ---------------------------------------------------------------------------
# Airglow / zodiacal brightness
# ---------------------------------------------------------------------------


def test_calc_airglow_brightness() -> None:
    """Airglow returns the module constant."""
    assert schaefer.calc_airglow_brightness() == schaefer.AIRGLOW_BRIGHTNESS_NL


def test_calc_zodiacal_brightness_too_close_to_sun() -> None:
    """Elongation < 20 is overwhelmed by twilight (zero contribution)."""
    assert schaefer.calc_zodiacal_brightness(0.0, 10.0) == 0.0


def test_calc_zodiacal_brightness_peak() -> None:
    """20 <= elongation < 90 is the peak region."""
    b = schaefer.calc_zodiacal_brightness(0.0, 45.0)
    expected = 150.0 * math.exp(0.0) * math.sin(math.radians(45.0))
    assert b == pytest.approx(expected)


def test_calc_zodiacal_brightness_anti_sun() -> None:
    """Elongation >= 90 decreases toward the anti-Sun point."""
    b = schaefer.calc_zodiacal_brightness(0.0, 120.0)
    expected = 150.0 * math.exp(0.0) * math.sin(math.radians(60.0))
    assert b == pytest.approx(expected)


def test_calc_zodiacal_brightness_latitude_factor() -> None:
    """Higher ecliptic latitude dims the zodiacal light."""
    near = schaefer.calc_zodiacal_brightness(0.0, 45.0)
    far = schaefer.calc_zodiacal_brightness(30.0, 45.0)
    assert far < near


# ---------------------------------------------------------------------------
# Total sky brightness
# ---------------------------------------------------------------------------


def test_calc_total_sky_brightness_with_zodiacal() -> None:
    """Sun < -12 includes the zodiacal contribution."""
    b = schaefer.calc_total_sky_brightness(
        sun_altitude_deg=-15.0,
        moon_altitude_deg=-90.0,
        moon_phase_deg=180.0,
        moon_distance_deg=90.0,
        ecliptic_latitude_deg=0.0,
        sun_elongation_deg=45.0,
    )
    expected = (
        schaefer.calc_twilight_brightness(-15.0)
        + 0.0
        + schaefer.calc_airglow_brightness()
        + schaefer.calc_zodiacal_brightness(0.0, 45.0)
    )
    assert b == pytest.approx(expected)


def test_calc_total_sky_brightness_without_zodiacal() -> None:
    """Sun >= -12 excludes the zodiacal contribution."""
    b = schaefer.calc_total_sky_brightness(
        sun_altitude_deg=-6.0,
        sun_elongation_deg=45.0,
    )
    expected = (
        schaefer.calc_twilight_brightness(-6.0)
        + 0.0
        + schaefer.calc_airglow_brightness()
        + 0.0
    )
    assert b == pytest.approx(expected)


# ---------------------------------------------------------------------------
# Brightness <-> mag/arcsec^2
# ---------------------------------------------------------------------------


def test_brightness_to_mag_arcsec2_dark() -> None:
    """Non-positive brightness returns the dark-sky limit."""
    assert schaefer.brightness_to_mag_arcsec2(0.0) == 22.0
    assert schaefer.brightness_to_mag_arcsec2(-5.0) == 22.0


def test_brightness_to_mag_arcsec2_positive() -> None:
    """Positive brightness uses the Garstang formula."""
    m = schaefer.brightness_to_mag_arcsec2(250.0)
    assert m == pytest.approx(26.33 - 2.5 * math.log10(250.0))


def test_mag_arcsec2_to_brightness_roundtrip() -> None:
    """mag/arcsec^2 -> brightness inverts the positive conversion."""
    b = schaefer.mag_arcsec2_to_brightness(21.0)
    assert b == pytest.approx(10 ** ((26.33 - 21.0) / 2.5))
    assert schaefer.brightness_to_mag_arcsec2(b) == pytest.approx(21.0)


# ---------------------------------------------------------------------------
# Limiting magnitude
# ---------------------------------------------------------------------------


def test_calc_limiting_magnitude_dark_young_normal() -> None:
    """Dark sky, age <= baseline, normal Snellen: only extinction reduces."""
    lim = schaefer.calc_limiting_magnitude(
        sky_brightness_nl=100.0,
        extinction_coeff=0.0,
        object_altitude_deg=90.0,
        observer_age=20.0,
        snellen_ratio=1.0,
    )
    assert lim == pytest.approx(schaefer.IDEAL_NAKED_EYE_LIMIT_MAG)


def test_calc_limiting_magnitude_bright_sky_old_poor_vision() -> None:
    """Bright sky, age > baseline, Snellen != 1 all reduce the limit."""
    lim = schaefer.calc_limiting_magnitude(
        sky_brightness_nl=2500.0,
        extinction_coeff=0.3,
        object_altitude_deg=30.0,
        observer_age=60.0,
        snellen_ratio=0.5,
    )
    base = schaefer.IDEAL_NAKED_EYE_LIMIT_MAG
    sky_factor = 2.5 * math.log10(2500.0 / 250.0)
    age_factor = schaefer.AGE_DEGRADATION_PER_DECADE * ((60.0 - 20.0) / 10.0)
    snellen_factor = 2.5 * math.log10(1.0 / 0.5)
    airmass = schaefer.calc_airmass(30.0)
    extinction_loss = 0.3 * (airmass - 1.0)
    expected = base - sky_factor - age_factor - snellen_factor - extinction_loss
    assert lim == pytest.approx(expected)


def test_calc_limiting_magnitude_snellen_zero() -> None:
    """snellen_ratio <= 0 takes the zero-factor branch."""
    lim = schaefer.calc_limiting_magnitude(
        sky_brightness_nl=100.0,
        extinction_coeff=0.0,
        object_altitude_deg=90.0,
        observer_age=20.0,
        snellen_ratio=0.0,
    )
    assert lim == pytest.approx(schaefer.IDEAL_NAKED_EYE_LIMIT_MAG)


# ---------------------------------------------------------------------------
# is_object_visible
# ---------------------------------------------------------------------------


def test_is_object_visible_below_horizon() -> None:
    """Object below the horizon is never visible."""
    vis, lim, app = schaefer.is_object_visible(
        object_magnitude=0.0,
        object_altitude_deg=-1.0,
        sun_altitude_deg=-15.0,
    )
    assert vis is False
    assert lim == 0.0
    assert app == 0.0


def test_is_object_visible_sun_up() -> None:
    """Sun above the horizon blocks heliacal visibility."""
    vis, lim, app = schaefer.is_object_visible(
        object_magnitude=0.0,
        object_altitude_deg=20.0,
        sun_altitude_deg=1.0,
    )
    assert vis is False
    assert lim == -5.0
    assert app == 0.0


def test_is_object_visible_bright_object_with_moon() -> None:
    """Bright object (mag < -2) at night, include_moon path, visible."""
    vis, lim, app = schaefer.is_object_visible(
        object_magnitude=-4.0,
        object_altitude_deg=20.0,
        sun_altitude_deg=-15.0,
        sun_elongation_deg=40.0,
        moon_altitude_deg=-90.0,
        include_moon=True,
    )
    assert vis is True
    assert lim > app


def test_is_object_visible_normal_object_no_moon() -> None:
    """Normal object (-2 <= mag < 2) with include_moon=False branch."""
    vis, lim, app = schaefer.is_object_visible(
        object_magnitude=0.0,
        object_altitude_deg=30.0,
        sun_altitude_deg=-15.0,
        sun_elongation_deg=40.0,
        include_moon=False,
    )
    assert isinstance(vis, bool)
    assert isinstance(lim, float)
    assert isinstance(app, float)


def test_is_object_visible_faint_object_elongation_fail() -> None:
    """Faint object (mag >= 2) below the minimum elongation fails the gate."""
    vis, _lim, _app = schaefer.is_object_visible(
        object_magnitude=4.0,
        object_altitude_deg=30.0,
        sun_altitude_deg=-15.0,
        sun_elongation_deg=5.0,
    )
    assert vis is False


def test_is_object_visible_bright_object_elongation_fail() -> None:
    """Bright object below MIN_ELONGATION_BRIGHT also fails the gate."""
    vis, _lim, _app = schaefer.is_object_visible(
        object_magnitude=-4.0,
        object_altitude_deg=30.0,
        sun_altitude_deg=-15.0,
        sun_elongation_deg=1.0,
    )
    assert vis is False


# ---------------------------------------------------------------------------
# Arcus visionis
# ---------------------------------------------------------------------------


def test_get_arcus_visionis_table_match() -> None:
    """A magnitude within a table band returns that band's value."""
    assert schaefer.get_arcus_visionis(0.5) == 9.0


def test_get_arcus_visionis_very_faint() -> None:
    """Magnitudes at/above 6.0 return the very-faint value."""
    assert schaefer.get_arcus_visionis(6.5) == 17.0


def test_get_arcus_visionis_extremely_bright() -> None:
    """Magnitudes below the brightest band return the bright fallback."""
    assert schaefer.get_arcus_visionis(-10.0) == 4.0


# ---------------------------------------------------------------------------
# Optimal sun altitude
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("mag", "expected"),
    [
        (-5.0, -8.0),
        (-1.0, -9.0),
        (1.0, -10.0),
        (3.0, -11.0),
        (5.0, -12.0),
    ],
)
def test_get_optimal_sun_altitude(mag: float, expected: float) -> None:
    """Each magnitude band maps to its optimal depression."""
    assert schaefer.get_optimal_sun_altitude(mag, True) == expected


# ---------------------------------------------------------------------------
# Heliacal altitude threshold
# ---------------------------------------------------------------------------


def test_calc_heliacal_altitude_threshold_deep_sun() -> None:
    """Sun < -12 uses the 0.8 factor."""
    th = schaefer.calc_heliacal_altitude_threshold(0.5, -15.0, 0.25)
    base_av = schaefer.get_arcus_visionis(0.5)
    expected = max(1.0, base_av * 1.0 * 0.8)
    assert th == pytest.approx(expected)


def test_calc_heliacal_altitude_threshold_mid_sun() -> None:
    """-12 <= Sun < -9 uses the 1.0 factor."""
    th = schaefer.calc_heliacal_altitude_threshold(0.5, -10.0, 0.25)
    base_av = schaefer.get_arcus_visionis(0.5)
    expected = max(1.0, base_av * 1.0 * 1.0)
    assert th == pytest.approx(expected)


def test_calc_heliacal_altitude_threshold_shallow_sun() -> None:
    """Sun >= -9 uses the linearly increasing factor."""
    th = schaefer.calc_heliacal_altitude_threshold(0.5, -6.0, 0.25)
    base_av = schaefer.get_arcus_visionis(0.5)
    sun_factor = 1.2 + (-6.0 + 9) * 0.1
    expected = max(1.0, base_av * 1.0 * sun_factor)
    assert th == pytest.approx(expected)


def test_calc_heliacal_altitude_threshold_floor() -> None:
    """The result is floored at 1.0 degree."""
    th = schaefer.calc_heliacal_altitude_threshold(-4.0, -15.0, 0.0)
    assert th >= 1.0


# ---------------------------------------------------------------------------
# Convenience dictionary
# ---------------------------------------------------------------------------


def test_get_visibility_conditions_keys_and_values() -> None:
    """The convenience function returns the full set of derived values."""
    conditions = schaefer.get_visibility_conditions(
        sun_altitude_deg=-15.0,
        object_altitude_deg=30.0,
        object_magnitude=1.0,
    )
    expected_keys = {
        "extinction_total",
        "extinction_rayleigh",
        "extinction_aerosol",
        "airmass",
        "extinction_magnitude",
        "sky_brightness_nl",
        "sky_brightness_mag_arcsec2",
        "limiting_magnitude",
        "apparent_magnitude",
        "magnitude_margin",
        "is_visible",
        "arcus_visionis",
        "optimal_sun_altitude",
    }
    assert set(conditions) == expected_keys
    assert isinstance(conditions["is_visible"], bool)
    assert conditions["magnitude_margin"] == pytest.approx(
        conditions["limiting_magnitude"] - conditions["apparent_magnitude"]
    )
