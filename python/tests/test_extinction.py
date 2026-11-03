"""Atmospheric models expose native mutable scalar records and physical limits."""

import inspect

import libephemeris as eph
import pytest


def test_airmass_models_and_horizon_limits():
    for method in ("secant", "kasten_young", "rozenberg"):
        assert eph.calc_airmass(90.0, method) == pytest.approx(1.0, abs=0.001)
        assert eph.calc_airmass(0.0, method) == 40.0
        assert eph.calc_airmass(10.0, method) > eph.calc_airmass(45.0, method)
    with pytest.raises(eph.InputValidationError):
        eph.calc_airmass(45.0, "unknown")


def test_extinction_components_and_plain_records():
    result = eph.calc_extinction_coefficient()
    assert isinstance(result, eph.ExtinctionCoefficients)
    assert result.k_total == result.k_rayleigh + result.k_aerosol + result.k_ozone + result.k_water
    assert eph.calc_extinction_magnitude(45.0) == result.k_total * eph.calc_airmass(45.0)
    assert eph.calc_extinction_magnitude(-1.0) == 99.0
    assert eph.apparent_magnitude_with_extinction(2.0, 45.0) == 2.0 + eph.calc_extinction_magnitude(45.0)
    record = eph.ExtinctionCoefficients(1.0, 2.0, 3.0, 4.0, 10.0)
    record.k_total = 12.0
    assert record.k_total == 12.0
    assert list(inspect.signature(eph.ExtinctionCoefficients).parameters) == ["k_rayleigh", "k_aerosol", "k_ozone", "k_water", "k_total"]


def test_twilight_phases_and_record_units():
    assert [eph.get_twilight_phase(value) for value in (1.0, -3.0, -9.0, -15.0, -20.0)] == ["day", "civil", "nautical", "astronomical", "night"]
    result = eph.calc_twilight_sky_brightness(-12.0)
    assert isinstance(result, eph.TwilightSkyBrightness)
    assert result.twilight_phase == "astronomical"
    assert 0 <= result.surface_brightness <= 22.5
    assert result.nanolamberts > 0
    assert result.sun_altitude_deg == -12.0
    assert eph.calc_twilight_brightness_simple(-12.0) == result.surface_brightness


def test_visibility_threshold_and_eye_adaptation():
    assert [eph.calc_eye_adaptation_state(value) for value in (10.0, 18.0, 21.5)] == ["photopic", "mesopic", "dark"]
    result = eph.calc_visibility_threshold(2.0, 21.5, observer_skill=2)
    assert isinstance(result, eph.VisibilityResult)
    assert result.is_visible == (result.object_magnitude < result.limiting_magnitude)
    assert result.visibility_margin == result.limiting_magnitude - result.object_magnitude
    assert result.is_visible == eph.is_object_visible(2.0, 21.5)
    assert result.limiting_magnitude == eph.calc_limiting_magnitude_for_sky(21.5)
    assert eph.calc_contrast_threshold(21.5, "scotopic") == eph.calc_contrast_threshold(21.5, "dark")
    with pytest.raises(eph.InputValidationError):
        eph.calc_contrast_threshold(21.5, "unknown")
