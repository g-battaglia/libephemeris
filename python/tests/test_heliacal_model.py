"""Typed VISLIMIT models retain the pinned public optical calculation methods."""

import inspect
import json
import math
import subprocess
from pathlib import Path

import pytest

CONFIGURATIONS = [
    {},
    {
        "pressure": 800,
        "temperature": 5,
        "humidity": 70,
        "altitude": 2000,
        "latitude": 42,
        "jd": 2460310.5,
        "observer_age": 50,
        "snellen": 0.8,
        "binocular": True,
        "telescope_mag": 10,
        "aperture": 70,
        "transmission": 0.85,
    },
    {
        "pressure": 800,
        "temperature": 5,
        "humidity": 70,
        "altitude": 2000,
        "latitude": 42,
        "jd": 2460310.5,
        "pressure_scales_extinction": True,
    },
    {"met_range": 0.2, "latitude": -35, "jd": 2451545.0},
    {"met_range": 20, "latitude": -35, "jd": 2451545.0},
]
METHOD_CALLS = [
    *[("airmass", [altitude]) for altitude in (-5, 0, 15, 80, 90)],
    *[("extinction", [altitude]) for altitude in (0, 15, 80)],
    ("sky_brightness_bl", [-18, -90, 0, 45, 180, 180]),
    ("sky_brightness_bl", [-9, 30, 0.5, 20, 60, 50]),
    ("limiting_magnitude", [-18]),
    ("limiting_magnitude", [-9, 30, 0.5, 20, 60, 50]),
    ("arcus_visionis_required", [-2]),
    ("arcus_visionis_required", [2, "star"]),
    ("is_visible", [-1, 2, -12, 30]),
    ("is_visible", [20, 2, -12, 30]),
    ("is_visible", [20, 2, -12, 30, 30, 0.5, 50, 0.2]),
    ("heliacal_visibility_angle", [-2]),
    ("heliacal_visibility_angle", [2, -10]),
]
COEFFICIENTS = ("k_rayleigh", "k_aerosol", "k_ozone", "k_water", "k_total")


@pytest.mark.parametrize("configuration", CONFIGURATIONS)
def test_season_updates_match_fresh_models_at_the_actual_requested_date(configuration):
    from libephemeris import heliacal as h

    model = h.create_schaefer_model(**configuration)
    start = configuration.get("jd", 2451545.0)
    # Include forward, backward, sub-day and repeated requests. Identical
    # atmosphere/observer/date inputs must not depend on sampling chronology.
    for step in (1.0, 0.25, -0.5, 2.0, 4.0, 0.25, 0.25):
        epoch = start + step
        model.update_season(epoch)
        fresh = h.create_schaefer_model(**{**configuration, "jd": epoch})
        assert model.jd == fresh.jd == epoch
        for name in COEFFICIENTS:
            assert getattr(model, name) == getattr(fresh, name)
        for name, arguments in METHOD_CALLS:
            assert getattr(model, name)(*arguments) == getattr(fresh, name)(*arguments)


def test_typed_record_constructor_and_factory_use_the_same_native_model():
    from libephemeris import heliacal as h

    atmosphere = h.AtmosphericConditions(800, 5, 70, 0, 2000)
    observer = h.ObserverParams(50, 0.8, True, 10, 70, 0.85)
    model = h.SchaeferModel(atmosphere, observer, latitude=42, jd=2460310.5)
    factory = h.create_schaefer_model(**CONFIGURATIONS[1])
    assert model.atmo == atmosphere and model.observer == observer
    for name in COEFFICIENTS:
        assert getattr(model, name) == getattr(factory, name)
    assert model.limiting_magnitude(-12) == factory.limiting_magnitude(-12)
    with pytest.raises(ValueError):
        h.create_schaefer_model(temperature=float("nan"))
