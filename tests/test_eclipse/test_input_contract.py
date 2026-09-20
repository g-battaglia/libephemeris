# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Public observer-position contract for eclipse and rise/set entry points."""

from __future__ import annotations

from typing import Any

import pytest

import libephemeris as ephem
from libephemeris import eclipse as eclipse_module
from libephemeris.exceptions import CoordinateError, InputValidationError


_JD = 2451545.0
_ENTRY_POINTS = (
    "sol_eclipse_when_loc",
    "sol_eclipse_how",
    "sol_eclipse_how_details",
    "lun_eclipse_when_loc",
    "lun_eclipse_how",
    "lun_occult_when_loc",
    "rise_trans",
    "rise_trans_true_hor",
    "sol_eclipse_magnitude_at_loc",
    "sol_eclipse_obscuration_at_loc",
)


def _call(name: str, geopos: Any):
    """Call one public entry point with the same observer argument."""
    function = getattr(ephem, name)
    if name == "lun_occult_when_loc":
        return function(_JD, ephem.SUN, geopos)
    if name in {"rise_trans", "rise_trans_true_hor"}:
        return function(_JD, ephem.SUN, ephem.CALC_RISE, geopos)
    return function(_JD, geopos)


@pytest.mark.parametrize("name", _ENTRY_POINTS)
@pytest.mark.parametrize(
    "geopos,error_type",
    [
        ((), ValueError),
        ((0.0, 0.0), ValueError),
        (None, ValueError),
        (("0", "0", "0"), InputValidationError),
        ((True, 0.0, 0.0), InputValidationError),
        ((None, 0.0, 0.0), InputValidationError),
        ((0.0, 91.0, 0.0), CoordinateError),
        ((0.0, float("nan"), 0.0), CoordinateError),
        ((400.0, 0.0, 0.0), CoordinateError),
        ((0.0, 0.0, float("inf")), InputValidationError),
    ],
    ids=(
        "empty",
        "missing-altitude",
        "not-a-sequence",
        "numeric-text",
        "boolean",
        "nonnumeric",
        "latitude-out-of-range",
        "latitude-nan",
        "longitude-out-of-range",
        "altitude-infinite",
    ),
)
def test_invalid_geopos_is_rejected_before_backend(
    monkeypatch: pytest.MonkeyPatch,
    name: str,
    geopos: Any,
    error_type: type[Exception],
) -> None:
    """Every invalid location fails before a calculation can start."""

    def unexpected_backend(*_args, **_kwargs):
        raise AssertionError("backend reached with an invalid observer position")

    monkeypatch.setattr(
        eclipse_module, "_call_with_leb_skyfield_fallback", unexpected_backend
    )
    with pytest.raises(error_type):
        _call(name, geopos)


def test_valid_geopos_is_normalized_to_native_floats(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Padding is ignored and valid numeric coordinates retain their values."""
    captured = None

    def capture(_impl, _jd, geopos, _flags):
        nonlocal captured
        captured = geopos
        return 0, (0.0,) * 20

    monkeypatch.setattr(eclipse_module, "_call_with_leb_skyfield_fallback", capture)
    ephem.sol_eclipse_how(_JD, [12, 45, 100, 999])
    assert isinstance(captured, tuple)
    assert captured == (12.0, 45.0, 100.0)
    assert all(type(value) is float for value in captured)


@pytest.mark.parametrize(
    "name,component,value,error_type",
    [
        ("sol_eclipse_max_time", "lat", 91.0, CoordinateError),
        ("sol_eclipse_max_time", "lon", 400.0, CoordinateError),
        ("sol_eclipse_max_time", "lat", float("nan"), CoordinateError),
        ("sol_eclipse_max_time", "altitude", float("inf"), InputValidationError),
        ("sol_eclipse_max_time", "lat", "41.9", InputValidationError),
        ("calc_eclipse_path_width", "lat", -91.0, CoordinateError),
        ("calc_eclipse_path_width", "lon", float("nan"), CoordinateError),
        ("calc_eclipse_path_width", "lon", 400.0, CoordinateError),
        ("calc_eclipse_path_width", "lat", True, InputValidationError),
        ("planet_occult_when_loc", "lat", 91.0, CoordinateError),
        ("planet_occult_when_loc", "lon", -181.0, CoordinateError),
        ("planet_occult_when_loc", "lon", float("inf"), CoordinateError),
        ("planet_occult_when_loc", "altitude", float("nan"), InputValidationError),
        ("planet_occult_when_loc", "lat", None, InputValidationError),
    ],
)
def test_invalid_scalar_observer_coordinate_is_rejected_before_calculation(
    monkeypatch: pytest.MonkeyPatch,
    name: str,
    component: str,
    value: Any,
    error_type: type[Exception],
) -> None:
    """Scalar observer coordinates use the same physical domain as geopos."""

    def unexpected_backend(*_args, **_kwargs):
        raise AssertionError("backend reached with an invalid observer position")

    monkeypatch.setattr(eclipse_module, "_get_leb_reader_safe", unexpected_backend)
    monkeypatch.setattr(
        eclipse_module, "_calc_local_eclipse_max_time", unexpected_backend
    )
    monkeypatch.setattr(
        eclipse_module, "_calc_eclipse_path_width_impl", unexpected_backend
    )
    monkeypatch.setattr(
        eclipse_module, "_calc_global_eclipse_max_time", unexpected_backend
    )
    monkeypatch.setattr("libephemeris.state.get_planets", unexpected_backend)

    coords = {"lat": 41.9, "lon": 12.5, "altitude": 0.0}
    coords[component] = value
    with pytest.raises(error_type):
        if name == "sol_eclipse_max_time":
            ephem.sol_eclipse_max_time(_JD, **coords)
        elif name == "calc_eclipse_path_width":
            ephem.calc_eclipse_path_width(_JD, lat=coords["lat"], lon=coords["lon"])
        else:
            ephem.planet_occult_when_loc(_JD, ephem.VENUS, ephem.MARS, "", **coords)


def test_scalar_observer_coordinates_are_native_floats(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Accepted integer coordinates reach local solvers as Python floats."""
    captured = None

    def capture_max(_jd, lat, lon, altitude, _search_range):
        nonlocal captured
        captured = (lat, lon, altitude)
        return _JD, 0.0

    monkeypatch.setattr(eclipse_module, "_calc_local_eclipse_max_time", capture_max)
    ephem.sol_eclipse_max_time(_JD, 42, 12, 100)
    assert captured == (42.0, 12.0, 100.0)
    assert all(type(part) is float for part in captured)

    def capture_width(_impl, _jd, lat, lon, _flags):
        nonlocal captured
        captured = (lat, lon)
        return 0.0

    monkeypatch.setattr(
        eclipse_module, "_call_with_leb_skyfield_fallback", capture_width
    )
    ephem.calc_eclipse_path_width(_JD, 42, 12)
    assert captured == (42.0, 12.0)
    ephem.calc_eclipse_path_width(_JD, 42)
    assert captured == (42.0, None)


def test_global_maximum_ignores_unused_altitude(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Without a site, altitude remains outside the global calculation."""
    monkeypatch.setattr(
        eclipse_module, "_calc_global_eclipse_max_time", lambda *_args: (_JD, 0.0)
    )
    assert ephem.sol_eclipse_max_time(_JD, altitude="unused") == (_JD, 0.0)


def test_planet_body_refusal_precedes_location_refusal(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The established illegal-body error remains the first public refusal."""

    def unexpected_reader():
        raise AssertionError("reader opened before invalid input was refused")

    monkeypatch.setattr(eclipse_module, "_get_leb_reader_safe", unexpected_reader)
    with pytest.raises(ValueError, match="The Sun cannot be the occulting body"):
        ephem.planet_occult_when_loc(_JD, ephem.SUN, ephem.MARS, lat=float("nan"))
