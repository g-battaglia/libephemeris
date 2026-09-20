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
