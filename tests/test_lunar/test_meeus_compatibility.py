# SPDX-License-Identifier: Apache-2.0
"""Compatibility checks for historical lunar warning symbols."""

from __future__ import annotations

import warnings

from libephemeris.lunar import (
    MEEUS_MAX_CENTURIES,
    MEEUS_OPTIMAL_CENTURIES,
    MEEUS_VALID_CENTURIES,
    MeeusPolynomialWarning,
    MeeusRangeError,
    calc_mean_lilith,
    calc_mean_lunar_node,
)


def test_historical_compatibility_symbols_remain_exported() -> None:
    assert MEEUS_OPTIMAL_CENTURIES == 2.0
    assert MEEUS_VALID_CENTURIES == 10.0
    assert MEEUS_MAX_CENTURIES == 20.0
    assert issubclass(MeeusPolynomialWarning, UserWarning)
    assert issubclass(MeeusRangeError, ValueError)
    assert str(MeeusRangeError("outside range")) == "outside range"


def test_clean_room_mean_model_does_not_emit_legacy_polynomial_warnings() -> None:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        calc_mean_lunar_node(2451545.0)
        calc_mean_lilith(2451545.0)

    assert not any(issubclass(item.category, MeeusPolynomialWarning) for item in caught)
