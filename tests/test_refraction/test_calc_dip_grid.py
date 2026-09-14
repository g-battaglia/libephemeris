"""Grid contracts for the public sea-horizon dip calculation."""

from __future__ import annotations

import math

import pytest

from libephemeris.refraction import calc_dip

pytestmark = pytest.mark.unit

HEIGHTS_M = (0.0, 1.0, 10.0, 100.0, 1_000.0, 5_000.0, 10_000.0)
PRESSURES_MBAR = (0.0, 500.0, 900.0, 1013.25, 1050.0)
TEMPERATURES_C = (-40.0, 0.0, 10.0, 40.0)
LAPSE_RATES_K_M = (0.0, 0.003, 0.0065, 0.01)


@pytest.mark.parametrize("pressure", PRESSURES_MBAR)
@pytest.mark.parametrize("temperature", TEMPERATURES_C)
@pytest.mark.parametrize("lapse_rate", LAPSE_RATES_K_M)
def test_height_grid_is_finite_nonpositive_and_monotone(
    pressure: float, temperature: float, lapse_rate: float
) -> None:
    """Ordinary atmospheres produce finite dips whose magnitude grows with height."""
    values = tuple(
        calc_dip(height, lapse_rate, pressure, temperature) for height in HEIGHTS_M
    )
    assert values[0] == 0.0
    assert math.copysign(1.0, values[0]) == 1.0
    assert all(type(value) is float and math.isfinite(value) for value in values)
    assert all(value <= 0.0 for value in values)
    magnitudes = tuple(-value for value in values)
    assert all(after >= before for before, after in zip(magnitudes, magnitudes[1:]))


@pytest.mark.parametrize("height", HEIGHTS_M[1:])
@pytest.mark.parametrize("temperature", TEMPERATURES_C)
@pytest.mark.parametrize("lapse_rate", LAPSE_RATES_K_M)
def test_zero_pressure_grid_is_exact_spherical_geometry(
    height: float, temperature: float, lapse_rate: float
) -> None:
    """With no atmosphere, temperature and lapse rate cannot alter geometry."""
    earth_radius_m = 6_378_136.6  # IERS Conventions 2010, Table 1.1.
    expected = -math.degrees(math.acos(earth_radius_m / (earth_radius_m + height)))
    assert calc_dip(height, lapse_rate, 0.0, temperature) == expected


@pytest.mark.parametrize("height", (10.0, 100.0, 1_000.0, 10_000.0))
@pytest.mark.parametrize("temperature", TEMPERATURES_C)
@pytest.mark.parametrize("lapse_rate", LAPSE_RATES_K_M)
def test_pressure_grid_bends_the_visible_horizon_upward(
    height: float, temperature: float, lapse_rate: float
) -> None:
    """For nonnegative pressure, refraction cannot deepen the geometric dip."""
    geometric = abs(calc_dip(height, lapse_rate, 0.0, temperature))
    magnitudes = tuple(
        abs(calc_dip(height, lapse_rate, pressure, temperature))
        for pressure in PRESSURES_MBAR
    )
    assert magnitudes[0] == geometric
    assert all(value <= geometric for value in magnitudes)
    assert all(after <= before for before, after in zip(magnitudes, magnitudes[1:]))


@pytest.mark.parametrize("height", (10.0, 1_000.0, 10_000.0))
@pytest.mark.parametrize("pressure", PRESSURES_MBAR[1:])
def test_temperature_and_lapse_rate_grid_follow_density_gradient(
    height: float, pressure: float
) -> None:
    """Warmer air deepens the dip, while a steeper lapse rate raises the horizon."""
    by_temperature = tuple(
        abs(calc_dip(height, 0.0065, pressure, temperature))
        for temperature in TEMPERATURES_C
    )
    assert all(
        after >= before for before, after in zip(by_temperature, by_temperature[1:])
    )

    by_lapse = tuple(
        abs(calc_dip(height, lapse_rate, pressure, 10.0))
        for lapse_rate in LAPSE_RATES_K_M
    )
    assert all(after <= before for before, after in zip(by_lapse, by_lapse[1:]))


def test_grid_covers_the_declared_axes() -> None:
    """The regression matrix retains each independent atmospheric dimension."""
    assert len(HEIGHTS_M) == 7
    assert len(PRESSURES_MBAR) == 5
    assert len(TEMPERATURES_C) == 4
    assert len(LAPSE_RATES_K_M) == 4
