# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Regression tests for issue #64: refraction refuses what it cannot model.

Every refraction formula in the module is a fit valid over a physical
domain, and none of them checked their arguments. Fed something outside it,
the arithmetic still produced a float:

    refrac(180, 1013.25, 15, TRUE_TO_APP)        -> 9.83e+43 degrees
    refrac(nan, ...)                             -> nan
    calc_dip(1000, 0.0065, -1013.25, -273.15)    -> -14629 degrees
    calc_refrac_compat_app_to_true(10, 1e9, -273.15) -> 1.68e+08 degrees

Those are not errors a caller can notice — they are answers, and they flow
onward into whatever they feed. The guards refuse them instead.

Zero pressure keeps its documented meaning at every entry point, and a
negative lapse rate stays legal: it is a temperature inversion, a real
atmospheric condition.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as eph
from libephemeris import refraction
from libephemeris.exceptions import Error, InputValidationError

TRUE_TO_APP = eph.TRUE_TO_APP
APP_TO_TRUE = eph.APP_TO_TRUE

_NON_FINITE = [math.nan, math.inf, -math.inf]


# ---------------------------------------------------------------------------
# The values from the issue
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_the_reported_altitude_no_longer_returns_1e43():
    with pytest.raises(InputValidationError):
        eph.refrac(180, 1013.25, 15, TRUE_TO_APP)


@pytest.mark.unit
def test_the_reported_dip_no_longer_returns_minus_14629_degrees():
    with pytest.raises(InputValidationError):
        refraction.calc_dip(1000, 0.0065, -1013.25, -273.15)


@pytest.mark.unit
def test_the_reported_compat_reverse_no_longer_returns_1e8():
    with pytest.raises(InputValidationError):
        refraction.calc_refrac_compat_app_to_true(10, 1e9, -273.15)


@pytest.mark.unit
def test_the_reported_lapse_rate_is_refused_at_the_point_of_setting():
    with pytest.raises(InputValidationError):
        eph.set_lapse_rate(float("nan"))
    # The rejected value must not have been stored.
    assert math.isfinite(eph.get_lapse_rate())


# ---------------------------------------------------------------------------
# Altitude
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("alt", [90.001, 180.0, -90.001, -360.0, 1e9])
@pytest.mark.parametrize("flag", [TRUE_TO_APP, APP_TO_TRUE])
def test_an_altitude_outside_the_sky_is_refused(alt, flag):
    """An altitude outside [-90, 90] names no direction in the sky."""
    with pytest.raises(InputValidationError):
        eph.refrac(alt, 1013.25, 15.0, flag)


@pytest.mark.unit
@pytest.mark.parametrize("alt", [-90.0, -45.0, -5.0, 0.0, 10.0, 45.0, 90.0])
@pytest.mark.parametrize("flag", [TRUE_TO_APP, APP_TO_TRUE])
def test_every_altitude_in_the_sky_is_accepted(alt, flag):
    """Narrowing the domain must not cost a single legal altitude."""
    result = eph.refrac(alt, 1013.25, 15.0, flag)
    assert math.isfinite(result)


@pytest.mark.unit
@pytest.mark.parametrize("alt", _NON_FINITE)
def test_a_non_finite_altitude_is_refused(alt):
    with pytest.raises(InputValidationError):
        eph.refrac(alt, 1013.25, 15.0, TRUE_TO_APP)


# ---------------------------------------------------------------------------
# Pressure and temperature
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("pressure", [-0.001, -1013.25, -1e9])
def test_a_negative_absolute_pressure_is_refused(pressure):
    with pytest.raises(InputValidationError):
        eph.refrac(10.0, pressure, 15.0, TRUE_TO_APP)


@pytest.mark.unit
def test_zero_pressure_keeps_its_documented_meaning():
    """Zero disables refraction; it is not an invalid value."""
    assert eph.refrac(0.0, 0.0, 15.0, TRUE_TO_APP) == 0.0
    assert refraction.calc_refraction_true_to_app(45.0, pressure=0.0) == 0.0
    assert refraction.calc_refraction_app_to_true(45.0, pressure=0.0) == 0.0


@pytest.mark.unit
@pytest.mark.parametrize("temperature", [-273.15, -273.16, -300.0, -1e6])
def test_a_temperature_at_or_below_absolute_zero_is_refused(temperature):
    with pytest.raises(InputValidationError):
        eph.refrac(10.0, 1013.25, temperature, TRUE_TO_APP)


@pytest.mark.unit
def test_the_pole_of_the_conventional_scale_is_refused():
    """-273 C passes the physical guard but is where the scale diverges."""
    with pytest.raises(InputValidationError):
        eph.refrac(10.0, 1.0, -273.0, TRUE_TO_APP)


@pytest.mark.unit
@pytest.mark.parametrize("value", _NON_FINITE)
def test_a_non_finite_pressure_or_temperature_is_refused(value):
    with pytest.raises(InputValidationError):
        eph.refrac(10.0, value, 15.0, TRUE_TO_APP)
    with pytest.raises(InputValidationError):
        eph.refrac(10.0, 1013.25, value, TRUE_TO_APP)


# ---------------------------------------------------------------------------
# The ICAO surfaces and the dip carry the same guards
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize(
    "call",
    [
        lambda **kw: refraction.calc_refraction_true_to_app(10.0, **kw),
        lambda **kw: refraction.calc_refraction_app_to_true(10.0, **kw),
    ],
)
@pytest.mark.parametrize(
    "kwargs",
    [
        {"pressure": -1.0},
        {"pressure": math.nan},
        {"temperature_C": -273.15},
        {"temperature_C": math.inf},
        {"obs_alt": math.nan},
        {"lapse_rate": math.nan},
    ],
)
def test_the_icao_surfaces_refuse_the_same_arguments(call, kwargs):
    with pytest.raises(InputValidationError):
        call(**kwargs)


@pytest.mark.unit
@pytest.mark.parametrize("alt", [90.001, -90.001, math.nan, math.inf])
def test_the_icao_surfaces_refuse_an_impossible_altitude(alt):
    with pytest.raises(InputValidationError):
        refraction.calc_refraction_true_to_app(alt)
    with pytest.raises(InputValidationError):
        refraction.calc_refraction_app_to_true(alt)


@pytest.mark.unit
@pytest.mark.parametrize(
    "kwargs",
    [
        {"atpress": -1.0},
        {"atpress": math.nan},
        {"attemp": -273.15},
        {"attemp": -math.inf},
        {"lapse_rate": math.inf},
    ],
)
def test_the_dip_refuses_the_same_arguments(kwargs):
    with pytest.raises(InputValidationError):
        refraction.calc_dip(1000.0, **kwargs)


@pytest.mark.unit
@pytest.mark.parametrize("obs_alt", _NON_FINITE)
def test_the_dip_refuses_a_non_finite_observer_altitude(obs_alt):
    with pytest.raises(InputValidationError):
        refraction.calc_dip(obs_alt)


@pytest.mark.unit
def test_the_dip_still_answers_for_a_real_observer():
    assert refraction.calc_dip(1000.0) == pytest.approx(-0.873, abs=1e-3)
    assert refraction.calc_dip(0.0) == 0.0


# ---------------------------------------------------------------------------
# set_lapse_rate
# ---------------------------------------------------------------------------


class TestLapseRate:
    def teardown_method(self):
        eph.set_lapse_rate(None)

    @pytest.mark.unit
    @pytest.mark.parametrize("value", _NON_FINITE)
    def test_a_non_finite_lapse_rate_is_refused(self, value):
        with pytest.raises(InputValidationError):
            eph.set_lapse_rate(value)

    @pytest.mark.unit
    def test_a_rejected_lapse_rate_is_not_stored(self):
        eph.set_lapse_rate(0.005)
        with pytest.raises(InputValidationError):
            eph.set_lapse_rate(math.nan)
        assert eph.get_lapse_rate() == 0.005

    @pytest.mark.unit
    @pytest.mark.parametrize("value", [0.0, 0.0034, 0.0065, 0.010, -0.002])
    def test_physical_lapse_rates_are_accepted(self, value):
        """A negative lapse rate is a temperature inversion, not an error."""
        eph.set_lapse_rate(value)
        assert eph.get_lapse_rate() == value

    @pytest.mark.unit
    def test_none_still_resets_to_the_default(self):
        eph.set_lapse_rate(0.008)
        eph.set_lapse_rate(None)
        assert eph.get_lapse_rate() == 0.0065


# ---------------------------------------------------------------------------
# The errors stay inside the documented hierarchy
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_the_errors_are_catchable_as_the_library_base_error():
    """The typed hierarchy is the contract; no new exception family."""
    assert issubclass(InputValidationError, Error)
    with pytest.raises(Error):
        eph.refrac(180, 1013.25, 15, TRUE_TO_APP)


@pytest.mark.unit
def test_refrac_extended_is_guarded_through_the_same_helpers():
    with pytest.raises(InputValidationError):
        eph.refrac_extended(180.0, 0.0, 1013.25, 15.0)
    with pytest.raises(InputValidationError):
        eph.refrac_extended(10.0, 0.0, -1.0, 15.0)
    altitude, details = eph.refrac_extended(0.0, 0.0)
    assert math.isfinite(altitude)
    assert len(details) == 4
