"""Independent invariants for the closed-form ``refrac()`` API.

The tests intentionally contain no numeric output captured from another
ephemeris. They exercise properties of the published Bennett/Saemundsson
formula family and the public function's input guards.
"""

from __future__ import annotations

import math

import pytest

from libephemeris import APP_TO_TRUE, TRUE_TO_APP, refrac


@pytest.mark.unit
@pytest.mark.parametrize("altitude", [-4.0, 0.0, 10.0, 45.0, 89.0])
def test_positive_pressure_lifts_true_altitude(altitude: float) -> None:
    apparent = refrac(altitude, 1013.25, 15.0, TRUE_TO_APP)
    if apparent != altitude:  # below-horizon guard may retain the input
        assert apparent > altitude


@pytest.mark.unit
@pytest.mark.parametrize("altitude", [1.0, 5.0, 10.0, 45.0, 85.0])
def test_positive_pressure_reduces_apparent_altitude(altitude: float) -> None:
    true = refrac(altitude, 1013.25, 15.0, APP_TO_TRUE)
    assert true < altitude


@pytest.mark.unit
def test_refraction_grows_with_pressure() -> None:
    altitude = 10.0
    low = refrac(altitude, 500.0, 15.0, TRUE_TO_APP) - altitude
    high = refrac(altitude, 1000.0, 15.0, TRUE_TO_APP) - altitude
    assert high > low > 0.0


@pytest.mark.unit
def test_refraction_decreases_with_temperature() -> None:
    altitude = 10.0
    cold = refrac(altitude, 1013.25, -20.0, TRUE_TO_APP) - altitude
    warm = refrac(altitude, 1013.25, 35.0, TRUE_TO_APP) - altitude
    assert cold > warm > 0.0


@pytest.mark.unit
def test_correction_decreases_toward_zenith() -> None:
    corrections = [
        refrac(alt, 1013.25, 15.0, TRUE_TO_APP) - alt
        for alt in (5.0, 10.0, 30.0, 60.0, 85.0)
    ]
    assert all(a > b > 0.0 for a, b in zip(corrections, corrections[1:]))


@pytest.mark.unit
@pytest.mark.parametrize("flag", [-1, 2, 999])
def test_every_nonzero_flag_selects_reverse_direction(flag: int) -> None:
    assert refrac(10.0, 1013.25, 15.0, flag) == refrac(10.0, 1013.25, 15.0, APP_TO_TRUE)


@pytest.mark.unit
def test_nonfinite_inputs_are_guarded() -> None:
    assert math.isnan(refrac(math.nan, 1013.25, 15.0, TRUE_TO_APP))
    assert refrac(math.inf, 1013.25, 15.0, TRUE_TO_APP) == math.inf
    assert refrac(-math.inf, 1013.25, 15.0, APP_TO_TRUE) == -math.inf
    assert refrac(10.0, math.nan, 15.0, TRUE_TO_APP) == 10.0


@pytest.mark.unit
def test_singular_temperature_does_not_raise() -> None:
    assert math.isinf(refrac(10.0, 1.0, -273.0, TRUE_TO_APP))


@pytest.mark.unit
def test_guarded_results_are_native_floats() -> None:
    assert isinstance(refrac(-5, 1013.25, 15.0, TRUE_TO_APP), float)
    assert isinstance(refrac(0, 1013.25, 15.0, APP_TO_TRUE), float)


@pytest.mark.unit
def test_domain_floor_is_not_the_horizon_clamp() -> None:
    """Below its floor the forward direction has no fit, whatever the air.

    An extreme pressure lifts an altitude just above the floor by thousands
    of degrees, so returning the input at the floor cannot be explained by
    the horizon clamp alone.
    """
    assert refrac(-5.0, 1e9, 15.0, TRUE_TO_APP) == -5.0
    assert refrac(-5.001, 1e9, 15.0, TRUE_TO_APP) == -5.001
    assert refrac(-4.99, 1e9, 15.0, TRUE_TO_APP) > 1000.0


@pytest.mark.unit
def test_reverse_stops_correcting_near_the_zenith() -> None:
    """Bennett's argument leaves the first quadrant just below 89.93 degrees."""
    assert refrac(89.92, 1013.25, 15.0, APP_TO_TRUE) < 89.92
    for altitude in (89.93, 89.95, 90.0):
        assert refrac(altitude, 1013.25, 15.0, APP_TO_TRUE) == altitude


@pytest.mark.unit
def test_forward_correction_vanishes_at_the_zenith() -> None:
    """The tangent of a zero zenith distance vanishes, so the zenith is exact."""
    assert refrac(90.0, 1013.25, 15.0, TRUE_TO_APP) == 90.0


@pytest.mark.unit
def test_forward_correction_steps_down_just_above_fifteen_degrees() -> None:
    """Two published fits meet at 15 degrees and do not join smoothly."""
    corrections = [
        refrac(altitude, 1013.25, 15.0, TRUE_TO_APP) - altitude
        for altitude in (14.99, 15.0, 15.01)
    ]
    smooth = corrections[0] - corrections[1]
    step = corrections[1] - corrections[2]
    assert smooth > 0.0
    assert step > 10.0 * smooth


@pytest.mark.unit
def test_negative_zero_keeps_its_sign_only_where_nothing_is_applied() -> None:
    """Negative zero is an altitude of zero, and the clamp hands it back as it came."""
    assert refrac(-0.0, 1013.25, 15.0, TRUE_TO_APP) == refrac(
        0.0, 1013.25, 15.0, TRUE_TO_APP
    )
    reverse = refrac(-0.0, 1013.25, 15.0, APP_TO_TRUE)
    assert math.copysign(1.0, reverse) == -1.0
    assert math.copysign(1.0, refrac(0.0, 1013.25, 15.0, APP_TO_TRUE)) == 1.0


@pytest.mark.unit
def test_pole_of_the_reverse_fit_raises_nothing() -> None:
    """At -4.4 degrees Bennett's argument is singular; the input comes back."""
    assert refrac(-4.4, 1013.25, 15.0, APP_TO_TRUE) == -4.4


@pytest.mark.unit
def test_horizon_lift_scales_with_pressure_and_inverse_temperature() -> None:
    """The atmosphere factor is linear in pressure and inverse in 273 + T."""
    lift_900 = refrac(0.0, 900.0, 15.0, TRUE_TO_APP)
    lift_1050 = refrac(0.0, 1050.0, 15.0, TRUE_TO_APP)
    assert lift_1050 / lift_900 == pytest.approx(1050.0 / 900.0, rel=1e-12)
    cold = refrac(0.0, 1013.25, -40.0, TRUE_TO_APP)
    warm = refrac(0.0, 1013.25, 40.0, TRUE_TO_APP)
    assert cold / warm == pytest.approx((273.0 + 40.0) / (273.0 - 40.0), rel=1e-12)


@pytest.mark.unit
def test_negative_pressure_reverses_both_directions() -> None:
    """A negative factor is a compatibility extrapolation, not an error."""
    assert refrac(10.0, -1013.25, 15.0, TRUE_TO_APP) < 10.0
    assert refrac(10.0, -1013.25, 15.0, APP_TO_TRUE) > 10.0


@pytest.mark.unit
def test_infinite_temperature_removes_the_refraction() -> None:
    """An infinite temperature drives the factor to zero in both directions."""
    for flag in (TRUE_TO_APP, APP_TO_TRUE):
        assert refrac(10.0, 1013.25, math.inf, flag) == 10.0
        assert refrac(10.0, 1013.25, -math.inf, flag) == 10.0


@pytest.mark.unit
def test_infinite_pressure_gives_an_infinite_lift() -> None:
    """An infinite factor is answered, not raised, where nothing clamps it."""
    assert refrac(10.0, math.inf, 15.0, TRUE_TO_APP) == math.inf
    assert refrac(10.0, math.inf, 15.0, APP_TO_TRUE) == 10.0


@pytest.mark.unit
def test_integer_arguments_give_a_native_float() -> None:
    """Integers are accepted for all three numbers and never leak into the result."""
    assert isinstance(refrac(10, 1013, 15, TRUE_TO_APP), float)
    assert isinstance(refrac(10, 1013, 15, APP_TO_TRUE), float)
