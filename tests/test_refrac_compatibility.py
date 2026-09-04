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
    """A non-finite argument is refused, not carried through the result."""
    from libephemeris.exceptions import InputValidationError

    with pytest.raises(InputValidationError):
        refrac(math.nan, 1013.25, 15.0, TRUE_TO_APP)
    with pytest.raises(InputValidationError):
        refrac(math.inf, 1013.25, 15.0, TRUE_TO_APP)
    with pytest.raises(InputValidationError):
        refrac(-math.inf, 1013.25, 15.0, APP_TO_TRUE)
    with pytest.raises(InputValidationError):
        refrac(10.0, math.nan, 15.0, TRUE_TO_APP)


@pytest.mark.unit
def test_singular_temperature_is_refused() -> None:
    """The pole of the conventional scale has no finite refraction.

    -273 C is a fraction of a kelvin above absolute zero, so it passes the
    physical guard, but it is exactly where the 273 + T denominator
    vanishes. Returning an infinity there is the extreme case of the answer
    this module must not produce.
    """
    from libephemeris.exceptions import InputValidationError

    with pytest.raises(InputValidationError):
        refrac(10.0, 1.0, -273.0, TRUE_TO_APP)


@pytest.mark.unit
def test_guarded_results_are_native_floats() -> None:
    assert isinstance(refrac(-5, 1013.25, 15.0, TRUE_TO_APP), float)
    assert isinstance(refrac(0, 1013.25, 15.0, APP_TO_TRUE), float)
