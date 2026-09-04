"""Termination guard for the eclipse path samplers.

``calc_eclipse_central_line``, ``calc_eclipse_northern_limit`` and
``calc_eclipse_southern_limit`` walk their window with the accumulator loop
``while jd <= jd_end: ...; jd += step_days``. These tests pin down the
arguments that stop that walk from ever finishing, which must be refused
before the loop starts, and the arguments that must keep sampling exactly as
they did before the guard existed.

Bounded by construction
-----------------------
A test that fed a bad argument straight to the real loop would hang instead of
failing if the guard were ever removed, weakened or moved into the loop, and a
hanging suite is worse than the bug. Two properties keep that from happening:

* every rejection test installs :func:`sealed_loop_body`, which replaces the
  functions the loop body calls with a tripwire, so the first iteration raises
  at once rather than spinning;
* every acceptance test uses a window and a step that admit at most a few
  dozen iterations, so they terminate on their own merits.

Neither relies on a timeout: a thread or alarm watchdog is not a real bound
here, because the runaway loop is pure Python and would starve the very
interpreter that has to notice the timeout.
"""

from __future__ import annotations

import math
from typing import Any, Callable

import pytest

import libephemeris as ephem
from libephemeris import (
    InputValidationError,
    calc_eclipse_central_line,
    calc_eclipse_northern_limit,
    calc_eclipse_southern_limit,
)
from libephemeris import eclipse as eclipse_module
from libephemeris.eclipse import _eclipse_sampling_step_days

# Around the total solar eclipse of 2024-04-08, close to greatest eclipse.
JD_ECLIPSE = ephem.julday(2024, 4, 8, 18.3)

# Spacing between consecutive doubles at that Julian Day: ~4.66e-10 days,
# i.e. ~40 microseconds. Any step at or below half of it is absorbed by
# ``jd += step_days`` and the accumulator stands still.
JD_SPACING_DAYS = math.nextafter(JD_ECLIPSE, math.inf) - JD_ECLIPSE
MINUTES_PER_DAY = 24.0 * 60.0

INF = float("inf")
NAN = float("nan")

SAMPLERS = pytest.mark.parametrize(
    "sampler",
    [
        calc_eclipse_central_line,
        calc_eclipse_northern_limit,
        calc_eclipse_southern_limit,
    ],
    ids=["central_line", "northern_limit", "southern_limit"],
)

# Everything the three loop bodies call on their first statement onwards.
_LOOP_BODY_CALLEES = (
    "sol_eclipse_where",
    "calc_besselian_x",
    "calc_besselian_y",
    "calc_besselian_d",
    "calc_besselian_mu",
    "calc_besselian_l2",
)


class _SamplingLoopEntered(AssertionError):
    """Raised when a call that must be rejected reaches the loop body."""


@pytest.fixture
def sealed_loop_body(monkeypatch: pytest.MonkeyPatch) -> None:
    """Turn the first sampling iteration into an immediate failure.

    Every argument checked in this module must be refused before the walk
    starts. Sealing the loop body makes that testable and keeps the tests
    bounded: without the guard the first iteration raises
    :class:`_SamplingLoopEntered` instead of looping for ever.

    Args:
        monkeypatch: pytest fixture used to swap the module globals.
    """

    def _tripwire(*args: Any, **kwargs: Any) -> float:
        raise _SamplingLoopEntered(
            "the sampling loop ran with arguments that must be rejected "
            f"before it starts: args={args!r}"
        )

    for name in _LOOP_BODY_CALLEES:
        monkeypatch.setattr(eclipse_module, name, _tripwire)


def _expect_rejected(
    sampler: Callable[..., Any], *args: Any, **kwargs: Any
) -> InputValidationError:
    """Assert that ``sampler`` refuses the given arguments.

    Args:
        sampler: One of the three eclipse path samplers.
        *args: Positional arguments for the sampler.
        **kwargs: Keyword arguments for the sampler.

    Returns:
        The raised exception, for further assertions.
    """
    with pytest.raises(InputValidationError) as excinfo:
        sampler(*args, **kwargs)
    return excinfo.value


def _assert_path_shape(
    result: Any, jd_start: float, jd_end: float
) -> tuple[float, ...]:
    """Assert a sampler returned a well formed path inside its window.

    Args:
        result: Return value of a sampler.
        jd_start: First Julian Day of the requested window.
        jd_end: Last Julian Day of the requested window.

    Returns:
        The tuple of sample times.
    """
    assert isinstance(result, tuple)
    assert len(result) == 3
    times, lats, lons = result
    assert len(times) == len(lats) == len(lons)
    for index, jd in enumerate(times):
        assert isinstance(jd, float)
        assert jd_start <= jd <= jd_end
        if index:
            assert jd > times[index - 1]
    for lat in lats:
        assert -90.0 <= lat <= 90.0
    for lon in lons:
        assert -360.0 <= lon <= 360.0
    return times


class TestStepsThatNeverTerminate:
    """Steps for which ``jd`` never crosses ``jd_end``."""

    @SAMPLERS
    def test_zero_step_is_rejected(self, sampler, sealed_loop_body):
        """A zero step never moves the accumulator (the reported hang)."""
        error = _expect_rejected(
            sampler, JD_ECLIPSE, JD_ECLIPSE + 1.0 / 1440.0, step_minutes=0.0
        )
        assert "step_minutes" in str(error)

    @SAMPLERS
    def test_negative_zero_step_is_rejected(self, sampler, sealed_loop_body):
        """``-0.0`` passes a ``< 0`` test yet still never advances ``jd``."""
        assert not (-0.0 < 0.0)
        _expect_rejected(
            sampler, JD_ECLIPSE, JD_ECLIPSE + 1.0 / 1440.0, step_minutes=-0.0
        )

    @SAMPLERS
    @pytest.mark.parametrize("step_minutes", [-1.0, -60.0, -1e-30])
    def test_negative_step_is_rejected(self, sampler, step_minutes, sealed_loop_body):
        """A negative step walks away from ``jd_end`` for ever."""
        _expect_rejected(
            sampler,
            JD_ECLIPSE,
            JD_ECLIPSE + 1.0 / 1440.0,
            step_minutes=step_minutes,
        )

    @SAMPLERS
    @pytest.mark.parametrize(
        "step_minutes",
        [1e-9, 1e-7, 5e-324, 1e-320],
        ids=["1e-9", "1e-7", "smallest_denormal", "denormal"],
    )
    def test_sub_spacing_step_is_rejected(
        self, sampler, step_minutes, sealed_loop_body
    ):
        """A positive step below half the double spacing is absorbed.

        These steps are greater than zero, so a plain sign test lets them
        through, yet ``jd + step_days == jd`` at Julian Day 2.46e6 and the
        walk never ends. The denormal cases additionally underflow to exactly
        zero once divided by 1440.
        """
        assert step_minutes > 0.0
        assert JD_ECLIPSE + step_minutes / MINUTES_PER_DAY == JD_ECLIPSE
        _expect_rejected(
            sampler,
            JD_ECLIPSE,
            JD_ECLIPSE + 1.0 / 1440.0,
            step_minutes=step_minutes,
        )

    @SAMPLERS
    def test_half_spacing_step_is_rejected(self, sampler, sealed_loop_body):
        """Exactly half the spacing is a round-half-to-even coin flip.

        It advances from samples whose mantissa is odd and stalls on the
        others, so it cannot be relied on to terminate.
        """
        step_minutes = (JD_SPACING_DAYS / 2.0) * MINUTES_PER_DAY
        _expect_rejected(
            sampler,
            JD_ECLIPSE,
            JD_ECLIPSE + 1.0 / 1440.0,
            step_minutes=step_minutes,
        )

    @SAMPLERS
    @pytest.mark.parametrize(
        "jd_start, jd_end",
        [(JD_ECLIPSE, INF), (-INF, JD_ECLIPSE), (-INF, INF)],
        ids=["infinite_end", "infinite_start", "infinite_both"],
    )
    def test_infinite_window_is_rejected(
        self, sampler, jd_start, jd_end, sealed_loop_body
    ):
        """An infinite bound breaks the exit test itself.

        ``jd`` never exceeds an infinite ``jd_end``, and ``-inf`` plus any
        finite step is still ``-inf``.
        """
        _expect_rejected(sampler, jd_start, jd_end, step_minutes=1.0)


class TestStepsThatTerminateButLie:
    """Arguments that stop, but only by producing a meaningless path."""

    @SAMPLERS
    @pytest.mark.parametrize(
        "step_minutes", [NAN, INF, -INF], ids=["nan", "inf", "-inf"]
    )
    def test_non_finite_step_is_rejected(self, sampler, step_minutes, sealed_loop_body):
        """``nan`` truncates the walk and an infinite step skips the window.

        Neither is caught by a sign test: ``nan <= 0.0`` is false and ``inf``
        is positive. Both leave the loop after a single sample and return a
        path that cannot be told apart from a genuine partial eclipse.
        """
        error = _expect_rejected(
            sampler,
            JD_ECLIPSE,
            JD_ECLIPSE + 1.0 / 1440.0,
            step_minutes=step_minutes,
        )
        assert "finite" in str(error)

    @SAMPLERS
    @pytest.mark.parametrize(
        "jd_start, jd_end",
        [(NAN, JD_ECLIPSE), (JD_ECLIPSE, NAN), (NAN, NAN)],
        ids=["nan_start", "nan_end", "nan_both"],
    )
    def test_nan_window_is_rejected(self, sampler, jd_start, jd_end, sealed_loop_body):
        """A ``nan`` bound silently yields the "no central path" answer."""
        error = _expect_rejected(sampler, jd_start, jd_end, step_minutes=1.0)
        assert "finite" in str(error)


class TestGuardContract:
    """Where the rejection comes from and what it looks like."""

    @SAMPLERS
    def test_error_is_a_library_error(self, sampler, sealed_loop_body):
        """The refusal uses the library's own typed error hierarchy."""
        error = _expect_rejected(
            sampler, JD_ECLIPSE, JD_ECLIPSE + 1.0 / 1440.0, step_minutes=0.0
        )
        assert isinstance(error, ephem.Error)
        assert isinstance(error, ephem.InputValidationError)

    @SAMPLERS
    def test_message_names_the_function_and_parameter(self, sampler, sealed_loop_body):
        """The message points at the offending call and argument."""
        error = _expect_rejected(
            sampler, JD_ECLIPSE, JD_ECLIPSE + 1.0 / 1440.0, step_minutes=-5.0
        )
        assert sampler.__name__ in str(error)
        assert "step_minutes" in str(error)

    @SAMPLERS
    def test_rejection_precedes_any_calculation(self, sampler, sealed_loop_body):
        """The step is validated before the walk, not inside it.

        The window here is also far outside every ephemeris tier, so reaching
        the first sample would raise something else entirely.
        """
        error = _expect_rejected(sampler, 1.0e9, 1.0e9 + 1.0, step_minutes=0.0)
        assert "step_minutes" in str(error)

    def test_helper_returns_the_step_in_days(self):
        """Accepted steps come back converted to days, unchanged in value."""
        step_days = _eclipse_sampling_step_days(
            JD_ECLIPSE, JD_ECLIPSE + 1.0, 30.0, "test"
        )
        assert step_days == 30.0 / 1440.0

    def test_helper_boundary_is_half_the_double_spacing(self):
        """The accept/reject boundary sits exactly at half the spacing.

        Half the spacing is refused (ties round to even, so it advances from
        some samples only); the next double above it is accepted, and it does
        advance the accumulator.
        """
        half_spacing_minutes = (JD_SPACING_DAYS / 2.0) * MINUTES_PER_DAY
        with pytest.raises(InputValidationError):
            _eclipse_sampling_step_days(
                JD_ECLIPSE, JD_ECLIPSE, half_spacing_minutes, "test"
            )
        smallest_ok = math.nextafter(half_spacing_minutes, math.inf)
        step_days = _eclipse_sampling_step_days(
            JD_ECLIPSE, JD_ECLIPSE, smallest_ok, "test"
        )
        assert JD_ECLIPSE + step_days > JD_ECLIPSE

    def test_helper_uses_the_largest_magnitude_of_the_window(self):
        """The spacing is taken where the walk is hardest to advance.

        A step fine enough near JD 0 is still absorbed at JD 2.46e6, so the
        bound has to come from the larger end of the window, not the first.
        """
        step_minutes = 1e-7
        assert 1.0 + step_minutes / MINUTES_PER_DAY > 1.0
        assert JD_ECLIPSE + step_minutes / MINUTES_PER_DAY == JD_ECLIPSE
        with pytest.raises(InputValidationError):
            _eclipse_sampling_step_days(1.0, JD_ECLIPSE, step_minutes, "test")


class TestInputsThatMustKeepWorking:
    """Behaviour that predates the guard and must be preserved."""

    @SAMPLERS
    def test_documented_default_step_is_accepted(self, sampler):
        """The documented default (1 minute) still samples the window."""
        jd_end = JD_ECLIPSE + 3.0 / 1440.0
        _assert_path_shape(sampler(JD_ECLIPSE, jd_end), JD_ECLIPSE, jd_end)

    @SAMPLERS
    @pytest.mark.parametrize("step_minutes", [1.0, 5.0, 10.0, 30.0])
    def test_documented_steps_are_accepted(self, sampler, step_minutes):
        """The step sizes used in the docstring examples still work."""
        jd_start = JD_ECLIPSE - 30.0 / 1440.0
        jd_end = JD_ECLIPSE + 30.0 / 1440.0
        result = sampler(jd_start, jd_end, step_minutes=step_minutes)
        _assert_path_shape(result, jd_start, jd_end)

    def test_central_line_still_finds_the_2024_path(self):
        """A coarse sampling of the 2024-04-08 eclipse still returns points."""
        jd_start = JD_ECLIPSE - 15.0 / 1440.0
        jd_end = JD_ECLIPSE + 15.0 / 1440.0
        times, lats, lons = calc_eclipse_central_line(
            jd_start, jd_end, step_minutes=15.0
        )
        assert len(times) >= 1
        # Greatest eclipse crossed Mexico and the United States.
        assert all(0.0 < lat < 60.0 for lat in lats)
        assert all(-140.0 < lon < -20.0 for lon in lons)

    @SAMPLERS
    def test_smallest_accepted_step_terminates(self, sampler):
        """One spacing unit is accepted and the walk still ends.

        The window is a single instant, so this samples at most once however
        small the accepted step is.
        """
        step_minutes = JD_SPACING_DAYS * MINUTES_PER_DAY
        result = sampler(JD_ECLIPSE, JD_ECLIPSE, step_minutes=step_minutes)
        _assert_path_shape(result, JD_ECLIPSE, JD_ECLIPSE)

    @SAMPLERS
    def test_reversed_window_still_returns_empty(self, sampler):
        """A window that ends before it starts keeps yielding empty tuples."""
        assert sampler(JD_ECLIPSE + 1.0, JD_ECLIPSE, step_minutes=1.0) == (
            (),
            (),
            (),
        )

    @SAMPLERS
    def test_step_larger_than_the_window_is_accepted(self, sampler):
        """A step past ``jd_end`` is fine: it samples the start only."""
        jd_end = JD_ECLIPSE + 1.0 / 1440.0
        result = sampler(JD_ECLIPSE, jd_end, step_minutes=1e6)
        assert len(_assert_path_shape(result, JD_ECLIPSE, jd_end)) <= 1
