# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Regression tests for issue #61: a path sampler must always terminate.

``calc_eclipse_central_line()``, ``calc_eclipse_northern_limit()`` and
``calc_eclipse_southern_limit()`` walk an instant from ``jd_start`` and stop
once it passes ``jd_end``. That exit depends on the instant actually
advancing, so a step that does not advance it spins forever:

* zero or negative moves the instant backwards or not at all;
* a step small enough that ``jd + step`` rounds back to ``jd`` is inert, and
  a sign test accepts it — Julian dates are ~2.4e6, so their representable
  spacing is around 5e-10 d, roughly 7e-7 minutes;
* two infinite bounds keep ``jd <= jd_end`` true with no possible progress.

Every case below is bounded: the guard is asserted to run *before* the loop
by making the path function's per-step dependency fail if it is ever reached. A
test that called the real sampler with a bad step would hang rather than
fail if the guard were ever moved into or after the loop — which is exactly
the defect being fixed.
"""

from __future__ import annotations

import math
from unittest.mock import patch

import pytest

from libephemeris import eclipse
from libephemeris.eclipse import (
    calc_eclipse_central_line,
    calc_eclipse_northern_limit,
    calc_eclipse_southern_limit,
)
from libephemeris.exceptions import Error, InputValidationError

# A short, well-formed window in the middle of the modern era. Nothing about
# these two instants matters except that they are ordered and finite: the
# step, not the window, is what every case below is about.
_WINDOW_OPEN = 2455000.5
_WINDOW_CLOSE = _WINDOW_OPEN + 0.25

_PATH_FUNCTIONS = [
    ("central_line", calc_eclipse_central_line),
    ("northern_limit", calc_eclipse_northern_limit),
    ("southern_limit", calc_eclipse_southern_limit),
]

# The per-step work each sampler performs. Patching these to raise turns
# "the guard runs first" into a bounded assertion instead of a hang.
_PER_STEP_DEPENDENCIES = ("sol_eclipse_where", "calc_besselian_x")


def _fail_if_sampled(*args, **kwargs):
    raise AssertionError("the sampling loop was entered despite an invalid step")


def _no_sampling_allowed():
    """Patch every per-step dependency so entering the loop fails fast."""
    from contextlib import ExitStack

    stack = ExitStack()
    for name in _PER_STEP_DEPENDENCIES:
        stack.enter_context(patch.object(eclipse, name, _fail_if_sampled))
    return stack


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "path_function"),
    _PATH_FUNCTIONS,
    ids=lambda v: v if isinstance(v, str) else "",
)
@pytest.mark.parametrize("step_minutes", [0.0, -0.0, -1.0, -0.5, -1e6])
def test_a_non_positive_step_is_refused_before_sampling(
    name, path_function, step_minutes
):
    """The reported case: these calls used to hang the interpreter."""
    with _no_sampling_allowed():
        with pytest.raises(InputValidationError):
            path_function(_WINDOW_OPEN, _WINDOW_CLOSE, step_minutes=step_minutes)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "path_function"),
    _PATH_FUNCTIONS,
    ids=lambda v: v if isinstance(v, str) else "",
)
@pytest.mark.parametrize("step_minutes", [1e-9, 1e-12, 5e-322])
def test_a_step_too_small_to_advance_the_date_is_refused(
    name, path_function, step_minutes
):
    """Positive, but inert: `jd + step` is the same float as `jd`."""
    assert _WINDOW_OPEN + step_minutes / (24.0 * 60.0) == _WINDOW_OPEN
    with _no_sampling_allowed():
        with pytest.raises(InputValidationError):
            path_function(_WINDOW_OPEN, _WINDOW_CLOSE, step_minutes=step_minutes)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "path_function"),
    _PATH_FUNCTIONS,
    ids=lambda v: v if isinstance(v, str) else "",
)
@pytest.mark.parametrize("step_minutes", [math.nan, math.inf, -math.inf])
def test_a_non_finite_step_is_refused(name, path_function, step_minutes):
    with _no_sampling_allowed():
        with pytest.raises(InputValidationError):
            path_function(_WINDOW_OPEN, _WINDOW_CLOSE, step_minutes=step_minutes)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "path_function"),
    _PATH_FUNCTIONS,
    ids=lambda v: v if isinstance(v, str) else "",
)
@pytest.mark.parametrize(
    ("jd_start", "jd_end"),
    [
        (math.nan, _WINDOW_CLOSE),
        (_WINDOW_OPEN, math.nan),
        (math.inf, math.inf),
        (-math.inf, math.inf),
    ],
)
def test_a_non_finite_window_is_refused(name, path_function, jd_start, jd_end):
    """Two infinite bounds keep the loop condition true with no progress."""
    with _no_sampling_allowed():
        with pytest.raises(InputValidationError):
            path_function(jd_start, jd_end, step_minutes=1.0)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "path_function"),
    _PATH_FUNCTIONS,
    ids=lambda v: v if isinstance(v, str) else "",
)
@pytest.mark.parametrize("step_minutes", [1e-3, 1.0, 30.0, 1440.0])
def test_a_step_that_advances_the_date_is_accepted(name, path_function, step_minutes):
    """Narrowing the input must not cost a step that can actually sample.

    The window is sized to three steps so the assertion is about acceptance,
    not about walking a real eclipse: a 1e-3 minute step over a tenth of a
    day would be 144 000 evaluations.
    """
    step_days = step_minutes / (24.0 * 60.0)
    jd_end = _WINDOW_OPEN + 3.0 * step_days

    times, latitudes, longitudes = path_function(
        _WINDOW_OPEN, jd_end, step_minutes=step_minutes
    )
    assert isinstance(times, tuple)
    assert len(times) == len(latitudes) == len(longitudes)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "path_function"),
    _PATH_FUNCTIONS,
    ids=lambda v: v if isinstance(v, str) else "",
)
def test_an_inverted_window_returns_nothing_rather_than_raising(name, path_function):
    """jd_end before jd_start is an empty window, not an invalid one."""
    assert path_function(_WINDOW_CLOSE, _WINDOW_OPEN, step_minutes=1.0) == ((), (), ())


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "path_function"),
    _PATH_FUNCTIONS,
    ids=lambda v: v if isinstance(v, str) else "",
)
def test_the_default_step_is_still_valid(name, path_function):
    """The documented default must not be rejected by its own guard."""
    times, _latitudes, _longitudes = path_function(_WINDOW_OPEN, _WINDOW_OPEN + 0.01)
    assert isinstance(times, tuple)


@pytest.mark.unit
def test_the_message_names_the_function_and_the_offending_value():
    with pytest.raises(InputValidationError) as excinfo:
        calc_eclipse_central_line(_WINDOW_OPEN, _WINDOW_CLOSE, step_minutes=0.0)
    message = str(excinfo.value)
    assert "calc_eclipse_central_line" in message
    assert "step_minutes" in message


@pytest.mark.unit
def test_the_error_stays_inside_the_documented_hierarchy():
    assert issubclass(InputValidationError, Error)
    with pytest.raises(Error):
        calc_eclipse_central_line(_WINDOW_OPEN, _WINDOW_CLOSE, step_minutes=-1.0)


@pytest.mark.unit
def test_the_advancement_rule_follows_the_date_magnitude():
    """The test is against the caller's own jd_start, not a fixed epsilon.

    Representable spacing grows with the magnitude of the Julian date, so a
    step that is fine near JD 0 can be inert at JD 2.4e6.
    """
    step_minutes = 1e-7
    step_days = step_minutes / (24.0 * 60.0)
    assert 1.0 + step_days > 1.0  # advances a small date
    assert _WINDOW_OPEN + step_days == _WINDOW_OPEN  # inert at a Julian date

    with _no_sampling_allowed():
        with pytest.raises(InputValidationError):
            calc_eclipse_central_line(
                _WINDOW_OPEN, _WINDOW_CLOSE, step_minutes=step_minutes
            )
