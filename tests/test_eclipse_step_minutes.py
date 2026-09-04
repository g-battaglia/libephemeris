"""Non-positive sampling steps must not hang the eclipse path helpers."""

from __future__ import annotations

import pytest

from libephemeris.eclipse import (
    calc_eclipse_central_line,
    calc_eclipse_northern_limit,
    calc_eclipse_southern_limit,
)
from libephemeris.exceptions import Error

_SAMPLERS = (
    calc_eclipse_central_line,
    calc_eclipse_northern_limit,
    calc_eclipse_southern_limit,
)

# A short window where jd_start <= jd_end. A zero or negative step never
# satisfies the loop exit condition and would spin forever without a guard.
_JD_START = 2451545.0
_JD_END = 2451545.1


@pytest.mark.unit
@pytest.mark.parametrize("sampler", _SAMPLERS, ids=lambda fn: fn.__name__)
@pytest.mark.parametrize("step_minutes", (0.0, -1.0, -0.5))
def test_non_positive_step_minutes_raises_before_loop(sampler, step_minutes):
    """Reject a non-positive step with the library Error, before sampling."""
    with pytest.raises(Error):
        sampler(_JD_START, _JD_END, step_minutes=step_minutes)


@pytest.mark.unit
@pytest.mark.parametrize("sampler", _SAMPLERS, ids=lambda fn: fn.__name__)
def test_positive_step_minutes_accepted_on_empty_window(sampler):
    """A positive step is valid; an inverted window yields no samples."""
    times, lats, lons = sampler(_JD_END, _JD_START, step_minutes=1.0)
    assert times == ()
    assert lats == ()
    assert lons == ()
