# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Regression tests for issue #62: sealed mode refuses on the date alone.

In sealed ``leb`` mode the interpolated apsides answered for every date
tested, down to year -13000 and up to +17000, while every other body raised
``EphemerisRangeError`` outside the loaded artefact's coverage.

The cause is not that these points skip a check but that they never reach
one. ``MEAN_NODE``, ``MEAN_APOG``, ``INTP_APOG`` and ``INTP_PERG`` are all
served from packaged runtime models rather than from a stored channel, so no
reader lookup fails for them and the sealed range policy — which is applied
where a reader exception is handled — was never invoked.

What made it look like only the apsides were affected is that the policy did
fire whenever the requested frame happened to need nutation, and the four
points need it under complementary conditions:

    MEAN_NODE / MEAN_APOG   refused of-date, answered under FLG_NONUT
    INTP_APOG / INTP_PERG   answered of-date, refused under FLG_NONUT

So all four leaked; the flag decided which. The refusal now depends on the
date alone.
"""

from __future__ import annotations

import pytest

import libephemeris as eph
from libephemeris.exceptions import EphemerisRangeError

# Well outside every canonical tier: year -3000, JD 0, and year 20000.
_OUT_OF_RANGE_JD = [625673.5, 0.0, 9000000.5]
_IN_RANGE_JD = 2451545.0

_MODEL_SERVED_POINTS = [
    ("MEAN_NODE", eph.MEAN_NODE),
    ("MEAN_APOG", eph.MEAN_APOG),
    ("INTP_APOG", eph.INTP_APOG),
    ("INTP_PERG", eph.INTP_PERG),
]

# The frames that decided the outcome before the guard became deterministic.
_FRAME_FLAGS = [
    ("of-date", 0),
    ("nonut", eph.FLG_NONUT),
    ("j2000", eph.FLG_J2000),
    ("speed", eph.FLG_SPEED),
    ("nonut-speed", eph.FLG_NONUT | eph.FLG_SPEED),
]


@pytest.fixture
def sealed_leb():
    """Run the body of the test in sealed leb mode, then restore the mode."""
    previous = eph.get_calc_mode()
    eph.set_calc_mode("leb")
    try:
        yield
    finally:
        eph.set_calc_mode(previous)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "body"),
    _MODEL_SERVED_POINTS,
    ids=lambda v: v if isinstance(v, str) else "",
)
@pytest.mark.parametrize(
    ("frame", "flags"), _FRAME_FLAGS, ids=lambda v: v if isinstance(v, str) else ""
)
@pytest.mark.parametrize("jd", _OUT_OF_RANGE_JD)
def test_model_served_points_refuse_outside_coverage(
    sealed_leb, name, body, frame, flags, jd
):
    """Sealed mode promises a refusal; the frame must not decide it."""
    with pytest.raises(EphemerisRangeError):
        eph.calc_ut(jd, body, flags)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "body"),
    _MODEL_SERVED_POINTS,
    ids=lambda v: v if isinstance(v, str) else "",
)
@pytest.mark.parametrize(
    ("frame", "flags"), _FRAME_FLAGS, ids=lambda v: v if isinstance(v, str) else ""
)
def test_model_served_points_still_answer_inside_coverage(
    sealed_leb, name, body, frame, flags
):
    """Narrowing the range must not disturb an in-coverage request."""
    position, _retflag = eph.calc_ut(_IN_RANGE_JD, body, flags)
    assert 0.0 <= position[0] < 360.0


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "body"),
    _MODEL_SERVED_POINTS,
    ids=lambda v: v if isinstance(v, str) else "",
)
def test_the_tt_entry_point_agrees_with_the_ut_one(sealed_leb, name, body):
    """calc() carried the same gap as calc_ut() and gets the same guard."""
    with pytest.raises(EphemerisRangeError):
        eph.calc(_OUT_OF_RANGE_JD[0], body, 0)
    position, _retflag = eph.calc(_IN_RANGE_JD, body, 0)
    assert 0.0 <= position[0] < 360.0


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "body"),
    _MODEL_SERVED_POINTS,
    ids=lambda v: v if isinstance(v, str) else "",
)
def test_the_error_describes_the_coverage_it_enforced(sealed_leb, name, body):
    """The typed error carries the bounds, as it does for every other body."""
    with pytest.raises(EphemerisRangeError) as excinfo:
        eph.calc_ut(_OUT_OF_RANGE_JD[0], body, 0)
    error = excinfo.value
    assert error.body_id == body
    assert error.start_jd is not None and error.end_jd is not None
    assert error.start_jd < error.end_jd


@pytest.mark.unit
def test_the_apsides_agree_with_the_other_bodies(sealed_leb):
    """The reported asymmetry, stated directly: same date, same outcome."""
    for body in (eph.SUN, eph.MOON, eph.MEAN_NODE, eph.INTP_APOG, eph.INTP_PERG):
        with pytest.raises(EphemerisRangeError):
            eph.calc_ut(_OUT_OF_RANGE_JD[0], body)


@pytest.mark.unit
@pytest.mark.parametrize("mode", ["auto", "skyfield"])
@pytest.mark.parametrize(
    ("name", "body"),
    [("INTP_APOG", eph.INTP_APOG), ("INTP_PERG", eph.INTP_PERG)],
    ids=lambda v: v if isinstance(v, str) else "",
)
def test_unsealed_modes_keep_the_documented_taper(mode, name, body):
    """Only sealed mode promises a refusal.

    Outside the fitted interval the model's residual table tapers to zero and
    the remaining trig-only longitude is documented with its own error bound.
    That degradation is deliberate and stays available in auto and skyfield.
    """
    previous = eph.get_calc_mode()
    eph.set_calc_mode(mode)
    try:
        position, _retflag = eph.calc_ut(_OUT_OF_RANGE_JD[0], body)
    finally:
        eph.set_calc_mode(previous)
    assert 0.0 <= position[0] < 360.0
