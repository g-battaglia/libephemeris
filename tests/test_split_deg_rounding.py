# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Regression tests for issue #56: rounding must not leak its own offset.

``split_deg()`` rounds by adding half of the requested unit and then
splitting the result. The half-unit has to be discarded once it has done its
job; leaving it in place put it straight into the output, so rounding 0
degrees to the nearest minute reported 0d00'30" — half a minute that came
from the algorithm, not from the input.

The error was in the returned values, not only in a rendering of them: both
the integer seconds and the fractional seconds carried it.
"""

from __future__ import annotations

import pytest

from libephemeris.utils import (
    SPLIT_DEG_KEEP_DEG,
    SPLIT_DEG_KEEP_SIGN,
    SPLIT_DEG_NAKSHATRA,
    SPLIT_DEG_ROUND_DEG,
    SPLIT_DEG_ROUND_MIN,
    SPLIT_DEG_ROUND_SEC,
    SPLIT_DEG_ZODIACAL,
    split_deg,
)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("degree", "expected"),
    [
        (0.0, (0, 0, 0, 0.0, 1)),
        (10.6, (10, 36, 0, 0.0, 1)),  # exactly 10d36'00"
        (45.5, (45, 30, 0, 0.0, 1)),
        (45.999, (46, 0, 0, 0.0, 1)),  # 45d59'56.4" rounds up a degree
    ],
)
def test_round_min_leaves_no_seconds(degree, expected):
    """Minute rounding zeroes the seconds and the sub-second fraction."""
    assert split_deg(degree, SPLIT_DEG_ROUND_MIN) == expected


@pytest.mark.unit
@pytest.mark.parametrize(
    ("degree", "expected"),
    [
        (0.0, (0, 0, 0, 0.0, 1)),
        (10.6, (11, 0, 0, 0.0, 1)),
        (10.4, (10, 0, 0, 0.0, 1)),
        (45.5, (46, 0, 0, 0.0, 1)),
    ],
)
def test_round_deg_leaves_no_minutes_or_seconds(degree, expected):
    """Degree rounding zeroes the minutes, seconds and fraction."""
    assert split_deg(degree, SPLIT_DEG_ROUND_DEG) == expected


@pytest.mark.unit
@pytest.mark.parametrize("degree", [0.0, 10.6, 45.5, 123.456789, 359.99])
def test_round_sec_leaves_no_sub_second_fraction(degree):
    """Second rounding leaves no fraction: there is nothing below the unit."""
    assert split_deg(degree, SPLIT_DEG_ROUND_SEC)[3] == 0.0


@pytest.mark.unit
def test_unrounded_split_is_unchanged():
    """Without a rounding flag the decomposition keeps full precision."""
    assert split_deg(0.0, 0) == (0, 0, 0, 0.0, 1)
    assert split_deg(45.5, 0) == (45, 30, 0, 0.0, 1)
    ideg, imin, isec, secfr, sign = split_deg(10.60001, 0)
    assert (ideg, imin, isec, sign) == (10, 36, 0, 1)
    assert 0.0 < secfr < 1.0


@pytest.mark.unit
def test_negative_values_round_the_same_way():
    """The sign is carried separately, so rounding behaves identically."""
    assert split_deg(-10.6, SPLIT_DEG_ROUND_MIN) == (10, 36, 0, 0.0, -1)
    assert split_deg(-10.6, SPLIT_DEG_ROUND_DEG) == (11, 0, 0, 0.0, -1)


@pytest.mark.unit
def test_zodiacal_rounding_leaves_no_offset():
    """The zodiacal split rounds within its sign without a leftover."""
    # 40.6 deg is 10d36' of Taurus (sign index 1).
    assert split_deg(40.6, SPLIT_DEG_ZODIACAL | SPLIT_DEG_ROUND_MIN) == (
        10,
        36,
        0,
        0.0,
        1,
    )


@pytest.mark.unit
def test_nakshatra_rounding_leaves_no_offset():
    """The nakshatra split rounds within its segment without a leftover."""
    assert split_deg(0.09, SPLIT_DEG_NAKSHATRA | SPLIT_DEG_ROUND_MIN) == (
        0,
        5,
        0,
        0.0,
        0,
    )


@pytest.mark.unit
def test_keep_flags_still_suppress_the_offset():
    """KEEP_DEG and KEEP_SIGN keep their meaning after the truncation."""
    # 45.6 rounds to 46 deg, so KEEP_DEG suppresses the offset entirely.
    assert split_deg(45.6, SPLIT_DEG_ROUND_DEG | SPLIT_DEG_KEEP_DEG) == (
        45,
        0,
        0,
        0.0,
        1,
    )
    # 29.9 would cross into the next sign, so KEEP_SIGN suppresses it.
    result = split_deg(
        29.9, SPLIT_DEG_ROUND_DEG | SPLIT_DEG_KEEP_SIGN | SPLIT_DEG_ZODIACAL
    )
    assert result == (29, 0, 0, 0.0, 0)


@pytest.mark.unit
@pytest.mark.parametrize(
    "flag", [SPLIT_DEG_ROUND_DEG, SPLIT_DEG_ROUND_MIN, SPLIT_DEG_ROUND_SEC]
)
@pytest.mark.parametrize("degree", [0.0, 0.09, 10.6, 45.5, 123.456789, 359.99])
def test_no_rounding_flag_ever_reports_a_sub_second_fraction(flag, degree):
    """A rounded value has, by definition, nothing left below its unit."""
    assert split_deg(degree, flag)[3] == 0.0
    assert split_deg(degree, flag | SPLIT_DEG_ZODIACAL)[3] == 0.0
