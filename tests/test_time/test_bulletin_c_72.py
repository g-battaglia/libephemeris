# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""UTC boundary facts announced by IERS Bulletin C 72 (6 July 2026).

Source: https://hpiers.obspm.fr/eoppc/bul/bulc/bulletinc.72
The bulletin excludes a leap second at the end of December 2026. It does not
supply UT1-UTC or an accuracy claim for predicted Earth orientation.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as ephem
from libephemeris.exceptions import Error


@pytest.mark.parametrize(
    "before,after,calendar",
    (
        ((2026, 12, 31), (2027, 1, 1), ephem.GREG_CAL),
        ((2026, 12, 18), (2026, 12, 19), ephem.JUL_CAL),
    ),
)
def test_no_leap_second_at_end_of_december_2026(
    before: tuple[int, int, int],
    after: tuple[int, int, int],
    calendar: int,
) -> None:
    """Refuse second 60 and advance TT by one SI second across midnight."""
    with pytest.raises(Error):
        ephem.utc_to_jd(*before, 23, 59, 60.0, calendar)

    tt_before, _ = ephem.utc_to_jd(*before, 23, 59, 59.0, calendar)
    tt_after, _ = ephem.utc_to_jd(*after, 0, 0, 0.0, calendar)
    elapsed_seconds = (tt_after - tt_before) * 86400.0
    jd_quantum_seconds = math.ulp(tt_after) * 86400.0
    assert abs(elapsed_seconds - 1.0) <= jd_quantum_seconds
