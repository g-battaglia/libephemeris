"""Output-mode coverage for restored historical fictitious bodies."""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    FLG_SPEED,
    FLG_HELCTR,
    FLG_EQUATORIAL,
    FLG_J2000,
    FLG_NOABERR,
    FLG_SIDEREAL,
)

FICTITIOUS_BODIES = [
    (40, "Cupido"),
    (41, "Hades"),
    (42, "Zeus"),
    (43, "Kronos"),
    (44, "Apollon"),
    (45, "Admetos"),
    (46, "Vulkanus"),
    (47, "Poseidon"),
    (48, "Transpluto"),
    (49, "Nibiru"),
    (50, "Harrington"),
    (51, "Leverrier"),
    (52, "Adams"),
    (53, "Lowell"),
    (54, "Pickering"),
    (55, "Vulcan"),
    (56, "White Moon"),
    (57, "Proserpina"),
    (58, "Waldemath"),
]


@pytest.mark.unit
@pytest.mark.parametrize("body_id,name", FICTITIOUS_BODIES)
@pytest.mark.parametrize(
    "flags",
    [
        0,
        FLG_SPEED,
        FLG_EQUATORIAL,
        FLG_J2000,
        FLG_NOABERR,
        FLG_HELCTR,
        FLG_SPEED | FLG_EQUATORIAL,
        FLG_HELCTR | FLG_SPEED,
        FLG_SIDEREAL,
    ],
)
def test_restored_fictitious_body_computes_for_all_output_modes(
    body_id: int, name: str, flags: int
) -> None:
    """Every historical ID returns a finite state across public output flags."""
    swe.set_calc_mode("skyfield")
    position, _ = swe.calc_ut(2451545.0, body_id, flags)
    assert len(position) == 6, name
    assert all(math.isfinite(value) for value in position), name
