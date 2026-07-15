"""Output-mode coverage for independently sourced fictitious bodies."""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    FLG_EQUATORIAL,
    FLG_HELCTR,
    FLG_J2000,
    FLG_NOABERR,
    FLG_SIDEREAL,
    FLG_SPEED,
)
from libephemeris.exceptions import UnknownBodyError

UNVERIFIED_FICTITIOUS_BODIES = [
    (48, "Transpluto"),
    (49, "Nibiru"),
    (54, "Pickering"),
    (55, "Vulcan"),
    (57, "Proserpina"),
    (58, "Waldemath"),
]

REVIEWED_FICTITIOUS_BODIES = [*range(40, 48), 50, 51, 52, 53, 56]


@pytest.mark.unit
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
def test_harrington_computes_for_all_output_modes(flags: int) -> None:
    """The primary-derived Harrington orbit supports public output flags."""
    swe.set_calc_mode("skyfield")
    position, _ = swe.calc_ut(2451545.0, swe.HARRINGTON, flags)
    assert len(position) == 6
    assert all(math.isfinite(value) for value in position)


@pytest.mark.unit
@pytest.mark.parametrize("body_id", REVIEWED_FICTITIOUS_BODIES)
def test_reviewed_fictitious_body_returns_finite_speed_state(body_id: int) -> None:
    """Every independently transcribed orbit enters the normal output pipeline."""
    swe.set_calc_mode("skyfield")
    position, _ = swe.calc_ut(2451545.0, body_id, FLG_SPEED)
    assert len(position) == 6
    assert all(math.isfinite(value) for value in position)


@pytest.mark.unit
@pytest.mark.parametrize("body_id,name", UNVERIFIED_FICTITIOUS_BODIES)
def test_unverified_fictitious_body_is_explicitly_unavailable(
    body_id: int, name: str
) -> None:
    """Recognised IDs without independently reviewed elements fail closed."""
    swe.set_calc_mode("skyfield")
    with pytest.raises(UnknownBodyError) as raised:
        swe.calc_ut(2451545.0, body_id, FLG_SPEED)
    assert raised.value.body_id == body_id, name
