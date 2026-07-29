"""Reference-free contracts for the restored Mardyks compatibility mode."""

from __future__ import annotations

from collections.abc import Iterator
import math
import warnings

import pytest

import libephemeris as le
from libephemeris.ayanamsha_definitions import MARDYKS_DEFINING
from libephemeris.planets import _iau2006_general_precession_deg
from libephemeris.sidereal_epoch import FIXED_EPOCH_LON_OFFSET, FIXED_EPOCH_T0


MODE = le.SIDM_GALALIGN_MARDYKS


@pytest.fixture(autouse=True)
def _reset_state() -> Iterator[None]:
    le.reset_session()
    yield
    le.reset_session()


def _signed_angle(angle: float) -> float:
    return (angle + 180.0) % 360.0 - 180.0


def test_mardyks_is_a_fixed_epoch_frame_mode() -> None:
    # Mode 34 is a full fixed-epoch frame request anchored at the published
    # defining epoch, with the defining 30 degree ayanamsha as a constant
    # ecliptic-longitude offset of the t0 frame.
    assert FIXED_EPOCH_T0[MODE] == MARDYKS_DEFINING[1]
    assert FIXED_EPOCH_LON_OFFSET[MODE] == MARDYKS_DEFINING[0]


@pytest.mark.parametrize("jd_ut", [2415020.0, 2451545.0, 2488070.0])
def test_mardyks_uses_its_restored_precession_anchor(jd_ut: float) -> None:
    le.set_sid_mode(MODE)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        actual = le.get_ayanamsa_ut(jd_ut)
    j2000_value = le.get_ayanamsa_ut(2451545.0)
    assert math.isfinite(actual)
    # Published anchor: exactly 30 degrees at the September equinox 1998
    # (Mardyks, Sacred Astronomy 1991), as a fixed-epoch precession frame.
    mardyks_value, mardyks_epoch = MARDYKS_DEFINING
    expected_j2000 = mardyks_value - _iau2006_general_precession_deg(mardyks_epoch)
    assert _signed_angle(j2000_value - expected_j2000) == pytest.approx(0.0, abs=1e-5)
    assert not caught


@pytest.mark.parametrize(
    "channel",
    [
        0,
        le.FLG_XYZ,
        le.FLG_EQUATORIAL,
        le.FLG_XYZ | le.FLG_EQUATORIAL,
    ],
)
def test_mardyks_is_distinct_from_fagan_bradley(channel: int) -> None:
    flags = le.FLG_SPEED | le.FLG_SIDEREAL | channel
    le.set_sid_mode(MODE)
    actual = le.calc_ut(2432000.5, le.MARS, flags)
    le.set_sid_mode(le.SIDM_FAGAN_BRADLEY)
    fagan = le.calc_ut(2432000.5, le.MARS, flags)
    # Every representation differs: the fixed-epoch frame projection
    # (equatorial of t0, offset t0 ecliptic) is not a scalar shift of the
    # Fagan/Bradley of-date output on any channel.
    assert actual[0] != pytest.approx(fagan[0], abs=1e-5)
    # Frame modes echo FLG_NONUT on top of the caller's flags.
    assert actual[1] == fagan[1] | le.FLG_NONUT


def test_context_and_module_use_the_same_definition() -> None:
    context = le.EphemerisContext()
    context.set_sid_mode(MODE)
    context_result = context.calc_ut(
        2460310.5, le.JUPITER, le.FLG_SPEED | le.FLG_SIDEREAL
    )
    le.set_sid_mode(MODE)
    module_result = le.calc_ut(2460310.5, le.JUPITER, le.FLG_SPEED | le.FLG_SIDEREAL)
    assert context_result[0] == pytest.approx(module_result[0], abs=2e-12)
    assert context_result[1] == module_result[1]


def test_leb_and_skyfield_backends_share_the_definition() -> None:
    le.set_sid_mode(MODE)
    flags = le.FLG_SPEED | le.FLG_SIDEREAL | le.FLG_XYZ
    try:
        le.set_calc_mode("skyfield")
        skyfield, _ = le.calc_ut(2451545.0, le.MARS, flags)
        le.set_calc_mode("leb")
        leb, _ = le.calc_ut(2451545.0, le.MARS, flags)
    finally:
        le.set_calc_mode("auto")
    assert leb == pytest.approx(skyfield, abs=5e-6)
