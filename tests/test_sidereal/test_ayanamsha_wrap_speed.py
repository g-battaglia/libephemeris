"""Regression: sidereal speed at the ayanamsha 0/360 wrap (star-anchored modes).

Star-anchored ayanamshas cross 0 deg on supported dates (SIDM_GALCENT_COCHRANE
does so around year 2225, inside the medium tier). The sidereal speed
correction central-differences the mod-360 ayanamsha at +/-1 s, so for ~2 s
around each crossing the two samples straddle the wrap; without shortest-arc
unwrapping the raw -360 deg delta reported ~1.6e7 deg/day (360 deg / 2 s).
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import (
    VENUS,
    FLG_SPEED,
    FLG_SIDEREAL,
    SIDM_GALCENT_COCHRANE,
)
from libephemeris.exceptions import EphemerisRangeError

# Zero crossing of the SIDM_GALCENT_COCHRANE true ayanamsha (the curve the
# Skyfield speed correction central-differences), bisected to sub-second
# precision against _get_ayanamsa_for_flags.
WRAP_JD = 2533851.698196034


class TestAyanamshaWrapSpeed:
    """No speed spike while the +/-1 s samples straddle the 0/360 wrap."""

    @pytest.mark.unit
    def test_no_speed_spike_at_wrap(self):
        ephem.set_calc_mode("skyfield")
        ephem.set_sid_mode(SIDM_GALCENT_COCHRANE, 2451545.0, 0.0)
        flags = FLG_SPEED | FLG_SIDEREAL
        step = 0.5 / 86400.0
        for k in range(-6, 7):
            jd = WRAP_JD + k * step
            try:
                pos, _ = ephem.calc_ut(jd, VENUS, flags)
            except EphemerisRangeError:
                pytest.skip("active ephemeris tier does not cover year ~2225")
            assert abs(pos[3]) < 10.0, (
                f"sidereal speed spike {pos[3]:+.1f} deg/day at JD {jd:.8f}"
            )
