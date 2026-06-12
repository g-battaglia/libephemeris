"""
Tests for planetary moon calculations.

Planetary moons require a registered satellite SPK kernel.  Without one,
``calc_ut`` raises ``Error`` exactly like the reference API does when its
sepm* satellite ephemeris file is missing.  With a kernel, positions are
apparent geocentric (light-time + annual aberration).

The position tests use the small Jupiter-excerpt kernel that ships with
Skyfield's own test suite (coverage 2015-03-02 .. 2015-03-04) and skip
when it is not present.
"""

from __future__ import annotations

import math
import os

import pytest

import libephemeris as swe
from libephemeris import planetary_moons
from libephemeris.constants import (
    FLG_EQUATORIAL,
    FLG_SPEED,
    FLG_TRUEPOS,
)


# Legacy alias ids (deprecated) and canonical ids (PLMOON_OFFSET + NAIF)
MOON_IO = 9001
MOON_EUROPA = 9002
MOON_GANYMEDE = 9003
MOON_CALLISTO = 9004

GALILEAN_MOONS = [
    (MOON_IO, 9501, "Io"),
    (MOON_EUROPA, 9502, "Europa"),
    (MOON_GANYMEDE, 9503, "Ganymede"),
    (MOON_CALLISTO, 9504, "Callisto"),
]

# Mid-coverage instant of the Skyfield test excerpt (2015-03-03 TT)
JD_EXCERPT = 2457084.5


def _jupiter_excerpt_path() -> str | None:
    """Skyfield's bundled jup310 test excerpt, if installed."""
    try:
        import skyfield
    except ImportError:
        return None
    path = os.path.join(
        os.path.dirname(skyfield.__file__),
        "tests",
        "data",
        "jup310-2015-03-02.bsp",
    )
    return path if os.path.exists(path) else None


@pytest.fixture
def jupiter_kernel():
    """Register the Jupiter excerpt kernel; unregister afterwards."""
    path = _jupiter_excerpt_path()
    if path is None:
        pytest.skip("Skyfield jup310 test excerpt not available")
    planetary_moons.register_moon_spk(path)
    yield path
    planetary_moons.close_moon_kernels()


class TestUnregisteredMoonsRaise:
    """Without a satellite SPK, moon ids raise Error (reference parity)."""

    def setup_method(self):
        planetary_moons.close_moon_kernels()

    @pytest.mark.unit
    @pytest.mark.parametrize("legacy_id,canonical_id,name", GALILEAN_MOONS)
    def test_unregistered_raises(self, legacy_id: int, canonical_id: int, name: str):
        jd = 2451545.0
        with pytest.raises(swe.Error):
            swe.calc_ut(jd, legacy_id, FLG_SPEED)
        with pytest.raises(swe.Error):
            swe.calc_ut(jd, canonical_id, FLG_SPEED)


class TestGalileanMoonPositions:
    """Real positions from the Skyfield test excerpt kernel."""

    @pytest.mark.unit
    @pytest.mark.parametrize("legacy_id,canonical_id,name", GALILEAN_MOONS)
    def test_position_shape_and_floats(
        self, jupiter_kernel, legacy_id: int, canonical_id: int, name: str
    ):
        """6-tuple of finite native floats, plausible Jupiter distance."""
        result, retflag = swe.calc_ut(JD_EXCERPT, canonical_id, FLG_SPEED)
        assert len(result) == 6
        for i, val in enumerate(result):
            assert type(val) is float, (
                f"{name}: result[{i}] is {type(val).__name__}, expected float"
            )
            assert math.isfinite(val)
        # Jupiter system was ~4.4 AU from Earth in March 2015
        assert 3.5 < result[2] < 7.0, f"{name} distance {result[2]:.2f} AU"

    @pytest.mark.unit
    @pytest.mark.parametrize("legacy_id,canonical_id,name", GALILEAN_MOONS)
    def test_legacy_alias_matches_canonical(
        self, jupiter_kernel, legacy_id: int, canonical_id: int, name: str
    ):
        """Deprecated MOON_* ids return the same position as canonical ids."""
        pos_legacy, _ = swe.calc_ut(JD_EXCERPT, legacy_id, FLG_SPEED)
        pos_canonical, _ = swe.calc_ut(JD_EXCERPT, canonical_id, FLG_SPEED)
        assert pos_legacy == pos_canonical

    @pytest.mark.unit
    def test_light_time_applied(self, jupiter_kernel):
        """FLG_TRUEPOS differs from apparent by the light-time swing."""
        apparent, _ = swe.calc_ut(JD_EXCERPT, 9501, 0)
        true_pos, _ = swe.calc_ut(JD_EXCERPT, 9501, FLG_TRUEPOS)
        delta_arcsec = abs((true_pos[0] - apparent[0] + 180) % 360 - 180) * 3600
        # Io's orbital motion over ~37 min of light time is arcsecond-scale
        assert 1.0 < delta_arcsec < 120.0

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "flags",
        [0, FLG_SPEED, FLG_EQUATORIAL, FLG_SPEED | FLG_EQUATORIAL],
    )
    def test_flag_combos(self, jupiter_kernel, flags: int):
        """Flag combinations return finite 6-tuples."""
        result, _ = swe.calc_ut(JD_EXCERPT, 9501, flags)
        assert len(result) == 6
        assert all(math.isfinite(v) for v in result)

    @pytest.mark.unit
    def test_out_of_coverage_raises_range_error(self, jupiter_kernel):
        """Dates outside the excerpt's window raise EphemerisRangeError."""
        with pytest.raises(swe.EphemerisRangeError):
            swe.calc_ut(2451545.0, 9501, 0)

    @pytest.mark.unit
    def test_arbitrary_naif_moon(self, jupiter_kernel):
        """Moons without named constants work via PLMOON_OFFSET + NAIF."""
        # Amalthea (NAIF 505) is in the excerpt but has no MOON_* constant
        result, _ = swe.calc_ut(JD_EXCERPT, 9505, 0)
        assert 3.5 < result[2] < 7.0
