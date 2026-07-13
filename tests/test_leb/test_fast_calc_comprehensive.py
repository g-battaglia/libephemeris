"""
Tests for fast_calc_ut function.

Verifies that fast_calc_ut produces results consistent with
calc_ut for supported bodies and flags.
"""

from __future__ import annotations

import math
import os

import pytest

from libephemeris.fast_calc import fast_calc_ut
from libephemeris.leb_reader import open_leb
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    FLG_SPEED,
    FLG_SIDEREAL,
    FLG_TOPOCTR,
    FLG_XYZ,
    SIDM_J1900,
)
import libephemeris as swe


LEB_BASE_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
    "data",
    "leb",
    "ephemeris_base.leb",
)

SKIP_NO_LEB = pytest.mark.skipif(
    not os.path.exists(LEB_BASE_PATH),
    reason="LEB base file not found",
)

CORE_BODIES = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
]


@SKIP_NO_LEB
class TestFastCalcBasic:
    """Basic fast_calc_ut functionality."""

    @pytest.mark.unit
    def test_returns_tuple_pair(self):
        """fast_calc_ut returns ((6 floats), iflag)."""
        with open_leb(LEB_BASE_PATH) as reader:
            result, iflag = fast_calc_ut(reader, 2451545.0, SUN, FLG_SPEED)
            assert len(result) == 6
            assert isinstance(iflag, int)

    @pytest.mark.unit
    def test_returns_native_floats(self):
        """All result values should be native Python floats."""
        with open_leb(LEB_BASE_PATH) as reader:
            result, _ = fast_calc_ut(reader, 2451545.0, MARS, FLG_SPEED)
            for i, v in enumerate(result):
                assert isinstance(v, (int, float)), f"result[{i}] is {type(v).__name__}"

    @pytest.mark.unit
    def test_all_values_finite(self):
        """All result values should be finite."""
        with open_leb(LEB_BASE_PATH) as reader:
            result, _ = fast_calc_ut(reader, 2451545.0, MARS, FLG_SPEED)
            for i, v in enumerate(result):
                assert math.isfinite(v), f"result[{i}] = {v}"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", CORE_BODIES)
    def test_various_bodies(self, body_id: int, name: str):
        """fast_calc_ut works for core bodies."""
        with open_leb(LEB_BASE_PATH) as reader:
            if reader.has_body(body_id):
                result, _ = fast_calc_ut(reader, 2451545.0, body_id, FLG_SPEED)
                lon = result[0]
                assert 0 <= lon < 360, f"{name}: lon={lon}"


@SKIP_NO_LEB
class TestFastCalcConsistency:
    """Test fast_calc_ut consistency with calc_ut."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", CORE_BODIES)
    def test_matches_calc_ut(self, body_id: int, name: str):
        """fast_calc_ut should match calc_ut within tolerance."""
        jd = 2451545.0
        flags = FLG_SPEED

        with open_leb(LEB_BASE_PATH) as reader:
            if not reader.has_body(body_id):
                pytest.skip(f"{name} not in LEB file")

            fast_result, _ = fast_calc_ut(reader, jd, body_id, flags)

            # Compare with calc_ut (LEB mode)
            swe.close()
            swe_result, _ = swe.calc_ut(jd, body_id, flags)

            # Longitude within 1 arcsecond
            lon_diff = abs(fast_result[0] - swe_result[0])
            if lon_diff > 180:
                lon_diff = 360 - lon_diff
            assert lon_diff < 1.0 / 3600, (
                f"{name}: lon diff {lon_diff * 3600:.3f} arcsec"
            )

    @pytest.mark.unit
    def test_speed_flag_produces_speeds(self):
        """With FLG_SPEED, speed components should be nonzero."""
        with open_leb(LEB_BASE_PATH) as reader:
            result, _ = fast_calc_ut(reader, 2451545.0, MARS, FLG_SPEED)
            speed_lon = result[3]
            # Mars speed should be nonzero
            assert abs(speed_lon) > 0.001, f"Mars speed: {speed_lon}"


@SKIP_NO_LEB
class TestFastCalcSidereal:
    """Test fast_calc_ut with sidereal mode."""

    @pytest.mark.unit
    def test_sidereal_differs_from_tropical(self):
        """Sidereal result should differ from tropical by ~ayanamsha."""
        with open_leb(LEB_BASE_PATH) as reader:
            trop, _ = fast_calc_ut(reader, 2451545.0, MARS, FLG_SPEED)
            sid, _ = fast_calc_ut(
                reader,
                2451545.0,
                MARS,
                FLG_SPEED | FLG_SIDEREAL,
                sid_mode=SIDM_J1900,
            )

            diff = abs(trop[0] - sid[0])
            if diff > 180:
                diff = 360 - diff
            assert diff > 1.0, f"Trop-Sid diff: {diff:.2f} deg"


@SKIP_NO_LEB
class TestFastCalcUnsupported:
    """Test fast_calc_ut with unsupported flags."""

    @pytest.mark.unit
    def test_topocentric_requires_topo(self):
        """FLG_TOPOCTR without set_topo() raises ValueError."""
        from libephemeris import state

        saved_topo = state._TOPO
        state._TOPO = None
        try:
            with open_leb(LEB_BASE_PATH) as reader:
                with pytest.raises(ValueError, match="set_topo"):
                    fast_calc_ut(reader, 2451545.0, SUN, FLG_TOPOCTR)
        finally:
            state._TOPO = saved_topo

    @pytest.mark.unit
    def test_xyz_returns_cartesian(self):
        """FLG_XYZ returns Cartesian coordinates (~1 AU for Sun)."""
        import math

        with open_leb(LEB_BASE_PATH) as reader:
            result, _ = fast_calc_ut(reader, 2451545.0, SUN, FLG_XYZ)
            r = math.sqrt(result[0] ** 2 + result[1] ** 2 + result[2] ** 2)
            assert 0.98 < r < 1.02

    @pytest.mark.unit
    def test_unknown_body_raises(self):
        """Unknown body ID should raise KeyError."""
        with open_leb(LEB_BASE_PATH) as reader:
            with pytest.raises(KeyError):
                fast_calc_ut(reader, 2451545.0, 99999, FLG_SPEED)

    @pytest.mark.unit
    def test_out_of_range_raises(self):
        """JD out of LEB range should raise ValueError."""
        with open_leb(LEB_BASE_PATH) as reader:
            with pytest.raises(ValueError):
                fast_calc_ut(reader, 1000000.0, SUN, FLG_SPEED)


@SKIP_NO_LEB
class TestFastCalcDateRange:
    """Test fast_calc_ut across dates."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "year",
        [1900, 1950, 2000, 2024, 2050, 2100],
    )
    def test_sun_across_years(self, year: int):
        """Sun position valid across years."""
        jd = swe.julday(year, 6, 21, 12.0)
        with open_leb(LEB_BASE_PATH) as reader:
            result, _ = fast_calc_ut(reader, jd, SUN, FLG_SPEED)
            assert 0 <= result[0] < 360
            assert math.isfinite(result[1])
            assert result[2] > 0  # distance positive
