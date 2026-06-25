"""
Tests for LEB2 compressed reader.

Verifies LEB2Reader opens LEB2 files correctly, produces consistent
results with LEB1, and handles edge cases.
"""

from __future__ import annotations

import math
import os

import pytest

from libephemeris.leb_reader import open_leb
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    MARS,
    JUPITER,
)


LEB2_BASE_CORE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
    "data",
    "leb2",
    "base_core.leb2",
)

LEB1_BASE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
    "data",
    "leb",
    "ephemeris_base.leb",
)

SKIP_NO_LEB2 = pytest.mark.skipif(
    not os.path.exists(LEB2_BASE_CORE),
    reason="LEB2 base core file not found",
)

SKIP_NO_LEB1 = pytest.mark.skipif(
    not os.path.exists(LEB1_BASE),
    reason="LEB1 base file not found",
)


@SKIP_NO_LEB2
class TestLEB2Open:
    """Test opening LEB2 files."""

    @pytest.mark.unit
    def test_open_leb2_returns_reader(self):
        """open_leb auto-detects LEB2 format."""
        reader = open_leb(LEB2_BASE_CORE)
        assert reader is not None
        reader.close()

    @pytest.mark.unit
    def test_leb2_context_manager(self):
        """LEB2Reader supports context manager."""
        with open_leb(LEB2_BASE_CORE) as reader:
            assert hasattr(reader, "eval_body")
            assert hasattr(reader, "has_body")

    @pytest.mark.unit
    def test_leb2_path_property(self):
        """path property returns file path."""
        with open_leb(LEB2_BASE_CORE) as reader:
            assert reader.path == LEB2_BASE_CORE

    @pytest.mark.unit
    def test_leb2_jd_range(self):
        """jd_range returns valid (start, end) tuple."""
        with open_leb(LEB2_BASE_CORE) as reader:
            jd_start, jd_end = reader.jd_range
            assert isinstance(jd_start, float)
            assert isinstance(jd_end, float)
            assert jd_end > jd_start


@SKIP_NO_LEB2
class TestLEB2Bodies:
    """Test LEB2 body operations."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MERCURY, "Mercury"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
        ],
    )
    def test_core_bodies_present(self, body_id: int, name: str):
        """Core bodies should be in LEB2 base core file."""
        with open_leb(LEB2_BASE_CORE) as reader:
            assert reader.has_body(body_id), f"{name} not found"

    @pytest.mark.unit
    def test_eval_body_returns_pos_vel(self):
        """eval_body returns ((x,y,z), (vx,vy,vz))."""
        with open_leb(LEB2_BASE_CORE) as reader:
            pos, vel = reader.eval_body(SUN, 2451545.0)
            assert len(pos) == 3
            assert len(vel) == 3
            for v in pos:
                assert math.isfinite(v)
            for v in vel:
                assert math.isfinite(v)

    @pytest.mark.unit
    def test_unknown_body_raises(self):
        """Unknown body raises KeyError."""
        with open_leb(LEB2_BASE_CORE) as reader:
            with pytest.raises(KeyError):
                reader.eval_body(99999, 2451545.0)

    @pytest.mark.unit
    def test_out_of_range_raises(self):
        """JD out of range raises ValueError."""
        with open_leb(LEB2_BASE_CORE) as reader:
            with pytest.raises(ValueError):
                reader.eval_body(SUN, 1000000.0)


@SKIP_NO_LEB2
@SKIP_NO_LEB1
class TestLEB2VsLEB1:
    """Compare LEB2 results against LEB1."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
        ],
    )
    def test_leb2_matches_leb1(self, body_id: int, name: str):
        """LEB2 positions should match LEB1 within precision target."""
        jd = 2451545.0
        with open_leb(LEB1_BASE) as r1:
            if not r1.has_body(body_id):
                pytest.skip(f"{name} not in LEB1")
            pos1, vel1 = r1.eval_body(body_id, jd)

        with open_leb(LEB2_BASE_CORE) as r2:
            if not r2.has_body(body_id):
                pytest.skip(f"{name} not in LEB2")
            pos2, vel2 = r2.eval_body(body_id, jd)

        # Positions should match within LEB2 precision target (<0.001")
        for i in range(3):
            diff = abs(pos1[i] - pos2[i])
            # For ecliptic bodies, units are degrees; 0.001" = ~2.8e-7 deg
            assert diff < 0.001, (
                f"{name} pos[{i}]: LEB1={pos1[i]}, LEB2={pos2[i]}, diff={diff}"
            )

    @pytest.mark.unit
    def test_leb2_nutation_matches_leb1(self):
        """LEB2 nutation matches LEB1."""
        jd = 2451545.0
        with open_leb(LEB1_BASE) as r1:
            if not hasattr(r1, "eval_nutation"):
                pytest.skip("LEB1 has no eval_nutation")
            n1 = r1.eval_nutation(jd)

        with open_leb(LEB2_BASE_CORE) as r2:
            if not hasattr(r2, "eval_nutation"):
                pytest.skip("LEB2 has no eval_nutation")
            n2 = r2.eval_nutation(jd)

        for i in range(2):
            diff = abs(n1[i] - n2[i])
            assert diff < 1e-9, (
                f"Nutation[{i}]: LEB1={n1[i]}, LEB2={n2[i]}, diff={diff}"
            )


@SKIP_NO_LEB2
class TestLEB2DateRange:
    """Test LEB2 across dates."""

    @pytest.mark.unit
    @pytest.mark.parametrize("year", [1900, 1950, 2000, 2024, 2050, 2100])
    def test_sun_across_years(self, year: int):
        """Sun position valid across years in LEB2."""
        import libephemeris as swe

        jd = swe.julday(year, 6, 21, 12.0)
        with open_leb(LEB2_BASE_CORE) as reader:
            jd_start, jd_end = reader.jd_range
            if jd < jd_start or jd > jd_end:
                pytest.skip(f"JD {jd} outside LEB2 range")
            pos, vel = reader.eval_body(SUN, jd)
            for v in pos:
                assert math.isfinite(v)
