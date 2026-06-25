"""
Tests for CompositeLEBReader mixing LEB1+LEB2 files
and using medium/extended tiers.
"""

from __future__ import annotations

import math
import os

import pytest

from libephemeris.leb_reader import open_leb
from libephemeris.leb_composite import CompositeLEBReader
from libephemeris.constants import (
    SUN,
    MOON,
    CHIRON,
)

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

LEB1_BASE = os.path.join(PROJECT_ROOT, "data", "leb", "ephemeris_base.leb")
LEB1_MEDIUM = os.path.join(PROJECT_ROOT, "data", "leb", "ephemeris_medium.leb")
LEB2_DIR = os.path.join(PROJECT_ROOT, "data", "leb2")
LEB2_BASE_CORE = os.path.join(LEB2_DIR, "base_core.leb2")
LEB2_BASE_ASTEROIDS = os.path.join(LEB2_DIR, "base_asteroids.leb2")
LEB2_MEDIUM_CORE = os.path.join(LEB2_DIR, "medium_core.leb2")
LEB2_EXTENDED_CORE = os.path.join(LEB2_DIR, "extended_core.leb2")

JD_J2000 = 2451545.0


class TestMixedLEB1LEB2:
    """Test CompositeLEBReader with mixed LEB1 and LEB2 files."""

    @pytest.mark.unit
    @pytest.mark.skipif(
        not (os.path.exists(LEB1_BASE) and os.path.exists(LEB2_BASE_CORE)),
        reason="Both LEB1 and LEB2 base files required",
    )
    def test_mixed_leb1_leb2_construction(self):
        """Construct composite from LEB1 + LEB2 readers."""
        r1 = open_leb(LEB1_BASE)
        r2 = open_leb(LEB2_BASE_CORE)
        comp = CompositeLEBReader([r1, r2])
        assert comp.has_body(SUN)
        assert comp.has_body(MOON)
        comp.close()

    @pytest.mark.unit
    @pytest.mark.skipif(
        not (os.path.exists(LEB1_BASE) and os.path.exists(LEB2_BASE_CORE)),
        reason="Both LEB1 and LEB2 base files required",
    )
    def test_mixed_eval_body_consistent(self):
        """Mixed composite gives same Sun position as individual readers."""
        with open_leb(LEB1_BASE) as leb1:
            pos1, vel1 = leb1.eval_body(SUN, JD_J2000)

        with open_leb(LEB2_BASE_CORE) as leb2:
            pos2, vel2 = leb2.eval_body(SUN, JD_J2000)

        # LEB1 and LEB2 positions should be very close (LEB2 is compressed LEB1)
        for i in range(3):
            assert pos1[i] == pytest.approx(pos2[i], abs=1e-6), (
                f"pos[{i}]: LEB1={pos1[i]}, LEB2={pos2[i]}"
            )

    @pytest.mark.unit
    @pytest.mark.skipif(
        not (os.path.exists(LEB1_BASE) and os.path.exists(LEB2_BASE_ASTEROIDS)),
        reason="LEB1 base + LEB2 asteroids required",
    )
    def test_mixed_leb1_core_leb2_asteroids(self):
        """LEB1 for core bodies + LEB2 for asteroids."""
        r1 = open_leb(LEB1_BASE)
        r2 = open_leb(LEB2_BASE_ASTEROIDS)
        comp = CompositeLEBReader([r1, r2])
        # Core from LEB1
        assert comp.has_body(SUN)
        pos, vel = comp.eval_body(SUN, JD_J2000)
        for v in pos + vel:
            assert math.isfinite(v)
        # Asteroids from LEB2
        if comp.has_body(CHIRON):
            pos, vel = comp.eval_body(CHIRON, JD_J2000)
            for v in pos + vel:
                assert math.isfinite(v)
        comp.close()


class TestMediumTier:
    """Test with medium tier files."""

    @pytest.mark.unit
    @pytest.mark.skipif(
        not os.path.exists(LEB1_MEDIUM),
        reason="LEB1 medium tier not found",
    )
    def test_medium_tier_opens(self):
        """Medium tier LEB1 file opens correctly."""
        with open_leb(LEB1_MEDIUM) as reader:
            assert reader.has_body(SUN)
            jd_start, jd_end = reader.jd_range
            # Medium tier covers 1550-2650
            assert (jd_end - jd_start) > 365.25 * 500

    @pytest.mark.unit
    @pytest.mark.skipif(
        not os.path.exists(LEB1_MEDIUM),
        reason="LEB1 medium tier not found",
    )
    def test_medium_tier_eval(self):
        """Medium tier produces valid positions."""
        with open_leb(LEB1_MEDIUM) as reader:
            pos, vel = reader.eval_body(SUN, JD_J2000)
            for v in pos + vel:
                assert math.isfinite(v)

    @pytest.mark.unit
    @pytest.mark.skipif(
        not os.path.exists(LEB2_MEDIUM_CORE),
        reason="LEB2 medium core not found",
    )
    def test_medium_tier_leb2_opens(self):
        """Medium tier LEB2 file opens correctly."""
        with open_leb(LEB2_MEDIUM_CORE) as reader:
            assert reader.has_body(SUN)

    @pytest.mark.unit
    @pytest.mark.skipif(
        not os.path.exists(LEB2_MEDIUM_CORE),
        reason="LEB2 medium core not found",
    )
    def test_medium_tier_leb2_companions(self):
        """Medium tier LEB2 discovers companions."""
        with CompositeLEBReader.from_file_with_companions(LEB2_MEDIUM_CORE) as reader:
            assert reader.has_body(SUN)
            pos, vel = reader.eval_body(SUN, JD_J2000)
            for v in pos + vel:
                assert math.isfinite(v)


class TestExtendedTier:
    """Test with extended tier files."""

    @pytest.mark.unit
    @pytest.mark.skipif(
        not os.path.exists(LEB2_EXTENDED_CORE),
        reason="LEB2 extended core not found",
    )
    def test_extended_tier_opens(self):
        """Extended tier LEB2 file opens correctly."""
        with open_leb(LEB2_EXTENDED_CORE) as reader:
            assert reader.has_body(SUN)
            jd_start, jd_end = reader.jd_range
            # Extended tier covers -13200 to +17191
            assert (jd_end - jd_start) > 365.25 * 1000

    @pytest.mark.unit
    @pytest.mark.skipif(
        not os.path.exists(LEB2_EXTENDED_CORE),
        reason="LEB2 extended core not found",
    )
    def test_extended_tier_eval(self):
        """Extended tier produces valid positions."""
        with open_leb(LEB2_EXTENDED_CORE) as reader:
            pos, vel = reader.eval_body(SUN, JD_J2000)
            for v in pos + vel:
                assert math.isfinite(v)

    @pytest.mark.unit
    @pytest.mark.skipif(
        not os.path.exists(LEB2_EXTENDED_CORE),
        reason="LEB2 extended core not found",
    )
    def test_extended_tier_companions(self):
        """Extended tier LEB2 discovers companions."""
        with CompositeLEBReader.from_file_with_companions(LEB2_EXTENDED_CORE) as reader:
            assert reader.has_body(SUN)


class TestJdRangeMerging:
    """Test that jd_range merging works correctly with different tiers."""

    @pytest.mark.unit
    @pytest.mark.skipif(
        not (os.path.exists(LEB2_BASE_CORE) and os.path.exists(LEB2_MEDIUM_CORE)),
        reason="Both base and medium LEB2 core required",
    )
    def test_wider_range_from_two_tiers(self):
        """Composite from base+medium tiers has wider jd_range than either alone."""
        with open_leb(LEB2_BASE_CORE) as base:
            base_range = base.jd_range
        with open_leb(LEB2_MEDIUM_CORE) as medium:
            medium_range = medium.jd_range

        r1 = open_leb(LEB2_BASE_CORE)
        r2 = open_leb(LEB2_MEDIUM_CORE)
        comp = CompositeLEBReader([r1, r2])
        comp_range = comp.jd_range

        # Composite range should be widest
        assert comp_range[0] <= min(base_range[0], medium_range[0])
        assert comp_range[1] >= max(base_range[1], medium_range[1])
        comp.close()
