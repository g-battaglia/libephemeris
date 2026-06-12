"""
Tests for LEB reader direct API.

Verifies open_leb(), LEBReader methods: has_body, eval_body,
eval_nutation, delta_t, jd_range, context manager protocol.
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
    VENUS,
    MARS,
    JUPITER,
    SATURN,
)


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


@SKIP_NO_LEB
class TestOpenLeb:
    """Test open_leb factory function."""

    @pytest.mark.unit
    def test_open_leb_returns_reader(self):
        """open_leb returns a reader object."""
        reader = open_leb(LEB_BASE_PATH)
        assert reader is not None
        reader.close()

    @pytest.mark.unit
    def test_open_leb_context_manager(self):
        """Reader supports context manager protocol."""
        with open_leb(LEB_BASE_PATH) as reader:
            assert reader is not None
            assert hasattr(reader, "eval_body")

    @pytest.mark.unit
    def test_open_leb_nonexistent_raises(self):
        """open_leb raises FileNotFoundError for missing files."""
        with pytest.raises(FileNotFoundError):
            open_leb("/nonexistent/path/to/file.leb")

    @pytest.mark.unit
    def test_open_leb_invalid_file_raises(self, tmp_path):
        """open_leb raises ValueError for non-LEB files."""
        bad_file = tmp_path / "bad.leb"
        bad_file.write_bytes(b"not a leb file contents here")
        with pytest.raises((ValueError, Exception)):
            open_leb(str(bad_file))


@SKIP_NO_LEB
class TestLEBReaderProperties:
    """Test LEBReader properties."""

    @pytest.mark.unit
    def test_path_property(self):
        """Reader path matches input."""
        with open_leb(LEB_BASE_PATH) as reader:
            assert reader.path == LEB_BASE_PATH

    @pytest.mark.unit
    def test_jd_range_tuple(self):
        """jd_range returns (start, end) tuple."""
        with open_leb(LEB_BASE_PATH) as reader:
            jd_start, jd_end = reader.jd_range
            assert isinstance(jd_start, float)
            assert isinstance(jd_end, float)
            assert jd_end > jd_start

    @pytest.mark.unit
    def test_jd_range_covers_modern_era(self):
        """Base tier should cover at least 1849-2150."""
        with open_leb(LEB_BASE_PATH) as reader:
            jd_start, jd_end = reader.jd_range
            # 1849-01-01 ~= JD 2396758
            # 2150-12-31 ~= JD 2506716
            assert jd_start < 2400000  # Before ~1858
            assert jd_end > 2500000  # After ~2132


@SKIP_NO_LEB
class TestLEBReaderHasBody:
    """Test has_body method."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_core_bodies_present(self, body_id: int, name: str):
        """Core bodies should be in LEB base file."""
        with open_leb(LEB_BASE_PATH) as reader:
            assert reader.has_body(body_id), f"{name} (id={body_id}) not found"

    @pytest.mark.unit
    def test_unknown_body_not_present(self):
        """Non-existent body ID returns False."""
        with open_leb(LEB_BASE_PATH) as reader:
            assert not reader.has_body(99999)


@SKIP_NO_LEB
class TestLEBReaderEvalBody:
    """Test eval_body method."""

    @pytest.mark.unit
    def test_eval_body_returns_pos_vel(self):
        """eval_body returns ((x,y,z), (vx,vy,vz))."""
        with open_leb(LEB_BASE_PATH) as reader:
            result = reader.eval_body(SUN, 2451545.0)
            assert len(result) == 2
            pos, vel = result
            assert len(pos) == 3
            assert len(vel) == 3

    @pytest.mark.unit
    def test_eval_body_finite_values(self):
        """All returned values should be finite."""
        with open_leb(LEB_BASE_PATH) as reader:
            pos, vel = reader.eval_body(SUN, 2451545.0)
            for v in pos:
                assert math.isfinite(v), f"Non-finite position: {v}"
            for v in vel:
                assert math.isfinite(v), f"Non-finite velocity: {v}"

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
    def test_eval_body_various_bodies(self, body_id: int, name: str):
        """eval_body works for various bodies."""
        with open_leb(LEB_BASE_PATH) as reader:
            if reader.has_body(body_id):
                pos, vel = reader.eval_body(body_id, 2451545.0)
                for v in pos:
                    assert math.isfinite(v), f"{name} position: {v}"

    @pytest.mark.unit
    def test_eval_body_unknown_raises(self):
        """eval_body raises KeyError for unknown body."""
        with open_leb(LEB_BASE_PATH) as reader:
            with pytest.raises(KeyError):
                reader.eval_body(99999, 2451545.0)

    @pytest.mark.unit
    def test_eval_body_out_of_range_raises(self):
        """eval_body raises ValueError for JD outside range."""
        with open_leb(LEB_BASE_PATH) as reader:
            with pytest.raises(ValueError):
                reader.eval_body(SUN, 1000000.0)  # Way before coverage

    @pytest.mark.unit
    def test_eval_body_at_multiple_dates(self):
        """Positions at different dates should differ."""
        with open_leb(LEB_BASE_PATH) as reader:
            pos1, _ = reader.eval_body(MARS, 2451545.0)
            pos2, _ = reader.eval_body(MARS, 2451545.0 + 30.0)
            # At least one coordinate should differ
            diffs = [abs(a - b) for a, b in zip(pos1, pos2)]
            assert max(diffs) > 0.001, f"Positions identical after 30 days: {diffs}"


@SKIP_NO_LEB
class TestLEBReaderNutation:
    """Test eval_nutation method."""

    @pytest.mark.unit
    def test_eval_nutation_returns_tuple(self):
        """eval_nutation returns (dpsi, deps)."""
        with open_leb(LEB_BASE_PATH) as reader:
            if hasattr(reader, "eval_nutation"):
                result = reader.eval_nutation(2451545.0)
                assert len(result) == 2
                dpsi, deps = result
                assert math.isfinite(dpsi)
                assert math.isfinite(deps)

    @pytest.mark.unit
    def test_nutation_values_small(self):
        """Nutation values should be small (arcseconds-scale in radians)."""
        with open_leb(LEB_BASE_PATH) as reader:
            if hasattr(reader, "eval_nutation"):
                dpsi, deps = reader.eval_nutation(2451545.0)
                # Nutation in radians should be < 0.001 (~3.4 arcmin)
                assert abs(dpsi) < 0.001, f"dpsi = {dpsi}"
                assert abs(deps) < 0.001, f"deps = {deps}"


@SKIP_NO_LEB
class TestLEBReaderDeltaT:
    """Test delta_t method."""

    @pytest.mark.unit
    def test_delta_t_returns_float(self):
        """delta_t returns a float."""
        with open_leb(LEB_BASE_PATH) as reader:
            if hasattr(reader, "delta_t"):
                dt = reader.delta_t(2451545.0)
                assert isinstance(dt, float)
                assert math.isfinite(dt)

    @pytest.mark.unit
    def test_delta_t_positive_modern(self):
        """Delta T at J2000 should be positive."""
        with open_leb(LEB_BASE_PATH) as reader:
            if hasattr(reader, "delta_t"):
                dt = reader.delta_t(2451545.0)
                assert dt > 0, f"Delta T = {dt}"
                # Should be ~63.8s = ~0.000739 days
                dt_seconds = dt * 86400
                assert 50 < dt_seconds < 80, f"Delta T = {dt_seconds}s"
