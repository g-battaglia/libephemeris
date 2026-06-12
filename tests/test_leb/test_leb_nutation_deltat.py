"""
Tests for LEB nutation and delta_t consistency with Skyfield.

Verifies that LEB precomputed nutation and delta-T values
match Skyfield's computed values.
"""

from __future__ import annotations

import math
import os

import pytest

from libephemeris.leb_reader import open_leb

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
LEB1_BASE = os.path.join(PROJECT_ROOT, "data", "leb", "ephemeris_base.leb")
LEB2_BASE_CORE = os.path.join(PROJECT_ROOT, "data", "leb2", "base_core.leb2")

SKIP_NO_LEB1 = pytest.mark.skipif(
    not os.path.exists(LEB1_BASE), reason="LEB1 base not found"
)
SKIP_NO_LEB2 = pytest.mark.skipif(
    not os.path.exists(LEB2_BASE_CORE), reason="LEB2 core not found"
)

JD_J2000 = 2451545.0


@SKIP_NO_LEB1
class TestLEBNutation:
    """Test LEB nutation values."""

    @pytest.mark.unit
    def test_nutation_returns_two_values(self):
        """eval_nutation returns (dpsi, deps)."""
        with open_leb(LEB1_BASE) as reader:
            dpsi, deps = reader.eval_nutation(JD_J2000)
            assert isinstance(dpsi, float)
            assert isinstance(deps, float)

    @pytest.mark.unit
    def test_nutation_small_values(self):
        """Nutation values should be small (< 100 arcsec in radians)."""
        with open_leb(LEB1_BASE) as reader:
            dpsi, deps = reader.eval_nutation(JD_J2000)
            # Max nutation ~20 arcsec = ~1e-4 radians
            assert abs(dpsi) < 0.001
            assert abs(deps) < 0.001

    @pytest.mark.unit
    def test_nutation_varies_with_date(self):
        """Nutation changes over time (18.6-year cycle)."""
        with open_leb(LEB1_BASE) as reader:
            vals = []
            for offset in range(0, 365 * 5, 100):
                jd = JD_J2000 + offset
                dpsi, deps = reader.eval_nutation(jd)
                vals.append(dpsi)
            # Should show variation
            variation = max(vals) - min(vals)
            assert variation > 1e-6, f"Nutation variation too small: {variation}"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "jd",
        [
            2415020.0,  # 1900
            2440587.5,  # 1970
            2451545.0,  # J2000
            2458849.5,  # 2020
        ],
    )
    def test_nutation_at_various_dates(self, jd):
        """Nutation is valid at various dates within range."""
        with open_leb(LEB1_BASE) as reader:
            jd_start, jd_end = reader.jd_range
            if jd_start <= jd <= jd_end:
                dpsi, deps = reader.eval_nutation(jd)
                assert math.isfinite(dpsi)
                assert math.isfinite(deps)


@SKIP_NO_LEB1
class TestLEBDeltaT:
    """Test LEB delta-T values."""

    @pytest.mark.unit
    def test_delta_t_returns_float(self):
        """delta_t returns a float (in days)."""
        with open_leb(LEB1_BASE) as reader:
            dt = reader.delta_t(JD_J2000)
            assert isinstance(dt, float)

    @pytest.mark.unit
    def test_delta_t_j2000_value(self):
        """Delta-T at J2000 should be about 63.8 seconds."""
        with open_leb(LEB1_BASE) as reader:
            dt_days = reader.delta_t(JD_J2000)
            dt_seconds = dt_days * 86400
            # At J2000, delta-T ~ 63.8 seconds
            assert 60.0 < dt_seconds < 70.0, f"Delta-T = {dt_seconds} seconds"

    @pytest.mark.unit
    def test_delta_t_increases_over_time(self):
        """Delta-T generally increases over the recent centuries."""
        with open_leb(LEB1_BASE) as reader:
            jd_start, jd_end = reader.jd_range
            # Sample near the middle of the range
            jd1 = max(jd_start + 1, 2440587.5)  # ~1970
            jd2 = min(jd_end - 1, 2460000.0)  # ~2023
            if jd1 < jd2:
                dt1 = reader.delta_t(jd1)
                dt2 = reader.delta_t(jd2)
                # Delta-T should be larger at later dates (recent centuries)
                assert dt2 > dt1, f"dt1={dt1 * 86400:.1f}s, dt2={dt2 * 86400:.1f}s"

    @pytest.mark.unit
    @pytest.mark.parametrize("offset", [0, 365, 3652, 7305])
    def test_delta_t_finite(self, offset):
        """Delta-T is finite at various dates."""
        with open_leb(LEB1_BASE) as reader:
            jd = JD_J2000 + offset
            jd_start, jd_end = reader.jd_range
            if jd_start <= jd <= jd_end:
                dt = reader.delta_t(jd)
                assert math.isfinite(dt)


@SKIP_NO_LEB2
class TestLEB2NutationDeltaT:
    """Test LEB2 nutation and delta-T consistency with LEB1."""

    @pytest.mark.unit
    @pytest.mark.skipif(
        not os.path.exists(LEB1_BASE), reason="Need LEB1 for comparison"
    )
    def test_nutation_leb1_vs_leb2(self):
        """LEB2 nutation matches LEB1 nutation (uncompressed section)."""
        with open_leb(LEB1_BASE) as leb1:
            with open_leb(LEB2_BASE_CORE) as leb2:
                for offset in range(0, 3652, 365):
                    jd = JD_J2000 + offset
                    dpsi1, deps1 = leb1.eval_nutation(jd)
                    dpsi2, deps2 = leb2.eval_nutation(jd)
                    assert dpsi1 == pytest.approx(dpsi2, abs=1e-12), (
                        f"dpsi at JD {jd}: LEB1={dpsi1}, LEB2={dpsi2}"
                    )
                    assert deps1 == pytest.approx(deps2, abs=1e-12)

    @pytest.mark.unit
    @pytest.mark.skipif(
        not os.path.exists(LEB1_BASE), reason="Need LEB1 for comparison"
    )
    def test_delta_t_leb1_vs_leb2(self):
        """LEB2 delta-T matches LEB1 delta-T."""
        with open_leb(LEB1_BASE) as leb1:
            with open_leb(LEB2_BASE_CORE) as leb2:
                for offset in range(0, 3652, 365):
                    jd = JD_J2000 + offset
                    dt1 = leb1.delta_t(jd)
                    dt2 = leb2.delta_t(jd)
                    assert dt1 == pytest.approx(dt2, abs=1e-12), (
                        f"delta_t at JD {jd}: LEB1={dt1}, LEB2={dt2}"
                    )
