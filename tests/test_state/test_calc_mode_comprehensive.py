"""
Tests for calc_mode switching and backend selection.

Verifies set_calc_mode/get_calc_mode behavior, backend switching,
and that results are consistent across backends.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.state import set_calc_mode, get_calc_mode
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    FLG_SPEED,
)


class TestCalcModeBasic:
    """Test calc_mode get/set basics."""

    @pytest.mark.unit
    def test_get_calc_mode_returns_string(self):
        """get_calc_mode returns a string."""
        mode = get_calc_mode()
        assert isinstance(mode, str)

    @pytest.mark.unit
    def test_get_calc_mode_valid_value(self):
        """get_calc_mode returns a valid mode string."""
        mode = get_calc_mode()
        assert mode in ("auto", "skyfield", "leb", "horizons"), f"mode={mode}"

    @pytest.mark.unit
    def test_set_calc_mode_skyfield(self):
        """set_calc_mode('skyfield') switches to skyfield."""
        set_calc_mode("skyfield")
        assert get_calc_mode() == "skyfield"

    @pytest.mark.unit
    def test_set_calc_mode_auto(self):
        """set_calc_mode('auto') resets to auto."""
        set_calc_mode("auto")
        mode = get_calc_mode()
        # Auto may resolve to leb if LEB files are available
        assert mode in ("auto", "leb"), f"mode={mode}"

    @pytest.mark.unit
    def test_set_calc_mode_invalid_raises(self):
        """Invalid mode string raises ValueError."""
        with pytest.raises(ValueError):
            set_calc_mode("invalid_mode")

    @pytest.mark.unit
    def test_set_calc_mode_none_resets(self):
        """set_calc_mode(None) resets to default."""
        set_calc_mode(None)
        mode = get_calc_mode()
        assert mode in ("auto", "leb", "skyfield"), f"mode={mode}"


class TestSkyfieldBackend:
    """Test calculations in skyfield mode."""

    @pytest.mark.unit
    def test_skyfield_sun_position(self):
        """Sun position computable in skyfield mode."""
        set_calc_mode("skyfield")
        swe.close()
        result, flags = swe.calc_ut(2451545.0, SUN, FLG_SPEED)
        lon = result[0]
        assert 0 <= lon < 360, f"Sun lon={lon}"
        assert math.isfinite(result[1])
        assert result[2] > 0  # distance positive

    @pytest.mark.unit
    def test_skyfield_moon_position(self):
        """Moon position computable in skyfield mode."""
        set_calc_mode("skyfield")
        swe.close()
        result, flags = swe.calc_ut(2451545.0, MOON, FLG_SPEED)
        lon = result[0]
        assert 0 <= lon < 360, f"Moon lon={lon}"
        # Moon distance should be ~0.0025 AU
        assert 0.002 < result[2] < 0.003, f"Moon dist={result[2]}"


class TestBackendConsistency:
    """Test that different backends produce similar results."""

    @pytest.mark.unit
    def test_leb_vs_skyfield_sun(self):
        """LEB and Skyfield Sun positions should agree within 1 arcsec."""
        jd = 2451545.0

        # LEB mode
        set_calc_mode("leb")
        swe.close()
        try:
            leb_result, _ = swe.calc_ut(jd, SUN, FLG_SPEED)
        except Exception:
            pytest.skip("LEB mode not available")

        # Skyfield mode
        set_calc_mode("skyfield")
        swe.close()
        sky_result, _ = swe.calc_ut(jd, SUN, FLG_SPEED)

        lon_diff = abs(leb_result[0] - sky_result[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff

        # Should agree within 1 arcsecond
        assert lon_diff < 1.0 / 3600, f"Sun lon diff: {lon_diff * 3600:.4f} arcsec"

    @pytest.mark.unit
    def test_leb_vs_skyfield_mars(self):
        """LEB and Skyfield Mars positions should agree within 2 arcsec."""
        jd = 2451545.0

        set_calc_mode("leb")
        swe.close()
        try:
            leb_result, _ = swe.calc_ut(jd, MARS, FLG_SPEED)
        except Exception:
            pytest.skip("LEB mode not available")

        set_calc_mode("skyfield")
        swe.close()
        sky_result, _ = swe.calc_ut(jd, MARS, FLG_SPEED)

        lon_diff = abs(leb_result[0] - sky_result[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff

        assert lon_diff < 2.0 / 3600, f"Mars lon diff: {lon_diff * 3600:.4f} arcsec"


class TestCloseAndReset:
    """Test close behavior."""

    @pytest.mark.unit
    def test_close_allows_recalc(self):
        """close followed by calc_ut should work."""
        swe.close()
        result, _ = swe.calc_ut(2451545.0, SUN, FLG_SPEED)
        assert 0 <= result[0] < 360

    @pytest.mark.unit
    def test_multiple_close_safe(self):
        """Multiple close calls should not crash."""
        swe.close()
        swe.close()
        swe.close()
        result, _ = swe.calc_ut(2451545.0, SUN, FLG_SPEED)
        assert 0 <= result[0] < 360
