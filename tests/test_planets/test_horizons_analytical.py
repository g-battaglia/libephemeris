"""Tests for Horizons backend analytical bodies and calc_mode switching.

These tests verify the offline analytical paths in the Horizons backend
(mean lunar points and reviewed fictitious orbits), which require no HTTP.
"""

from __future__ import annotations

import math

import pytest
import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MEAN_NODE,
    MEAN_APOG,
    FLG_SWIEPH,
    FLG_SPEED,
    FLG_HELCTR,
)
from libephemeris.horizons_backend import horizons_calc_ut
from libephemeris.exceptions import UnknownBodyError

JD_J2000 = 2451545.0


@pytest.mark.unit
class TestCalcModeSwitch:
    """Test set_calc_mode / get_calc_mode state management."""

    def test_default_mode_is_auto(self):
        """Default calc mode should be 'auto'."""
        mode = swe.get_calc_mode()
        assert mode in ("auto", "skyfield", "leb"), f"Unexpected default mode: {mode}"

    def test_set_skyfield_mode(self):
        """Setting mode to 'skyfield' should work."""
        original = swe.get_calc_mode()
        try:
            swe.set_calc_mode("skyfield")
            assert swe.get_calc_mode() == "skyfield"
        finally:
            swe.set_calc_mode(original)

    def test_set_auto_mode(self):
        """Setting mode to 'auto' should work."""
        original = swe.get_calc_mode()
        try:
            swe.set_calc_mode("auto")
            assert swe.get_calc_mode() == "auto"
        finally:
            swe.set_calc_mode(original)

    def test_invalid_mode_raises(self):
        """Setting an invalid mode should raise ValueError."""
        with pytest.raises(ValueError):
            swe.set_calc_mode("invalid_mode_xyz")

    def test_mode_switch_preserves_results(self):
        """Switching between auto and skyfield should give same results."""
        original = swe.get_calc_mode()
        try:
            swe.set_calc_mode("skyfield")
            res_sky, _ = swe.calc_ut(JD_J2000, SUN, FLG_SWIEPH | FLG_SPEED)

            swe.set_calc_mode("auto")
            res_auto, _ = swe.calc_ut(JD_J2000, SUN, FLG_SWIEPH | FLG_SPEED)

            for i in range(3):
                assert abs(res_sky[i] - res_auto[i]) < 0.001, (
                    f"Component {i}: sky={res_sky[i]}, auto={res_auto[i]}"
                )
        finally:
            swe.set_calc_mode(original)


@pytest.mark.unit
class TestHorizonsAnalyticalBodies:
    """Test Horizons backend analytical bodies (no HTTP required)."""

    def test_mean_node_analytical(self):
        """Mean Node via horizons_calc_ut should match skyfield calc_ut."""
        # horizons_calc_ut with analytical body should work without a client
        result = horizons_calc_ut(None, JD_J2000, 10, FLG_SWIEPH | FLG_SPEED)
        assert len(result) == 2  # (data, retflag)
        data = result[0]
        assert len(data) == 6

        # Should match the Skyfield-based result closely
        ref, _ = swe.calc_ut(JD_J2000, MEAN_NODE, FLG_SWIEPH | FLG_SPEED)
        assert abs(data[0] - ref[0]) < 0.01, (
            f"Mean Node lon: horizons={data[0]}, skyfield={ref[0]}"
        )

    def test_mean_apogee_analytical(self):
        """Mean Apogee (Lilith) via horizons_calc_ut should match calc_ut."""
        result = horizons_calc_ut(None, JD_J2000, 12, FLG_SWIEPH | FLG_SPEED)
        data = result[0]
        assert len(data) == 6

        ref, _ = swe.calc_ut(JD_J2000, MEAN_APOG, FLG_SWIEPH | FLG_SPEED)
        assert abs(data[0] - ref[0]) < 0.01, (
            f"Mean Apogee lon: horizons={data[0]}, skyfield={ref[0]}"
        )

    def test_mean_node_range(self):
        """Mean Node analytical should be in [0, 360) over many dates."""
        for i in range(20):
            jd = JD_J2000 + i * 180.0
            result = horizons_calc_ut(None, jd, 10, FLG_SWIEPH)
            lon = result[0][0]
            assert 0.0 <= lon < 360.0, f"Mean Node at JD {jd}: lon={lon}"

    def test_mean_apogee_range(self):
        """Mean Apogee analytical should be in [0, 360) over many dates."""
        for i in range(20):
            jd = JD_J2000 + i * 180.0
            result = horizons_calc_ut(None, jd, 12, FLG_SWIEPH)
            lon = result[0][0]
            assert 0.0 <= lon < 360.0, f"Mean Apogee at JD {jd}: lon={lon}"

    def test_mean_node_speed_negative(self):
        """Mean Node speed should be negative (retrograde) via analytical path."""
        result = horizons_calc_ut(None, JD_J2000, 10, FLG_SWIEPH | FLG_SPEED)
        speed = result[0][3]
        assert speed < 0, f"Mean Node speed {speed} should be negative (retrograde)"


@pytest.mark.unit
class TestHorizonsFictitiousAnalytical:
    """Exercise the no-HTTP native-center hypothetical paths."""

    @pytest.mark.parametrize("body_id", [*range(40, 56), 57])
    def test_heliocentric_historical_models_follow_provenance_boundary(self, body_id):
        """Reviewed models compute; unverified heliocentric models fail closed."""
        if body_id in {*range(40, 49), *range(50, 56), 57}:
            data, _ = horizons_calc_ut(None, JD_J2000, body_id, FLG_SWIEPH | FLG_HELCTR)
            assert len(data) == 6
            assert 0.0 <= data[0] < 360.0
            return

        with pytest.raises(UnknownBodyError) as raised:
            horizons_calc_ut(None, JD_J2000, body_id, FLG_SWIEPH | FLG_HELCTR)
        assert raised.value.body_id == body_id

    def test_reviewed_white_moon_native_geocentric_model(self):
        """White Moon uses its local published seven-year convention."""
        data, _ = horizons_calc_ut(None, JD_J2000, 56, FLG_SWIEPH | FLG_SPEED)
        assert len(data) == 6
        assert all(math.isfinite(value) for value in data)

    def test_unverified_waldemath_native_geocentric_fails_closed(self):
        """Waldemath has no reviewed built-in model."""
        with pytest.raises(UnknownBodyError) as raised:
            horizons_calc_ut(None, JD_J2000, 58, FLG_SWIEPH)
        assert raised.value.body_id == 58


@pytest.mark.unit
class TestHorizonsUnsupportedFallback:
    """Test that unsupported flags/bodies raise appropriate errors."""

    def test_topocentric_raises(self):
        """FLG_TOPOCTR should raise KeyError in horizons_calc_ut."""
        from libephemeris.constants import FLG_TOPOCTR

        with pytest.raises(KeyError):
            horizons_calc_ut(None, JD_J2000, 0, FLG_TOPOCTR)
