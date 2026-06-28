# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for libephemeris.extinction band/edge selections."""

from __future__ import annotations

from libephemeris import extinction as ext


class TestAirmass:
    def test_secant_near_horizon_caps(self):
        # zenith >= 85 deg -> capped airmass
        assert ext.calc_airmass(4.0, method="secant") == 40.0

    def test_rozenberg_method(self):
        am = ext.calc_airmass(30.0, method="rozenberg")
        assert am > 1.0

    def test_unknown_method_defaults_to_kasten_young(self):
        # The else branch delegates to the kasten_young computation.
        assert ext.calc_airmass(30.0, method="not-a-real-method") == ext.calc_airmass(
            30.0, method="kasten_young"
        )


class TestWavelengthBands:
    def test_ozone_blue_and_red_bands(self):
        assert ext.calc_ozone_coefficient(450.0) == 0.04  # blue band
        assert ext.calc_ozone_coefficient(800.0) == 0.008  # red/IR band

    def test_water_vapor_red_band(self):
        # 600 <= wavelength < 800 -> 0.005 * humidity_fraction
        val = ext.calc_water_vapor_coefficient(
            humidity_percent=50.0, wavelength_nm=700.0
        )
        assert val == 0.005 * 0.5


class TestSimpleExtinctionEdges:
    def test_near_horizon_practical_maximum(self):
        # altitude > 0 but sin(alt) <= 0.001 -> k * 40
        assert ext.calc_simple_extinction(0.01, k=0.28) == 0.28 * 40.0


class TestAltitudeFactor:
    def test_low_altitude_returns_horizon_factor(self):
        # zenith >= 85 deg -> brighter near-horizon factor
        assert ext._calc_altitude_factor(3.0) == 1.5


class TestContrastThreshold:
    def test_photopic_mid_brightness_branch(self):
        # photopic adaptation with 5 <= sky_brightness < 10 exercises the
        # middle limiting_mag arm.
        assert ext.calc_contrast_threshold(7.0, eye_adaptation="photopic") > 0.0

    def test_photopic_very_bright_branch(self):
        # photopic adaptation with sky_brightness < 5 exercises the daylight arm.
        assert ext.calc_contrast_threshold(2.0, eye_adaptation="photopic") > 0.0


class TestVisibilityThreshold:
    def test_explicit_eye_adaptation_skips_auto_detect(self):
        # Passing eye_adaptation skips the `is None` auto-detection branch.
        result = ext.calc_visibility_threshold(
            object_magnitude=3.0,
            sky_brightness=20.0,
            eye_adaptation="dark",
        )
        assert result.eye_adaptation == "dark"

    def test_zero_sky_flux_gives_infinite_contrast(self, monkeypatch):
        # A large sky_brightness underflows sky_flux (10**-400) to 0.0 -> the
        # contrast falls into the inf branch. Pin the threshold + adaptation so
        # only the sky_flux underflow path is exercised (the auto-detect /
        # threshold helpers would overflow at this brightness).
        monkeypatch.setattr(ext, "calc_contrast_threshold", lambda *a, **k: 0.1)
        result = ext.calc_visibility_threshold(
            object_magnitude=0.0,
            sky_brightness=1000.0,
            eye_adaptation="dark",
        )
        assert result.contrast == float("inf")

    def test_nonpositive_threshold_uses_faint_limit(self, monkeypatch):
        # Force a non-positive contrast threshold to take the faint-limit arm.
        monkeypatch.setattr(ext, "calc_contrast_threshold", lambda *a, **k: 0.0)
        result = ext.calc_visibility_threshold(
            object_magnitude=3.0,
            sky_brightness=20.0,
            eye_adaptation="dark",
        )
        assert result.limiting_magnitude == 20.0 + 10.0
