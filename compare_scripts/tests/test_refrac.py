"""
Tests for the refrac atmospheric refraction function.

Tests verify that refrac correctly calculates atmospheric refraction,
matching pyswisseph's swe.refrac() behavior.
"""

import pytest
import math
import swisseph as swe
import libephemeris as ephem


class TestRefracBasic:
    """Basic functionality tests for refrac."""

    def test_refrac_exported(self):
        """Test that refrac is exported from the package."""
        assert hasattr(ephem, "refrac")
        assert callable(ephem.refrac)

    def test_refrac_constants_exported(self):
        """Test that refrac constants are exported."""
        assert hasattr(ephem, "TRUE_TO_APP")
        assert hasattr(ephem, "APP_TO_TRUE")
        assert ephem.TRUE_TO_APP == 0
        assert ephem.APP_TO_TRUE == 1

    def test_refrac_returns_float(self):
        """Test that refrac returns a float."""
        result = ephem.refrac(30.0, 1013.25, 15.0, ephem.TRUE_TO_APP)
        assert isinstance(result, float)

    def test_refrac_zero_pressure_returns_input(self):
        """Test that refrac returns input altitude when pressure is 0."""
        result = ephem.refrac(30.0, 0.0, 15.0, ephem.TRUE_TO_APP)
        assert result == 30.0

    def test_refrac_negative_pressure_returns_input(self):
        """Test that refrac returns input altitude when pressure is negative."""
        result = ephem.refrac(30.0, -10.0, 15.0, ephem.TRUE_TO_APP)
        assert result == 30.0


class TestRefracTrueToApp:
    """Tests for TRUE_TO_APP mode (true altitude to apparent altitude)."""

    def test_refrac_true_to_app_returns_higher_altitude(self):
        """Test that apparent altitude is higher than true altitude."""
        true_alt = 0.0
        apparent_alt = ephem.refrac(true_alt, 1013.25, 15.0, ephem.TRUE_TO_APP)
        assert apparent_alt > true_alt, "Apparent altitude should be higher than true"

    def test_refrac_true_to_app_horizon_value(self):
        """Test that at 0° true altitude, apparent altitude is ~0.5°."""
        apparent_alt = ephem.refrac(0.0, 1013.25, 15.0, ephem.TRUE_TO_APP)
        # At horizon, refraction is ~34 arcmin = ~0.57 degrees
        assert 0.4 < apparent_alt < 0.7, (
            f"Apparent altitude at 0° true should be ~0.5°, got {apparent_alt}"
        )

    def test_refrac_increases_with_altitude(self):
        """Test that apparent altitude increases with true altitude."""
        app_0 = ephem.refrac(0.0, 1013.25, 15.0, ephem.TRUE_TO_APP)
        app_15 = ephem.refrac(15.0, 1013.25, 15.0, ephem.TRUE_TO_APP)
        app_45 = ephem.refrac(45.0, 1013.25, 15.0, ephem.TRUE_TO_APP)
        app_90 = ephem.refrac(90.0, 1013.25, 15.0, ephem.TRUE_TO_APP)

        assert app_0 < app_15 < app_45 < app_90

    def test_refrac_at_zenith_near_input(self):
        """Test that refraction at zenith (90°) adds very little."""
        apparent_alt = ephem.refrac(90.0, 1013.25, 15.0, ephem.TRUE_TO_APP)
        # At zenith, refraction should be essentially 0
        assert abs(apparent_alt - 90.0) < 0.01, (
            f"Apparent altitude at zenith should be ~90°, got {apparent_alt}"
        )


class TestRefracAppToTrue:
    """Tests for APP_TO_TRUE mode (apparent altitude to true altitude)."""

    def test_refrac_app_to_true_returns_lower_altitude(self):
        """Test that true altitude is lower than apparent altitude."""
        apparent_alt = 0.5  # Just above horizon
        true_alt = ephem.refrac(apparent_alt, 1013.25, 15.0, ephem.APP_TO_TRUE)
        assert true_alt < apparent_alt, "True altitude should be lower than apparent"

    def test_refrac_app_to_true_near_horizon(self):
        """Test that apparent altitude near horizon returns valid true altitude."""
        # pyswisseph returns input unchanged for apparent altitudes below the
        # horizon refraction threshold (~0.48°)
        true_alt = ephem.refrac(0.0, 1013.25, 15.0, ephem.APP_TO_TRUE)
        # At 0° apparent, pyswisseph returns 0° (input unchanged)
        assert true_alt == 0.0, f"Expected 0.0, got {true_alt}"

        # Above the threshold, we get a smaller result
        true_alt = ephem.refrac(0.5, 1013.25, 15.0, ephem.APP_TO_TRUE)
        assert true_alt < 0.5, "True altitude should be less than apparent"

    def test_refrac_roundtrip(self):
        """Test that true->app->true gives approximately the original altitude."""
        true_alt = 30.0
        # Get apparent altitude
        apparent_alt = ephem.refrac(true_alt, 1013.25, 15.0, ephem.TRUE_TO_APP)

        # Get back to true altitude
        recovered_true = ephem.refrac(apparent_alt, 1013.25, 15.0, ephem.APP_TO_TRUE)

        # Should be close to original (within 0.01 degree)
        assert abs(recovered_true - true_alt) < 0.01, (
            f"Roundtrip: {true_alt} -> {apparent_alt} -> {recovered_true}"
        )


class TestRefracAtmosphericConditions:
    """Tests for refraction under different atmospheric conditions."""

    def test_refrac_higher_pressure_more_refraction(self):
        """Test that higher pressure causes more refraction."""
        app_low_p = ephem.refrac(10.0, 900.0, 15.0, ephem.TRUE_TO_APP)
        app_high_p = ephem.refrac(10.0, 1100.0, 15.0, ephem.TRUE_TO_APP)

        # Higher pressure -> more refraction -> higher apparent altitude
        assert app_high_p > app_low_p, (
            f"Higher pressure should cause higher apparent alt: {app_high_p} vs {app_low_p}"
        )

    def test_refrac_lower_temp_more_refraction(self):
        """Test that lower temperature causes more refraction."""
        app_cold = ephem.refrac(10.0, 1013.25, -10.0, ephem.TRUE_TO_APP)
        app_warm = ephem.refrac(10.0, 1013.25, 30.0, ephem.TRUE_TO_APP)

        # Colder air -> more refraction -> higher apparent altitude
        assert app_cold > app_warm, (
            f"Colder air should cause higher apparent alt: {app_cold} vs {app_warm}"
        )

    def test_refrac_standard_conditions(self):
        """Test refraction at standard conditions matches expected values."""
        # At 10° altitude, standard conditions, refraction is about 5.3 arcmin
        apparent_alt = ephem.refrac(10.0, 1013.25, 10.0, ephem.TRUE_TO_APP)
        refraction = apparent_alt - 10.0
        # 5.3 arcmin ≈ 0.088 degrees
        assert 0.05 < refraction < 0.15, (
            f"Refraction at 10° should be ~0.09 deg, got {refraction}"
        )


class TestRefracVsSwisseph:
    """Comparison tests with pyswisseph's swe.refrac()."""

    @pytest.mark.parametrize(
        "altitude",
        [0.0, 5.0, 10.0, 15.0, 20.0, 30.0, 45.0, 60.0, 75.0, 90.0],
    )
    def test_refrac_true_to_app_matches_swisseph(self, altitude):
        """Test TRUE_TO_APP matches pyswisseph."""
        pressure = 1013.25
        temp = 15.0

        result_lib = ephem.refrac(altitude, pressure, temp, ephem.TRUE_TO_APP)
        result_swe = swe.refrac(altitude, pressure, temp, swe.TRUE_TO_APP)

        # Allow 10% tolerance for small refraction values, absolute for larger
        if result_swe > 0.01:
            tolerance = result_swe * 0.1
        else:
            tolerance = 0.01

        assert abs(result_lib - result_swe) < tolerance, (
            f"At alt={altitude}°: lib={result_lib:.6f}, swe={result_swe:.6f}"
        )

    @pytest.mark.parametrize(
        "altitude",
        [0.0, 5.0, 10.0, 15.0, 20.0, 30.0, 45.0, 60.0, 75.0, 90.0],
    )
    def test_refrac_app_to_true_matches_swisseph(self, altitude):
        """Test APP_TO_TRUE matches pyswisseph."""
        pressure = 1013.25
        temp = 15.0

        result_lib = ephem.refrac(altitude, pressure, temp, ephem.APP_TO_TRUE)
        result_swe = swe.refrac(altitude, pressure, temp, swe.APP_TO_TRUE)

        # Allow tolerance based on magnitude
        if abs(result_swe) > 0.01:
            tolerance = abs(result_swe) * 0.1
        else:
            tolerance = 0.01

        assert abs(result_lib - result_swe) < tolerance, (
            f"At alt={altitude}°: lib={result_lib:.6f}, swe={result_swe:.6f}"
        )

    @pytest.mark.parametrize(
        "pressure,temp",
        [
            (1013.25, 15.0),  # Standard
            (900.0, 10.0),  # Low pressure
            (1100.0, 20.0),  # High pressure
            (1013.25, -10.0),  # Cold
            (1013.25, 35.0),  # Hot
        ],
    )
    def test_refrac_conditions_match_swisseph(self, pressure, temp):
        """Test refraction under various conditions matches pyswisseph."""
        altitude = 10.0

        result_lib = ephem.refrac(altitude, pressure, temp, ephem.TRUE_TO_APP)
        result_swe = swe.refrac(altitude, pressure, temp, swe.TRUE_TO_APP)

        tolerance = result_swe * 0.1 if result_swe > 0.01 else 0.01

        assert abs(result_lib - result_swe) < tolerance, (
            f"Mismatch at P={pressure}, T={temp}: lib={result_lib:.6f}, swe={result_swe:.6f}"
        )


class TestRefracEdgeCases:
    """Edge case tests for refrac."""

    def test_refrac_negative_altitude(self):
        """Test refraction for negative altitudes (below horizon).

        For true altitudes where refraction cannot bring the object above
        the horizon, the input altitude is returned unchanged (matching
        pyswisseph behavior). For altitudes just below 0° (e.g. -0.5°),
        refraction CAN bring the object above the horizon, so a corrected
        value is returned.
        """
        # At -1°, refraction (~0.5°) can't bring it above 0° -> unchanged
        apparent_alt = ephem.refrac(-1.0, 1013.25, 15.0, ephem.TRUE_TO_APP)
        assert isinstance(apparent_alt, float)
        assert apparent_alt == -1.0

        # At -0.5°, refraction (~0.55°) CAN bring it above 0° -> corrected
        apparent_half = ephem.refrac(-0.5, 1013.25, 15.0, ephem.TRUE_TO_APP)
        assert apparent_half > -0.5

    def test_refrac_very_negative_altitude(self):
        """Test refraction for very negative altitudes."""
        # At -5 degrees, we're in extrapolation territory
        apparent_alt = ephem.refrac(-5.0, 1013.25, 15.0, ephem.TRUE_TO_APP)
        assert isinstance(apparent_alt, float)

    def test_refrac_extreme_pressure(self):
        """Test refraction at extreme pressure values."""
        # Very high altitude (low pressure)
        app_high_alt = ephem.refrac(10.0, 500.0, 0.0, ephem.TRUE_TO_APP)
        # Sea level
        app_sea = ephem.refrac(10.0, 1013.25, 15.0, ephem.TRUE_TO_APP)

        # Low pressure -> less refraction -> lower apparent altitude
        assert app_high_alt < app_sea, "Low pressure should give less refraction"

    def test_refrac_extreme_temperature(self):
        """Test refraction at extreme temperatures."""
        # Very cold
        app_cold = ephem.refrac(10.0, 1013.25, -40.0, ephem.TRUE_TO_APP)
        # Very hot
        app_hot = ephem.refrac(10.0, 1013.25, 50.0, ephem.TRUE_TO_APP)

        # Cold -> more refraction -> higher apparent altitude
        assert app_cold > app_hot, "Cold air should give more refraction"
        assert isinstance(app_cold, float)
        assert isinstance(app_hot, float)

    def test_refrac_default_parameters(self):
        """Test that default parameters work correctly."""
        # Only altitude is required
        result = ephem.refrac(30.0)
        assert isinstance(result, float)
        # Should return apparent altitude > 30 (TRUE_TO_APP by default)
        assert result > 30.0

    def test_refrac_handles_near_horizon(self):
        """Test refraction near the horizon (most sensitive region)."""
        altitudes = [-0.5, 0.0, 0.5, 1.0, 2.0]
        previous = None
        for alt in altitudes:
            apparent = ephem.refrac(alt, 1013.25, 15.0, ephem.TRUE_TO_APP)
            assert isinstance(apparent, float)
            assert not math.isnan(apparent), f"Apparent altitude at {alt}° is NaN"
            assert not math.isinf(apparent), f"Apparent altitude at {alt}° is infinite"
            if previous is not None:
                # Apparent altitude should increase with true altitude
                assert apparent > previous, (
                    f"Apparent should increase: at {alt}°={apparent}, previous={previous}"
                )
            previous = apparent
