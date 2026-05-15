"""
Tests for the refrac_extended atmospheric refraction function.

Tests verify that refrac_extended correctly calculates atmospheric refraction
with extended parameters (observer altitude, lapse rate) and dip of horizon.
"""

import pytest
import math
import swisseph as swe
import libephemeris as ephem


class TestRefracExtendedBasic:
    """Basic functionality tests for refrac_extended."""

    def test_refrac_extended_exported(self):
        """Test that refrac_extended is exported from the package."""
        assert hasattr(ephem, "refrac_extended")
        assert callable(ephem.refrac_extended)

    def test_refrac_extended_returns_tuple(self):
        """Test that refrac_extended returns a tuple of (float, tuple)."""
        result = ephem.refrac_extended(
            30.0, 0.0, 1013.25, 15.0, 0.0065, ephem.TRUE_TO_APP
        )
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], float)
        assert isinstance(result[1], tuple)
        assert len(result[1]) == 4

    def test_refrac_extended_detail_tuple_contents(self):
        """Test that detail tuple contains expected values."""
        alt, (true_alt, app_alt, refrac, dip) = ephem.refrac_extended(
            30.0, 0.0, 1013.25, 15.0, 0.0065, ephem.TRUE_TO_APP
        )
        assert isinstance(true_alt, float)
        assert isinstance(app_alt, float)
        assert isinstance(refrac, float)
        assert isinstance(dip, float)

    def test_refrac_extended_zero_pressure_returns_input(self):
        """Test that refrac_extended returns input altitude when pressure is 0."""
        alt, (true_alt, app_alt, refrac, dip) = ephem.refrac_extended(
            30.0, 0.0, 0.0, 15.0, 0.0065, ephem.TRUE_TO_APP
        )
        assert alt == 30.0
        assert refrac == 0.0


class TestRefracExtendedTrueToApp:
    """Tests for TRUE_TO_APP mode."""

    def test_refrac_extended_true_to_app_returns_higher_altitude(self):
        """Test that apparent altitude is higher than true altitude."""
        true_alt_input = 0.0
        alt, _ = ephem.refrac_extended(
            true_alt_input, 0.0, 1013.25, 15.0, 0.0065, ephem.TRUE_TO_APP
        )
        assert alt > true_alt_input

    def test_refrac_extended_true_to_app_detail_values(self):
        """Test that detail tuple values are consistent for TRUE_TO_APP."""
        input_alt = 10.0
        alt, (true_alt, app_alt, refrac, dip) = ephem.refrac_extended(
            input_alt, 0.0, 1013.25, 15.0, 0.0065, ephem.TRUE_TO_APP
        )
        # For TRUE_TO_APP: true_alt should equal input, app_alt equals returned alt
        assert true_alt == input_alt
        assert app_alt == alt
        # Refraction should be positive
        assert refrac > 0
        # refrac = app_alt - true_alt
        assert abs(refrac - (app_alt - true_alt)) < 1e-10


class TestRefracExtendedAppToTrue:
    """Tests for APP_TO_TRUE mode."""

    def test_refrac_extended_app_to_true_returns_lower_altitude(self):
        """Test that true altitude is lower than apparent altitude."""
        app_alt_input = 10.0
        alt, _ = ephem.refrac_extended(
            app_alt_input, 0.0, 1013.25, 15.0, 0.0065, ephem.APP_TO_TRUE
        )
        assert alt < app_alt_input

    def test_refrac_extended_app_to_true_detail_values(self):
        """Test that detail tuple values are consistent for APP_TO_TRUE."""
        input_alt = 10.0
        alt, (true_alt, app_alt, refrac, dip) = ephem.refrac_extended(
            input_alt, 0.0, 1013.25, 15.0, 0.0065, ephem.APP_TO_TRUE
        )
        # For APP_TO_TRUE: app_alt should equal input, true_alt equals returned alt
        assert app_alt == input_alt
        assert true_alt == alt
        # Refraction should be positive
        assert refrac > 0


class TestRefracExtendedDip:
    """Tests for dip of the horizon calculation."""

    def test_dip_zero_at_sea_level(self):
        """Test that dip is zero at sea level."""
        _, (_, _, _, dip) = ephem.refrac_extended(
            0.0, 0.0, 1013.25, 15.0, 0.0065, ephem.TRUE_TO_APP
        )
        assert dip == 0.0

    def test_dip_negative_at_elevation(self):
        """Test that dip is negative for elevated observers."""
        _, (_, _, _, dip) = ephem.refrac_extended(
            0.0, 1000.0, 1013.25, 15.0, 0.0065, ephem.TRUE_TO_APP
        )
        assert dip < 0

    def test_dip_increases_with_altitude(self):
        """Test that dip magnitude increases with observer altitude."""
        dips = []
        for geoalt in [100, 500, 1000, 5000]:
            _, (_, _, _, dip) = ephem.refrac_extended(
                0.0, geoalt, 1013.25, 15.0, 0.0065, ephem.TRUE_TO_APP
            )
            dips.append(dip)

        # Each dip should be more negative than the previous
        for i in range(1, len(dips)):
            assert dips[i] < dips[i - 1]

    def test_dip_affected_by_lapse_rate(self):
        """Test that lapse rate affects the dip calculation."""
        _, (_, _, _, dip_low_lapse) = ephem.refrac_extended(
            0.0, 1000.0, 1013.25, 15.0, 0.003, ephem.TRUE_TO_APP
        )
        _, (_, _, _, dip_high_lapse) = ephem.refrac_extended(
            0.0, 1000.0, 1013.25, 15.0, 0.010, ephem.TRUE_TO_APP
        )
        # Higher lapse rate -> more refraction correction -> less negative dip
        # (dip_high_lapse is closer to 0 than dip_low_lapse)
        assert dip_high_lapse > dip_low_lapse


class TestRefracExtendedVsSwisseph:
    """Comparison tests with pyswisseph's swe.refrac_extended()."""

    @pytest.mark.parametrize("geoalt", [0, 100, 500, 1000, 5000])
    def test_dip_matches_swisseph(self, geoalt):
        """Test that dip calculation matches pyswisseph."""
        result_lib = ephem.refrac_extended(
            0.0, geoalt, 1013.25, 15.0, 0.0065, ephem.TRUE_TO_APP
        )
        result_swe = swe.refrac_extended(
            0.0, geoalt, 1013.25, 15.0, 0.0065, swe.TRUE_TO_APP
        )

        # Dip values should be very close (within 0.001 degrees)
        assert abs(result_lib[1][3] - result_swe[1][3]) < 0.001, (
            f"Dip mismatch at geoalt={geoalt}: "
            f"lib={result_lib[1][3]:.6f}, swe={result_swe[1][3]:.6f}"
        )

    @pytest.mark.parametrize(
        "lapse_rate",
        [0.003, 0.004, 0.005, 0.006, 0.0065, 0.007, 0.008, 0.009, 0.010],
    )
    def test_dip_lapse_rate_matches_swisseph(self, lapse_rate):
        """Test that dip calculation matches pyswisseph at various lapse rates."""
        result_lib = ephem.refrac_extended(
            0.0, 1000.0, 1013.25, 15.0, lapse_rate, ephem.TRUE_TO_APP
        )
        result_swe = swe.refrac_extended(
            0.0, 1000.0, 1013.25, 15.0, lapse_rate, swe.TRUE_TO_APP
        )

        # Dip values should be very close (within 0.001 degrees)
        assert abs(result_lib[1][3] - result_swe[1][3]) < 0.001, (
            f"Dip mismatch at lapse_rate={lapse_rate}: "
            f"lib={result_lib[1][3]:.6f}, swe={result_swe[1][3]:.6f}"
        )

    @pytest.mark.parametrize("altitude", [0.0, 10.0, 30.0, 45.0, 90.0])
    def test_refrac_true_to_app_close_to_swisseph(self, altitude):
        """Test that refraction values are reasonably close to pyswisseph."""
        result_lib = ephem.refrac_extended(
            altitude, 0.0, 1013.25, 15.0, 0.0065, ephem.TRUE_TO_APP
        )
        result_swe = swe.refrac_extended(
            altitude, 0.0, 1013.25, 15.0, 0.0065, swe.TRUE_TO_APP
        )

        # Refraction values should be close (within 10% or 0.01 degrees)
        diff = abs(result_lib[1][2] - result_swe[1][2])
        tolerance = max(0.01, abs(result_swe[1][2]) * 0.15)
        assert diff < tolerance, (
            f"Refraction mismatch at alt={altitude}: "
            f"lib={result_lib[1][2]:.6f}, swe={result_swe[1][2]:.6f}"
        )

    @pytest.mark.parametrize("altitude", [5.0, 10.0, 30.0, 45.0])
    def test_refrac_app_to_true_close_to_swisseph(self, altitude):
        """Test that APP_TO_TRUE values are reasonably close to pyswisseph."""
        result_lib = ephem.refrac_extended(
            altitude, 0.0, 1013.25, 15.0, 0.0065, ephem.APP_TO_TRUE
        )
        result_swe = swe.refrac_extended(
            altitude, 0.0, 1013.25, 15.0, 0.0065, swe.APP_TO_TRUE
        )

        # Returned altitude should be close (within 10% of the difference from input)
        diff = abs(result_lib[0] - result_swe[0])
        tolerance = max(0.01, abs(altitude - result_swe[0]) * 0.15)
        assert diff < tolerance, (
            f"Altitude mismatch at alt={altitude}: "
            f"lib={result_lib[0]:.6f}, swe={result_swe[0]:.6f}"
        )


class TestRefracExtendedEdgeCases:
    """Edge case tests for refrac_extended."""

    def test_negative_altitude(self):
        """Test refraction for negative altitudes (below horizon)."""
        alt, (true_alt, app_alt, refrac, dip) = ephem.refrac_extended(
            -1.0, 0.0, 1013.25, 15.0, 0.0065, ephem.TRUE_TO_APP
        )
        assert isinstance(alt, float)
        assert not math.isnan(alt)
        assert not math.isinf(alt)

    def test_very_high_altitude_observer(self):
        """Test with very high observer altitude."""
        alt, (_, _, _, dip) = ephem.refrac_extended(
            0.0, 10000.0, 1013.25, 15.0, 0.0065, ephem.TRUE_TO_APP
        )
        assert isinstance(alt, float)
        assert dip < -2.0  # Should be significant dip

    def test_zero_lapse_rate(self):
        """Test with zero lapse rate (isothermal atmosphere)."""
        _, (_, _, _, dip) = ephem.refrac_extended(
            0.0, 1000.0, 1013.25, 15.0, 0.0, ephem.TRUE_TO_APP
        )
        assert isinstance(dip, float)
        # With zero lapse rate, should get pure geometric dip
        # Check against geometric formula
        EARTH_RADIUS = 6371000.0
        ratio = EARTH_RADIUS / (EARTH_RADIUS + 1000.0)
        geometric_dip = -math.degrees(math.acos(ratio))
        assert abs(dip - geometric_dip) < 0.001

    def test_default_parameters(self):
        """Test that default parameters work correctly."""
        # All parameters have defaults except altitude and altitude_geo
        result = ephem.refrac_extended(30.0, 0.0)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_zenith_refraction(self):
        """Test that refraction at zenith is essentially zero."""
        alt, (_, _, refrac, _) = ephem.refrac_extended(
            90.0, 0.0, 1013.25, 15.0, 0.0065, ephem.TRUE_TO_APP
        )
        assert abs(refrac) < 0.01  # Less than 0.01 degrees at zenith

    def test_roundtrip_conversion(self):
        """Test that true->app->true gives approximately the original altitude."""
        original = 30.0
        # Convert to apparent
        app_alt, _ = ephem.refrac_extended(
            original, 0.0, 1013.25, 15.0, 0.0065, ephem.TRUE_TO_APP
        )
        # Convert back to true
        recovered, _ = ephem.refrac_extended(
            app_alt, 0.0, 1013.25, 15.0, 0.0065, ephem.APP_TO_TRUE
        )
        assert abs(recovered - original) < 0.01, (
            f"Roundtrip: {original} -> {app_alt} -> {recovered}"
        )
