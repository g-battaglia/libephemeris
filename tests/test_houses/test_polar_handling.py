"""
Comprehensive tests for improved Placidus polar latitude handling.

Tests the PolarCircleError exception, get_polar_latitude_threshold(),
houses_with_fallback(), and related functionality.

Also includes tests for extreme latitude (>80°) edge cases.
"""

import pytest
import math
import libephemeris as ephem
from libephemeris.exceptions import PolarCircleError


class TestPolarCircleError:
    """Test the PolarCircleError exception class."""

    @pytest.mark.unit
    def test_polar_circle_error_is_error_subclass(self):
        """PolarCircleError should be a subclass of Error."""
        assert issubclass(PolarCircleError, ephem.Error)

    @pytest.mark.unit
    def test_polar_circle_error_attributes(self):
        """PolarCircleError should have expected attributes."""
        err = PolarCircleError(
            message="test message",
            latitude=70.0,
            threshold=66.56,
            obliquity=23.44,
            house_system="P",
        )
        assert err.latitude == 70.0
        assert err.threshold == 66.56
        assert err.obliquity == 23.44
        assert err.house_system == "P"
        assert str(err) == "test message"

    @pytest.mark.unit
    def test_polar_circle_error_repr(self):
        """PolarCircleError repr should be informative."""
        err = PolarCircleError(
            message="test",
            latitude=70.0,
            threshold=66.56,
            obliquity=23.44,
            house_system="P",
        )
        repr_str = repr(err)
        assert "PolarCircleError" in repr_str
        assert "70.0" in repr_str
        assert "66.56" in repr_str

    @pytest.mark.unit
    def test_polar_circle_error_can_be_caught_as_error(self):
        """PolarCircleError should be catchable as ephem.Error."""
        jd = 2451545.0

        try:
            ephem.houses(jd, 70.0, 0.0, ord("P"))
            pytest.fail("Should have raised PolarCircleError")
        except ephem.Error as e:
            # Should be caught as Error
            assert isinstance(e, PolarCircleError)


class TestGetPolarLatitudeThreshold:
    """Test the get_polar_latitude_threshold() helper function."""

    @pytest.mark.unit
    def test_default_obliquity(self):
        """Default obliquity should give approximately 66.56° threshold."""
        threshold = ephem.get_polar_latitude_threshold()
        assert abs(threshold - 66.56) < 0.01

    @pytest.mark.unit
    def test_custom_obliquity(self):
        """Custom obliquity should give correct threshold."""
        # For obliquity of 23.5°, threshold should be 66.5°
        threshold = ephem.get_polar_latitude_threshold(23.5)
        assert threshold == 66.5

    @pytest.mark.unit
    def test_zero_obliquity(self):
        """Zero obliquity should give 90° threshold."""
        threshold = ephem.get_polar_latitude_threshold(0.0)
        assert threshold == 90.0

    @pytest.mark.unit
    def test_max_obliquity(self):
        """Obliquity of 90° should give 0° threshold."""
        threshold = ephem.get_polar_latitude_threshold(90.0)
        assert threshold == 0.0


class TestPolarCircleErrorDetails:
    """Test that PolarCircleError provides helpful information."""

    @pytest.mark.unit
    def test_error_includes_latitude(self):
        """Error message should include the problematic latitude."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.houses(jd, 70.0, 0.0, ord("P"))

        assert exc_info.value.latitude == 70.0
        assert "70" in str(exc_info.value)

    @pytest.mark.unit
    def test_error_includes_threshold(self):
        """Error message should include the threshold."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.houses(jd, 70.0, 0.0, ord("P"))

        assert exc_info.value.threshold is not None
        assert exc_info.value.threshold > 66.0
        assert exc_info.value.threshold < 67.0

    @pytest.mark.unit
    def test_error_includes_house_system(self):
        """Error message should include the house system."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.houses(jd, 70.0, 0.0, ord("P"))

        assert exc_info.value.house_system == "P"
        assert "Placidus" in str(exc_info.value)

    @pytest.mark.unit
    def test_error_suggests_alternatives(self):
        """Error message should suggest alternative house systems."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.houses(jd, 70.0, 0.0, ord("P"))

        msg = str(exc_info.value).lower()
        assert "porphyry" in msg or "equal" in msg or "whole sign" in msg

    @pytest.mark.unit
    def test_error_mentions_fallback_function(self):
        """Error message should mention houses_with_fallback."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.houses(jd, 70.0, 0.0, ord("P"))

        assert "fallback" in str(exc_info.value).lower()

    @pytest.mark.unit
    @pytest.mark.parametrize("lat", [67.0, 70.0, 80.0, 89.0])
    def test_northern_polar_latitudes(self, lat):
        """Northern polar latitudes should raise detailed error."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.houses(jd, lat, 0.0, ord("P"))

        assert exc_info.value.latitude == lat
        # Should mention Northern
        assert "N" in str(exc_info.value) or "Northern" in str(exc_info.value)

    @pytest.mark.unit
    @pytest.mark.parametrize("lat", [-67.0, -70.0, -80.0, -89.0])
    def test_southern_polar_latitudes(self, lat):
        """Southern polar latitudes should raise detailed error."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.houses(jd, lat, 0.0, ord("P"))

        assert exc_info.value.latitude == lat
        # Should mention Southern
        assert "S" in str(exc_info.value) or "Southern" in str(exc_info.value)


class TestSweHousesWithFallback:
    """Test the houses_with_fallback() convenience function."""

    @pytest.mark.unit
    def test_returns_normal_result_for_non_polar(self):
        """Should return normal result without fallback for non-polar latitudes."""
        jd = 2451545.0

        cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
            jd, 45.0, 0.0, ord("P")
        )

        assert len(cusps) == 12
        assert len(ascmc) == 8
        assert used_fallback is False
        assert warning is None

    @pytest.mark.unit
    def test_falls_back_for_polar_latitude(self):
        """Should fall back to Porphyry for polar latitudes."""
        jd = 2451545.0

        cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
            jd, 70.0, 0.0, ord("P")
        )

        assert len(cusps) == 12
        assert len(ascmc) == 8
        assert used_fallback is True
        assert warning is not None
        assert "Placidus" in warning
        assert "Porphyry" in warning

    @pytest.mark.unit
    def test_warning_includes_threshold(self):
        """Warning message should include the threshold."""
        jd = 2451545.0

        _, _, _, warning = ephem.houses_with_fallback(jd, 70.0, 0.0, ord("P"))

        assert "threshold" in warning.lower()

    @pytest.mark.unit
    def test_custom_fallback_system(self):
        """Should use custom fallback system when specified."""
        jd = 2451545.0

        cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
            jd,
            70.0,
            0.0,
            ord("P"),
            fallback_hsys=ord("E"),  # Equal houses
        )

        assert used_fallback is True
        assert "Equal" in warning

    @pytest.mark.unit
    def test_all_polar_affected_systems(self):
        """Should handle all polar-affected systems."""
        jd = 2451545.0
        lat = 70.0

        for hsys in [ord("P"), ord("K"), ord("G")]:  # Placidus, Koch, Gauquelin
            cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
                jd, lat, 0.0, hsys
            )
            assert used_fallback is True
            assert len(cusps) == 12

    @pytest.mark.unit
    def test_polar_safe_systems_no_fallback(self):
        """Polar-safe systems should not trigger fallback."""
        jd = 2451545.0
        lat = 70.0

        # These systems should work at polar latitudes
        for hsys in [ord("O"), ord("E"), ord("W"), ord("M")]:
            cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
                jd, lat, 0.0, hsys
            )
            assert used_fallback is False
            assert warning is None

    @pytest.mark.unit
    def test_matches_porphyry_when_fallback_used(self):
        """Fallback cusps should match direct Porphyry calculation."""
        jd = 2451545.0
        lat = 70.0

        # Get fallback result
        cusps_fallback, ascmc_fallback, used_fallback, _ = ephem.houses_with_fallback(
            jd, lat, 0.0, ord("P")
        )

        # Get direct Porphyry
        cusps_porphyry, ascmc_porphyry = ephem.houses(jd, lat, 0.0, ord("O"))

        # Should match
        assert used_fallback is True
        for i in range(12):
            assert abs(cusps_fallback[i] - cusps_porphyry[i]) < 0.001

    @pytest.mark.unit
    def test_returns_valid_cusps_at_extreme_latitudes(self):
        """Should return valid cusps even at extreme polar latitudes."""
        jd = 2451545.0

        for lat in [80.0, 85.0, 89.0, -80.0, -85.0, -89.0]:
            cusps, ascmc, _, _ = ephem.houses_with_fallback(jd, lat, 0.0, ord("P"))
            for cusp in cusps:
                assert 0 <= cusp < 360


class TestSweHousesArmcWithFallback:
    """Test the houses_armc_with_fallback() convenience function."""

    @pytest.mark.unit
    def test_returns_normal_result_for_non_polar(self):
        """Should return normal result without fallback for non-polar latitudes."""
        armc = 280.0
        eps = 23.44

        cusps, ascmc, used_fallback, warning = ephem.houses_armc_with_fallback(
            armc, 45.0, eps, ord("P")
        )

        assert len(cusps) == 12
        assert len(ascmc) == 8
        assert used_fallback is False
        assert warning is None

    @pytest.mark.unit
    def test_falls_back_for_polar_latitude(self):
        """Should fall back for polar latitudes."""
        armc = 280.0
        eps = 23.44

        cusps, ascmc, used_fallback, warning = ephem.houses_armc_with_fallback(
            armc, 70.0, eps, ord("P")
        )

        assert used_fallback is True
        assert warning is not None
        assert "Placidus" in warning

    @pytest.mark.unit
    def test_custom_fallback(self):
        """Should use custom fallback when specified."""
        armc = 280.0
        eps = 23.44

        _, _, used_fallback, warning = ephem.houses_armc_with_fallback(
            armc, 70.0, eps, ord("P"), fallback_hsys=ord("W")
        )

        assert used_fallback is True
        assert "Whole Sign" in warning


class TestPolarCircleEdgeCases:
    """Test edge cases at and near the polar circle threshold."""

    @pytest.mark.unit
    def test_just_below_threshold(self):
        """Latitude just below threshold should work normally."""
        jd = 2451545.0

        # Get the actual threshold for the current obliquity
        # At J2000, obliquity is ~23.44°, so threshold is ~66.56°
        # Test at 66.0° which should be safely below
        cusps, ascmc = ephem.houses(jd, 66.0, 0.0, ord("P"))
        assert len(cusps) == 12

    @pytest.mark.unit
    def test_just_above_threshold(self):
        """Latitude just above threshold should raise error."""
        jd = 2451545.0

        # At 67° it should definitely fail
        with pytest.raises(PolarCircleError):
            ephem.houses(jd, 67.0, 0.0, ord("P"))

    @pytest.mark.unit
    def test_different_times_different_thresholds(self):
        """Different Julian dates have slightly different obliquities."""
        # J2000.0
        jd_2000 = 2451545.0
        # J1900.0
        jd_1900 = 2415020.0

        # Both should raise error at 70°
        with pytest.raises(PolarCircleError) as exc1:
            ephem.houses(jd_2000, 70.0, 0.0, ord("P"))

        with pytest.raises(PolarCircleError) as exc2:
            ephem.houses(jd_1900, 70.0, 0.0, ord("P"))

        # Thresholds should be slightly different due to different obliquities
        # The difference is very small (~0.01° per century)
        assert exc1.value.threshold is not None
        assert exc2.value.threshold is not None


class TestBackwardCompatibility:
    """Test backward compatibility with existing code."""

    @pytest.mark.unit
    def test_still_raises_error_type(self):
        """PolarCircleError should still be catchable as the base Error."""
        jd = 2451545.0

        # Old code catching Error should still work
        try:
            ephem.houses(jd, 70.0, 0.0, ord("P"))
            pytest.fail("Should have raised an error")
        except ephem.Error:
            pass  # Old code would catch this

    @pytest.mark.unit
    def test_error_message_still_mentions_polar_circle(self):
        """Error message should still contain 'polar circle' for compatibility."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.houses(jd, 70.0, 0.0, ord("P"))

        # Tests may check for this string
        assert "polar" in str(exc_info.value).lower()

    @pytest.mark.unit
    def test_houses_armc_also_raises_polar_error(self):
        """houses_armc should also raise PolarCircleError."""
        armc = 280.0
        eps = 23.44

        with pytest.raises(PolarCircleError):
            ephem.houses_armc(armc, 70.0, eps, ord("P"))


class TestKochAndGauquelinPolar:
    """Test polar handling for Koch and Gauquelin systems."""

    @pytest.mark.unit
    def test_koch_raises_polar_error(self):
        """Koch should raise PolarCircleError at polar latitudes."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.houses(jd, 70.0, 0.0, ord("K"))

        assert exc_info.value.house_system == "K"
        assert "Koch" in str(exc_info.value)

    @pytest.mark.unit
    def test_gauquelin_raises_polar_error(self):
        """Gauquelin should raise PolarCircleError at polar latitudes."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.houses(jd, 70.0, 0.0, ord("G"))

        assert exc_info.value.house_system == "G"
        assert "Gauquelin" in str(exc_info.value)

    @pytest.mark.unit
    def test_koch_fallback_works(self):
        """houses_with_fallback should work for Koch."""
        jd = 2451545.0

        cusps, _, used_fallback, warning = ephem.houses_with_fallback(
            jd, 70.0, 0.0, ord("K")
        )

        assert used_fallback is True
        assert "Koch" in warning
        assert len(cusps) == 12

    @pytest.mark.unit
    def test_gauquelin_fallback_works(self):
        """houses_with_fallback should work for Gauquelin."""
        jd = 2451545.0

        cusps, _, used_fallback, warning = ephem.houses_with_fallback(
            jd, 70.0, 0.0, ord("G")
        )

        assert used_fallback is True
        assert "Gauquelin" in warning
        assert len(cusps) == 12


class TestKochPolarCircleErrorDetails:
    """Test that Koch PolarCircleError provides helpful information (similar to Placidus)."""

    @pytest.mark.unit
    def test_koch_error_includes_latitude(self):
        """Koch error message should include the problematic latitude."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.houses(jd, 70.0, 0.0, ord("K"))

        assert exc_info.value.latitude == 70.0
        assert "70" in str(exc_info.value)

    @pytest.mark.unit
    def test_koch_error_includes_threshold(self):
        """Koch error message should include the threshold."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.houses(jd, 70.0, 0.0, ord("K"))

        assert exc_info.value.threshold is not None
        assert exc_info.value.threshold > 66.0
        assert exc_info.value.threshold < 67.0

    @pytest.mark.unit
    def test_koch_error_includes_house_system(self):
        """Koch error message should include the house system."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.houses(jd, 70.0, 0.0, ord("K"))

        assert exc_info.value.house_system == "K"
        assert "Koch" in str(exc_info.value)

    @pytest.mark.unit
    def test_koch_error_suggests_alternatives(self):
        """Koch error message should suggest alternative house systems."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.houses(jd, 70.0, 0.0, ord("K"))

        msg = str(exc_info.value).lower()
        assert "porphyry" in msg or "equal" in msg or "whole sign" in msg

    @pytest.mark.unit
    def test_koch_error_mentions_fallback_function(self):
        """Koch error message should mention houses_with_fallback."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.houses(jd, 70.0, 0.0, ord("K"))

        assert "fallback" in str(exc_info.value).lower()

    @pytest.mark.unit
    @pytest.mark.parametrize("lat", [67.0, 70.0, 80.0, 89.0])
    def test_koch_northern_polar_latitudes(self, lat):
        """Koch at northern polar latitudes should raise detailed error."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.houses(jd, lat, 0.0, ord("K"))

        assert exc_info.value.latitude == lat
        # Should mention Northern
        assert "N" in str(exc_info.value) or "Northern" in str(exc_info.value)

    @pytest.mark.unit
    @pytest.mark.parametrize("lat", [-67.0, -70.0, -80.0, -89.0])
    def test_koch_southern_polar_latitudes(self, lat):
        """Koch at southern polar latitudes should raise detailed error."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.houses(jd, lat, 0.0, ord("K"))

        assert exc_info.value.latitude == lat
        # Should mention Southern
        assert "S" in str(exc_info.value) or "Southern" in str(exc_info.value)


class TestKochPolarEdgeCases:
    """Test Koch edge cases at and near the polar circle threshold."""

    @pytest.mark.unit
    def test_koch_just_below_threshold(self):
        """Koch at latitude just below threshold should work normally."""
        jd = 2451545.0

        # At J2000, obliquity is ~23.44deg, so threshold is ~66.56deg
        # Test at 66.0deg which should be safely below
        cusps, ascmc = ephem.houses(jd, 66.0, 0.0, ord("K"))
        assert len(cusps) == 12

    @pytest.mark.unit
    def test_koch_just_above_threshold(self):
        """Koch at latitude just above threshold should raise error."""
        jd = 2451545.0

        # At 67deg it should definitely fail
        with pytest.raises(PolarCircleError):
            ephem.houses(jd, 67.0, 0.0, ord("K"))

    @pytest.mark.unit
    def test_koch_matches_porphyry_when_fallback_used(self):
        """Koch fallback cusps should match direct Porphyry calculation."""
        jd = 2451545.0
        lat = 70.0

        # Get fallback result
        cusps_fallback, ascmc_fallback, used_fallback, _ = ephem.houses_with_fallback(
            jd, lat, 0.0, ord("K")
        )

        # Get direct Porphyry
        cusps_porphyry, ascmc_porphyry = ephem.houses(jd, lat, 0.0, ord("O"))

        # Should match
        assert used_fallback is True
        for i in range(12):
            assert abs(cusps_fallback[i] - cusps_porphyry[i]) < 0.001

    @pytest.mark.unit
    def test_koch_returns_valid_cusps_at_extreme_latitudes(self):
        """Koch fallback should return valid cusps even at extreme polar latitudes."""
        jd = 2451545.0

        for lat in [80.0, 85.0, 89.0, -80.0, -85.0, -89.0]:
            cusps, ascmc, _, _ = ephem.houses_with_fallback(jd, lat, 0.0, ord("K"))
            for cusp in cusps:
                assert 0 <= cusp < 360

    @pytest.mark.unit
    def test_koch_custom_fallback_system(self):
        """Koch should use custom fallback system when specified."""
        jd = 2451545.0

        cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
            jd,
            70.0,
            0.0,
            ord("K"),
            fallback_hsys=ord("E"),  # Equal houses
        )

        assert used_fallback is True
        assert "Equal" in warning

    @pytest.mark.unit
    def test_koch_warning_includes_threshold(self):
        """Koch warning message should include the threshold."""
        jd = 2451545.0

        _, _, _, warning = ephem.houses_with_fallback(jd, 70.0, 0.0, ord("K"))

        assert "threshold" in warning.lower()


class TestKochBackwardCompatibility:
    """Test Koch backward compatibility with existing code."""

    @pytest.mark.unit
    def test_koch_still_raises_error_type(self):
        """Koch PolarCircleError should still be catchable as the base Error."""
        jd = 2451545.0

        # Old code catching Error should still work
        try:
            ephem.houses(jd, 70.0, 0.0, ord("K"))
            pytest.fail("Should have raised an error")
        except ephem.Error:
            pass  # Old code would catch this

    @pytest.mark.unit
    def test_koch_error_message_still_mentions_polar_circle(self):
        """Koch error message should still contain 'polar circle' for compatibility."""
        jd = 2451545.0

        with pytest.raises(PolarCircleError) as exc_info:
            ephem.houses(jd, 70.0, 0.0, ord("K"))

        # Tests may check for this string
        assert "polar" in str(exc_info.value).lower()

    @pytest.mark.unit
    def test_koch_houses_armc_also_raises_polar_error(self):
        """Koch houses_armc should also raise PolarCircleError."""
        armc = 280.0
        eps = 23.44

        with pytest.raises(PolarCircleError):
            ephem.houses_armc(armc, 70.0, eps, ord("K"))

    @pytest.mark.unit
    def test_koch_armc_fallback_works(self):
        """houses_armc_with_fallback should work for Koch."""
        armc = 280.0
        eps = 23.44

        cusps, ascmc, used_fallback, warning = ephem.houses_armc_with_fallback(
            armc, 70.0, eps, ord("K")
        )

        assert used_fallback is True
        assert "Koch" in warning
        assert len(cusps) == 12

    @pytest.mark.unit
    def test_koch_armc_custom_fallback(self):
        """Koch houses_armc_with_fallback should use custom fallback when specified."""
        armc = 280.0
        eps = 23.44

        _, _, used_fallback, warning = ephem.houses_armc_with_fallback(
            armc, 70.0, eps, ord("K"), fallback_hsys=ord("W")
        )

        assert used_fallback is True
        assert "Whole Sign" in warning


class TestExtremeLatitudeInfo:
    """Test the get_extreme_latitude_info() function."""

    @pytest.mark.unit
    def test_extreme_latitude_info_at_normal_latitude(self):
        """Should return correct info for normal latitudes."""
        info = ephem.get_extreme_latitude_info(45.0)

        assert info["latitude"] == 45.0
        assert info["is_extreme"] is False
        assert info["is_polar_circle"] is False
        assert info["hemisphere"] == "N"
        assert info["affected_systems"] == []
        assert info["unstable_systems"] == []
        assert "E" in info["stable_systems"]
        assert "W" in info["stable_systems"]
        assert "O" in info["stable_systems"]

    @pytest.mark.unit
    def test_extreme_latitude_info_at_polar_circle(self):
        """Should correctly identify polar circle locations."""
        info = ephem.get_extreme_latitude_info(70.0)

        assert info["latitude"] == 70.0
        assert info["is_extreme"] is False  # 70 < 80
        assert info["is_polar_circle"] is True  # 70 > ~66.56
        assert "P" in info["affected_systems"]
        assert "K" in info["affected_systems"]
        assert "G" in info["affected_systems"]

    @pytest.mark.unit
    def test_extreme_latitude_info_at_extreme_latitude(self):
        """Should correctly identify extreme latitude locations."""
        info = ephem.get_extreme_latitude_info(85.0)

        assert info["latitude"] == 85.0
        assert info["is_extreme"] is True
        assert info["is_polar_circle"] is True
        assert "C" in info["unstable_systems"]  # Campanus
        assert "R" in info["unstable_systems"]  # Regiomontanus
        assert "T" in info["unstable_systems"]  # Topocentric
        assert info["extreme_threshold"] == 80.0

    @pytest.mark.unit
    def test_extreme_latitude_info_southern_hemisphere(self):
        """Should correctly identify southern hemisphere."""
        info = ephem.get_extreme_latitude_info(-85.0)

        assert info["hemisphere"] == "S"
        assert info["is_extreme"] is True

    @pytest.mark.unit
    def test_extreme_latitude_threshold_constant(self):
        """EXTREME_LATITUDE_THRESHOLD should be 80."""
        assert ephem.EXTREME_LATITUDE_THRESHOLD == 80.0


class TestExtremeLatitudeWarnings:
    """Test that extreme latitude warnings are generated correctly."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys,name",
        [
            (ord("C"), "Campanus"),
            (ord("R"), "Regiomontanus"),
            (ord("T"), "Topocentric"),
        ],
    )
    def test_warning_for_unstable_systems_at_extreme_latitude(self, hsys, name):
        """Should generate warning for potentially unstable systems at >80°."""
        jd = 2451545.0

        cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
            jd, 85.0, 0.0, hsys
        )

        # Should NOT use fallback - these systems still work
        assert used_fallback is False
        # But should have a warning
        assert warning is not None
        assert "reduced accuracy" in warning.lower() or "may have" in warning.lower()
        assert name in warning

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys",
        [
            ord("E"),  # Equal
            ord("W"),  # Whole Sign
            ord("O"),  # Porphyry
            ord("M"),  # Morinus
        ],
    )
    def test_no_warning_for_stable_systems_at_extreme_latitude(self, hsys):
        """Stable systems should have no warning at extreme latitudes."""
        jd = 2451545.0

        cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
            jd, 85.0, 0.0, hsys
        )

        assert used_fallback is False
        assert warning is None

    @pytest.mark.unit
    def test_polar_failing_systems_use_fallback_at_extreme_latitude(self):
        """Placidus/Koch/Gauquelin should use fallback at extreme latitudes."""
        jd = 2451545.0

        for hsys in [ord("P"), ord("K"), ord("G")]:
            cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
                jd, 85.0, 0.0, hsys
            )

            assert used_fallback is True
            assert warning is not None
            assert len(cusps) == 12


class TestCuspValidation:
    """Test cusp validation at extreme latitudes."""

    @pytest.mark.unit
    def test_valid_cusps_at_normal_latitude(self):
        """Cusps should be valid at normal latitudes."""
        jd = 2451545.0

        cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
            jd, 45.0, 0.0, ord("C"), validate_cusps=True
        )

        for cusp in cusps:
            assert 0 <= cusp < 360
            assert not math.isnan(cusp)
            assert not math.isinf(cusp)

    @pytest.mark.unit
    def test_valid_cusps_at_extreme_latitude(self):
        """Cusps should still be valid at extreme latitudes (possibly with fallback)."""
        jd = 2451545.0

        for lat in [80.0, 85.0, 89.0, -85.0]:
            cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
                jd, lat, 0.0, ord("C"), validate_cusps=True
            )

            # All cusps should be valid regardless of fallback
            for cusp in cusps:
                assert 0 <= cusp < 360
                assert not math.isnan(cusp)
                assert not math.isinf(cusp)

    @pytest.mark.unit
    def test_validation_can_be_disabled(self):
        """Cusp validation should be disableable."""
        jd = 2451545.0

        # Just verify it doesn't crash with validation disabled
        cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
            jd, 85.0, 0.0, ord("C"), validate_cusps=False
        )

        assert len(cusps) == 12


class TestAllHouseSystemsAtExtremeLatitudes:
    """Test all house systems at various extreme latitudes."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "lat", [80.0, 85.0, 89.0, 89.5, -80.0, -85.0, -89.0, -89.5]
    )
    def test_all_systems_produce_valid_output_with_fallback(self, lat):
        """All house systems should produce valid output with fallback enabled."""
        jd = 2451545.0

        # All house system codes
        all_systems = [
            ord("P"),
            ord("K"),
            ord("R"),
            ord("C"),
            ord("E"),
            ord("W"),
            ord("O"),
            ord("B"),
            ord("T"),
            ord("M"),
            ord("X"),
            ord("V"),
            ord("H"),
            ord("G"),
            ord("U"),
            ord("F"),
            ord("Y"),
            ord("N"),
        ]

        for hsys in all_systems:
            cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
                jd, lat, 0.0, hsys
            )

            # Should always get 12 valid cusps
            assert len(cusps) == 12, f"System {chr(hsys)} at lat {lat}"
            for i, cusp in enumerate(cusps):
                assert 0 <= cusp < 360, (
                    f"System {chr(hsys)} at lat {lat}: cusp {i + 1}={cusp}"
                )
                assert not math.isnan(cusp), (
                    f"System {chr(hsys)} at lat {lat}: cusp {i + 1} is NaN"
                )

    @pytest.mark.unit
    def test_89_5_degrees_latitude(self):
        """Test at 89.5° (near pole) with fallback."""
        jd = 2451545.0
        lat = 89.5

        # Placidus should fall back
        cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
            jd, lat, 0.0, ord("P")
        )

        assert used_fallback is True
        assert len(cusps) == 12
        for cusp in cusps:
            assert 0 <= cusp < 360

    @pytest.mark.unit
    def test_multiple_julian_days_at_extreme_latitude(self):
        """Test multiple Julian days at extreme latitudes."""
        lat = 85.0

        julian_days = [
            2451545.0,  # J2000.0
            2455000.0,  # 2009
            2460000.0,  # 2023
        ]

        for jd in julian_days:
            cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
                jd,
                lat,
                0.0,
                ord("O"),  # Porphyry should always work
            )

            assert used_fallback is False
            assert warning is None
            assert len(cusps) == 12


class TestArmcWithFallbackAtExtremeLatitudes:
    """Test houses_armc_with_fallback at extreme latitudes."""

    @pytest.mark.unit
    def test_armc_fallback_at_extreme_latitude(self):
        """Should provide fallback for ARMC-based calculations at extreme latitudes."""
        armc = 280.0
        eps = 23.44
        lat = 85.0

        # Placidus should fail
        cusps, ascmc, used_fallback, warning = ephem.houses_armc_with_fallback(
            armc, lat, eps, ord("P")
        )

        assert used_fallback is True
        assert len(cusps) == 12

    @pytest.mark.unit
    def test_armc_warning_for_unstable_systems(self):
        """ARMC-based calculations should warn about unstable systems."""
        armc = 280.0
        eps = 23.44
        lat = 85.0

        # Campanus should warn but not fail
        cusps, ascmc, used_fallback, warning = ephem.houses_armc_with_fallback(
            armc, lat, eps, ord("C")
        )

        assert used_fallback is False
        assert warning is not None
        assert "Campanus" in warning

    @pytest.mark.unit
    def test_armc_cusp_validation(self):
        """ARMC-based calculations should validate cusps."""
        armc = 280.0
        eps = 23.44
        lat = 85.0

        cusps, ascmc, used_fallback, warning = ephem.houses_armc_with_fallback(
            armc, lat, eps, ord("R"), validate_cusps=True
        )

        for cusp in cusps:
            assert 0 <= cusp < 360


class TestEdgeCasesNearPoles:
    """Test edge cases very near the geographic poles."""

    @pytest.mark.unit
    def test_latitude_89_9(self):
        """Test at 89.9° latitude."""
        jd = 2451545.0

        cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
            jd, 89.9, 0.0, ord("P")
        )

        assert used_fallback is True  # Placidus fails at polar circle
        for cusp in cusps:
            assert 0 <= cusp < 360

    @pytest.mark.unit
    def test_latitude_exactly_90_raises_error(self):
        """Latitude of exactly ±90 should raise CoordinateError."""
        jd = 2451545.0

        # 90.1 degrees should be out of range
        with pytest.raises(ephem.CoordinateError):
            ephem.houses(jd, 90.1, 0.0, ord("O"))

    @pytest.mark.unit
    def test_negative_extreme_latitude(self):
        """Should handle extreme southern latitudes."""
        jd = 2451545.0

        cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
            jd, -89.0, 0.0, ord("O")
        )

        assert len(cusps) == 12
        for cusp in cusps:
            assert 0 <= cusp < 360

    @pytest.mark.unit
    def test_varying_longitudes_at_extreme_latitude(self):
        """Test multiple longitudes at extreme latitude."""
        jd = 2451545.0
        lat = 85.0

        # Use longitudes within valid range (-180 to 180)
        for lon in [0.0, 90.0, 180.0, -90.0, -45.0, -120.0]:
            cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
                jd, lat, lon, ord("O")
            )

            assert len(cusps) == 12
            for cusp in cusps:
                assert 0 <= cusp < 360
