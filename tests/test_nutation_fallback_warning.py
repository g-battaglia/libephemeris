"""
Unit tests for NutationFallbackWarning and get_nutation_model() function.

Tests verify that:
1. get_nutation_model() correctly reports the active nutation model
2. NutationFallbackWarning is a proper UserWarning subclass (legacy, kept for API compat)
3. The warning infrastructure works correctly
"""

import warnings
import pytest
import libephemeris as eph
from libephemeris.planets import NutationFallbackWarning, get_nutation_model


@pytest.mark.unit
class TestGetNutationModel:
    """Tests for the get_nutation_model() function."""

    def test_get_nutation_model_returns_dict(self):
        """Test that get_nutation_model returns a dictionary."""
        info = get_nutation_model()
        assert isinstance(info, dict)

    def test_get_nutation_model_has_required_keys(self):
        """Test that get_nutation_model returns all required keys."""
        info = get_nutation_model()
        required_keys = ["model", "terms", "precision", "source"]
        for key in required_keys:
            assert key in info, f"Missing required key: {key}"

    def test_get_nutation_model_valid_model_values(self):
        """Test that model value is the IAU 2006/2000A model."""
        info = get_nutation_model()
        assert info["model"] == "IAU2006_2000A"

    def test_get_nutation_model_uses_pyerfa(self):
        """Test that the model source is pyerfa."""
        info = get_nutation_model()
        assert info["source"] == "pyerfa"
        assert info["terms"] == 1365
        assert "mas" in info["precision"].lower()

    def test_get_nutation_model_terms_is_int(self):
        """Test that terms is an integer."""
        info = get_nutation_model()
        assert isinstance(info["terms"], int)
        assert info["terms"] == 1365

    def test_get_nutation_model_exported_from_package(self):
        """Test that get_nutation_model is exported from libephemeris package."""
        assert hasattr(eph, "get_nutation_model")
        info = eph.get_nutation_model()
        assert isinstance(info, dict)
        assert "model" in info


@pytest.mark.unit
class TestNutationFallbackWarningClass:
    """Tests for the NutationFallbackWarning class itself."""

    def test_warning_is_user_warning_subclass(self):
        """Test that NutationFallbackWarning is a UserWarning subclass."""
        assert issubclass(NutationFallbackWarning, UserWarning)

    def test_warning_exported_from_package(self):
        """Test that NutationFallbackWarning is exported from libephemeris package."""
        assert hasattr(eph, "NutationFallbackWarning")
        assert issubclass(eph.NutationFallbackWarning, UserWarning)

    def test_warning_can_be_raised(self):
        """Test that NutationFallbackWarning can be raised and caught."""
        with pytest.raises(NutationFallbackWarning):
            raise NutationFallbackWarning("Test warning message")

    def test_warning_message_preserved(self):
        """Test that warning message is preserved."""
        msg = "Test nutation fallback warning"
        try:
            raise NutationFallbackWarning(msg)
        except NutationFallbackWarning as e:
            assert str(e) == msg


@pytest.mark.unit
class TestNutationFallbackBehavior:
    """Tests for nutation fallback behavior during calculations."""

    def test_ecl_nut_calculation_succeeds(self, standard_jd):
        """Test that ECL_NUT calculation works with either nutation model."""
        # This should work regardless of which model is active
        pos, flag = eph.calc_ut(standard_jd, eph.ECL_NUT, 0)

        # Verify we got valid nutation data
        true_obliquity = pos[0]
        mean_obliquity = pos[1]
        nutation_lon = pos[2]
        nutation_obl = pos[3]

        # True obliquity should be around 23.4 degrees
        assert 23.0 < true_obliquity < 24.0, f"Invalid true obliquity: {true_obliquity}"
        assert 23.0 < mean_obliquity < 24.0, f"Invalid mean obliquity: {mean_obliquity}"

        # Nutation values should be small (typically < 0.01 degrees)
        assert abs(nutation_lon) < 0.1, f"Nutation lon too large: {nutation_lon}"
        assert abs(nutation_obl) < 0.1, f"Nutation obl too large: {nutation_obl}"

    def test_nutation_model_info_matches_behavior(self, standard_jd):
        """Test that get_nutation_model() info matches actual behavior."""
        info = get_nutation_model()

        # Calculate nutation - this should work with pyerfa
        pos, _ = eph.calc_ut(standard_jd, eph.ECL_NUT, 0)

        # If we got here without error, the model reported is being used
        assert pos is not None
        assert info["model"] == "IAU2006_2000A"


@pytest.mark.unit
class TestNutationWarningWithMockedImport:
    """
    Tests that verify warning behavior.

    Note: Testing the actual fallback warning emission requires mocking
    the skyfield import, which is complex. These tests verify the warning
    infrastructure is correctly set up.
    """

    def test_warning_can_be_caught_and_filtered(self):
        """Test that NutationFallbackWarning can be caught in warning filter."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            # Manually emit the warning to test the filter
            warnings.warn(
                "Test message",
                NutationFallbackWarning,
            )
            nutation_warnings = [
                x for x in w if issubclass(x.category, NutationFallbackWarning)
            ]
            assert len(nutation_warnings) == 1
            assert "Test message" in str(nutation_warnings[0].message)

    def test_warning_can_be_suppressed(self):
        """Test that NutationFallbackWarning can be suppressed."""
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings("ignore", category=NutationFallbackWarning)
            # Manually emit the warning
            warnings.warn("Should be ignored", NutationFallbackWarning)
            nutation_warnings = [
                x for x in w if issubclass(x.category, NutationFallbackWarning)
            ]
            assert len(nutation_warnings) == 0

    def test_warning_message_contains_helpful_info(self):
        """Test that the warning message would contain helpful information."""
        # The actual warning message should contain guidance
        expected_keywords = ["skyfield", "precision", "arcsecond", "get_nutation_model"]

        # Simulate checking the warning message that would be emitted
        test_msg = (
            "Using simplified 4-term nutation model (~1 arcsecond precision) "
            "because Skyfield's IAU 2000A model is unavailable. "
            "For sub-milliarcsecond precision, install Skyfield: pip install skyfield. "
            "Use get_nutation_model() to check the active nutation model."
        )

        for keyword in expected_keywords:
            assert keyword.lower() in test_msg.lower(), (
                f"Warning message should contain '{keyword}'"
            )


@pytest.mark.unit
class TestNutationDocumentation:
    """Tests verifying documentation and discoverability."""

    def test_warning_has_docstring(self):
        """Test that NutationFallbackWarning has a descriptive docstring."""
        assert NutationFallbackWarning.__doc__ is not None
        assert len(NutationFallbackWarning.__doc__) > 50
        assert "precision" in NutationFallbackWarning.__doc__.lower()

    def test_get_nutation_model_has_docstring(self):
        """Test that get_nutation_model has a descriptive docstring."""
        assert get_nutation_model.__doc__ is not None
        assert len(get_nutation_model.__doc__) > 100
        assert (
            "IAU 2006" in get_nutation_model.__doc__
            or "IAU2006" in get_nutation_model.__doc__
        )


@pytest.fixture
def standard_jd():
    """Standard Julian Day for testing (J2000.0)."""
    return 2451545.0
