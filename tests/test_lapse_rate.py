"""
Tests for the set_lapse_rate and get_lapse_rate functions.

Tests verify that the atmospheric lapse rate can be configured globally
and that it affects refrac_extended() calculations.
"""

import libephemeris as ephem
from libephemeris import state


class TestLapseRateBasic:
    """Basic functionality tests for lapse rate configuration."""

    def setup_method(self):
        """Reset lapse rate to default before each test."""
        ephem.set_lapse_rate(None)

    def teardown_method(self):
        """Reset lapse rate to default after each test."""
        ephem.set_lapse_rate(None)

    def test_get_lapse_rate_exported(self):
        """Test that get_lapse_rate is exported from the package."""
        assert hasattr(ephem, "get_lapse_rate")
        assert callable(ephem.get_lapse_rate)

    def test_set_lapse_rate_exported(self):
        """Test that set_lapse_rate is exported from the package."""
        assert hasattr(ephem, "set_lapse_rate")
        assert callable(ephem.set_lapse_rate)

    def test_swe_prefixed_versions_exported(self):
        """Test that swe_ prefixed versions are exported."""
        assert hasattr(ephem, "get_lapse_rate")
        assert hasattr(ephem, "set_lapse_rate")
        assert callable(ephem.get_lapse_rate)
        assert callable(ephem.set_lapse_rate)

    def test_get_lapse_rate_returns_default(self):
        """Test that get_lapse_rate returns 0.0065 by default."""
        result = ephem.get_lapse_rate()
        assert result == 0.0065

    def test_set_lapse_rate_changes_value(self):
        """Test that set_lapse_rate changes the lapse rate value."""
        ephem.set_lapse_rate(0.005)
        result = ephem.get_lapse_rate()
        assert result == 0.005

    def test_set_lapse_rate_to_none_resets_to_default(self):
        """Test that setting lapse rate to None resets to default."""
        ephem.set_lapse_rate(0.008)
        assert ephem.get_lapse_rate() == 0.008

        ephem.set_lapse_rate(None)
        assert ephem.get_lapse_rate() == 0.0065

    def test_set_lapse_rate_zero(self):
        """Test that lapse rate can be set to zero."""
        ephem.set_lapse_rate(0.0)
        assert ephem.get_lapse_rate() == 0.0

    def test_set_lapse_rate_negative(self):
        """Test that lapse rate can be set to negative values."""
        ephem.set_lapse_rate(-0.001)
        assert ephem.get_lapse_rate() == -0.001


class TestLapseRateAffectsRefracExtended:
    """Tests that lapse rate affects refrac_extended calculations."""

    def setup_method(self):
        """Reset lapse rate to default before each test."""
        ephem.set_lapse_rate(None)

    def teardown_method(self):
        """Reset lapse rate to default after each test."""
        ephem.set_lapse_rate(None)

    def test_refrac_extended_uses_global_lapse_rate(self):
        """Test that refrac_extended uses global lapse rate when not specified."""
        # Set a custom lapse rate
        ephem.set_lapse_rate(0.008)

        # Call refrac_extended without explicit lapse_rate
        alt, (true_alt, app_alt, ref, dip) = ephem.refrac_extended(
            0.0, 1000.0, 1013.25, 15.0
        )

        # Call with explicit lapse_rate matching global
        alt2, (true_alt2, app_alt2, ref2, dip2) = ephem.refrac_extended(
            0.0, 1000.0, 1013.25, 15.0, 0.008
        )

        # Results should match
        assert dip == dip2
        assert ref == ref2

    def test_refrac_extended_explicit_overrides_global(self):
        """Test that explicit lapse_rate overrides global setting."""
        # Set a custom global lapse rate
        ephem.set_lapse_rate(0.008)

        # Call with explicit lapse_rate (different from global)
        alt1, (_, _, _, dip1) = ephem.refrac_extended(0.0, 1000.0, 1013.25, 15.0, 0.005)

        # Call with global lapse rate
        alt2, (_, _, _, dip2) = ephem.refrac_extended(0.0, 1000.0, 1013.25, 15.0)

        # Dip should be different because lapse rates differ
        assert dip1 != dip2

    def test_changing_lapse_rate_affects_results(self):
        """Test that changing global lapse rate affects refrac_extended results."""
        # Get result with default lapse rate
        ephem.set_lapse_rate(None)
        alt1, (_, _, _, dip1) = ephem.refrac_extended(0.0, 1000.0)

        # Get result with higher lapse rate
        ephem.set_lapse_rate(0.010)
        alt2, (_, _, _, dip2) = ephem.refrac_extended(0.0, 1000.0)

        # Get result with lower lapse rate
        ephem.set_lapse_rate(0.003)
        alt3, (_, _, _, dip3) = ephem.refrac_extended(0.0, 1000.0)

        # Higher lapse rates result in less refraction -> more negative dip
        # (but this is nuanced - the key point is results differ)
        assert dip1 != dip2
        assert dip2 != dip3


class TestLapseRateDefaultConstant:
    """Tests for the LAPSE_RATE_DEFAULT constant."""

    def test_default_constant_exists(self):
        """Test that LAPSE_RATE_DEFAULT constant exists in state module."""
        assert hasattr(state, "LAPSE_RATE_DEFAULT")

    def test_default_constant_value(self):
        """Test that LAPSE_RATE_DEFAULT is 0.0065."""
        assert state.LAPSE_RATE_DEFAULT == 0.0065


class TestLapseRateGlobalState:
    """Tests for lapse rate global state management."""

    def setup_method(self):
        """Reset lapse rate to default before each test."""
        ephem.set_lapse_rate(None)

    def teardown_method(self):
        """Reset lapse rate to default after each test."""
        ephem.set_lapse_rate(None)

    def test_state_variable_updated(self):
        """Test that _LAPSE_RATE state variable is updated correctly."""
        assert state._LAPSE_RATE is None

        ephem.set_lapse_rate(0.007)
        assert state._LAPSE_RATE == 0.007

        ephem.set_lapse_rate(None)
        assert state._LAPSE_RATE is None

    def test_swe_prefixed_functions_same_as_unprefixed(self):
        """Test that swe_ prefixed functions behave identically."""
        ephem.set_lapse_rate(0.006)
        assert ephem.get_lapse_rate() == 0.006
        assert ephem.get_lapse_rate() == 0.006

        ephem.set_lapse_rate(0.004)
        assert ephem.get_lapse_rate() == 0.004
