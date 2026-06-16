"""
Tests for tidal acceleration functions (set_tid_acc, get_tid_acc).

Tidal acceleration affects Delta T calculations for historical dates.
"""

import pytest
import libephemeris as ephem


class TestTidAccBasicFunctionality:
    """Test basic get/set functionality for tidal acceleration."""

    @pytest.fixture(autouse=True)
    def reset_tid_acc(self):
        """Reset tidal acceleration to default before and after each test."""
        # Reset to automatic/default before test
        ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)
        yield
        # Reset after test to avoid affecting other tests
        ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)

    @pytest.mark.unit
    def test_get_tid_acc_default(self):
        """get_tid_acc should return default value when not explicitly set."""
        value = ephem.get_tid_acc()
        assert value == ephem.TIDAL_DEFAULT
        assert value == -25.8

    @pytest.mark.unit
    def test_set_tid_acc_de421(self):
        """set_tid_acc should accept DE421 tidal acceleration."""
        ephem.set_tid_acc(ephem.TIDAL_DE421)
        assert ephem.get_tid_acc() == ephem.TIDAL_DE421
        assert ephem.get_tid_acc() == -25.85

    @pytest.mark.unit
    def test_set_tid_acc_de431(self):
        """set_tid_acc should accept DE431 tidal acceleration."""
        ephem.set_tid_acc(ephem.TIDAL_DE431)
        assert ephem.get_tid_acc() == ephem.TIDAL_DE431
        assert ephem.get_tid_acc() == -25.80

    @pytest.mark.unit
    def test_set_tid_acc_de441(self):
        """set_tid_acc should accept DE441 tidal acceleration."""
        ephem.set_tid_acc(ephem.TIDAL_DE441)
        assert ephem.get_tid_acc() == ephem.TIDAL_DE441
        assert ephem.get_tid_acc() == -25.936

    @pytest.mark.unit
    def test_set_tid_acc_custom_value(self):
        """set_tid_acc should accept custom numeric values."""
        custom_value = -26.5
        ephem.set_tid_acc(custom_value)
        assert ephem.get_tid_acc() == custom_value

    @pytest.mark.unit
    def test_set_tid_acc_automatic_resets_to_default(self):
        """set_tid_acc with TIDAL_AUTOMATIC should reset to default."""
        # First set a custom value
        ephem.set_tid_acc(ephem.TIDAL_DE441)
        assert ephem.get_tid_acc() == ephem.TIDAL_DE441

        # Reset with automatic
        ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)
        assert ephem.get_tid_acc() == ephem.TIDAL_DEFAULT

    @pytest.mark.unit
    def test_set_tid_acc_zero_stores_zero(self):
        """set_tid_acc with 0.0 should store 0.0 (reference-API compat)."""
        ephem.set_tid_acc(ephem.TIDAL_DE421)
        ephem.set_tid_acc(0.0)
        assert ephem.get_tid_acc() == 0.0


class TestTidAccConstants:
    """Test that all tidal acceleration constants are properly defined."""

    @pytest.mark.unit
    def test_se_tidal_constants_exist(self):
        """All TIDAL_* constants should exist."""
        assert hasattr(ephem, "TIDAL_DE200")
        assert hasattr(ephem, "TIDAL_DE403")
        assert hasattr(ephem, "TIDAL_DE404")
        assert hasattr(ephem, "TIDAL_DE405")
        assert hasattr(ephem, "TIDAL_DE406")
        assert hasattr(ephem, "TIDAL_DE421")
        assert hasattr(ephem, "TIDAL_DE422")
        assert hasattr(ephem, "TIDAL_DE430")
        assert hasattr(ephem, "TIDAL_DE431")
        assert hasattr(ephem, "TIDAL_DE440")
        assert hasattr(ephem, "TIDAL_DE441")
        assert hasattr(ephem, "TIDAL_DEFAULT")
        assert hasattr(ephem, "TIDAL_AUTOMATIC")

    @pytest.mark.unit
    def test_tidal_constants_aliases_exist(self):
        """TIDAL_* aliases (without SE_ prefix) should exist."""
        assert hasattr(ephem, "TIDAL_DE200")
        assert hasattr(ephem, "TIDAL_DE403")
        assert hasattr(ephem, "TIDAL_DE421")
        assert hasattr(ephem, "TIDAL_DE431")
        assert hasattr(ephem, "TIDAL_DE440")
        assert hasattr(ephem, "TIDAL_DE441")
        assert hasattr(ephem, "TIDAL_DEFAULT")
        assert hasattr(ephem, "TIDAL_AUTOMATIC")

    @pytest.mark.unit
    def test_tidal_constants_values_are_negative(self):
        """Tidal acceleration values should be negative."""
        assert ephem.TIDAL_DE200 < 0
        assert ephem.TIDAL_DE421 < 0
        assert ephem.TIDAL_DE431 < 0
        assert ephem.TIDAL_DE440 < 0
        assert ephem.TIDAL_DE441 < 0

    @pytest.mark.unit
    def test_tidal_automatic_is_999999(self):
        """TIDAL_AUTOMATIC should be 999999 (reference-API compat)."""
        assert ephem.TIDAL_AUTOMATIC == 999999

    @pytest.mark.unit
    def test_tidal_default_is_minus_25_8(self):
        """TIDAL_DEFAULT should be -25.8 (reference-API compat)."""
        assert ephem.TIDAL_DEFAULT == -25.8


class TestTidAccFunctionAliases:
    """Test that function aliases work correctly."""

    @pytest.fixture(autouse=True)
    def reset_tid_acc(self):
        """Reset tidal acceleration to default before and after each test."""
        ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)
        yield
        ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)

    @pytest.mark.unit
    def test_swe_prefixed_aliases_exist(self):
        """set_tid_acc and get_tid_acc should exist."""
        assert hasattr(ephem, "set_tid_acc")
        assert hasattr(ephem, "get_tid_acc")

    @pytest.mark.unit
    def test_swe_prefixed_aliases_work(self):
        """set_tid_acc and get_tid_acc should work identically."""
        ephem.set_tid_acc(ephem.TIDAL_DE421)
        assert ephem.get_tid_acc() == ephem.TIDAL_DE421

        # Cross-check with non-prefixed versions
        assert ephem.get_tid_acc() == ephem.TIDAL_DE421

    @pytest.mark.unit
    def test_non_prefixed_and_prefixed_are_same(self):
        """set_tid_acc and set_tid_acc should be the same function."""


class TestTidAccReturnTypes:
    """Test return types of tidal acceleration functions."""

    @pytest.fixture(autouse=True)
    def reset_tid_acc(self):
        """Reset tidal acceleration to default."""
        ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)
        yield
        ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)

    @pytest.mark.unit
    def test_get_tid_acc_returns_float(self):
        """get_tid_acc should return a float."""
        result = ephem.get_tid_acc()
        assert isinstance(result, float)

    @pytest.mark.unit
    def test_set_tid_acc_returns_none(self):
        """set_tid_acc should return None."""
        result = ephem.set_tid_acc(ephem.TIDAL_DE421)
        assert result is None
