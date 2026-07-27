"""
Tests for type safety and mypy compliance.

These tests verify that the type annotations in libephemeris are correct
and that functions return values with the expected types.
"""

from pathlib import Path

import libephemeris as ephem
from libephemeris import EphemerisContext
from libephemeris.state import get_sid_mode, set_sid_mode, close
from libephemeris.houses import houses, houses_ex, house_name
from libephemeris.constants import SIDM_LAHIRI


class TestPEP561Compliance:
    """Tests for PEP 561 py.typed marker file."""

    def test_py_typed_exists_in_package(self):
        """py.typed marker file must exist in the package directory for PEP 561 compliance."""
        package_dir = Path(ephem.__file__).parent
        py_typed_path = package_dir / "py.typed"
        assert py_typed_path.exists(), (
            f"py.typed marker file not found at {py_typed_path}. "
            "This file is required for PEP 561 compliance so type checkers "
            "can use libephemeris type hints."
        )

    def test_py_typed_is_empty_or_minimal(self):
        """py.typed should be an empty marker file."""
        package_dir = Path(ephem.__file__).parent
        py_typed_path = package_dir / "py.typed"
        content = py_typed_path.read_text()
        # py.typed should be empty or contain only whitespace/comments
        assert content.strip() == "" or content.startswith("#"), (
            "py.typed should be empty or contain only comments"
        )


class TestGetSidMode:
    """Tests for get_sid_mode type safety."""

    def setup_method(self):
        """Reset state before each test."""
        close()

    def teardown_method(self):
        """Reset state after each test."""
        close()

    def test_get_sid_mode_returns_int_when_not_set(self):
        """get_sid_mode() should return int even when mode is not set."""
        # After close(), _SIDEREAL_MODE is None; the default mirrors the
        # reference API: SIDM_FAGAN_BRADLEY (0)
        close()
        mode = get_sid_mode()
        assert isinstance(mode, int)
        assert mode == 0  # Default to SIDM_FAGAN_BRADLEY

    def test_get_sid_mode_returns_int_when_set(self):
        """get_sid_mode() should return int when mode is explicitly set."""
        set_sid_mode(SIDM_LAHIRI)
        mode = get_sid_mode()
        assert isinstance(mode, int)
        assert mode == SIDM_LAHIRI

    def test_get_sid_mode_full_returns_tuple(self):
        """get_sid_mode(full=True) should return a tuple of (int, float, float)."""
        set_sid_mode(SIDM_LAHIRI, 2451545.0, 23.5)
        result = get_sid_mode(full=True)
        assert isinstance(result, tuple)
        assert len(result) == 3
        mode, t0, ayan_t0 = result
        assert isinstance(mode, int)
        assert isinstance(t0, float)
        assert isinstance(ayan_t0, float)

    def test_get_sid_mode_full_returns_default_int_when_not_set(self):
        """get_sid_mode(full=True) returns int (not None) for mode even when unset."""
        close()
        result = get_sid_mode(full=True)
        assert isinstance(result, tuple)
        mode, t0, ayan_t0 = result
        # Mode should be an int (default 0 = Fagan/Bradley), not None
        assert isinstance(mode, int)
        assert mode == 0


class TestHousesReturnTypes:
    """Tests for house calculation return types."""

    def test_swe_houses_returns_tuples(self):
        """houses should return tuple of tuples, not lists."""
        cusps, ascmc = houses(2451545.0, 51.5, -0.12, ord("P"))
        assert isinstance(cusps, tuple)
        assert isinstance(ascmc, tuple)
        assert len(cusps) == 12
        assert len(ascmc) == 8

    def test_swe_houses_ex_returns_tuples(self):
        """houses_ex should return tuple of tuples."""
        cusps, ascmc = houses_ex(2451545.0, 51.5, -0.12, ord("P"), 0)
        assert isinstance(cusps, tuple)
        assert isinstance(ascmc, tuple)

    def test_swe_house_name_with_int(self):
        """house_name should accept int and return str."""
        name = house_name(ord("P"))
        assert isinstance(name, str)
        assert name == "Placidus"

    def test_swe_house_name_with_unknown(self):
        """house_name returns '' for an unknown selector (measured contract)."""
        name = house_name(ord("Z"))  # Unknown system
        assert isinstance(name, str)
        assert name == ""


class TestEphemerisContextHouses:
    """Tests for EphemerisContext.houses return types."""

    def test_context_houses_returns_tuples(self):
        """EphemerisContext.houses should return tuple of tuples."""
        ctx = EphemerisContext()
        cusps, ascmc = ctx.houses(2451545.0, 51.5, -0.12, ord("P"))
        assert isinstance(cusps, tuple)
        assert isinstance(ascmc, tuple)
        assert len(cusps) == 12
        assert len(ascmc) == 8


class TestAyanamsaTypes:
    """Tests for ayanamsa calculation types."""

    def test_get_ayanamsa_ut_returns_float(self):
        """get_ayanamsa_ut should return float."""
        set_sid_mode(SIDM_LAHIRI)
        ayanamsa = ephem.get_ayanamsa_ut(2451545.0)
        assert isinstance(ayanamsa, float)

    def test_get_ayanamsa_ut_with_unset_mode(self):
        """get_ayanamsa_ut works even when mode was never explicitly set."""
        close()
        # Should use default mode (1 = Lahiri)
        ayanamsa = ephem.get_ayanamsa_ut(2451545.0)
        assert isinstance(ayanamsa, float)
