"""
Unit tests for HIP prefix search functionality in star lookup.

Tests the ability to search for stars using "HIP NNNNN" format,
which is required for pyswisseph fixstar2 compatibility.

Examples:
    - "HIP 49669" -> Regulus
    - "HIP 65474" -> Spica
    - "HIP65474" -> Spica (no space)
"""

import pytest
from libephemeris.fixed_stars import _resolve_star2, fixstar2_ut


@pytest.mark.unit
class TestHipPrefixSearch:
    """Tests for HIP prefix search in _resolve_star2."""

    def test_hip_prefix_with_space_regulus(self):
        """Test HIP 49669 resolves to Regulus."""
        entry, err = _resolve_star2("HIP 49669")
        assert err is None, f"Unexpected error: {err}"
        assert entry is not None
        assert entry.name == "Regulus"
        assert entry.hip_number == 49669

    def test_hip_prefix_with_space_spica(self):
        """Test HIP 65474 resolves to Spica."""
        entry, err = _resolve_star2("HIP 65474")
        assert err is None, f"Unexpected error: {err}"
        assert entry is not None
        assert entry.name == "Spica"
        assert entry.hip_number == 65474

    def test_hip_prefix_without_space(self):
        """Test HIP65474 (no space) resolves to Spica."""
        entry, err = _resolve_star2("HIP65474")
        assert err is None, f"Unexpected error: {err}"
        assert entry is not None
        assert entry.name == "Spica"
        assert entry.hip_number == 65474

    def test_hip_prefix_lowercase(self):
        """Test hip 49669 (lowercase) resolves to Regulus."""
        entry, err = _resolve_star2("hip 49669")
        assert err is None, f"Unexpected error: {err}"
        assert entry is not None
        assert entry.name == "Regulus"
        assert entry.hip_number == 49669

    def test_hip_prefix_mixed_case(self):
        """Test Hip 49669 (mixed case) resolves to Regulus."""
        entry, err = _resolve_star2("Hip 49669")
        assert err is None, f"Unexpected error: {err}"
        assert entry is not None
        assert entry.name == "Regulus"
        assert entry.hip_number == 49669

    def test_hip_prefix_with_whitespace(self):
        """Test HIP prefix with leading/trailing whitespace."""
        entry, err = _resolve_star2("  HIP 49669  ")
        assert err is None, f"Unexpected error: {err}"
        assert entry is not None
        assert entry.name == "Regulus"

    def test_hip_prefix_multiple_spaces(self):
        """Test HIP prefix with multiple spaces between HIP and number."""
        entry, err = _resolve_star2("HIP   49669")
        assert err is None, f"Unexpected error: {err}"
        assert entry is not None
        assert entry.name == "Regulus"

    def test_hip_prefix_unknown_number(self):
        """Test HIP prefix with unknown HIP number returns error."""
        entry, err = _resolve_star2("HIP 999999999")
        assert entry is None
        assert err is not None
        assert "999999999" in err

    def test_hip_prefix_with_invalid_number(self):
        """Test HIP prefix with non-numeric value doesn't match HIP pattern."""
        # "HIP abc" should not match the HIP pattern and fall through to other searches
        entry, err = _resolve_star2("HIP abc")
        # This should attempt to find a star named "HIP abc" and fail
        assert entry is None
        assert err is not None


@pytest.mark.unit
class TestHipPrefixWithFixstar2:
    """Tests for HIP prefix search via fixstar2_ut."""

    def test_swe_fixstar2_ut_hip_prefix_regulus(self):
        """Test fixstar2_ut with HIP 49669 returns Regulus."""
        jd = 2451545.0  # J2000.0
        pos, name, retflag = fixstar2_ut("HIP 49669", jd, 0)
        assert "Regulus" in name
        assert pos[0] > 0  # Longitude should be positive

    def test_swe_fixstar2_ut_hip_prefix_spica(self):
        """Test fixstar2_ut with HIP 65474 returns Spica."""
        jd = 2451545.0  # J2000.0
        pos, name, retflag = fixstar2_ut("HIP 65474", jd, 0)
        assert "Spica" in name
        assert pos[0] > 0

    def test_swe_fixstar2_ut_hip_prefix_no_space(self):
        """Test fixstar2_ut with HIP65474 (no space) returns Spica."""
        jd = 2451545.0  # J2000.0
        pos, name, retflag = fixstar2_ut("HIP65474", jd, 0)
        assert "Spica" in name


@pytest.mark.unit
class TestHipPrefixCatalogStars:
    """Tests for HIP prefix with various catalog stars."""

    @pytest.mark.parametrize(
        "hip_number,expected_name",
        [
            (49669, "Regulus"),
            (65474, "Spica"),
            (21421, "Aldebaran"),
            (80763, "Antares"),
            (32349, "Sirius"),
            (91262, "Vega"),
            (69673, "Arcturus"),
            (24608, "Capella"),
            (37279, "Procyon"),
            (37826, "Pollux"),
            (11767, "Polaris"),
            (27989, "Betelgeuse"),
            (24436, "Rigel"),
        ],
    )
    def test_hip_prefix_major_stars(self, hip_number, expected_name):
        """Test HIP prefix for major catalog stars."""
        entry, err = _resolve_star2(f"HIP {hip_number}")
        assert err is None, f"Error for HIP {hip_number}: {err}"
        assert entry is not None, f"No entry for HIP {hip_number}"
        assert entry.name == expected_name
        assert entry.hip_number == hip_number


@pytest.mark.unit
class TestHipPrefixConsistency:
    """Tests for consistency between HIP prefix and numeric search."""

    def test_hip_prefix_matches_numeric_search(self):
        """Test that HIP prefix search returns same result as numeric search."""
        # Search with HIP prefix
        entry_hip, err_hip = _resolve_star2("HIP 49669")
        # Search with just number
        entry_num, err_num = _resolve_star2("49669")

        assert err_hip is None
        assert err_num is None
        assert entry_hip is not None
        assert entry_num is not None
        assert entry_hip.name == entry_num.name
        assert entry_hip.hip_number == entry_num.hip_number

    def test_hip_prefix_matches_comma_prefix_search(self):
        """Test that HIP prefix search returns same result as comma prefix search."""
        # Search with HIP prefix
        entry_hip, err_hip = _resolve_star2("HIP 65474")
        # Search with comma prefix (pyswisseph HIP format)
        entry_comma, err_comma = _resolve_star2(",65474")

        assert err_hip is None
        assert err_comma is None
        assert entry_hip is not None
        assert entry_comma is not None
        assert entry_hip.name == entry_comma.name
        assert entry_hip.hip_number == entry_comma.hip_number
