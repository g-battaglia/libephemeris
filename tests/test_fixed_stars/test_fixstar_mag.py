"""
Unit tests for fixstar_mag and fixstar2_mag functions.

Tests the magnitude lookup functions that return only the visual magnitude
of a fixed star without calculating position.
"""

import pytest
import libephemeris as ephem
from libephemeris.exceptions import Error


@pytest.mark.unit
class TestFixstarMag:
    """Tests for fixstar_mag()."""

    def test_fixstar_mag_regulus(self):
        """Test magnitude lookup for Regulus."""
        mag, name = ephem.fixstar_mag("Regulus")

        assert mag == pytest.approx(1.40, abs=0.01), f"Regulus magnitude: {mag}"
        assert name == "Regulus,alLeo"

    def test_fixstar_mag_spica(self):
        """Test magnitude lookup for Spica."""
        mag, name = ephem.fixstar_mag("Spica")

        assert mag == pytest.approx(1.04, abs=0.01), f"Spica magnitude: {mag}"
        assert name == "Spica,alVir"

    def test_fixstar_mag_case_insensitive(self):
        """Test that star name lookup is case insensitive."""
        mag_upper, _ = ephem.fixstar_mag("REGULUS")
        mag_lower, _ = ephem.fixstar_mag("regulus")
        mag_mixed, _ = ephem.fixstar_mag("ReGuLuS")

        assert mag_upper == mag_lower == mag_mixed

    def test_fixstar_mag_unknown_star(self):
        """Test handling of unknown star name."""
        with pytest.raises(Error):
            ephem.fixstar_mag("UnknownStar")

    def test_fixstar_mag_with_comma_in_name(self):
        """Test star name with comma (common in catalogs)."""
        mag, name = ephem.fixstar_mag("Regulus,alLeo")

        assert mag == pytest.approx(1.40, abs=0.01)

    def test_fixstar_mag_whitespace_in_name(self):
        """Test star name with extra whitespace."""
        mag1, _ = ephem.fixstar_mag("Regulus")
        mag2, _ = ephem.fixstar_mag("  Regulus  ")

        assert mag1 == mag2, "Should handle whitespace in name"

    def test_fixstar_mag_return_structure(self):
        """Test the return structure of fixstar_mag.

        The reference ephemeris fixstar_mag returns (magnitude, star_name) tuple.
        """
        result = ephem.fixstar_mag("Regulus")

        assert isinstance(result, tuple), "Should return a tuple"
        assert len(result) == 2, "Should return 2-tuple (magnitude, star_name)"

        mag, name = result
        assert isinstance(mag, float), "Magnitude should be float"
        assert isinstance(name, str), "Name should be string"

    def test_fixstar_mag_alias(self):
        """Test that fixstar_mag is an alias for fixstar_mag."""
        mag1, name1 = ephem.fixstar_mag("Regulus")
        mag2, name2 = ephem.fixstar_mag("Regulus")

        assert mag1 == mag2
        assert name1 == name2


@pytest.mark.unit
class TestFixstar2Mag:
    """Tests for fixstar2_mag()."""

    def test_fixstar2_mag_exact_name_regulus(self):
        """Test exact name lookup for Regulus magnitude."""
        mag, name = ephem.fixstar2_mag("Regulus")

        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"
        assert mag == pytest.approx(1.40, abs=0.01), f"Regulus magnitude: {mag}"

    def test_fixstar2_mag_exact_name_spica(self):
        """Test exact name lookup for Spica magnitude."""
        mag, name = ephem.fixstar2_mag("Spica")

        assert name == "Spica,alVir", f"Expected 'Spica,alVir', got '{name}'"
        assert mag == pytest.approx(1.04, abs=0.01), f"Spica magnitude: {mag}"

    def test_fixstar2_mag_case_insensitive(self):
        """Test that star name lookup is case insensitive."""
        mag_upper, name_upper = ephem.fixstar2_mag("REGULUS")
        mag_lower, name_lower = ephem.fixstar2_mag("regulus")
        mag_mixed, name_mixed = ephem.fixstar2_mag("ReGuLuS")

        assert name_upper == name_lower == name_mixed == "Regulus,alLeo"
        assert mag_upper == mag_lower == mag_mixed

    def test_fixstar2_mag_no_hip_number(self):
        """A bare number is a sequential index, never a HIP lookup."""
        with pytest.raises(Error, match="sequential fixed star number"):
            ephem.fixstar2_mag("49669")

    def test_fixstar2_mag_no_hip_number_with_comma(self):
        """A comma form keys on the nomenclature, never on a HIP number."""
        with pytest.raises(
            Error, match="no fixed star matches the search string ',65474'"
        ):
            ephem.fixstar2_mag(",65474")

    def test_fixstar2_mag_nomenclature(self):
        """Nomenclature resolves after a comma; the bare form errors."""
        mag, name = ephem.fixstar2_mag(",alLeo")

        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"
        assert mag == pytest.approx(1.40, abs=0.01)

        with pytest.raises(
            Error, match="no fixed star matches the search string 'alleo'"
        ):
            ephem.fixstar2_mag("alLeo")

    def test_fixstar2_mag_nomenclature_case_sensitive(self):
        """The nomenclature key is matched case-sensitively."""
        mag1, name1 = ephem.fixstar2_mag(",alVir")
        assert name1 == "Spica,alVir"

        for wrong_case in (",ALVIR", ",AlViR"):
            with pytest.raises(Error, match="no fixed star matches the search string"):
                ephem.fixstar2_mag(wrong_case)

    def test_fixstar2_mag_no_implicit_partial_name(self):
        """A partial name errors; the explicit '%' wildcard resolves."""
        with pytest.raises(
            Error, match="no fixed star matches the search string 'reg'"
        ):
            ephem.fixstar2_mag("Reg")

        mag, name = ephem.fixstar2_mag("Regu%")
        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"

    def test_fixstar2_mag_wildcard_spica(self):
        """Trailing-'%' prefix lookup for Spica."""
        mag, name = ephem.fixstar2_mag("Spi%")

        assert name == "Spica,alVir", f"Expected 'Spica,alVir', got '{name}'"

    def test_fixstar2_mag_unknown_star(self):
        """Test handling of unknown star name."""
        with pytest.raises(Error, match="no fixed star matches the search string"):
            ephem.fixstar2_mag("UnknownStar")

    def test_fixstar2_mag_out_of_range_number(self):
        """An out-of-range sequential number errors."""
        with pytest.raises(Error, match="sequential fixed star number"):
            ephem.fixstar2_mag("99999")

    def test_fixstar2_mag_empty_string(self):
        """Test handling of empty star name."""
        with pytest.raises(Error, match="the star name is empty"):
            ephem.fixstar2_mag("")

    def test_fixstar2_mag_whitespace_handling(self):
        """Test star name with extra whitespace."""
        mag1, name1 = ephem.fixstar2_mag("Regulus")
        mag2, name2 = ephem.fixstar2_mag("  Regulus  ")

        assert name1 == name2, "Should handle whitespace in name"
        assert mag1 == mag2, "Magnitudes should match"

    def test_fixstar2_mag_comma_format_input(self):
        """Test input with comma (catalog format)."""
        mag, name = ephem.fixstar2_mag("Regulus,alLeo")

        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"

    def test_fixstar2_mag_return_structure(self):
        """Test the return structure of fixstar2_mag."""
        result = ephem.fixstar2_mag("Regulus")

        # Should return a 2-tuple (magnitude, star_name)
        assert len(result) == 2, "Should return 2-tuple (magnitude, star_name)"

        mag, name = result

        # Magnitude should be float
        assert isinstance(mag, float), "Magnitude should be float"

        # Name should be string
        assert isinstance(name, str), "Name should be string"
        assert name != "", "Name should not be empty on success"

    def test_fixstar2_mag_vs_fixstar_mag_consistency(self):
        """Test that fixstar2_mag and fixstar_mag give same magnitude for same star."""
        mag2, name = ephem.fixstar2_mag("Regulus")
        mag1, _ = ephem.fixstar_mag("Regulus")

        assert mag1 == mag2, "Magnitude should match"

    def test_fixstar2_mag_alias(self):
        """Test that fixstar2_mag is an alias for fixstar2_mag."""
        mag1, name1 = ephem.fixstar2_mag("Regulus")
        mag2, name2 = ephem.fixstar2_mag("Regulus")

        assert name1 == name2
        assert mag1 == mag2


@pytest.mark.unit
class TestMagnitudeValues:
    """Tests to validate magnitude values are astronomically correct."""

    def test_regulus_magnitude_value(self):
        """Regulus is a 1st magnitude star (V=1.40)."""
        mag, _ = ephem.fixstar_mag("Regulus")
        # Regulus V magnitude is approximately 1.40
        assert 1.35 < mag < 1.45, f"Regulus magnitude {mag} not in expected range"

    def test_spica_magnitude_value(self):
        """Spica is a 1st magnitude star (V=1.04)."""
        mag, _ = ephem.fixstar_mag("Spica")
        # Spica V magnitude is approximately 1.04
        assert 0.99 < mag < 1.09, f"Spica magnitude {mag} not in expected range"

    def test_spica_brighter_than_regulus(self):
        """Spica (V=1.04) should be brighter than Regulus (V=1.40)."""
        mag_spica, _ = ephem.fixstar_mag("Spica")
        mag_regulus, _ = ephem.fixstar_mag("Regulus")

        # Lower magnitude means brighter
        assert mag_spica < mag_regulus, "Spica should be brighter than Regulus"
