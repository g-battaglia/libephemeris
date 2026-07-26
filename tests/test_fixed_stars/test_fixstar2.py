"""
Unit tests for fixstar2 and fixstar2_ut functions.

Tests the v2-family fixed star lookup functions that support:
- Exact traditional star name lookup (case- and whitespace-insensitive)
- Nomenclature after a comma (",alLeo"), exact and case-sensitive
- Sequential catalog numbers ("12")
- Trailing-'%' prefix wildcard ("Alde%")
- Returns the full star name with the position

Bare nomenclature, bare HIP numbers, Bayer word forms and implicit
partial matching do NOT resolve (measured reference behavior).
"""

import pytest
import libephemeris as ephem
from libephemeris.exceptions import Error


@pytest.mark.unit
class TestFixstar2TT:
    """Tests for fixstar2() with Terrestrial Time."""

    def test_fixstar2_exact_name_regulus(self, standard_jd):
        """Test exact name lookup for Regulus."""
        pos, name, retflag = ephem.fixstar2("Regulus", standard_jd, 0)

        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"
        assert 149 < pos[0] < 151, f"Regulus lon: {pos[0]:.2f} out of range"
        assert -1 < pos[1] < 2, f"Regulus lat: {pos[1]:.2f} out of range"
        assert pos[2] > 1000, "Fixed stars should be very distant"

    def test_fixstar2_exact_name_spica(self, standard_jd):
        """Test exact name lookup for Spica."""
        pos, name, retflag = ephem.fixstar2("Spica", standard_jd, 0)

        assert name == "Spica,alVir", f"Expected 'Spica,alVir', got '{name}'"
        assert 203 < pos[0] < 205, f"Spica lon: {pos[0]:.2f} out of range"
        assert -3 < pos[1] < -1, f"Spica lat: {pos[1]:.2f} out of range"

    def test_fixstar2_case_insensitive(self, standard_jd):
        """Test that star name lookup is case insensitive."""
        pos_upper, name_upper, _ = ephem.fixstar2("REGULUS", standard_jd, 0)
        pos_lower, name_lower, _ = ephem.fixstar2("regulus", standard_jd, 0)
        pos_mixed, name_mixed, _ = ephem.fixstar2("ReGuLuS", standard_jd, 0)

        assert name_upper == name_lower == name_mixed == "Regulus,alLeo"
        assert pos_upper[0] == pos_lower[0] == pos_mixed[0]

    def test_fixstar2_no_bare_hip_number(self, standard_jd):
        """A bare number is a sequential index, never a HIP lookup."""
        # HIP 49669 (Regulus) is far beyond the catalog size, so the
        # sequential interpretation fails loudly instead of resolving.
        with pytest.raises(Error, match="sequential fixed star number"):
            ephem.fixstar2("49669", standard_jd, 0)

    def test_fixstar2_no_hip_number_with_comma(self, standard_jd):
        """A comma form keys on the nomenclature, never on a HIP number."""
        with pytest.raises(Error, match="could not find star name ,65474"):
            ephem.fixstar2(",65474", standard_jd, 0)

    def test_fixstar2_nomenclature(self, standard_jd):
        """Nomenclature resolves only AFTER a comma; bare nomen errors."""
        pos, name, retflag = ephem.fixstar2(",alLeo", standard_jd, 0)
        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"

        with pytest.raises(Error, match="could not find star name alleo"):
            ephem.fixstar2("alLeo", standard_jd, 0)

    def test_fixstar2_nomenclature_case_sensitive(self, standard_jd):
        """The nomenclature key is matched case-sensitively."""
        _, name1, _ = ephem.fixstar2(",alVir", standard_jd, 0)
        assert name1 == "Spica,alVir"

        for wrong_case in (",ALVIR", ",AlViR"):
            with pytest.raises(Error, match="could not find star name"):
                ephem.fixstar2(wrong_case, standard_jd, 0)

    def test_fixstar2_no_implicit_partial_name(self, standard_jd):
        """A partial name errors; the explicit '%' wildcard resolves."""
        with pytest.raises(Error, match="could not find star name reg"):
            ephem.fixstar2("Reg", standard_jd, 0)

        pos, name, retflag = ephem.fixstar2("Regu%", standard_jd, 0)
        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"

    def test_fixstar2_wildcard_spica(self, standard_jd):
        """Trailing-'%' prefix lookup for Spica."""
        pos, name, retflag = ephem.fixstar2("Spi%", standard_jd, 0)

        assert name == "Spica,alVir", f"Expected 'Spica,alVir', got '{name}'"

    def test_fixstar2_unknown_star(self, standard_jd):
        """Test handling of unknown star name."""
        with pytest.raises(Error):
            ephem.fixstar2("UnknownStar", standard_jd, 0)

    def test_fixstar2_unknown_hip_number(self, standard_jd):
        """Test handling of unknown HIP number."""
        with pytest.raises(Error):
            ephem.fixstar2("99999", standard_jd, 0)

    def test_fixstar2_empty_string(self, standard_jd):
        """Test handling of empty star name."""
        with pytest.raises(Error):
            ephem.fixstar2("", standard_jd, 0)

    def test_fixstar2_whitespace_handling(self, standard_jd):
        """Test star name with extra whitespace."""
        pos1, name1, _ = ephem.fixstar2("Regulus", standard_jd, 0)
        pos2, name2, _ = ephem.fixstar2("  Regulus  ", standard_jd, 0)

        assert name1 == name2, "Should handle whitespace in name"
        assert pos1[0] == pos2[0], "Positions should match"

    def test_fixstar2_comma_format_input(self, standard_jd):
        """Test input with comma (catalog format)."""
        pos, name, retflag = ephem.fixstar2("Regulus,alLeo", standard_jd, 0)

        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"

    def test_fixstar2_return_structure(self, standard_jd):
        """Test the return structure of fixstar2."""
        result = ephem.fixstar2("Regulus", standard_jd, 0)

        # Should return a 3-tuple (pos, name, retflag)
        assert len(result) == 3, "Should return 3-tuple (pos, name, retflag)"

        pos, name, retflag = result

        # Position should be 6-tuple
        assert len(pos) == 6, "Position should be 6-tuple"
        assert all(isinstance(p, float) for p in pos), "All positions should be floats"

        # Name should be string
        assert isinstance(name, str), "Name should be string"
        assert name != "", "Name should not be empty on success"

        # retflag should be int
        assert isinstance(retflag, int), "retflag should be int"

    def test_fixstar2_vs_fixstar_consistency(self, standard_jd):
        """Test that fixstar2 and fixstar give same position for same star."""
        pos2, name2, _ = ephem.fixstar2("Regulus", standard_jd, 0)
        pos1, _, _ = ephem.fixstar("Regulus", standard_jd, 0)

        assert pos1[0] == pos2[0], "Longitude should match"
        assert pos1[1] == pos2[1], "Latitude should match"
        assert pos1[2] == pos2[2], "Distance should match"

    def test_fixstar2_alias(self, standard_jd):
        """Test that fixstar2 is an alias for fixstar2."""
        pos1, name1, flag1 = ephem.fixstar2("Regulus", standard_jd, 0)
        pos2, name2, flag2 = ephem.fixstar2("Regulus", standard_jd, 0)

        assert name1 == name2
        assert pos1 == pos2
        assert flag1 == flag2


@pytest.mark.unit
class TestFixstar2UT:
    """Tests for fixstar2_ut() with Universal Time."""

    def test_fixstar2_ut_basic(self, standard_jd):
        """Test basic fixstar2_ut functionality."""
        pos, name, retflag = ephem.fixstar2_ut("Regulus", standard_jd, 0)

        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"
        assert 149 < pos[0] < 151, f"Regulus lon: {pos[0]:.2f} out of range"

    def test_fixstar2_ut_nomenclature(self, standard_jd):
        """Test fixstar2_ut with comma-nomenclature lookup."""
        pos, name, retflag = ephem.fixstar2_ut(",alVir", standard_jd, 0)

        assert name == "Spica,alVir", f"Expected 'Spica,alVir', got '{name}'"
        assert 203 < pos[0] < 205, f"Spica lon: {pos[0]:.2f} out of range"

    def test_fixstar2_ut_wildcard_name(self, standard_jd):
        """Test fixstar2_ut with the explicit '%' prefix wildcard."""
        pos, name, retflag = ephem.fixstar2_ut("Spi%", standard_jd, 0)

        assert name == "Spica,alVir"

    def test_fixstar2_ut_unknown_star(self, standard_jd):
        """Test fixstar2_ut handling of unknown star."""
        with pytest.raises(Error):
            ephem.fixstar2_ut("UnknownStar", standard_jd, 0)

    def test_fixstar2_ut_alias(self, standard_jd):
        """Test that fixstar2_ut is an alias for fixstar2_ut."""
        pos1, name1, flag1 = ephem.fixstar2_ut("Spica", standard_jd, 0)
        pos2, name2, flag2 = ephem.fixstar2_ut("Spica", standard_jd, 0)

        assert name1 == name2
        assert pos1 == pos2
        assert flag1 == flag2

    def test_fixstar2_ut_vs_fixstar_ut_consistency(self, standard_jd):
        """Test that fixstar2_ut and fixstar_ut give same position for same star."""
        pos2, name2, _ = ephem.fixstar2_ut("Spica", standard_jd, 0)
        pos1, _, _ = ephem.fixstar_ut("Spica", standard_jd, 0)

        assert pos1[0] == pos2[0], "Longitude should match"
        assert pos1[1] == pos2[1], "Latitude should match"

    def test_fixstar2_ut_vs_tt_relationship(self, standard_jd):
        """
        Test the relationship between fixstar2 (TT) and fixstar2_ut (UT).

        When called with the same JD value, they should give slightly different
        results because of Delta T.
        """
        jd = standard_jd

        pos_tt, name_tt, _ = ephem.fixstar2("Regulus", jd, 0)
        pos_ut, name_ut, _ = ephem.fixstar2_ut("Regulus", jd, 0)

        # Names should be the same
        assert name_tt == name_ut == "Regulus,alLeo"

        # Positions should be slightly different due to Delta T
        diff = abs(pos_tt[0] - pos_ut[0])
        assert diff < 0.01, f"Difference too large: {diff}"


@pytest.mark.unit
class TestFixstar2AndFixstar2UtConsistency:
    """Tests for consistency between fixstar2 and fixstar2_ut."""

    def test_deltat_effect(self):
        """
        Test that the Delta T effect is correctly applied.

        If we calculate a JD in UT, then convert to TT by adding Delta T,
        calling fixstar2(jd_tt) should give the same result as fixstar2_ut(jd_ut).
        """
        jd_ut = ephem.julday(2000, 1, 1, 12.0)

        # Get Delta T for this date
        delta_t = ephem.deltat(jd_ut)  # Returns Delta T in days
        jd_tt = jd_ut + delta_t

        # Calculate positions
        pos_from_ut, name_ut, _ = ephem.fixstar2_ut("Regulus", jd_ut, 0)
        pos_from_tt, name_tt, _ = ephem.fixstar2("Regulus", jd_tt, 0)

        # Names should match
        assert name_ut == name_tt == "Regulus,alLeo"

        # Positions should be essentially identical (allowing for floating point)
        diff = abs(pos_from_ut[0] - pos_from_tt[0])
        assert diff < 1e-6, f"Positions should match: diff={diff}"

    def test_both_functions_return_names(self, standard_jd):
        """Test that both functions return proper star names."""
        pos_tt, name_tt, _ = ephem.fixstar2("Spica", standard_jd, 0)
        pos_ut, name_ut, _ = ephem.fixstar2_ut("Spica", standard_jd, 0)

        assert name_tt == name_ut == "Spica,alVir"
        assert 203 < pos_tt[0] < 205
        assert 203 < pos_ut[0] < 205

    def test_multiple_dates_proper_motion(self):
        """Test that positions change with proper motion over time."""
        jd_2000 = ephem.julday(2000, 1, 1, 12.0)
        jd_2050 = ephem.julday(2050, 1, 1, 12.0)

        pos_2000, name_2000, _ = ephem.fixstar2("Regulus", jd_2000, 0)
        pos_2050, name_2050, _ = ephem.fixstar2("Regulus", jd_2050, 0)

        assert name_2000 == name_2050 == "Regulus,alLeo"

        # Star should have moved due to proper motion + precession
        diff = abs(pos_2050[0] - pos_2000[0])
        assert 0.001 < diff < 5.0, f"Proper motion over 50 years: {diff}"
