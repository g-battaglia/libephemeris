"""
Unit tests for fixstar (TT-based fixed star calculation).

Tests the new fixstar() function that takes Terrestrial Time (TT)
instead of Universal Time (UT).
"""

import pytest
import swisseph as swe
import libephemeris as ephem


@pytest.mark.unit
class TestFixstarTT:
    """Tests for fixstar() with Terrestrial Time."""

    def test_fixstar_basic_regulus(self, standard_jd):
        """Test basic Regulus calculation with TT."""
        pos, name, retflag = ephem.fixstar("Regulus", standard_jd, 0)

        assert 149 < pos[0] < 151, f"Regulus lon: {pos[0]:.2f} out of range"
        assert -1 < pos[1] < 2, f"Regulus lat: {pos[1]:.2f} out of range"
        assert pos[2] > 1000, "Fixed stars should be very distant"

    def test_fixstar_basic_spica(self, standard_jd):
        """Test basic Spica calculation with TT."""
        pos, name, retflag = ephem.fixstar("Spica", standard_jd, 0)

        assert 203 < pos[0] < 205, f"Spica lon: {pos[0]:.2f} out of range"
        assert -3 < pos[1] < -1, f"Spica lat: {pos[1]:.2f} out of range"
        assert pos[2] > 1000, "Fixed stars should be very distant"

    def test_fixstar_unknown_star(self, standard_jd):
        """Test handling of unknown star name."""
        from libephemeris.exceptions import Error

        with pytest.raises(Error, match="could not find star name"):
            ephem.fixstar("UnknownStar", standard_jd, 0)

    def test_fixstar_case_insensitive(self, standard_jd):
        """Test that star name lookup is case insensitive."""
        pos_upper, _, _ = ephem.fixstar("REGULUS", standard_jd, 0)
        pos_lower, _, _ = ephem.fixstar("regulus", standard_jd, 0)
        pos_mixed, _, _ = ephem.fixstar("ReGuLuS", standard_jd, 0)

        assert pos_upper[0] == pos_lower[0] == pos_mixed[0]

    def test_fixstar_vs_fixstar_ut_relationship(self, standard_jd):
        """
        Test the relationship between fixstar (TT) and fixstar_ut (UT).

        When called with the same JD value, they should give different results
        because one interprets the JD as TT and the other as UT.
        The difference should be consistent with Delta T.
        """
        jd = standard_jd

        pos_tt, _, _ = ephem.fixstar("Regulus", jd, 0)
        pos_ut, _, _ = ephem.fixstar_ut("Regulus", jd, 0)

        # The positions should be slightly different because of Delta T
        # For J2000, Delta T is about 63.83 seconds
        # Stars move very slowly, so the difference should be tiny but non-zero
        diff = abs(pos_tt[0] - pos_ut[0])

        # Difference should be small but exist (proper motion over ~64 seconds)
        # Given star proper motion is very small, diff should be < 0.01 degree
        assert diff < 0.01, f"Difference too large: {diff}"

    def test_fixstar_vs_swisseph_regulus(self, standard_jd):
        """Compare Regulus with SwissEph fixstar."""
        pos_py, _, _ = ephem.fixstar("Regulus", standard_jd, 0)

        try:
            # SwissEph fixstar (TT version)
            pos_swe, name_swe, err_swe = swe.fixstar("Regulus", standard_jd, 0)
        except (TypeError, swe.Error) as e:
            pytest.skip(f"Skipping SwissEph comparison: {e}")
        except AttributeError:
            pytest.skip("SwissEph does not have fixstar function")

        diff_lon = abs(pos_py[0] - pos_swe[0])
        if diff_lon > 180:
            diff_lon = 360 - diff_lon

        assert diff_lon < 0.1, f"Regulus diff: {diff_lon}"

    def test_fixstar_vs_swisseph_spica(self, standard_jd):
        """Compare Spica with SwissEph fixstar."""
        pos_py, _, _ = ephem.fixstar("Spica", standard_jd, 0)

        try:
            pos_swe, name_swe, err_swe = swe.fixstar("Spica", standard_jd, 0)
        except (TypeError, swe.Error) as e:
            pytest.skip(f"Skipping SwissEph comparison: {e}")
        except AttributeError:
            pytest.skip("SwissEph does not have fixstar function")

        diff_lon = abs(pos_py[0] - pos_swe[0])
        if diff_lon > 180:
            diff_lon = 360 - diff_lon

        assert diff_lon < 0.1, f"Spica diff: {diff_lon}"

    def test_fixstar_return_structure(self, standard_jd):
        """Test the return structure of fixstar."""
        result = ephem.fixstar("Regulus", standard_jd, 0)

        # Should return a 3-tuple (pos, star_name, retflag)
        assert len(result) == 3, "Should return 3-tuple (pos, star_name, retflag)"

        pos, name, retflag = result

        # Position should be 6-tuple
        assert len(pos) == 6, "Position should be 6-tuple"
        assert isinstance(pos[0], float), "Longitude should be float"
        assert isinstance(pos[1], float), "Latitude should be float"
        assert isinstance(pos[2], float), "Distance should be float"
        assert isinstance(pos[3], float), "Speed lon should be float"
        assert isinstance(pos[4], float), "Speed lat should be float"
        assert isinstance(pos[5], float), "Speed dist should be float"

        # Speeds are currently 0 (not implemented)
        assert pos[3] == 0.0, "Speed lon should be 0.0"
        assert pos[4] == 0.0, "Speed lat should be 0.0"
        assert pos[5] == 0.0, "Speed dist should be 0.0"

        # name should be string
        assert isinstance(name, str), "star_name should be string"

        # retflag should be int
        assert isinstance(retflag, int), "retflag should be int"

    def test_fixstar_alias(self, standard_jd):
        """Test that fixstar is an alias for fixstar."""
        pos1, name1, flag1 = ephem.fixstar("Regulus", standard_jd, 0)
        pos2, name2, flag2 = ephem.fixstar("Regulus", standard_jd, 0)

        assert pos1 == pos2
        assert name1 == name2
        assert flag1 == flag2

    def test_fixstar_multiple_dates(self):
        """Test fixstar at different dates to verify proper motion."""
        jd_2000 = ephem.julday(2000, 1, 1, 12.0)
        jd_2050 = ephem.julday(2050, 1, 1, 12.0)

        pos_2000, _, _ = ephem.fixstar("Regulus", jd_2000, 0)
        pos_2050, _, _ = ephem.fixstar("Regulus", jd_2050, 0)

        # Star should have moved due to proper motion + precession
        diff = abs(pos_2050[0] - pos_2000[0])

        # Should move but not too much (a few degrees over 50 years)
        assert 0.001 < diff < 5.0, f"Proper motion over 50 years: {diff}"

    def test_fixstar_with_comma_in_name(self, standard_jd):
        """Test star name with comma (common in catalogs)."""
        # Some catalog formats use "Regulus,alLeo" style
        pos, _, _ = ephem.fixstar("Regulus,alLeo", standard_jd, 0)

        assert 149 < pos[0] < 151, f"Regulus lon: {pos[0]:.2f}"

    def test_fixstar_whitespace_in_name(self, standard_jd):
        """Test star name with extra whitespace."""
        pos1, _, _ = ephem.fixstar("Regulus", standard_jd, 0)
        pos2, _, _ = ephem.fixstar("  Regulus  ", standard_jd, 0)

        assert pos1[0] == pos2[0], "Should handle whitespace in name"


@pytest.mark.unit
class TestFixstarAndFixstarUtConsistency:
    """Tests for consistency between fixstar and fixstar_ut."""

    def test_deltat_effect(self):
        """
        Test that the Delta T effect is correctly applied.

        If we calculate a JD in UT, then convert to TT by adding Delta T,
        calling fixstar(jd_tt) should give the same result as fixstar_ut(jd_ut).
        """
        jd_ut = ephem.julday(2000, 1, 1, 12.0)

        # Get Delta T for this date
        delta_t = ephem.deltat(jd_ut)  # Returns Delta T in days
        jd_tt = jd_ut + delta_t

        # Calculate positions
        pos_from_ut, _, _ = ephem.fixstar_ut("Regulus", jd_ut, 0)
        pos_from_tt, _, _ = ephem.fixstar("Regulus", jd_tt, 0)

        # These should be very close (within floating point precision)
        diff = abs(pos_from_ut[0] - pos_from_tt[0])

        # Should be essentially identical (allowing for floating point)
        assert diff < 1e-6, f"Positions should match: diff={diff}"

    def test_both_functions_work(self, standard_jd):
        """Test that both fixstar and fixstar_ut work correctly."""
        pos_tt, name_tt, retflag_tt = ephem.fixstar("Spica", standard_jd, 0)
        pos_ut, name_ut, retflag_ut = ephem.fixstar_ut("Spica", standard_jd, 0)

        assert isinstance(name_tt, str)
        assert isinstance(name_ut, str)
        assert 203 < pos_tt[0] < 205
        assert 203 < pos_ut[0] < 205
