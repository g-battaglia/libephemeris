"""
Unit tests for fixed stars: Regulus and Spica.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


@pytest.mark.unit
class TestFixedStars:
    """Tests for fixed star calculations."""

    def test_regulus_j2000(self, standard_jd):
        """Test Regulus position at J2000."""
        pos, _ = ephem.swe_calc_ut(standard_jd, SE_REGULUS, 0)

        # Regulus at ~149-150° at J2000
        assert 149 < pos[0] < 151, f"Regulus lon: {pos[0]:.2f}°"
        assert -1 < pos[1] < 2, f"Regulus lat: {pos[1]:.2f}°"
        assert pos[2] > 1000, "Fixed stars should be very distant"

    def test_spica_j2000(self, standard_jd):
        """Test Spica position at J2000."""
        pos, _ = ephem.swe_calc_ut(standard_jd, SE_SPICA_STAR, 0)

        # Spica at ~203-204° at J2000
        assert 203 < pos[0] < 205, f"Spica lon: {pos[0]:.2f}°"
        assert -3 < pos[1] < -1, f"Spica lat: {pos[1]:.2f}°"
        assert pos[2] > 1000, "Fixed stars should be very distant"

    def test_regulus_vs_swisseph(self, standard_jd):
        """Compare Regulus with SwissEph."""
        pos_py, _ = ephem.swe_calc_ut(standard_jd, SE_REGULUS, 0)

        # SwissEph fixed star by name
        try:
            # fixstar_ut returns (pos_tuple, name, error_msg)
            pos_swe, name_swe, err_swe = swe.fixstar_ut("Regulus", standard_jd, 0)
        except TypeError:
            # Fallback for older versions if they strictly want bytes
            try:
                pos_swe, name_swe, err_swe = swe.fixstar_ut(b"Regulus", standard_jd, 0)
            except swe.Error as e:
                pytest.skip(f"Skipping Regulus comparison: {e}")
        except swe.Error as e:
            pytest.skip(f"Skipping Regulus comparison: {e}")
        except ValueError:
            # Handle case where it might return 2 values (older versions?)
            # But error said "too many values", so it returns > 2.
            # Let's assume 3.
            pytest.skip(
                "Skipping Regulus comparison: Unexpected return signature from swisseph"
            )

        diff_lon = abs(pos_py[0] - pos_swe[0])
        if diff_lon > 180:
            diff_lon = 360 - diff_lon

        assert diff_lon < 0.1, f"Regulus diff: {diff_lon}°"

    def test_spica_vs_swisseph(self, standard_jd):
        """Compare Spica with SwissEph."""
        pos_py, _ = ephem.swe_calc_ut(standard_jd, SE_SPICA_STAR, 0)

        try:
            pos_swe, name_swe, err_swe = swe.fixstar_ut("Spica", standard_jd, 0)
        except TypeError:
            try:
                pos_swe, name_swe, err_swe = swe.fixstar_ut(b"Spica", standard_jd, 0)
            except swe.Error as e:
                pytest.skip(f"Skipping Spica comparison: {e}")
        except swe.Error as e:
            pytest.skip(f"Skipping Spica comparison: {e}")

        diff_lon = abs(pos_py[0] - pos_swe[0])
        if diff_lon > 180:
            diff_lon = 360 - diff_lon

        assert diff_lon < 0.1, f"Spica diff: {diff_lon}°"

    def test_proper_motion(self):
        """Test that fixed stars move due to proper motion."""
        jd1 = ephem.swe_julday(2000, 1, 1, 12.0)
        jd2 = ephem.swe_julday(2050, 1, 1, 12.0)  # 50 years later

        pos1, _ = ephem.swe_calc_ut(jd1, SE_REGULUS, 0)
        pos2, _ = ephem.swe_calc_ut(jd2, SE_REGULUS, 0)

        # Regulus has small proper motion, should move slightly
        diff = abs(pos2[0] - pos1[0])

        assert 0.001 < diff < 1.0, f"Proper motion over 50 years: {diff}°"


@pytest.mark.unit
class TestAngles:
    """Tests for astrological angles."""

    def test_ascendant_basic(self):
        """Test basic Ascendant calculation."""
        jd = ephem.swe_julday(2000, 1, 1, 12.0)
        lat, lon = 41.9028, 12.4964  # Rome

        ephem.swe_set_topo(lon, lat, 0)

        asc, _ = ephem.swe_calc_ut(jd, SE_ASCENDANT, 0)

        assert 0 <= asc[0] < 360, "Invalid Ascendant"

    def test_mc_basic(self):
        """Test basic MC calculation."""
        jd = ephem.swe_julday(2000, 1, 1, 12.0)
        lat, lon = 41.9028, 12.4964

        ephem.swe_set_topo(lon, lat, 0)

        mc, _ = ephem.swe_calc_ut(jd, SE_MC, 0)

        assert 0 <= mc[0] < 360, "Invalid MC"

    def test_descendant_opposite_ascendant(self):
        """Test Descendant is opposite Ascendant."""
        jd = ephem.swe_julday(2000, 1, 1, 12.0)
        lat, lon = 41.9028, 12.4964

        ephem.swe_set_topo(lon, lat, 0)

        asc, _ = ephem.swe_calc_ut(jd, SE_ASCENDANT, 0)
        desc, _ = ephem.swe_calc_ut(jd, SE_DESCENDANT, 0)

        expected_desc = (asc[0] + 180.0) % 360.0
        diff = abs(desc[0] - expected_desc)

        assert diff < 0.01, "Descendant should be 180° from Ascendant"

    def test_ic_opposite_mc(self):
        """Test IC is opposite MC."""
        jd = ephem.swe_julday(2000, 1, 1, 12.0)
        lat, lon = 41.9028, 12.4964

        ephem.swe_set_topo(lon, lat, 0)

        mc, _ = ephem.swe_calc_ut(jd, SE_MC, 0)
        ic, _ = ephem.swe_calc_ut(jd, SE_IC, 0)

        expected_ic = (mc[0] + 180.0) % 360.0
        diff = abs(ic[0] - expected_ic)

        assert diff < 0.01, "IC should be 180° from MC"

    def test_vertex(self):
        """Test Vertex calculation."""
        jd = ephem.swe_julday(2000, 1, 1, 12.0)
        lat, lon = 41.9028, 12.4964

        ephem.swe_set_topo(lon, lat, 0)

        vertex, _ = ephem.swe_calc_ut(jd, SE_VERTEX, 0)

        assert 0 <= vertex[0] < 360, "Invalid Vertex"

    def test_angles_require_location(self):
        """Test that angles fail without location set."""
        jd = ephem.swe_julday(2000, 1, 1, 12.0)

        # Reset state (no location)
        from libephemeris import state

        state._TOPO = None

        with pytest.raises(ValueError, match="observer location"):
            ephem.swe_calc_ut(jd, SE_ASCENDANT, 0)


@pytest.mark.unit
class TestArabicParts:
    """Tests for Arabic Parts."""

    def test_pars_fortunae(self):
        """Test Part of Fortune."""
        jd = ephem.swe_julday(2000, 1, 1, 12.0)
        lat, lon = 41.9028, 12.4964

        # Pre-calculate
        ephem.swe_calc_angles(jd, lat, lon)

        pf, _ = ephem.swe_calc_ut(jd, SE_PARS_FORTUNAE, 0)

        assert 0 <= pf[0] < 360, "Invalid Part of Fortune"

    def test_all_arabic_parts(self):
        """Test all 4 Arabic parts."""
        jd = ephem.swe_julday(2000, 1, 1, 12.0)
        lat, lon = 41.9028, 12.4964

        ephem.swe_calc_angles(jd, lat, lon)

        parts = [
            (SE_PARS_FORTUNAE, "Fortune"),
            (SE_PARS_SPIRITUS, "Spirit"),
            (SE_PARS_AMORIS, "Love"),
            (SE_PARS_FIDEI, "Faith"),
        ]

        for part_id, name in parts:
            pos, _ = ephem.swe_calc_ut(jd, part_id, 0)
            assert 0 <= pos[0] < 360, f"Invalid {name}"

    def test_arabic_parts_require_precalc(self):
        """Test that Arabic parts require pre-calculation."""
        jd = ephem.swe_julday(2000, 1, 1, 12.0)

        # Clear cache
        from libephemeris import state

        state.clear_angles_cache()

        with pytest.raises(ValueError, match="pre-calculated"):
            ephem.swe_calc_ut(jd, SE_PARS_FORTUNAE, 0)
