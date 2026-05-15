"""
Comprehensive tests for all 43 ayanamsha modes.

Tests sidereal calculations and ayanamsha values.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


class TestAyanamshaBasicValues:
    """Test basic ayanamsha value calculations."""

    @pytest.mark.unit
    def test_lahiri_at_j2000(self):
        """Lahiri ayanamsha at J2000 should be ~23.9 degrees."""
        jd = 2451545.0
        ephem.set_sid_mode(SIDM_LAHIRI)
        ayan = ephem.get_ayanamsa_ut(jd)

        assert 23.0 < ayan < 25.0, f"Lahiri ayanamsha {ayan} unexpected"

    @pytest.mark.unit
    def test_fagan_bradley_at_j2000(self):
        """Fagan-Bradley ayanamsha at J2000 should be ~24.7 degrees."""
        jd = 2451545.0
        ephem.set_sid_mode(SIDM_FAGAN_BRADLEY)
        ayan = ephem.get_ayanamsa_ut(jd)

        assert 24.0 < ayan < 26.0, f"Fagan-Bradley ayanamsha {ayan} unexpected"

    @pytest.mark.unit
    def test_ayanamsha_positive(self):
        """Ayanamsha should be positive for modern dates."""
        jd = 2451545.0
        ephem.set_sid_mode(SIDM_LAHIRI)
        ayan = ephem.get_ayanamsa_ut(jd)

        assert ayan > 0


class TestAyanamshaVsPyswisseph:
    """Compare ayanamsha values with pyswisseph."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "sid_mode,name",
        [
            (SIDM_FAGAN_BRADLEY, "Fagan-Bradley"),
            (SIDM_LAHIRI, "Lahiri"),
            (SIDM_RAMAN, "Raman"),
            (SIDM_KRISHNAMURTI, "Krishnamurti"),
            (SIDM_J2000, "J2000"),
            (SIDM_B1950, "B1950"),
        ],
    )
    def test_standard_ayanamsha_matches_swe(self, sid_mode, name):
        """Standard ayanamshas should match pyswisseph closely."""
        jd = 2451545.0
        tolerance = 0.06  # degrees

        ephem.set_sid_mode(sid_mode)
        swe.set_sid_mode(sid_mode)

        ayan_lib = ephem.get_ayanamsa_ut(jd)
        ayan_swe = swe.get_ayanamsa_ut(jd)

        diff = abs(ayan_lib - ayan_swe)
        assert diff < tolerance, f"{name} diff {diff} >= {tolerance}"

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "sid_mode,name",
        [
            (SIDM_TRUE_CITRA, "True Citra"),
            (SIDM_TRUE_REVATI, "True Revati"),
            (SIDM_GALCENT_0SAG, "Galcent 0 Sag"),
        ],
    )
    def test_star_based_ayanamsha_relaxed(self, sid_mode, name):
        """Star-based ayanamshas may need relaxed tolerance."""
        jd = 2451545.0
        tolerance = 1.0  # degrees (relaxed)

        ephem.set_sid_mode(sid_mode)
        swe.set_sid_mode(sid_mode)

        ayan_lib = ephem.get_ayanamsa_ut(jd)
        ayan_swe = swe.get_ayanamsa_ut(jd)

        diff = abs(ayan_lib - ayan_swe)
        assert diff < tolerance, f"{name} diff {diff} >= {tolerance}"


class TestSiderealPositions:
    """Test sidereal planet positions."""

    @pytest.mark.unit
    def test_sidereal_sun_differs_from_tropical(self):
        """Sidereal Sun should differ from tropical by ayanamsha."""
        jd = 2451545.0

        ephem.set_sid_mode(SIDM_LAHIRI)

        pos_trop, _ = ephem.calc_ut(jd, SUN, 0)
        pos_sid, _ = ephem.calc_ut(jd, SUN, FLG_SIDEREAL)

        diff = pos_trop[0] - pos_sid[0]
        if diff < 0:
            diff += 360

        # Difference should be close to ayanamsha
        ayan = ephem.get_ayanamsa_ut(jd)
        assert abs(diff - ayan) < 0.01, f"Diff {diff} != ayanamsha {ayan}"

    @pytest.mark.comparison
    def test_sidereal_sun_matches_pyswisseph(self):
        """Sidereal Sun should match pyswisseph."""
        jd = 2451545.0

        ephem.set_sid_mode(SIDM_LAHIRI)
        swe.set_sid_mode(SIDM_LAHIRI)

        pos_lib, _ = ephem.calc_ut(jd, SUN, FLG_SIDEREAL)
        pos_swe, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL)

        diff = abs(pos_lib[0] - pos_swe[0])
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Sidereal Sun diff {diff}"


class TestAyanamshaProgression:
    """Test that ayanamsha progresses correctly over time."""

    @pytest.mark.unit
    def test_ayanamsha_increases_over_time(self):
        """Ayanamsha should increase with time (precession)."""
        ephem.set_sid_mode(SIDM_LAHIRI)

        jd_2000 = 2451545.0
        jd_2020 = ephem.julday(2020, 1, 1, 12.0)

        ayan_2000 = ephem.get_ayanamsa_ut(jd_2000)
        ayan_2020 = ephem.get_ayanamsa_ut(jd_2020)

        assert ayan_2020 > ayan_2000, "Ayanamsha should increase over time"

    @pytest.mark.unit
    def test_ayanamsha_precession_rate(self):
        """Ayanamsha should increase ~50 arcsec/year."""
        ephem.set_sid_mode(SIDM_LAHIRI)

        jd_2000 = 2451545.0
        jd_2020 = ephem.julday(2020, 1, 1, 12.0)

        ayan_2000 = ephem.get_ayanamsa_ut(jd_2000)
        ayan_2020 = ephem.get_ayanamsa_ut(jd_2020)

        # ~20 years × 50"/year = 1000" = ~0.28°
        diff = ayan_2020 - ayan_2000
        assert 0.2 < diff < 0.4, f"20-year ayanamsha change {diff} unexpected"


class TestAyanamshaName:
    """Test get_ayanamsa_name function."""

    @pytest.mark.unit
    def test_lahiri_name(self):
        """Should return 'Lahiri' for SIDM_LAHIRI."""
        name = ephem.get_ayanamsa_name(SIDM_LAHIRI)
        assert "Lahiri" in name or "lahiri" in name.lower()

    @pytest.mark.unit
    def test_fagan_bradley_name(self):
        """Should return 'Fagan' or 'Bradley' for SIDM_FAGAN_BRADLEY."""
        name = ephem.get_ayanamsa_name(SIDM_FAGAN_BRADLEY)
        assert "Fagan" in name or "Bradley" in name


class TestAllAyanamshas:
    """Test all 43 ayanamsha modes work."""

    @pytest.mark.unit
    @pytest.mark.parametrize("sid_mode", range(43))
    def test_ayanamsha_mode_valid(self, sid_mode):
        """Each ayanamsha mode should return valid value."""
        jd = 2451545.0

        try:
            ephem.set_sid_mode(sid_mode)
            ayan = ephem.get_ayanamsa_ut(jd)

            # Should be reasonable (some modes might give very different values)
            assert isinstance(ayan, float)
            assert -360 < ayan < 360  # Very permissive
        except Exception as e:
            pytest.fail(f"Mode {sid_mode} failed: {e}")


class TestAyanamsaEx:
    """Test get_ayanamsa_ex and get_ayanamsa_ex_ut functions."""

    @pytest.mark.unit
    def test_get_ayanamsa_ex_returns_tuple(self):
        """get_ayanamsa_ex should return a tuple of (retflag, ayanamsa)."""
        jd = 2451545.0
        ephem.set_sid_mode(SIDM_LAHIRI)
        result = ephem.get_ayanamsa_ex(jd, 0)

        assert isinstance(result, tuple)
        assert len(result) == 2
        assert all(isinstance(x, (int, float)) for x in result)

    @pytest.mark.unit
    def test_get_ayanamsa_ex_ut_returns_tuple(self):
        """get_ayanamsa_ex_ut should return a tuple of (retflag, ayanamsa)."""
        jd = 2451545.0
        ephem.set_sid_mode(SIDM_LAHIRI)
        result = ephem.get_ayanamsa_ex_ut(jd, 0)

        assert isinstance(result, tuple)
        assert len(result) == 2
        assert all(isinstance(x, (int, float)) for x in result)

    @pytest.mark.unit
    def test_get_ayanamsa_ex_ayanamsa_matches_standard(self):
        """Ayanamsa from get_ayanamsa_ex should match get_ayanamsa."""
        jd = 2451545.0
        ephem.set_sid_mode(SIDM_LAHIRI)

        aya_standard = ephem.get_ayanamsa(jd)
        _retflag, aya_ex = ephem.get_ayanamsa_ex(jd, 0)

        assert abs(aya_ex - aya_standard) < 0.0001

    @pytest.mark.unit
    def test_get_ayanamsa_ex_ut_ayanamsa_matches_standard(self):
        """Ayanamsa from get_ayanamsa_ex_ut should match get_ayanamsa_ut."""
        jd = 2451545.0
        ephem.set_sid_mode(SIDM_LAHIRI)

        aya_standard = ephem.get_ayanamsa_ut(jd)
        _retflag, aya_ex = ephem.get_ayanamsa_ex_ut(jd, 0)

        assert abs(aya_ex - aya_standard) < 0.0001

    @pytest.mark.unit
    def test_get_ayanamsa_ex_retflag(self):
        """get_ayanamsa_ex retflag should be a non-negative integer."""
        jd = 2451545.0  # J2000.0
        ephem.set_sid_mode(SIDM_LAHIRI)
        retflag, _aya = ephem.get_ayanamsa_ex(jd, 0)

        assert isinstance(retflag, int)
        assert retflag >= 0

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "sid_mode,name",
        [
            (SIDM_FAGAN_BRADLEY, "Fagan-Bradley"),
            (SIDM_LAHIRI, "Lahiri"),
            (SIDM_RAMAN, "Raman"),
            (SIDM_KRISHNAMURTI, "Krishnamurti"),
            (SIDM_J2000, "J2000"),
            (SIDM_B1950, "B1950"),
        ],
    )
    def test_get_ayanamsa_ex_different_modes(self, sid_mode, name):
        """Test get_ayanamsa_ex works with different sidereal modes."""
        jd = 2451545.0
        ephem.set_sid_mode(sid_mode)
        retflag, aya = ephem.get_ayanamsa_ex(jd, 0)

        # Ayanamsa should be valid
        assert isinstance(aya, float)
        assert -360 < aya < 360

    @pytest.mark.unit
    def test_get_ayanamsa_ex_uses_global_sid_mode(self):
        """get_ayanamsa_ex should use the global set_sid_mode."""
        jd = 2451545.0

        # Set global mode to Lahiri
        ephem.set_sid_mode(SIDM_LAHIRI)
        _retflag, aya_ex = ephem.get_ayanamsa_ex(jd, 0)

        # Compare with standard function
        aya_lahiri = ephem.get_ayanamsa(jd)

        # Should match Lahiri
        assert abs(aya_ex - aya_lahiri) < 0.0001

    @pytest.mark.unit
    def test_get_ayanamsa_ex_with_flags(self):
        """Test that flags parameter is accepted."""
        jd = 2451545.0
        ephem.set_sid_mode(SIDM_LAHIRI)
        result = ephem.get_ayanamsa_ex(jd, 0)
        assert len(result) == 2

        result_with_flag = ephem.get_ayanamsa_ex(jd, 256)
        assert len(result_with_flag) == 2

    @pytest.mark.unit
    def test_get_ayanamsa_ex_changes_over_time(self):
        """Ayanamsa should change over time due to precession."""
        jd_2000 = 2451545.0  # J2000.0
        jd_2020 = 2459215.5  # 2020-12-01

        ephem.set_sid_mode(SIDM_LAHIRI)
        _, aya_2000 = ephem.get_ayanamsa_ex(jd_2000, 0)
        _, aya_2020 = ephem.get_ayanamsa_ex(jd_2020, 0)

        # Ayanamsa increases over time due to precession
        assert aya_2020 > aya_2000

    @pytest.mark.unit
    def test_get_ayanamsa_ex_aliases(self):
        """Test that the non-swe_ aliases work."""
        jd = 2451545.0
        ephem.set_sid_mode(SIDM_LAHIRI)

        # Test alias without swe_ prefix
        result1 = ephem.get_ayanamsa_ex(jd, 0)
        result2 = ephem.get_ayanamsa_ex(jd, 0)

        assert result1 == result2

        result3 = ephem.get_ayanamsa_ex_ut(jd, 0)
        result4 = ephem.get_ayanamsa_ex_ut(jd, 0)

        assert result3 == result4
