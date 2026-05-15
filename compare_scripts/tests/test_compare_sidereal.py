"""
Sidereal/Ayanamsha Comparison Tests.

Compares all ayanamsha (sidereal) modes between pyswisseph and libephemeris.
Tests all 43 ayanamsha systems across different dates and planets.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    JUPITER,
    SATURN,
    SIDM_USER,
    SIDM_FAGAN_BRADLEY,
    SIDM_LAHIRI,
    SIDM_DELUCE,
    SIDM_RAMAN,
    SIDM_USHASHASHI,
    SIDM_KRISHNAMURTI,
    SIDM_DJWHAL_KHUL,
    SIDM_YUKTESHWAR,
    SIDM_JN_BHASIN,
    SIDM_BABYL_KUGLER1,
    SIDM_BABYL_KUGLER2,
    SIDM_BABYL_KUGLER3,
    SIDM_BABYL_HUBER,
    SIDM_BABYL_ETPSC,
    SIDM_ALDEBARAN_15TAU,
    SIDM_HIPPARCHOS,
    SIDM_SASSANIAN,
    SIDM_GALCENT_0SAG,
    SIDM_J2000,
    SIDM_J1900,
    SIDM_B1950,
    SIDM_SURYASIDDHANTA,
    SIDM_SURYASIDDHANTA_MSUN,
    SIDM_ARYABHATA,
    SIDM_ARYABHATA_MSUN,
    SIDM_SS_REVATI,
    SIDM_SS_CITRA,
    SIDM_TRUE_CITRA,
    SIDM_TRUE_REVATI,
    SIDM_TRUE_PUSHYA,
    SIDM_GALCENT_RGILBRAND,
    SIDM_GALEQU_IAU1958,
    SIDM_GALEQU_TRUE,
    SIDM_GALEQU_MULA,
    SIDM_GALALIGN_MARDYKS,
    SIDM_TRUE_MULA,
    SIDM_GALCENT_MULA_WILHELM,
    SIDM_ARYABHATA_522,
    SIDM_BABYL_BRITTON,
    SIDM_TRUE_SHEORAN,
    SIDM_GALCENT_COCHRANE,
    SIDM_GALEQU_FIORENZA,
    SIDM_VALENS_MOON,
    FLG_SIDEREAL,
)


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# AYANAMSHA MODES
# ============================================================================

FORMULA_BASED_AYANAMSHAS = [
    (SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
    (SIDM_LAHIRI, "Lahiri"),
    (SIDM_DELUCE, "De Luce"),
    (SIDM_RAMAN, "Raman"),
    (SIDM_USHASHASHI, "Ushashashi"),
    (SIDM_KRISHNAMURTI, "Krishnamurti"),
    (SIDM_DJWHAL_KHUL, "Djwhal Khul"),
    (SIDM_YUKTESHWAR, "Yukteshwar"),
    (SIDM_JN_BHASIN, "JN Bhasin"),
    (SIDM_BABYL_KUGLER1, "Babylonian Kugler 1"),
    (SIDM_BABYL_KUGLER2, "Babylonian Kugler 2"),
    (SIDM_BABYL_KUGLER3, "Babylonian Kugler 3"),
    (SIDM_BABYL_HUBER, "Babylonian Huber"),
    (SIDM_BABYL_ETPSC, "Babylonian ETPSC"),
    (SIDM_ALDEBARAN_15TAU, "Aldebaran 15 Tau"),
    (SIDM_HIPPARCHOS, "Hipparchos"),
    (SIDM_SASSANIAN, "Sassanian"),
    (SIDM_J2000, "J2000"),
    (SIDM_J1900, "J1900"),
    (SIDM_B1950, "B1950"),
    (SIDM_SURYASIDDHANTA, "Suryasiddhanta"),
    (SIDM_SURYASIDDHANTA_MSUN, "Suryasiddhanta Mean Sun"),
    (SIDM_ARYABHATA, "Aryabhata"),
    (SIDM_ARYABHATA_MSUN, "Aryabhata Mean Sun"),
    (SIDM_SS_REVATI, "SS Revati"),
    (SIDM_SS_CITRA, "SS Citra"),
    (SIDM_ARYABHATA_522, "Aryabhata 522"),
    (SIDM_BABYL_BRITTON, "Babylonian Britton"),
]

STAR_BASED_AYANAMSHAS = [
    (SIDM_TRUE_CITRA, "True Citra"),
    (SIDM_TRUE_REVATI, "True Revati"),
    (SIDM_TRUE_PUSHYA, "True Pushya"),
    (SIDM_TRUE_MULA, "True Mula"),
    (SIDM_TRUE_SHEORAN, "True Sheoran"),
    (SIDM_GALCENT_0SAG, "Galactic Center 0 Sag"),
    (SIDM_GALCENT_RGILBRAND, "Galactic Center Rgilbrand"),
    (SIDM_GALCENT_MULA_WILHELM, "Galactic Center Mula Wilhelm"),
    (SIDM_GALCENT_COCHRANE, "Galactic Center Cochrane"),
    (SIDM_GALEQU_IAU1958, "Galactic Equator IAU 1958"),
    (SIDM_GALEQU_TRUE, "Galactic Equator True"),
    (SIDM_GALEQU_MULA, "Galactic Equator Mula"),
    (SIDM_GALEQU_FIORENZA, "Galactic Equator Fiorenza"),
    (SIDM_GALALIGN_MARDYKS, "Galactic Alignment Mardyks"),
    (SIDM_VALENS_MOON, "Valens Moon"),
]

ALL_AYANAMSHAS = FORMULA_BASED_AYANAMSHAS + STAR_BASED_AYANAMSHAS

TEST_DATES = [
    (2451545.0, "J2000.0"),
    (2415020.5, "1900-01-01"),
    (2433282.5, "1950-01-01"),
    (2469807.5, "2050-01-01"),
]

TEST_PLANETS = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
]

# Tolerances
STRICT_TOLERANCE = 0.001  # 3.6 arcseconds - for formula-based
RELAXED_TOLERANCE = 0.1  # 6 arcminutes - for star-based


def get_tolerance(sid_mode: int) -> float:
    """Get appropriate tolerance for ayanamsha mode."""
    star_modes = {m[0] for m in STAR_BASED_AYANAMSHAS}
    if sid_mode in star_modes:
        return RELAXED_TOLERANCE
    return STRICT_TOLERANCE


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestAyanamshaValues:
    """Compare ayanamsha values between implementations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("sid_mode,sid_name", FORMULA_BASED_AYANAMSHAS)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    def test_formula_based_ayanamsha(self, sid_mode, sid_name, jd, date_desc):
        """Test formula-based ayanamsha values match within strict tolerance."""
        swe.set_sid_mode(sid_mode)
        ephem.set_sid_mode(sid_mode)

        ayan_swe = swe.get_ayanamsa_ut(jd)
        ayan_py = ephem.get_ayanamsa_ut(jd)

        diff = angular_diff(ayan_swe, ayan_py)

        assert diff < STRICT_TOLERANCE, (
            f"{sid_name} at {date_desc}: diff {diff:.6f}° exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("sid_mode,sid_name", STAR_BASED_AYANAMSHAS)
    @pytest.mark.parametrize(
        "jd,date_desc", TEST_DATES[:2]
    )  # Fewer dates for star-based
    def test_star_based_ayanamsha(self, sid_mode, sid_name, jd, date_desc):
        """Test star-based ayanamsha values match within relaxed tolerance."""
        swe.set_sid_mode(sid_mode)
        ephem.set_sid_mode(sid_mode)

        ayan_swe = swe.get_ayanamsa_ut(jd)
        ayan_py = ephem.get_ayanamsa_ut(jd)

        diff = angular_diff(ayan_swe, ayan_py)

        assert diff < RELAXED_TOLERANCE, (
            f"{sid_name} at {date_desc}: diff {diff:.6f}° exceeds tolerance"
        )


class TestSiderealPlanets:
    """Compare sidereal planetary positions."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("sid_mode,sid_name", FORMULA_BASED_AYANAMSHAS[:10])
    @pytest.mark.parametrize("planet_id,planet_name", TEST_PLANETS)
    def test_sidereal_planet_positions(
        self, sid_mode, sid_name, planet_id, planet_name
    ):
        """Test sidereal planetary positions match."""
        jd = 2451545.0  # J2000

        swe.set_sid_mode(sid_mode)
        ephem.set_sid_mode(sid_mode)

        pos_swe, _ = swe.calc_ut(jd, planet_id, FLG_SIDEREAL)
        pos_py, _ = ephem.calc_ut(jd, planet_id, FLG_SIDEREAL)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])

        tolerance = get_tolerance(sid_mode)

        assert diff_lon < tolerance, (
            f"{planet_name} sidereal ({sid_name}): diff {diff_lon:.6f}° exceeds tolerance"
        )


class TestMajorAyanamshas:
    """Detailed tests for most commonly used ayanamshas."""

    @pytest.mark.comparison
    def test_lahiri_at_j2000(self):
        """Test Lahiri ayanamsha at J2000 epoch."""
        jd = 2451545.0

        swe.set_sid_mode(SIDM_LAHIRI)
        ephem.set_sid_mode(SIDM_LAHIRI)

        ayan_swe = swe.get_ayanamsa_ut(jd)
        ayan_py = ephem.get_ayanamsa_ut(jd)

        diff = abs(ayan_swe - ayan_py)

        # Lahiri should be approximately 23.86 degrees at J2000
        assert 23.5 < ayan_py < 24.5, f"Lahiri value {ayan_py}° out of expected range"
        assert diff < 0.001, f"Lahiri diff {diff:.6f}° exceeds tight tolerance"

    @pytest.mark.comparison
    def test_fagan_bradley_at_j2000(self):
        """Test Fagan-Bradley ayanamsha at J2000 epoch."""
        jd = 2451545.0

        swe.set_sid_mode(SIDM_FAGAN_BRADLEY)
        ephem.set_sid_mode(SIDM_FAGAN_BRADLEY)

        ayan_swe = swe.get_ayanamsa_ut(jd)
        ayan_py = ephem.get_ayanamsa_ut(jd)

        diff = abs(ayan_swe - ayan_py)

        # Fagan-Bradley should be approximately 24.74 degrees at J2000
        assert 24.0 < ayan_py < 25.5, (
            f"Fagan-Bradley value {ayan_py}° out of expected range"
        )
        assert diff < 0.001, f"Fagan-Bradley diff {diff:.6f}° exceeds tight tolerance"

    @pytest.mark.comparison
    def test_raman_at_j2000(self):
        """Test Raman ayanamsha at J2000 epoch."""
        jd = 2451545.0

        swe.set_sid_mode(SIDM_RAMAN)
        ephem.set_sid_mode(SIDM_RAMAN)

        ayan_swe = swe.get_ayanamsa_ut(jd)
        ayan_py = ephem.get_ayanamsa_ut(jd)

        diff = abs(ayan_swe - ayan_py)

        # Raman should be approximately 22.4 degrees at J2000
        assert 22.0 < ayan_py < 23.0, f"Raman value {ayan_py}° out of expected range"
        assert diff < 0.001, f"Raman diff {diff:.6f}° exceeds tight tolerance"


class TestAyanamshaProgression:
    """Test that ayanamsha changes correctly over time."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("sid_mode,sid_name", FORMULA_BASED_AYANAMSHAS[:5])
    def test_ayanamsha_increases_over_century(self, sid_mode, sid_name):
        """Test that ayanamsha increases over a century (precession)."""
        jd_1900 = 2415020.5
        jd_2000 = 2451545.0

        swe.set_sid_mode(sid_mode)
        ephem.set_sid_mode(sid_mode)

        # Get ayanamsha at both epochs from both libraries
        ayan_1900_swe = swe.get_ayanamsa_ut(jd_1900)
        ayan_2000_swe = swe.get_ayanamsa_ut(jd_2000)
        ayan_1900_py = ephem.get_ayanamsa_ut(jd_1900)
        ayan_2000_py = ephem.get_ayanamsa_ut(jd_2000)

        # Change should be positive (precession increases ayanamsha)
        change_swe = ayan_2000_swe - ayan_1900_swe
        change_py = ayan_2000_py - ayan_1900_py

        # Skip J2000 mode which is 0 at J2000 by definition
        if sid_mode != SIDM_J2000:
            assert change_swe > 0.5, f"{sid_name}: expected precession over century"
            assert change_py > 0.5, f"{sid_name}: expected precession over century"

        # Changes should be similar
        change_diff = abs(change_swe - change_py)
        assert change_diff < 0.01, (
            f"{sid_name}: century change differs by {change_diff:.6f}°"
        )


class TestJ2000Mode:
    """Test J2000 ayanamsha mode specifically."""

    @pytest.mark.comparison
    def test_j2000_is_zero_at_epoch(self):
        """J2000 ayanamsha should be exactly 0 at J2000 epoch."""
        jd = 2451545.0

        swe.set_sid_mode(SIDM_J2000)
        ephem.set_sid_mode(SIDM_J2000)

        ayan_swe = swe.get_ayanamsa_ut(jd)
        ayan_py = ephem.get_ayanamsa_ut(jd)

        assert abs(ayan_swe) < 0.0001, f"SWE J2000 at J2000 = {ayan_swe}, expected 0"
        assert abs(ayan_py) < 0.0001, f"PY J2000 at J2000 = {ayan_py}, expected 0"

    @pytest.mark.comparison
    def test_j2000_changes_over_time(self):
        """J2000 ayanamsha should change from the J2000 epoch."""
        swe.set_sid_mode(SIDM_J2000)
        ephem.set_sid_mode(SIDM_J2000)

        ayan_1900 = ephem.get_ayanamsa_ut(2415020.5)
        ayan_2000 = ephem.get_ayanamsa_ut(2451545.0)
        ayan_2100 = ephem.get_ayanamsa_ut(2488069.5)

        ayan_1900_swe = swe.get_ayanamsa_ut(2415020.5)
        ayan_2100_swe = swe.get_ayanamsa_ut(2488069.5)

        assert abs(ayan_1900 - ayan_1900_swe) < 0.001, (
            f"J2000 at 1900: lib={ayan_1900}, swe={ayan_1900_swe}"
        )
        assert abs(ayan_2000) < 0.001, f"J2000 at 2000 should be ~0, got {ayan_2000}"
        assert abs(ayan_2100 - ayan_2100_swe) < 0.001, (
            f"J2000 at 2100: lib={ayan_2100}, swe={ayan_2100_swe}"
        )


class TestUserDefinedAyanamsha:
    """
    Test SIDM_USER (255) custom ayanamsha mode.

    This mode allows users to define custom ayanamsha with t0 (reference time)
    and ayan_t0 (ayanamsha value at reference time) parameters.
    Many Vedic astrologers use custom ayanamsha values.
    """

    @pytest.mark.comparison
    def test_user_ayanamsha_value_at_j2000(self):
        """Test SIDM_USER ayanamsha value matches pyswisseph at J2000."""
        jd = 2451545.0  # J2000
        t0 = jd
        ayan_t0 = 23.5  # Custom ayanamsha value

        # Set user-defined mode with same parameters
        swe.set_sid_mode(SIDM_USER, t0, ayan_t0)
        ephem.set_sid_mode(SIDM_USER, t0=t0, ayan_t0=ayan_t0)

        ayan_swe = swe.get_ayanamsa_ut(jd)
        ayan_py = ephem.get_ayanamsa_ut(jd)

        diff = abs(ayan_swe - ayan_py)

        # At t0, both should return exactly ayan_t0
        assert abs(ayan_py - ayan_t0) < 0.001, f"Expected {ayan_t0}, got {ayan_py}"
        assert diff < STRICT_TOLERANCE, (
            f"SIDM_USER at J2000: diff {diff:.6f}° exceeds tolerance"
        )

    @pytest.mark.comparison
    def test_user_ayanamsha_at_multiple_dates(self):
        """Test SIDM_USER ayanamsha matches pyswisseph across multiple dates."""
        t0 = 2451545.0  # J2000
        ayan_t0 = 23.5

        swe.set_sid_mode(SIDM_USER, t0, ayan_t0)
        ephem.set_sid_mode(SIDM_USER, t0=t0, ayan_t0=ayan_t0)

        for jd, date_desc in TEST_DATES:
            ayan_swe = swe.get_ayanamsa_ut(jd)
            ayan_py = ephem.get_ayanamsa_ut(jd)

            diff = abs(ayan_swe - ayan_py)

            assert diff < STRICT_TOLERANCE, (
                f"SIDM_USER at {date_desc}: diff {diff:.6f}° exceeds tolerance "
                f"(swe={ayan_swe:.6f}, py={ayan_py:.6f})"
            )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", TEST_PLANETS)
    def test_user_sidereal_planet_positions(self, planet_id, planet_name):
        """Test sidereal planetary positions with SIDM_USER match pyswisseph."""
        jd = 2451545.0  # J2000
        t0 = jd
        ayan_t0 = 23.5

        swe.set_sid_mode(SIDM_USER, t0, ayan_t0)
        ephem.set_sid_mode(SIDM_USER, t0=t0, ayan_t0=ayan_t0)

        pos_swe, _ = swe.calc_ut(jd, planet_id, FLG_SIDEREAL)
        pos_py, _ = ephem.calc_ut(jd, planet_id, FLG_SIDEREAL)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])

        assert diff_lon < STRICT_TOLERANCE, (
            f"{planet_name} sidereal (SIDM_USER): diff {diff_lon:.6f}° exceeds tolerance"
        )

    @pytest.mark.comparison
    def test_changing_parameters_changes_results(self):
        """Verify that changing t0/ayan_t0 parameters changes results."""
        jd = 2451545.0  # J2000

        # First set of parameters
        t0_1 = jd
        ayan_t0_1 = 23.5

        ephem.set_sid_mode(SIDM_USER, t0=t0_1, ayan_t0=ayan_t0_1)
        ayan_1 = ephem.get_ayanamsa_ut(jd)

        # Second set with different ayan_t0
        t0_2 = jd
        ayan_t0_2 = 25.0

        ephem.set_sid_mode(SIDM_USER, t0=t0_2, ayan_t0=ayan_t0_2)
        ayan_2 = ephem.get_ayanamsa_ut(jd)

        # Results should be different
        diff = abs(ayan_1 - ayan_2)
        assert diff > 1.0, (
            f"Changing ayan_t0 should change result: got {ayan_1:.6f} vs {ayan_2:.6f}"
        )

        # Third set with different t0
        t0_3 = 2415020.5  # J1900
        ayan_t0_3 = 23.5

        ephem.set_sid_mode(SIDM_USER, t0=t0_3, ayan_t0=ayan_t0_3)
        ayan_3 = ephem.get_ayanamsa_ut(jd)

        # Should differ from first (same ayan_t0 but different t0)
        diff_t0 = abs(ayan_1 - ayan_3)
        assert diff_t0 > 0.1, (
            f"Changing t0 should change result: got {ayan_1:.6f} vs {ayan_3:.6f}"
        )

    @pytest.mark.comparison
    def test_user_ayanamsha_different_reference_epochs(self):
        """Test SIDM_USER with different reference epochs matches pyswisseph."""
        # Test with t0 at J1900
        t0 = 2415020.5  # J1900
        ayan_t0 = 22.5
        jd = 2451545.0  # Test at J2000

        swe.set_sid_mode(SIDM_USER, t0, ayan_t0)
        ephem.set_sid_mode(SIDM_USER, t0=t0, ayan_t0=ayan_t0)

        ayan_swe = swe.get_ayanamsa_ut(jd)
        ayan_py = ephem.get_ayanamsa_ut(jd)

        diff = abs(ayan_swe - ayan_py)

        assert diff < STRICT_TOLERANCE, (
            f"SIDM_USER with J1900 reference: diff {diff:.6f}° exceeds tolerance"
        )

    @pytest.mark.comparison
    def test_user_sidereal_positions_at_future_date(self):
        """Test SIDM_USER sidereal positions at future date match pyswisseph."""
        t0 = 2451545.0  # J2000
        ayan_t0 = 23.5
        jd = 2469807.5  # 2050-01-01

        swe.set_sid_mode(SIDM_USER, t0, ayan_t0)
        ephem.set_sid_mode(SIDM_USER, t0=t0, ayan_t0=ayan_t0)

        for planet_id, planet_name in TEST_PLANETS[:2]:  # Sun and Moon
            pos_swe, _ = swe.calc_ut(jd, planet_id, FLG_SIDEREAL)
            pos_py, _ = ephem.calc_ut(jd, planet_id, FLG_SIDEREAL)

            diff_lon = angular_diff(pos_swe[0], pos_py[0])

            assert diff_lon < STRICT_TOLERANCE, (
                f"{planet_name} sidereal (SIDM_USER) at 2050: "
                f"diff {diff_lon:.6f}° exceeds tolerance"
            )
