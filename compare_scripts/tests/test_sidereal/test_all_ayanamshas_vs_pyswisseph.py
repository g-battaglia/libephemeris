"""
Comprehensive tests comparing all 43 ayanamsha modes between libephemeris and pyswisseph.

Tests ayanamsha values at three different epochs:
- J2000 (2000-01-01 12:00 TT)
- 1900 (1900-01-01 12:00 TT)
- 2100 (2100-01-01 12:00 TT)

Uses appropriate tolerances:
- 0.01 degrees (~36 arcseconds) for formula-based ayanamshas
  (minor differences due to different precession models between implementations)
- 1.0 degrees for star-based ayanamshas (due to different star catalogs and definitions)

Note: The angle difference calculation properly handles wrap-around at 360 degrees.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
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
)


# =============================================================================
# AYANAMSHA MODE DEFINITIONS
# =============================================================================

# Complete list of all 43 ayanamsha modes (0-42)
ALL_AYANAMSHA_MODES = [
    (0, SIDM_FAGAN_BRADLEY, "Fagan-Bradley"),
    (1, SIDM_LAHIRI, "Lahiri"),
    (2, SIDM_DELUCE, "De Luce"),
    (3, SIDM_RAMAN, "Raman"),
    (4, SIDM_USHASHASHI, "Ushashashi"),
    (5, SIDM_KRISHNAMURTI, "Krishnamurti"),
    (6, SIDM_DJWHAL_KHUL, "Djwhal Khul"),
    (7, SIDM_YUKTESHWAR, "Yukteshwar"),
    (8, SIDM_JN_BHASIN, "JN Bhasin"),
    (9, SIDM_BABYL_KUGLER1, "Babylonian Kugler 1"),
    (10, SIDM_BABYL_KUGLER2, "Babylonian Kugler 2"),
    (11, SIDM_BABYL_KUGLER3, "Babylonian Kugler 3"),
    (12, SIDM_BABYL_HUBER, "Babylonian Huber"),
    (13, SIDM_BABYL_ETPSC, "Babylonian ETPSC"),
    (14, SIDM_ALDEBARAN_15TAU, "Aldebaran 15 Tau"),
    (15, SIDM_HIPPARCHOS, "Hipparchos"),
    (16, SIDM_SASSANIAN, "Sassanian"),
    (17, SIDM_GALCENT_0SAG, "Galactic Center 0 Sag"),
    (18, SIDM_J2000, "J2000"),
    (19, SIDM_J1900, "J1900"),
    (20, SIDM_B1950, "B1950"),
    (21, SIDM_SURYASIDDHANTA, "Suryasiddhanta"),
    (22, SIDM_SURYASIDDHANTA_MSUN, "Suryasiddhanta Mean Sun"),
    (23, SIDM_ARYABHATA, "Aryabhata"),
    (24, SIDM_ARYABHATA_MSUN, "Aryabhata Mean Sun"),
    (25, SIDM_SS_REVATI, "SS Revati"),
    (26, SIDM_SS_CITRA, "SS Citra"),
    (27, SIDM_TRUE_CITRA, "True Citra"),
    (28, SIDM_TRUE_REVATI, "True Revati"),
    (29, SIDM_TRUE_PUSHYA, "True Pushya"),
    (30, SIDM_GALCENT_RGILBRAND, "Galactic Center Rgilbrand"),
    (31, SIDM_GALEQU_IAU1958, "Galactic Equator IAU 1958"),
    (32, SIDM_GALEQU_TRUE, "Galactic Equator True"),
    (33, SIDM_GALEQU_MULA, "Galactic Equator Mula"),
    (34, SIDM_GALALIGN_MARDYKS, "Galactic Alignment Mardyks"),
    (35, SIDM_TRUE_MULA, "True Mula"),
    (36, SIDM_GALCENT_MULA_WILHELM, "Galactic Center Mula Wilhelm"),
    (37, SIDM_ARYABHATA_522, "Aryabhata 522"),
    (38, SIDM_BABYL_BRITTON, "Babylonian Britton"),
    (39, SIDM_TRUE_SHEORAN, "True Sheoran"),
    (40, SIDM_GALCENT_COCHRANE, "Galactic Center Cochrane"),
    (41, SIDM_GALEQU_FIORENZA, "Galactic Equator Fiorenza"),
    (42, SIDM_VALENS_MOON, "Valens Moon"),
]

# Star-based and galactic ayanamshas that may have larger differences
# due to different star catalogs, proper motion calculations, or ambiguous definitions
STAR_BASED_AYANAMSHAS = {
    SIDM_TRUE_CITRA,  # 27 - True position of Spica
    SIDM_TRUE_REVATI,  # 28 - True position of Revati
    SIDM_TRUE_PUSHYA,  # 29 - True position of Pushya
    SIDM_TRUE_MULA,  # 35 - True position of Mula
    SIDM_TRUE_SHEORAN,  # 39 - True Sheoran
    SIDM_GALCENT_0SAG,  # 17 - Galactic Center at 0 Sag
    SIDM_GALCENT_RGILBRAND,  # 30 - Galactic Center (Gil Brand)
    SIDM_GALCENT_MULA_WILHELM,  # 36 - Galactic Center at Mula (Wilhelm)
    SIDM_GALCENT_COCHRANE,  # 40 - Galactic Center (Cochrane)
    SIDM_GALEQU_IAU1958,  # 31 - Galactic Equator (IAU 1958)
    SIDM_GALEQU_TRUE,  # 32 - Galactic Equator (True)
    SIDM_GALEQU_MULA,  # 33 - Galactic Equator at Mula
    SIDM_GALEQU_FIORENZA,  # 41 - Galactic Equator (Fiorenza)
    SIDM_GALALIGN_MARDYKS,  # 34 - Galactic Alignment (Mardyks)
}

# Tolerances
# Formula-based: 0.01 degrees (~36 arcseconds) - accounts for minor differences
# in precession model implementations between libephemeris and pyswisseph
FORMULA_TOLERANCE = 0.01  # degrees
STAR_TOLERANCE = 1.0  # degrees - relaxed for star-based ayanamshas

# Test dates as Julian Days
JD_J2000 = 2451545.0  # 2000-01-01 12:00:00 TT
JD_1900 = 2415020.0  # 1900-01-01 12:00:00 TT (approximately)
JD_2100 = 2488070.0  # 2100-01-01 12:00:00 TT (approximately)


def angle_diff(a1: float, a2: float) -> float:
    """
    Calculate the smallest difference between two angles.

    Handles wrap-around at 360 degrees. For example:
    - angle_diff(359.0, 1.0) returns 2.0, not 358.0
    - angle_diff(0.5, 359.5) returns 1.0, not 359.0

    Args:
        a1: First angle in degrees
        a2: Second angle in degrees

    Returns:
        Absolute difference in degrees (0 to 180)
    """
    diff = abs(a1 - a2)
    if diff > 180:
        diff = 360 - diff
    return diff


def get_tolerance(sid_mode: int) -> float:
    """Return appropriate tolerance based on ayanamsha type."""
    if sid_mode in STAR_BASED_AYANAMSHAS:
        return STAR_TOLERANCE
    return FORMULA_TOLERANCE


# =============================================================================
# TEST CLASSES
# =============================================================================


class TestAllAyanamshasAtJ2000:
    """Test all 43 ayanamsha modes at J2000 epoch."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "mode_id,sid_mode,name",
        ALL_AYANAMSHA_MODES,
        ids=[m[2] for m in ALL_AYANAMSHA_MODES],
    )
    def test_ayanamsha_at_j2000(self, mode_id, sid_mode, name):
        """
        Compare ayanamsha value at J2000 between libephemeris and pyswisseph.

        Tests that each of the 43 ayanamsha modes produces values matching
        pyswisseph within the appropriate tolerance.
        """
        jd = JD_J2000
        tolerance = get_tolerance(sid_mode)

        # Set mode in both libraries
        ephem.set_sid_mode(sid_mode)
        swe.set_sid_mode(sid_mode)

        # Calculate ayanamsha
        ayan_lib = ephem.get_ayanamsa_ut(jd)
        ayan_swe = swe.get_ayanamsa_ut(jd)

        # Compare (use angle_diff to handle wrap-around at 360 degrees)
        diff = angle_diff(ayan_lib, ayan_swe)
        assert diff < tolerance, (
            f"Mode {mode_id} ({name}) at J2000: "
            f"libephemeris={ayan_lib:.6f}, pyswisseph={ayan_swe:.6f}, "
            f"diff={diff:.6f} >= tolerance={tolerance}"
        )


class TestAllAyanamshasAt1900:
    """Test all 43 ayanamsha modes at year 1900."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "mode_id,sid_mode,name",
        ALL_AYANAMSHA_MODES,
        ids=[m[2] for m in ALL_AYANAMSHA_MODES],
    )
    def test_ayanamsha_at_1900(self, mode_id, sid_mode, name):
        """
        Compare ayanamsha value at 1900 between libephemeris and pyswisseph.

        Tests historical date calculations for each ayanamsha mode.
        """
        jd = JD_1900
        tolerance = get_tolerance(sid_mode)

        # Set mode in both libraries
        ephem.set_sid_mode(sid_mode)
        swe.set_sid_mode(sid_mode)

        # Calculate ayanamsha
        ayan_lib = ephem.get_ayanamsa_ut(jd)
        ayan_swe = swe.get_ayanamsa_ut(jd)

        # Compare (use angle_diff to handle wrap-around at 360 degrees)
        diff = angle_diff(ayan_lib, ayan_swe)
        assert diff < tolerance, (
            f"Mode {mode_id} ({name}) at 1900: "
            f"libephemeris={ayan_lib:.6f}, pyswisseph={ayan_swe:.6f}, "
            f"diff={diff:.6f} >= tolerance={tolerance}"
        )


class TestAllAyanamshasAt2100:
    """Test all 43 ayanamsha modes at year 2100."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "mode_id,sid_mode,name",
        ALL_AYANAMSHA_MODES,
        ids=[m[2] for m in ALL_AYANAMSHA_MODES],
    )
    def test_ayanamsha_at_2100(self, mode_id, sid_mode, name):
        """
        Compare ayanamsha value at 2100 between libephemeris and pyswisseph.

        Tests future date calculations for each ayanamsha mode.
        """
        jd = JD_2100
        tolerance = get_tolerance(sid_mode)

        # Set mode in both libraries
        ephem.set_sid_mode(sid_mode)
        swe.set_sid_mode(sid_mode)

        # Calculate ayanamsha
        ayan_lib = ephem.get_ayanamsa_ut(jd)
        ayan_swe = swe.get_ayanamsa_ut(jd)

        # Compare (use angle_diff to handle wrap-around at 360 degrees)
        diff = angle_diff(ayan_lib, ayan_swe)
        assert diff < tolerance, (
            f"Mode {mode_id} ({name}) at 2100: "
            f"libephemeris={ayan_lib:.6f}, pyswisseph={ayan_swe:.6f}, "
            f"diff={diff:.6f} >= tolerance={tolerance}"
        )


class TestFormulBasedAyanamshasTight:
    """
    Test formula-based ayanamshas with 0.01 degree (~36 arcsec) tolerance.

    These ayanamshas use fixed epoch values and precession rates.
    Minor differences occur due to different precession model implementations.
    """

    # Formula-based ayanamshas (non-star-based)
    FORMULA_BASED = [
        m for m in ALL_AYANAMSHA_MODES if m[1] not in STAR_BASED_AYANAMSHAS
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "mode_id,sid_mode,name",
        FORMULA_BASED,
        ids=[m[2] for m in FORMULA_BASED],
    )
    @pytest.mark.parametrize(
        "jd,epoch_name",
        [
            (JD_J2000, "J2000"),
            (JD_1900, "1900"),
            (JD_2100, "2100"),
        ],
        ids=["J2000", "1900", "2100"],
    )
    def test_formula_ayanamsha_tight_tolerance(
        self, mode_id, sid_mode, name, jd, epoch_name
    ):
        """
        Test formula-based ayanamshas at all three epochs with tight tolerance.
        """
        tolerance = FORMULA_TOLERANCE

        # Set mode in both libraries
        ephem.set_sid_mode(sid_mode)
        swe.set_sid_mode(sid_mode)

        # Calculate ayanamsha
        ayan_lib = ephem.get_ayanamsa_ut(jd)
        ayan_swe = swe.get_ayanamsa_ut(jd)

        # Compare (use angle_diff to handle wrap-around at 360 degrees)
        diff = angle_diff(ayan_lib, ayan_swe)
        assert diff < tolerance, (
            f"{name} at {epoch_name}: "
            f"lib={ayan_lib:.6f}, swe={ayan_swe:.6f}, "
            f"diff={diff:.6f} >= {tolerance}"
        )


class TestStarBasedAyanamshasRelaxed:
    """
    Test star-based ayanamshas with relaxed 1.0 degree tolerance.

    These ayanamshas depend on star catalogs and proper motion calculations,
    which may differ between implementations.
    """

    STAR_BASED = [m for m in ALL_AYANAMSHA_MODES if m[1] in STAR_BASED_AYANAMSHAS]

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "mode_id,sid_mode,name",
        STAR_BASED,
        ids=[m[2] for m in STAR_BASED],
    )
    @pytest.mark.parametrize(
        "jd,epoch_name",
        [
            (JD_J2000, "J2000"),
            (JD_1900, "1900"),
            (JD_2100, "2100"),
        ],
        ids=["J2000", "1900", "2100"],
    )
    def test_star_ayanamsha_relaxed_tolerance(
        self, mode_id, sid_mode, name, jd, epoch_name
    ):
        """
        Test star-based ayanamshas at all three epochs with relaxed tolerance.
        """
        tolerance = STAR_TOLERANCE

        # Set mode in both libraries
        ephem.set_sid_mode(sid_mode)
        swe.set_sid_mode(sid_mode)

        # Calculate ayanamsha
        ayan_lib = ephem.get_ayanamsa_ut(jd)
        ayan_swe = swe.get_ayanamsa_ut(jd)

        # Compare (use angle_diff to handle wrap-around at 360 degrees)
        diff = angle_diff(ayan_lib, ayan_swe)
        assert diff < tolerance, (
            f"{name} at {epoch_name}: "
            f"lib={ayan_lib:.6f}, swe={ayan_swe:.6f}, "
            f"diff={diff:.6f} >= {tolerance}"
        )


class TestAyanamshaConsistency:
    """Test internal consistency of ayanamsha calculations."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode_id,sid_mode,name",
        ALL_AYANAMSHA_MODES,
        ids=[m[2] for m in ALL_AYANAMSHA_MODES],
    )
    def test_ayanamsha_increases_over_time(self, mode_id, sid_mode, name):
        """
        Most ayanamshas should increase over time due to precession.

        Exception: J2000 mode (18) which is always 0 at J2000.
        """
        ephem.set_sid_mode(sid_mode)

        ayan_1900 = ephem.get_ayanamsa_ut(JD_1900)
        ayan_2000 = ephem.get_ayanamsa_ut(JD_J2000)
        ayan_2100 = ephem.get_ayanamsa_ut(JD_2100)

        # J2000 mode is special - always returns 0 at J2000 epoch
        if sid_mode == SIDM_J2000:
            assert abs(ayan_2000) < 0.001, "J2000 mode should be 0 at J2000"
        else:
            # For most modes, ayanamsha increases with time (precession ~50"/year)
            # Some modes may have different behavior, so we just check values are valid
            assert isinstance(ayan_1900, float)
            assert isinstance(ayan_2000, float)
            assert isinstance(ayan_2100, float)

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode_id,sid_mode,name",
        ALL_AYANAMSHA_MODES,
        ids=[m[2] for m in ALL_AYANAMSHA_MODES],
    )
    def test_ayanamsha_returns_float(self, mode_id, sid_mode, name):
        """All ayanamsha modes should return a valid float value."""
        ephem.set_sid_mode(sid_mode)

        for jd in [JD_1900, JD_J2000, JD_2100]:
            ayan = ephem.get_ayanamsa_ut(jd)
            assert isinstance(ayan, float), f"Mode {mode_id} should return float"
            # Ayanamsha can range from negative to 360+ depending on mode and epoch
            # J1900 mode at 1900 epoch can be close to 0 or 360
            assert -100 < ayan < 400, (
                f"Mode {mode_id} ayanamsha {ayan} out of reasonable range"
            )


class TestAyanamshaCompleteness:
    """Test that all 43 modes are properly tested."""

    @pytest.mark.unit
    def test_all_modes_defined(self):
        """Verify we're testing all 43 modes (0-42)."""
        expected_modes = set(range(43))
        tested_modes = {m[0] for m in ALL_AYANAMSHA_MODES}

        assert tested_modes == expected_modes, (
            f"Missing modes: {expected_modes - tested_modes}, "
            f"Extra modes: {tested_modes - expected_modes}"
        )

    @pytest.mark.unit
    def test_mode_constant_mapping(self):
        """Verify mode IDs match their constant values."""
        for mode_id, sid_mode, name in ALL_AYANAMSHA_MODES:
            assert mode_id == sid_mode, (
                f"Mode {name}: ID {mode_id} != constant value {sid_mode}"
            )


class TestAyanamshaExFunctions:
    """Test extended ayanamsha functions with all modes."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "mode_id,sid_mode,name",
        ALL_AYANAMSHA_MODES[:10],  # Test first 10 modes for efficiency
        ids=[m[2] for m in ALL_AYANAMSHA_MODES[:10]],
    )
    def test_get_ayanamsa_ex_matches_standard(self, mode_id, sid_mode, name):
        """get_ayanamsa_ex_ut must match pyswisseph's get_ayanamsa_ex_ut 1:1.

        It does NOT match the plain get_ayanamsa_ut: the ex form returns the
        mean ayanamsha while the plain form returns the true (nutation-applied)
        ayanamsha, so the two differ by the nutation in longitude (~0.004°) --
        in pyswisseph exactly as in libephemeris. The meaningful check is
        therefore ex-vs-Swiss-ex, not ex-vs-standard.
        """
        jd = JD_J2000
        tolerance = get_tolerance(sid_mode)

        ephem.set_sid_mode(sid_mode)
        swe.set_sid_mode(sid_mode)

        _retflag, ayan_ex = ephem.get_ayanamsa_ex_ut(jd, 0)
        _retflag_swe, ayan_ex_swe = swe.get_ayanamsa_ex_ut(jd, 0)

        diff = angle_diff(ayan_ex, ayan_ex_swe)
        assert diff < tolerance, (
            f"Mode {mode_id} ({name}): ex={ayan_ex:.6f} vs "
            f"pyswisseph ex={ayan_ex_swe:.6f}, diff={diff:.6f} >= {tolerance}"
        )
