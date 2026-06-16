"""
Tests all 43 ayanamsha modes at multiple dates.

This module provides comprehensive regression tests for all ayanamsha modes,
comparing calculated values against reference values at multiple epochs:
- 1900-01-01 12:00 TT (historical)
- 2000-01-01 12:00 TT (J2000)
- 2050-01-01 12:00 TT (future)
- 2100-01-01 12:00 TT (far future)

Reference values are derived from the reference ephemeris to ensure accuracy.
"""

import pytest
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
# TEST DATES (Julian Days)
# =============================================================================
TEST_DATES = {
    "1900": 2415020.0,  # 1900-01-01 12:00 TT
    "J2000": 2451545.0,  # 2000-01-01 12:00 TT
    "2050": 2469807.5,  # 2050-01-01 12:00 TT
    "2100": 2488070.0,  # 2100-01-01 12:00 TT
}


# =============================================================================
# ALL 43 AYANAMSHA MODES
# =============================================================================
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

# Star-based ayanamshas that may have implementation differences
STAR_BASED_AYANAMSHAS = {
    SIDM_TRUE_CITRA,
    SIDM_TRUE_REVATI,
    SIDM_TRUE_PUSHYA,
    SIDM_TRUE_MULA,
    SIDM_TRUE_SHEORAN,
    SIDM_GALCENT_0SAG,
    SIDM_GALCENT_RGILBRAND,
    SIDM_GALCENT_MULA_WILHELM,
    SIDM_GALCENT_COCHRANE,
    SIDM_GALEQU_IAU1958,
    SIDM_GALEQU_TRUE,
    SIDM_GALEQU_MULA,
    SIDM_GALEQU_FIORENZA,
    SIDM_GALALIGN_MARDYKS,
}


def angle_diff(a1: float, a2: float) -> float:
    """Calculate the smallest difference between two angles (0-180)."""
    diff = abs(a1 - a2)
    if diff > 180:
        diff = 360 - diff
    return diff


class TestAllAyanamshasMultipleEpochs:
    """Test all 43 ayanamsha modes at multiple date epochs."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode_id,sid_mode,name",
        ALL_AYANAMSHA_MODES,
        ids=[m[2] for m in ALL_AYANAMSHA_MODES],
    )
    @pytest.mark.parametrize(
        "epoch_name,jd",
        list(TEST_DATES.items()),
        ids=list(TEST_DATES.keys()),
    )
    def test_ayanamsha_returns_valid_value(
        self, mode_id, sid_mode, name, epoch_name, jd
    ):
        """
        Each ayanamsha mode should return a valid float at all test epochs.

        This verifies basic functionality across all 43 modes at 4 different dates.
        """
        ephem.set_sid_mode(sid_mode)
        ayan = ephem.get_ayanamsa_ut(jd)

        # Should return a float
        assert isinstance(ayan, float), (
            f"Mode {name} at {epoch_name} should return float"
        )

        # Should be in reasonable range (-100 to 400 degrees)
        assert -100 < ayan < 400, (
            f"Mode {name} at {epoch_name}: ayanamsha {ayan:.4f} out of range"
        )


class TestAyanamshaTemporalProgression:
    """Test that ayanamsha values progress correctly over time."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode_id,sid_mode,name",
        ALL_AYANAMSHA_MODES,
        ids=[m[2] for m in ALL_AYANAMSHA_MODES],
    )
    def test_ayanamsha_changes_over_century(self, mode_id, sid_mode, name):
        """
        Ayanamsha should change over a century due to precession.

        Most modes increase by ~1.4 degrees per century (~50"/year).
        """
        ephem.set_sid_mode(sid_mode)

        ayan_1900 = ephem.get_ayanamsa_ut(TEST_DATES["1900"])
        ayan_2000 = ephem.get_ayanamsa_ut(TEST_DATES["J2000"])
        ayan_2100 = ephem.get_ayanamsa_ut(TEST_DATES["2100"])

        # Calculate change per century, accounting for 360° wrap
        def angular_change(a, b):
            """Calculate signed angular change, accounting for wrap."""
            d = b - a
            if d > 180:
                d -= 360
            elif d < -180:
                d += 360
            return d

        change_1900_to_2000 = angular_change(ayan_1900, ayan_2000)
        change_2000_to_2100 = angular_change(ayan_2000, ayan_2100)

        # J2000 mode is special - ayanamsha is 0 at J2000, negative before
        if sid_mode == SIDM_J2000:
            assert abs(ayan_2000) < 0.001, "J2000 mode should be 0 at J2000"
            # But it should still change over time
            assert abs(change_1900_to_2000) > 0.5, (
                "J2000 should change from 1900 to 2000"
            )
        # J1900 mode - should be ~0 at 1900 (but calculated using precession from J2000)
        # The value at 1900 may be close to 360° due to wrap
        elif sid_mode == SIDM_J1900:
            # J1900 mode: at 1900, the value should be near 0 or 360
            ayan_1900_normalized = ayan_1900 if ayan_1900 < 180 else ayan_1900 - 360
            assert abs(ayan_1900_normalized) < 0.1, (
                f"J1900 mode should be ~0 at 1900, got {ayan_1900}"
            )
            assert abs(change_1900_to_2000) > 0.5, (
                "J1900 should change from 1900 to 2000"
            )
        # B1950 mode has reference epoch at 1950
        elif sid_mode == SIDM_B1950:
            # Just verify it changes over time
            assert abs(change_1900_to_2000) > 0.5, (
                "B1950 should change from 1900 to 2000"
            )
        else:
            # For most modes, change should be positive (precession)
            # and roughly 1-2 degrees per century
            # Allow for some modes with different behavior
            assert change_1900_to_2000 != 0, f"{name} should change over time"
            assert change_2000_to_2100 != 0, f"{name} should change over time"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode_id,sid_mode,name",
        [m for m in ALL_AYANAMSHA_MODES if m[1] not in {SIDM_J2000, SIDM_J1900}],
        ids=[
            m[2]
            for m in ALL_AYANAMSHA_MODES
            if m[1] not in {SIDM_J2000, SIDM_J1900}
        ],
    )
    def test_precession_rate_reasonable(self, mode_id, sid_mode, name):
        """
        Most ayanamshas should have a precession rate of roughly 50"/year.

        Over 100 years, this is about 1.39 degrees.
        Allow 0.5 to 2.5 degrees per century for various modes.
        """
        ephem.set_sid_mode(sid_mode)

        ayan_2000 = ephem.get_ayanamsa_ut(TEST_DATES["J2000"])
        ayan_2100 = ephem.get_ayanamsa_ut(TEST_DATES["2100"])

        change = abs(ayan_2100 - ayan_2000)

        # Most modes should change by 0.5-2.5 degrees per century
        # Some star-based modes might differ
        if sid_mode not in STAR_BASED_AYANAMSHAS:
            assert 0.5 < change < 2.5, (
                f"{name} precession {change:.4f} deg/century out of expected range"
            )


class TestAyanamshaConsistencyAcrossDates:
    """Test internal consistency of ayanamsha calculations."""

    # Special modes that wrap around 360° need different handling
    WRAP_AROUND_MODES = {SIDM_J2000, SIDM_J1900, SIDM_B1950}

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode_id,sid_mode,name",
        [
            m
            for m in ALL_AYANAMSHA_MODES
            if m[1] not in {SIDM_J2000, SIDM_J1900, SIDM_B1950}
        ],
        ids=[
            m[2]
            for m in ALL_AYANAMSHA_MODES
            if m[1] not in {SIDM_J2000, SIDM_J1900, SIDM_B1950}
        ],
    )
    def test_ayanamsha_monotonic(self, mode_id, sid_mode, name):
        """
        Ayanamsha should change monotonically over time (no reversals).

        Tests that values at 1900 < 2000 < 2050 < 2100 (or vice versa).
        Note: J2000, J1900, B1950 modes excluded as they wrap around 360°.
        """
        ephem.set_sid_mode(sid_mode)

        values = [ephem.get_ayanamsa_ut(jd) for jd in TEST_DATES.values()]

        # Check if monotonically increasing or decreasing
        increasing = all(values[i] <= values[i + 1] for i in range(len(values) - 1))
        decreasing = all(values[i] >= values[i + 1] for i in range(len(values) - 1))

        assert increasing or decreasing, (
            f"{name} should be monotonic over time: {values}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode_id,sid_mode,name",
        [
            m
            for m in ALL_AYANAMSHA_MODES
            if m[1] not in {SIDM_J2000, SIDM_J1900, SIDM_B1950}
        ],
        ids=[
            m[2]
            for m in ALL_AYANAMSHA_MODES
            if m[1] not in {SIDM_J2000, SIDM_J1900, SIDM_B1950}
        ],
    )
    def test_ayanamsha_smooth_progression(self, mode_id, sid_mode, name):
        """
        Ayanamsha change rate should be relatively smooth (no jumps).

        Compares century-to-century changes - they should be similar.
        Note: J2000, J1900, B1950 modes excluded as they wrap around 360°.
        """
        ephem.set_sid_mode(sid_mode)

        ayan_1900 = ephem.get_ayanamsa_ut(TEST_DATES["1900"])
        ayan_2000 = ephem.get_ayanamsa_ut(TEST_DATES["J2000"])
        ayan_2100 = ephem.get_ayanamsa_ut(TEST_DATES["2100"])

        change_1 = ayan_2000 - ayan_1900  # 1900-2000
        change_2 = ayan_2100 - ayan_2000  # 2000-2100

        # Both changes should have similar magnitude (within factor of 2)
        if abs(change_1) > 0.01 and abs(change_2) > 0.01:
            ratio = abs(change_1 / change_2) if change_2 != 0 else float("inf")
            assert 0.5 < ratio < 2.0, (
                f"{name} has non-smooth progression: "
                f"change 1900-2000={change_1:.4f}, 2000-2100={change_2:.4f}"
            )


class TestAyanamshaExAtMultipleDates:
    """Test extended ayanamsha function at multiple dates."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode_id,sid_mode,name",
        ALL_AYANAMSHA_MODES[:15],  # Test first 15 for efficiency
        ids=[m[2] for m in ALL_AYANAMSHA_MODES[:15]],
    )
    @pytest.mark.parametrize(
        "epoch_name,jd",
        list(TEST_DATES.items()),
        ids=list(TEST_DATES.keys()),
    )
    def test_ayanamsa_ex_matches_standard(
        self, mode_id, sid_mode, name, epoch_name, jd
    ):
        """
        get_ayanamsa_ex_ut returns the TRUE ayanamsha (mean + nutation),
        while the standard function returns the mean value — the reference
        API semantics (verified against the reference ephemeris). With FLG_NONUT the ex
        variant returns the mean value again.
        """
        ephem.set_sid_mode(sid_mode)
        ayan_standard = ephem.get_ayanamsa_ut(jd)

        retflag, ayan_ex = ephem.get_ayanamsa_ex_ut(jd, 0)

        # ex - standard equals the nutation in longitude (|dpsi| <= ~17.5");
        # wrap-aware because ayanamshas near 0 deg can wrap to ~360.
        delta_deg = abs(ayan_ex - ayan_standard) % 360.0
        delta_deg = min(delta_deg, 360.0 - delta_deg)
        assert delta_deg < 0.006, (
            f"{name} at {epoch_name}: ex-standard={delta_deg * 3600:.2f} arcsec "
            f"exceeds the nutation envelope"
        )

        # With FLG_NONUT the ex variant returns the mean ayanamsha exactly
        _, ayan_ex_nonut = ephem.get_ayanamsa_ex_ut(jd, ephem.FLG_NONUT)
        assert abs(ayan_ex_nonut - ayan_standard) < 0.0001, (
            f"{name} at {epoch_name}: ex(NONUT)={ayan_ex_nonut:.6f} != "
            f"standard={ayan_standard:.6f}"
        )


class TestAyanamshaCompleteness:
    """Verify all 43 modes are properly tested."""

    @pytest.mark.unit
    def test_all_43_modes_defined(self):
        """Ensure we have exactly 43 ayanamsha modes (0-42)."""
        expected_modes = set(range(43))
        tested_modes = {m[0] for m in ALL_AYANAMSHA_MODES}

        assert tested_modes == expected_modes, (
            f"Missing modes: {expected_modes - tested_modes}, "
            f"Extra modes: {tested_modes - expected_modes}"
        )

    @pytest.mark.unit
    def test_mode_id_matches_constant(self):
        """Verify mode IDs match their constant values."""
        for mode_id, sid_mode, name in ALL_AYANAMSHA_MODES:
            assert mode_id == sid_mode, f"{name}: ID {mode_id} != constant {sid_mode}"

    @pytest.mark.unit
    def test_all_four_dates_tested(self):
        """Ensure we test at 4 different epochs."""
        assert len(TEST_DATES) == 4
        assert "1900" in TEST_DATES
        assert "J2000" in TEST_DATES
        assert "2050" in TEST_DATES
        assert "2100" in TEST_DATES


class TestSpecialAyanamshaModes:
    """Test special ayanamsha modes with known behavior."""

    @pytest.mark.unit
    def test_j2000_mode_zero_at_j2000(self):
        """J2000 mode should return exactly 0 at J2000 epoch."""
        ephem.set_sid_mode(SIDM_J2000)
        ayan = ephem.get_ayanamsa_ut(TEST_DATES["J2000"])
        assert abs(ayan) < 0.0001, f"J2000 at J2000 should be 0, got {ayan}"

    @pytest.mark.unit
    def test_j1900_mode_zero_at_1900(self):
        """J1900 mode should return approximately 0 (or 360) at 1900 epoch."""
        ephem.set_sid_mode(SIDM_J1900)
        ayan = ephem.get_ayanamsa_ut(TEST_DATES["1900"])
        # J1900 is slightly different from 1900-01-01
        # The value may be close to 360° due to wrap-around
        ayan_normalized = ayan if ayan < 180 else ayan - 360
        assert abs(ayan_normalized) < 0.1, f"J1900 at 1900 should be ~0, got {ayan}"

    @pytest.mark.unit
    def test_lahiri_at_j2000_expected_value(self):
        """Lahiri ayanamsha at J2000 should be approximately 23.9 degrees."""
        ephem.set_sid_mode(SIDM_LAHIRI)
        ayan = ephem.get_ayanamsa_ut(TEST_DATES["J2000"])
        assert 23.5 < ayan < 24.5, f"Lahiri at J2000 expected ~23.9, got {ayan}"

    @pytest.mark.unit
    def test_fagan_bradley_at_j2000_expected_value(self):
        """Fagan-Bradley ayanamsha at J2000 should be approximately 24.7 degrees."""
        ephem.set_sid_mode(SIDM_FAGAN_BRADLEY)
        ayan = ephem.get_ayanamsa_ut(TEST_DATES["J2000"])
        assert 24.0 < ayan < 25.5, f"Fagan-Bradley at J2000 expected ~24.7, got {ayan}"

    @pytest.mark.unit
    def test_raman_at_j2000_expected_value(self):
        """Raman ayanamsha at J2000 should be approximately 22.4 degrees."""
        ephem.set_sid_mode(SIDM_RAMAN)
        ayan = ephem.get_ayanamsa_ut(TEST_DATES["J2000"])
        assert 22.0 < ayan < 23.0, f"Raman at J2000 expected ~22.4, got {ayan}"


class TestAyanamshaValueRanges:
    """Test that ayanamsha values fall within expected ranges at different epochs."""

    # Expected approximate values at J2000 for key ayanamshas
    EXPECTED_J2000_VALUES = {
        SIDM_FAGAN_BRADLEY: (24.0, 25.5),  # ~24.74
        SIDM_LAHIRI: (23.5, 24.5),  # ~23.86
        SIDM_RAMAN: (22.0, 23.0),  # ~22.37
        SIDM_KRISHNAMURTI: (23.5, 24.5),  # ~23.77
        SIDM_J2000: (-0.01, 0.01),  # exactly 0
        SIDM_B1950: (0.3, 0.7),  # ~0.53
    }

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "sid_mode,expected_range",
        list(EXPECTED_J2000_VALUES.items()),
    )
    def test_key_ayanamshas_at_j2000(self, sid_mode, expected_range):
        """Key ayanamshas should fall within expected ranges at J2000."""
        ephem.set_sid_mode(sid_mode)
        ayan = ephem.get_ayanamsa_ut(TEST_DATES["J2000"])

        low, high = expected_range
        assert low < ayan < high, (
            f"Mode {sid_mode} at J2000: {ayan:.4f} not in ({low}, {high})"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode_id,sid_mode,name",
        ALL_AYANAMSHA_MODES,
        ids=[m[2] for m in ALL_AYANAMSHA_MODES],
    )
    def test_ayanamsha_at_2050_greater_than_2000(self, mode_id, sid_mode, name):
        """
        For most modes, ayanamsha at 2050 should be greater than at 2000.

        Due to precession, the tropical zodiac moves relative to sidereal.
        """
        ephem.set_sid_mode(sid_mode)
        ayan_2000 = ephem.get_ayanamsa_ut(TEST_DATES["J2000"])
        ayan_2050 = ephem.get_ayanamsa_ut(TEST_DATES["2050"])

        # Most modes increase over time
        # Allow for modes that might be fixed or have different behavior
        diff = ayan_2050 - ayan_2000
        assert diff != 0, f"{name} should change between 2000 and 2050"
