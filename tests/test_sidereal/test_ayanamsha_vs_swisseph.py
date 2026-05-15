"""
Ayanamsha comparison tests vs pyswisseph.

Compares all 47 ayanamsha modes (0-46) at J2000 against pyswisseph,
verifies get_ayanamsa_name returns non-empty strings, and checks
sidereal positions for selected modes x bodies x dates.
"""

from __future__ import annotations

import math

import pytest
import swisseph as swe_ref

import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    FLG_SIDEREAL,
    FLG_SPEED,
    SIDM_FAGAN_BRADLEY,
    SIDM_LAHIRI,
    SIDM_RAMAN,
    SIDM_KRISHNAMURTI,
    SIDM_YUKTESHWAR,
    SIDM_TRUE_CITRA,
    SIDM_TRUE_REVATI,
    SIDM_SURYASIDDHANTA,
    SIDM_J2000,
    SIDM_HIPPARCHOS,
    SIDM_LAHIRI_1940,
    SIDM_LAHIRI_VP285,
    SIDM_KRISHNAMURTI_VP291,
    SIDM_LAHIRI_ICRC,
)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# All valid ayanamsha modes (0-46); mode 47+ returns empty name in pyswisseph
ALL_MODES = list(range(47))

J2000 = 2451545.0

# "True star" modes compute ayanamsha from actual star positions, which differ
# by ~10 arcseconds between DE440 (libephemeris) and the Swiss Ephemeris
# analytical fit.  They also tend to return numpy float64 instead of native
# Python float because of the stellar position pipeline.
TRUE_STAR_MODES = {17, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 39, 40, 42}

# Modes that use galactic-center based offsets and may coincide at J2000
# (modes 31/GALEQU_IAU1958 and 34/GALALIGN_MARDYKS share the same offset)
GALACTIC_COINCIDENT_PAIRS = {(31, 34)}

# Tolerance for ayanamsha comparison at J2000.
# Standard modes: implementations share the same polynomial coefficients
AYAN_TOL_DEG = 0.001  # 0.001 degrees = 3.6 arcseconds
# True star modes: differ because of different star catalogs / ephemerides
AYAN_TRUE_TOL_DEG = 0.005  # 0.005 degrees = 18 arcseconds

# Modes where libephemeris and pyswisseph are known to agree well
SELECTED_MODES = [
    (SIDM_FAGAN_BRADLEY, "Fagan-Bradley"),
    (SIDM_LAHIRI, "Lahiri"),
    (SIDM_RAMAN, "Raman"),
    (SIDM_KRISHNAMURTI, "Krishnamurti"),
    (SIDM_YUKTESHWAR, "Yukteshwar"),
]

POSITION_TEST_BODIES = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MARS, "Mars"),
]

POSITION_TEST_DATES = [
    2415020.0,  # 1900-01-01
    2451545.0,  # J2000
    2460310.5,  # 2024-01-18
]

# Tolerance for sidereal position comparison: 2 arcseconds
POS_TOL_DEG = 2.0 / 3600.0


def _angle_diff(a: float, b: float) -> float:
    """Signed angular difference, handling wrap-around at 360."""
    d = float(a) - float(b)
    while d > 180:
        d -= 360
    while d < -180:
        d += 360
    return d


# ---------------------------------------------------------------------------
# Tests: Ayanamsha values at J2000
# ---------------------------------------------------------------------------


class TestAyanamshaAtJ2000:
    """Compare ayanamsha values at J2000 against pyswisseph."""

    @pytest.mark.unit
    @pytest.mark.parametrize("mode", ALL_MODES)
    def test_ayanamsha_matches_swisseph(self, mode: int):
        """Ayanamsha at J2000 matches pyswisseph within tolerance.

        Standard modes use 0.001 degree; true-star modes use 0.005 degree
        to account for different star position catalogs.
        """
        swe.set_sid_mode(mode)
        swe_ref.set_sid_mode(mode)

        lib_ayan = float(swe.get_ayanamsa_ut(J2000))
        ref_ayan = float(swe_ref.get_ayanamsa_ut(J2000))

        tol = AYAN_TRUE_TOL_DEG if mode in TRUE_STAR_MODES else AYAN_TOL_DEG
        diff = abs(lib_ayan - ref_ayan)
        assert diff < tol, (
            f"Mode {mode}: lib={lib_ayan:.8f}, ref={ref_ayan:.8f}, "
            f'diff={diff:.6f} deg ({diff * 3600:.2f}"), tol={tol}'
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("mode", ALL_MODES)
    def test_ayanamsha_finite(self, mode: int):
        """Each ayanamsha mode returns a finite value at J2000."""
        swe.set_sid_mode(mode)
        ayan = swe.get_ayanamsa_ut(J2000)
        assert math.isfinite(float(ayan)), f"Mode {mode}: ayanamsha={ayan} not finite"

    @pytest.mark.unit
    @pytest.mark.parametrize("mode", ALL_MODES)
    def test_ayanamsha_returns_numeric(self, mode: int):
        """get_ayanamsa_ut returns a Python float or numpy float64.

        Some true-star modes return numpy float64 from the stellar position
        pipeline; both are acceptable numeric types.
        """
        swe.set_sid_mode(mode)
        ayan = swe.get_ayanamsa_ut(J2000)
        assert isinstance(ayan, (int, float)) or hasattr(ayan, "__float__"), (
            f"Mode {mode}: type={type(ayan)} is not numeric"
        )
        # Verify it can be converted to float without error
        assert math.isfinite(float(ayan))


# ---------------------------------------------------------------------------
# Tests: Ayanamsha names
# ---------------------------------------------------------------------------


class TestAyanamshaNames:
    """Verify get_ayanamsa_name returns non-empty strings."""

    @pytest.mark.unit
    @pytest.mark.parametrize("mode", ALL_MODES)
    def test_name_non_empty(self, mode: int):
        """Each valid mode returns a non-empty name string."""
        name = swe.get_ayanamsa_name(mode)
        assert isinstance(name, str), f"Mode {mode}: type={type(name)}"
        assert len(name) > 0, f"Mode {mode}: name is empty"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode,expected_fragment",
        [
            (SIDM_FAGAN_BRADLEY, "Fagan"),
            (SIDM_LAHIRI, "Lahiri"),
            (SIDM_RAMAN, "Raman"),
            (SIDM_KRISHNAMURTI, "Krishnamurti"),
            (SIDM_YUKTESHWAR, "Yukteshwar"),
            (SIDM_TRUE_CITRA, "Citra"),
            (SIDM_TRUE_REVATI, "Revati"),
            (SIDM_SURYASIDDHANTA, "Surya"),
            (SIDM_J2000, "J2000"),
            (SIDM_HIPPARCHOS, "Hipparchos"),
        ],
    )
    def test_name_contains_expected(self, mode: int, expected_fragment: str):
        """Known mode names contain expected keyword."""
        name = swe.get_ayanamsa_name(mode)
        assert expected_fragment.lower() in name.lower(), (
            f"Mode {mode}: name='{name}' does not contain '{expected_fragment}'"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("mode", ALL_MODES)
    def test_name_matches_swisseph(self, mode: int):
        """Mode names match pyswisseph (case-insensitive substring)."""
        lib_name = swe.get_ayanamsa_name(mode)
        ref_name = swe_ref.get_ayanamsa_name(mode)
        # Some names may differ slightly in formatting, so we compare
        # lowercased first few characters
        assert lib_name.lower()[:8] == ref_name.lower()[:8], (
            f"Mode {mode}: lib='{lib_name}', ref='{ref_name}'"
        )


# ---------------------------------------------------------------------------
# Tests: Ayanamsha variation over time
# ---------------------------------------------------------------------------


class TestAyanamshaTimeVariation:
    """Verify ayanamsha changes over time (precession)."""

    @pytest.mark.unit
    @pytest.mark.parametrize("mode", [0, 1, 3, 5, 7, 27, 28])
    def test_ayanamsha_increases_with_time(self, mode: int):
        """Most ayanamsha modes increase over time due to precession.

        For most modes, ayanamsha in 2000 > ayanamsha in 1900.
        (Exception: mode 40/Cochrane wraps around, but we test safe modes.)
        """
        swe.set_sid_mode(mode)
        ayan_1900 = float(swe.get_ayanamsa_ut(2415020.0))
        ayan_2000 = float(swe.get_ayanamsa_ut(2451545.0))
        ayan_2100 = float(swe.get_ayanamsa_ut(2488070.0))

        assert ayan_2000 > ayan_1900, (
            f"Mode {mode}: ayanamsha should increase 1900->2000 "
            f"({ayan_1900:.4f} -> {ayan_2000:.4f})"
        )
        assert ayan_2100 > ayan_2000, (
            f"Mode {mode}: ayanamsha should increase 2000->2100 "
            f"({ayan_2000:.4f} -> {ayan_2100:.4f})"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("mode", [0, 1, 3, 5, 7])
    def test_precession_rate_approximately_50_arcsec_per_year(self, mode: int):
        """Annual precession rate is approximately 50\"/year for standard modes."""
        swe.set_sid_mode(mode)
        ayan_2000 = float(swe.get_ayanamsa_ut(2451545.0))
        ayan_2100 = float(swe.get_ayanamsa_ut(2488070.0))

        rate_arcsec = (ayan_2100 - ayan_2000) * 3600 / 100.0  # arcsec/year
        assert 49 < rate_arcsec < 51, (
            f'Mode {mode}: precession rate={rate_arcsec:.2f}"/yr, expected ~50'
        )


# ---------------------------------------------------------------------------
# Tests: Sidereal positions for selected modes
# ---------------------------------------------------------------------------


class TestSiderealPositions:
    """Compare sidereal positions against pyswisseph for selected modes."""

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,mode_name", SELECTED_MODES)
    @pytest.mark.parametrize("body_id,body_name", POSITION_TEST_BODIES)
    @pytest.mark.parametrize("jd", POSITION_TEST_DATES)
    def test_sidereal_position_matches(
        self,
        mode: int,
        mode_name: str,
        body_id: int,
        body_name: str,
        jd: float,
    ):
        """Sidereal longitude matches pyswisseph within 2 arcseconds."""
        swe.set_sid_mode(mode)
        swe_ref.set_sid_mode(mode)

        flags = FLG_SIDEREAL | FLG_SPEED
        lib_vals, _ = swe.calc_ut(jd, body_id, flags)
        ref_vals, _ = swe_ref.calc_ut(jd, body_id, flags)

        dlon = abs(_angle_diff(lib_vals[0], ref_vals[0]))
        assert dlon < POS_TOL_DEG, (
            f'{mode_name} {body_name} jd={jd}: lon diff={dlon * 3600:.4f}" '
            f"(lib={lib_vals[0]:.8f}, ref={ref_vals[0]:.8f})"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,mode_name", SELECTED_MODES)
    @pytest.mark.parametrize("body_id,body_name", POSITION_TEST_BODIES)
    def test_sidereal_longitude_in_range(
        self, mode: int, mode_name: str, body_id: int, body_name: str
    ):
        """Sidereal longitude is in [0, 360) range."""
        swe.set_sid_mode(mode)
        vals, _ = swe.calc_ut(J2000, body_id, FLG_SIDEREAL)
        assert 0 <= vals[0] < 360, (
            f"{mode_name} {body_name}: lon={vals[0]} out of range"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,mode_name", SELECTED_MODES)
    def test_sidereal_vs_tropical_offset(self, mode: int, mode_name: str):
        """Sidereal longitude = tropical longitude - ayanamsha (mod 360)."""
        swe.set_sid_mode(mode)
        jd = J2000

        trop_vals, _ = swe.calc_ut(jd, SUN, FLG_SPEED)
        sid_vals, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL | FLG_SPEED)
        ayan = float(swe.get_ayanamsa_ut(jd))

        expected_sid = (trop_vals[0] - ayan) % 360
        actual_sid = sid_vals[0]

        diff = abs(_angle_diff(actual_sid, expected_sid))
        # Allow 0.01 degree because ayanamsha is applied internally
        # with slightly different precision
        assert diff < 0.01, (
            f"{mode_name}: sid={actual_sid:.8f}, "
            f"trop-ayan={expected_sid:.8f}, diff={diff:.6f}"
        )


# ---------------------------------------------------------------------------
# Tests: Extended mode set
# ---------------------------------------------------------------------------


class TestExtendedModes:
    """Test newer/less common ayanamsha modes (43-46)."""

    EXTENDED_MODES = [
        (SIDM_LAHIRI_1940, "Lahiri 1940"),
        (SIDM_LAHIRI_VP285, "Lahiri VP285"),
        (SIDM_KRISHNAMURTI_VP291, "Krishnamurti VP291"),
        (SIDM_LAHIRI_ICRC, "Lahiri ICRC"),
    ]

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,name", EXTENDED_MODES)
    def test_extended_mode_finite(self, mode: int, name: str):
        """Extended mode returns a finite ayanamsha."""
        swe.set_sid_mode(mode)
        ayan = swe.get_ayanamsa_ut(J2000)
        assert math.isfinite(float(ayan)), f"{name}: ayan={ayan}"

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,name", EXTENDED_MODES)
    def test_extended_mode_matches_swisseph(self, mode: int, name: str):
        """Extended mode matches pyswisseph within tolerance."""
        swe.set_sid_mode(mode)
        swe_ref.set_sid_mode(mode)

        lib_ayan = float(swe.get_ayanamsa_ut(J2000))
        ref_ayan = float(swe_ref.get_ayanamsa_ut(J2000))

        diff = abs(lib_ayan - ref_ayan)
        assert diff < AYAN_TOL_DEG, (
            f"{name}: lib={lib_ayan:.8f}, ref={ref_ayan:.8f}, diff={diff:.6f}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,name", EXTENDED_MODES)
    def test_extended_mode_calc_works(self, mode: int, name: str):
        """Extended mode works with calc_ut + FLG_SIDEREAL."""
        swe.set_sid_mode(mode)
        vals, _ = swe.calc_ut(J2000, SUN, FLG_SIDEREAL)
        assert 0 <= vals[0] < 360, f"{name}: lon={vals[0]}"


class TestDistinctModes:
    """Verify that different modes produce different ayanamsha values."""

    @pytest.mark.unit
    def test_lahiri_vs_fagan_bradley(self):
        """Lahiri and Fagan-Bradley produce different ayanamshas."""
        swe.set_sid_mode(SIDM_LAHIRI)
        ayan_lahiri = float(swe.get_ayanamsa_ut(J2000))

        swe.set_sid_mode(SIDM_FAGAN_BRADLEY)
        ayan_fb = float(swe.get_ayanamsa_ut(J2000))

        assert abs(ayan_lahiri - ayan_fb) > 0.1, (
            f"Lahiri={ayan_lahiri:.6f}, Fagan-Bradley={ayan_fb:.6f} "
            "should differ by more than 0.1 degree"
        )

    @pytest.mark.unit
    def test_j2000_mode_is_near_zero(self):
        """J2000 mode should have ayanamsha very close to 0.0."""
        swe.set_sid_mode(SIDM_J2000)
        ayan = float(swe.get_ayanamsa_ut(J2000))
        # Allow 1e-6 degree (about 0.004 arcsecond) for numerical noise
        assert abs(ayan) < 1e-6, f"J2000 mode ayanamsha={ayan}, expected ~0.0"

    @pytest.mark.unit
    def test_all_standard_modes_mostly_distinct(self):
        """Standard modes (0-42) produce mostly distinct values at J2000.

        Known coincidences: modes 31 (GALEQU_IAU1958) and 34
        (GALALIGN_MARDYKS) share the same galactic-plane offset.
        """
        values = {}
        for mode in range(43):
            swe.set_sid_mode(mode)
            ayan = float(swe.get_ayanamsa_ut(J2000))
            values[mode] = ayan

        # Check that most pairs differ by at least 0.001 degree
        modes = sorted(values.keys())
        duplicate_count = 0
        for i in range(len(modes)):
            for j in range(i + 1, len(modes)):
                m1, m2 = modes[i], modes[j]
                diff = abs(values[m1] - values[m2])
                if diff < 1e-6:
                    pair = (min(m1, m2), max(m1, m2))
                    if pair not in GALACTIC_COINCIDENT_PAIRS:
                        duplicate_count += 1

        # At most 0 unexpected duplicates
        assert duplicate_count == 0, (
            f"Found {duplicate_count} unexpected duplicate pairs"
        )
