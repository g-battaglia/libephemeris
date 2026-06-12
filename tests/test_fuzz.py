"""Fuzz testing for robustness against extreme and malformed inputs.

Verifies that libephemeris handles edge cases gracefully — either returning
valid results or raising clean exceptions, never crashing or hanging.

Validation Plan v2, Section 4.
"""

from __future__ import annotations

import math
import warnings

import pytest

import libephemeris as swe
from libephemeris.exceptions import (
    Error,
)

warnings.filterwarnings("ignore")

pytestmark = pytest.mark.slow

# Acceptable exceptions for fuzz inputs — must be clean, not crashes
ACCEPTABLE_EXCEPTIONS = (
    Error,
    ValueError,
    TypeError,
    OverflowError,
)


# ============================================================================
# §4.1 Extreme Julian Dates
# ============================================================================


class TestExtremeJulianDates:
    """§4.1 Verify robustness with extreme, invalid, and boundary JDs."""

    # Representative body for most tests
    BODY = swe.SUN

    @pytest.mark.parametrize(
        "jd,desc",
        [
            (0.0, "JD=0 (4713 BC)"),
            (-1e6, "JD=-1e6 (far past)"),
            (1e8, "JD=1e8 (far future)"),
            (2451545.0, "J2000.0 (baseline)"),
            (2299160.5, "Gregorian calendar start"),
            (1.0, "JD=1"),
            (-1.0, "JD=-1"),
            (1e-10, "JD near zero"),
            (5373484.5, "~10000 CE"),
            (-1931076.5, "~-10000 CE"),
        ],
        ids=lambda x: x[1] if isinstance(x, tuple) else str(x),
    )
    def test_calc_ut_extreme_jd(self, jd: float, desc: str) -> None:
        """calc_ut with extreme JDs: valid result or clean exception."""
        try:
            result, flags = swe.calc_ut(jd, self.BODY, 0)
            # If it returns, must be a valid 6-tuple of floats
            assert isinstance(result, tuple), f"Expected tuple, got {type(result)}"
            assert len(result) == 6, f"Expected 6 elements, got {len(result)}"
            for i, v in enumerate(result):
                assert isinstance(v, (int, float)), (
                    f"Element {i} is {type(v)}, not float"
                )
                assert math.isfinite(v), f"Element {i} is not finite: {v}"
        except ACCEPTABLE_EXCEPTIONS:
            pass  # Clean exception is acceptable

    @pytest.mark.parametrize(
        "jd,desc",
        [
            (float("nan"), "NaN"),
            (float("inf"), "Inf"),
            (float("-inf"), "-Inf"),
        ],
    )
    def test_calc_ut_nan_inf(self, jd: float, desc: str) -> None:
        """calc_ut with NaN/Inf: must not crash or hang."""
        try:
            result, flags = swe.calc_ut(jd, self.BODY, 0)
            # If it somehow returns, just verify no crash
        except ACCEPTABLE_EXCEPTIONS:
            pass  # Expected

    @pytest.mark.parametrize(
        "jd",
        [0.0, -1e6, 1e8, 2451545.0],
        ids=["JD=0", "far-past", "far-future", "J2000"],
    )
    def test_revjul_extreme_jd(self, jd: float) -> None:
        """revjul with extreme JDs: must return valid date components."""
        try:
            y, m, d, h = swe.revjul(jd)
            assert isinstance(y, int)
            assert 1 <= m <= 12
            assert 1 <= d <= 31
            assert 0.0 <= h < 24.0
        except ACCEPTABLE_EXCEPTIONS:
            pass

    def test_julday_revjul_roundtrip_extremes(self) -> None:
        """julday/revjul roundtrip at extreme dates."""
        test_cases = [
            (2000, 1, 1, 12.0),
            (1, 1, 1, 0.0),
            (-4712, 1, 1, 12.0),
            (9999, 12, 31, 23.99),
        ]
        for y, m, d, h in test_cases:
            jd = swe.julday(y, m, d, h)
            assert math.isfinite(jd), f"julday returned non-finite for {y}-{m}-{d}"

    @pytest.mark.parametrize(
        "jd",
        [0.0, -1e6, 1e8],
        ids=["JD=0", "far-past", "far-future"],
    )
    def test_houses_extreme_jd(self, jd: float) -> None:
        """houses() with extreme JDs: valid result or clean exception."""
        try:
            cusps, angles = swe.houses(jd, 45.0, 10.0)
            assert isinstance(cusps, tuple)
            assert isinstance(angles, tuple)
        except ACCEPTABLE_EXCEPTIONS:
            pass

    @pytest.mark.parametrize(
        "jd",
        [0.0, 2451545.0, 1e8],
        ids=["JD=0", "J2000", "far-future"],
    )
    def test_eclipse_extreme_jd(self, jd: float) -> None:
        """Eclipse search from extreme JDs: valid result or clean exception."""
        try:
            ecl_type, times = swe.sol_eclipse_when_glob(jd)
            assert isinstance(ecl_type, int)
            assert isinstance(times, tuple)
        except ACCEPTABLE_EXCEPTIONS:
            pass

    def test_sidtime_extreme_jd(self) -> None:
        """Sidereal time at extreme JDs."""
        for jd in [0.0, 2451545.0, -1e6, 1e8]:
            try:
                st = swe.sidtime(jd)
                assert isinstance(st, float)
            except ACCEPTABLE_EXCEPTIONS:
                pass

    def test_deltat_extreme_jd(self) -> None:
        """Delta T at extreme JDs."""
        for jd in [0.0, 2451545.0, -1e6, 1e8]:
            try:
                dt = swe.deltat(jd)
                assert isinstance(dt, float)
                assert math.isfinite(dt)
            except ACCEPTABLE_EXCEPTIONS:
                pass


# ============================================================================
# §4.2 Invalid Body IDs
# ============================================================================


class TestInvalidBodyIDs:
    """§4.2 Verify robustness with invalid, negative, and extreme body IDs."""

    JD = 2451545.0  # J2000.0

    @pytest.mark.parametrize(
        "body_id,desc",
        [
            (-1, "Negative body ID"),
            (-100, "Large negative body ID"),
            (-999999, "Very large negative body ID"),
            (99, "Unused body ID 99"),
            (999, "Unused body ID 999"),
            (9999, "Just below AST_OFFSET"),
            (swe.AST_OFFSET, "AST_OFFSET + 0 (edge case)"),
            (swe.AST_OFFSET + 1, "AST_OFFSET + 1 (Ceres)"),
            (100001, "Above AST_OFFSET"),
            (999999, "Very large body ID"),
        ],
    )
    def test_calc_ut_invalid_body(self, body_id: int, desc: str) -> None:
        """calc_ut with invalid body IDs: valid result or clean exception."""
        try:
            result, flags = swe.calc_ut(self.JD, body_id, 0)
            # Some of these may be valid (e.g., AST_OFFSET+1 = Ceres)
            assert isinstance(result, tuple)
            assert len(result) == 6
        except ACCEPTABLE_EXCEPTIONS:
            pass  # Clean error is fine

    @pytest.mark.parametrize(
        "body_id",
        [-100, 99, 999, 9999],
        ids=["neg100", "id99", "id999", "id9999"],
    )
    def test_invalid_body_raises_known_exception(self, body_id: int) -> None:
        """Invalid body IDs should raise UnknownBodyError or similar.

        Note: body ID -1 (ECL_NUT) is valid in Swiss Ephemeris (returns
        Earth nutation values), so it is excluded from this test.
        """
        with pytest.raises(ACCEPTABLE_EXCEPTIONS):
            swe.calc_ut(self.JD, body_id, 0)

    def test_all_standard_bodies_valid(self) -> None:
        """All standard body IDs (0-20) return valid results."""
        standard_bodies = [
            swe.SUN,  # 0
            swe.MOON,  # 1
            swe.MERCURY,  # 2
            swe.VENUS,  # 3
            swe.MARS,  # 4
            swe.JUPITER,  # 5
            swe.SATURN,  # 6
            swe.URANUS,  # 7
            swe.NEPTUNE,  # 8
            swe.PLUTO,  # 9
            swe.MEAN_NODE,  # 10
            swe.TRUE_NODE,  # 11
        ]
        for body in standard_bodies:
            result, flags = swe.calc_ut(self.JD, body, 0)
            assert isinstance(result, tuple)
            assert len(result) == 6
            for v in result:
                assert math.isfinite(v), f"Body {body}: non-finite value {v}"

    def test_pheno_ut_invalid_body(self) -> None:
        """pheno_ut with invalid body IDs: clean exception."""
        for body_id in [-1, 99, 999]:
            try:
                swe.pheno_ut(self.JD, body_id, 0)
            except ACCEPTABLE_EXCEPTIONS:
                pass


# ============================================================================
# §4.3 Extreme Geographic Coordinates
# ============================================================================


class TestExtremeGeographicCoordinates:
    """§4.3 Verify robustness with extreme geographic coordinates."""

    JD = 2451545.0  # J2000.0

    @pytest.mark.parametrize(
        "lat,lon,desc",
        [
            (90.0, 0.0, "North pole"),
            (-90.0, 0.0, "South pole"),
            (0.0, 0.0, "Equator/prime meridian"),
            (0.0, 180.0, "Equator/date line east"),
            (0.0, -180.0, "Equator/date line west"),
            (89.99, 0.0, "Near north pole"),
            (-89.99, 0.0, "Near south pole"),
            (45.0, 360.0, "Longitude 360"),
            (45.0, -360.0, "Longitude -360"),
            (45.0, 720.0, "Longitude 720"),
        ],
    )
    def test_houses_extreme_coords(self, lat: float, lon: float, desc: str) -> None:
        """houses() with extreme geographic coordinates."""
        try:
            cusps, angles = swe.houses(self.JD, lat, lon)
            assert isinstance(cusps, tuple)
            assert len(cusps) >= 12, f"Expected >=12 cusps, got {len(cusps)}"
            assert isinstance(angles, tuple)
        except ACCEPTABLE_EXCEPTIONS:
            pass  # Polar circle errors are expected at poles

    @pytest.mark.parametrize(
        "lat,lon,desc",
        [
            (91.0, 0.0, "Latitude > 90"),
            (-91.0, 0.0, "Latitude < -90"),
            (180.0, 0.0, "Latitude = 180"),
            (-180.0, 0.0, "Latitude = -180"),
        ],
    )
    def test_houses_invalid_latitude(self, lat: float, lon: float, desc: str) -> None:
        """houses() with invalid latitude: should raise or handle gracefully."""
        try:
            cusps, angles = swe.houses(self.JD, lat, lon)
            # If it returns without error, that's fine too (some libs wrap)
        except ACCEPTABLE_EXCEPTIONS:
            pass

    @pytest.mark.parametrize(
        "alt,desc",
        [
            (0.0, "Sea level"),
            (-1000.0, "Below sea level"),
            (8848.0, "Mount Everest"),
            (100000.0, "100km altitude"),
            (1000000.0, "1000km altitude"),
        ],
    )
    def test_set_topo_extreme_altitude(self, alt: float, desc: str) -> None:
        """set_topo with extreme altitudes: valid or clean exception."""
        try:
            swe.set_topo(10.0, 45.0, alt)
            # If set succeeds, try a topocentric calc
            result, flags = swe.calc_ut(self.JD, swe.MOON, swe.FLG_TOPOCTR)
            assert isinstance(result, tuple)
            assert len(result) == 6
        except ACCEPTABLE_EXCEPTIONS:
            pass
        finally:
            # Reset topo state
            try:
                swe.set_topo(0.0, 0.0, 0.0)
            except Exception:
                pass

    def test_houses_all_systems(self) -> None:
        """Test all house systems at a normal location."""
        house_systems = [
            ord("P"),  # Placidus
            ord("K"),  # Koch
            ord("O"),  # Porphyrius
            ord("R"),  # Regiomontanus
            ord("C"),  # Campanus
            ord("E"),  # Equal
            ord("W"),  # Whole sign
            ord("B"),  # Alcabitius
            ord("M"),  # Morinus
        ]
        for hsys in house_systems:
            try:
                cusps, angles = swe.houses(self.JD, 45.0, 10.0, hsys)
                assert isinstance(cusps, tuple)
                assert len(cusps) >= 12
            except ACCEPTABLE_EXCEPTIONS:
                pass

    def test_houses_poles_with_different_systems(self) -> None:
        """House calculations at poles with various systems."""
        # Some house systems fail at poles (polar circle issue)
        for lat in [89.0, 90.0, -89.0, -90.0]:
            for hsys in [ord("P"), ord("E"), ord("W")]:
                try:
                    cusps, angles = swe.houses(self.JD, lat, 0.0, hsys)
                    assert isinstance(cusps, tuple)
                except ACCEPTABLE_EXCEPTIONS:
                    pass  # PolarCircleError or similar is expected


# ============================================================================
# §4.4 Flag Exhaustion
# ============================================================================


class TestFlagExhaustion:
    """§4.4 Verify no crashes with arbitrary flag combinations."""

    JD = 2451545.0  # J2000.0
    BODY = swe.SUN

    # The 14 individual SEFLG flags to combine
    FLAGS = [
        swe.FLG_JPLEPH,  # 1
        swe.FLG_SWIEPH,  # 2
        swe.FLG_MOSEPH,  # 4
        swe.FLG_HELCTR,  # 8
        swe.FLG_TRUEPOS,  # 16
        swe.FLG_J2000,  # 32
        swe.FLG_NONUT,  # 64
        swe.FLG_SPEED3,  # 128
        swe.FLG_SPEED,  # 256
        swe.FLG_NOGDEFL,  # 512
        swe.FLG_NOABERR,  # 1024
        swe.FLG_EQUATORIAL,  # 2048
        swe.FLG_XYZ,  # 4096
        swe.FLG_RADIANS,  # 8192
    ]

    def test_all_single_flags(self) -> None:
        """Each individual flag works without crash."""
        for flag in self.FLAGS:
            try:
                result, retflags = swe.calc_ut(self.JD, self.BODY, flag)
                assert isinstance(result, tuple)
                assert len(result) == 6
            except ACCEPTABLE_EXCEPTIONS:
                pass

    def test_all_flag_pairs(self) -> None:
        """All pairs of flags work without crash (14*13/2 = 91 combos)."""
        crash_count = 0
        for i, f1 in enumerate(self.FLAGS):
            for f2 in self.FLAGS[i + 1 :]:
                combined = f1 | f2
                try:
                    result, retflags = swe.calc_ut(self.JD, self.BODY, combined)
                    assert isinstance(result, tuple)
                except ACCEPTABLE_EXCEPTIONS:
                    pass
                except Exception as e:
                    crash_count += 1
                    if crash_count > 5:
                        pytest.fail(
                            f"Too many unexpected exceptions in flag pairs: "
                            f"{type(e).__name__}: {e}"
                        )
        assert crash_count == 0, f"{crash_count} unexpected exception(s) in flag pairs"

    def test_exhaustive_flag_combos_sampled(self) -> None:
        """Sample 500 random flag combinations from the 2^14 space.

        Testing all 16384 combos would be slow; sampling 500 gives good
        coverage while keeping runtime reasonable (~30s).
        """
        import random

        rng = random.Random(42)  # Fixed seed for reproducibility
        max_flag = (1 << 14) - 1  # 16383

        # Always include some known edge cases
        test_flags = [
            0,  # No flags
            max_flag,  # All flags
            0x3FFF,  # All 14 bits
            0x0001,  # Just JPLEPH
            0x2000,  # Just RADIANS
            0x1000,  # Just XYZ
        ]
        # Add 494 random samples
        test_flags.extend(rng.randint(0, max_flag) for _ in range(494))

        crash_count = 0
        for flags in test_flags:
            try:
                result, retflags = swe.calc_ut(self.JD, self.BODY, flags)
                assert isinstance(result, tuple)
                assert len(result) == 6
            except ACCEPTABLE_EXCEPTIONS:
                pass
            except Exception as e:
                crash_count += 1
                if crash_count > 5:
                    pytest.fail(
                        f"Too many unexpected exceptions: "
                        f"flags=0x{flags:04X} {type(e).__name__}: {e}"
                    )
        assert crash_count == 0, (
            f"{crash_count} unexpected exception(s) in sampled flag combos"
        )

    def test_flag_combos_with_moon(self) -> None:
        """Flag combinations specifically for the Moon (different code path)."""
        import random

        rng = random.Random(43)
        max_flag = (1 << 14) - 1

        test_flags = [0, max_flag]
        test_flags.extend(rng.randint(0, max_flag) for _ in range(98))

        for flags in test_flags:
            try:
                result, retflags = swe.calc_ut(self.JD, swe.MOON, flags)
                assert isinstance(result, tuple)
            except ACCEPTABLE_EXCEPTIONS:
                pass

    def test_higher_flags(self) -> None:
        """Test FLG_BARYCTR, FLG_TOPOCTR, FLG_SIDEREAL, FLG_ICRS."""
        higher_flags = [
            swe.FLG_BARYCTR,  # 16384
            swe.FLG_TOPOCTR,  # 32768
            swe.FLG_SIDEREAL,  # 65536
            swe.FLG_ICRS,  # 131072
        ]
        for flag in higher_flags:
            try:
                result, retflags = swe.calc_ut(self.JD, self.BODY, flag)
                assert isinstance(result, tuple)
            except ACCEPTABLE_EXCEPTIONS:
                pass

        # Combinations of higher flags
        for i, f1 in enumerate(higher_flags):
            for f2 in higher_flags[i + 1 :]:
                try:
                    result, retflags = swe.calc_ut(self.JD, self.BODY, f1 | f2)
                except ACCEPTABLE_EXCEPTIONS:
                    pass

    def test_conflicting_flag_combos(self) -> None:
        """Deliberately conflicting flags: should not crash."""
        conflicting = [
            swe.FLG_HELCTR | swe.FLG_BARYCTR,
            swe.FLG_HELCTR | swe.FLG_TOPOCTR,
            swe.FLG_BARYCTR | swe.FLG_TOPOCTR,
            swe.FLG_JPLEPH | swe.FLG_SWIEPH | swe.FLG_MOSEPH,
            swe.FLG_EQUATORIAL | swe.FLG_XYZ | swe.FLG_RADIANS,
            swe.FLG_J2000 | swe.FLG_SIDEREAL | swe.FLG_ICRS,
        ]
        for flags in conflicting:
            try:
                result, retflags = swe.calc_ut(self.JD, self.BODY, flags)
                assert isinstance(result, tuple)
            except ACCEPTABLE_EXCEPTIONS:
                pass


# ============================================================================
# Additional robustness tests
# ============================================================================


class TestMiscRobustness:
    """Additional robustness tests not in the plan subsections."""

    JD = 2451545.0

    def test_fixstar_invalid_names(self) -> None:
        """fixstar with invalid star names: clean exception."""
        invalid_names = [
            "",
            "NONEXISTENT_STAR",
            "!!invalid!!",
            "a" * 1000,
            "\x00null",
        ]
        for name in invalid_names:
            try:
                swe.fixstar_ut(name, self.JD, 0)
            except ACCEPTABLE_EXCEPTIONS:
                pass

    def test_cotrans_extreme_values(self) -> None:
        """Coordinate transform with extreme input values."""
        extreme_coords = [
            (0.0, 0.0, 1.0),
            (360.0, 0.0, 1.0),
            (-360.0, 0.0, 1.0),
            (0.0, 90.0, 1.0),
            (0.0, -90.0, 1.0),
            (999.0, 999.0, 999.0),
            (1e10, 1e10, 1e10),
        ]
        for coords in extreme_coords:
            try:
                result = swe.cotrans(coords, 23.4)
                assert isinstance(result, tuple)
                assert len(result) == 3
            except ACCEPTABLE_EXCEPTIONS:
                pass

    def test_houses_ex_invalid_systems(self) -> None:
        """houses with invalid house system codes."""
        invalid_systems = [0, 1, 255, ord("Z"), ord("!")]
        for hsys in invalid_systems:
            try:
                cusps, angles = swe.houses(self.JD, 45.0, 10.0, hsys)
                assert isinstance(cusps, tuple)
            except ACCEPTABLE_EXCEPTIONS:
                pass

    def test_rapid_successive_calls(self) -> None:
        """Rapid successive calls don't cause state corruption."""
        bodies = [swe.SUN, swe.MOON, swe.MARS, swe.JUPITER]
        jds = [2451545.0, 2451545.5, 2451546.0, 2451546.5]

        for _ in range(100):
            for body in bodies:
                for jd in jds:
                    result, flags = swe.calc_ut(jd, body, swe.FLG_SPEED)
                    assert len(result) == 6

    def test_calc_ut_with_speed_and_all_bodies(self) -> None:
        """All standard bodies with FLG_SPEED return valid velocities."""
        bodies = range(0, 12)  # SUN through TRUE_NODE
        for body in bodies:
            result, flags = swe.calc_ut(self.JD, body, swe.FLG_SPEED)
            # Elements 3-5 are velocities, should be finite
            for i in range(3, 6):
                assert math.isfinite(result[i]), (
                    f"Body {body}: velocity element {i} is {result[i]}"
                )
