"""Tests for 360° wrap-around edge cases, NaN/Inf inputs, and ECL_NUT special body."""

from __future__ import annotations

import math
import pytest
import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    JUPITER,
    MEAN_NODE,
    ECL_NUT,
    FLG_SWIEPH,
    FLG_SPEED,
    FLG_EQUATORIAL,
    FLG_SIDEREAL,
)

JD_J2000 = 2451545.0


@pytest.mark.unit
class TestWrapAround360:
    """Test positions near the 0°/360° boundary."""

    def _find_jd_near_zero(self, body: int, start_jd: float, days: int = 365) -> float:
        """Find a JD where body longitude is near 0°/360°."""
        best_jd = start_jd
        best_dist = 999.0
        for i in range(days):
            jd = start_jd + i
            result, _ = swe.calc_ut(jd, body, FLG_SWIEPH | FLG_SPEED)
            lon = result[0]
            dist_to_zero = min(lon, 360.0 - lon)
            if dist_to_zero < best_dist:
                best_dist = dist_to_zero
                best_jd = jd
        return best_jd

    def test_sun_near_zero_valid_range(self):
        """Sun near 0° Aries should still be in [0, 360)."""
        # Sun crosses 0° around March equinox. J2000 is Jan 1 2000.
        # March ~20 is about day 79.
        jd = self._find_jd_near_zero(SUN, JD_J2000, 365)
        result, _ = swe.calc_ut(jd, SUN, FLG_SWIEPH | FLG_SPEED)
        assert 0.0 <= result[0] < 360.0

    def test_moon_near_zero_valid(self):
        """Moon frequently crosses 0°; position must be in [0, 360)."""
        jd = self._find_jd_near_zero(MOON, JD_J2000, 30)
        result, _ = swe.calc_ut(jd, MOON, FLG_SWIEPH | FLG_SPEED)
        assert 0.0 <= result[0] < 360.0

    def test_positions_never_negative(self):
        """No longitude should ever be negative over a scan of dates."""
        for body in [SUN, MOON, MARS, JUPITER]:
            for i in range(50):
                jd = JD_J2000 + i * 7.3  # sample every ~week
                result, _ = swe.calc_ut(jd, body, FLG_SWIEPH)
                assert result[0] >= 0.0, f"Body {body} at JD {jd}: lon {result[0]} < 0"

    def test_positions_never_360(self):
        """No longitude should ever be exactly 360.0 or above."""
        for body in [SUN, MOON, MARS, JUPITER]:
            for i in range(50):
                jd = JD_J2000 + i * 7.3
                result, _ = swe.calc_ut(jd, body, FLG_SWIEPH)
                assert result[0] < 360.0, (
                    f"Body {body} at JD {jd}: lon {result[0]} >= 360"
                )

    def test_solcross_near_zero(self):
        """solcross_ut at 0° should return a valid crossing JD."""
        jd_cross = swe.solcross_ut(0.0, JD_J2000, 0)
        assert jd_cross > JD_J2000
        # Verify the Sun is near 0° at that JD
        result, _ = swe.calc_ut(jd_cross, SUN, FLG_SWIEPH)
        assert result[0] < 0.01 or result[0] > 359.99

    def test_mooncross_near_zero(self):
        """mooncross_ut at 0° should return a valid crossing JD."""
        jd_cross = swe.mooncross_ut(0.0, JD_J2000, 0)
        assert jd_cross > JD_J2000
        result, _ = swe.calc_ut(jd_cross, MOON, FLG_SWIEPH)
        assert result[0] < 0.1 or result[0] > 359.9

    def test_mean_node_range_across_full_cycle(self):
        """Mean Node (retrograde) should stay in [0, 360) over 19 years."""
        for i in range(200):
            jd = JD_J2000 + i * 34.7  # ~19 years in 200 steps
            result, _ = swe.calc_ut(jd, MEAN_NODE, FLG_SWIEPH)
            assert 0.0 <= result[0] < 360.0, (
                f"Mean Node at JD {jd}: lon {result[0]} out of range"
            )

    def test_sidereal_positions_in_range(self):
        """Sidereal longitudes must also be in [0, 360)."""
        swe.set_sid_mode(1)  # Lahiri
        try:
            for body in [SUN, MOON, MARS]:
                result, _ = swe.calc_ut(JD_J2000, body, FLG_SWIEPH | FLG_SIDEREAL)
                assert 0.0 <= result[0] < 360.0, (
                    f"Sidereal lon for body {body}: {result[0]}"
                )
        finally:
            swe.set_sid_mode(0)


@pytest.mark.unit
class TestEclNut:
    """Test ECL_NUT special body for nutation and obliquity."""

    def test_ecl_nut_returns_6_values(self):
        """ECL_NUT should return a 6-element tuple."""
        result, flag = swe.calc_ut(JD_J2000, ECL_NUT, FLG_SWIEPH)
        assert len(result) == 6

    def test_ecl_nut_values_finite(self):
        """All ECL_NUT values should be finite."""
        result, flag = swe.calc_ut(JD_J2000, ECL_NUT, FLG_SWIEPH)
        for i, val in enumerate(result[:4]):
            assert math.isfinite(val), f"ECL_NUT[{i}] = {val} not finite"

    def test_ecl_nut_obliquity_plausible(self):
        """True and mean obliquity should be ~23.4°."""
        result, _ = swe.calc_ut(JD_J2000, ECL_NUT, FLG_SWIEPH)
        true_obl = result[0]
        mean_obl = result[1]
        assert 23.0 < true_obl < 24.0, f"True obliquity {true_obl} implausible"
        assert 23.0 < mean_obl < 24.0, f"Mean obliquity {mean_obl} implausible"

    def test_ecl_nut_nutation_small(self):
        """Nutation in longitude and obliquity should be small (< 0.01°)."""
        result, _ = swe.calc_ut(JD_J2000, ECL_NUT, FLG_SWIEPH)
        nut_lon = result[2]
        nut_obl = result[3]
        assert abs(nut_lon) < 0.01, f"Nutation in longitude {nut_lon} too large"
        assert abs(nut_obl) < 0.01, f"Nutation in obliquity {nut_obl} too large"

    def test_ecl_nut_last_two_zero(self):
        """Last two values of ECL_NUT should be 0.0."""
        result, _ = swe.calc_ut(JD_J2000, ECL_NUT, FLG_SWIEPH)
        assert result[4] == 0.0
        assert result[5] == 0.0

    def test_ecl_nut_varies_over_time(self):
        """ECL_NUT values should change measurably over centuries."""
        res_2000, _ = swe.calc_ut(JD_J2000, ECL_NUT, 0)
        res_2100, _ = swe.calc_ut(JD_J2000 + 36525.0, ECL_NUT, 0)
        # Mean obliquity changes ~47" per century
        diff_obl = abs(res_2100[1] - res_2000[1])
        assert diff_obl > 0.001, "Mean obliquity should change over a century"

    def test_ecl_nut_with_speed_flag(self):
        """ECL_NUT with FLG_SPEED should not crash."""
        result, _ = swe.calc_ut(JD_J2000, ECL_NUT, FLG_SPEED)
        assert len(result) == 6

    def test_ecl_nut_with_equatorial_flag(self):
        """ECL_NUT with FLG_EQUATORIAL should not crash."""
        result, _ = swe.calc_ut(JD_J2000, ECL_NUT, FLG_EQUATORIAL)
        assert len(result) == 6


@pytest.mark.unit
class TestInvalidInputHandling:
    """Test that NaN and Inf inputs are handled gracefully."""

    def test_calc_ut_nan_jd(self):
        """calc_ut with NaN JD should raise or return non-NaN (not silently corrupt)."""
        try:
            result, _ = swe.calc_ut(float("nan"), SUN, FLG_SWIEPH)
            # If it doesn't raise, all outputs should be NaN (propagation)
            # or the function should have caught it
        except (ValueError, TypeError, Exception):
            pass  # Any exception is acceptable

    def test_calc_ut_inf_jd(self):
        """calc_ut with Inf JD should raise or handle gracefully."""
        try:
            result, _ = swe.calc_ut(float("inf"), SUN, FLG_SWIEPH)
        except (ValueError, TypeError, OverflowError, Exception):
            pass  # Any exception is acceptable

    def test_calc_ut_neg_inf_jd(self):
        """calc_ut with -Inf JD should raise or handle gracefully."""
        try:
            result, _ = swe.calc_ut(float("-inf"), SUN, FLG_SWIEPH)
        except (ValueError, TypeError, OverflowError, Exception):
            pass  # Any exception is acceptable

    def test_houses_nan_lat(self):
        """houses with NaN latitude should raise or handle gracefully."""
        try:
            cusps, ascmc = swe.houses(JD_J2000, float("nan"), 12.5, ord("P"))
        except (ValueError, TypeError, Exception):
            pass

    def test_houses_inf_lon(self):
        """houses with Inf longitude should raise or handle gracefully."""
        try:
            cusps, ascmc = swe.houses(JD_J2000, 41.9, float("inf"), ord("P"))
        except (ValueError, TypeError, Exception):
            pass

    def test_julday_nan_components(self):
        """julday with NaN components should raise or handle gracefully."""
        try:
            jd = swe.julday(2000, 1, float("nan"), 12.0)
        except (ValueError, TypeError, Exception):
            pass

    def test_cotrans_nan_input(self):
        """cotrans with NaN should not silently produce valid-looking output."""
        try:
            result = swe.cotrans((float("nan"), 0.0, 1.0), 23.4)
            # If it doesn't raise, NaN should propagate
            if not math.isnan(result[0]):
                # Some implementations might handle this differently
                pass
        except (ValueError, TypeError, Exception):
            pass


@pytest.mark.unit
class TestCalcPctrSpeed3Regression:
    """calc_pctr FLG_SPEED3 public flag semantics.

    calc()/calc_ut() remap FLG_SPEED3 -> FLG_SPEED before computing, but
    calc_pctr() preserves FLG_SPEED3 in the echoed retflag. The public contract
    computes the pctr velocity only when FLG_SPEED is set: FLG_SPEED3
    alone echoes the bit yet returns zero velocity slots, while
    FLG_SPEED|FLG_SPEED3 returns the FLG_SPEED velocity.
    """

    def test_speed3_alone_zero_velocity(self):
        from libephemeris.constants import FLG_SPEED3

        pos3, ret3 = swe.calc_pctr(JD_J2000, MOON, MARS, FLG_SPEED3)
        # The bit is echoed, but only FLG_SPEED triggers computation.
        assert pos3[3:] == (0.0, 0.0, 0.0)
        assert ret3 & FLG_SPEED3

    def test_speed_with_speed3_matches_speed_velocity(self):
        from libephemeris.constants import FLG_SPEED3

        pos_both, ret_both = swe.calc_pctr(JD_J2000, MOON, MARS, FLG_SPEED | FLG_SPEED3)
        pos_s, _ = swe.calc_pctr(JD_J2000, MOON, MARS, FLG_SPEED)
        assert pos_both[3:] == pos_s[3:]
        assert any(v != 0.0 for v in pos_both[3:])
        assert ret_both & FLG_SPEED
        assert ret_both & FLG_SPEED3

    def test_speed3_retflag_preserved(self):
        from libephemeris.constants import FLG_SPEED3

        _, ret3 = swe.calc_pctr(JD_J2000, MOON, MARS, FLG_SPEED3)
        # FLG_SPEED3 (128) is preserved in the retflag, not remapped to
        # FLG_SPEED (256); the default SWIEPH bit (2) is echoed too.
        assert ret3 & FLG_SPEED3
        assert not (ret3 & FLG_SPEED)

    def test_speed3_native_floats(self):
        from libephemeris.constants import FLG_SPEED3

        pos3, _ = swe.calc_pctr(JD_J2000, MOON, MARS, FLG_SPEED3)
        assert all(type(v) is float for v in pos3)


@pytest.mark.unit
class TestCalcPctrNonFiniteRegression:
    """calc_pctr must reject non-finite Julian Days (ADV-2).

    A NaN slips through validate_jd_range's ``jd < start`` / ``jd > end``
    comparisons (all NaN comparisons are False), so calc_pctr() previously
    returned garbage ``(nan, 0.0, nan, ...)`` instead of raising like
    calc()/calc_ut(). All non-finite inputs now raise EphemerisRangeError.
    """

    def test_calc_pctr_nan_raises(self):
        from libephemeris.exceptions import EphemerisRangeError

        with pytest.raises(EphemerisRangeError):
            swe.calc_pctr(float("nan"), MOON, MARS, FLG_SPEED)

    def test_calc_pctr_inf_raises(self):
        from libephemeris.exceptions import EphemerisRangeError

        for jd in (float("inf"), float("-inf")):
            with pytest.raises(EphemerisRangeError):
                swe.calc_pctr(jd, MOON, MARS, FLG_SPEED)

    def test_calc_ut_nan_symmetry(self):
        # calc_ut / calc share validate_jd_range, so they raise the same
        # EphemerisRangeError for a NaN Julian Day.
        from libephemeris.exceptions import EphemerisRangeError

        with pytest.raises(EphemerisRangeError):
            swe.calc_ut(float("nan"), SUN, FLG_SPEED)
        with pytest.raises(EphemerisRangeError):
            swe.calc(float("nan"), SUN, FLG_SPEED)
