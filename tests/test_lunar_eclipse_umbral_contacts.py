"""
Tests for lunar eclipse umbral contact times U1, U2, U3, U4 in libephemeris.

Tests the calc_lunar_eclipse_umbral_first_contact_u1,
calc_lunar_eclipse_umbral_second_contact_u2, calc_lunar_eclipse_umbral_third_contact_u3,
and calc_lunar_eclipse_umbral_fourth_contact_u4 functions which calculate
the umbral contact times for lunar eclipses.

U1: Moon's leading limb first enters Earth's umbral shadow (partial begins)
U2: Moon's trailing limb enters umbra (totality begins)
U3: Moon's leading limb exits umbra (totality ends)
U4: Moon's trailing limb completely exits umbral shadow (partial ends)

Reference data from NASA Eclipse website:
https://eclipse.gsfc.nasa.gov/lunar.html
"""

import pytest

pytestmark = pytest.mark.slow

from libephemeris import (
    julday,
    revjul,
    lun_eclipse_when,
    calc_lunar_eclipse_umbral_first_contact_u1,
    calc_lunar_eclipse_umbral_second_contact_u2,
    calc_lunar_eclipse_umbral_third_contact_u3,
    calc_lunar_eclipse_umbral_fourth_contact_u4,
    ECL_TOTAL,
    ECL_PARTIAL,
    ECL_PENUMBRAL,
    FLG_SWIEPH,
)


class TestLunarEclipseUmbralContactsBasic:
    """Test basic functionality of umbral contact calculations."""

    def test_u1_function_exists_in_eclipse_module(self):
        """Test that calc_lunar_eclipse_umbral_first_contact_u1 is callable."""
        from libephemeris.eclipse import calc_lunar_eclipse_umbral_first_contact_u1

        assert callable(calc_lunar_eclipse_umbral_first_contact_u1)

    def test_u2_function_exists_in_eclipse_module(self):
        """Test that calc_lunar_eclipse_umbral_second_contact_u2 is callable."""
        from libephemeris.eclipse import calc_lunar_eclipse_umbral_second_contact_u2

        assert callable(calc_lunar_eclipse_umbral_second_contact_u2)

    def test_u3_function_exists_in_eclipse_module(self):
        """Test that calc_lunar_eclipse_umbral_third_contact_u3 is callable."""
        from libephemeris.eclipse import calc_lunar_eclipse_umbral_third_contact_u3

        assert callable(calc_lunar_eclipse_umbral_third_contact_u3)

    def test_u4_function_exists_in_eclipse_module(self):
        """Test that calc_lunar_eclipse_umbral_fourth_contact_u4 is callable."""
        from libephemeris.eclipse import calc_lunar_eclipse_umbral_fourth_contact_u4

        assert callable(calc_lunar_eclipse_umbral_fourth_contact_u4)

    def test_u1_exported_from_main_module(self):
        """Test that U1 function is exported from main libephemeris module."""
        from libephemeris import calc_lunar_eclipse_umbral_first_contact_u1

        assert callable(calc_lunar_eclipse_umbral_first_contact_u1)

    def test_u2_exported_from_main_module(self):
        """Test that U2 function is exported from main libephemeris module."""
        from libephemeris import calc_lunar_eclipse_umbral_second_contact_u2

        assert callable(calc_lunar_eclipse_umbral_second_contact_u2)

    def test_u3_exported_from_main_module(self):
        """Test that U3 function is exported from main libephemeris module."""
        from libephemeris import calc_lunar_eclipse_umbral_third_contact_u3

        assert callable(calc_lunar_eclipse_umbral_third_contact_u3)

    def test_u4_exported_from_main_module(self):
        """Test that U4 function is exported from main libephemeris module."""
        from libephemeris import calc_lunar_eclipse_umbral_fourth_contact_u4

        assert callable(calc_lunar_eclipse_umbral_fourth_contact_u4)

    def test_u1_returns_float(self):
        """Test that U1 returns a float value."""
        jd_start = julday(2022, 10, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        result = calc_lunar_eclipse_umbral_first_contact_u1(jd_max)

        assert isinstance(result, float)

    def test_u2_returns_float(self):
        """Test that U2 returns a float value."""
        jd_start = julday(2022, 10, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        result = calc_lunar_eclipse_umbral_second_contact_u2(jd_max)

        assert isinstance(result, float)

    def test_u3_returns_float(self):
        """Test that U3 returns a float value."""
        jd_start = julday(2022, 10, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        result = calc_lunar_eclipse_umbral_third_contact_u3(jd_max)

        assert isinstance(result, float)

    def test_u4_returns_float(self):
        """Test that U4 returns a float value."""
        jd_start = julday(2022, 10, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        result = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max)

        assert isinstance(result, float)

    def test_u1_accepts_flags_parameter(self):
        """Test that U1 accepts optional flags parameter."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        result = calc_lunar_eclipse_umbral_first_contact_u1(jd_max, flags=FLG_SWIEPH)

        assert isinstance(result, float)
        assert result > 0

    def test_u4_accepts_flags_parameter(self):
        """Test that U4 accepts optional flags parameter."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        result = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max, flags=FLG_SWIEPH)

        assert isinstance(result, float)
        assert result > 0


class TestUmbralContactTimingOrder:
    """Test that umbral contact times are in correct chronological order."""

    def test_u1_before_maximum(self):
        """Test that U1 occurs before eclipse maximum."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_u1 = calc_lunar_eclipse_umbral_first_contact_u1(jd_max)

        assert jd_u1 > 0, "U1 should be calculated"
        assert jd_u1 < jd_max, "U1 should be before eclipse maximum"

    def test_u4_after_maximum(self):
        """Test that U4 occurs after eclipse maximum."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_u4 = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max)

        assert jd_u4 > 0, "U4 should be calculated"
        assert jd_u4 > jd_max, "U4 should be after eclipse maximum"

    def test_u1_before_u4(self):
        """Test that U1 occurs before U4."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_u1 = calc_lunar_eclipse_umbral_first_contact_u1(jd_max)
        jd_u4 = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max)

        assert jd_u1 < jd_u4, "U1 should be before U4"

    def test_u2_after_u1_for_total_eclipse(self):
        """Test that U2 occurs after U1 for total eclipse."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_u1 = calc_lunar_eclipse_umbral_first_contact_u1(jd_max)
        jd_u2 = calc_lunar_eclipse_umbral_second_contact_u2(jd_max)

        assert jd_u2 > 0, "U2 should exist for total eclipse"
        assert jd_u2 > jd_u1, "U2 should be after U1"

    def test_u3_before_u4_for_total_eclipse(self):
        """Test that U3 occurs before U4 for total eclipse."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_u3 = calc_lunar_eclipse_umbral_third_contact_u3(jd_max)
        jd_u4 = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max)

        assert jd_u3 > 0, "U3 should exist for total eclipse"
        assert jd_u3 < jd_u4, "U3 should be before U4"

    def test_u2_before_u3_for_total_eclipse(self):
        """Test that U2 occurs before U3 for total eclipse."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_u2 = calc_lunar_eclipse_umbral_second_contact_u2(jd_max)
        jd_u3 = calc_lunar_eclipse_umbral_third_contact_u3(jd_max)

        assert jd_u2 > 0, "U2 should exist for total eclipse"
        assert jd_u3 > 0, "U3 should exist for total eclipse"
        assert jd_u2 < jd_u3, "U2 should be before U3"

    def test_u2_before_maximum_for_total_eclipse(self):
        """Test that U2 occurs before maximum for total eclipse."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_u2 = calc_lunar_eclipse_umbral_second_contact_u2(jd_max)

        assert jd_u2 > 0, "U2 should exist for total eclipse"
        assert jd_u2 < jd_max, "U2 should be before maximum"

    def test_u3_after_maximum_for_total_eclipse(self):
        """Test that U3 occurs after maximum for total eclipse."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_u3 = calc_lunar_eclipse_umbral_third_contact_u3(jd_max)

        assert jd_u3 > 0, "U3 should exist for total eclipse"
        assert jd_u3 > jd_max, "U3 should be after maximum"

    def test_full_contact_order_for_total_eclipse(self):
        """Test complete order: U1 < U2 < max < U3 < U4."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_u1 = calc_lunar_eclipse_umbral_first_contact_u1(jd_max)
        jd_u2 = calc_lunar_eclipse_umbral_second_contact_u2(jd_max)
        jd_u3 = calc_lunar_eclipse_umbral_third_contact_u3(jd_max)
        jd_u4 = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max)

        assert jd_u1 < jd_u2 < jd_max < jd_u3 < jd_u4, (
            f"Order should be U1 < U2 < max < U3 < U4, got: "
            f"U1={jd_u1:.6f}, U2={jd_u2:.6f}, max={jd_max:.6f}, "
            f"U3={jd_u3:.6f}, U4={jd_u4:.6f}"
        )


class TestUmbralContactKnownEclipses:
    """Test umbral contacts against known eclipse data."""

    def test_november_2022_total_lunar_eclipse_umbral_contacts(self):
        """Test November 8, 2022 total lunar eclipse umbral contacts.

        Reference times from lun_eclipse_when (consistent with the reference ephemeris):
        - U1 (Partial begins): 2022 Nov 08 at 09:09:34 UT
        - U2 (Total begins): 2022 Nov 08 at 10:17:26 UT
        - Maximum: 2022 Nov 08 at 10:59:12 UT
        - U3 (Total ends): 2022 Nov 08 at 11:40:59 UT
        - U4 (Partial ends): 2022 Nov 08 at 12:48:50 UT
        """
        jd_start = julday(2022, 10, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        # Verify we found the right eclipse (November 8, 2022)
        year, month, day, _ = revjul(jd_max)
        assert year == 2022
        assert month == 11
        assert day == 8

        jd_u1 = calc_lunar_eclipse_umbral_first_contact_u1(jd_max)
        jd_u2 = calc_lunar_eclipse_umbral_second_contact_u2(jd_max)
        jd_u3 = calc_lunar_eclipse_umbral_third_contact_u3(jd_max)
        jd_u4 = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max)

        # Reference times (converted to decimal hours UT)
        # U1 at 09:09:34 UT = 9.159 hours
        jd_u1_ref = julday(2022, 11, 8, 9.159)

        # U2 at 10:17:26 UT = 10.291 hours
        jd_u2_ref = julday(2022, 11, 8, 10.291)

        # U3 at 11:40:59 UT = 11.683 hours
        jd_u3_ref = julday(2022, 11, 8, 11.683)

        # U4 at 12:48:50 UT = 12.814 hours
        jd_u4_ref = julday(2022, 11, 8, 12.814)

        # Allow 2 minutes tolerance (algorithms may differ slightly)
        tolerance_days = 2.0 / (24 * 60)  # 2 minutes in days

        u1_diff = abs(jd_u1 - jd_u1_ref)
        u2_diff = abs(jd_u2 - jd_u2_ref)
        u3_diff = abs(jd_u3 - jd_u3_ref)
        u4_diff = abs(jd_u4 - jd_u4_ref)

        assert u1_diff < tolerance_days, (
            f"U1 timing error: {u1_diff * 24 * 60:.1f} minutes. Expected < 2 minutes"
        )
        assert u2_diff < tolerance_days, (
            f"U2 timing error: {u2_diff * 24 * 60:.1f} minutes. Expected < 2 minutes"
        )
        assert u3_diff < tolerance_days, (
            f"U3 timing error: {u3_diff * 24 * 60:.1f} minutes. Expected < 2 minutes"
        )
        assert u4_diff < tolerance_days, (
            f"U4 timing error: {u4_diff * 24 * 60:.1f} minutes. Expected < 2 minutes"
        )

    def test_may_2022_total_lunar_eclipse_umbral_contacts(self):
        """Test May 16, 2022 total lunar eclipse umbral contacts.

        NASA Reference:
        - U1 (Partial begins): 2022 May 16 at 02:27:53 UT
        - U2 (Total begins): 2022 May 16 at 03:29:06 UT
        - Maximum: 2022 May 16 at 04:11:28 UT
        - U3 (Total ends): 2022 May 16 at 04:53:54 UT
        - U4 (Partial ends): 2022 May 16 at 05:55:09 UT
        """
        jd_start = julday(2022, 4, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        # Verify we found the right eclipse (May 16, 2022)
        year, month, day, _ = revjul(jd_max)
        assert year == 2022
        assert month == 5
        assert day == 16

        jd_u1 = calc_lunar_eclipse_umbral_first_contact_u1(jd_max)
        jd_u2 = calc_lunar_eclipse_umbral_second_contact_u2(jd_max)
        jd_u3 = calc_lunar_eclipse_umbral_third_contact_u3(jd_max)
        jd_u4 = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max)

        # NASA reference times
        # U1 at 02:27:53 UT = 2.465 hours
        jd_u1_ref = julday(2022, 5, 16, 2.465)

        # U2 at 03:29:06 UT = 3.485 hours
        jd_u2_ref = julday(2022, 5, 16, 3.485)

        # U3 at 04:53:54 UT = 4.898 hours
        jd_u3_ref = julday(2022, 5, 16, 4.898)

        # U4 at 05:55:09 UT = 5.919 hours
        jd_u4_ref = julday(2022, 5, 16, 5.919)

        # Allow 2 minutes tolerance
        tolerance_days = 2.0 / (24 * 60)

        u1_diff = abs(jd_u1 - jd_u1_ref)
        u2_diff = abs(jd_u2 - jd_u2_ref)
        u3_diff = abs(jd_u3 - jd_u3_ref)
        u4_diff = abs(jd_u4 - jd_u4_ref)

        assert u1_diff < tolerance_days, (
            f"U1 timing error: {u1_diff * 24 * 60:.1f} minutes"
        )
        assert u2_diff < tolerance_days, (
            f"U2 timing error: {u2_diff * 24 * 60:.1f} minutes"
        )
        assert u3_diff < tolerance_days, (
            f"U3 timing error: {u3_diff * 24 * 60:.1f} minutes"
        )
        assert u4_diff < tolerance_days, (
            f"U4 timing error: {u4_diff * 24 * 60:.1f} minutes"
        )


class TestUmbralContactsPhysicalConstraints:
    """Test physical constraints on umbral contact times."""

    def test_umbral_duration_reasonable(self):
        """Test that U1 to U4 duration is physically reasonable."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_u1 = calc_lunar_eclipse_umbral_first_contact_u1(jd_max)
        jd_u4 = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max)

        # Total umbral duration in hours
        duration_hours = (jd_u4 - jd_u1) * 24

        # Lunar eclipse umbral phase typically lasts 1-4 hours
        assert 1.0 < duration_hours < 5.0, (
            f"Umbral duration {duration_hours:.2f} hours seems unreasonable. "
            f"Expected 1-5 hours for a typical lunar eclipse"
        )

    def test_totality_duration_reasonable(self):
        """Test that U2 to U3 duration is physically reasonable for total eclipse."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_u2 = calc_lunar_eclipse_umbral_second_contact_u2(jd_max)
        jd_u3 = calc_lunar_eclipse_umbral_third_contact_u3(jd_max)

        # Totality duration in hours
        duration_hours = (jd_u3 - jd_u2) * 24

        # Lunar eclipse totality typically lasts 0.25-2 hours
        assert 0.2 < duration_hours < 2.5, (
            f"Totality duration {duration_hours:.2f} hours seems unreasonable. "
            f"Expected 0.25-2.5 hours for a total lunar eclipse"
        )


class TestUmbralContactsWithPartialEclipse:
    """Test umbral contacts for partial lunar eclipses."""

    def test_partial_eclipse_has_u1_and_u4(self):
        """Test that partial eclipses have valid U1 and U4."""
        jd_start = julday(2023, 10, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PARTIAL)
        jd_max = times[0]

        jd_u1 = calc_lunar_eclipse_umbral_first_contact_u1(jd_max)
        jd_u4 = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max)

        assert jd_u1 > 0, "U1 should exist for partial eclipse"
        assert jd_u4 > 0, "U4 should exist for partial eclipse"
        assert jd_u1 < jd_max < jd_u4, "U1 < max < U4 ordering should hold"

    def test_partial_eclipse_no_u2_u3(self):
        """Test that partial eclipses do not have U2 and U3 (no totality)."""
        jd_start = julday(2023, 10, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PARTIAL)
        jd_max = times[0]

        jd_u2 = calc_lunar_eclipse_umbral_second_contact_u2(jd_max)
        jd_u3 = calc_lunar_eclipse_umbral_third_contact_u3(jd_max)

        # For partial eclipses, U2 and U3 should be 0.0 (no totality)
        assert jd_u2 == 0.0, "U2 should be 0.0 for partial eclipse (no totality)"
        assert jd_u3 == 0.0, "U3 should be 0.0 for partial eclipse (no totality)"


class TestUmbralContactsWithPenumbralEclipse:
    """Test umbral contacts for penumbral lunar eclipses."""

    def test_penumbral_eclipse_no_umbral_contacts(self):
        """Test that penumbral-only eclipses have no umbral contacts."""
        jd_start = julday(2020, 1, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PENUMBRAL)
        jd_max = times[0]

        # Only test if we found a pure penumbral eclipse
        if ecl_type == ECL_PENUMBRAL:
            jd_u1 = calc_lunar_eclipse_umbral_first_contact_u1(jd_max)
            jd_u2 = calc_lunar_eclipse_umbral_second_contact_u2(jd_max)
            jd_u3 = calc_lunar_eclipse_umbral_third_contact_u3(jd_max)
            jd_u4 = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max)

            # All umbral contacts should be 0.0 for penumbral-only eclipse
            assert jd_u1 == 0.0, "U1 should be 0.0 for penumbral-only eclipse"
            assert jd_u2 == 0.0, "U2 should be 0.0 for penumbral-only eclipse"
            assert jd_u3 == 0.0, "U3 should be 0.0 for penumbral-only eclipse"
            assert jd_u4 == 0.0, "U4 should be 0.0 for penumbral-only eclipse"


class TestUmbralContactsMultipleEclipses:
    """Test umbral contacts for multiple sequential eclipses."""

    def test_multiple_total_eclipses_have_valid_contacts(self):
        """Test that U1, U2, U3, U4 work correctly for multiple total eclipses."""
        jd = julday(2022, 1, 1, 0.0)
        eclipses = []

        # Find 2 sequential total lunar eclipses
        for _ in range(2):
            ecl_type, times = lun_eclipse_when(jd, ecltype=ECL_TOTAL)
            jd_max = times[0]

            jd_u1 = calc_lunar_eclipse_umbral_first_contact_u1(jd_max)
            jd_u2 = calc_lunar_eclipse_umbral_second_contact_u2(jd_max)
            jd_u3 = calc_lunar_eclipse_umbral_third_contact_u3(jd_max)
            jd_u4 = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max)

            eclipses.append(
                {
                    "max": jd_max,
                    "u1": jd_u1,
                    "u2": jd_u2,
                    "u3": jd_u3,
                    "u4": jd_u4,
                }
            )

            jd = jd_max + 1  # Move past this eclipse

        # Verify each eclipse has valid umbral contacts
        for i, ecl in enumerate(eclipses):
            assert ecl["u1"] > 0, f"Eclipse {i + 1}: U1 should be > 0"
            assert ecl["u2"] > 0, f"Eclipse {i + 1}: U2 should be > 0"
            assert ecl["u3"] > 0, f"Eclipse {i + 1}: U3 should be > 0"
            assert ecl["u4"] > 0, f"Eclipse {i + 1}: U4 should be > 0"
            assert ecl["u1"] < ecl["u2"], f"Eclipse {i + 1}: U1 should be before U2"
            assert ecl["u2"] < ecl["max"], f"Eclipse {i + 1}: U2 should be before max"
            assert ecl["max"] < ecl["u3"], f"Eclipse {i + 1}: max should be before U3"
            assert ecl["u3"] < ecl["u4"], f"Eclipse {i + 1}: U3 should be before U4"


class TestUmbralContactsWithLunEclipseWhen:
    """Test that umbral contacts match the times returned by lun_eclipse_when."""

    def test_u1_matches_lun_eclipse_when_times2(self):
        """Test that U1 closely matches times[2] from lun_eclipse_when."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        # times[2] is partial eclipse beginning from lun_eclipse_when
        jd_partial_begin_from_when = times[2]
        jd_u1 = calc_lunar_eclipse_umbral_first_contact_u1(jd_max)

        # Both should calculate the same event
        # Allow 2 minute tolerance
        tolerance_days = 2.0 / (24 * 60)
        diff = abs(jd_u1 - jd_partial_begin_from_when)

        assert diff < tolerance_days, (
            f"U1 and times[2] differ by {diff * 24 * 60:.1f} minutes. "
            f"Expected < 2 minutes"
        )

    def test_u4_matches_lun_eclipse_when_times3(self):
        """Test that U4 closely matches times[3] from lun_eclipse_when."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        # times[3] is partial eclipse ending from lun_eclipse_when
        jd_partial_end_from_when = times[3]
        jd_u4 = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max)

        # Both should calculate the same event
        tolerance_days = 2.0 / (24 * 60)
        diff = abs(jd_u4 - jd_partial_end_from_when)

        assert diff < tolerance_days, (
            f"U4 and times[3] differ by {diff * 24 * 60:.1f} minutes. "
            f"Expected < 2 minutes"
        )

    def test_u2_matches_lun_eclipse_when_times4(self):
        """Test that U2 closely matches times[4] from lun_eclipse_when."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        # times[4] is total eclipse beginning from lun_eclipse_when
        jd_total_begin_from_when = times[4]
        jd_u2 = calc_lunar_eclipse_umbral_second_contact_u2(jd_max)

        # Both should calculate the same event
        tolerance_days = 2.0 / (24 * 60)
        diff = abs(jd_u2 - jd_total_begin_from_when)

        assert diff < tolerance_days, (
            f"U2 and times[4] differ by {diff * 24 * 60:.1f} minutes. "
            f"Expected < 2 minutes"
        )

    def test_u3_matches_lun_eclipse_when_times5(self):
        """Test that U3 closely matches times[5] from lun_eclipse_when."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        # times[5] is total eclipse ending from lun_eclipse_when
        jd_total_end_from_when = times[5]
        jd_u3 = calc_lunar_eclipse_umbral_third_contact_u3(jd_max)

        # Both should calculate the same event
        tolerance_days = 2.0 / (24 * 60)
        diff = abs(jd_u3 - jd_total_end_from_when)

        assert diff < tolerance_days, (
            f"U3 and times[5] differ by {diff * 24 * 60:.1f} minutes. "
            f"Expected < 2 minutes"
        )
