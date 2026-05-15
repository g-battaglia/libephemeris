"""
Tests for lunar eclipse penumbral contact times P1 and P4 in libephemeris.

Tests the calc_lunar_eclipse_penumbral_first_contact_p1 and
calc_lunar_eclipse_penumbral_fourth_contact_p4 functions which calculate
the penumbral contact times for lunar eclipses.

P1: Moon's leading limb first enters Earth's penumbral shadow
P4: Moon's trailing limb completely exits Earth's penumbral shadow

Reference data from NASA Eclipse website:
https://eclipse.gsfc.nasa.gov/lunar.html
"""

import pytest

pytestmark = pytest.mark.slow

from libephemeris import (
    julday,
    revjul,
    lun_eclipse_when,
    calc_lunar_eclipse_penumbral_first_contact_p1,
    calc_lunar_eclipse_penumbral_fourth_contact_p4,
    ECL_TOTAL,
    ECL_PARTIAL,
    ECL_PENUMBRAL,
    FLG_SWIEPH,
)


class TestLunarEclipsePenumbralContactsBasic:
    """Test basic functionality of penumbral contact calculations."""

    def test_p1_function_exists_in_eclipse_module(self):
        """Test that calc_lunar_eclipse_penumbral_first_contact_p1 is callable."""
        from libephemeris.eclipse import calc_lunar_eclipse_penumbral_first_contact_p1

        assert callable(calc_lunar_eclipse_penumbral_first_contact_p1)

    def test_p4_function_exists_in_eclipse_module(self):
        """Test that calc_lunar_eclipse_penumbral_fourth_contact_p4 is callable."""
        from libephemeris.eclipse import calc_lunar_eclipse_penumbral_fourth_contact_p4

        assert callable(calc_lunar_eclipse_penumbral_fourth_contact_p4)

    def test_p1_exported_from_main_module(self):
        """Test that P1 function is exported from main libephemeris module."""
        from libephemeris import calc_lunar_eclipse_penumbral_first_contact_p1

        assert callable(calc_lunar_eclipse_penumbral_first_contact_p1)

    def test_p4_exported_from_main_module(self):
        """Test that P4 function is exported from main libephemeris module."""
        from libephemeris import calc_lunar_eclipse_penumbral_fourth_contact_p4

        assert callable(calc_lunar_eclipse_penumbral_fourth_contact_p4)

    def test_p1_returns_float(self):
        """Test that P1 returns a float value."""
        jd_start = julday(2022, 10, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        result = calc_lunar_eclipse_penumbral_first_contact_p1(jd_max)

        assert isinstance(result, float)

    def test_p4_returns_float(self):
        """Test that P4 returns a float value."""
        jd_start = julday(2022, 10, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        result = calc_lunar_eclipse_penumbral_fourth_contact_p4(jd_max)

        assert isinstance(result, float)

    def test_p1_accepts_flags_parameter(self):
        """Test that P1 accepts optional flags parameter."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        result = calc_lunar_eclipse_penumbral_first_contact_p1(
            jd_max, flags=FLG_SWIEPH
        )

        assert isinstance(result, float)
        assert result > 0

    def test_p4_accepts_flags_parameter(self):
        """Test that P4 accepts optional flags parameter."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        result = calc_lunar_eclipse_penumbral_fourth_contact_p4(
            jd_max, flags=FLG_SWIEPH
        )

        assert isinstance(result, float)
        assert result > 0


class TestPenumbralContactTimingOrder:
    """Test that penumbral contact times are in correct chronological order."""

    def test_p1_before_maximum(self):
        """Test that P1 occurs before eclipse maximum."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_p1 = calc_lunar_eclipse_penumbral_first_contact_p1(jd_max)

        assert jd_p1 > 0, "P1 should be calculated"
        assert jd_p1 < jd_max, "P1 should be before eclipse maximum"

    def test_p4_after_maximum(self):
        """Test that P4 occurs after eclipse maximum."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_p4 = calc_lunar_eclipse_penumbral_fourth_contact_p4(jd_max)

        assert jd_p4 > 0, "P4 should be calculated"
        assert jd_p4 > jd_max, "P4 should be after eclipse maximum"

    def test_p1_before_p4(self):
        """Test that P1 occurs before P4."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_p1 = calc_lunar_eclipse_penumbral_first_contact_p1(jd_max)
        jd_p4 = calc_lunar_eclipse_penumbral_fourth_contact_p4(jd_max)

        assert jd_p1 < jd_p4, "P1 should be before P4"

    def test_p1_p4_symmetric_around_maximum(self):
        """Test that P1 and P4 are roughly symmetric around maximum."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_p1 = calc_lunar_eclipse_penumbral_first_contact_p1(jd_max)
        jd_p4 = calc_lunar_eclipse_penumbral_fourth_contact_p4(jd_max)

        # Time from P1 to max and max to P4 should be similar (within 10%)
        t_p1_to_max = jd_max - jd_p1
        t_max_to_p4 = jd_p4 - jd_max

        ratio = t_p1_to_max / t_max_to_p4
        assert 0.8 < ratio < 1.2, (
            f"P1 and P4 should be symmetric around maximum (ratio={ratio:.3f})"
        )


class TestPenumbralContactKnownEclipses:
    """Test penumbral contacts against known eclipse data."""

    def test_november_2022_total_lunar_eclipse_penumbral_contacts(self):
        """Test November 8, 2022 total lunar eclipse penumbral contacts.

        NASA Reference (eclipse.gsfc.nasa.gov):
        - P1 (Penumbral begins): 2022 Nov 08 at 08:02:17 UT
        - Maximum: 2022 Nov 08 at 10:59:08 UT
        - P4 (Penumbral ends): 2022 Nov 08 at 13:56:00 UT
        """
        jd_start = julday(2022, 10, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        # Verify we found the right eclipse (November 8, 2022)
        year, month, day, _ = revjul(jd_max)
        assert year == 2022
        assert month == 11
        assert day == 8

        jd_p1 = calc_lunar_eclipse_penumbral_first_contact_p1(jd_max)
        jd_p4 = calc_lunar_eclipse_penumbral_fourth_contact_p4(jd_max)

        # NASA reference: P1 at 08:02:17 UT
        jd_p1_ref = julday(2022, 11, 8, 8.038)  # 08:02:17 UT = 8.038 hours

        # NASA reference: P4 at 13:56:00 UT
        jd_p4_ref = julday(2022, 11, 8, 13.933)  # 13:56:00 UT = 13.933 hours

        # Allow 2 minutes tolerance (NASA times may use slightly different algorithms)
        tolerance_days = 2.0 / (24 * 60)  # 2 minutes in days

        p1_diff = abs(jd_p1 - jd_p1_ref)
        p4_diff = abs(jd_p4 - jd_p4_ref)

        assert p1_diff < tolerance_days, (
            f"P1 timing error: {p1_diff * 24 * 60:.1f} minutes. Expected < 2 minutes"
        )
        assert p4_diff < tolerance_days, (
            f"P4 timing error: {p4_diff * 24 * 60:.1f} minutes. Expected < 2 minutes"
        )

    def test_may_2022_total_lunar_eclipse_penumbral_contacts(self):
        """Test May 16, 2022 total lunar eclipse penumbral contacts.

        NASA Reference:
        - P1 (Penumbral begins): 2022 May 16 at 01:32:05 UT
        - Maximum: 2022 May 16 at 04:11:28 UT
        - P4 (Penumbral ends): 2022 May 16 at 06:50:54 UT
        """
        jd_start = julday(2022, 4, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        # Verify we found the right eclipse (May 16, 2022)
        year, month, day, _ = revjul(jd_max)
        assert year == 2022
        assert month == 5
        assert day == 16

        jd_p1 = calc_lunar_eclipse_penumbral_first_contact_p1(jd_max)
        jd_p4 = calc_lunar_eclipse_penumbral_fourth_contact_p4(jd_max)

        # NASA reference: P1 at 01:32:05 UT
        jd_p1_ref = julday(2022, 5, 16, 1.535)  # 01:32:05 UT

        # NASA reference: P4 at 06:50:54 UT
        jd_p4_ref = julday(2022, 5, 16, 6.848)  # 06:50:54 UT

        # Allow 2 minutes tolerance
        tolerance_days = 2.0 / (24 * 60)

        p1_diff = abs(jd_p1 - jd_p1_ref)
        p4_diff = abs(jd_p4 - jd_p4_ref)

        assert p1_diff < tolerance_days, (
            f"P1 timing error: {p1_diff * 24 * 60:.1f} minutes"
        )
        assert p4_diff < tolerance_days, (
            f"P4 timing error: {p4_diff * 24 * 60:.1f} minutes"
        )


class TestPenumbralContactsWithLunEclipseWhen:
    """Test that P1 and P4 match the times returned by lun_eclipse_when."""

    def test_p1_matches_lun_eclipse_when_times6(self):
        """Test that P1 closely matches times[6] from lun_eclipse_when."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        # times[6] is penumbral eclipse beginning from lun_eclipse_when
        jd_pen_begin_from_when = times[6]
        jd_p1 = calc_lunar_eclipse_penumbral_first_contact_p1(jd_max)

        # Both should calculate the same event
        # Allow 2 minute tolerance (different algorithms may give slightly different results)
        tolerance_days = 2.0 / (24 * 60)
        diff = abs(jd_p1 - jd_pen_begin_from_when)

        assert diff < tolerance_days, (
            f"P1 and times[6] differ by {diff * 24 * 60:.1f} minutes. "
            f"Expected < 2 minutes"
        )

    def test_p4_matches_lun_eclipse_when_times7(self):
        """Test that P4 closely matches times[7] from lun_eclipse_when."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        # times[7] is penumbral eclipse ending from lun_eclipse_when
        jd_pen_end_from_when = times[7]
        jd_p4 = calc_lunar_eclipse_penumbral_fourth_contact_p4(jd_max)

        # Both should calculate the same event
        tolerance_days = 2.0 / (24 * 60)
        diff = abs(jd_p4 - jd_pen_end_from_when)

        assert diff < tolerance_days, (
            f"P4 and times[7] differ by {diff * 24 * 60:.1f} minutes. "
            f"Expected < 2 minutes"
        )


class TestPenumbralContactsPhysicalConstraints:
    """Test physical constraints on penumbral contact times."""

    def test_penumbral_duration_reasonable(self):
        """Test that P1 to P4 duration is physically reasonable."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_p1 = calc_lunar_eclipse_penumbral_first_contact_p1(jd_max)
        jd_p4 = calc_lunar_eclipse_penumbral_fourth_contact_p4(jd_max)

        # Total penumbral duration in hours
        duration_hours = (jd_p4 - jd_p1) * 24

        # Lunar eclipse penumbral phase typically lasts 4-6 hours
        assert 3.0 < duration_hours < 7.0, (
            f"Penumbral duration {duration_hours:.2f} hours seems unreasonable. "
            f"Expected 3-7 hours for a typical lunar eclipse"
        )

    def test_p1_to_max_reasonable(self):
        """Test that P1 to maximum duration is reasonable."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_p1 = calc_lunar_eclipse_penumbral_first_contact_p1(jd_max)

        # Duration from P1 to maximum in hours
        duration_hours = (jd_max - jd_p1) * 24

        # Typically 1.5-3.5 hours from P1 to maximum
        assert 1.5 < duration_hours < 3.5, (
            f"P1 to maximum: {duration_hours:.2f} hours seems unreasonable"
        )

    def test_max_to_p4_reasonable(self):
        """Test that maximum to P4 duration is reasonable."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_p4 = calc_lunar_eclipse_penumbral_fourth_contact_p4(jd_max)

        # Duration from maximum to P4 in hours
        duration_hours = (jd_p4 - jd_max) * 24

        # Typically 1.5-3.5 hours from maximum to P4
        assert 1.5 < duration_hours < 3.5, (
            f"Maximum to P4: {duration_hours:.2f} hours seems unreasonable"
        )


class TestPenumbralContactsDifferentEclipseTypes:
    """Test penumbral contacts for different eclipse types."""

    def test_partial_eclipse_has_penumbral_contacts(self):
        """Test that partial eclipses have valid P1 and P4."""
        jd_start = julday(2023, 10, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PARTIAL)
        jd_max = times[0]

        jd_p1 = calc_lunar_eclipse_penumbral_first_contact_p1(jd_max)
        jd_p4 = calc_lunar_eclipse_penumbral_fourth_contact_p4(jd_max)

        assert jd_p1 > 0, "P1 should exist for partial eclipse"
        assert jd_p4 > 0, "P4 should exist for partial eclipse"
        assert jd_p1 < jd_max < jd_p4, "P1 < max < P4 ordering should hold"

    def test_penumbral_only_eclipse_has_penumbral_contacts(self):
        """Test that penumbral-only eclipses have valid P1 and P4."""
        jd_start = julday(2020, 1, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PENUMBRAL)
        jd_max = times[0]

        # Only test if we found a pure penumbral eclipse
        if ecl_type == ECL_PENUMBRAL:
            jd_p1 = calc_lunar_eclipse_penumbral_first_contact_p1(jd_max)
            jd_p4 = calc_lunar_eclipse_penumbral_fourth_contact_p4(jd_max)

            assert jd_p1 > 0, "P1 should exist for penumbral eclipse"
            assert jd_p4 > 0, "P4 should exist for penumbral eclipse"
            assert jd_p1 < jd_max < jd_p4, "P1 < max < P4 ordering should hold"


class TestPenumbralContactsMultipleEclipses:
    """Test penumbral contacts for multiple sequential eclipses."""

    def test_multiple_eclipses_have_valid_contacts(self):
        """Test that P1 and P4 work correctly for multiple eclipses."""
        jd = julday(2022, 1, 1, 0.0)
        eclipses = []

        # Find 3 sequential lunar eclipses
        for _ in range(3):
            ecl_type, times = lun_eclipse_when(jd)
            jd_max = times[0]

            jd_p1 = calc_lunar_eclipse_penumbral_first_contact_p1(jd_max)
            jd_p4 = calc_lunar_eclipse_penumbral_fourth_contact_p4(jd_max)

            eclipses.append(
                {
                    "max": jd_max,
                    "p1": jd_p1,
                    "p4": jd_p4,
                }
            )

            jd = jd_max + 1  # Move past this eclipse

        # Verify each eclipse has valid penumbral contacts
        for i, ecl in enumerate(eclipses):
            assert ecl["p1"] > 0, f"Eclipse {i + 1}: P1 should be > 0"
            assert ecl["p4"] > 0, f"Eclipse {i + 1}: P4 should be > 0"
            assert ecl["p1"] < ecl["max"], f"Eclipse {i + 1}: P1 should be before max"
            assert ecl["p4"] > ecl["max"], f"Eclipse {i + 1}: P4 should be after max"
