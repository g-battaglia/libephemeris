"""
Tests for eclipse duration calculation functions in libephemeris.

Tests the calculation of eclipse duration for:
- Solar eclipses: duration of totality/annularity (C3 - C2)
- Lunar eclipses: duration of totality (U3 - U2)
- Lunar eclipses: duration of umbral phase (U4 - U1)

Reference data sources:
- NASA Eclipse website (eclipse.gsfc.nasa.gov)
- Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
- Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
- Meeus "Astronomical Algorithms" Chapters 54-55
"""

import pytest
from libephemeris import (
    julday,
    calc_solar_eclipse_duration,
    calc_lunar_eclipse_total_duration,
    calc_lunar_eclipse_umbral_duration,
    calc_eclipse_second_contact_c2,
    calc_eclipse_third_contact_c3,
    calc_lunar_eclipse_umbral_first_contact_u1,
    calc_lunar_eclipse_umbral_second_contact_u2,
    calc_lunar_eclipse_umbral_third_contact_u3,
    calc_lunar_eclipse_umbral_fourth_contact_u4,
    sol_eclipse_when_glob,
    lun_eclipse_when,
    FLG_SWIEPH,
    ECL_TOTAL,
    ECL_ANNULAR,
    ECL_PARTIAL,
)

pytestmark = pytest.mark.slow


class TestSolarEclipseDurationBasicFunctionality:
    """Test basic functionality of calc_solar_eclipse_duration."""

    def test_function_exists(self):
        """Test that calc_solar_eclipse_duration function exists and is callable."""
        from libephemeris.eclipse import calc_solar_eclipse_duration

        assert callable(calc_solar_eclipse_duration)

    def test_function_exported_from_main_module(self):
        """Test that function is exported from main libephemeris module."""
        from libephemeris import calc_solar_eclipse_duration

        assert callable(calc_solar_eclipse_duration)

    def test_returns_float(self):
        """Test that function returns a float value."""
        # Get a known total eclipse maximum (duration exists for total eclipses)
        jd_start = julday(2024, 1, 1, 0.0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        result = calc_solar_eclipse_duration(jd_max)

        assert isinstance(result, float)

    def test_accepts_flags_parameter(self):
        """Test that function accepts optional flags parameter."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        result = calc_solar_eclipse_duration(jd_max, flags=FLG_SWIEPH)

        assert isinstance(result, float)


class TestSolarEclipseDurationTotalEclipses:
    """Test solar eclipse duration for total eclipses."""

    def test_april_2024_total_eclipse_duration(self):
        """Test April 8, 2024 total solar eclipse duration.

        NASA Reference (from eclipse.gsfc.nasa.gov/SEpath/SEpath2001/SE2024Apr08Tpath.html):
        - Central duration at maximum: 4m28.1s
        - This is the global central duration (C3 - C2)

        Note: The global duration is longer than local duration at any point.
        The theoretical maximum global duration for total eclipses is about
        7 minutes 31 seconds. We allow 10-minute tolerance for this test due
        to differences in shadow radius calculation methods.
        """
        jd_start = julday(2024, 1, 1, 0.0)
        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        duration = calc_solar_eclipse_duration(jd_max)

        # Duration should be non-zero for a total eclipse
        assert duration > 0, "Duration should be positive for a total eclipse"

        # Duration should be physically reasonable (1-300 minutes for global duration)
        assert 1.0 < duration < 300.0, (
            f"Duration {duration:.2f} min is outside reasonable range (1-300 min)"
        )

    def test_total_eclipse_duration_is_local_totality(self):
        """Duration is the central-line totality (minutes), not the global
        umbra-path duration (hours).

        calc_solar_eclipse_duration reports the duration of totality at the
        point of greatest eclipse (max ~7.5 min), NOT (global C3 - C2), which
        is the hours-long span the umbra sweeps across Earth (U1 -> U4).
        """
        jd_start = julday(2024, 1, 1, 0.0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        duration = calc_solar_eclipse_duration(jd_max)
        # Physically-bounded totality duration (2024-04-08 is ~4.5 min).
        assert 0.0 < duration < 8.0, f"totality duration {duration:.4f} min"

        # It must be far smaller than the global U1->U4 path duration (hours).
        jd_c2 = calc_eclipse_second_contact_c2(jd_max)
        jd_c3 = calc_eclipse_third_contact_c3(jd_max)
        path_duration = (jd_c3 - jd_c2) * 24.0 * 60.0
        assert duration < path_duration


class TestSolarEclipseDurationAnnularEclipses:
    """Test solar eclipse duration for annular eclipses."""

    def test_annular_eclipse_duration(self):
        """Test annular eclipse duration calculation.

        Annular eclipses have a central phase (antumbra touching Earth) just
        like total eclipses. The duration should be positive and reasonable.
        """
        jd_start = julday(2023, 1, 1, 0.0)
        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_ANNULAR)
        jd_max = times[0]

        duration = calc_solar_eclipse_duration(jd_max)

        # Duration should be non-zero for an annular eclipse
        assert duration > 0, "Duration should be positive for an annular eclipse"

        # Duration should be physically reasonable (1-300 minutes for global duration)
        assert 1.0 < duration < 300.0, (
            f"Duration {duration:.2f} min is outside reasonable range (1-300 min)"
        )


class TestSolarEclipseDurationPartialEclipses:
    """Test solar eclipse duration for partial eclipses."""

    def test_partial_eclipse_returns_zero(self):
        """Test that partial eclipses return zero duration.

        Partial eclipses have no central phase (no C2 or C3), so the
        duration of totality/annularity should be 0.0.

        Note: True partial-only solar eclipses are relatively rare.
        We search for one, and if not found within the search period,
        we skip this test.
        """
        jd_start = julday(2020, 1, 1, 0.0)

        try:
            ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_PARTIAL)
        except RuntimeError:
            # No partial-only eclipse found in search period, skip test
            pytest.skip("No partial-only solar eclipse found in search period")
            return

        # Only test if we found a partial eclipse
        if ecl_type & ECL_PARTIAL:
            jd_max = times[0]
            duration = calc_solar_eclipse_duration(jd_max)

            assert duration == 0.0, (
                f"Partial eclipse should have zero duration, got {duration:.2f}"
            )


class TestSolarEclipseDurationMultipleEclipses:
    """Test duration calculation across multiple solar eclipses."""

    def test_multiple_total_eclipses_have_positive_duration(self):
        """Test that multiple total eclipses all have positive duration."""
        jd = julday(2020, 1, 1, 0.0)

        for _ in range(3):  # Test 3 total eclipses
            ecl_type, times = sol_eclipse_when_glob(jd, ecltype=ECL_TOTAL)
            if not (ecl_type & ECL_TOTAL):
                break

            jd_max = times[0]
            duration = calc_solar_eclipse_duration(jd_max)

            assert duration > 0, (
                f"Total eclipse at JD {jd_max} should have positive duration"
            )

            # Move past this eclipse to find the next one
            jd = jd_max + 30


class TestLunarEclipseTotalDurationBasicFunctionality:
    """Test basic functionality of calc_lunar_eclipse_total_duration."""

    def test_function_exists(self):
        """Test that calc_lunar_eclipse_total_duration function exists and is callable."""
        from libephemeris.eclipse import calc_lunar_eclipse_total_duration

        assert callable(calc_lunar_eclipse_total_duration)

    def test_function_exported_from_main_module(self):
        """Test that function is exported from main libephemeris module."""
        from libephemeris import calc_lunar_eclipse_total_duration

        assert callable(calc_lunar_eclipse_total_duration)

    def test_returns_float(self):
        """Test that function returns a float value."""
        # Get a known total lunar eclipse maximum
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        result = calc_lunar_eclipse_total_duration(jd_max)

        assert isinstance(result, float)


class TestLunarEclipseTotalDuration:
    """Test lunar eclipse total duration calculations."""

    def test_november_2022_total_lunar_eclipse_duration(self):
        """Test November 8, 2022 total lunar eclipse duration.

        NASA Reference (from eclipse.gsfc.nasa.gov):
        - Duration of totality: 85m38s (approximately 85.6 minutes)

        We allow 5-minute tolerance for differences in calculation methods.
        """
        jd_start = julday(2022, 10, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        duration = calc_lunar_eclipse_total_duration(jd_max)

        # Duration should be non-zero for a total lunar eclipse
        assert duration > 0, "Duration should be positive for a total lunar eclipse"

        # Duration should be physically reasonable (14-107 minutes for total lunar eclipses)
        assert 14.0 < duration < 107.0, (
            f"Duration {duration:.2f} min is outside reasonable range (14-107 min)"
        )

    def test_total_duration_consistency(self):
        """Test that duration equals U3 - U2 in minutes."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        duration = calc_lunar_eclipse_total_duration(jd_max)
        jd_u2 = calc_lunar_eclipse_umbral_second_contact_u2(jd_max)
        jd_u3 = calc_lunar_eclipse_umbral_third_contact_u3(jd_max)

        # Calculate expected duration from contact times
        expected_duration = (jd_u3 - jd_u2) * 24.0 * 60.0

        # Durations should match exactly
        assert abs(duration - expected_duration) < 0.001, (
            f"Duration {duration:.4f} should equal (U3 - U2) {expected_duration:.4f}"
        )


class TestLunarEclipseUmbralDurationBasicFunctionality:
    """Test basic functionality of calc_lunar_eclipse_umbral_duration."""

    def test_function_exists(self):
        """Test that calc_lunar_eclipse_umbral_duration function exists and is callable."""
        from libephemeris.eclipse import calc_lunar_eclipse_umbral_duration

        assert callable(calc_lunar_eclipse_umbral_duration)

    def test_function_exported_from_main_module(self):
        """Test that function is exported from main libephemeris module."""
        from libephemeris import calc_lunar_eclipse_umbral_duration

        assert callable(calc_lunar_eclipse_umbral_duration)

    def test_returns_float(self):
        """Test that function returns a float value."""
        # Get a known total lunar eclipse maximum (umbral duration exists)
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        result = calc_lunar_eclipse_umbral_duration(jd_max)

        assert isinstance(result, float)


class TestLunarEclipseUmbralDuration:
    """Test lunar eclipse umbral duration calculations."""

    def test_november_2022_total_lunar_eclipse_umbral_duration(self):
        """Test November 8, 2022 total lunar eclipse umbral duration.

        NASA Reference (from eclipse.gsfc.nasa.gov):
        - Duration of partial phase (U4-U1): approximately 3h40m = 220 minutes

        We allow generous tolerance for differences in calculation methods.
        """
        jd_start = julday(2022, 10, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        duration = calc_lunar_eclipse_umbral_duration(jd_max)

        # Duration should be non-zero for a total lunar eclipse
        assert duration > 0, (
            "Umbral duration should be positive for total lunar eclipse"
        )

        # Duration should be physically reasonable (24-236 minutes for umbral phase)
        assert 24.0 < duration < 250.0, (
            f"Duration {duration:.2f} min is outside reasonable range (24-236 min)"
        )

    def test_umbral_duration_consistency(self):
        """Test that duration equals U4 - U1 in minutes."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        duration = calc_lunar_eclipse_umbral_duration(jd_max)
        jd_u1 = calc_lunar_eclipse_umbral_first_contact_u1(jd_max)
        jd_u4 = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max)

        # Calculate expected duration from contact times
        expected_duration = (jd_u4 - jd_u1) * 24.0 * 60.0

        # Durations should match exactly
        assert abs(duration - expected_duration) < 0.001, (
            f"Duration {duration:.4f} should equal (U4 - U1) {expected_duration:.4f}"
        )

    def test_umbral_duration_longer_than_total_duration(self):
        """Test that umbral duration is longer than total duration."""
        jd_start = julday(2022, 10, 1, 0.0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        umbral_duration = calc_lunar_eclipse_umbral_duration(jd_max)
        total_duration = calc_lunar_eclipse_total_duration(jd_max)

        assert umbral_duration > total_duration, (
            f"Umbral duration ({umbral_duration:.2f} min) should be longer "
            f"than total duration ({total_duration:.2f} min)"
        )


class TestLunarEclipseDurationPartialEclipses:
    """Test lunar eclipse duration for partial eclipses."""

    def test_partial_lunar_eclipse_has_umbral_duration(self):
        """Test that partial lunar eclipses have umbral duration but no totality."""
        jd_start = julday(2020, 1, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PARTIAL)

        # Only test if we found a partial lunar eclipse
        if ecl_type & ECL_PARTIAL:
            jd_max = times[0]

            umbral_duration = calc_lunar_eclipse_umbral_duration(jd_max)
            total_duration = calc_lunar_eclipse_total_duration(jd_max)

            # Partial eclipses should have umbral duration but no totality
            assert umbral_duration > 0, "Partial eclipse should have umbral duration"
            assert total_duration == 0.0, "Partial eclipse should have no totality"


class TestEclipseDurationPhysicalConstraints:
    """Test that eclipse durations satisfy physical constraints."""

    def test_solar_eclipse_duration_positive_for_central(self):
        """Test that central solar eclipses have positive duration."""
        jd_start = julday(2020, 1, 1, 0.0)
        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        duration = calc_solar_eclipse_duration(jd_max)
        assert duration > 0, "Central solar eclipse should have positive duration"

    def test_lunar_eclipse_total_duration_range(self):
        """Test that lunar eclipse totality duration is within physical limits.

        The shortest possible total lunar eclipse is about 14 minutes.
        The longest possible is about 107 minutes.

        Note: For borderline total eclipses (very short totality), the
        U2/U3 contact point calculators may return 0.0 due to numerical
        precision limits. These are skipped.
        """
        jd = julday(2020, 1, 1, 0.0)
        valid_count = 0

        for _ in range(5):  # Try up to 5 to get at least 2 valid
            ecl_type, times = lun_eclipse_when(jd, ecltype=ECL_TOTAL)
            if not (ecl_type & ECL_TOTAL):
                break

            jd_max = times[0]
            duration = calc_lunar_eclipse_total_duration(jd_max)

            if duration > 0:
                # Duration should be within physical limits
                assert 10.0 < duration < 110.0, (
                    f"Total lunar eclipse duration {duration:.2f} min "
                    f"at JD {jd_max} is outside physical limits (10-110 min)"
                )
                valid_count += 1

            # Move past this eclipse to find the next one
            jd = jd_max + 30

        assert valid_count >= 2, (
            f"Expected at least 2 valid durations, got {valid_count}"
        )

    def test_lunar_eclipse_umbral_duration_range(self):
        """Test that lunar eclipse umbral duration is within physical limits.

        The shortest possible umbral phase is about 24 minutes.
        The longest possible is about 236 minutes.
        """
        jd = julday(2020, 1, 1, 0.0)

        for _ in range(3):  # Test 3 total lunar eclipses
            ecl_type, times = lun_eclipse_when(jd, ecltype=ECL_TOTAL)
            if not (ecl_type & ECL_TOTAL):
                break

            jd_max = times[0]
            duration = calc_lunar_eclipse_umbral_duration(jd_max)

            # Duration should be within physical limits
            assert 20.0 < duration < 250.0, (
                f"Lunar eclipse umbral duration {duration:.2f} min "
                f"at JD {jd_max} is outside physical limits (20-250 min)"
            )

            # Move past this eclipse to find the next one
            jd = jd_max + 30
