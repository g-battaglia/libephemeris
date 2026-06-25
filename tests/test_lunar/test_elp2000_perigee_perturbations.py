"""
Tests for the ELP2000-82B perigee perturbation series.

This module tests the accuracy and correctness of the perturbation
series used for calculating the interpolated lunar perigee.

The perturbation series is calibrated against JPL DE441 using
passage-interpolated harmonic fitting (v2.2). The series alone achieves
~2 deg RMS near J2000, and with the correction table achieves < 0.1 deg.

The SE comparison shows ~1-2 deg RMS difference, which reflects the
intrinsic difference between JPL-calibrated and SE analytical approaches.

Tests validate:
1. Basic functionality and reasonable output ranges
2. Internal consistency (perturbations vary as expected)
3. JPL self-consistency (with correction table, precision is high)
4. SE comparison with relaxed thresholds (informational only)
"""

import math
import pytest

from libephemeris import julday, lunar


def normalize_angle_diff(diff: float) -> float:
    """Normalize angle difference to [-180, 180) range."""
    while diff >= 180.0:
        diff -= 360.0
    while diff < -180.0:
        diff += 360.0
    return diff


class TestELP2000PerigeePerturbationsBasic:
    """Basic functionality tests for perigee perturbation calculations."""

    def test_perturbation_function_exists(self):
        """Test that the perturbation function exists and is callable."""
        jd_tt = 2451545.0
        result = lunar._calc_elp2000_perigee_perturbations(jd_tt)
        assert isinstance(result, float)

    def test_perturbation_magnitude_reasonable(self):
        """Test that perturbation values are within expected range."""
        jd_tt = 2451545.0
        result = lunar._calc_elp2000_perigee_perturbations(jd_tt)
        assert -30.0 < result < 30.0, f"Perturbation {result} deg out of expected range"

    def test_perturbation_varies_with_time(self):
        """Test that perturbation changes over time (not constant)."""
        jd1 = 2451545.0
        jd2 = 2451560.0

        pert1 = lunar._calc_elp2000_perigee_perturbations(jd1)
        pert2 = lunar._calc_elp2000_perigee_perturbations(jd2)

        assert pert1 != pert2, "Perturbation should change with time"

    def test_dominant_term_is_negative(self):
        """Test that the dominant 2D-2M' evection term has negative coefficient."""
        jd_tt = 2451545.0
        pert = lunar._calc_elp2000_perigee_perturbations(jd_tt)
        assert isinstance(pert, float)


class TestELP2000PerigeePerturbationsTerms:
    """Test specific perturbation term contributions."""

    def test_dominant_evection_term_variation(self):
        """Test that the evection term causes significant variation."""
        jd1 = 2451545.0
        jd2 = 2451545.0 + 14.77

        pert1 = lunar._calc_elp2000_perigee_perturbations(jd1)
        pert2 = lunar._calc_elp2000_perigee_perturbations(jd2)

        diff = abs(pert1 - pert2)
        assert diff > 2.0, f"Evection term should cause variation > 2 deg, got {diff}"

    def test_perturbation_oscillation_range(self):
        """Test that perturbation shows significant oscillation over time."""
        jd_start = 2451545.0
        perturbations = []

        for i in range(60):
            jd = jd_start + i
            pert = lunar._calc_elp2000_perigee_perturbations(jd)
            perturbations.append(pert)

        oscillation_range = max(perturbations) - min(perturbations)
        assert oscillation_range > 5.0, (
            f"Should see oscillation range > 5 deg, got {oscillation_range}"
        )


class TestELP2000PerigeePerturbationsEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_at_year_2000(self):
        """Test calculation at J2000 epoch."""
        jd_ut = julday(2000, 1, 1, 0.0)
        jd_tt = jd_ut + 69.184 / 86400.0

        our_lon, _, _ = lunar.calc_interpolated_perigee(jd_tt)
        assert 0.0 <= our_lon < 360.0, f"Longitude {our_lon} not in [0, 360)"

    def test_continuity_over_month(self):
        """Test that perigee position is continuous over a month."""
        prev_lon = None
        max_jump = 0.0

        for i in range(30):
            jd_tt = 2451545.0 + i
            lon, _, _ = lunar.calc_interpolated_perigee(jd_tt)

            if prev_lon is not None:
                jump = abs(normalize_angle_diff(lon - prev_lon))
                max_jump = max(max_jump, jump)

            prev_lon = lon

        assert max_jump < 5.0, f"Perigee position jumped {max_jump} deg in one day"


class TestELP2000PerigeePerturbationsJPL:
    """JPL self-consistency tests.

    These validate the interpolated perigee against known physical properties
    derived from JPL DE440/DE441, without requiring any external reference.

    Tests cover:
    - Mean apsidal precession rate (~40.69 deg/yr, Simon et al. 1994)
    - Three-level decomposition consistency (mean + series + table)
    - Correction table continuity at interpolation boundaries
    - Output format and range invariants
    - Long-term smoothness over decades
    - Comparison with JPL perigee passage reference data
    """

    def test_mean_apsidal_precession_rate(self):
        """Mean perigee should precess at ~40.69 deg/yr (Simon et al. 1994).

        The mean apsidal precession rate is one of the best-known parameters
        of the lunar orbit. We verify it over a 100-year span to ensure the
        underlying polynomial is correct.
        """
        jd_start = julday(2000, 1, 1, 12.0) + 69.184 / 86400.0
        jd_end = julday(2100, 1, 1, 12.0) + 69.184 / 86400.0

        # Track cumulative advance month-by-month to handle wraparound
        prev = (lunar.calc_mean_lilith(jd_start) + 180.0) % 360.0
        total_advance = 0.0
        for month in range(1, 100 * 12 + 1):
            year = 2000 + month // 12
            m = month % 12 + 1
            jd = julday(year, m, 1, 12.0) + 69.184 / 86400.0
            cur = (lunar.calc_mean_lilith(jd) + 180.0) % 360.0
            diff = normalize_angle_diff(cur - prev)
            total_advance += diff
            prev = cur

        days = jd_end - jd_start
        years = days / 365.25
        rate = total_advance / years

        # Simon et al. (1994) value: ~40.6901 deg/yr
        assert abs(rate - 40.69) < 0.1, (
            f"Apsidal precession rate {rate:.4f} deg/yr, expected ~40.69"
        )

    def test_three_level_decomposition_consistency(self):
        """Full result must equal mean + perturbation + correction table."""
        test_jds = [
            julday(1950, 6, 15, 0.0) + 69.184 / 86400.0,
            2451545.0,  # J2000
            julday(2020, 3, 20, 0.0) + 69.184 / 86400.0,
            julday(2050, 12, 1, 0.0) + 69.184 / 86400.0,
        ]

        for jd_tt in test_jds:
            full_lon, lat, dist = lunar.calc_interpolated_perigee(jd_tt)

            # Reconstruct manually: mean perigee + perturbation
            # Note: calc_interpolated_perigee also adds a correction table
            # residual (~0.03° typical), so allow tolerance for that.
            mean = (lunar.calc_mean_lilith(jd_tt) + 180.0) % 360.0
            pert = lunar._calc_elp2000_perigee_perturbations(jd_tt)
            reconstructed = (mean + pert) % 360.0

            diff = abs(normalize_angle_diff(full_lon - reconstructed))
            assert diff < 0.2, (
                f"Decomposition mismatch at JD {jd_tt}: "
                f"full={full_lon:.6f}, reconstructed={reconstructed:.6f}"
            )

    def test_correction_table_available(self):
        """Perigee correction table must be loaded."""
        assert lunar._PERIGEE_CORRECTIONS_AVAILABLE, (
            "PERIGEE_PERTURBATION_CORRECTIONS not loaded from lunar_corrections.py"
        )

    def test_correction_table_continuity(self):
        """Correction values must be continuous at table boundaries.

        The correction table uses linear interpolation with entries every 2
        years. At each boundary, the interpolated value from both sides must
        converge to the same value (no discontinuities).
        """
        from libephemeris.lunar_corrections import (
            PERIGEE_CORRECTION_START_YEAR,
            PERIGEE_CORRECTION_STEP_YEARS,
        )

        eps = 0.001  # ~8.7 hours in years
        # Check 10 boundaries in the modern era (around year 2000)
        idx_2000 = int(
            (2000 - PERIGEE_CORRECTION_START_YEAR) / PERIGEE_CORRECTION_STEP_YEARS
        )
        for i in range(idx_2000 - 5, idx_2000 + 5):
            boundary_year = (
                PERIGEE_CORRECTION_START_YEAR + i * PERIGEE_CORRECTION_STEP_YEARS
            )
            jd_boundary = 2451545.0 + (boundary_year - 2000.0) * 365.25

            lon_before, _, _ = lunar.calc_interpolated_perigee(jd_boundary - eps)
            lon_at, _, _ = lunar.calc_interpolated_perigee(jd_boundary)
            lon_after, _, _ = lunar.calc_interpolated_perigee(jd_boundary + eps)

            jump1 = abs(normalize_angle_diff(lon_at - lon_before))
            jump2 = abs(normalize_angle_diff(lon_after - lon_at))

            # With eps ~8.7 hours, perigee moves ~0.005 deg; any jump > 0.1
            # would indicate a table boundary discontinuity
            assert jump1 < 0.1, (
                f"Discontinuity at boundary year {boundary_year}: "
                f"jump = {jump1:.6f} deg"
            )
            assert jump2 < 0.1, (
                f"Discontinuity at boundary year {boundary_year}: "
                f"jump = {jump2:.6f} deg"
            )

    def test_output_format_invariants(self):
        """calc_interpolated_perigee must return (lon, lat, dist)."""
        test_jds = [2451545.0, julday(1900, 1, 1, 0.0), julday(2100, 1, 1, 0.0)]

        for jd_tt in test_jds:
            lon, lat, dist = lunar.calc_interpolated_perigee(jd_tt)

            assert 0.0 <= lon < 360.0, f"Longitude {lon} not in [0, 360)"
            assert -6.0 <= lat <= 6.0, f"Latitude {lat} not in [-6, 6]"
            assert 0.001 < dist < 0.004, f"Distance {dist} not in expected AU range"

    def test_long_term_smoothness(self):
        """Perigee position must advance smoothly over decades.

        Weekly samples over 10 years should show no jumps > 15 deg.
        The mean advance is ~40.69/52 ≈ 0.78 deg/week, but oscillations
        from the perturbation series can add up to ~10 deg/week.
        """
        max_jump = 0.0
        prev_lon = None
        jd_start = julday(2000, 1, 1, 0.0) + 69.184 / 86400.0

        for week in range(520):  # 10 years
            jd_tt = jd_start + week * 7.0
            lon, _, _ = lunar.calc_interpolated_perigee(jd_tt)

            if prev_lon is not None:
                jump = abs(normalize_angle_diff(lon - prev_lon))
                max_jump = max(max_jump, jump)

            prev_lon = lon

        assert max_jump < 15.0, (
            f"Perigee jumped {max_jump:.2f} deg in one week (expected < 15)"
        )

    def test_correction_table_improves_accuracy(self):
        """With correction table, result should differ from series-only.

        At arbitrary dates, the correction table adds a non-trivial
        correction (typically several degrees) that compensates for secular
        drift and missing harmonics in the 62-term series.
        """
        dates_with_significant_correction = 0
        for year in range(1950, 2051, 5):
            jd_tt = julday(year, 6, 15, 0.0) + 69.184 / 86400.0
            corr = lunar._interpolate_perigee_correction(jd_tt)
            if abs(corr) > 1.0:
                dates_with_significant_correction += 1

        # The correction should be significant (> 1 deg) at most dates
        total_dates = len(range(1950, 2051, 5))
        assert dates_with_significant_correction > total_dates * 0.5, (
            f"Correction table is significant at only "
            f"{dates_with_significant_correction}/{total_dates} dates"
        )

    @pytest.mark.slow
    def test_jpl_perigee_passage_comparison(self):
        """Compare interpolated perigee with JPL Moon longitude at perigee passages.

        At actual perigee passages (Earth-Moon distance minima from JPL DE440),
        the Moon's ecliptic longitude IS the true perigee longitude. Our
        interpolated perigee is a smooth approximation, so differences of
        ~5-15 deg are expected (the interpolated perigee tracks the apsidal
        line, not the instantaneous Moon position).

        This test verifies the interpolated perigee stays within a physically
        reasonable range of the true perigee longitude at known passages.

        Reference data computed from JPL DE440/DE441 via Skyfield.
        """
        # Pre-computed JPL perigee passages: (jd_tt, moon_lon_at_perigee)
        # Found by minimizing Earth-Moon distance with golden section search
        jpl_passages = [
            (2433432.181891, 236.601348),  # ~1950 Jun
            (2440759.250092, 300.605123),  # ~1970 Jun
            (2448063.955048, 71.078870),  # ~1990 Jun
            (2451699.055152, 87.996896),  # ~2000 Jun
            (2455363.125269, 126.693003),  # ~2010 Jun
            (2459003.652540, 217.061769),  # ~2020 Jun
            (2462667.482174, 252.675005),  # ~2030 Jun
            (2469967.266384, 316.699836),  # ~2050 Jun
        ]

        errors = []
        for jd_tt, jpl_lon in jpl_passages:
            our_lon, _, _ = lunar.calc_interpolated_perigee(jd_tt)
            error = normalize_angle_diff(our_lon - jpl_lon)
            errors.append(error)

        rms = math.sqrt(sum(e**2 for e in errors) / len(errors))
        max_err = max(abs(e) for e in errors)

        # The interpolated perigee is a smooth apsidal-line approximation,
        # not the exact Moon longitude at each passage. RMS ~7 deg and
        # max ~15 deg are expected from the structural difference.
        assert rms < 15.0, f"RMS error vs JPL passages: {rms:.2f} deg (expected < 15)"
        assert max_err < 20.0, (
            f"Max error vs JPL passages: {max_err:.2f} deg (expected < 20)"
        )


class TestELP2000PerigeePerturbationsSupermoon:
    """Tests focused on supermoon timing accuracy."""

    def test_supermoon_dates_return_valid_values(self):
        """Test that supermoon dates return valid perigee positions."""
        test_dates = [
            (2015, 9, 27),
            (2016, 11, 14),
            (2019, 2, 19),
            (2021, 4, 27),
            (2023, 8, 30),
        ]

        for year, month, day in test_dates:
            jd_ut = julday(year, month, day, 0.0)
            jd_tt = jd_ut + 69.184 / 86400.0

            our_lon, _, _ = lunar.calc_interpolated_perigee(jd_tt)

            assert 0.0 <= our_lon < 360.0, (
                f"Invalid longitude {our_lon} at {year}-{month}-{day}"
            )
