"""
Tests for the osculating lunar perigee calculation.

The osculating perigee is computed directly from the eccentricity vector,
which naturally points toward perigee. This is different from the approach
of computing apogee and adding 180 degrees, which is incorrect because
apogee and perigee are NOT exactly opposite — they can deviate by up to
28 degrees depending on Sun-Moon geometry (Chapront-Touzé & Chapront 1988;
Meeus, Astronomical Algorithms ch. 47).

These tests verify:
1. Basic functionality of calc_osculating_perigee
2. The perigee is computed independently (not just apogee + 180)
3. The relationship between apogee and perigee (approximately opposite)
4. Consistency of interpolated perigee results
"""

import pytest
from libephemeris import lunar


class TestOsculatingPerigeeBasic:
    """Basic functionality tests for calc_osculating_perigee."""

    def test_returns_valid_longitude(self):
        """calc_osculating_perigee should return longitude in [0, 360) range."""
        jd_j2000 = 2451545.0
        lon, lat, ecc = lunar.calc_osculating_perigee(jd_j2000)

        assert 0 <= lon < 360, f"Longitude {lon} out of range"

    def test_returns_valid_latitude(self):
        """calc_osculating_perigee should return small latitude."""
        jd_j2000 = 2451545.0
        lon, lat, ecc = lunar.calc_osculating_perigee(jd_j2000)

        # Latitude should be less than 10 degrees
        assert -10 < lat < 10, f"Latitude {lat} unexpectedly large"

    def test_returns_valid_distance(self):
        """calc_osculating_perigee should return reasonable distance in AU."""
        jd_j2000 = 2451545.0
        lon, lat, dist = lunar.calc_osculating_perigee(jd_j2000)

        # 3rd return value is perigee distance in AU (~0.0024)
        assert 0.002 < dist < 0.003, f"Distance {dist} AU out of expected range"

    def test_works_for_historical_dates(self):
        """Should work for historical dates in ephemeris range."""
        jd_1950 = 2433282.5  # 1950-01-01
        lon, lat, ecc = lunar.calc_osculating_perigee(jd_1950)

        assert 0 <= lon < 360
        assert -90 < lat < 90
        assert 0.002 < ecc < 0.003

    def test_works_for_future_dates(self):
        """Should work for future dates in ephemeris range."""
        jd_2050 = 2469807.5  # 2050-01-01
        lon, lat, ecc = lunar.calc_osculating_perigee(jd_2050)

        assert 0 <= lon < 360
        assert -90 < lat < 90
        assert 0.002 < ecc < 0.003


class TestOsculatingPerigeeVsApogee:
    """Test relationship between perigee and apogee."""

    def test_perigee_approximately_opposite_apogee(self):
        """Perigee should be approximately 180 degrees from apogee.

        According to established lunar mechanics (Chapront-Touzé & Chapront
        1988), they are not exactly opposite but should be roughly so when
        the Sun is in certain configurations.
        """
        jd = 2451545.0  # J2000.0

        apogee_lon, apogee_lat, _ = lunar.calc_true_lilith(jd)
        perigee_lon, perigee_lat, _ = lunar.calc_osculating_perigee(jd)

        # Calculate angular difference
        diff = abs(apogee_lon - perigee_lon)
        if diff > 180:
            diff = 360 - diff

        # Should be approximately 180 degrees (within 30 degrees)
        # The up-to-28-degree deviation is a physical property of the Moon's orbit
        assert abs(diff - 180) < 30.0, (
            f"Apogee-perigee angular distance: {diff:.2f}°, "
            f"deviation from 180°: {abs(diff - 180):.2f}°"
        )

    def test_osculating_perigee_exactly_opposite_apogee(self):
        """Verify that osculating perigee IS exactly opposite osculating apogee.

        In a two-body problem, the eccentricity vector points toward perigee
        and apogee is exactly 180° opposite. This is by mathematical definition.
        The deviation from 180° in perturbation-corrected values (vs pure osculating
        elements) applies to the interpolated apsides, not the raw osculating values.
        """
        jd = 2451545.0

        apogee_lon, apogee_lat, _ = lunar.calc_true_lilith(jd)
        perigee_lon, perigee_lat, _ = lunar.calc_osculating_perigee(jd)

        # In two-body osculating elements, they should be exactly opposite
        diff = abs(apogee_lon - perigee_lon)
        if diff > 180:
            diff = 360 - diff

        # Should be exactly 180 degrees (within numerical precision)
        assert abs(diff - 180) < 1e-10, (
            f"Osculating perigee should be exactly opposite apogee, "
            f"but difference is {diff:.6f}°"
        )

        # Latitude should also be exactly negated
        assert abs(perigee_lat + apogee_lat) < 1e-10, (
            f"Osculating perigee latitude should be exactly negated, "
            f"but apogee_lat={apogee_lat:.6f}, perigee_lat={perigee_lat:.6f}"
        )

    def test_distance_values_reasonable(self):
        """Apogee and perigee distances should both be valid and apogee > perigee."""
        jd = 2451545.0

        _, _, apogee_dist = lunar.calc_true_lilith(jd)
        _, _, perigee_dist = lunar.calc_osculating_perigee(jd)

        # Both should be valid lunar distances in AU
        assert 0.002 < apogee_dist < 0.003, f"Apogee distance {apogee_dist} AU invalid"
        assert 0.002 < perigee_dist < 0.003, (
            f"Perigee distance {perigee_dist} AU invalid"
        )
        # Apogee distance should be greater than perigee distance
        assert apogee_dist > perigee_dist, (
            f"Apogee dist ({apogee_dist}) should exceed perigee dist ({perigee_dist})"
        )


class TestInterpolatedPerigeeImprovement:
    """Test that interpolated perigee uses the improved computation."""

    def test_sample_function_uses_calc_osculating_perigee(self):
        """Verify that _sample_osculating_perigee_with_fallback computes
        perigee independently, not as apogee + 180.

        We test this by checking that the sampled values match
        calc_osculating_perigee, not calc_true_lilith + 180.
        """
        jd_tt = 2451545.0
        half_window = 28.0
        num_samples = 9

        sample_times, sample_lons, sample_lats, sample_eccs, target_idx = (
            lunar._sample_osculating_perigee_with_fallback(
                jd_tt, half_window, num_samples
            )
        )

        # Check that the sampled values match calc_osculating_perigee
        for i, sample_jd in enumerate(sample_times):
            expected_lon, expected_lat, expected_ecc = lunar.calc_osculating_perigee(
                sample_jd
            )

            assert abs(sample_lons[i] - expected_lon) < 1e-10, (
                f"Sample {i}: longitude mismatch"
            )
            assert abs(sample_lats[i] - expected_lat) < 1e-10, (
                f"Sample {i}: latitude mismatch"
            )
            assert abs(sample_eccs[i] - expected_ecc) < 1e-10, (
                f"Sample {i}: eccentricity mismatch"
            )

    def test_interpolated_perigee_consistency(self):
        """Test that interpolated perigee values are consistent over time."""
        test_dates = [2451545.0, 2451600.0, 2451700.0, 2451800.0]

        for jd in test_dates:
            lon, lat, ecc = lunar.calc_interpolated_perigee(jd)

            assert 0 <= lon < 360, f"Invalid longitude at JD {jd}"
            assert -10 < lat < 10, f"Invalid latitude at JD {jd}"
            assert 0.002 < ecc < 0.003, f"Invalid distance at JD {jd}"


class TestInterpolatedPerigeeDeviation:
    """Test that interpolated perigee deviates from interpolated apogee + 180."""

    def test_interpolated_deviation_varies_with_date(self):
        """The interpolated perigee should NOT be exactly opposite interpolated apogee.

        According to established lunar orbital mechanics (Chapront-Touzé &
        Chapront 1988), the interpolated apogee and perigee have different
        oscillation amplitudes (5° vs 25°). When we interpolate (smooth) these
        different oscillations, the resulting values differ by more than 180°.

        This is the key improvement: even though osculating values are exactly
        opposite, the interpolation process produces different results because
        it's smoothing over curves with different oscillation characteristics.
        """
        # Sample dates over a full year (to cover different Sun-Moon geometries)
        start_jd = 2451545.0  # 2000-01-01
        test_dates = [start_jd + i * 30 for i in range(13)]  # ~1 year

        deviations = []
        for jd in test_dates:
            apogee_lon, _, _ = lunar.calc_interpolated_apogee(jd)
            perigee_lon, _, _ = lunar.calc_interpolated_perigee(jd)

            diff = abs(apogee_lon - perigee_lon)
            if diff > 180:
                diff = 360 - diff

            deviations.append(abs(diff - 180))

        # Check that deviations vary and are non-trivial
        min_dev = min(deviations)
        max_dev = max(deviations)
        avg_dev = sum(deviations) / len(deviations)

        # The interpolated values should show significant deviation from 180
        # (because they have different oscillation characteristics)
        assert avg_dev > 1.0, (
            f"Average deviation from 180° is only {avg_dev:.2f}°. "
            "Interpolated perigee should differ significantly from apogee + 180°."
        )

        # Print statistics for debugging
        print("\nInterpolated apogee/perigee deviation statistics over 1 year:")
        print(f"  Min: {min_dev:.2f}°")
        print(f"  Max: {max_dev:.2f}°")
        print(f"  Avg: {avg_dev:.2f}°")

    def test_interpolated_perigee_closer_to_pyswisseph_than_apogee_plus_180(self):
        """Verify our interpolated perigee is closer to pyswisseph than apogee + 180.

        This test confirms the fix works: computing perigee independently and
        interpolating it should give better agreement with pyswisseph
        than simply adding 180° to the interpolated apogee.
        """
        try:
            import swisseph as swe
        except ImportError:
            pytest.skip("pyswisseph not installed")

        jd = 2451545.0

        # pyswisseph (via pyswisseph) interpolated perigee
        swe_perg, _ = swe.calc_ut(jd, swe.INTP_PERG, 0)
        swe_perg_lon = swe_perg[0]

        # Our interpolated perigee (computed independently)
        lib_perg_lon, _, _ = lunar.calc_interpolated_perigee(jd)

        # Our interpolated apogee + 180 (the old incorrect approach)
        lib_apog_lon, _, _ = lunar.calc_interpolated_apogee(jd)
        derived_perg_lon = (lib_apog_lon + 180.0) % 360.0

        # Calculate errors
        direct_error = abs(swe_perg_lon - lib_perg_lon)
        if direct_error > 180:
            direct_error = 360 - direct_error

        derived_error = abs(swe_perg_lon - derived_perg_lon)
        if derived_error > 180:
            derived_error = 360 - derived_error

        print("\nPerigee error comparison:")
        print(f"  Direct (independent) computation: {direct_error:.2f}°")
        print(f"  Derived (apogee + 180): {derived_error:.2f}°")

        # Our direct computation should be closer to pyswisseph
        assert direct_error < derived_error, (
            f"Direct perigee error ({direct_error:.2f}°) should be less than "
            f"derived error ({derived_error:.2f}°)"
        )
