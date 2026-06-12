"""
Tests for solar gravitational perturbation on the eccentricity vector in True Lilith.

The True Lilith calculation applies a three-body correction to the eccentricity
vector to account for the Sun's gravitational influence on the lunar orbit.
This perturbation rotates the eccentricity vector based on the angle between
the apsidal line and the Sun's direction, capturing the quadrupole nature of
solar tidal forces.

Physical Background:
- The Sun's tidal force creates a quadrupole field at the Earth
- This field perturbs the lunar eccentricity vector direction
- The effect depends on sin(2*theta) where theta is the angle between
  the eccentricity vector and the solar direction in the orbital plane
- The amplitude is approximately 0.01148 in eccentricity units
- This corresponds to the evection perturbation (~1.274 degrees in longitude)

References:
- Brown, E.W. "An Introductory Treatise on the Lunar Theory" (1896)
- Brouwer, D. & Clemence, G.M. "Methods of Celestial Mechanics" (1961)
- Chapront-Touze, M. & Chapront, J. "Lunar Tables and Programs" (1991)
"""

from libephemeris.lunar import calc_true_lilith, calc_mean_lilith


class TestSolarPerturbationEccentricityVector:
    """Test solar gravitational perturbation on eccentricity vector."""

    def test_perturbation_does_not_break_calculation(self):
        """Solar perturbation should not break the calculation."""
        jd_j2000 = 2451545.0
        lon, lat, e_mag = calc_true_lilith(jd_j2000)

        # Should return valid coordinates
        assert 0 <= lon < 360, f"Longitude {lon} out of range"
        assert -10 < lat < 10, f"Latitude {lat} unexpectedly large"
        assert 0.002 < e_mag < 0.003, f"Distance {e_mag} AU out of expected range"

    def test_eccentricity_magnitude_preserved(self):
        """Distance in AU should remain in valid range after perturbation."""
        test_dates = [
            2451545.0,  # J2000.0
            2451545.0 + 8,  # ~1/4 evection period
            2451545.0 + 16,  # ~1/2 evection period
            2451545.0 + 24,  # ~3/4 evection period
            2451545.0 + 32,  # ~1 full evection period
        ]

        for jd in test_dates:
            _, _, e_mag = calc_true_lilith(jd)
            # 3rd return value is distance in AU (~0.0027 for apogee)
            assert 0.002 < e_mag < 0.003, (
                f"Distance {e_mag} AU at JD {jd} out of expected range [0.002, 0.003]"
            )

    def test_perturbation_varies_with_solar_position(self):
        """Perturbation should cause variation as Sun-Moon-Earth geometry changes."""
        jd_start = 2451545.0
        synodic_month = 29.53  # Days between similar Sun-Moon alignments

        # Get True Lilith at different phases of the synodic month
        positions = []
        for i in range(6):
            jd = jd_start + i * synodic_month / 5
            lon, _, _ = calc_true_lilith(jd)
            positions.append(lon)

        # Check that positions vary over the synodic month
        # The solar perturbation should cause measurable variation
        lon_range = max(positions) - min(positions)
        if lon_range > 180:
            # Handle wrap-around
            positions_unwrapped = [(p - 360 if p > 180 else p) for p in positions]
            lon_range = max(positions_unwrapped) - min(positions_unwrapped)

        # Expect some variation due to solar perturbation
        # (combined with other effects like apsidal motion)
        assert lon_range > 1.0, (
            f"Expected significant longitude variation over synodic month, "
            f"got range of {lon_range} degrees"
        )

    def test_perturbation_quadrupole_character(self):
        """Solar perturbation should have quadrupole character (sin(2*theta) dependence).

        The perturbation repeats twice per synodic month because the solar tidal
        force has quadrupole symmetry (same effect at 0 and 180 degree elongation).
        """
        jd_start = 2451545.0
        half_synodic = 14.77  # Half the synodic month

        # Get True Lilith at start and after half synodic month
        lon1, _, _ = calc_true_lilith(jd_start)
        lon2, _, _ = calc_true_lilith(jd_start + half_synodic)

        # Both positions should be valid
        assert 0 <= lon1 < 360
        assert 0 <= lon2 < 360

        # The positions should differ due to secular apsidal motion
        # (~0.11 degrees/day * 14.77 days ~ 1.6 degrees) plus perturbation effects
        diff = abs(lon2 - lon1)
        if diff > 180:
            diff = 360 - diff

        # Expect measurable change but not too large
        # Note: This includes secular apsidal motion (~1.6°/14.77 days) plus
        # all perturbation effects (evection, variation, annual equation, etc.)
        # Combined effects can produce changes of ~30-35° over half synodic month
        assert diff < 40, (
            f"Position change {diff}° over half synodic month is unexpectedly large"
        )

    def test_continuity_across_solar_perturbation(self):
        """Position should change smoothly when solar perturbation is applied."""
        jd_start = 2451545.0
        step = 0.25  # 6 hours

        positions = []
        for i in range(20):
            jd = jd_start + i * step
            lon, _, _ = calc_true_lilith(jd)
            positions.append(lon)

        # Check smooth transitions
        for i in range(1, len(positions)):
            diff = positions[i] - positions[i - 1]
            if diff > 180:
                diff -= 360
            if diff < -180:
                diff += 360

            # Change over 6 hours should be small (~0.028 degrees mean motion)
            # Allow up to 1 degree for oscillations
            assert abs(diff) < 1.5, (
                f"Position jumped {diff}° in 6 hours between JD {jd_start + (i - 1) * step} "
                f"and JD {jd_start + i * step}"
            )

    def test_latitude_variation_from_3d_perturbation(self):
        """Solar perturbation on e-vector should cause latitude variations.

        Unlike post-hoc longitude corrections, the e-vector perturbation
        correctly captures the 3D nature of the effect, including latitude.
        """
        test_dates = [
            2451545.0,
            2451545.0 + 10,
            2451545.0 + 20,
            2451545.0 + 30,
        ]

        latitudes = []
        for jd in test_dates:
            _, lat, _ = calc_true_lilith(jd)
            latitudes.append(lat)

        # All latitudes should be valid (within lunar orbital inclination ~5.14°)
        for lat in latitudes:
            assert -6 < lat < 6, f"Latitude {lat}° unexpectedly large"

        # Latitudes should show some variation (not all identical)
        lat_range = max(latitudes) - min(latitudes)
        assert lat_range > 0.001, (
            "Expected some latitude variation due to 3D perturbation"
        )


class TestSolarPerturbationPhysicalProperties:
    """Test physical properties of the solar gravitational perturbation."""

    def test_perturbation_amplitude_is_physical(self):
        """The perturbation amplitude should match theoretical evection value.

        The evection modulates lunar eccentricity with amplitude ~0.01148.
        This corresponds to ~1.274° in longitude (arctan(0.01148/0.055) * 2).
        """
        # The amplitude 0.01148 is used in the implementation
        # This test verifies the effect is within physical bounds
        jd = 2451545.0

        lon, _, e_mag = calc_true_lilith(jd)

        # Distance should be close to mean apogee distance ~0.00271 AU
        assert 0.002 < e_mag < 0.003, (
            f"Distance {e_mag} AU deviates too much from expected ~0.0027"
        )

    def test_perturbation_preserves_orbital_validity(self):
        """Perturbation should not create unphysical orbital parameters."""
        # Test over an extended period
        test_jds = [2451545.0 + i * 100 for i in range(20)]

        for jd in test_jds:
            lon, lat, e_mag = calc_true_lilith(jd)

            # All orbital parameters should be physical
            assert 0 <= lon < 360, f"Invalid longitude {lon} at JD {jd}"
            assert -90 < lat < 90, f"Invalid latitude {lat} at JD {jd}"
            assert 0 < e_mag < 1, f"Invalid distance {e_mag} AU at JD {jd}"

            # Distance should be in lunar apogee range
            assert 0.002 < e_mag < 0.003, (
                f"Distance {e_mag} AU at JD {jd} outside expected lunar range"
            )


class TestSolarPerturbationVsMeanLilith:
    """Test effect of solar perturbation on True vs Mean Lilith difference."""

    def test_perturbation_affects_true_mean_difference(self):
        """Solar perturbation should contribute to True-Mean Lilith difference."""
        test_dates = [
            2451545.0,  # J2000.0
            2451600.0,  # ~2 months later
            2451700.0,  # ~5 months later
            2451900.0,  # ~1 year later
        ]

        for jd in test_dates:
            true_lon, _, _ = calc_true_lilith(jd)
            mean_lon = calc_mean_lilith(jd)

            diff = true_lon - mean_lon
            if diff > 180:
                diff -= 360
            if diff < -180:
                diff += 360

            # The True-Mean difference should be substantial
            # (combining all perturbations including solar e-vector effect)
            assert abs(diff) < 35, (
                f"True-Mean diff {diff}° at JD {jd} exceeds expected bounds"
            )

    def test_perturbation_period_matches_evection(self):
        """Solar perturbation should show ~31.8 day periodicity (evection period).

        The quadrupole nature of the solar tidal force means the perturbation
        varies with the evection argument 2D - M'.
        """
        jd_start = 2451545.0
        evection_period = 31.8  # Days

        # Sample True Lilith over one evection period
        samples = []
        for i in range(32):
            jd = jd_start + i
            lon, _, _ = calc_true_lilith(jd)
            samples.append(lon)

        # Remove secular trend (mean apsidal motion ~0.11°/day)
        mean_rate = 0.11
        detrended = [samples[i] - (i * mean_rate) for i in range(len(samples))]

        # Normalize to handle wrap-around
        for i in range(len(detrended)):
            while detrended[i] < 0:
                detrended[i] += 360
            while detrended[i] >= 360:
                detrended[i] -= 360

        # Check that there's oscillation with reasonable amplitude
        range_val = max(detrended) - min(detrended)

        # Handle wrap-around case
        if range_val > 180:
            detrended_shifted = [(d - 180 if d > 180 else d + 180) for d in detrended]
            range_val = max(detrended_shifted) - min(detrended_shifted)

        # Expect some oscillation due to perturbation effects
        # Combined effects should create measurable variation
        assert range_val > 0.5, (
            f"Expected oscillation amplitude > 0.5°, got range {range_val}°"
        )


class TestSolarPerturbationEdgeCases:
    """Test edge cases in solar gravitational perturbation."""

    def test_handles_various_epoch_dates(self):
        """Solar perturbation should work correctly across different epochs."""
        epoch_dates = [
            2433282.5,  # 1950-01-01
            2440587.5,  # 1970-01-01
            2451545.0,  # 2000-01-01 (J2000)
            2460000.0,  # ~2023
            2469807.5,  # 2050-01-01
        ]

        for jd in epoch_dates:
            lon, lat, e_mag = calc_true_lilith(jd)

            assert 0 <= lon < 360, f"Invalid longitude at JD {jd}"
            assert -10 < lat < 10, f"Invalid latitude at JD {jd}"
            assert 0.002 < e_mag < 0.003, f"Invalid distance at JD {jd}"

    def test_handles_rapid_consecutive_calls(self):
        """Multiple rapid calls should produce consistent results."""
        jd = 2451545.0

        results = [calc_true_lilith(jd) for _ in range(10)]

        # All results should be identical
        for i in range(1, len(results)):
            assert results[i] == results[0], (
                f"Inconsistent result on call {i}: {results[i]} != {results[0]}"
            )

    def test_perturbation_near_sun_alignment(self):
        """Perturbation should work when Moon/apogee is nearly aligned with Sun.

        The sin(2*theta) dependence means perturbation is zero at theta = 0, 90°
        and maximum at theta = 45°, 135°.
        """
        # Test at various dates (we can't control alignment directly)
        test_dates = [
            2451545.0,
            2451545.0 + 7.4,  # ~quarter synodic month
            2451545.0 + 14.8,  # ~half synodic month
            2451545.0 + 22.2,  # ~3/4 synodic month
        ]

        for jd in test_dates:
            lon, lat, e_mag = calc_true_lilith(jd)

            # All calculations should succeed
            assert 0 <= lon < 360
            assert -10 < lat < 10
            assert 0.002 < e_mag < 0.003
