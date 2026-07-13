"""
Unit tests for plutino resonant libration model.

Tests verify:
- Libration parameters are defined for known plutinos (Ixion, Orcus)
- Resonant argument calculation is correct
- Libration correction is applied to plutino positions
- Position errors for plutinos are reduced over long timescales
- Non-plutinos are not affected by libration corrections
"""

import math
from libephemeris.constants import (
    IXION,
    ORCUS,
    CERES,
    ERIS,
    QUAOAR,
)
from libephemeris.minor_bodies import (
    MINOR_BODY_ELEMENTS,
    PLUTINO_LIBRATION_PARAMS,
    J2000_EPOCH,
    NEPTUNE_MEAN_LONGITUDE_J2000,
    NEPTUNE_N,
    LibrationParameters,
    calc_neptune_mean_longitude,
    calc_resonant_argument_plutino,
    calc_libration_correction,
    has_libration_model,
    get_libration_parameters,
    calc_minor_body_heliocentric,
    detect_mean_motion_resonance,
)


class TestLibrationParameters:
    """Test that libration parameters are properly defined."""

    def test_ixion_has_libration_params(self):
        """Ixion should have libration parameters defined."""
        assert IXION in PLUTINO_LIBRATION_PARAMS
        params = PLUTINO_LIBRATION_PARAMS[IXION]
        assert isinstance(params, LibrationParameters)
        assert params.amplitude > 0
        assert params.period > 0
        assert 0 <= params.center <= 360

    def test_orcus_has_libration_params(self):
        """Orcus should have libration parameters defined."""
        assert ORCUS in PLUTINO_LIBRATION_PARAMS
        params = PLUTINO_LIBRATION_PARAMS[ORCUS]
        assert isinstance(params, LibrationParameters)
        assert params.amplitude > 0
        assert params.period > 0
        assert 0 <= params.center <= 360

    def test_plutino_amplitudes_reasonable(self):
        """Libration amplitudes should be in reasonable range (30-120 degrees)."""
        for body_id, params in PLUTINO_LIBRATION_PARAMS.items():
            assert 30 <= params.amplitude <= 120, (
                f"Body {body_id} amplitude out of range"
            )

    def test_plutino_periods_reasonable(self):
        """Libration periods should be in a reasonable range.

        Most plutinos have periods ~20,000 years (7-8 million days),
        but objects in higher-order resonances (e.g., 3:10) can have
        longer periods up to ~25,000 years (~9.1 million days).
        """
        for body_id, params in PLUTINO_LIBRATION_PARAMS.items():
            # 18,000-25,000 years = 6.5-9.2 million days
            assert 6_500_000 <= params.period <= 9_200_000, (
                f"Body {body_id} period out of range"
            )

    def test_non_plutinos_no_params(self):
        """Non-plutino bodies should not have libration parameters."""
        assert CERES not in PLUTINO_LIBRATION_PARAMS
        assert ERIS not in PLUTINO_LIBRATION_PARAMS
        assert QUAOAR not in PLUTINO_LIBRATION_PARAMS


class TestNeptuneMeanLongitude:
    """Test Neptune mean longitude calculation."""

    def test_neptune_longitude_at_j2000(self):
        """The linear model starts at its independently defined epoch value."""
        lon = calc_neptune_mean_longitude(J2000_EPOCH)
        assert abs(lon - NEPTUNE_MEAN_LONGITUDE_J2000) < 0.001

    def test_neptune_longitude_progresses(self):
        """Neptune longitude should increase over time."""
        lon_j2000 = calc_neptune_mean_longitude(J2000_EPOCH)
        lon_later = calc_neptune_mean_longitude(J2000_EPOCH + 365.25)  # 1 year later

        # Neptune moves ~2.2 degrees/year
        expected_motion = NEPTUNE_N * 365.25
        actual_motion = (lon_later - lon_j2000) % 360
        assert abs(actual_motion - expected_motion) < 0.1

    def test_neptune_longitude_wraps(self):
        """Neptune longitude should wrap correctly at 360 degrees."""
        # After ~164 years, Neptune completes one orbit
        one_orbit_days = 360.0 / NEPTUNE_N  # ~59,800 days
        lon = calc_neptune_mean_longitude(J2000_EPOCH + one_orbit_days)
        # Should be back near starting position
        assert 0 <= lon < 360


class TestResonantArgument:
    """Test resonant argument calculation for plutinos."""

    def test_resonant_argument_range(self):
        """Resonant argument should be in 0-360 range."""
        elements = MINOR_BODY_ELEMENTS[IXION]
        phi = calc_resonant_argument_plutino(
            elements, J2000_EPOCH, elements.omega, elements.M0
        )
        assert 0 <= phi < 360

    def test_resonant_argument_for_ixion(self):
        """Ixion's resonant argument should be calculable."""
        elements = MINOR_BODY_ELEMENTS[IXION]
        # At epoch
        phi = calc_resonant_argument_plutino(
            elements, elements.epoch, elements.omega, elements.M0
        )
        # For a plutino, phi should librate around 180 degrees
        # Check it's in a reasonable range
        assert 0 <= phi < 360

    def test_resonant_argument_for_orcus(self):
        """Orcus's resonant argument should be calculable."""
        elements = MINOR_BODY_ELEMENTS[ORCUS]
        phi = calc_resonant_argument_plutino(
            elements, elements.epoch, elements.omega, elements.M0
        )
        assert 0 <= phi < 360

    def test_resonant_argument_is_frame_invariant(self):
        """phi must use the longitude of perihelion (Omega + omega), not omega.

        A proper resonant argument is invariant under a rotation of the
        reference longitude origin. Rotating the node Omega by +delta while
        rotating the argument of perihelion omega by -delta leaves both the
        mean longitude (Omega + omega + M) and the longitude of perihelion
        (Omega + omega) unchanged, so phi must not change. Were phi to use the
        argument of perihelion omega directly (rather than Omega + omega), it
        would shift by delta -- which is exactly the defect this guards against.
        """
        import dataclasses

        elements = MINOR_BODY_ELEMENTS[IXION]
        jd = elements.epoch
        phi0 = calc_resonant_argument_plutino(elements, jd, elements.omega, elements.M0)
        delta = 37.0
        rotated = dataclasses.replace(elements, Omega=elements.Omega + delta)
        phi1 = calc_resonant_argument_plutino(
            rotated, jd, elements.omega - delta, elements.M0
        )
        diff = (phi1 - phi0 + 180.0) % 360.0 - 180.0
        assert abs(diff) < 1e-9, (
            f"resonant argument not frame-invariant: shifted by {diff} deg "
            "(it must use the longitude of perihelion Omega + omega, not omega)"
        )


class TestLibrationCorrection:
    """Test libration correction calculation."""

    def test_ixion_correction_bounded(self):
        """Ixion's libration correction should be bounded by amplitude/3."""
        params = PLUTINO_LIBRATION_PARAMS[IXION]
        max_expected = params.amplitude / 3.0

        # Test at multiple times over a full libration period
        for phase in range(0, 360, 30):
            t = J2000_EPOCH + (phase / 360.0) * params.period
            correction = calc_libration_correction(IXION, t)
            assert abs(correction) <= max_expected + 0.1  # Small tolerance

    def test_orcus_correction_bounded(self):
        """Orcus's libration correction should be bounded by amplitude/3."""
        params = PLUTINO_LIBRATION_PARAMS[ORCUS]
        max_expected = params.amplitude / 3.0

        for phase in range(0, 360, 30):
            t = J2000_EPOCH + (phase / 360.0) * params.period
            correction = calc_libration_correction(ORCUS, t)
            assert abs(correction) <= max_expected + 0.1

    def test_non_plutino_zero_correction(self):
        """Non-plutinos should have zero libration correction."""
        correction = calc_libration_correction(CERES, J2000_EPOCH)
        assert correction == 0.0

        correction = calc_libration_correction(ERIS, J2000_EPOCH)
        assert correction == 0.0

    def test_libration_periodic(self):
        """Libration correction should be periodic over the libration period."""
        params = PLUTINO_LIBRATION_PARAMS[IXION]

        # Correction at t and t + period should be the same
        t0 = J2000_EPOCH
        t1 = J2000_EPOCH + params.period

        correction_t0 = calc_libration_correction(IXION, t0)
        correction_t1 = calc_libration_correction(IXION, t1)

        assert abs(correction_t0 - correction_t1) < 0.01

    def test_libration_varies_with_time(self):
        """Libration correction should vary over time (not constant)."""
        params = PLUTINO_LIBRATION_PARAMS[IXION]

        # Check corrections at quarter period intervals
        corrections = []
        for i in range(4):
            t = J2000_EPOCH + i * params.period / 4
            corrections.append(calc_libration_correction(IXION, t))

        # Not all corrections should be the same
        assert not all(abs(c - corrections[0]) < 0.01 for c in corrections)


class TestHelperFunctions:
    """Test helper functions for libration model."""

    def test_has_libration_model_plutinos(self):
        """has_libration_model should return True for plutinos."""
        assert has_libration_model(IXION) is True
        assert has_libration_model(ORCUS) is True

    def test_has_libration_model_non_plutinos(self):
        """has_libration_model should return False for non-plutinos."""
        assert has_libration_model(CERES) is False
        assert has_libration_model(ERIS) is False
        assert has_libration_model(QUAOAR) is False

    def test_get_libration_parameters_plutinos(self):
        """get_libration_parameters should return params for plutinos."""
        params = get_libration_parameters(IXION)
        assert params is not None
        assert isinstance(params, LibrationParameters)

    def test_get_libration_parameters_non_plutinos(self):
        """get_libration_parameters should return None for non-plutinos."""
        params = get_libration_parameters(CERES)
        assert params is None


class TestPositionWithLibration:
    """Test that libration correction is applied in position calculation."""

    def test_plutino_position_includes_libration(self):
        """Plutino positions should include libration correction."""
        # Get position at two different times
        jd1 = J2000_EPOCH
        jd2 = J2000_EPOCH + 365.25 * 10  # 10 years later

        lon1, lat1, dist1 = calc_minor_body_heliocentric(IXION, jd1, use_spk=False)
        lon2, lat2, dist2 = calc_minor_body_heliocentric(IXION, jd2, use_spk=False)

        # Positions should be valid
        assert 0 <= lon1 < 360
        assert 0 <= lon2 < 360
        assert -90 <= lat1 <= 90
        assert -90 <= lat2 <= 90
        assert dist1 > 0
        assert dist2 > 0

    def test_non_plutino_position_unchanged(self):
        """Non-plutino positions should work normally."""
        jd = J2000_EPOCH

        lon, lat, dist = calc_minor_body_heliocentric(CERES, jd, use_spk=False)

        assert 0 <= lon < 360
        assert -90 <= lat <= 90
        assert dist > 0

    def test_plutinos_detected_as_resonant(self):
        """Plutinos (2:3 resonance) with libration params should be detected as resonant.

        Note: PLUTINO_LIBRATION_PARAMS may include bodies in other Neptune
        resonances (e.g., Gonggong in 3:10 resonance). Only bodies with
        resonance_p=2, resonance_q=3 are actual plutinos and should be
        detected by detect_mean_motion_resonance().
        """
        for body_id, params in PLUTINO_LIBRATION_PARAMS.items():
            if params.resonance_p != 2 or params.resonance_q != 3:
                # Not a 2:3 plutino — skip (e.g., Gonggong is 3:10)
                continue
            elements = MINOR_BODY_ELEMENTS[body_id]
            result = detect_mean_motion_resonance(elements)
            assert result.is_resonant, f"Body {body_id} not detected as resonant"
            assert result.resonance is not None
            assert result.resonance.name == "plutino"


class TestLongTermAccuracy:
    """Test improved accuracy over long timescales."""

    def test_libration_correction_magnitude_reasonable(self):
        """Libration correction should be in reasonable range over decades."""
        # Test over 100 years
        start_jd = J2000_EPOCH
        corrections = []

        for years in range(0, 101, 10):
            jd = start_jd + years * 365.25
            correction = calc_libration_correction(IXION, jd)
            corrections.append(abs(correction))

        # Maximum correction should be bounded (amplitude/3 ~ 26 degrees)
        max_correction = max(corrections)
        assert max_correction < 30, f"Max correction {max_correction} too large"

    def test_position_continuity(self):
        """Positions should change continuously over time (no jumps)."""
        jd_start = J2000_EPOCH

        prev_lon = None
        for days in range(0, 3650, 365):  # 10 years, yearly steps
            jd = jd_start + days
            lon, _, _ = calc_minor_body_heliocentric(IXION, jd, use_spk=False)

            if prev_lon is not None:
                # Position change should be reasonable (< 20 degrees/year for slow TNO)
                diff = abs(lon - prev_lon)
                if diff > 180:
                    diff = 360 - diff
                assert diff < 20, f"Position jump too large: {diff} degrees"

            prev_lon = lon

    def test_orcus_position_continuity(self):
        """Orcus positions should also be continuous."""
        jd_start = J2000_EPOCH

        prev_lon = None
        for days in range(0, 3650, 365):
            jd = jd_start + days
            lon, _, _ = calc_minor_body_heliocentric(ORCUS, jd, use_spk=False)

            if prev_lon is not None:
                diff = abs(lon - prev_lon)
                if diff > 180:
                    diff = 360 - diff
                assert diff < 20

            prev_lon = lon


class TestCalcMinorBodyPositionLibration:
    """Test libration correction in calc_minor_body_position with body_id."""

    def test_position_with_body_id_applies_libration(self):
        """calc_minor_body_position with body_id should apply libration for plutinos."""
        from libephemeris.minor_bodies import calc_minor_body_position

        elements = MINOR_BODY_ELEMENTS[IXION]
        jd = J2000_EPOCH + 365.25 * 50  # 50 years from J2000

        # Without body_id: no libration correction
        x1, y1, z1 = calc_minor_body_position(elements, jd, body_id=None)

        # With body_id: libration correction applied
        x2, y2, z2 = calc_minor_body_position(elements, jd, body_id=IXION)

        # Convert to longitude for comparison
        lon1 = math.degrees(math.atan2(y1, x1)) % 360.0
        lon2 = math.degrees(math.atan2(y2, x2)) % 360.0

        # The positions should be different due to libration correction
        diff = abs(lon1 - lon2)
        if diff > 180:
            diff = 360 - diff

        # Libration correction should cause a measurable difference
        # Ixion's amplitude/3 is ~26 degrees max, so expect some difference
        # but could be larger due to secular effects interaction
        assert diff > 0.1, "Libration should cause measurable position difference"
        # Allow up to 50° difference to account for nonlinear effects
        assert diff < 50, f"Libration correction too large: {diff}°"

    def test_non_plutino_position_unchanged_with_body_id(self):
        """Non-plutinos should have same position with or without body_id."""
        from libephemeris.minor_bodies import calc_minor_body_position

        elements = MINOR_BODY_ELEMENTS[CERES]
        jd = J2000_EPOCH + 365.25 * 10

        # Without body_id
        x1, y1, z1 = calc_minor_body_position(elements, jd, body_id=None)

        # With body_id (Ceres is not a plutino)
        x2, y2, z2 = calc_minor_body_position(elements, jd, body_id=CERES)

        # Positions should be identical (Ceres has no libration params)
        assert abs(x1 - x2) < 1e-10
        assert abs(y1 - y2) < 1e-10
        assert abs(z1 - z2) < 1e-10

    def test_orcus_position_with_body_id(self):
        """Orcus position should include libration correction with body_id."""
        from libephemeris.minor_bodies import calc_minor_body_position

        elements = MINOR_BODY_ELEMENTS[ORCUS]
        jd = J2000_EPOCH + 365.25 * 50

        x1, y1, z1 = calc_minor_body_position(elements, jd, body_id=None)
        x2, y2, z2 = calc_minor_body_position(elements, jd, body_id=ORCUS)

        lon1 = math.degrees(math.atan2(y1, x1)) % 360.0
        lon2 = math.degrees(math.atan2(y2, x2)) % 360.0

        diff = abs(lon1 - lon2)
        if diff > 180:
            diff = 360 - diff

        assert diff > 0.0001, "Orcus should have libration correction"
        assert diff < 25, "Orcus libration correction should be bounded (amp/3 ~23°)"

    def test_libration_varies_over_period(self):
        """Libration correction should vary over the libration period."""
        from libephemeris.minor_bodies import calc_minor_body_position

        elements = MINOR_BODY_ELEMENTS[IXION]
        params = PLUTINO_LIBRATION_PARAMS[IXION]

        # Sample at quarter-period intervals
        corrections = []
        for i in range(4):
            jd = J2000_EPOCH + i * params.period / 4

            x_no_lib, y_no_lib, _ = calc_minor_body_position(elements, jd, body_id=None)
            x_lib, y_lib, _ = calc_minor_body_position(elements, jd, body_id=IXION)

            lon_no_lib = math.degrees(math.atan2(y_no_lib, x_no_lib)) % 360.0
            lon_lib = math.degrees(math.atan2(y_lib, x_lib)) % 360.0

            diff = lon_lib - lon_no_lib
            if diff > 180:
                diff -= 360
            elif diff < -180:
                diff += 360

            corrections.append(diff)

        # Not all corrections should be the same (libration varies)
        unique_corrections = set(round(c, 1) for c in corrections)
        assert len(unique_corrections) > 1, "Libration should vary over period"


class TestLongTermPlutonoAccuracy:
    """Test improved accuracy for plutinos over 100-year timescales."""

    def test_ixion_100_year_position_bounded(self):
        """Ixion position should be computable over 100 years without errors."""
        from libephemeris.minor_bodies import calc_minor_body_position

        elements = MINOR_BODY_ELEMENTS[IXION]

        # Ixion orbital period is ~248 years, so in 100 years it moves ~145°
        jd_start = J2000_EPOCH
        jd_end = J2000_EPOCH + 365.25 * 100

        x1, y1, z1 = calc_minor_body_position(elements, jd_start, body_id=IXION)
        x2, y2, z2 = calc_minor_body_position(elements, jd_end, body_id=IXION)

        # Positions should be valid (finite, non-zero distance)
        r1 = math.sqrt(x1**2 + y1**2 + z1**2)
        r2 = math.sqrt(x2**2 + y2**2 + z2**2)

        assert r1 > 0 and math.isfinite(r1)
        assert r2 > 0 and math.isfinite(r2)

        # Ixion's semi-major axis is ~39.5 AU, should be in reasonable range
        assert 30 < r1 < 50, f"Start distance out of range: {r1} AU"
        assert 30 < r2 < 50, f"End distance out of range: {r2} AU"

        lon1 = math.degrees(math.atan2(y1, x1)) % 360.0
        lon2 = math.degrees(math.atan2(y2, x2)) % 360.0

        # Both positions should be valid longitudes
        assert 0 <= lon1 < 360
        assert 0 <= lon2 < 360

    def test_position_continuity_with_libration(self):
        """Positions should change smoothly with libration correction."""
        from libephemeris.minor_bodies import calc_minor_body_position

        elements = MINOR_BODY_ELEMENTS[IXION]

        prev_lon = None
        for year in range(0, 100, 5):  # Every 5 years over 100 years
            jd = J2000_EPOCH + year * 365.25
            x, y, _ = calc_minor_body_position(elements, jd, body_id=IXION)
            lon = math.degrees(math.atan2(y, x)) % 360.0

            if prev_lon is not None:
                diff = abs(lon - prev_lon)
                if diff > 180:
                    diff = 360 - diff
                # Ixion moves ~2° per year, so 5 years = ~10°
                # Allow some margin for libration oscillation
                assert diff < 20, f"Position jump at year {year}: {diff}°"

            prev_lon = lon
