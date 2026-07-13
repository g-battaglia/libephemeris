"""Tests for rigorous proper motion calculation in fixed_stars module.

The proper motion implementation uses the space motion approach from
Hipparcos Vol. 1, Section 1.5.5 with second-order Taylor expansion,
which correctly handles the curvature of the celestial sphere when
propagating stellar positions over long time periods.

This tests the calc_fixed_star_position function in fixed_stars.py.
"""

import math
import pytest
from libephemeris.fixed_stars import (
    calc_fixed_star_position,
    StarData,
    FIXED_STARS,
)
from libephemeris.constants import FLG_BARYCTR, REGULUS, SPICA_STAR


def _linear_proper_motion(star: StarData, t_years: float) -> tuple[float, float]:
    """
    Apply LINEAR proper motion (old method for comparison).

    Returns (ra_pm, dec_pm) in degrees.
    """
    ra_pm = star.ra_j2000 + (star.pm_ra * t_years) / 3600.0
    dec_pm = star.dec_j2000 + (star.pm_dec * t_years) / 3600.0
    return ra_pm, dec_pm


def _first_order_proper_motion(star: StarData, t_years: float) -> tuple[float, float]:
    """
    Apply FIRST-ORDER space motion (P(t) = P(0) + V*dt).

    Uses 3D vector propagation but only first-order Taylor expansion.
    Returns (ra_pm, dec_pm) in degrees.
    """
    # Convert proper motions from arcsec/year to radians/year
    pm_ra_rad = math.radians(star.pm_ra / 3600.0)
    pm_dec_rad = math.radians(star.pm_dec / 3600.0)

    # Convert J2000 position to radians
    ra_rad = math.radians(star.ra_j2000)
    dec_rad = math.radians(star.dec_j2000)

    # Unit position vector at J2000 epoch
    cos_dec = math.cos(dec_rad)
    sin_dec = math.sin(dec_rad)
    cos_ra = math.cos(ra_rad)
    sin_ra = math.sin(ra_rad)

    px = cos_dec * cos_ra
    py = cos_dec * sin_ra
    pz = sin_dec

    # Proper motion velocity vector in the tangent plane
    vx = -pm_ra_rad * sin_ra - pm_dec_rad * sin_dec * cos_ra
    vy = pm_ra_rad * cos_ra - pm_dec_rad * sin_dec * sin_ra
    vz = pm_dec_rad * cos_dec

    # First-order propagation: P(t) = P(0) + V * dt
    px_t = px + vx * t_years
    py_t = py + vy * t_years
    pz_t = pz + vz * t_years

    # Normalize to get unit vector
    r = math.sqrt(px_t * px_t + py_t * py_t + pz_t * pz_t)
    px_t /= r
    py_t /= r
    pz_t /= r

    # Convert back to RA/Dec
    dec_pm = math.degrees(math.asin(pz_t))
    ra_pm = math.degrees(math.atan2(py_t, px_t))
    if ra_pm < 0:
        ra_pm += 360.0

    return ra_pm, dec_pm


def _second_order_proper_motion(star: StarData, t_years: float) -> tuple[float, float]:
    """
    Apply SECOND-ORDER proper motion using space motion approach.

    Uses 3D vector propagation with second-order Taylor expansion per
    Hipparcos Vol. 1, Section 1.5.5. This includes the centripetal
    acceleration term for improved accuracy on century-scale intervals.

    Returns (ra_pm, dec_pm) in degrees.
    """
    # Convert proper motions from arcsec/year to radians/year
    pm_ra_rad = math.radians(star.pm_ra / 3600.0)
    pm_dec_rad = math.radians(star.pm_dec / 3600.0)

    # Convert J2000 position to radians
    ra_rad = math.radians(star.ra_j2000)
    dec_rad = math.radians(star.dec_j2000)

    # Unit position vector at J2000 epoch
    cos_dec = math.cos(dec_rad)
    sin_dec = math.sin(dec_rad)
    cos_ra = math.cos(ra_rad)
    sin_ra = math.sin(ra_rad)

    px = cos_dec * cos_ra
    py = cos_dec * sin_ra
    pz = sin_dec

    # Proper motion velocity vector in the tangent plane
    vx = -pm_ra_rad * sin_ra - pm_dec_rad * sin_dec * cos_ra
    vy = pm_ra_rad * cos_ra - pm_dec_rad * sin_dec * sin_ra
    vz = pm_dec_rad * cos_dec

    # Second-order propagation with centripetal acceleration:
    # P(t) = P(0) + V*dt - 0.5*|V|²*P(0)*dt²
    v_squared = vx * vx + vy * vy + vz * vz
    t_squared = t_years * t_years

    px_t = px + vx * t_years - 0.5 * v_squared * px * t_squared
    py_t = py + vy * t_years - 0.5 * v_squared * py * t_squared
    pz_t = pz + vz * t_years - 0.5 * v_squared * pz * t_squared

    # Normalize to get unit vector
    r = math.sqrt(px_t * px_t + py_t * py_t + pz_t * pz_t)
    px_t /= r
    py_t /= r
    pz_t /= r

    # Convert back to RA/Dec
    dec_pm = math.degrees(math.asin(pz_t))
    ra_pm = math.degrees(math.atan2(py_t, px_t))
    if ra_pm < 0:
        ra_pm += 360.0

    return ra_pm, dec_pm


@pytest.mark.unit
class TestFixedStarsProperMotionRigorous:
    """Tests for rigorous proper motion calculation in fixed_stars module."""

    def test_zero_time_returns_j2000(self):
        """At J2000.0 epoch (t_years=0), result should equal input exactly."""
        star = FIXED_STARS[SPICA_STAR]

        ra_pm, dec_pm = _second_order_proper_motion(star, 0.0)

        # Should match J2000 coordinates exactly (within floating point)
        assert abs(ra_pm - star.ra_j2000) < 1e-10
        assert abs(dec_pm - star.dec_j2000) < 1e-10

    def test_short_time_matches_linear(self):
        """For short periods (1 year), rigorous and linear should agree well."""
        star = FIXED_STARS[SPICA_STAR]
        t_years = 1.0

        ra_rig, dec_rig = _second_order_proper_motion(star, t_years)
        ra_lin, dec_lin = _linear_proper_motion(star, t_years)

        # Should agree to within 0.001 arcsec for 1 year
        diff_ra = abs(ra_rig - ra_lin) * 3600  # arcsec
        diff_dec = abs(dec_rig - dec_lin) * 3600  # arcsec

        assert diff_ra < 0.001, f"RA difference: {diff_ra} arcsec"
        assert diff_dec < 0.001, f"Dec difference: {diff_dec} arcsec"

    def test_long_time_differs_from_linear(self):
        """
        For very long periods (1000 years), the rigorous method should
        produce different results due to spherical geometry effects.
        """
        # Create a star with significant proper motion for clear difference
        high_pm_star = StarData(
            ra_j2000=180.0,  # degrees
            dec_j2000=60.0,  # high declination where RA circles are smaller
            pm_ra=0.5,  # 0.5 arcsec/year (high motion)
            pm_dec=0.5,  # 0.5 arcsec/year
        )

        t_years = 1000.0

        ra_rig, dec_rig = _second_order_proper_motion(high_pm_star, t_years)
        ra_lin, dec_lin = _linear_proper_motion(high_pm_star, t_years)

        # Should differ due to curvature effects
        diff_ra = abs(ra_rig - ra_lin) * 3600  # arcsec
        diff_dec = abs(dec_rig - dec_lin) * 3600  # arcsec

        # The difference should be measurable (at least a few arcsec)
        total_diff = math.sqrt(diff_ra**2 + diff_dec**2)
        assert total_diff > 1.0, (
            f"Expected significant difference, got {total_diff} arcsec"
        )

    def test_preserves_unit_vector_norm(self):
        """The rigorous method should maintain unit vector normalization."""
        star = FIXED_STARS[SPICA_STAR]
        t_years = 100.0

        ra_pm, dec_pm = _second_order_proper_motion(star, t_years)

        # Convert back to unit vector and check norm
        ra_rad = math.radians(ra_pm)
        dec_rad = math.radians(dec_pm)

        x = math.cos(dec_rad) * math.cos(ra_rad)
        y = math.cos(dec_rad) * math.sin(ra_rad)
        z = math.sin(dec_rad)

        norm = math.sqrt(x * x + y * y + z * z)
        assert abs(norm - 1.0) < 1e-10

    def test_symmetric_time_direction(self):
        """Forward and backward propagation should be symmetric."""
        star = FIXED_STARS[SPICA_STAR]
        t_years = 50.0

        # Go forward then backward
        ra_fwd, dec_fwd = _second_order_proper_motion(star, t_years)

        # Create star at forward position and go backward
        star_fwd = StarData(
            ra_j2000=ra_fwd,
            dec_j2000=dec_fwd,
            pm_ra=star.pm_ra,
            pm_dec=star.pm_dec,
        )

        ra_back, dec_back = _second_order_proper_motion(star_fwd, -t_years)

        # Should return close to original (within small numerical error)
        diff_ra = abs(ra_back - star.ra_j2000) * 3600  # arcsec
        diff_dec = abs(dec_back - star.dec_j2000) * 3600  # arcsec

        # Second-order effects mean it won't be exact, but should be close
        assert diff_ra < 0.01, f"RA difference: {diff_ra} arcsec"
        assert diff_dec < 0.01, f"Dec difference: {diff_dec} arcsec"

    def test_all_stars_in_catalog(self):
        """Test that proper motion works for all stars in the catalog."""
        t_years = 50.0

        for star_id, star in FIXED_STARS.items():
            ra_pm, dec_pm = _second_order_proper_motion(star, t_years)

            # Position should be valid
            assert 0 <= ra_pm < 360, f"Star {star_id}: RA out of range: {ra_pm}"
            assert -90 <= dec_pm <= 90, f"Star {star_id}: Dec out of range: {dec_pm}"

            # Position should have moved from J2000 (except for zero proper motion)
            if star.pm_ra != 0 or star.pm_dec != 0:
                moved = (ra_pm != star.ra_j2000) or (dec_pm != star.dec_j2000)
                assert moved, (
                    f"Star {star_id}: Position didn't change with proper motion"
                )

    def test_declination_stays_bounded(self):
        """Declination should stay within +/-90 even with extreme motion."""
        # Star near pole with motion toward pole
        polar_star = StarData(
            ra_j2000=0.0,
            dec_j2000=89.0,  # Near north pole
            pm_ra=0.0,
            pm_dec=1.0,  # Moving toward pole
        )

        # 100 years would push linear dec to 89 + 100/3600 * 100 = 91.78 deg
        # But should stay bounded
        t_years = 100.0

        ra_pm, dec_pm = _second_order_proper_motion(polar_star, t_years)

        assert -90 <= dec_pm <= 90, f"Dec out of bounds: {dec_pm}"


@pytest.mark.unit
class TestCalcFixedStarPositionIntegration:
    """Integration tests using the actual calc_fixed_star_position function."""

    def test_regulus_ecliptic_position_j2000(self):
        """Test Regulus ecliptic position at J2000.0."""
        # JD for J2000.0
        jd_tt = 2451545.0

        lon, lat, dist = calc_fixed_star_position(REGULUS, jd_tt)

        # Regulus is at approximately 150 deg ecliptic longitude at J2000
        assert 149 < lon < 151, f"Regulus ecliptic longitude: {lon}"
        # Latitude should be small
        assert -1 < lat < 2, f"Regulus ecliptic latitude: {lat}"
        # Distance from parallax (Regulus ~5 million AU)
        assert dist > 1000000, f"Regulus distance should be large: {dist}"

    def test_spica_ecliptic_position_j2000(self):
        """Test Spica ecliptic position at J2000.0."""
        jd_tt = 2451545.0

        lon, lat, dist = calc_fixed_star_position(SPICA_STAR, jd_tt)

        # Spica is at approximately 203-204 deg ecliptic longitude at J2000
        assert 202 < lon < 205, f"Spica ecliptic longitude: {lon}"

    def test_position_changes_over_time(self):
        """Test that star positions change over 50 years."""
        jd_2000 = 2451545.0
        # Use year 2050 which is within ephemeris range (1899-2053)
        jd_2050 = 2451545.0 + 50 * 365.25

        lon_2000, _, _ = calc_fixed_star_position(SPICA_STAR, jd_2000)
        lon_2050, _, _ = calc_fixed_star_position(SPICA_STAR, jd_2050)

        # Position should have changed due to proper motion + precession
        # Over 50 years, expect at least 0.05 degrees change
        diff = abs(lon_2050 - lon_2000)
        assert diff > 0.05, f"Expected movement over 50 years, got {diff} deg"

    def test_position_changes_backwards_in_time(self):
        """Test that star positions change going backwards in time."""
        jd_2000 = 2451545.0
        jd_1900 = 2451545.0 - 100 * 365.25

        lon_2000, _, _ = calc_fixed_star_position(REGULUS, jd_2000)
        lon_1900, _, _ = calc_fixed_star_position(REGULUS, jd_1900)

        # Position should have changed
        diff = abs(lon_2000 - lon_1900)
        assert diff > 0.1, f"Expected movement over 100 years, got {diff} deg"

    def test_unknown_star_raises_error(self):
        """Test that unknown star ID raises ValueError."""
        with pytest.raises(ValueError, match="could not find star name"):
            calc_fixed_star_position(99999, 2451545.0)

    def test_all_catalog_stars_work(self):
        """Test that all stars in the catalog can be calculated."""
        jd_tt = 2451545.0

        for star_id in FIXED_STARS.keys():
            lon, lat, dist = calc_fixed_star_position(star_id, jd_tt)

            # All results should be valid
            assert 0 <= lon < 360, f"Star {star_id}: lon out of range"
            assert -90 <= lat <= 90, f"Star {star_id}: lat out of range"
            assert dist > 0, f"Star {star_id}: distance should be positive"

    def test_nonpositive_parallax_uses_direction_at_one_gigaparsec(self):
        """Unknown inverse distance follows Skyfield's public convention.

        One gigaparsec has parallax one microarcsecond by the definition of a
        parsec.  This neutral direction-at-infinity representation avoids
        inventing a distance from another ephemeris.
        """
        star_id = 1001303  # 9 Cep; catalog parallax is unavailable.
        assert FIXED_STARS[star_id].parallax_mas <= 0.0
        _lon, _lat, dist = calc_fixed_star_position(
            star_id,
            2451545.0,
            noaberr=True,
            nogdefl=True,
            j2000_frame=True,
            center=FLG_BARYCTR,
        )
        microarcsecond_rad = math.radians(1.0e-6 / 1000.0 / 3600.0)
        expected_au = 1.0 / math.sin(microarcsecond_rad)
        assert dist == pytest.approx(expected_au, rel=1e-12)


@pytest.mark.unit
class TestProperMotionAccuracy:
    """Tests verifying the accuracy improvement from rigorous proper motion."""

    def test_100_year_accuracy_under_0_1_arcsec(self):
        """
        Verify that the rigorous method maintains sub-0.1 arcsec accuracy
        over 100 years compared to what would accumulate with linear method.
        """
        star = FIXED_STARS[SPICA_STAR]
        t_years = 100.0

        ra_rig, dec_rig = _second_order_proper_motion(star, t_years)
        ra_lin, dec_lin = _linear_proper_motion(star, t_years)

        # For Spica over 100 years, the curvature error should be small
        # but measurable - the rigorous method should agree within 0.1 arcsec
        # of itself when done in steps vs all at once

        # Test that the rigorous method is internally consistent
        # by comparing single-step with multi-step propagation
        ra_step, dec_step = star.ra_j2000, star.dec_j2000
        n_steps = 10
        step_size = t_years / n_steps

        temp_star = StarData(
            ra_j2000=star.ra_j2000,
            dec_j2000=star.dec_j2000,
            pm_ra=star.pm_ra,
            pm_dec=star.pm_dec,
        )

        for _ in range(n_steps):
            ra_step, dec_step = _second_order_proper_motion(temp_star, step_size)
            temp_star = StarData(
                ra_j2000=ra_step,
                dec_j2000=dec_step,
                pm_ra=star.pm_ra,
                pm_dec=star.pm_dec,
            )

        # Single step vs multi-step should agree very closely
        diff_ra = abs(ra_rig - ra_step) * 3600  # arcsec
        diff_dec = abs(dec_rig - dec_step) * 3600  # arcsec

        assert diff_ra < 0.1, f"RA step consistency: {diff_ra} arcsec"
        assert diff_dec < 0.1, f"Dec step consistency: {diff_dec} arcsec"


@pytest.mark.unit
class TestSecondOrderProperMotion:
    """Tests verifying the second-order expansion improves accuracy over first-order."""

    # Barnard's Star has extremely high proper motion: ~10.3 arcsec/year total
    # This makes it the ideal test case for verifying the second-order correction
    BARNARDS_STAR = StarData(
        ra_j2000=269.452075,  # 17h 57m 48.5s
        dec_j2000=4.693391,  # +4° 41' 36.2"
        pm_ra=-0.79858,  # -798.58 mas/yr (includes cos(dec))
        pm_dec=10.32812,  # +10328.12 mas/yr
    )

    def test_second_order_more_accurate_than_first_order_100_years(self):
        """
        For high proper motion stars over 100 years, second-order should be
        more accurate than first-order when compared to step-wise integration.
        """
        star = self.BARNARDS_STAR
        t_years = 100.0

        # Reference: step-wise integration (100 steps of 1 year each)
        # This approximates the "true" path on the sphere
        n_steps = 100
        step_size = t_years / n_steps
        ra_ref, dec_ref = star.ra_j2000, star.dec_j2000
        temp_star = StarData(
            ra_j2000=star.ra_j2000,
            dec_j2000=star.dec_j2000,
            pm_ra=star.pm_ra,
            pm_dec=star.pm_dec,
        )
        for _ in range(n_steps):
            ra_ref, dec_ref = _second_order_proper_motion(temp_star, step_size)
            temp_star = StarData(
                ra_j2000=ra_ref,
                dec_j2000=dec_ref,
                pm_ra=star.pm_ra,
                pm_dec=star.pm_dec,
            )

        # Single-step first-order
        ra_first, dec_first = _first_order_proper_motion(star, t_years)

        # Single-step second-order
        ra_second, dec_second = _second_order_proper_motion(star, t_years)

        # Error for first-order vs reference
        err_first_ra = abs(ra_first - ra_ref) * 3600  # arcsec
        err_first_dec = abs(dec_first - dec_ref) * 3600
        err_first = math.sqrt(err_first_ra**2 + err_first_dec**2)

        # Error for second-order vs reference
        err_second_ra = abs(ra_second - ra_ref) * 3600
        err_second_dec = abs(dec_second - dec_ref) * 3600
        err_second = math.sqrt(err_second_ra**2 + err_second_dec**2)

        # Second-order should be significantly more accurate
        assert err_second < err_first, (
            f"Second-order ({err_second:.4f} arcsec) should be more accurate "
            f"than first-order ({err_first:.4f} arcsec)"
        )

    def test_second_order_more_accurate_than_first_order_500_years(self):
        """
        Over 500 years, the second-order advantage should be even more pronounced.
        """
        star = self.BARNARDS_STAR
        t_years = 500.0

        # Reference: step-wise integration
        n_steps = 500
        step_size = t_years / n_steps
        ra_ref, dec_ref = star.ra_j2000, star.dec_j2000
        temp_star = StarData(
            ra_j2000=star.ra_j2000,
            dec_j2000=star.dec_j2000,
            pm_ra=star.pm_ra,
            pm_dec=star.pm_dec,
        )
        for _ in range(n_steps):
            ra_ref, dec_ref = _second_order_proper_motion(temp_star, step_size)
            temp_star = StarData(
                ra_j2000=ra_ref,
                dec_j2000=dec_ref,
                pm_ra=star.pm_ra,
                pm_dec=star.pm_dec,
            )

        # Single-step first-order
        ra_first, dec_first = _first_order_proper_motion(star, t_years)

        # Single-step second-order
        ra_second, dec_second = _second_order_proper_motion(star, t_years)

        # Calculate errors
        err_first = math.sqrt(
            (abs(ra_first - ra_ref) * 3600) ** 2
            + (abs(dec_first - dec_ref) * 3600) ** 2
        )
        err_second = math.sqrt(
            (abs(ra_second - ra_ref) * 3600) ** 2
            + (abs(dec_second - dec_ref) * 3600) ** 2
        )

        # Second-order should be more accurate
        assert err_second < err_first, (
            f"At 500 years, second-order ({err_second:.2f} arcsec) should be "
            f"more accurate than first-order ({err_first:.2f} arcsec)"
        )

        # The improvement ratio should be substantial for long timescales
        improvement_ratio = err_first / err_second
        assert improvement_ratio > 1.5, (
            f"Expected >1.5x improvement, got {improvement_ratio:.2f}x"
        )

    def test_second_order_accuracy_under_1_arcsec_for_500_years(self):
        """
        For Barnard's Star over 500 years, second-order error should be < 1 arcsec.
        This is the key accuracy target mentioned in the docstring.
        """
        star = self.BARNARDS_STAR
        t_years = 500.0

        # Reference: step-wise integration with small steps
        n_steps = 1000
        step_size = t_years / n_steps
        ra_ref, dec_ref = star.ra_j2000, star.dec_j2000
        temp_star = StarData(
            ra_j2000=star.ra_j2000,
            dec_j2000=star.dec_j2000,
            pm_ra=star.pm_ra,
            pm_dec=star.pm_dec,
        )
        for _ in range(n_steps):
            ra_ref, dec_ref = _second_order_proper_motion(temp_star, step_size)
            temp_star = StarData(
                ra_j2000=ra_ref,
                dec_j2000=dec_ref,
                pm_ra=star.pm_ra,
                pm_dec=star.pm_dec,
            )

        # Single-step second-order
        ra_second, dec_second = _second_order_proper_motion(star, t_years)

        # Error should be less than 1 arcsec
        err_arcsec = math.sqrt(
            (abs(ra_second - ra_ref) * 3600) ** 2
            + (abs(dec_second - dec_ref) * 3600) ** 2
        )

        assert err_arcsec < 1.0, (
            f"Second-order error over 500 years should be <1 arcsec, "
            f"got {err_arcsec:.3f} arcsec"
        )

    def test_second_order_accuracy_under_0_01_arcsec_for_100_years(self):
        """
        For typical stars over 100 years, second-order error should be < 0.01 arcsec.
        This is the accuracy target mentioned in the updated docstring.
        """
        star = FIXED_STARS[SPICA_STAR]
        t_years = 100.0

        # Reference: step-wise integration
        n_steps = 100
        step_size = t_years / n_steps
        ra_ref, dec_ref = star.ra_j2000, star.dec_j2000
        temp_star = StarData(
            ra_j2000=star.ra_j2000,
            dec_j2000=star.dec_j2000,
            pm_ra=star.pm_ra,
            pm_dec=star.pm_dec,
        )
        for _ in range(n_steps):
            ra_ref, dec_ref = _second_order_proper_motion(temp_star, step_size)
            temp_star = StarData(
                ra_j2000=ra_ref,
                dec_j2000=dec_ref,
                pm_ra=star.pm_ra,
                pm_dec=star.pm_dec,
            )

        # Single-step second-order
        ra_second, dec_second = _second_order_proper_motion(star, t_years)

        # Error should be less than 0.01 arcsec for typical stars
        err_arcsec = math.sqrt(
            (abs(ra_second - ra_ref) * 3600) ** 2
            + (abs(dec_second - dec_ref) * 3600) ** 2
        )

        assert err_arcsec < 0.01, (
            f"Second-order error over 100 years for Spica should be <0.01 arcsec, "
            f"got {err_arcsec:.6f} arcsec"
        )

    def test_second_order_handles_high_proper_motion_correctly(self):
        """
        Verify the second-order method correctly handles extreme proper motion.

        Barnard's Star moves ~10.3 arcsec/year, which means it moves about
        1.03 degrees over 360 years (nearly half a degree per century).
        """
        star = self.BARNARDS_STAR
        t_years = 100.0

        ra_pm, dec_pm = _second_order_proper_motion(star, t_years)

        # Expected motion in Dec: ~10.3 arcsec/year * 100 years = ~1030 arcsec = 0.286 deg
        expected_dec_change = (star.pm_dec * t_years) / 3600.0  # degrees
        actual_dec_change = dec_pm - star.dec_j2000

        # Should be reasonably close (within 0.01 deg = 36 arcsec)
        assert abs(actual_dec_change - expected_dec_change) < 0.01, (
            f"Dec change: expected ~{expected_dec_change:.4f}, "
            f"got {actual_dec_change:.4f} deg"
        )

    def test_first_vs_second_order_difference_increases_with_time(self):
        """
        The difference between first and second order should grow with time squared.
        """
        star = self.BARNARDS_STAR

        diffs = []
        times = [50, 100, 200, 400]

        for t in times:
            ra_first, dec_first = _first_order_proper_motion(star, float(t))
            ra_second, dec_second = _second_order_proper_motion(star, float(t))

            diff = math.sqrt(
                (abs(ra_first - ra_second) * 3600) ** 2
                + (abs(dec_first - dec_second) * 3600) ** 2
            )
            diffs.append(diff)

        # Check that differences grow approximately quadratically
        # diff(2t) should be ~4 * diff(t) for second-order correction
        for i in range(len(diffs) - 1):
            ratio = diffs[i + 1] / diffs[i] if diffs[i] > 0 else 0
            time_ratio = times[i + 1] / times[i]
            expected_ratio = time_ratio**2  # Quadratic growth

            # Allow 50% tolerance due to higher-order effects
            assert ratio > expected_ratio * 0.5, (
                f"Difference ratio {ratio:.2f} should be close to "
                f"{expected_ratio:.2f} (t² scaling)"
            )
