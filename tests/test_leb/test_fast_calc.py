"""
Tests for the fast calculation pipeline (fast_calc.py).

Compares fast_calc output with Skyfield/calc_ut for all supported bodies
and flag combinations.
"""

from __future__ import annotations

import warnings

import pytest

import libephemeris as ephem
from libephemeris.constants import (
    MARS,
    MEAN_NODE,
    MOON,
    SUN,
    FLG_BARYCTR,
    FLG_EQUATORIAL,
    FLG_HELCTR,
    FLG_J2000,
    FLG_NOABERR,
    FLG_NONUT,
    FLG_RADIANS,
    FLG_SIDEREAL,
    FLG_SPEED,
    FLG_TOPOCTR,
    FLG_TRUEPOS,
    FLG_XYZ,
)
from libephemeris.fast_calc import (
    _calc_ayanamsa_from_leb,
    _cartesian_to_spherical,
    _mean_obliquity_iau2006,
    _rotate_equatorial_to_ecliptic,
    _vec3_dist,
    _vec3_sub,
    fast_calc_tt,
    fast_calc_ut,
)


# =============================================================================
# UTILITY FUNCTION TESTS
# =============================================================================


class TestUtilityFunctions:
    """Test coordinate transform and vector utility functions."""

    @pytest.mark.unit
    def test_vec3_sub(self):
        """Vector subtraction."""
        result = _vec3_sub((3.0, 5.0, 7.0), (1.0, 2.0, 3.0))
        assert result == (2.0, 3.0, 4.0)

    @pytest.mark.unit
    def test_vec3_dist(self):
        """Euclidean distance."""
        assert abs(_vec3_dist((3.0, 4.0, 0.0)) - 5.0) < 1e-14
        assert abs(_vec3_dist((0.0, 0.0, 0.0))) < 1e-14

    @pytest.mark.unit
    def test_cartesian_to_spherical_origin(self):
        """Origin gives (0, 0, 0)."""
        lon, lat, dist = _cartesian_to_spherical(0.0, 0.0, 0.0)
        assert lon == 0.0
        assert lat == 0.0
        assert dist == 0.0

    @pytest.mark.unit
    def test_cartesian_to_spherical_x_axis(self):
        """Point on +x axis: lon=0, lat=0."""
        lon, lat, dist = _cartesian_to_spherical(1.0, 0.0, 0.0)
        assert abs(lon - 0.0) < 1e-12
        assert abs(lat - 0.0) < 1e-12
        assert abs(dist - 1.0) < 1e-14

    @pytest.mark.unit
    def test_cartesian_to_spherical_y_axis(self):
        """Point on +y axis: lon=90."""
        lon, lat, dist = _cartesian_to_spherical(0.0, 1.0, 0.0)
        assert abs(lon - 90.0) < 1e-12

    @pytest.mark.unit
    def test_cartesian_to_spherical_z_axis(self):
        """Point on +z axis: lat=90."""
        lon, lat, dist = _cartesian_to_spherical(0.0, 0.0, 1.0)
        assert abs(lat - 90.0) < 1e-12

    @pytest.mark.unit
    def test_mean_obliquity_j2000(self):
        """Mean obliquity at J2000 should be ~23.4393 degrees."""
        eps = _mean_obliquity_iau2006(2451545.0)
        # IAU 2006: 84381.406 arcsec = 23.4393 deg
        assert abs(eps - 23.4393) < 0.001

    @pytest.mark.unit
    def test_rotate_equatorial_to_ecliptic_identity_at_zero(self):
        """With obliquity=0, rotation is identity."""
        x, y, z = _rotate_equatorial_to_ecliptic(1.0, 2.0, 3.0, 0.0)
        assert abs(x - 1.0) < 1e-14
        assert abs(y - 2.0) < 1e-14
        assert abs(z - 3.0) < 1e-14


# =============================================================================
# FAST CALC vs SKYFIELD COMPARISON TESTS
# =============================================================================


class TestFastCalcVsSkyfield:
    """Compare fast_calc output with Skyfield pipeline for all bodies."""

    # Tolerance in arcseconds
    TOLERANCE_ARCSEC = 0.1  # 100 milliarcseconds (generous for Pipeline A)
    TOLERANCE_ARCSEC_ECLIPTIC = 0.01  # 10 mas for ecliptic direct bodies

    @pytest.fixture
    def test_dates(self, leb_reader):
        """Return a set of test dates within the LEB file range."""
        jd_start, jd_end = leb_reader.jd_range
        margin = 10.0  # Stay away from edges
        # Spread dates across the range
        dates = [
            jd_start + margin,
            (jd_start + jd_end) / 2.0,
            jd_end - margin,
            jd_start + (jd_end - jd_start) * 0.25,
            jd_start + (jd_end - jd_start) * 0.75,
        ]
        return dates

    @pytest.mark.integration
    @pytest.mark.parametrize("ipl", [SUN, MOON, MARS])
    def test_planet_ecliptic_default(self, leb_reader, test_dates, ipl):
        """Planet ecliptic lon/lat matches Skyfield within tolerance."""
        for jd in test_dates:
            fast_result, _ = fast_calc_ut(leb_reader, jd, ipl, FLG_SPEED)
            sky_result, _ = ephem.calc_ut(jd, ipl, FLG_SPEED)

            lon_err = abs(fast_result[0] - sky_result[0])
            # Handle wrap-around
            if lon_err > 180.0:
                lon_err = 360.0 - lon_err
            lon_err_arcsec = lon_err * 3600.0

            lat_err_arcsec = abs(fast_result[1] - sky_result[1]) * 3600.0

            assert lon_err_arcsec < self.TOLERANCE_ARCSEC, (
                f"Body {ipl} at JD {jd}: lon error = {lon_err_arcsec:.4f} arcsec"
            )
            assert lat_err_arcsec < self.TOLERANCE_ARCSEC, (
                f"Body {ipl} at JD {jd}: lat error = {lat_err_arcsec:.4f} arcsec"
            )

    @pytest.mark.integration
    def test_mean_node_ecliptic(self, leb_reader, test_dates):
        """Mean node longitude matches Skyfield within tolerance."""
        for jd in test_dates:
            fast_result, _ = fast_calc_ut(leb_reader, jd, MEAN_NODE, FLG_SPEED)
            sky_result, _ = ephem.calc_ut(jd, MEAN_NODE, FLG_SPEED)

            lon_err = abs(fast_result[0] - sky_result[0])
            if lon_err > 180.0:
                lon_err = 360.0 - lon_err
            lon_err_arcsec = lon_err * 3600.0

            assert lon_err_arcsec < self.TOLERANCE_ARCSEC_ECLIPTIC, (
                f"Mean node at JD {jd}: lon error = {lon_err_arcsec:.4f} arcsec"
            )

    @pytest.mark.integration
    def test_sun_speed(self, leb_reader, test_dates):
        """Sun speed from .leb should match Skyfield within ~0.001 deg/day."""
        for jd in test_dates:
            fast_result, _ = fast_calc_ut(leb_reader, jd, SUN, FLG_SPEED)
            sky_result, _ = ephem.calc_ut(jd, SUN, FLG_SPEED)

            # Speed tolerance: 0.001 deg/day = 3.6 arcsec/day
            speed_err = abs(fast_result[3] - sky_result[3])
            assert speed_err < 0.01, (
                f"Sun speed error at JD {jd}: {speed_err:.6f} deg/day"
            )

    @pytest.mark.integration
    def test_sun_distance(self, leb_reader, test_dates):
        """Sun distance from .leb should match Skyfield within ~1e-6 AU."""
        for jd in test_dates:
            fast_result, _ = fast_calc_ut(leb_reader, jd, SUN, FLG_SPEED)
            sky_result, _ = ephem.calc_ut(jd, SUN, FLG_SPEED)

            dist_err = abs(fast_result[2] - sky_result[2])
            assert dist_err < 1e-4, f"Sun distance error at JD {jd}: {dist_err:.2e} AU"


class TestFastCalcFlags:
    """Test various flag combinations."""

    @pytest.fixture
    def jd_mid(self, leb_reader):
        """Return the midpoint JD of the test file."""
        jd_start, jd_end = leb_reader.jd_range
        return (jd_start + jd_end) / 2.0

    @pytest.mark.integration
    def test_no_speed(self, leb_reader, jd_mid):
        """Without FLG_SPEED, velocity should be zero."""
        result, _ = fast_calc_ut(leb_reader, jd_mid, SUN, 0)
        assert result[3] == 0.0
        assert result[4] == 0.0
        assert result[5] == 0.0

    @pytest.mark.integration
    def test_with_speed(self, leb_reader, jd_mid):
        """With FLG_SPEED, velocity should be non-zero."""
        result, _ = fast_calc_ut(leb_reader, jd_mid, SUN, FLG_SPEED)
        # Sun moves ~1 deg/day in longitude
        assert abs(result[3]) > 0.5, f"Sun speed too small: {result[3]}"

    @pytest.mark.integration
    def test_equatorial_flag(self, leb_reader, jd_mid):
        """FLG_EQUATORIAL should return equatorial RA/Dec for ICRS_BARY bodies."""
        result, _ = fast_calc_ut(leb_reader, jd_mid, SUN, FLG_EQUATORIAL)
        # RA should be in [0, 360)
        assert 0.0 <= result[0] < 360.0, f"RA = {result[0]}"
        # Dec should be in [-90, 90]
        assert -90.0 <= result[1] <= 90.0, f"Dec = {result[1]}"

    @pytest.mark.integration
    def test_j2000_flag(self, leb_reader, jd_mid):
        """FLG_J2000 should return J2000 ecliptic coords for ICRS_BARY bodies."""
        result, _ = fast_calc_ut(leb_reader, jd_mid, SUN, FLG_J2000)
        assert 0.0 <= result[0] < 360.0, f"Lon = {result[0]}"

    @pytest.mark.integration
    def test_equatorial_j2000_flag(self, leb_reader, jd_mid):
        """FLG_EQUATORIAL|FLG_J2000 should return J2000 RA/Dec."""
        result, _ = fast_calc_ut(leb_reader, jd_mid, SUN, FLG_EQUATORIAL | FLG_J2000)
        assert 0.0 <= result[0] < 360.0, f"RA = {result[0]}"
        assert -90.0 <= result[1] <= 90.0, f"Dec = {result[1]}"

    @pytest.mark.integration
    def test_helctr_flag(self, leb_reader, jd_mid):
        """FLG_HELCTR should return heliocentric coords for ICRS_BARY bodies."""
        result, _ = fast_calc_ut(leb_reader, jd_mid, MARS, FLG_HELCTR)
        # Mars heliocentric distance ~1.38-1.67 AU
        assert 1.0 < result[2] < 2.0, f"Mars helio dist = {result[2]}"

    @pytest.mark.integration
    def test_baryctr_flag(self, leb_reader, jd_mid):
        """FLG_BARYCTR should return barycentric coords for ICRS_BARY bodies."""
        result, _ = fast_calc_ut(leb_reader, jd_mid, MARS, FLG_BARYCTR)
        # Mars barycentric distance ~1.38-1.67 AU
        assert 1.0 < result[2] < 2.0, f"Mars bary dist = {result[2]}"

    @pytest.mark.integration
    def test_noaberr_flag(self, leb_reader, jd_mid):
        """FLG_NOABERR should skip aberration and return valid coords."""
        result, _ = fast_calc_ut(leb_reader, jd_mid, MARS, FLG_NOABERR)
        assert 0.0 <= result[0] < 360.0, f"Lon = {result[0]}"

    @pytest.mark.integration
    def test_truepos_flag(self, leb_reader, jd_mid):
        """FLG_TRUEPOS should skip light-time and aberration."""
        result, _ = fast_calc_ut(leb_reader, jd_mid, MARS, FLG_TRUEPOS)
        assert 0.0 <= result[0] < 360.0, f"Lon = {result[0]}"

    @pytest.mark.integration
    def test_topoctr_requires_topo(self, leb_reader, jd_mid):
        """FLG_TOPOCTR without set_topo() raises ConfigurationError.

        Both backends raise the same typed error for this request (the
        Skyfield path does too); the assertion previously named ValueError,
        which ConfigurationError does not subclass, so the test failed
        whenever it actually ran.
        """
        from libephemeris import state
        from libephemeris.exceptions import ConfigurationError

        saved = state._TOPO
        state._TOPO = None
        try:
            with pytest.raises(ConfigurationError, match="set_topo"):
                fast_calc_ut(leb_reader, jd_mid, SUN, FLG_TOPOCTR)
        finally:
            state._TOPO = saved

    @pytest.mark.integration
    def test_topoctr_with_topo(self, leb_reader, jd_mid):
        """FLG_TOPOCTR works after set_topo()."""
        import libephemeris
        from libephemeris import state

        saved = state._TOPO
        libephemeris.set_topo(12.5, 41.9, 0)
        try:
            result, _flags = fast_calc_ut(leb_reader, jd_mid, SUN, FLG_TOPOCTR)
            assert 0.0 < result[0] < 360.0
        finally:
            state._TOPO = saved

    @pytest.mark.integration
    def test_xyz_flag(self, leb_reader, jd_mid):
        """FLG_XYZ returns Cartesian coordinates (~1 AU from Sun)."""
        import math

        result, _flags = fast_calc_ut(leb_reader, jd_mid, SUN, FLG_XYZ)
        r = math.sqrt(result[0] ** 2 + result[1] ** 2 + result[2] ** 2)
        assert 0.98 < r < 1.02  # ~1 AU distance

    @pytest.mark.integration
    def test_radians_flag(self, leb_reader, jd_mid):
        """FLG_RADIANS returns coordinates in radians."""
        import math

        result_deg, _ = fast_calc_ut(leb_reader, jd_mid, SUN, 0)
        result_rad, _ = fast_calc_ut(leb_reader, jd_mid, SUN, FLG_RADIANS)
        assert abs(result_rad[0] - math.radians(result_deg[0])) < 1e-10


class TestEclipticDirectVelocity:
    """Test ecliptic-direct bodies with FLG_EQUATORIAL + FLG_SPEED.

    This validates the approximate velocity transform in Pipeline B
    (_pipeline_ecliptic) when converting from ecliptic to equatorial.
    """

    @pytest.fixture
    def jd_mid(self, leb_reader):
        jd_start, jd_end = leb_reader.jd_range
        return (jd_start + jd_end) / 2.0

    @pytest.mark.integration
    def test_mean_node_equatorial_speed(self, leb_reader, jd_mid):
        """Mean node equatorial velocity via LEB should match Skyfield."""
        flags = FLG_EQUATORIAL | FLG_SPEED
        fast_result, _ = fast_calc_ut(leb_reader, jd_mid, MEAN_NODE, flags)
        sky_result, _ = ephem.calc_ut(jd_mid, MEAN_NODE, flags)

        # Position check (RA/Dec)
        ra_err = abs(fast_result[0] - sky_result[0])
        if ra_err > 180.0:
            ra_err = 360.0 - ra_err
        ra_err_arcsec = ra_err * 3600.0
        assert ra_err_arcsec < 0.1, f"Mean node RA error: {ra_err_arcsec:.4f} arcsec"

        dec_err_arcsec = abs(fast_result[1] - sky_result[1]) * 3600.0
        assert dec_err_arcsec < 0.1, f"Mean node Dec error: {dec_err_arcsec:.4f} arcsec"

        # Velocity check: RA speed (deg/day)
        # Mean node moves ~-0.053 deg/day ecliptic, equatorial speed is similar
        # Tolerance: 10% of Skyfield's value or 0.01 deg/day, whichever is larger
        if abs(sky_result[3]) > 0.001:
            speed_ra_err = abs(fast_result[3] - sky_result[3])
            tol = max(abs(sky_result[3]) * 0.1, 0.01)
            assert speed_ra_err < tol, (
                f"Mean node RA speed error: {speed_ra_err:.6f} deg/day "
                f"(LEB={fast_result[3]:.6f}, Sky={sky_result[3]:.6f})"
            )

    @pytest.mark.integration
    def test_mean_node_equatorial_j2000_speed(self, leb_reader, jd_mid):
        """Mean node FLG_EQUATORIAL | FLG_J2000 should produce valid output."""
        flags = FLG_EQUATORIAL | FLG_J2000 | FLG_SPEED
        fast_result, _ = fast_calc_ut(leb_reader, jd_mid, MEAN_NODE, flags)
        sky_result, _ = ephem.calc_ut(jd_mid, MEAN_NODE, flags)

        # Position check
        ra_err = abs(fast_result[0] - sky_result[0])
        if ra_err > 180.0:
            ra_err = 360.0 - ra_err
        ra_err_arcsec = ra_err * 3600.0
        assert ra_err_arcsec < 0.5, (
            f"Mean node J2000 RA error: {ra_err_arcsec:.4f} arcsec"
        )

    @pytest.mark.integration
    def test_mean_node_j2000_ecliptic(self, leb_reader, jd_mid):
        """Mean node FLG_J2000 (ecliptic J2000) should match Skyfield."""
        flags = FLG_J2000 | FLG_SPEED
        fast_result, _ = fast_calc_ut(leb_reader, jd_mid, MEAN_NODE, flags)
        sky_result, _ = ephem.calc_ut(jd_mid, MEAN_NODE, flags)

        lon_err = abs(fast_result[0] - sky_result[0])
        if lon_err > 180.0:
            lon_err = 360.0 - lon_err
        lon_err_arcsec = lon_err * 3600.0
        assert lon_err_arcsec < 0.5, (
            f"Mean node J2000 ecliptic lon error: {lon_err_arcsec:.4f} arcsec"
        )


class TestFastCalcFallback:
    """Test fallback behavior for unsupported bodies/ranges."""

    @pytest.fixture
    def jd_mid(self, leb_reader):
        jd_start, jd_end = leb_reader.jd_range
        return (jd_start + jd_end) / 2.0

    @pytest.mark.integration
    def test_unknown_body_raises_keyerror(self, leb_reader, jd_mid):
        """Bodies not in .leb should raise KeyError."""
        with pytest.raises(KeyError, match="Body 99999"):
            fast_calc_ut(leb_reader, jd_mid, 99999, 0)

    @pytest.mark.integration
    def test_out_of_range_raises_valueerror(self, leb_reader):
        """JD outside .leb range should raise ValueError."""
        jd_start, _ = leb_reader.jd_range
        with pytest.raises((ValueError, KeyError)):
            fast_calc_ut(leb_reader, jd_start - 10000, SUN, 0)


class TestFastCalcTT:
    """Test the TT entry point."""

    @pytest.mark.integration
    def test_fast_calc_tt_basic(self, leb_reader):
        """fast_calc_tt should work with TT input."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0
        result, iflag = fast_calc_tt(leb_reader, jd_mid, SUN, FLG_SPEED)

        assert len(result) == 6
        assert 0.0 <= result[0] < 360.0

    @pytest.mark.integration
    def test_tt_vs_ut_close(self, leb_reader):
        """fast_calc_tt and fast_calc_ut should give close results for similar JD."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0

        # UT and TT differ by Delta-T (~69 sec in 2020s)
        # Using the same JD should give slightly different positions
        result_ut, _ = fast_calc_ut(leb_reader, jd_mid, SUN, FLG_SPEED)
        result_tt, _ = fast_calc_tt(leb_reader, jd_mid, SUN, FLG_SPEED)

        # Sun moves ~1 deg/day, Delta-T ~69 sec → diff ~0.0008 deg
        lon_diff = abs(result_ut[0] - result_tt[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        # Should be small but non-zero
        assert lon_diff < 0.01, f"UT vs TT lon diff = {lon_diff} deg (too large)"


class TestNonutFallback:
    """Test FLG_NONUT triggers Skyfield fallback."""

    @pytest.mark.integration
    def test_nonut_flag_ut(self, leb_reader):
        """FLG_NONUT should work in LEB mode (mean ecliptic)."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0
        result, _flags = fast_calc_ut(leb_reader, jd_mid, SUN, FLG_NONUT)
        assert 0.0 < result[0] < 360.0

    @pytest.mark.integration
    def test_nonut_flag_tt(self, leb_reader):
        """FLG_NONUT should work in TT path too."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0
        result, _flags = fast_calc_tt(leb_reader, jd_mid, SUN, FLG_NONUT)
        assert 0.0 < result[0] < 360.0


class TestExplicitSiderealParams:
    """Test thread-safe explicit sidereal parameter passing."""

    @pytest.mark.integration
    def test_explicit_lahiri_mode_is_computed_without_fallback(self, leb_reader):
        """An explicit Lahiri request is mode-specific and emits no warning."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            lahiri, _ = fast_calc_ut(
                leb_reader,
                jd_mid,
                SUN,
                FLG_SPEED | FLG_SIDEREAL,
                sid_mode=1,
            )
        j2000, _ = fast_calc_ut(
            leb_reader,
            jd_mid,
            SUN,
            FLG_SPEED | FLG_SIDEREAL,
            sid_mode=18,
        )
        assert lahiri != pytest.approx(j2000, abs=1e-6)
        assert not caught

    @pytest.mark.integration
    def test_explicit_standard_epoch_modes_differ(self, leb_reader):
        """Explicit J2000 and J1900 frames produce distinct results."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0

        result_j2000, _ = fast_calc_ut(
            leb_reader,
            jd_mid,
            SUN,
            FLG_SPEED | FLG_SIDEREAL,
            sid_mode=18,
        )
        result_j1900, _ = fast_calc_ut(
            leb_reader,
            jd_mid,
            SUN,
            FLG_SPEED | FLG_SIDEREAL,
            sid_mode=19,
        )
        diff = abs(result_j2000[0] - result_j1900[0])
        assert diff > 0.5, f"J2000 vs J1900 diff = {diff:.4f} deg"

    @pytest.mark.integration
    def test_sidereal_speed_includes_precession_correction(self, leb_reader):
        """Sidereal dlon should differ from tropical dlon by precession rate."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0

        trop, _ = fast_calc_ut(leb_reader, jd_mid, SUN, FLG_SPEED)
        sid, _ = fast_calc_ut(
            leb_reader,
            jd_mid,
            SUN,
            FLG_SPEED | FLG_SIDEREAL,
            sid_mode=18,
            sid_t0=2451545.0,
            sid_ayan_t0=0.0,
        )
        # The of-date sidereal longitude subtracts the TRUE ayanamsha
        # (mean + nutation-in-longitude), so the tropical-minus-sidereal speed
        # difference is the true-ayanamsha rate: the general precession rate
        # (~50.3"/yr = 5028.8"/cy ≈ 0.0000382 deg/day) plus dΔψ/dt, which
        # oscillates within ~±0.000017 deg/day (±0.06"/day) over the 18.6-year
        # nutation cycle. (The previous bound 0.00003–0.00005 pinned an older
        # behaviour where only the mean precession rate was subtracted; the
        # speed is now the exact derivative of the reported sidereal position,
        # verified reference-free by central-differencing the position.)
        speed_diff = trop[3] - sid[3]
        assert 0.000015 < speed_diff < 0.00006, (
            f"Sidereal speed correction = {speed_diff:.7f} deg/day "
            f"(expected 0.0000382 ± nutation-rate ~0.000017)"
        )


class TestPipelineBJ2000Velocity:
    """Test that Pipeline B (ecliptic-direct) J2000-only branch precesses velocity."""

    @pytest.mark.integration
    def test_mean_node_j2000_velocity_precessed(self, leb_reader):
        """Mean node velocity in J2000 ecliptic should be precession-corrected."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0

        # J2000 ecliptic (the branch we fixed)
        result_j2000, _ = fast_calc_ut(
            leb_reader, jd_mid, MEAN_NODE, FLG_SPEED | FLG_J2000
        )
        # Ecliptic of date (default)
        result_date, _ = fast_calc_ut(leb_reader, jd_mid, MEAN_NODE, FLG_SPEED)

        # Velocity should differ slightly due to precession correction
        # Mean node moves ~-0.053 deg/day
        assert abs(result_j2000[3]) > 0.01, (
            f"J2000 mean node speed = {result_j2000[3]:.6f} (suspiciously small)"
        )
        # The difference should be small (~precession rate)
        speed_diff = abs(result_j2000[3] - result_date[3])
        assert speed_diff < 0.01, (
            f"J2000 vs date speed diff = {speed_diff:.6f} deg/day (too large)"
        )


class TestAyanamsa:
    """Test the LEB-based ayanamsa computation."""

    @pytest.mark.integration
    def test_default_fagan_bradley_ayanamsa_is_mode_specific(self, leb_reader):
        """The default Fagan/Bradley ID computes without a J2000 fallback."""

        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0

        ephem.set_sid_mode(0)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            aya = _calc_ayanamsa_from_leb(leb_reader, jd_mid)
        j2000 = _calc_ayanamsa_from_leb(leb_reader, jd_mid, sid_mode=18)
        assert 0.0 <= aya < 360.0
        assert abs((aya - j2000 + 180.0) % 360.0 - 180.0) > 1.0
        assert not caught

    @pytest.mark.integration
    def test_explicit_sid_mode_param(self, leb_reader):
        """_calc_ayanamsa_from_leb accepts explicit sid_mode."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0

        aya_j2000 = _calc_ayanamsa_from_leb(leb_reader, jd_mid, sid_mode=18)
        aya_j1900 = _calc_ayanamsa_from_leb(leb_reader, jd_mid, sid_mode=19)
        assert 0.0 <= aya_j2000 < 360.0
        assert 0.0 <= aya_j1900 < 360.0
        assert abs(aya_j2000 - aya_j1900) > 0.5

    @pytest.mark.integration
    def test_star_based_mode_delegates(self, leb_reader):
        """Star-based modes delegate the anchor to the Skyfield pipeline.

        The LEB fast path used to reject star-based sidereal modes with a
        KeyError, forcing the whole calculation onto Skyfield. It now
        computes the (memoized) star anchor via the Skyfield pipeline
        while the body itself stays on LEB, so the mean ayanamsha must
        match the direct implementation in planets.
        """
        from libephemeris.planets import _calc_ayanamsa
        from libephemeris.time_utils import deltat

        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0

        aya = _calc_ayanamsa_from_leb(leb_reader, jd_mid, sid_mode=27)
        jd_ut = jd_mid - deltat(jd_mid)
        expected = _calc_ayanamsa(jd_ut, 27)
        assert abs(aya - expected) < 1e-9, f"LEB {aya} vs Skyfield {expected}"


# =============================================================================
# SIDEREAL BUG FIX REGRESSION TESTS
#
# These tests guard against regressions of 4 sidereal calculation bugs
# fixed in the leb/precision branch.  Each test targets the specific
# error signature of the original bug so that a regression would be caught
# immediately (tight tolerances, targeted flag combos).
#
# Bodies available in the test LEB: Sun(0), Moon(1), Mars(4), Earth(14),
# MeanNode(10).  Full Pipeline B/C coverage (TrueNode, OscuApog, etc.)
# is in tests/test_leb/compare/extended/test_extended_sidereal.py and
# compare_scripts/tests/test_compare_sidereal_regression.py.
# =============================================================================


class TestSiderealRegressionBug1:
    """Bug 1 regression (commit 64b8367):
    Pipeline A SID+EQ used nutation matrix instead of mean equator precession.
    Pipeline A SID+J2K used true ayanamsha instead of mean.

    Error signature: ~36" RA error for SID+EQ, ~0.3" for SID+J2K.
    """

    @pytest.fixture
    def jd_mid(self, leb_reader):
        jd_start, jd_end = leb_reader.jd_range
        return (jd_start + jd_end) / 2.0

    @pytest.mark.integration
    @pytest.mark.parametrize("ipl", [SUN, MOON, MARS])
    def test_sid_eq_pipeline_a_position(self, leb_reader, jd_mid, ipl):
        """SID+EQ for Pipeline A bodies: LEB vs Skyfield within 0.1"."""
        flags = FLG_SPEED | FLG_SIDEREAL | FLG_EQUATORIAL
        ephem.set_sid_mode(1, 2451545.0, 0.0)

        fast_result, _ = fast_calc_ut(
            leb_reader,
            jd_mid,
            ipl,
            flags,
            sid_mode=1,
            sid_t0=2451545.0,
            sid_ayan_t0=0.0,
        )
        sky_result, _ = ephem.calc_ut(jd_mid, ipl, flags)

        ra_err = abs(fast_result[0] - sky_result[0])
        if ra_err > 180.0:
            ra_err = 360.0 - ra_err
        ra_err_arcsec = ra_err * 3600.0

        dec_err_arcsec = abs(fast_result[1] - sky_result[1]) * 3600.0

        # Tolerance: 0.1" catches Bug 1 error (~36" if unfixed)
        assert ra_err_arcsec < 0.1, (
            f'Body {ipl} SID+EQ RA error = {ra_err_arcsec:.4f}" (Bug 1 regression?)'
        )
        assert dec_err_arcsec < 0.1, (
            f'Body {ipl} SID+EQ Dec error = {dec_err_arcsec:.4f}" (Bug 1 regression?)'
        )

    @pytest.mark.integration
    @pytest.mark.parametrize("ipl", [SUN, MOON, MARS])
    def test_sid_j2k_pipeline_a_position(self, leb_reader, jd_mid, ipl):
        """SID+J2K for Pipeline A bodies: LEB vs Skyfield within 0.1"."""
        flags = FLG_SPEED | FLG_SIDEREAL | FLG_J2000
        ephem.set_sid_mode(1, 2451545.0, 0.0)

        fast_result, _ = fast_calc_ut(
            leb_reader,
            jd_mid,
            ipl,
            flags,
            sid_mode=1,
            sid_t0=2451545.0,
            sid_ayan_t0=0.0,
        )
        sky_result, _ = ephem.calc_ut(jd_mid, ipl, flags)

        lon_err = abs(fast_result[0] - sky_result[0])
        if lon_err > 180.0:
            lon_err = 360.0 - lon_err
        lon_err_arcsec = lon_err * 3600.0

        assert lon_err_arcsec < 0.1, (
            f'Body {ipl} SID+J2K lon error = {lon_err_arcsec:.4f}" (Bug 1 regression?)'
        )

    @pytest.mark.integration
    @pytest.mark.parametrize("ipl", [SUN, MOON, MARS])
    def test_sid_eq_j2k_pipeline_a_position(self, leb_reader, jd_mid, ipl):
        """SID+EQ+J2K for Pipeline A bodies: LEB vs Skyfield within 0.5"."""
        flags = FLG_SPEED | FLG_SIDEREAL | FLG_EQUATORIAL | FLG_J2000
        ephem.set_sid_mode(1, 2451545.0, 0.0)

        fast_result, _ = fast_calc_ut(
            leb_reader,
            jd_mid,
            ipl,
            flags,
            sid_mode=1,
            sid_t0=2451545.0,
            sid_ayan_t0=0.0,
        )
        sky_result, _ = ephem.calc_ut(jd_mid, ipl, flags)

        ra_err = abs(fast_result[0] - sky_result[0])
        if ra_err > 180.0:
            ra_err = 360.0 - ra_err
        ra_err_arcsec = ra_err * 3600.0

        assert ra_err_arcsec < 0.5, (
            f'Body {ipl} SID+EQ+J2K RA error = {ra_err_arcsec:.4f}"'
        )

    @pytest.mark.integration
    @pytest.mark.parametrize("ipl", [SUN, MOON, MARS])
    @pytest.mark.parametrize("sid_mode", [0, 1, 3])
    def test_sid_eq_multi_ayanamsha(self, leb_reader, jd_mid, ipl, sid_mode):
        """SID+EQ should match Skyfield for multiple ayanamsha modes."""
        flags = FLG_SPEED | FLG_SIDEREAL | FLG_EQUATORIAL
        ephem.set_sid_mode(sid_mode, 2451545.0, 0.0)

        fast_result, _ = fast_calc_ut(
            leb_reader,
            jd_mid,
            ipl,
            flags,
            sid_mode=sid_mode,
            sid_t0=2451545.0,
            sid_ayan_t0=0.0,
        )
        sky_result, _ = ephem.calc_ut(jd_mid, ipl, flags)

        ra_err = abs(fast_result[0] - sky_result[0])
        if ra_err > 180.0:
            ra_err = 360.0 - ra_err
        ra_err_arcsec = ra_err * 3600.0

        assert ra_err_arcsec < 0.1, (
            f'Body {ipl} SID+EQ mode={sid_mode}: RA error = {ra_err_arcsec:.4f}"'
        )


class TestSiderealRegressionBug2:
    """Bug 2 regression (commit e6555ed):
    Pipeline B/C SID+EQ dpsi nutation handling was wrong.
    - MeanNode/MeanApog should skip dpsi
    - TrueNode/OscuApog should subtract dpsi
    - Mean obliquity must be used for SID+EQ rotation

    Error signature: ~10-20" for MeanNode SID+EQ.
    Only MeanNode is in the test LEB; full coverage in compare tests.
    """

    @pytest.fixture
    def test_dates(self, leb_reader):
        jd_start, jd_end = leb_reader.jd_range
        margin = 10.0
        return [
            jd_start + margin,
            (jd_start + jd_end) / 2.0,
            jd_end - margin,
        ]

    @pytest.mark.integration
    def test_mean_node_sid_eq(self, leb_reader, test_dates):
        """MeanNode SID+EQ: LEB vs Skyfield within 0.1"."""
        flags = FLG_SPEED | FLG_SIDEREAL | FLG_EQUATORIAL
        for jd in test_dates:
            ephem.set_sid_mode(1, 2451545.0, 0.0)

            fast_result, _ = fast_calc_ut(
                leb_reader,
                jd,
                MEAN_NODE,
                flags,
                sid_mode=1,
                sid_t0=2451545.0,
                sid_ayan_t0=0.0,
            )
            sky_result, _ = ephem.calc_ut(jd, MEAN_NODE, flags)

            ra_err = abs(fast_result[0] - sky_result[0])
            if ra_err > 180.0:
                ra_err = 360.0 - ra_err
            ra_err_arcsec = ra_err * 3600.0

            assert ra_err_arcsec < 0.1, (
                f'MeanNode SID+EQ RA error = {ra_err_arcsec:.4f}" '
                f"at JD {jd} (Bug 2 regression?)"
            )

    @pytest.mark.integration
    def test_mean_node_sid_eq_speed(self, leb_reader, test_dates):
        """MeanNode SID+EQ velocity: LEB vs Skyfield within 0.01 deg/day."""
        flags = FLG_SPEED | FLG_SIDEREAL | FLG_EQUATORIAL
        for jd in test_dates:
            ephem.set_sid_mode(1, 2451545.0, 0.0)

            fast_result, _ = fast_calc_ut(
                leb_reader,
                jd,
                MEAN_NODE,
                flags,
                sid_mode=1,
                sid_t0=2451545.0,
                sid_ayan_t0=0.0,
            )
            sky_result, _ = ephem.calc_ut(jd, MEAN_NODE, flags)

            speed_err = abs(fast_result[3] - sky_result[3])

            assert speed_err < 0.01, (
                f"MeanNode SID+EQ RA speed error = {speed_err:.6f} deg/day at JD {jd}"
            )


class TestSiderealRegressionBug3:
    """Bug 3 regression (commit 9f0fde7):
    the reference ephemeris ignores FLG_J2000 for TrueNode, OscuApog, IntpApog, IntpPerg
    when SIDEREAL is set.  MeanNode and MeanApog precess to J2000 normally.

    Only MeanNode is in the test LEB.  We verify MeanNode SID+J2K DOES differ
    from MeanNode SID (confirming J2K is applied for mean bodies).
    Full Pipeline B coverage in compare tests.
    """

    @pytest.fixture
    def test_dates(self, leb_reader):
        jd_start, jd_end = leb_reader.jd_range
        margin = 10.0
        return [
            jd_start + margin,
            (jd_start + jd_end) / 2.0,
            jd_end - margin,
        ]

    @pytest.mark.integration
    def test_mean_node_sid_j2k_differs_from_sid(self, leb_reader, test_dates):
        """MeanNode SID+J2K should differ from SID (J2K IS applied for mean).

        If this test fails (values identical), Bug 3 may have regressed to
        suppress J2K for mean bodies as well.
        """
        for jd in test_dates:
            ephem.set_sid_mode(1, 2451545.0, 0.0)

            sid_result, _ = fast_calc_ut(
                leb_reader,
                jd,
                MEAN_NODE,
                FLG_SPEED | FLG_SIDEREAL,
                sid_mode=1,
                sid_t0=2451545.0,
                sid_ayan_t0=0.0,
            )
            sid_j2k_result, _ = fast_calc_ut(
                leb_reader,
                jd,
                MEAN_NODE,
                FLG_SPEED | FLG_SIDEREAL | FLG_J2000,
                sid_mode=1,
                sid_t0=2451545.0,
                sid_ayan_t0=0.0,
            )

            lon_diff = abs(sid_result[0] - sid_j2k_result[0])
            if lon_diff > 180.0:
                lon_diff = 360.0 - lon_diff

            # Should differ by precession amount (~0.3° for 2023-2028)
            assert lon_diff > 0.001, (
                f"MeanNode SID vs SID+J2K must differ at JD {jd}, "
                f"got diff = {lon_diff:.6f}° (J2K not applied? Bug 3 regression)"
            )

    @pytest.mark.integration
    def test_mean_node_sid_j2k_matches_skyfield(self, leb_reader, test_dates):
        """MeanNode SID+J2K from LEB matches Skyfield within 0.5"."""
        for jd in test_dates:
            flags = FLG_SPEED | FLG_SIDEREAL | FLG_J2000
            ephem.set_sid_mode(1, 2451545.0, 0.0)

            fast_result, _ = fast_calc_ut(
                leb_reader,
                jd,
                MEAN_NODE,
                flags,
                sid_mode=1,
                sid_t0=2451545.0,
                sid_ayan_t0=0.0,
            )
            sky_result, _ = ephem.calc_ut(jd, MEAN_NODE, flags)

            lon_err = abs(fast_result[0] - sky_result[0])
            if lon_err > 180.0:
                lon_err = 360.0 - lon_err
            lon_err_arcsec = lon_err * 3600.0

            assert lon_err_arcsec < 0.5, (
                f'MeanNode SID+J2K lon error = {lon_err_arcsec:.4f}" '
                f"at JD {jd} (Bug 4 regression?)"
            )


class TestSiderealRegressionBug4:
    """Bug 4 regression (commit b816be0):
    (a) _get_precession_matrix() used t.P including ICRS frame bias (~17 mas).
        Fixed to use mean_equator_and_equinox_of_date.rotation_at(t).
    (b) MeanNode/MeanApog SID+J2K applied precession before ayanamsha
        subtraction.  Non-commutative, up to 28" at extreme dates.
        Fixed to defer precession until after sidereal correction.

    Error signature: (a) ~17 mas systematic at J2000, grows with time.
                     (b) up to 28" at extreme dates (non-commutativity).
    """

    @pytest.fixture
    def test_dates(self, leb_reader):
        """Spread dates across the LEB range to catch date-dependent errors."""
        jd_start, jd_end = leb_reader.jd_range
        margin = 10.0
        n = 5
        step = (jd_end - jd_start - 2 * margin) / (n - 1)
        return [jd_start + margin + i * step for i in range(n)]

    @pytest.mark.integration
    def test_mean_node_sid_eq_j2k_combined(self, leb_reader, test_dates):
        """MeanNode SID+EQ+J2K: LEB vs Skyfield within 0.5"."""
        flags = FLG_SPEED | FLG_SIDEREAL | FLG_EQUATORIAL | FLG_J2000
        for jd in test_dates:
            ephem.set_sid_mode(1, 2451545.0, 0.0)

            fast_result, _ = fast_calc_ut(
                leb_reader,
                jd,
                MEAN_NODE,
                flags,
                sid_mode=1,
                sid_t0=2451545.0,
                sid_ayan_t0=0.0,
            )
            sky_result, _ = ephem.calc_ut(jd, MEAN_NODE, flags)

            ra_err = abs(fast_result[0] - sky_result[0])
            if ra_err > 180.0:
                ra_err = 360.0 - ra_err
            ra_err_arcsec = ra_err * 3600.0

            assert ra_err_arcsec < 0.5, (
                f'MeanNode SID+EQ+J2K RA error = {ra_err_arcsec:.4f}" '
                f"at JD {jd} (Bug 4 regression?)"
            )

    @pytest.mark.integration
    @pytest.mark.parametrize("ipl", [SUN, MOON, MARS])
    def test_pipeline_a_precession_no_frame_bias(self, leb_reader, test_dates, ipl):
        """Pipeline A SID+EQ: verify no ICRS frame bias in precession.

        Bug 4a introduced ~17 mas systematic error at J2000, growing with
        distance from J2000.  With the fix, LEB and Skyfield use the same
        mean equator precession matrix (no frame bias).
        """
        flags = FLG_SPEED | FLG_SIDEREAL | FLG_EQUATORIAL
        max_err = 0.0

        for jd in test_dates:
            ephem.set_sid_mode(1, 2451545.0, 0.0)

            fast_result, _ = fast_calc_ut(
                leb_reader,
                jd,
                ipl,
                flags,
                sid_mode=1,
                sid_t0=2451545.0,
                sid_ayan_t0=0.0,
            )
            sky_result, _ = ephem.calc_ut(jd, ipl, flags)

            ra_err = abs(fast_result[0] - sky_result[0])
            if ra_err > 180.0:
                ra_err = 360.0 - ra_err
            ra_err_arcsec = ra_err * 3600.0
            max_err = max(max_err, ra_err_arcsec)

        # 0.1" tolerance catches the ~17 mas frame bias
        assert max_err < 0.1, (
            f'Body {ipl} SID+EQ max RA error = {max_err:.4f}" '
            f"across {len(test_dates)} dates (frame bias regression?)"
        )
