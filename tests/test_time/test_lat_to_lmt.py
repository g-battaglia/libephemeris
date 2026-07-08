"""
Comprehensive tests for LAT to LMT conversion (lat_to_lmt).

LAT (Local Apparent Time) is the true solar time as shown by a sundial.
LMT (Local Mean Time) is the mean solar time for a specific longitude.
The difference between them is the Equation of Time.
"""

import pytest
import libephemeris as ephem


class TestLatToLmtBasicValues:
    """Test LAT to LMT conversion returns reasonable values."""

    @pytest.mark.unit
    def test_lat_to_lmt_returns_float(self):
        """lat_to_lmt should return a float."""
        jd_lat = 2451545.0  # J2000
        jd_lmt = ephem.lat_to_lmt(jd_lat, 0.0)
        assert isinstance(jd_lmt, float)

    @pytest.mark.unit
    def test_lat_to_lmt_at_greenwich(self):
        """Test LAT to LMT conversion at Greenwich (0° longitude)."""
        jd_lat = ephem.julday(2000, 6, 15, 12.0)
        jd_lmt = ephem.lat_to_lmt(jd_lat, 0.0)
        # The difference should be approximately the Equation of Time
        # In mid-June, EoT is near zero
        diff_minutes = (jd_lat - jd_lmt) * 1440
        assert -3 < diff_minutes < 3, (
            f"LAT-LMT difference = {diff_minutes:.2f} minutes, expected ~0 in June"
        )

    @pytest.mark.unit
    def test_lat_to_lmt_february_minimum(self):
        """Test LAT to LMT at February minimum (EoT ~-14 min)."""
        # Feb 11, 2000 - EoT minimum
        jd_lat = ephem.julday(2000, 2, 11, 12.0)
        jd_lmt = ephem.lat_to_lmt(jd_lat, 0.0)
        # EoT is negative (~-14 min), so LMT = LAT - EoT means LMT > LAT
        diff_minutes = (jd_lat - jd_lmt) * 1440
        # EoT = LAT - LMT (when negative, LAT < LMT)
        assert -16 < diff_minutes < -12, (
            f"LAT-LMT difference = {diff_minutes:.2f} minutes, expected ~-14 in Feb"
        )

    @pytest.mark.unit
    def test_lat_to_lmt_november_maximum(self):
        """Test LAT to LMT at November maximum (EoT ~+16 min)."""
        # Nov 3, 2000 - EoT maximum
        jd_lat = ephem.julday(2000, 11, 3, 12.0)
        jd_lmt = ephem.lat_to_lmt(jd_lat, 0.0)
        # EoT is positive (~+16 min), so LMT = LAT - EoT means LMT < LAT
        diff_minutes = (jd_lat - jd_lmt) * 1440
        assert 14 < diff_minutes < 18, (
            f"LAT-LMT difference = {diff_minutes:.2f} minutes, expected ~+16 in Nov"
        )


class TestLatToLmtProperties:
    """Test properties of LAT to LMT conversion."""

    @pytest.mark.unit
    def test_difference_in_valid_range(self):
        """LAT-LMT difference should be within ±17 minutes."""
        for month in range(1, 13):
            jd_lat = ephem.julday(2000, month, 15, 12.0)
            jd_lmt = ephem.lat_to_lmt(jd_lat, 0.0)
            diff_minutes = abs(jd_lat - jd_lmt) * 1440
            assert diff_minutes < 17, (
                f"Month {month}: LAT-LMT = {diff_minutes:.2f} minutes, exceeds 17"
            )

    @pytest.mark.unit
    def test_continuous_over_year(self):
        """LAT-LMT difference should change smoothly over the year."""
        jd_start = ephem.julday(2000, 1, 1, 12.0)
        prev_lmt = ephem.lat_to_lmt(jd_start, 0.0)
        prev_diff = jd_start - prev_lmt

        for i in range(1, 366):
            jd_lat = jd_start + i
            jd_lmt = ephem.lat_to_lmt(jd_lat, 0.0)
            diff = jd_lat - jd_lmt
            # Change in EoT should be small (less than 0.5 minutes per day)
            change = abs(diff - prev_diff) * 1440
            assert change < 0.5, (
                f"LAT-LMT changed by {change:.3f} minutes in one day at day {i}"
            )
            prev_diff = diff

    @pytest.mark.unit
    def test_periodic_over_years(self):
        """LAT-LMT difference should be similar after one year."""
        jd_lat_2000 = ephem.julday(2000, 6, 15, 12.0)
        jd_lat_2001 = ephem.julday(2001, 6, 15, 12.0)

        diff_2000 = jd_lat_2000 - ephem.lat_to_lmt(jd_lat_2000, 0.0)
        diff_2001 = jd_lat_2001 - ephem.lat_to_lmt(jd_lat_2001, 0.0)

        diff_between_years = abs(diff_2000 - diff_2001) * 1440
        # Should be within 0.5 minutes after exactly one year
        assert diff_between_years < 0.5, (
            f"LAT-LMT differs by {diff_between_years:.3f} minutes between years"
        )


class TestLatToLmtDifferentLongitudes:
    """Test LAT to LMT conversion at different longitudes."""

    @pytest.mark.unit
    def test_east_longitude(self):
        """Test LAT to LMT at eastern longitude (Rome, 12.5°E)."""
        jd_lat = ephem.julday(2000, 6, 15, 12.0)
        jd_lmt = ephem.lat_to_lmt(jd_lat, 12.5)
        # The difference should still be approximately the Equation of Time
        diff_minutes = (jd_lat - jd_lmt) * 1440
        # In June, EoT is near zero
        assert -3 < diff_minutes < 3, (
            f"LAT-LMT at Rome = {diff_minutes:.2f} minutes, expected ~0 in June"
        )

    @pytest.mark.unit
    def test_west_longitude(self):
        """Test LAT to LMT at western longitude (New York, -74°W)."""
        jd_lat = ephem.julday(2000, 6, 15, 12.0)
        jd_lmt = ephem.lat_to_lmt(jd_lat, -74.0)
        # The difference should still be approximately the Equation of Time
        diff_minutes = (jd_lat - jd_lmt) * 1440
        # In June, EoT is near zero
        assert -3 < diff_minutes < 3, (
            f"LAT-LMT at NYC = {diff_minutes:.2f} minutes, expected ~0 in June"
        )

    @pytest.mark.unit
    def test_180_longitude(self):
        """Test LAT to LMT at date line (180° longitude)."""
        jd_lat = ephem.julday(2000, 6, 15, 12.0)
        jd_lmt = ephem.lat_to_lmt(jd_lat, 180.0)
        diff_minutes = (jd_lat - jd_lmt) * 1440
        # In June, EoT is near zero
        assert -3 < diff_minutes < 3

    @pytest.mark.unit
    def test_minus_180_longitude(self):
        """Test LAT to LMT at date line (-180° longitude)."""
        jd_lat = ephem.julday(2000, 6, 15, 12.0)
        jd_lmt = ephem.lat_to_lmt(jd_lat, -180.0)
        diff_minutes = (jd_lat - jd_lmt) * 1440
        # In June, EoT is near zero
        assert -3 < diff_minutes < 3

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "longitude",
        [
            0.0,
            15.0,
            30.0,
            45.0,
            60.0,
            90.0,
            120.0,
            150.0,
            180.0,
            -15.0,
            -30.0,
            -45.0,
            -60.0,
            -90.0,
            -120.0,
            -150.0,
            -180.0,
        ],
    )
    def test_various_longitudes_same_eot(self, longitude):
        """LAT-LMT difference should be similar regardless of longitude."""
        # The Equation of Time depends on the date, not on the observer's longitude
        jd_lat = ephem.julday(2000, 11, 3, 12.0)  # Nov 3 - EoT max
        jd_lmt = ephem.lat_to_lmt(jd_lat, longitude)
        diff_minutes = (jd_lat - jd_lmt) * 1440
        # EoT at Nov 3 is ~+16 minutes, with some variation due to the
        # longitude adjustment in the EoT calculation
        assert 13 < diff_minutes < 19, (
            f"Longitude {longitude}°: LAT-LMT = {diff_minutes:.2f} min"
        )


class TestLatToLmtConsistencyWithTimeEqu:
    """Test that LAT to LMT is consistent with time_equ."""

    @pytest.mark.unit
    def test_difference_equals_time_equ_at_greenwich(self):
        """At Greenwich (0°), LAT-LMT should equal time_equ."""
        jd_lat = ephem.julday(2000, 6, 15, 12.0)
        jd_lmt = ephem.lat_to_lmt(jd_lat, 0.0)

        # LAT - LMT = EoT, taken at the UT of the LMT instant. lat_to_lmt is the
        # exact iterative inverse of lmt_to_lat, so the Equation of Time must be
        # evaluated at the mean-time (LMT) instant, not at LAT; at Greenwich
        # (0 deg longitude) UT == LMT.
        lat_lmt_diff = jd_lat - jd_lmt
        eot = ephem.time_equ(jd_lmt)

        # Should be very close
        assert lat_lmt_diff == pytest.approx(eot, abs=1e-8), (
            f"LAT-LMT = {lat_lmt_diff}, EoT = {eot}"
        )

    @pytest.mark.unit
    def test_relationship_with_time_equ_various_dates(self):
        """Test LAT-LMT = EoT relationship at various dates."""
        dates = [
            (2000, 1, 1),  # New Year
            (2000, 2, 11),  # EoT minimum
            (2000, 4, 15),  # EoT zero crossing
            (2000, 6, 13),  # EoT zero crossing
            (2000, 9, 1),  # EoT zero crossing
            (2000, 11, 3),  # EoT maximum
            (2000, 12, 25),  # Christmas
        ]
        for year, month, day in dates:
            jd_lat = ephem.julday(year, month, day, 12.0)
            jd_lmt = ephem.lat_to_lmt(jd_lat, 0.0)
            lat_lmt_diff = jd_lat - jd_lmt
            # EoT at the LMT instant (exact inverse of lmt_to_lat); UT == LMT
            # at Greenwich.
            eot = ephem.time_equ(jd_lmt)
            assert lat_lmt_diff == pytest.approx(eot, abs=1e-8), (
                f"{year}-{month}-{day}: LAT-LMT = {lat_lmt_diff}, EoT = {eot}"
            )


class TestLatToLmtEdgeCases:
    """Test edge cases for LAT to LMT conversion."""

    @pytest.mark.edge_case
    def test_de421_range_start(self):
        """Test at DE421 range start (1900)."""
        jd_lat = ephem.julday(1900, 1, 1, 12.0)
        jd_lmt = ephem.lat_to_lmt(jd_lat, 0.0)
        diff_minutes = abs(jd_lat - jd_lmt) * 1440
        # Should be within valid range
        assert diff_minutes < 17

    @pytest.mark.edge_case
    def test_de421_range_end(self):
        """Test at DE421 range end (2050)."""
        jd_lat = ephem.julday(2050, 1, 1, 12.0)
        jd_lmt = ephem.lat_to_lmt(jd_lat, 0.0)
        diff_minutes = abs(jd_lat - jd_lmt) * 1440
        # Should be within valid range
        assert diff_minutes < 17

    @pytest.mark.edge_case
    def test_various_years(self):
        """LAT to LMT should work for various years."""
        for year in [1920, 1960, 2000, 2020, 2040]:
            jd_lat = ephem.julday(year, 6, 15, 12.0)
            jd_lmt = ephem.lat_to_lmt(jd_lat, 0.0)
            diff_minutes = abs(jd_lat - jd_lmt) * 1440
            assert diff_minutes < 17, (
                f"Year {year}: LAT-LMT = {diff_minutes:.2f} minutes"
            )

    @pytest.mark.edge_case
    def test_extreme_longitudes(self):
        """Test with extreme longitude values."""
        jd_lat = ephem.julday(2000, 6, 15, 12.0)

        # 180 degrees
        jd_lmt_180 = ephem.lat_to_lmt(jd_lat, 180.0)
        assert isinstance(jd_lmt_180, float)

        # -180 degrees
        jd_lmt_neg180 = ephem.lat_to_lmt(jd_lat, -180.0)
        assert isinstance(jd_lmt_neg180, float)


class TestLatToLmtOutputFormat:
    """Test output format of lat_to_lmt."""

    @pytest.mark.unit
    def test_output_is_julian_day(self):
        """lat_to_lmt should return a Julian Day number."""
        jd_lat = 2451545.0  # J2000
        jd_lmt = ephem.lat_to_lmt(jd_lat, 0.0)
        # Output should be close to input (difference is EoT, ~minutes)
        assert abs(jd_lmt - jd_lat) < 0.02  # Less than ~30 minutes

    @pytest.mark.unit
    def test_output_preserves_scale(self):
        """Output JD should be in the same scale as input."""
        jd_lat = ephem.julday(2000, 6, 15, 12.0)
        jd_lmt = ephem.lat_to_lmt(jd_lat, 0.0)
        # The difference should be the EoT in days
        diff_days = jd_lat - jd_lmt
        # EoT is at most ~17 minutes = ~0.012 days
        assert abs(diff_days) < 0.015


class TestLatToLmtIterativeInversion:
    """Regression (I2-F3): lat_to_lmt is the exact iterative inverse of
    lmt_to_lat.

    Recovering LMT from LAT needs the Equation of Time at the UT of the LMT
    instant. A single evaluation at UT ~ LAT is off by up to ~0.2s in late
    November, where dE/dt is steepest; three fixed-point iterations bring the
    round trip back to well under a millisecond.
    """

    # Late November (steep dE/dt) plus the Nov / Feb EoT extrema.
    CRITICAL_DATES = [
        (2024, 11, 29),
        (2024, 11, 3),
        (2024, 2, 11),
    ]

    @pytest.mark.unit
    @pytest.mark.parametrize("y,m,d", CRITICAL_DATES)
    @pytest.mark.parametrize("longitude", [0.0, 12.5, -75.0])
    def test_lmt_lat_roundtrip_under_millisecond(self, y, m, d, longitude):
        """lmt -> lat -> lmt must recover the input to < 0.001 s."""
        jd_lmt = ephem.julday(y, m, d, 12.0)
        jd_lat = ephem.lmt_to_lat(jd_lmt, longitude)
        jd_back = ephem.lat_to_lmt(jd_lat, longitude)
        err_seconds = abs(jd_back - jd_lmt) * 86400.0
        assert err_seconds < 0.001, (
            f"round-trip error {err_seconds:.6f}s at "
            f"{y}-{m:02d}-{d:02d} lon={longitude}"
        )

    @pytest.mark.unit
    def test_lat_to_lmt_inverts_late_november_offset(self):
        """A direct forward/back check at the steep-slope date and a nonzero
        longitude, where the single-evaluation inverse was ~0.17s off."""
        jd_lmt = ephem.julday(2024, 11, 29, 0.0)
        jd_lat = ephem.lmt_to_lat(jd_lmt, 0.0)
        jd_back = ephem.lat_to_lmt(jd_lat, 0.0)
        err_seconds = abs(jd_back - jd_lmt) * 86400.0
        assert err_seconds < 0.001, f"round-trip error {err_seconds:.6f}s"
