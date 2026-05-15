"""
Comprehensive tests for LMT to LAT conversion (lmt_to_lat).

LMT (Local Mean Time) is the mean solar time for a specific longitude.
LAT (Local Apparent Time) is the true solar time as shown by a sundial.
The difference between them is the Equation of Time.

This is the inverse operation of lat_to_lmt.
"""

import pytest
import libephemeris as ephem


class TestLmtToLatBasicValues:
    """Test LMT to LAT conversion returns reasonable values."""

    @pytest.mark.unit
    def test_lmt_to_lat_returns_float(self):
        """lmt_to_lat should return a float."""
        jd_lmt = 2451545.0  # J2000
        jd_lat = ephem.lmt_to_lat(jd_lmt, 0.0)
        assert isinstance(jd_lat, float)

    @pytest.mark.unit
    def test_lmt_to_lat_at_greenwich(self):
        """Test LMT to LAT conversion at Greenwich (0 longitude)."""
        jd_lmt = ephem.julday(2000, 6, 15, 12.0)
        jd_lat = ephem.lmt_to_lat(jd_lmt, 0.0)
        # The difference should be approximately the Equation of Time
        # In mid-June, EoT is near zero
        diff_minutes = (jd_lat - jd_lmt) * 1440
        assert -3 < diff_minutes < 3, (
            f"LAT-LMT difference = {diff_minutes:.2f} minutes, expected ~0 in June"
        )

    @pytest.mark.unit
    def test_lmt_to_lat_february_minimum(self):
        """Test LMT to LAT at February minimum (EoT ~-14 min)."""
        # Feb 11, 2000 - EoT minimum
        jd_lmt = ephem.julday(2000, 2, 11, 12.0)
        jd_lat = ephem.lmt_to_lat(jd_lmt, 0.0)
        # EoT is negative (~-14 min), so LAT = LMT + EoT means LAT < LMT
        diff_minutes = (jd_lat - jd_lmt) * 1440
        # EoT = LAT - LMT (when negative, LAT < LMT)
        assert -16 < diff_minutes < -12, (
            f"LAT-LMT difference = {diff_minutes:.2f} minutes, expected ~-14 in Feb"
        )

    @pytest.mark.unit
    def test_lmt_to_lat_november_maximum(self):
        """Test LMT to LAT at November maximum (EoT ~+16 min)."""
        # Nov 3, 2000 - EoT maximum
        jd_lmt = ephem.julday(2000, 11, 3, 12.0)
        jd_lat = ephem.lmt_to_lat(jd_lmt, 0.0)
        # EoT is positive (~+16 min), so LAT = LMT + EoT means LAT > LMT
        diff_minutes = (jd_lat - jd_lmt) * 1440
        assert 14 < diff_minutes < 18, (
            f"LAT-LMT difference = {diff_minutes:.2f} minutes, expected ~+16 in Nov"
        )


class TestLmtToLatProperties:
    """Test properties of LMT to LAT conversion."""

    @pytest.mark.unit
    def test_difference_in_valid_range(self):
        """LAT-LMT difference should be within +/-17 minutes."""
        for month in range(1, 13):
            jd_lmt = ephem.julday(2000, month, 15, 12.0)
            jd_lat = ephem.lmt_to_lat(jd_lmt, 0.0)
            diff_minutes = abs(jd_lat - jd_lmt) * 1440
            assert diff_minutes < 17, (
                f"Month {month}: LAT-LMT = {diff_minutes:.2f} minutes, exceeds 17"
            )

    @pytest.mark.unit
    def test_continuous_over_year(self):
        """LAT-LMT difference should change smoothly over the year."""
        jd_start = ephem.julday(2000, 1, 1, 12.0)
        prev_lat = ephem.lmt_to_lat(jd_start, 0.0)
        prev_diff = prev_lat - jd_start

        for i in range(1, 366):
            jd_lmt = jd_start + i
            jd_lat = ephem.lmt_to_lat(jd_lmt, 0.0)
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
        jd_lmt_2000 = ephem.julday(2000, 6, 15, 12.0)
        jd_lmt_2001 = ephem.julday(2001, 6, 15, 12.0)

        diff_2000 = ephem.lmt_to_lat(jd_lmt_2000, 0.0) - jd_lmt_2000
        diff_2001 = ephem.lmt_to_lat(jd_lmt_2001, 0.0) - jd_lmt_2001

        diff_between_years = abs(diff_2000 - diff_2001) * 1440
        # Should be within 0.5 minutes after exactly one year
        assert diff_between_years < 0.5, (
            f"LAT-LMT differs by {diff_between_years:.3f} minutes between years"
        )


class TestLmtToLatDifferentLongitudes:
    """Test LMT to LAT conversion at different longitudes."""

    @pytest.mark.unit
    def test_east_longitude(self):
        """Test LMT to LAT at eastern longitude (Rome, 12.5E)."""
        jd_lmt = ephem.julday(2000, 6, 15, 12.0)
        jd_lat = ephem.lmt_to_lat(jd_lmt, 12.5)
        # The difference should still be approximately the Equation of Time
        diff_minutes = (jd_lat - jd_lmt) * 1440
        # In June, EoT is near zero
        assert -3 < diff_minutes < 3, (
            f"LAT-LMT at Rome = {diff_minutes:.2f} minutes, expected ~0 in June"
        )

    @pytest.mark.unit
    def test_west_longitude(self):
        """Test LMT to LAT at western longitude (New York, -74W)."""
        jd_lmt = ephem.julday(2000, 6, 15, 12.0)
        jd_lat = ephem.lmt_to_lat(jd_lmt, -74.0)
        # The difference should still be approximately the Equation of Time
        diff_minutes = (jd_lat - jd_lmt) * 1440
        # In June, EoT is near zero
        assert -3 < diff_minutes < 3, (
            f"LAT-LMT at NYC = {diff_minutes:.2f} minutes, expected ~0 in June"
        )

    @pytest.mark.unit
    def test_180_longitude(self):
        """Test LMT to LAT at date line (180 longitude)."""
        jd_lmt = ephem.julday(2000, 6, 15, 12.0)
        jd_lat = ephem.lmt_to_lat(jd_lmt, 180.0)
        diff_minutes = (jd_lat - jd_lmt) * 1440
        # In June, EoT is near zero
        assert -3 < diff_minutes < 3

    @pytest.mark.unit
    def test_minus_180_longitude(self):
        """Test LMT to LAT at date line (-180 longitude)."""
        jd_lmt = ephem.julday(2000, 6, 15, 12.0)
        jd_lat = ephem.lmt_to_lat(jd_lmt, -180.0)
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
        jd_lmt = ephem.julday(2000, 11, 3, 12.0)  # Nov 3 - EoT max
        jd_lat = ephem.lmt_to_lat(jd_lmt, longitude)
        diff_minutes = (jd_lat - jd_lmt) * 1440
        # EoT at Nov 3 is ~+16 minutes, with some variation due to the
        # longitude adjustment in the EoT calculation
        assert 13 < diff_minutes < 19, (
            f"Longitude {longitude}: LAT-LMT = {diff_minutes:.2f} min"
        )


class TestLmtToLatConsistencyWithTimeEqu:
    """Test that LMT to LAT is consistent with time_equ."""

    @pytest.mark.unit
    def test_difference_equals_time_equ_at_greenwich(self):
        """At Greenwich (0), LAT-LMT should equal time_equ."""
        jd_lmt = ephem.julday(2000, 6, 15, 12.0)
        jd_lat = ephem.lmt_to_lat(jd_lmt, 0.0)

        # LAT - LMT = EoT
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
            jd_lmt = ephem.julday(year, month, day, 12.0)
            jd_lat = ephem.lmt_to_lat(jd_lmt, 0.0)
            lat_lmt_diff = jd_lat - jd_lmt
            eot = ephem.time_equ(jd_lmt)
            assert lat_lmt_diff == pytest.approx(eot, abs=1e-8), (
                f"{year}-{month}-{day}: LAT-LMT = {lat_lmt_diff}, EoT = {eot}"
            )


class TestLmtToLatInverseOfLatToLmt:
    """Test that lmt_to_lat is the inverse of lat_to_lmt."""

    @pytest.mark.unit
    def test_roundtrip_lat_to_lmt_to_lat(self):
        """Converting LAT -> LMT -> LAT should return the original value."""
        jd_lat_original = ephem.julday(2000, 6, 15, 12.0)
        longitude = 12.5  # Rome

        # LAT -> LMT
        jd_lmt = ephem.lat_to_lmt(jd_lat_original, longitude)
        # LMT -> LAT
        jd_lat_recovered = ephem.lmt_to_lat(jd_lmt, longitude)

        # Should be very close to original (within ~1 second)
        diff_seconds = abs(jd_lat_recovered - jd_lat_original) * 86400
        assert diff_seconds < 1.0, f"Roundtrip error = {diff_seconds:.3f} seconds"

    @pytest.mark.unit
    def test_roundtrip_lmt_to_lat_to_lmt(self):
        """Converting LMT -> LAT -> LMT should return the original value."""
        jd_lmt_original = ephem.julday(2000, 6, 15, 12.0)
        longitude = -74.0  # New York

        # LMT -> LAT
        jd_lat = ephem.lmt_to_lat(jd_lmt_original, longitude)
        # LAT -> LMT
        jd_lmt_recovered = ephem.lat_to_lmt(jd_lat, longitude)

        # Should be very close to original (within ~1 second)
        diff_seconds = abs(jd_lmt_recovered - jd_lmt_original) * 86400
        assert diff_seconds < 1.0, f"Roundtrip error = {diff_seconds:.3f} seconds"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "date,longitude",
        [
            ((2000, 1, 1, 12.0), 0.0),
            ((2000, 2, 11, 12.0), 12.5),  # EoT minimum
            ((2000, 6, 15, 12.0), -74.0),
            ((2000, 11, 3, 12.0), 180.0),  # EoT maximum
            ((2020, 7, 4, 6.0), -120.0),
        ],
    )
    def test_roundtrip_various_dates_and_longitudes(self, date, longitude):
        """Roundtrip should work for various dates and longitudes."""
        year, month, day, hour = date
        jd_lat_original = ephem.julday(year, month, day, hour)

        # LAT -> LMT -> LAT
        jd_lmt = ephem.lat_to_lmt(jd_lat_original, longitude)
        jd_lat_recovered = ephem.lmt_to_lat(jd_lmt, longitude)

        diff_seconds = abs(jd_lat_recovered - jd_lat_original) * 86400
        assert diff_seconds < 1.0, (
            f"Roundtrip error for {date} @ {longitude} = {diff_seconds:.3f}s"
        )


class TestLmtToLatEdgeCases:
    """Test edge cases for LMT to LAT conversion."""

    @pytest.mark.edge_case
    def test_de421_range_start(self):
        """Test at DE421 range start (1900)."""
        jd_lmt = ephem.julday(1900, 1, 1, 12.0)
        jd_lat = ephem.lmt_to_lat(jd_lmt, 0.0)
        diff_minutes = abs(jd_lat - jd_lmt) * 1440
        # Should be within valid range
        assert diff_minutes < 17

    @pytest.mark.edge_case
    def test_de421_range_end(self):
        """Test at DE421 range end (2050)."""
        jd_lmt = ephem.julday(2050, 1, 1, 12.0)
        jd_lat = ephem.lmt_to_lat(jd_lmt, 0.0)
        diff_minutes = abs(jd_lat - jd_lmt) * 1440
        # Should be within valid range
        assert diff_minutes < 17

    @pytest.mark.edge_case
    def test_various_years(self):
        """LMT to LAT should work for various years."""
        for year in [1920, 1960, 2000, 2020, 2040]:
            jd_lmt = ephem.julday(year, 6, 15, 12.0)
            jd_lat = ephem.lmt_to_lat(jd_lmt, 0.0)
            diff_minutes = abs(jd_lat - jd_lmt) * 1440
            assert diff_minutes < 17, (
                f"Year {year}: LAT-LMT = {diff_minutes:.2f} minutes"
            )

    @pytest.mark.edge_case
    def test_extreme_longitudes(self):
        """Test with extreme longitude values."""
        jd_lmt = ephem.julday(2000, 6, 15, 12.0)

        # 180 degrees
        jd_lat_180 = ephem.lmt_to_lat(jd_lmt, 180.0)
        assert isinstance(jd_lat_180, float)

        # -180 degrees
        jd_lat_neg180 = ephem.lmt_to_lat(jd_lmt, -180.0)
        assert isinstance(jd_lat_neg180, float)


class TestLmtToLatOutputFormat:
    """Test output format of lmt_to_lat."""

    @pytest.mark.unit
    def test_output_is_julian_day(self):
        """lmt_to_lat should return a Julian Day number."""
        jd_lmt = 2451545.0  # J2000
        jd_lat = ephem.lmt_to_lat(jd_lmt, 0.0)
        # Output should be close to input (difference is EoT, ~minutes)
        assert abs(jd_lat - jd_lmt) < 0.02  # Less than ~30 minutes

    @pytest.mark.unit
    def test_output_preserves_scale(self):
        """Output JD should be in the same scale as input."""
        jd_lmt = ephem.julday(2000, 6, 15, 12.0)
        jd_lat = ephem.lmt_to_lat(jd_lmt, 0.0)
        # The difference should be the EoT in days
        diff_days = jd_lat - jd_lmt
        # EoT is at most ~17 minutes = ~0.012 days
        assert abs(diff_days) < 0.015
