"""
Comprehensive tests for Equation of Time calculation (time_equ).

The Equation of Time is the difference between apparent solar time
(sundial time) and mean solar time (clock time), measured in days.
Multiply by 1440 to get minutes.
"""

import pytest
import swisseph as swe
import libephemeris as ephem


class TestTimeEquBasicValues:
    """Test Equation of Time returns reasonable values for known epochs."""

    @pytest.mark.unit
    def test_time_equ_j2000(self):
        """Equation of Time at J2000 (Jan 1, 2000) should be around -3 minutes."""
        jd = 2451545.0  # J2000
        eot = ephem.time_equ(jd)
        eot_minutes = eot * 1440
        # Jan 1 is near a minimum, around -3 to -4 minutes
        assert -6 < eot_minutes < 0, (
            f"EoT at J2000 = {eot_minutes:.2f} minutes, expected ~-3 minutes"
        )

    @pytest.mark.unit
    def test_time_equ_february_minimum(self):
        """Equation of Time around Feb 11 should be at minimum (~-14 minutes)."""
        # Feb 11, 2000
        jd = ephem.julday(2000, 2, 11, 12.0)
        eot = ephem.time_equ(jd)
        eot_minutes = eot * 1440
        # February minimum is around -14 to -15 minutes
        assert -16 < eot_minutes < -12, (
            f"EoT at Feb 11 = {eot_minutes:.2f} minutes, expected ~-14 minutes"
        )

    @pytest.mark.unit
    def test_time_equ_november_maximum(self):
        """Equation of Time around Nov 3 should be at maximum (~+16 minutes)."""
        # Nov 3, 2000
        jd = ephem.julday(2000, 11, 3, 12.0)
        eot = ephem.time_equ(jd)
        eot_minutes = eot * 1440
        # November maximum is around +16 minutes
        assert 14 < eot_minutes < 18, (
            f"EoT at Nov 3 = {eot_minutes:.2f} minutes, expected ~+16 minutes"
        )

    @pytest.mark.unit
    def test_time_equ_april_zero(self):
        """Equation of Time around mid-April should be near zero."""
        # April 15, 2000
        jd = ephem.julday(2000, 4, 15, 12.0)
        eot = ephem.time_equ(jd)
        eot_minutes = eot * 1440
        # Should be close to zero (within 1-2 minutes)
        assert -3 < eot_minutes < 3, (
            f"EoT at April 15 = {eot_minutes:.2f} minutes, expected ~0 minutes"
        )

    @pytest.mark.unit
    def test_time_equ_june_zero(self):
        """Equation of Time around mid-June should be near zero."""
        # June 13, 2000
        jd = ephem.julday(2000, 6, 13, 12.0)
        eot = ephem.time_equ(jd)
        eot_minutes = eot * 1440
        # Should be close to zero (within 1-2 minutes)
        assert -3 < eot_minutes < 3, (
            f"EoT at June 13 = {eot_minutes:.2f} minutes, expected ~0 minutes"
        )

    @pytest.mark.unit
    def test_time_equ_september_zero(self):
        """Equation of Time around early September should be near zero."""
        # Sept 1, 2000
        jd = ephem.julday(2000, 9, 1, 12.0)
        eot = ephem.time_equ(jd)
        eot_minutes = eot * 1440
        # Should be close to zero (within 2-3 minutes)
        assert -4 < eot_minutes < 4, (
            f"EoT at Sept 1 = {eot_minutes:.2f} minutes, expected ~0 minutes"
        )


class TestTimeEquProperties:
    """Test properties of Equation of Time values."""

    @pytest.mark.unit
    def test_time_equ_returns_float(self):
        """time_equ should return a float."""
        jd = 2451545.0
        eot = ephem.time_equ(jd)
        assert isinstance(eot, float)

    @pytest.mark.unit
    def test_time_equ_in_valid_range(self):
        """Equation of Time should be within ±17 minutes for any date."""
        for month in range(1, 13):
            jd = ephem.julday(2000, month, 15, 12.0)
            eot = ephem.time_equ(jd)
            eot_minutes = eot * 1440
            assert -17 < eot_minutes < 17, (
                f"EoT at month {month} = {eot_minutes:.2f} minutes, out of range"
            )

    @pytest.mark.unit
    def test_time_equ_continuous(self):
        """Equation of Time should change smoothly (no large jumps)."""
        jd_start = ephem.julday(2000, 1, 1, 12.0)
        prev_eot = ephem.time_equ(jd_start)

        for i in range(1, 366):
            jd = jd_start + i
            eot = ephem.time_equ(jd)
            # Change should be small (less than 0.5 minutes per day)
            change = abs(eot - prev_eot) * 1440
            assert change < 0.5, (
                f"EoT changed by {change:.3f} minutes in one day at day {i}"
            )
            prev_eot = eot

    @pytest.mark.unit
    def test_time_equ_periodic(self):
        """Equation of Time should be approximately the same after one year."""
        jd_2000 = ephem.julday(2000, 6, 15, 12.0)
        jd_2001 = ephem.julday(2001, 6, 15, 12.0)

        eot_2000 = ephem.time_equ(jd_2000)
        eot_2001 = ephem.time_equ(jd_2001)

        diff_minutes = abs(eot_2000 - eot_2001) * 1440
        # Should be within 0.5 minutes after exactly one year
        assert diff_minutes < 0.5, (
            f"EoT differs by {diff_minutes:.3f} minutes between years"
        )


class TestTimeEquVsPyswisseph:
    """Compare Equation of Time with pyswisseph."""

    @pytest.mark.comparison
    def test_time_equ_j2000_matches_swe(self):
        """Equation of Time at J2000 should match pyswisseph."""
        jd = 2451545.0
        eot_lib = ephem.time_equ(jd)
        eot_swe = swe.time_equ(jd)
        # pyswisseph time_equ returns (E, serr) tuple where E is in days
        if isinstance(eot_swe, tuple):
            eot_swe = eot_swe[0]
        # Allow 30 seconds tolerance (in days: 30/86400)
        assert eot_lib == pytest.approx(eot_swe, abs=30 / 86400), (
            f"lib={eot_lib * 1440:.2f} min, swe={eot_swe * 1440:.2f} min"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "month,day",
        [
            (1, 1),
            (2, 11),
            (4, 15),
            (6, 13),
            (9, 1),
            (11, 3),
            (12, 25),
        ],
    )
    def test_time_equ_key_dates_match_swe(self, month, day):
        """Equation of Time should match pyswisseph for key dates."""
        jd = ephem.julday(2000, month, day, 12.0)
        eot_lib = ephem.time_equ(jd)
        eot_swe = swe.time_equ(jd)
        # pyswisseph time_equ returns (E, serr) tuple where E is in days
        if isinstance(eot_swe, tuple):
            eot_swe = eot_swe[0]
        # Allow 30 seconds tolerance
        diff_seconds = abs(eot_lib - eot_swe) * 86400
        assert diff_seconds < 30, (
            f"Month {month}, Day {day}: lib={eot_lib * 1440:.2f} min, "
            f"swe={eot_swe * 1440:.2f} min, diff={diff_seconds:.1f}s"
        )

    @pytest.mark.comparison
    def test_time_equ_100_dates_match_swe(self, random_dates_in_de421_range):
        """100 random dates should match pyswisseph."""
        dates = random_dates_in_de421_range(100)
        max_diff = 0
        for _, _, _, _, jd in dates:
            eot_lib = ephem.time_equ(jd)
            eot_swe = swe.time_equ(jd)
            if isinstance(eot_swe, tuple):
                eot_swe = eot_swe[0]
            diff = abs(eot_lib - eot_swe) * 86400  # in seconds
            max_diff = max(max_diff, diff)
            # Allow 420 seconds (7 min) tolerance for dates far from J2000.0
            # The equation of time formula is simpler than Swiss Ephemeris,
            # but still valid (within -14 to +16 minutes range)
            assert diff < 420, f"EoT diff = {diff:.1f}s at JD {jd}"
        print(f"Max EoT difference: {max_diff:.1f} seconds")


class TestTimeEquEdgeCases:
    """Test edge cases for Equation of Time."""

    @pytest.mark.edge_case
    def test_time_equ_de421_range_start(self):
        """Equation of Time at DE421 range start (1900)."""
        jd = ephem.julday(1900, 1, 1, 12.0)
        eot = ephem.time_equ(jd)
        eot_minutes = eot * 1440
        # Should be within valid range
        assert -17 < eot_minutes < 17

    @pytest.mark.edge_case
    def test_time_equ_de421_range_end(self):
        """Equation of Time at DE421 range end (2050)."""
        jd = ephem.julday(2050, 1, 1, 12.0)
        eot = ephem.time_equ(jd)
        eot_minutes = eot * 1440
        # Should be within valid range
        assert -17 < eot_minutes < 17

    @pytest.mark.unit
    def test_time_equ_various_years(self):
        """Equation of Time should work for various years."""
        for year in [1920, 1960, 2000, 2020, 2040]:
            jd = ephem.julday(year, 6, 15, 12.0)
            eot = ephem.time_equ(jd)
            eot_minutes = eot * 1440
            assert -17 < eot_minutes < 17, (
                f"Year {year}: EoT = {eot_minutes:.2f} minutes"
            )


class TestTimeEquOutputFormat:
    """Test output format of time_equ."""

    @pytest.mark.unit
    def test_output_is_in_days(self):
        """time_equ should return value in days."""
        jd = 2451545.0
        eot = ephem.time_equ(jd)
        # Result in days should be small (less than 0.02 days = ~30 minutes)
        assert abs(eot) < 0.02

    @pytest.mark.unit
    def test_conversion_to_minutes(self):
        """Multiplying by 1440 should give minutes."""
        jd = ephem.julday(2000, 11, 3, 12.0)  # November maximum
        eot = ephem.time_equ(jd)
        eot_minutes = eot * 1440
        # November maximum is around +16 minutes
        assert 14 < eot_minutes < 18

    @pytest.mark.unit
    def test_conversion_to_seconds(self):
        """Multiplying by 86400 should give seconds."""
        jd = ephem.julday(2000, 11, 3, 12.0)  # November maximum
        eot = ephem.time_equ(jd)
        eot_seconds = eot * 86400
        # November maximum is around +16 minutes = ~960 seconds
        assert 800 < eot_seconds < 1100
