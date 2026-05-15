"""
Comprehensive tests for Local Sidereal Time calculation (sidtime).

Sidereal time is the hour angle of the vernal equinox and is used to
determine which celestial objects are visible at a given time and location.
The function returns local apparent sidereal time in hours (0-24).
"""

import pytest
import swisseph as swe
import libephemeris as ephem


class TestSidtimeBasicValues:
    """Test sidereal time returns correct values for known epochs."""

    @pytest.mark.unit
    def test_sidtime_j2000_greenwich(self):
        """Sidereal time at J2000.0 (Jan 1, 2000 noon) at Greenwich should be ~18.7h."""
        jd = 2451545.0  # J2000
        obliquity = 23.4393
        lst = ephem.sidtime(jd, 0.0, obliquity, 0.0)
        # GMST at J2000.0 is about 18.6974 hours
        assert 18.5 < lst < 19.0, f"LST at J2000 = {lst:.4f} hours, expected ~18.7h"

    @pytest.mark.unit
    def test_sidtime_returns_float(self):
        """sidtime should return a float."""
        jd = 2451545.0
        lst = ephem.sidtime(jd, 0.0, 23.44, 0.0)
        assert isinstance(lst, float)

    @pytest.mark.unit
    def test_sidtime_range_0_to_24(self):
        """sidtime should always return a value between 0 and 24."""
        for offset in range(0, 365, 10):
            jd = 2451545.0 + offset
            lst = ephem.sidtime(jd, 0.0, 23.44, 0.0)
            assert 0 <= lst < 24, f"LST = {lst} hours at JD {jd}, out of range"

    @pytest.mark.unit
    def test_sidtime_progresses_throughout_day(self):
        """Sidereal time should progress by ~24h in a solar day."""
        jd_start = 2451545.0
        jd_end = 2451545.0 + 1.0  # One solar day later

        lst_start = ephem.sidtime(jd_start, 0.0, 23.44, 0.0)
        lst_end = ephem.sidtime(jd_end, 0.0, 23.44, 0.0)

        # Sidereal day is shorter than solar day, so after 1 solar day,
        # sidereal time advances by ~24.066 hours (~4 minutes extra)
        # This means the difference should be about 0.066 hours
        diff = (lst_end - lst_start) % 24
        if diff > 12:  # Handle wrap-around
            diff = 24 - diff
        assert abs(diff - 0.066) < 0.01, f"Expected ~0.066h difference, got {diff:.4f}h"


class TestSidtimeLongitudeOffset:
    """Test that longitude correctly affects local sidereal time."""

    @pytest.mark.unit
    def test_positive_longitude_advances_time(self):
        """Positive (East) longitude should advance local sidereal time."""
        jd = 2451545.0
        lst_greenwich = ephem.sidtime(jd, 0.0, 23.44, 0.0)
        lst_rome = ephem.sidtime(jd, 12.5, 23.44, 0.0)  # Rome ~12.5°E

        # Rome is 12.5°E = 12.5/15 = 0.833 hours ahead
        expected_diff = 12.5 / 15.0
        actual_diff = (lst_rome - lst_greenwich) % 24
        if actual_diff > 12:
            actual_diff = actual_diff - 24

        assert abs(actual_diff - expected_diff) < 0.001, (
            f"Expected diff = {expected_diff:.4f}h, got {actual_diff:.4f}h"
        )

    @pytest.mark.unit
    def test_negative_longitude_retards_time(self):
        """Negative (West) longitude should retard local sidereal time."""
        jd = 2451545.0
        lst_greenwich = ephem.sidtime(jd, 0.0, 23.44, 0.0)
        lst_nyc = ephem.sidtime(jd, -74.0, 23.44, 0.0)  # NYC ~74°W

        # NYC is 74°W = -74/15 = -4.933 hours behind
        expected_diff = -74.0 / 15.0
        actual_diff = lst_nyc - lst_greenwich
        if actual_diff > 12:
            actual_diff = actual_diff - 24
        if actual_diff < -12:
            actual_diff = actual_diff + 24

        assert abs(actual_diff - expected_diff) < 0.001, (
            f"Expected diff = {expected_diff:.4f}h, got {actual_diff:.4f}h"
        )

    @pytest.mark.unit
    def test_180_degrees_longitude(self):
        """180° longitude should be exactly 12 hours different from Greenwich."""
        jd = 2451545.0
        lst_greenwich = ephem.sidtime(jd, 0.0, 23.44, 0.0)
        lst_pacific = ephem.sidtime(jd, 180.0, 23.44, 0.0)

        expected_diff = 12.0  # 180°/15 = 12 hours
        actual_diff = (lst_pacific - lst_greenwich) % 24
        if actual_diff > 12:
            actual_diff = actual_diff - 24

        assert abs(actual_diff - expected_diff) < 0.001, (
            f"Expected diff = {expected_diff:.4f}h, got {actual_diff:.4f}h"
        )


class TestSidtimeNutationEffect:
    """Test that nutation correctly affects sidereal time."""

    @pytest.mark.unit
    def test_nutation_affects_result(self):
        """Non-zero nutation should change the result."""
        jd = 2451545.0
        lst_no_nut = ephem.sidtime(jd, 0.0, 23.44, 0.0)
        lst_with_nut = ephem.sidtime(jd, 0.0, 23.44, 0.005)  # Typical nutation

        # Nutation effect should be small but measurable
        assert lst_no_nut != lst_with_nut
        diff = abs(lst_with_nut - lst_no_nut) * 3600  # in seconds
        # Typical nutation effect is a few seconds
        assert diff < 10, f"Nutation effect = {diff:.2f}s, seems too large"

    @pytest.mark.unit
    def test_nutation_sign(self):
        """Positive nutation should advance sidereal time (for typical obliquity)."""
        jd = 2451545.0
        obliquity = 23.44
        lst_neg = ephem.sidtime(jd, 0.0, obliquity, -0.005)
        lst_pos = ephem.sidtime(jd, 0.0, obliquity, 0.005)

        # With cos(23.44°) > 0, positive nutation adds to sidereal time
        assert lst_pos > lst_neg


class TestSidtimeVsPyswisseph:
    """Compare sidereal time with pyswisseph sidtime0 function."""

    @pytest.mark.comparison
    def test_sidtime_j2000_matches_swe(self):
        """Sidereal time at J2000 should match pyswisseph sidtime0."""
        jd = 2451545.0
        obliquity = 23.4393
        nutation = 0.0

        lst_lib = ephem.sidtime(jd, 0.0, obliquity, nutation)
        gst_swe = swe.sidtime0(jd, obliquity, nutation)

        # At Greenwich (lon=0), our LST equals GST
        diff_seconds = abs(lst_lib - gst_swe) * 3600
        assert diff_seconds < 1.0, (
            f"lib={lst_lib:.6f}h, swe={gst_swe:.6f}h, diff={diff_seconds:.2f}s"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("offset", [0, 0.25, 0.5, 0.75, 1.0, 10, 100, 365])
    def test_sidtime_various_times_match_swe(self, offset):
        """Sidereal time at various times should match pyswisseph."""
        jd = 2451545.0 + offset
        obliquity = 23.44
        nutation = 0.0

        lst_lib = ephem.sidtime(jd, 0.0, obliquity, nutation)
        gst_swe = swe.sidtime0(jd, obliquity, nutation)

        diff_seconds = abs(lst_lib - gst_swe) * 3600
        assert diff_seconds < 1.0, (
            f"JD {jd}: lib={lst_lib:.6f}h, swe={gst_swe:.6f}h, diff={diff_seconds:.2f}s"
        )

    @pytest.mark.comparison
    def test_sidtime_with_nutation_matches_swe(self):
        """Sidereal time with nutation should match pyswisseph."""
        jd = 2451545.0
        obliquity = 23.4393
        nutation = 0.00478  # Typical nutation value

        lst_lib = ephem.sidtime(jd, 0.0, obliquity, nutation)
        gst_swe = swe.sidtime0(jd, obliquity, nutation)

        diff_seconds = abs(lst_lib - gst_swe) * 3600
        assert diff_seconds < 1.0, (
            f"lib={lst_lib:.6f}h, swe={gst_swe:.6f}h, diff={diff_seconds:.2f}s"
        )

    @pytest.mark.comparison
    def test_sidtime_100_dates_match_swe(self, random_dates_in_de421_range):
        """100 random dates should match pyswisseph sidtime0."""
        dates = random_dates_in_de421_range(100)
        max_diff = 0
        obliquity = 23.44
        nutation = 0.0

        for _, _, _, _, jd in dates:
            lst_lib = ephem.sidtime(jd, 0.0, obliquity, nutation)
            gst_swe = swe.sidtime0(jd, obliquity, nutation)
            diff = abs(lst_lib - gst_swe) * 3600  # in seconds
            max_diff = max(max_diff, diff)
            assert diff < 1.0, f"Diff = {diff:.2f}s at JD {jd}"

        print(f"Max sidereal time difference: {max_diff:.3f} seconds")


class TestSidtimeLocalVsGreenwich:
    """Test the relationship between local and Greenwich sidereal time."""

    @pytest.mark.unit
    def test_lst_equals_gst_plus_longitude(self):
        """LST should equal GST + longitude/15."""
        jd = 2451545.0
        longitude = 45.0  # 45°E = 3 hours
        obliquity = 23.44

        lst = ephem.sidtime(jd, longitude, obliquity, 0.0)
        gst = ephem.sidtime(jd, 0.0, obliquity, 0.0)

        expected_lst = (gst + longitude / 15.0) % 24
        assert abs(lst - expected_lst) < 0.0001, (
            f"LST = {lst:.6f}h, expected {expected_lst:.6f}h"
        )


class TestSidtimeEdgeCases:
    """Test edge cases for sidereal time calculation."""

    @pytest.mark.edge_case
    def test_sidtime_de421_range_start(self):
        """Sidereal time at DE421 range start (1900)."""
        jd = ephem.julday(1900, 1, 1, 12.0)
        lst = ephem.sidtime(jd, 0.0, 23.44, 0.0)
        assert 0 <= lst < 24

    @pytest.mark.edge_case
    def test_sidtime_de421_range_end(self):
        """Sidereal time at DE421 range end (2050)."""
        jd = ephem.julday(2050, 1, 1, 12.0)
        lst = ephem.sidtime(jd, 0.0, 23.44, 0.0)
        assert 0 <= lst < 24

    @pytest.mark.edge_case
    def test_sidtime_extreme_longitude_east(self):
        """Test with extreme east longitude (+180°)."""
        jd = 2451545.0
        lst = ephem.sidtime(jd, 180.0, 23.44, 0.0)
        assert 0 <= lst < 24

    @pytest.mark.edge_case
    def test_sidtime_extreme_longitude_west(self):
        """Test with extreme west longitude (-180°)."""
        jd = 2451545.0
        lst = ephem.sidtime(jd, -180.0, 23.44, 0.0)
        assert 0 <= lst < 24

    @pytest.mark.edge_case
    def test_sidtime_midnight_ut(self):
        """Test at midnight UT (JD ending in .5)."""
        jd = 2451544.5  # Jan 1, 2000, 0h UT
        lst = ephem.sidtime(jd, 0.0, 23.44, 0.0)
        gst = swe.sidtime0(jd, 23.44, 0.0)
        diff_seconds = abs(lst - gst) * 3600
        assert diff_seconds < 1.0


class TestSidtimeOutputFormat:
    """Test output format and units of sidtime."""

    @pytest.mark.unit
    def test_output_is_in_hours(self):
        """sidtime should return value in hours, not radians or degrees."""
        jd = 2451545.0
        lst = ephem.sidtime(jd, 0.0, 23.44, 0.0)
        # Hours should be between 0 and 24
        assert 0 <= lst < 24

    @pytest.mark.unit
    def test_conversion_to_degrees(self):
        """Multiplying by 15 should give degrees (0-360)."""
        jd = 2451545.0
        lst = ephem.sidtime(jd, 0.0, 23.44, 0.0)
        lst_degrees = lst * 15.0
        assert 0 <= lst_degrees < 360

    @pytest.mark.unit
    def test_conversion_to_hms(self):
        """LST can be converted to hours, minutes, seconds."""
        jd = 2451545.0
        lst = ephem.sidtime(jd, 0.0, 23.44, 0.0)

        hours = int(lst)
        minutes_frac = (lst - hours) * 60
        minutes = int(minutes_frac)
        seconds = (minutes_frac - minutes) * 60

        assert 0 <= hours < 24
        assert 0 <= minutes < 60
        assert 0 <= seconds < 60


class TestSidtime0BasicValues:
    """Test sidtime0 returns correct Greenwich Sidereal Time values."""

    @pytest.mark.unit
    def test_sidtime0_j2000_greenwich(self):
        """Sidereal time at J2000.0 should be ~18.7h."""
        jd = 2451545.0  # J2000
        obliquity = 23.4393
        gst = ephem.sidtime0(jd, obliquity, 0.0)
        # GMST at J2000.0 is about 18.6974 hours
        assert 18.5 < gst < 19.0, f"GST at J2000 = {gst:.4f} hours, expected ~18.7h"

    @pytest.mark.unit
    def test_sidtime0_returns_float(self):
        """sidtime0 should return a float."""
        jd = 2451545.0
        gst = ephem.sidtime0(jd, 23.44, 0.0)
        assert isinstance(gst, float)

    @pytest.mark.unit
    def test_sidtime0_range_0_to_24(self):
        """sidtime0 should always return a value between 0 and 24."""
        for offset in range(0, 365, 10):
            jd = 2451545.0 + offset
            gst = ephem.sidtime0(jd, 23.44, 0.0)
            assert 0 <= gst < 24, f"GST = {gst} hours at JD {jd}, out of range"

    @pytest.mark.unit
    def test_sidtime0_equals_sidtime_at_greenwich(self):
        """sidtime0(jd, obl, nut) should equal sidtime(jd, 0.0, obl, nut)."""
        jd = 2451545.0
        obliquity = 23.44
        nutation = 0.005

        gst = ephem.sidtime0(jd, obliquity, nutation)
        lst_greenwich = ephem.sidtime(jd, 0.0, obliquity, nutation)

        assert gst == lst_greenwich, f"sidtime0={gst}, sidtime(lon=0)={lst_greenwich}"


class TestSidtime0VsPyswisseph:
    """Compare sidtime0 with pyswisseph sidtime0 function."""

    @pytest.mark.comparison
    def test_sidtime0_j2000_matches_swe(self):
        """GST at J2000 should match pyswisseph sidtime0."""
        jd = 2451545.0
        obliquity = 23.4393
        nutation = 0.0

        gst_lib = ephem.sidtime0(jd, obliquity, nutation)
        gst_swe = swe.sidtime0(jd, obliquity, nutation)

        diff_seconds = abs(gst_lib - gst_swe) * 3600
        assert diff_seconds < 1.0, (
            f"lib={gst_lib:.6f}h, swe={gst_swe:.6f}h, diff={diff_seconds:.2f}s"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("offset", [0, 0.25, 0.5, 0.75, 1.0, 10, 100, 365])
    def test_sidtime0_various_times_match_swe(self, offset):
        """GST at various times should match pyswisseph."""
        jd = 2451545.0 + offset
        obliquity = 23.44
        nutation = 0.0

        gst_lib = ephem.sidtime0(jd, obliquity, nutation)
        gst_swe = swe.sidtime0(jd, obliquity, nutation)

        diff_seconds = abs(gst_lib - gst_swe) * 3600
        assert diff_seconds < 1.0, (
            f"JD {jd}: lib={gst_lib:.6f}h, swe={gst_swe:.6f}h, diff={diff_seconds:.2f}s"
        )

    @pytest.mark.comparison
    def test_sidtime0_with_nutation_matches_swe(self):
        """GST with nutation should match pyswisseph."""
        jd = 2451545.0
        obliquity = 23.4393
        nutation = 0.00478  # Typical nutation value

        gst_lib = ephem.sidtime0(jd, obliquity, nutation)
        gst_swe = swe.sidtime0(jd, obliquity, nutation)

        diff_seconds = abs(gst_lib - gst_swe) * 3600
        assert diff_seconds < 1.0, (
            f"lib={gst_lib:.6f}h, swe={gst_swe:.6f}h, diff={diff_seconds:.2f}s"
        )

    @pytest.mark.comparison
    def test_sidtime0_100_dates_match_swe(self, random_dates_in_de421_range):
        """100 random dates should match pyswisseph sidtime0."""
        dates = random_dates_in_de421_range(100)
        max_diff = 0
        obliquity = 23.44
        nutation = 0.0

        for _, _, _, _, jd in dates:
            gst_lib = ephem.sidtime0(jd, obliquity, nutation)
            gst_swe = swe.sidtime0(jd, obliquity, nutation)
            diff = abs(gst_lib - gst_swe) * 3600  # in seconds
            max_diff = max(max_diff, diff)
            assert diff < 1.0, f"Diff = {diff:.2f}s at JD {jd}"

        print(f"Max sidtime0 difference: {max_diff:.3f} seconds")


class TestSidtime0EdgeCases:
    """Test edge cases for sidtime0 calculation."""

    @pytest.mark.edge_case
    def test_sidtime0_de421_range_start(self):
        """GST at DE421 range start (1900)."""
        jd = ephem.julday(1900, 1, 1, 12.0)
        gst = ephem.sidtime0(jd, 23.44, 0.0)
        assert 0 <= gst < 24

    @pytest.mark.edge_case
    def test_sidtime0_de421_range_end(self):
        """GST at DE421 range end (2050)."""
        jd = ephem.julday(2050, 1, 1, 12.0)
        gst = ephem.sidtime0(jd, 23.44, 0.0)
        assert 0 <= gst < 24

    @pytest.mark.edge_case
    def test_sidtime0_midnight_ut(self):
        """Test at midnight UT (JD ending in .5)."""
        jd = 2451544.5  # Jan 1, 2000, 0h UT
        gst_lib = ephem.sidtime0(jd, 23.44, 0.0)
        gst_swe = swe.sidtime0(jd, 23.44, 0.0)
        diff_seconds = abs(gst_lib - gst_swe) * 3600
        assert diff_seconds < 1.0
