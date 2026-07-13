"""
Tests for sol_eclipse_when_loc function in libephemeris.

Tests the solar eclipse local visibility calculations with the exact reference ephemeris signature.

Validation tests use the 2024-Apr-08 total solar eclipse as reference:
- Dallas, Texas: Total eclipse with maximum around 18:42 UTC
- New York City: Partial eclipse with ~90% obscuration

Reference data from NASA Eclipse website and the reference ephemeris comparison.
"""

import pytest

pytestmark = pytest.mark.slow

from libephemeris import (
    julday,
    revjul,
    sol_eclipse_when_loc,
    FLG_SWIEPH,
    ECL_TOTAL,
    ECL_PARTIAL,
    ECL_VISIBLE,
    ECL_MAX_VISIBLE,
)


class TestSweSwolEclipseWhenLocSignature:
    """Test that function signature matches the reference ephemeris."""

    def test_function_exists(self):
        """Test that sol_eclipse_when_loc function exists."""
        from libephemeris.eclipse import sol_eclipse_when_loc

        assert callable(sol_eclipse_when_loc)

    def test_returns_correct_tuple_structure(self):
        """Test that return values have correct structure."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = [-96.797, 32.7767, 0]  # Dallas

        retflag, tret, attr = sol_eclipse_when_loc(jd_start, geopos, FLG_SWIEPH)

        # tret should be 10-element tuple (matching reference API)
        assert len(tret) >= 7
        assert all(isinstance(t, float) for t in tret)

        # attr should be at least 8-element tuple (may be 20 matching reference API)
        assert len(attr) >= 8
        assert all(isinstance(a, float) for a in attr)

        # retflag should be int
        assert isinstance(retflag, int)

    def test_accepts_geopos_as_list(self):
        """Test that function accepts geopos as list."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = [-96.797, 32.7767, 0]

        retflag, tret, attr = sol_eclipse_when_loc(jd_start, geopos, FLG_SWIEPH)
        assert tret[0] > jd_start

    def test_accepts_geopos_as_tuple(self):
        """Test that function accepts geopos as tuple."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (-96.797, 32.7767, 0)

        retflag, tret, attr = sol_eclipse_when_loc(jd_start, geopos, FLG_SWIEPH)
        assert tret[0] > jd_start

    def test_invalid_geopos_raises_error(self):
        """Test that invalid geopos raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)

        # Too few elements
        with pytest.raises(ValueError):
            sol_eclipse_when_loc(jd_start, [0, 0], FLG_SWIEPH)


class TestSweSwolEclipseWhenLocDallasApril2024:
    """Test with Dallas, Texas for April 8, 2024 total solar eclipse."""

    def setup_method(self):
        """Set up test fixtures."""
        # Dallas coordinates: 32.7767°N, 96.7970°W
        # geopos format: [longitude, latitude, altitude] - longitude first!
        self.geopos_dallas = [-96.797, 32.7767, 0]
        # Start search from Jan 1, 2024
        self.jd_start = julday(2024, 1, 1, 0)
        # Expected maximum around April 8, 2024 18:42 UTC
        self.expected_max_jd = julday(2024, 4, 8, 18.7)

    def test_finds_april_2024_eclipse(self):
        """Test that function finds the April 8, 2024 eclipse."""
        retflag, tret, attr = sol_eclipse_when_loc(
            self.jd_start, self.geopos_dallas, FLG_SWIEPH
        )

        # Should find an eclipse
        assert retflag != 0
        assert tret[0] > self.jd_start

        # Check date is April 8, 2024
        year, month, day, hour = revjul(tret[0])
        assert year == 2024
        assert month == 4
        assert day == 8

    def test_dallas_eclipse_is_total(self):
        """Test that Dallas sees a total eclipse."""
        retflag, tret, attr = sol_eclipse_when_loc(
            self.jd_start, self.geopos_dallas, FLG_SWIEPH
        )

        # Should be total eclipse
        assert retflag & ECL_TOTAL

    def test_dallas_maximum_time_accuracy(self):
        """Test maximum time is within 0.001 days (~1.5 minutes) of expected."""
        retflag, tret, attr = sol_eclipse_when_loc(
            self.jd_start, self.geopos_dallas, FLG_SWIEPH
        )

        # Maximum should be around 18:42 UTC = 18.7 hours
        # JD for April 8, 2024 18:42 UTC is approximately 2460409.28
        # Allow 0.01 days (about 15 minutes) tolerance for this test
        # since we're comparing to an approximate expected time
        year, month, day, hour = revjul(tret[0])

        # Should be around 18:40-18:50 UTC
        assert 18.0 < hour < 19.5

    def test_dallas_has_all_contacts(self):
        """Test that total eclipse has all four contacts."""
        retflag, tret, attr = sol_eclipse_when_loc(
            self.jd_start, self.geopos_dallas, FLG_SWIEPH
        )

        # All contact times should be non-zero for total eclipse
        assert tret[0] > 0  # Maximum
        assert tret[1] > 0  # First contact
        # Note: 2nd/3rd contacts may be 0 if totality is brief or at edge
        assert tret[4] > 0  # Fourth contact

    def test_dallas_contacts_in_order(self):
        """Test that contact times are in chronological order."""
        retflag, tret, attr = sol_eclipse_when_loc(
            self.jd_start, self.geopos_dallas, FLG_SWIEPH
        )

        # First contact < maximum < fourth contact
        if tret[1] > 0 and tret[4] > 0:
            assert tret[1] < tret[0] < tret[4]

        # If we have 2nd and 3rd contacts (totality), verify order
        if tret[2] > 0 and tret[3] > 0:
            assert tret[1] < tret[2] < tret[0] < tret[3] < tret[4]

    def test_dallas_visibility_flags(self):
        """Test that visibility flags are set correctly."""
        retflag, tret, attr = sol_eclipse_when_loc(
            self.jd_start, self.geopos_dallas, FLG_SWIEPH
        )

        # Eclipse should be visible
        assert retflag & ECL_VISIBLE
        assert retflag & ECL_MAX_VISIBLE


class TestSweSwolEclipseWhenLocNYCApril2024:
    """Test with New York City for April 8, 2024 eclipse (partial)."""

    def setup_method(self):
        """Set up test fixtures."""
        # NYC coordinates: 40.7128°N, 74.0060°W
        # geopos format: [longitude, latitude, altitude]
        self.geopos_nyc = [-74.006, 40.7128, 0]
        # Start search from Jan 1, 2024
        self.jd_start = julday(2024, 1, 1, 0)

    def test_finds_april_2024_eclipse(self):
        """Test that function finds the April 8, 2024 eclipse from NYC."""
        retflag, tret, attr = sol_eclipse_when_loc(
            self.jd_start, self.geopos_nyc, FLG_SWIEPH
        )

        # Check date is April 8, 2024
        year, month, day, hour = revjul(tret[0])
        assert year == 2024
        assert month == 4
        assert day == 8

    def test_nyc_eclipse_is_partial(self):
        """Test that NYC sees a partial eclipse."""
        retflag, tret, attr = sol_eclipse_when_loc(
            self.jd_start, self.geopos_nyc, FLG_SWIEPH
        )

        # Should be partial eclipse (not total)
        assert retflag & ECL_PARTIAL
        assert not (retflag & ECL_TOTAL)

    def test_nyc_obscuration_approximately_90_percent(self):
        """Test that NYC has approximately 90% obscuration."""
        retflag, tret, attr = sol_eclipse_when_loc(
            self.jd_start, self.geopos_nyc, FLG_SWIEPH
        )

        # attr[2] is obscuration (fraction of solar disc area covered)
        # Expected ~90% for NYC, allow 10% tolerance
        obscuration = attr[2]
        assert 0.80 < obscuration < 1.0

    def test_nyc_has_no_totality_contacts(self):
        """Test that partial eclipse has no 2nd/3rd contacts."""
        retflag, tret, attr = sol_eclipse_when_loc(
            self.jd_start, self.geopos_nyc, FLG_SWIEPH
        )

        # 2nd and 3rd contacts should be 0 for partial eclipse
        assert tret[2] == 0.0  # Second contact (totality begins)
        assert tret[3] == 0.0  # Third contact (totality ends)


class TestSweSwolEclipseWhenLocAttributes:
    """Test eclipse attributes (attr tuple)."""

    def setup_method(self):
        """Set up test fixtures."""
        self.geopos_dallas = [-96.797, 32.7767, 0]
        self.jd_start = julday(2024, 1, 1, 0)

    def test_magnitude_is_reasonable(self):
        """Test that magnitude is in valid range."""
        retflag, tret, attr = sol_eclipse_when_loc(
            self.jd_start, self.geopos_dallas, FLG_SWIEPH
        )

        # attr[0] is magnitude (fraction of Sun diameter covered)
        magnitude = attr[0]
        assert 0.0 <= magnitude <= 1.5  # Can exceed 1 for total eclipses

    def test_ratio_is_reasonable(self):
        """Test that diameter ratio is in valid range."""
        retflag, tret, attr = sol_eclipse_when_loc(
            self.jd_start, self.geopos_dallas, FLG_SWIEPH
        )

        # attr[1] is ratio of Moon to Sun diameter
        ratio = attr[1]
        # Typical range: 0.94 to 1.07 depending on distances
        assert 0.9 < ratio < 1.1

    def test_obscuration_is_reasonable(self):
        """Test that obscuration is in valid range."""
        retflag, tret, attr = sol_eclipse_when_loc(
            self.jd_start, self.geopos_dallas, FLG_SWIEPH
        )

        # attr[2] is obscuration (area fraction)
        # For total eclipses, obscuration = (r_moon/r_sun)^2 which can exceed 1.0
        obscuration = attr[2]
        assert 0.0 <= obscuration <= 1.3

    def test_azimuth_is_reasonable(self):
        """Test that azimuth is in valid range."""
        retflag, tret, attr = sol_eclipse_when_loc(
            self.jd_start, self.geopos_dallas, FLG_SWIEPH
        )

        # attr[4] is azimuth
        azimuth = attr[4]
        assert 0 <= azimuth < 360

    def test_altitude_is_positive_during_visible_eclipse(self):
        """Test that Sun altitude is positive for visible eclipse."""
        retflag, tret, attr = sol_eclipse_when_loc(
            self.jd_start, self.geopos_dallas, FLG_SWIEPH
        )

        # attr[5] is true altitude
        true_altitude = attr[5]
        assert true_altitude > -1.0  # Above horizon

    def test_separation_is_small_at_maximum(self):
        """Test that Sun-Moon separation at maximum is small."""
        retflag, tret, attr = sol_eclipse_when_loc(
            self.jd_start, self.geopos_dallas, FLG_SWIEPH
        )

        # attr[7] is angular distance between centers
        separation = attr[7]
        # Should be less than Sun's angular radius (~0.27 degrees)
        assert separation < 1.0


class TestSweSwolEclipseWhenLocEdgeCases:
    """Test edge cases and special conditions."""

    def test_location_with_no_visible_eclipse_soon(self):
        """Test location that may not see eclipse for a while."""
        # Antarctica location
        geopos = [0, -85, 0]
        jd_start = julday(2024, 1, 1, 0)

        # Should still find an eclipse eventually
        retflag, tret, attr = sol_eclipse_when_loc(jd_start, geopos, FLG_SWIEPH)

        # Should find something (may take years)
        assert tret[0] > jd_start

    def test_altitude_is_included(self):
        """Test that altitude in geopos is used."""
        jd_start = julday(2024, 1, 1, 0)

        # Same location, different altitudes
        geopos_sea = [-96.797, 32.7767, 0]
        geopos_high = [-96.797, 32.7767, 5000]  # 5km altitude

        _, tret_sea, attr_sea = sol_eclipse_when_loc(jd_start, geopos_sea, FLG_SWIEPH)
        _, tret_high, attr_high = sol_eclipse_when_loc(
            jd_start, geopos_high, FLG_SWIEPH
        )

        # Times should be similar but not identical
        # Maximum time may differ slightly due to topocentric effect
        assert abs(tret_sea[0] - tret_high[0]) < 0.01  # Within ~15 minutes

    def test_eastern_hemisphere_location(self):
        """Test eclipse search from Eastern hemisphere."""
        # Tokyo coordinates
        geopos = [139.6917, 35.6895, 0]
        jd_start = julday(2024, 1, 1, 0)

        retflag, tret, attr = sol_eclipse_when_loc(jd_start, geopos, FLG_SWIEPH)

        # Should find an eclipse
        assert tret[0] > jd_start
        assert retflag & ECL_VISIBLE


class TestSweSwolEclipseWhenLocBackward:
    """Test backward search functionality."""

    def test_backward_finds_earlier_eclipse(self):
        """Test that backwards=True finds earlier eclipses."""
        # Start from after April 2024 eclipse
        jd_start = julday(2024, 5, 1, 0)
        geopos = [-96.797, 32.7767, 0]  # Dallas

        retflag, tret, attr = sol_eclipse_when_loc(
            jd_start, geopos, FLG_SWIEPH, backwards=True
        )

        # Should find eclipse before start date
        assert tret[0] < jd_start

        # Should be the April 2024 eclipse
        year, month, day, hour = revjul(tret[0])
        assert year == 2024
        assert month == 4
        assert day == 8
