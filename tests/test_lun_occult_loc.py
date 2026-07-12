"""
Tests for lunar occultation location-specific calculations in libephemeris.

Tests the lun_occult_when_loc function which finds lunar occultations
of planets and fixed stars visible from a specific geographic location.

Lunar occultations occur when the Moon passes in front of a planet or star.
This function adds location-specific visibility checking to ensure the
occultation is observable from the given coordinates.

The function matches the reference ephemeris lun_occult_when_loc() API:
    lun_occult_when_loc(tjdut, body, geopos, flags, backwards)
    -> (retflags, tret, attr)
"""

import pytest

pytestmark = pytest.mark.slow

from libephemeris import (
    julday,
    lun_occult_when_loc,
    FLG_SWIEPH,
    ECL_TOTAL,
    ECL_PARTIAL,
    ECL_VISIBLE,
    ECL_MAX_VISIBLE,
    ECL_1ST_VISIBLE,
    ECL_4TH_VISIBLE,
)


class TestSweLunOccultWhenLoc:
    """Test suite for lun_occult_when_loc function (reference ephemeris API)."""

    def test_finds_star_occultation_regulus_from_location(self):
        """Test finding a lunar occultation of Regulus visible from a location.

        Regulus (Alpha Leonis) is frequently occulted by the Moon
        because it lies very close to the ecliptic.
        """
        # Start from Jan 1, 2017 - known occultation series in 2017
        jd_start = julday(2017, 1, 1, 0)

        # Search from a location in Europe (central location for good visibility)
        # geopos = [lon, lat, alt] - reference ephemeris convention
        geopos = [10.0, 45.0, 0.0]  # Northern Italy (lon=10°E, lat=45°N)

        retflags, tret, attr = lun_occult_when_loc(
            jd_start, "Regulus", geopos, FLG_SWIEPH
        )

        # Should find an occultation
        assert retflags != 0
        assert tret[0] > jd_start  # Maximum should be after start
        assert tret[0] > 0  # Should have valid maximum time

    def test_returns_correct_tuple_structure(self):
        """Test that return values have correct structure (reference-API compatible)."""
        jd_start = julday(2017, 1, 1, 0)
        geopos = [10.0, 45.0, 0.0]

        retflags, tret, attr = lun_occult_when_loc(
            jd_start, "Regulus", geopos, FLG_SWIEPH
        )

        # retflags should be int
        assert isinstance(retflags, int)

        # tret should be 10-element tuple
        assert len(tret) == 10
        assert all(isinstance(t, float) for t in tret)

        # attr should be 20-element tuple (reference-API compatible)
        assert len(attr) == 20
        assert all(isinstance(a, float) for a in attr)

    def test_contact_times_ordering(self):
        """Test that contact times are in correct chronological order."""
        jd_start = julday(2017, 1, 1, 0)
        geopos = [10.0, 45.0, 0.0]

        retflags, tret, attr = lun_occult_when_loc(
            jd_start, "Regulus", geopos, FLG_SWIEPH
        )

        jd_max = tret[0]
        jd_first = tret[1]
        jd_second = tret[2]
        jd_third = tret[3]
        jd_fourth = tret[4]

        # First contact should be before maximum
        if jd_first > 0:
            assert jd_first < jd_max

        # Fourth contact should be after maximum
        if jd_fourth > 0:
            assert jd_fourth > jd_max

        # If we have second/third contacts (total occultation)
        if jd_second > 0 and jd_third > 0:
            assert jd_first <= jd_second <= jd_max
            assert jd_max <= jd_third <= jd_fourth

    def test_occultation_type_flags(self):
        """Test that occultation type flags are set correctly."""
        jd_start = julday(2017, 1, 1, 0)
        geopos = [10.0, 45.0, 0.0]

        retflags, tret, attr = lun_occult_when_loc(
            jd_start, "Regulus", geopos, FLG_SWIEPH
        )

        # Should be either total or partial
        assert (retflags & ECL_TOTAL) or (retflags & ECL_PARTIAL)

        # Should have visibility flags
        assert retflags & ECL_VISIBLE

    def test_visibility_flags_set(self):
        """Test that visibility flags are properly set."""
        jd_start = julday(2017, 1, 1, 0)
        geopos = [10.0, 45.0, 0.0]

        retflags, tret, attr = lun_occult_when_loc(
            jd_start, "Regulus", geopos, FLG_SWIEPH
        )

        # At least one visibility flag should be set
        visibility_flags = (
            ECL_VISIBLE | ECL_MAX_VISIBLE | ECL_1ST_VISIBLE | ECL_4TH_VISIBLE
        )
        assert retflags & visibility_flags

    def test_moon_attributes_returned(self):
        """Test that occultation attributes are returned.

        Uses reference-API compatible format.
        """
        jd_start = julday(2017, 1, 1, 0)
        geopos = [10.0, 45.0, 0.0]

        retflags, tret, attr = lun_occult_when_loc(
            jd_start, "Regulus", geopos, FLG_SWIEPH
        )

        # Fraction covered should be between 0 and 1
        assert 0.0 <= attr[0] <= 1.0

        # Target azimuth should be between 0 and 360
        assert 0.0 <= attr[4] < 360.0

        # Target altitude should be between -90 and 90
        assert -90.0 <= attr[5] <= 90.0

        # Angular separation (elongation) at maximum
        # For geocentric occultations, this is very small (<0.5 deg)
        # For global occultations visible only due to parallax, this can be
        # up to ~1.5 degrees (Moon radius + parallax margin)
        assert 0.0 <= attr[7] < 2.0

    def test_body_zero_is_the_sun(self):
        """Body 0 is the Sun: the search finds a solar eclipse.

        The reference treats ipl=0 (Sun) as a legitimate occultation
        target - an occultation of the Sun by the Moon is a solar
        eclipse (verified against the reference ephemeris 2.10.03).
        """
        jd_start = julday(2024, 1, 1, 0)
        geopos = [-96.797, 32.7767, 0.0]  # Dallas, 2024-04-08 path

        retflags, tret, attr = lun_occult_when_loc(jd_start, 0, geopos, FLG_SWIEPH)
        assert retflags != 0
        assert tret[0] > jd_start

    def test_raises_error_for_unknown_star(self):
        """Test that function raises error for unknown star name."""
        jd_start = julday(2017, 1, 1, 0)
        geopos = [10.0, 45.0, 0.0]

        from libephemeris.exceptions import Error

        with pytest.raises(Error):
            lun_occult_when_loc(jd_start, "UnknownStar123", geopos, FLG_SWIEPH)

    def test_raises_error_for_invalid_geopos(self):
        """Test that function raises error for invalid geopos."""
        jd_start = julday(2017, 1, 1, 0)
        geopos = [10.0, 45.0]  # Missing altitude

        with pytest.raises(ValueError):
            lun_occult_when_loc(jd_start, "Regulus", geopos, FLG_SWIEPH)

    def test_occultation_duration_reasonable(self):
        """Test that occultation duration is reasonable.

        Lunar occultations typically last a few seconds to about an hour,
        depending on the geometry.
        """
        jd_start = julday(2017, 1, 1, 0)
        geopos = [10.0, 45.0, 0.0]

        retflags, tret, attr = lun_occult_when_loc(
            jd_start, "Regulus", geopos, FLG_SWIEPH
        )

        jd_first = tret[1]
        jd_fourth = tret[4]

        if jd_first > 0 and jd_fourth > 0:
            duration_hours = (jd_fourth - jd_first) * 24
            # Duration should be between 0.01 hours (~30 seconds) and 2 hours
            assert 0.001 < duration_hours < 2.0


class TestSweLunOccultLocDifferentLocations:
    """Test occultation visibility from different locations."""

    def test_same_occultation_different_visibility(self):
        """Test that the same occultation may be found but with different attributes."""
        jd_start = julday(2017, 1, 1, 0)

        # Two different locations (reference ephemeris geopos = [lon, lat, alt])
        geopos1 = [-0.1, 51.5, 0.0]  # London
        geopos2 = [139.7, 35.7, 0.0]  # Tokyo

        # Both should find occultations
        retflags1, tret1, attr1 = lun_occult_when_loc(
            jd_start, "Regulus", geopos1, FLG_SWIEPH
        )
        retflags2, tret2, attr2 = lun_occult_when_loc(
            jd_start, "Regulus", geopos2, FLG_SWIEPH
        )

        # Both should find occultations
        assert tret1[0] > 0
        assert tret2[0] > 0

        # Target altitude will be different at different locations
        # (even for the same event, if visible from both)
        # Just verify they're valid altitudes
        assert -90.0 <= attr1[5] <= 90.0
        assert -90.0 <= attr2[5] <= 90.0


class TestSweLunOccultLocEdgeCases:
    """Test edge cases for lunar occultation location calculations."""

    def test_altitude_parameter(self):
        """Test that altitude parameter is accepted."""
        jd_start = julday(2017, 1, 1, 0)
        geopos = [10.0, 45.0, 1000.0]  # 1000m above sea level

        # Should not raise an error
        retflags, tret, attr = lun_occult_when_loc(
            jd_start, "Regulus", geopos, FLG_SWIEPH
        )

        assert tret[0] > 0

    def test_polar_region(self):
        """Test occultation search from polar region."""
        jd_start = julday(2017, 1, 1, 0)
        geopos = [25.0, 70.0, 0.0]  # Northern Norway

        # Should find an occultation or raise RuntimeError if none visible
        try:
            retflags, tret, attr = lun_occult_when_loc(
                jd_start, "Regulus", geopos, FLG_SWIEPH
            )
            assert tret[0] > jd_start
        except RuntimeError as e:
            # It's valid if no occultation is visible from this location
            assert "No lunar occultation" in str(e)

    def test_southern_hemisphere(self):
        """Test occultation search from southern hemisphere."""
        jd_start = julday(2017, 1, 1, 0)
        geopos = [18.4, -33.9, 0.0]  # Cape Town, South Africa

        # Should find an occultation or raise RuntimeError if none visible
        try:
            retflags, tret, attr = lun_occult_when_loc(
                jd_start, "Regulus", geopos, FLG_SWIEPH
            )
            assert tret[0] > jd_start
        except RuntimeError as e:
            assert "No lunar occultation" in str(e)

    def test_moonrise_moonset_times(self):
        """Test that moonrise/moonset times are provided when applicable."""
        jd_start = julday(2017, 1, 1, 0)
        geopos = [10.0, 45.0, 0.0]

        retflags, tret, attr = lun_occult_when_loc(
            jd_start, "Regulus", geopos, FLG_SWIEPH
        )

        # Moonrise time is at index 7
        moonrise = tret[7]
        # Moonset time is at index 8
        moonset = tret[8]

        # These should be 0 or valid JD times
        assert moonrise >= 0
        assert moonset >= 0

        # If moonrise is set, it should be within the occultation window
        if moonrise > 0 and tret[1] > 0 and tret[4] > 0:
            assert tret[1] <= moonrise <= tret[4]

        # If moonset is set, it should be within the occultation window
        if moonset > 0 and tret[1] > 0 and tret[4] > 0:
            assert tret[1] <= moonset <= tret[4]

    def test_angular_separation_at_maximum(self):
        """Test that angular separation at maximum is provided."""
        jd_start = julday(2017, 1, 1, 0)
        geopos = [10.0, 45.0, 0.0]

        retflags, tret, attr = lun_occult_when_loc(
            jd_start, "Regulus", geopos, FLG_SWIEPH
        )

        # Angular separation at maximum is at index 7
        min_separation = attr[7]

        # Should be a reasonable value during an occultation
        # Moon's angular radius is about 0.25 degrees, but for global
        # occultations visible due to parallax, separation can be up to
        # ~1.5 degrees (Moon radius + parallax margin)
        assert 0 <= min_separation < 2.0


class TestLunOccultWhenLocLegacy:
    """Test the legacy lun_occult_when_loc function for backward compatibility."""

    def test_legacy_signature_works(self):
        """Test that the legacy function signature still works."""
        jd_start = julday(2017, 1, 1, 0)
        lat, lon = 45.0, 10.0

        ocl_type, times, attr = lun_occult_when_loc(
            jd_start, "Regulus", (lon, lat, 0.0)
        )

        # Should find an occultation
        assert ocl_type != 0
        assert times[0] > jd_start

    def test_legacy_returns_correct_structure(self):
        """Test that legacy function returns (ocl_type, times, attr) order."""
        jd_start = julday(2017, 1, 1, 0)
        lat, lon = 45.0, 10.0

        ocl_type, times, attr = lun_occult_when_loc(
            jd_start, "Regulus", (lon, lat, 0.0)
        )

        # ocl_type is first (int)
        assert isinstance(ocl_type, int)
        # times is second
        assert len(times) == 10
        # attr is third
        assert len(attr) == 20

    def test_legacy_with_altitude(self):
        """Test that legacy function with altitude parameter works."""
        jd_start = julday(2017, 1, 1, 0)
        lat, lon, alt = 45.0, 10.0, 1000.0

        ocl_type, times, attr = lun_occult_when_loc(
            jd_start, "Regulus", (lon, lat, alt)
        )

        assert times[0] > jd_start


class TestIntegerStarIdRoutesToStarPath:
    """Integer fixed-star ids (e.g. REGULUS=1000001) must behave exactly
    like the star name in when_loc/where — they crashed with a raw
    KeyError while when_glob/calc_ut accepted them."""

    def test_when_loc_int_equals_str(self):
        import libephemeris as le

        jd = le.julday(2024, 1, 1, 0)
        r_int = le.lun_occult_when_loc(jd, le.REGULUS, [12.5, 41.9, 0])
        r_str = le.lun_occult_when_loc(jd, "Regulus", [12.5, 41.9, 0])
        assert r_int[0] == r_str[0]
        assert all(abs(a - b) < 1e-9 for a, b in zip(r_int[1], r_str[1]))

    def test_where_int_equals_str(self):
        import libephemeris as le

        w_int = le.lun_occult_where(2460883.3934, le.REGULUS)
        w_str = le.lun_occult_where(2460883.3934, "Regulus")
        assert w_int[0] == w_str[0]
        assert all(abs(a - b) < 1e-9 for a, b in zip(w_int[1], w_str[1]))
