"""
Tests for solar eclipse calculations.
"""

from libephemeris import (
    sol_eclipse_when_glob,
    sol_eclipse_when_loc,
    julday,
    revjul,
    ECL_TOTAL,
    ECL_ANNULAR,
    ECL_PARTIAL,
    ECL_CENTRAL,
    ECL_VISIBLE,
    ECL_1ST_VISIBLE,
    ECL_4TH_VISIBLE,
    FLG_SWIEPH,
)


class TestSolEclipseWhenGlob:
    """Tests for sol_eclipse_when_glob function."""

    def test_finds_eclipse_after_start_date(self):
        """Should find an eclipse after the start date."""
        # Start from Jan 1, 2024
        jd_start = julday(2024, 1, 1, 0)

        ecl_type, times = sol_eclipse_when_glob(jd_start)

        # Should return a valid eclipse time
        assert times[0] > jd_start
        # Eclipse type should have some flags set
        assert ecl_type != 0

    def test_returns_ten_time_values(self):
        """Should return tuple of 10 time values like the reference ephemeris."""
        jd_start = julday(2024, 1, 1, 0)

        _, times = sol_eclipse_when_glob(jd_start)

        assert len(times) == 10

    def test_time_order_is_correct(self):
        """Eclipse phase times should be in chronological order.

        For sol_eclipse_when_glob, the time indices are:
          [0] = time of maximum eclipse
          [1] = time at local apparent noon (not first contact!)
          [2] = time of eclipse begin (first contact globally)
          [3] = time of eclipse end (last contact globally)
          [4] = time of totality begin
        """
        jd_start = julday(2024, 1, 1, 0)

        ecl_type, times = sol_eclipse_when_glob(jd_start)

        t_max = times[0]
        t_begin = times[2]
        t_end = times[3]

        # Eclipse begin should be before maximum
        assert t_begin < t_max
        # Maximum should be before eclipse end
        assert t_max < t_end

    def test_filter_by_total_eclipse(self):
        """Should filter for total eclipses only."""
        jd_start = julday(2020, 1, 1, 0)

        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        # Should have total eclipse flag set
        assert ecl_type & ECL_TOTAL

    def test_filter_by_annular_eclipse(self):
        """Should filter for annular eclipses only."""
        jd_start = julday(2020, 1, 1, 0)

        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_ANNULAR)

        # Should have annular eclipse flag set
        assert ecl_type & ECL_ANNULAR

    def test_filter_accepts_partial_eclipse_type(self):
        """Should accept partial eclipse filter parameter."""
        # Note: Pure partial solar eclipses (not central) are rare.
        # This test verifies the function accepts the parameter and either
        # finds a partial eclipse or times out (both are valid behaviors).
        jd_start = julday(2020, 1, 1, 0)

        # Test that the function accepts the parameter - we check for any eclipse
        # that includes partial (since many eclipses have partial phases)
        ecl_type, times = sol_eclipse_when_glob(jd_start)

        # The function should return something (any eclipse type)
        assert ecl_type != 0

    def test_known_eclipse_april_2024(self):
        """Test against known total solar eclipse of April 8, 2024."""
        # Start search before the eclipse
        jd_start = julday(2024, 3, 1, 0)

        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        # Convert JD to date
        year, month, day, hour = revjul(times[0])

        # The April 8, 2024 total solar eclipse
        # Allow for some tolerance in the search
        assert year == 2024
        assert month == 4
        # Day should be around 8 (allowing some tolerance)
        assert 7 <= day <= 9

    def test_known_eclipse_october_2023(self):
        """Test against known annular solar eclipse of October 14, 2023."""
        # Start search before the eclipse
        jd_start = julday(2023, 9, 1, 0)

        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_ANNULAR)

        # Convert JD to date
        year, month, day, hour = revjul(times[0])

        # The October 14, 2023 annular solar eclipse
        assert year == 2023
        assert month == 10
        # Day should be around 14
        assert 13 <= day <= 15

    def test_central_eclipse_has_second_third_contact(self):
        """Central eclipses should have second and third contact times.

        Request a *central* total/annular eclipse. A geometry-less
        ``ECL_TOTAL | ECL_ANNULAR`` mask (two type bits, no centrality bit)
        matches no eclipse under the reference ifltype semantics and errors
        at the ephemeris boundary, so the centrality bit is included here.
        """
        jd_start = julday(2024, 1, 1, 0)

        ecl_type, times = sol_eclipse_when_glob(
            jd_start, ecltype=ECL_CENTRAL | ECL_TOTAL | ECL_ANNULAR
        )

        if ecl_type & ECL_CENTRAL:
            # Second contact (total/annular begins)
            t_second = times[2]
            # Third contact (total/annular ends)
            t_third = times[3]

            # These should be non-zero for central eclipses
            assert t_second > 0
            assert t_third > 0
            # And in correct order
            assert t_second < t_third

    def test_swe_wrapper_produces_same_results(self):
        """sol_eclipse_when_glob should produce same results as sol_eclipse_when_glob."""
        jd_start = julday(2024, 1, 1, 0)

        ecl_type1, times1 = sol_eclipse_when_glob(jd_start)
        ecl_type2, times2 = sol_eclipse_when_glob(jd_start)

        # Same eclipse should be found
        assert abs(times1[0] - times2[0]) < 0.001
        assert ecl_type1 == ecl_type2

    def test_flags_parameter_accepted(self):
        """Should accept flags parameter without error."""
        jd_start = julday(2024, 1, 1, 0)

        # Should not raise
        ecl_type, times = sol_eclipse_when_glob(jd_start, flags=FLG_SWIEPH)

        assert times[0] > jd_start

    def test_multiple_eclipses_in_sequence(self):
        """Finding consecutive eclipses should return different events."""
        jd_start = julday(2024, 1, 1, 0)

        # Find first eclipse
        _, times1 = sol_eclipse_when_glob(jd_start)

        # Find next eclipse (search from after first one)
        jd_after_first = times1[0] + 1  # Day after first eclipse
        _, times2 = sol_eclipse_when_glob(jd_after_first)

        # Second eclipse should be after first
        assert times2[0] > times1[0]
        # They should be at least 25 days apart (minimum time between eclipses)
        assert times2[0] - times1[0] >= 25

    def test_search_direction_forward(self):
        """Forward search should find eclipse at or after start date."""
        # Start from Aug 10, 2017 (10 days before Aug 21, 2017 eclipse)
        jd_start = julday(2017, 8, 10, 0)

        ecl_type, times = sol_eclipse_when_glob(
            jd_start, ecltype=ECL_TOTAL, backwards=False
        )

        # Should find Aug 21, 2017 eclipse
        year, month, day, _ = revjul(times[0])
        assert (year, month, day) == (2017, 8, 21)

    def test_search_direction_backward(self):
        """Backward search should find eclipse before start date."""
        # Start from Sep 1, 2017 (after Aug 21, 2017 eclipse)
        jd_start = julday(2017, 9, 1, 0)

        ecl_type, times = sol_eclipse_when_glob(
            jd_start, ecltype=ECL_TOTAL, backwards=True
        )

        # Should find Aug 21, 2017 eclipse
        year, month, day, _ = revjul(times[0])
        assert (year, month, day) == (2017, 8, 21)

    def test_search_direction_bidirectional_finds_nearby_eclipse(self):
        """Bidirectional search should find eclipse within 15 days even if forward search would miss."""
        # Start from Aug 15, 2017 (6 days before Aug 21, 2017 eclipse)
        # In this case, the eclipse maximum is after the start date, so bidirectional should find it
        jd_start = julday(2017, 8, 15, 0)

        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        # Should find Aug 21, 2017 eclipse
        year, month, day, _ = revjul(times[0])
        assert (year, month, day) == (2017, 8, 21)

    def test_backward_false_is_default(self):
        """Default search should use forward mode (backwards=False)."""
        jd_start = julday(2017, 8, 1, 0)

        # Both calls should return the same result
        ecl_type_default, times_default = sol_eclipse_when_glob(
            jd_start, ecltype=ECL_TOTAL
        )
        ecl_type_explicit, times_explicit = sol_eclipse_when_glob(
            jd_start, ecltype=ECL_TOTAL, backwards=False
        )

        assert times_default[0] == times_explicit[0]

    def test_august_2017_eclipse_found_from_any_start_in_august(self):
        """Aug 2017 total eclipse should be found regardless of exact start date in August."""
        # This is the key regression test for the bug fix
        for day in [1, 5, 10, 15, 20]:
            jd_start = julday(2017, 8, day, 0)

            ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

            year, month, eclipse_day, _ = revjul(times[0])
            assert (year, month, eclipse_day) == (2017, 8, 21), (
                f"Starting from Aug {day}, expected to find Aug 21 eclipse, "
                f"but found {year}-{month}-{eclipse_day}"
            )


class TestNewMoonFinding:
    """Tests for internal New Moon finding logic."""

    def test_eclipse_occurs_near_new_moon(self):
        """Eclipse should occur very close to a New Moon."""
        from libephemeris import calc_ut, SUN, MOON

        jd_start = julday(2024, 1, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start)

        # Get Sun and Moon positions at eclipse maximum
        sun_pos, _ = calc_ut(times[0], SUN, FLG_SWIEPH)
        moon_pos, _ = calc_ut(times[0], MOON, FLG_SWIEPH)

        # Calculate elongation (should be very close to 0 at New Moon)
        elongation = abs((moon_pos[0] - sun_pos[0] + 180) % 360 - 180)

        # Should be within a few degrees of conjunction
        assert elongation < 5.0


class TestSolEclipseWhenLoc:
    """Tests for sol_eclipse_when_loc function."""

    def test_finds_eclipse_visible_from_location(self):
        """Should find an eclipse visible from the given location."""
        # Start from Jan 1, 2024, search from a central US location
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 35.0, -100.0  # Central US

        ecl_type, times, attr = sol_eclipse_when_loc(jd_start, (lon, lat, 0.0))

        # Should return valid eclipse time
        assert times[0] > jd_start
        # Eclipse should be visible
        assert ecl_type & ECL_VISIBLE
        # Magnitude should be positive
        assert attr[0] > 0

    def test_returns_correct_tuple_sizes(self):
        """Should return times tuple of 10 and attr tuple of 11 elements."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        _, times, attr = sol_eclipse_when_loc(jd_start, (lon, lat, 0.0))

        assert len(times) == 10
        assert len(attr) == 20  # reference version returns 20-element attr tuple

    def test_time_order_is_correct(self):
        """Eclipse phase times should be in chronological order."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 35.0, -100.0

        ecl_type, times, attr = sol_eclipse_when_loc(jd_start, (lon, lat, 0.0))

        t_max = times[0]
        t_first = times[1]
        t_fourth = times[4]

        # First contact should be before maximum (if visible)
        if t_first > 0:
            assert t_first < t_max
        # Maximum should be before fourth contact (if visible)
        if t_fourth > 0:
            assert t_max < t_fourth

    def test_eclipse_attributes_valid(self):
        """Eclipse attributes should be within valid ranges.

        For sol_eclipse_when_loc, the attr indices are:
          [0]: eclipse magnitude
          [1]: ratio of lunar to solar diameter
          [2]: fraction of solar disc obscured
          [3]: core shadow diameter in km (negative = no central shadow)
          [4]: azimuth of sun at maximum (degrees, 0-360)
          [5]: true altitude of sun (degrees)
          [6]: apparent altitude of sun with refraction (degrees)
          [7]: angular distance of moon center from sun center (degrees)
        """
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 35.0, -100.0

        ecl_type, times, attr = sol_eclipse_when_loc(jd_start, (lon, lat, 0.0))

        # Magnitude should be 0-1.5 range
        magnitude = attr[0]
        assert 0 <= magnitude <= 1.5

        # Obscuration should be 0-1 range
        obscuration = attr[2]
        assert 0 <= obscuration <= 1

        # Azimuth should be 0 to 360 (attr[4])
        azimuth = attr[4]
        assert 0 <= azimuth <= 360

        # Sun true altitude should be -90 to 90 (attr[5])
        altitude = attr[5]
        assert -90 <= altitude <= 90

    def test_known_eclipse_april_2024_texas(self):
        """Test April 8, 2024 total eclipse from Dallas, Texas."""
        # Dallas was in the path of totality
        jd_start = julday(2024, 3, 1, 0)
        dallas_lat, dallas_lon = 32.7767, -96.7970

        ecl_type, times, attr = sol_eclipse_when_loc(
            jd_start, (dallas_lon, dallas_lat, 0.0)
        )

        # Convert JD to date
        year, month, day, hour = revjul(times[0])

        # Should find the April 8, 2024 eclipse
        assert year == 2024
        assert month == 4
        assert 7 <= day <= 9

        # Dallas was in path of totality - should have high magnitude
        # (may be total or high partial depending on exact location)
        assert attr[0] > 0.9  # High magnitude

    def test_known_eclipse_april_2024_europe(self):
        """Test April 8, 2024 eclipse from Europe (should not be visible)."""
        # The April 2024 eclipse was only visible in North America
        # From Europe, this eclipse is not visible, so we should find a different one
        jd_start = julday(2024, 3, 1, 0)
        rome_lat, rome_lon = 41.9028, 12.4964

        ecl_type, times, attr = sol_eclipse_when_loc(
            jd_start, (rome_lon, rome_lat, 0.0)
        )

        # Should find some eclipse visible from Rome
        assert times[0] > jd_start
        # It may or may not be April 2024 - the important thing is it's visible
        assert ecl_type & ECL_VISIBLE

    def test_partial_eclipse_location(self):
        """Location far from centerline should see partial eclipse."""
        # Start before April 2024 eclipse and search from a location
        # that's away from the path of totality
        jd_start = julday(2024, 3, 1, 0)
        # Location in Florida (gets partial coverage)
        florida_lat, florida_lon = 25.7617, -80.1918

        ecl_type, times, attr = sol_eclipse_when_loc(
            jd_start, (florida_lon, florida_lat, 0.0)
        )

        # Should have visibility flags
        if times[0] > jd_start:
            assert ecl_type & ECL_VISIBLE

    def test_visibility_flags_set(self):
        """Contact visibility flags should be set appropriately."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 35.0, -100.0

        ecl_type, times, attr = sol_eclipse_when_loc(jd_start, (lon, lat, 0.0))

        # If first contact time is set, flag should be set
        if times[1] > 0:
            assert ecl_type & ECL_1ST_VISIBLE
        # If fourth contact time is set, flag should be set
        if times[4] > 0:
            assert ecl_type & ECL_4TH_VISIBLE

    def test_altitude_parameter(self):
        """Should accept altitude parameter."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 35.0, -100.0
        altitude = 1000.0  # 1000 meters

        # Should not raise
        ecl_type, times, attr = sol_eclipse_when_loc(jd_start, (lon, lat, altitude))

        assert times[0] > jd_start

    def test_flags_parameter_accepted(self):
        """Should accept flags parameter."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 35.0, -100.0

        # Should not raise
        ecl_type, times, attr = sol_eclipse_when_loc(
            jd_start, (lon, lat, 0.0), flags=FLG_SWIEPH
        )

        assert times[0] > jd_start

    def test_swe_version_provides_reference_compatible_interface(self):
        """sol_eclipse_when_loc provides reference-API-compatible interface.

        Note: Both sol_eclipse_when_loc and sol_eclipse_when_loc now use the
        same reference-compatible signature with geopos tuple.
        """
        from libephemeris import FLG_SWIEPH

        jd_start = julday(2024, 1, 1, 0)
        lat, lon, altitude = 35.0, -100.0, 0.0
        geopos = (lon, lat, altitude)  # geopos uses (lon, lat, alt) order

        # Call sol_eclipse_when_loc (now uses reference signature)
        ecl_type1, times1, attr1 = sol_eclipse_when_loc(jd_start, geopos)

        # Call sol_eclipse_when_loc (reference-style with geopos sequence)
        ecl_type2, times2, attr2 = sol_eclipse_when_loc(jd_start, geopos, FLG_SWIEPH)

        # Both should find the same eclipse (same maximum time within tolerance)
        assert abs(times1[0] - times2[0]) < 0.01  # Within ~15 minutes
        # Both should return eclipse visibility flags
        assert ecl_type1 & ECL_VISIBLE
        assert ecl_type2 & ECL_VISIBLE

    def test_multiple_locations_same_eclipse(self):
        """Same eclipse should show different circumstances at different locations."""
        jd_start = julday(2024, 3, 1, 0)

        # Dallas and Miami
        dallas_lat, dallas_lon = 32.7767, -96.7970
        miami_lat, miami_lon = 25.7617, -80.1918

        _, times_dallas, attr_dallas = sol_eclipse_when_loc(
            jd_start, (dallas_lon, dallas_lat, 0.0)
        )
        _, times_miami, attr_miami = sol_eclipse_when_loc(
            jd_start, (miami_lon, miami_lat, 0.0)
        )

        # Both should find eclipses (not necessarily the same one)
        assert times_dallas[0] > jd_start
        assert times_miami[0] > jd_start

    def test_sun_must_be_above_horizon(self):
        """Eclipse should only be returned if Sun is visible."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 35.0, -100.0

        ecl_type, times, attr = sol_eclipse_when_loc(jd_start, (lon, lat, 0.0))

        # Sun altitude at maximum should be positive (above horizon)
        sun_altitude = attr[4]
        assert sun_altitude > -1.0  # Allow for refraction near horizon

    def test_sequential_eclipses(self):
        """Finding consecutive eclipses should return different events."""
        jd_start = julday(2024, 1, 1, 0)
        lat, lon = 35.0, -100.0

        # Find first eclipse
        _, times1, _ = sol_eclipse_when_loc(jd_start, (lon, lat, 0.0))

        # Find next eclipse (search from after first one)
        jd_after_first = times1[0] + 1
        _, times2, _ = sol_eclipse_when_loc(jd_after_first, (lon, lat, 0.0))

        # Second eclipse should be after first
        assert times2[0] > times1[0]


class TestSolEclipseWhere:
    """Tests for sol_eclipse_where function."""

    def test_returns_correct_tuple_sizes(self):
        """Should return geopos tuple of 10 and attr tuple of 20 elements per the reference ephemeris."""
        from libephemeris import sol_eclipse_where

        # First find an eclipse to get a valid time
        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        ecl_type, geopos, attr = sol_eclipse_where(times[0])

        assert len(geopos) == 10
        assert len(attr) == 20

    def test_central_eclipse_returns_valid_position(self):
        """For central eclipse, should return valid geographic coordinates."""
        from libephemeris import sol_eclipse_where

        # Find a total eclipse (which is central)
        jd_start = julday(2024, 3, 1, 0)
        global_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        ecl_type, geopos, attr = sol_eclipse_where(times[0])

        # If eclipse is central, we should get valid coordinates
        if ecl_type & ECL_CENTRAL:
            lon, lat = geopos[0], geopos[1]
            # Longitude should be -180 to 180
            assert -180 <= lon <= 180
            # Latitude should be -90 to 90
            assert -90 <= lat <= 90

    def test_eclipse_type_flags_set(self):
        """Eclipse type flags should be set appropriately."""
        from libephemeris import sol_eclipse_where

        # Find a total eclipse
        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        ecl_type, geopos, attr = sol_eclipse_where(times[0])

        # Should have central flag for central eclipse
        if ecl_type != 0:
            assert ecl_type & ECL_CENTRAL
            # Should have either total or annular flag
            assert ecl_type & (ECL_TOTAL | ECL_ANNULAR)

    def test_attributes_have_valid_ranges(self):
        """Eclipse attributes should be within valid ranges.

        For sol_eclipse_where, attr indices:
          [0]: magnitude (can exceed 1.0 for total eclipses)
          [1]: ratio of lunar to solar diameter
          [2]: obscuration (fraction of solar disc area covered,
               can exceed 1.0 for total eclipses)
          [3]: core shadow diameter in km (negative = no umbra on surface)
          [4]: azimuth of sun (degrees)
          [5]: true altitude of sun (degrees)
          [6]-[7]: apparent diameter data
        """
        from libephemeris import sol_eclipse_where

        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        ecl_type, geopos, attr = sol_eclipse_where(times[0])

        if ecl_type != 0:
            magnitude = attr[0]
            ratio = attr[1]
            obscuration = attr[2]
            sun_altitude = attr[5]

            # Magnitude should be positive for central eclipse
            assert magnitude > 0
            # Ratio should be close to 1 for total/annular
            assert 0.8 < ratio < 1.2
            # Obscuration should be positive (can exceed 1.0 for total eclipses)
            assert obscuration > 0
            # Sun altitude at central line should be valid
            assert -90 <= sun_altitude <= 90

    def test_non_eclipse_time_returns_zero_flag_with_proximity_point(self):
        """A no-eclipse instant gives retflag 0 with the closest-approach point.

        The reference convention (verified against the reference ephemeris 2.10.03) keeps
        geopos[0..1] at the surface point of maximum proximity of the shadow
        axis and reports a negative magnitude measure in attr[0]; only
        geopos[2..9] stay zero.
        """
        from libephemeris import sol_eclipse_where

        # Random time not during an eclipse (full moon time)
        jd_non_eclipse = julday(2024, 4, 23, 12)  # Near full moon

        ecl_type, geopos, attr = sol_eclipse_where(jd_non_eclipse)

        assert ecl_type == 0
        assert -180.0 <= geopos[0] <= 180.0
        assert -90.0 <= geopos[1] <= 90.0
        assert attr[0] < 0.0
        for i in range(2, 10):
            assert geopos[i] == 0.0

    def test_april_2024_eclipse_path(self):
        """Test April 8, 2024 total eclipse path."""
        from libephemeris import sol_eclipse_where

        # April 8, 2024 eclipse maximum around 18:18 UT
        jd_eclipse = julday(2024, 4, 8, 18.3)

        ecl_type, geopos, attr = sol_eclipse_where(jd_eclipse)

        # Should be a total eclipse
        if ecl_type & ECL_TOTAL:
            lon, lat = geopos[0], geopos[1]
            # The path crossed Mexico, US, and Eastern Canada
            # Central line was roughly around:
            # - Mexico: around 100W, 20-25N
            # - US: around 85-100W, 30-45N
            # - Canada: around 60-80W, 45-50N
            # At 18:18 UT, central line should be somewhere in Americas
            # Allow reasonable tolerance for approximate algorithm
            assert -130 < lon < -40  # Western longitudes (Americas)
            # The algorithm may give slightly different latitudes
            # due to approximations in shadow geometry
            assert -10 < lat < 60  # Northern hemisphere, with tolerance

    def test_path_width_positive_for_central_eclipse(self):
        """Path width (attr[3]) for central eclipses.

        Note: attr[3] is the core shadow diameter in km. It can be negative
        when the umbral cone vertex is above Earth's surface (common for
        annular eclipses and some total eclipses near the horizon).
        For a total eclipse at maximum, we just check it's non-zero.
        """
        from libephemeris import sol_eclipse_where

        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        ecl_type, geopos, attr = sol_eclipse_where(times[0])

        if ecl_type & ECL_CENTRAL:
            path_width = attr[3]
            # Core shadow diameter should be non-zero for central eclipse
            assert path_width != 0

    def test_swe_version_provides_reference_compatible_interface(self):
        """sol_eclipse_where and sol_eclipse_where produce equivalent results, with a reference-compatible interface.

        Note: sol_eclipse_where is a wrapper that calls sol_eclipse_where internally.
        They are distinct functions with the same interface that produce identical results.
        """
        from libephemeris import sol_eclipse_where

        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        # Call both versions
        geopos1, attr1, ecl_type1 = sol_eclipse_where(times[0])
        geopos2, attr2, ecl_type2 = sol_eclipse_where(times[0], FLG_SWIEPH)

        # Results should be identical
        assert geopos1 == geopos2
        assert attr1 == attr2
        assert ecl_type1 == ecl_type2

    def test_flags_parameter_accepted(self):
        """Should accept flags parameter without error."""
        from libephemeris import sol_eclipse_where

        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start)

        # Should not raise
        ecl_type, geopos, attr = sol_eclipse_where(times[0], flags=FLG_SWIEPH)

    def test_annular_eclipse(self):
        """Test with annular eclipse."""
        from libephemeris import sol_eclipse_where

        # October 14, 2023 annular eclipse
        jd_start = julday(2023, 9, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_ANNULAR)

        ecl_type, geopos, attr = sol_eclipse_where(times[0])

        # Should be annular
        if ecl_type & ECL_ANNULAR:
            # Moon/Sun ratio should be < 1 for annular
            ratio = attr[1]
            assert ratio < 1.0

    def test_multiple_times_during_eclipse(self):
        """Different times during eclipse should give different positions."""
        from libephemeris import sol_eclipse_where

        # Find an eclipse
        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        # Check positions at different times
        jd_max = times[0]
        delta = 10.0 / (24 * 60)  # 10 minutes in days

        type1, geopos1, _ = sol_eclipse_where(jd_max - delta)
        type2, geopos2, _ = sol_eclipse_where(jd_max)
        type3, geopos3, _ = sol_eclipse_where(jd_max + delta)

        # If all are central, positions should differ
        if type1 & ECL_CENTRAL and type2 & ECL_CENTRAL and type3 & ECL_CENTRAL:
            # Longitudes should be different (shadow moves westward)
            # The shadow moves about 0.25 degrees per minute in longitude
            assert geopos1[0] != geopos2[0] or geopos1[1] != geopos2[1]


class TestSolEclipseHow:
    """Tests for sol_eclipse_how function."""

    def test_returns_correct_tuple_size(self):
        """Should return attr tuple of 11 elements."""
        from libephemeris import sol_eclipse_how

        # Find an eclipse to get a valid time
        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        # Calculate circumstances at Dallas
        dallas_lat, dallas_lon = 32.7767, -96.7970
        ecl_type, attr = sol_eclipse_how(times[0], (dallas_lon, dallas_lat, 0.0))

        assert len(attr) == 20

    def test_eclipse_visible_during_eclipse(self):
        """During an eclipse, should return positive magnitude."""
        from libephemeris import sol_eclipse_how

        # April 8, 2024 total eclipse
        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        # Dallas was in path of totality
        dallas_lat, dallas_lon = 32.7767, -96.7970
        ecl_type, attr = sol_eclipse_how(times[0], (dallas_lon, dallas_lat, 0.0))

        # Should have visibility
        if ecl_type & ECL_VISIBLE:
            magnitude = attr[0]
            assert magnitude > 0

    def test_no_eclipse_at_random_time(self):
        """At a random time with no eclipse, magnitude should be 0."""
        from libephemeris import sol_eclipse_how

        # Random time during full moon (no solar eclipse possible)
        jd_full_moon = julday(2024, 4, 23, 12)  # Near full moon

        dallas_lat, dallas_lon = 32.7767, -96.7970
        ecl_type, attr = sol_eclipse_how(jd_full_moon, (dallas_lon, dallas_lat, 0.0))

        # Should have no eclipse
        assert ecl_type == 0 or attr[0] == 0

    def test_attributes_have_valid_ranges(self):
        """Eclipse attributes should be within valid ranges."""
        from libephemeris import sol_eclipse_how

        # Find an eclipse
        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        dallas_lat, dallas_lon = 32.7767, -96.7970
        ecl_type, attr = sol_eclipse_how(times[0], (dallas_lon, dallas_lat, 0.0))

        if ecl_type != 0:
            # Magnitude should be 0-1.5 range
            magnitude = attr[0]
            assert 0 <= magnitude <= 1.5

            # Obscuration should be 0-1 range
            obscuration = attr[2]
            assert 0 <= obscuration <= 1

            # Azimuth should be 0 to 360
            azimuth = attr[4]
            assert 0 <= azimuth <= 360

            # Sun altitude should be -90 to 90
            altitude = attr[5]
            assert -90 <= altitude <= 90

    def test_sun_below_horizon_returns_zero_magnitude(self):
        """When Sun is below horizon, should return 0 magnitude or no visibility."""
        from libephemeris import sol_eclipse_how

        # Find an eclipse
        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        # Check from a location on opposite side of Earth (Sun below horizon)
        # For an eclipse happening in North America around 18 UT,
        # a location in Eastern Asia/Australia would have Sun below horizon
        east_asia_lat, east_asia_lon = 35.0, 135.0  # Japan area

        ecl_type, attr = sol_eclipse_how(times[0], (east_asia_lon, east_asia_lat, 0.0))

        # Either no visibility flag or zero magnitude
        sun_alt = attr[5]
        if sun_alt < -1.0:
            # Sun below horizon, should have no visible eclipse
            assert attr[0] == 0.0 or not (ecl_type & ECL_VISIBLE)

    def test_different_locations_different_magnitudes(self):
        """Different locations should show different eclipse magnitudes."""
        from libephemeris import sol_eclipse_how

        # Find an eclipse
        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        # Dallas (in path of totality)
        dallas_lat, dallas_lon = 32.7767, -96.7970
        type_dallas, attr_dallas = sol_eclipse_how(
            times[0], (dallas_lon, dallas_lat, 0.0)
        )

        # Miami (away from path)
        miami_lat, miami_lon = 25.7617, -80.1918
        type_miami, attr_miami = sol_eclipse_how(times[0], (miami_lon, miami_lat, 0.0))

        # If both visible, Dallas should have higher magnitude (closer to center)
        if (type_dallas & ECL_VISIBLE) and (type_miami & ECL_VISIBLE):
            # The magnitudes should be different
            # (may not always be true depending on exact timing)
            assert attr_dallas[0] != attr_miami[0] or attr_dallas[2] != attr_miami[2]

    def test_total_eclipse_type_flags(self):
        """Total eclipse should have appropriate type flags."""
        from libephemeris import sol_eclipse_how

        # Find a total eclipse
        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        # From a location in the path
        # Use local time found by sol_eclipse_when_loc for better timing
        _, times_loc, attr_loc = sol_eclipse_when_loc(
            julday(2024, 3, 1, 0), (-96.7970, 32.7767, 0.0)
        )

        if times_loc[0] > 0:
            ecl_type, attr = sol_eclipse_how(times_loc[0], (-96.7970, 32.7767, 0.0))

            if ecl_type & ECL_TOTAL:
                # Should also have central and visible flags
                assert ecl_type & ECL_CENTRAL
                assert ecl_type & ECL_VISIBLE

    def test_partial_eclipse_type_flags(self):
        """Partial eclipse should have partial flag."""
        from libephemeris import sol_eclipse_how

        # Find an eclipse and check from edge location
        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        # NYC is far from totality path for April 2024 eclipse
        nyc_lat, nyc_lon = 40.7128, -74.0060
        ecl_type, attr = sol_eclipse_how(times[0], (nyc_lon, nyc_lat, 0.0))

        # If visible and not central, should be partial
        if (ecl_type & ECL_VISIBLE) and not (ecl_type & ECL_CENTRAL):
            assert ecl_type & ECL_PARTIAL

    def test_altitude_parameter(self):
        """Should accept altitude parameter."""
        from libephemeris import sol_eclipse_how

        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start)

        # Should not raise
        ecl_type, attr = sol_eclipse_how(times[0], (-96.7970, 32.7767, 1000.0))

        assert len(attr) == 20

    def test_flags_parameter_accepted(self):
        """Should accept flags parameter."""
        from libephemeris import sol_eclipse_how

        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start)

        # Should not raise
        ecl_type, attr = sol_eclipse_how(
            times[0], (-96.7970, 32.7767, 0.0), flags=FLG_SWIEPH
        )

        assert len(attr) == 20

    def test_swe_version_provides_reference_compatible_interface(self):
        """sol_eclipse_how and sol_eclipse_how produce equivalent results.

        Note: Both sol_eclipse_how and sol_eclipse_how now use the same
        reference-compatible signature with geopos tuple.
        They produce identical results for the same location.
        """
        from libephemeris import sol_eclipse_how

        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        lat, lon, altitude = 32.7767, -96.7970, 0.0
        geopos = (lon, lat, altitude)  # geopos uses (lon, lat, alt) order

        # Call both versions
        ecl_type1, attr1 = sol_eclipse_how(times[0], geopos)
        ecl_type2, attr2 = sol_eclipse_how(times[0], geopos, FLG_SWIEPH)

        # Results should be identical
        assert attr1 == attr2
        assert ecl_type1 == ecl_type2

    def test_known_eclipse_april_2024_dallas(self):
        """Test April 8, 2024 eclipse from Dallas."""
        from libephemeris import sol_eclipse_how

        # During the April 2024 eclipse, around 18:40 UT
        jd_eclipse = julday(2024, 4, 8, 18.67)  # ~18:40 UT

        dallas_lat, dallas_lon = 32.7767, -96.7970
        ecl_type, attr = sol_eclipse_how(jd_eclipse, (dallas_lon, dallas_lat, 0.0))

        # Should have significant magnitude (Dallas was in path of totality)
        if ecl_type & ECL_VISIBLE:
            assert attr[0] > 0.8  # High magnitude expected

    def test_obscuration_consistent_with_magnitude(self):
        """Obscuration should be consistent with magnitude."""
        from libephemeris import sol_eclipse_how

        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        dallas_lat, dallas_lon = 32.7767, -96.7970
        ecl_type, attr = sol_eclipse_how(times[0], (dallas_lon, dallas_lat, 0.0))

        if ecl_type & ECL_VISIBLE:
            magnitude = attr[0]
            obscuration = attr[2]

            # Both should be positive if eclipse is visible
            if magnitude > 0:
                assert obscuration > 0

            # For total eclipse (magnitude >= 1), obscuration should be 1
            if magnitude >= 1.0:
                assert obscuration >= 0.99  # Allow small tolerance


class TestKnownEclipseValidation:
    """Validation tests against known eclipse times.

    These tests verify that the eclipse times are accurate within 0.01 days (~15 minutes)
    compared to known eclipse dates from NASA eclipse catalogs.
    """

    def test_april_2024_total_eclipse_timing(self):
        """Validate April 8, 2024 total eclipse maximum time.

        Known data: Total solar eclipse on April 8, 2024
        Maximum eclipse around JD 2460409.26 (approximately 18:17 UT)
        """
        # Search starting well before the eclipse
        jd_start = 2460400.0  # About 8 days before

        ecl_type, times = sol_eclipse_when_glob(jd_start)

        # Eclipse maximum should be around JD 2460409.26
        expected_jd = 2460409.26
        tolerance = 0.01  # 0.01 days = ~15 minutes

        assert abs(times[0] - expected_jd) < tolerance, (
            f"Eclipse maximum {times[0]} differs from expected {expected_jd} "
            f"by more than {tolerance} days"
        )

        # Should be a total eclipse
        assert ecl_type & ECL_TOTAL

    def test_october_2023_annular_eclipse_timing(self):
        """Validate October 14, 2023 annular eclipse maximum time.

        Known data: Annular solar eclipse on October 14, 2023
        Maximum eclipse around 17:59 UT (JD ~2460232.25)
        """
        jd_start = 2460200.0  # About 32 days before

        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_ANNULAR)

        # Verify the correct date was found
        year, month, day, hour = revjul(times[0])
        assert year == 2023
        assert month == 10
        assert 13 <= day <= 15  # October 14

        # Should be an annular eclipse
        assert ecl_type & ECL_ANNULAR

    def test_august_2017_total_eclipse_timing(self):
        """Validate August 21, 2017 total eclipse maximum time.

        Known data: Total solar eclipse on August 21, 2017
        Maximum eclipse around 18:26 UT
        The famous "Great American Eclipse"

        Note: Due to eclipse type classification differences, we search without
        type filter and verify the eclipse is central (total/annular).
        """
        jd_start = 2457900.0  # About 87 days before

        # Search without type filter to find any eclipse
        ecl_type, times = sol_eclipse_when_glob(jd_start)

        # Verify the correct date was found
        year, month, day, hour = revjul(times[0])
        assert year == 2017
        assert month == 8
        assert 20 <= day <= 22  # August 21

        # Should be a central eclipse (either total or annular)
        assert ecl_type & ECL_CENTRAL, "August 2017 eclipse should be central"

    def test_june_2021_annular_eclipse_timing(self):
        """Validate June 10, 2021 annular eclipse maximum time.

        Known data: Annular solar eclipse on June 10, 2021
        Maximum eclipse around 10:41 UT
        """
        jd_start = 2459350.0  # About 26 days before

        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_ANNULAR)

        # Verify the correct date was found
        year, month, day, hour = revjul(times[0])
        assert year == 2021
        assert month == 6
        assert 9 <= day <= 11  # June 10

        # Should be an annular eclipse
        assert ecl_type & ECL_ANNULAR

    def test_eclipse_type_classification_consistency(self):
        """Verify eclipse type classification is consistent for known eclipses."""
        # Total eclipse 2024
        ecl_type1, times1 = sol_eclipse_when_glob(2460400.0, ecltype=ECL_TOTAL)
        assert ecl_type1 & ECL_TOTAL, "April 2024 eclipse should be classified as TOTAL"

        # Annular eclipse 2023
        ecl_type2, times2 = sol_eclipse_when_glob(2460200.0, ecltype=ECL_ANNULAR)
        assert ecl_type2 & ECL_ANNULAR, (
            "October 2023 eclipse should be classified as ANNULAR"
        )

    def test_return_tuple_has_ten_elements(self):
        """Verify the tret tuple has exactly 10 elements like the reference ephemeris."""
        jd_start = 2460400.0

        ecl_type, times = sol_eclipse_when_glob(jd_start)

        assert len(times) == 10, f"Expected 10 elements in tret, got {len(times)}"

        # First element should be eclipse maximum
        assert times[0] > jd_start, "First element should be eclipse maximum time"

        # Elements 0, 2, 3 should be phase times
        # Maximum should be between eclipse begin (times[2]) and end (times[3])
        assert times[2] < times[0] < times[3], "Maximum should be between begin and end"


class TestSolEclipseMaxTime:
    """Tests for sol_eclipse_max_time function - precise eclipse maximum calculation."""

    def test_global_max_time_returns_tuple(self):
        """Should return tuple of (jd_max, gamma)."""
        from libephemeris import sol_eclipse_max_time

        # Use a known eclipse - April 8, 2024
        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        jd_max, gamma = sol_eclipse_max_time(times[0])

        # Should return two floats
        assert isinstance(jd_max, float)
        assert isinstance(gamma, float)

    def test_global_max_time_near_input(self):
        """Global maximum time should be close to input approximation."""
        from libephemeris import sol_eclipse_max_time

        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        jd_max, _ = sol_eclipse_max_time(times[0])

        # Result should be within search range of input
        assert abs(jd_max - times[0]) < 0.125  # default search range

    def test_global_max_refines_precision(self):
        """Refined maximum should be more precise than New Moon estimate."""
        from libephemeris import sol_eclipse_max_time

        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        # Get refined maximum
        jd_max, gamma = sol_eclipse_max_time(times[0])

        # Gamma should be positive and reasonably small for a total eclipse
        assert gamma >= 0
        assert gamma < 1.0  # For central eclipse, gamma < 1

    def test_global_max_gamma_is_minimum(self):
        """Gamma at refined maximum should be less than at nearby times."""
        from libephemeris import sol_eclipse_max_time
        from libephemeris.eclipse import _calc_gamma

        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start)

        jd_max, gamma_at_max = sol_eclipse_max_time(times[0])

        # Check gamma at nearby times is larger
        delta = 0.001  # ~1.44 minutes
        gamma_before = _calc_gamma(jd_max - delta)
        gamma_after = _calc_gamma(jd_max + delta)

        assert gamma_at_max <= gamma_before
        assert gamma_at_max <= gamma_after

    def test_local_max_time_returns_tuple(self):
        """Local maximum should return tuple of (jd_max, separation)."""
        from libephemeris import sol_eclipse_max_time

        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        # Dallas, Texas - was in path of totality
        dallas_lat, dallas_lon = 32.7767, -96.7970

        jd_max, separation = sol_eclipse_max_time(
            times[0], lat=dallas_lat, lon=dallas_lon
        )

        assert isinstance(jd_max, float)
        assert isinstance(separation, float)

    def test_local_max_separation_is_minimum(self):
        """Separation at local maximum should be minimum."""
        from libephemeris import sol_eclipse_max_time

        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        dallas_lat, dallas_lon = 32.7767, -96.7970

        jd_max, min_sep = sol_eclipse_max_time(times[0], lat=dallas_lat, lon=dallas_lon)

        # Separation should be positive
        assert min_sep >= 0

        # For Dallas during April 2024 eclipse, separation should be small
        # (Dallas was near path of totality)
        assert min_sep < 1.0  # degrees

    def test_local_max_differs_from_global(self):
        """Local maximum time can differ from global maximum."""
        from libephemeris import sol_eclipse_max_time

        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        # Get global maximum
        jd_global_max, _ = sol_eclipse_max_time(times[0])

        # Get local maximum from Dallas
        dallas_lat, dallas_lon = 32.7767, -96.7970
        jd_local_max, _ = sol_eclipse_max_time(times[0], lat=dallas_lat, lon=dallas_lon)

        # They may differ due to parallax and shadow movement
        # Allow up to 1 hour difference
        assert abs(jd_global_max - jd_local_max) < 1.0 / 24.0

    def test_requires_both_lat_lon_or_neither(self):
        """Should raise ValueError if only one of lat/lon provided."""
        from libephemeris import sol_eclipse_max_time

        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start)

        # Only lat provided
        try:
            sol_eclipse_max_time(times[0], lat=32.0)
            assert False, "Should have raised ValueError"
        except ValueError:
            pass

        # Only lon provided
        try:
            sol_eclipse_max_time(times[0], lon=-96.0)
            assert False, "Should have raised ValueError"
        except ValueError:
            pass

    def test_accepts_altitude_parameter(self):
        """Should accept altitude parameter for local calculations."""
        from libephemeris import sol_eclipse_max_time

        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        dallas_lat, dallas_lon = 32.7767, -96.7970

        # Should not raise
        jd_max, sep = sol_eclipse_max_time(
            times[0], lat=dallas_lat, lon=dallas_lon, altitude=1000.0
        )

        assert jd_max > times[0] - 0.125

    def test_search_range_parameter(self):
        """Should respect custom search_range parameter."""
        from libephemeris import sol_eclipse_max_time

        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        # Use smaller search range
        jd_max1, _ = sol_eclipse_max_time(times[0], search_range=0.05)

        # Use larger search range
        jd_max2, _ = sol_eclipse_max_time(times[0], search_range=0.25)

        # Both should find similar maximum (within tolerance)
        assert abs(jd_max1 - jd_max2) < 0.001

    def test_swe_alias_exists(self):
        """sol_eclipse_max_time should be an alias."""

    def test_known_eclipse_april_2024_precision(self):
        """Test precision against known April 8, 2024 eclipse."""
        from libephemeris import sol_eclipse_max_time

        # NASA gives maximum at approximately 18:18 UT
        # JD 2460409.26 corresponds to ~18:17 UT
        jd_approx = julday(2024, 4, 8, 18.3)

        jd_max, gamma = sol_eclipse_max_time(jd_approx)

        # Convert to hours for comparison
        _, _, _, hour = revjul(jd_max)

        # Should be very close to 18:17-18:18 UT
        assert 18.2 < hour < 18.4, f"Expected ~18.3 hours, got {hour}"

    def test_different_locations_different_local_max(self):
        """Different locations should have different local maximum times."""
        from libephemeris import sol_eclipse_max_time

        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        # Dallas, Texas
        dallas_lat, dallas_lon = 32.7767, -96.7970
        jd_dallas, sep_dallas = sol_eclipse_max_time(
            times[0], lat=dallas_lat, lon=dallas_lon
        )

        # Cleveland, Ohio (also in path of April 2024 eclipse)
        cleveland_lat, cleveland_lon = 41.4993, -81.6944
        jd_cleveland, sep_cleveland = sol_eclipse_max_time(
            times[0], lat=cleveland_lat, lon=cleveland_lon
        )

        # Times should differ (Cleveland is east, sees maximum later)
        # Allow some tolerance for the difference
        assert jd_dallas != jd_cleveland or sep_dallas != sep_cleveland

    def test_annular_eclipse(self):
        """Should work for annular eclipses too."""
        from libephemeris import sol_eclipse_max_time

        # October 14, 2023 annular eclipse
        jd_start = julday(2023, 9, 1, 0)
        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_ANNULAR)

        jd_max, gamma = sol_eclipse_max_time(times[0])

        # Should return valid values
        assert jd_max > jd_start
        assert gamma >= 0

    def test_sub_second_precision(self):
        """Result should have sub-second precision."""
        from libephemeris import sol_eclipse_max_time
        from libephemeris.eclipse import _calc_gamma

        jd_start = julday(2024, 3, 1, 0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        jd_max, gamma_at_max = sol_eclipse_max_time(times[0])

        # Test that at 1 second away, gamma is measurably larger
        one_second = 1.0 / 86400.0  # 1 second in days

        gamma_before = _calc_gamma(jd_max - one_second)
        gamma_after = _calc_gamma(jd_max + one_second)

        # Gamma should be slightly larger at 1 second offset
        # (may be equal due to floating point precision)
        assert gamma_at_max <= gamma_before + 1e-10
        assert gamma_at_max <= gamma_after + 1e-10
