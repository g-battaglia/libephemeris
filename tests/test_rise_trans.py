"""
Tests for rise_trans function in libephemeris.

Tests the calculation of rise, set, and transit times for celestial bodies.

Reference data from USNO and timeanddate.com for verification.
"""

import pytest

pytestmark = pytest.mark.slow

from libephemeris import (
    julday,
    revjul,
    rise_trans,
    SUN,
    MOON,
    MARS,
    JUPITER,
    CALC_RISE,
    CALC_SET,
    CALC_MTRANSIT,
    CALC_ITRANSIT,
    BIT_DISC_CENTER,
    BIT_NO_REFRACTION,
    BIT_CIVIL_TWILIGHT,
    BIT_NAUTIC_TWILIGHT,
    BIT_ASTRO_TWILIGHT,
    BIT_GEOCTR_NO_ECL_LAT,
    BIT_HINDU_RISING,
)


class TestRiseTransBasic:
    """Basic tests for rise_trans function."""

    def test_sunrise_returns_valid_time(self):
        """Test that sunrise calculation returns a valid time."""
        # June 21, 2024 - summer solstice
        jd_start = julday(2024, 6, 21, 0)
        # Rome, Italy
        lat, lon = 41.9028, 12.4964

        retflag, tret = rise_trans(jd_start, SUN, CALC_RISE, [lon, lat, 0])
        jd_rise = tret[0]

        # Should find sunrise on June 21
        assert jd_rise > jd_start
        assert jd_rise < jd_start + 1  # Within 24 hours
        assert retflag == 0

        # Verify it's in the morning (before noon local time)
        year, month, day, hour = revjul(jd_rise)
        assert year == 2024
        assert month == 6
        assert day == 21
        # Rome is UTC+1 (summer: UTC+2), sunrise should be ~5:30 local = ~3:30 UTC
        assert 3.0 <= hour <= 6.0

    def test_sunset_returns_valid_time(self):
        """Test that sunset calculation returns a valid time."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        retflag, tret = rise_trans(jd_start, SUN, CALC_SET, [lon, lat, 0])
        jd_set = tret[0]

        # Should find sunset on June 21
        assert jd_set > jd_start
        assert jd_set < jd_start + 1
        assert retflag == 0

        # Verify it's in the evening (after noon)
        year, month, day, hour = revjul(jd_set)
        assert year == 2024
        assert month == 6
        # Sunset should be ~20:50 local = ~18:50 UTC
        assert 18.0 <= hour <= 21.0

    def test_transit_returns_valid_time(self):
        """Test that meridian transit (noon) calculation works."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        retflag, tret = rise_trans(jd_start, SUN, CALC_MTRANSIT, [lon, lat, 0])
        jd_transit = tret[0]

        # Should find transit on June 21
        assert jd_transit > jd_start
        assert jd_transit < jd_start + 1
        assert retflag == 0

        # Verify it's around solar noon (12:00-13:00 local)
        year, month, day, hour = revjul(jd_transit)
        assert year == 2024
        # Solar noon at Rome (~12.5°E) is about 12:10 - 0:50 (longitude offset) = ~11:20 UTC
        assert 10.5 <= hour <= 13.0

    def test_lower_transit_returns_valid_time(self):
        """Test that lower transit (midnight) calculation works."""
        jd_start = julday(2024, 6, 21, 12)  # Start at noon
        lat, lon = 41.9028, 12.4964  # Rome

        retflag, tret = rise_trans(jd_start, SUN, CALC_ITRANSIT, [lon, lat, 0])
        jd_itransit = tret[0]

        # Should find lower transit around midnight
        assert jd_itransit > jd_start
        assert jd_itransit < jd_start + 1
        assert retflag == 0

    def test_swe_alias_works(self):
        """Test that rise_trans is an alias for rise_trans."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964

        result1 = rise_trans(jd_start, SUN, CALC_RISE, [lon, lat, 0])
        result2 = rise_trans(jd_start, SUN, CALC_RISE, [lon, lat, 0])

        assert result1 == result2


class TestRiseTransMoon:
    """Tests for Moon rise/set calculations."""

    def test_moonrise_valid(self):
        """Test moonrise calculation."""
        jd_start = julday(2024, 1, 15, 0)
        lat, lon = 51.5074, -0.1278  # London

        retflag, tret = rise_trans(jd_start, MOON, CALC_RISE, [lon, lat, 0])
        jd_rise = tret[0]

        assert jd_rise > jd_start
        assert jd_rise < jd_start + 2  # Moon can rise any time
        assert retflag == 0

    def test_moonset_valid(self):
        """Test moonset calculation."""
        jd_start = julday(2024, 1, 15, 0)
        lat, lon = 51.5074, -0.1278  # London

        retflag, tret = rise_trans(jd_start, MOON, CALC_SET, [lon, lat, 0])
        jd_set = tret[0]

        assert jd_set > jd_start
        assert jd_set < jd_start + 2
        assert retflag == 0

    def test_moon_transit_valid(self):
        """Test Moon transit calculation."""
        jd_start = julday(2024, 1, 15, 0)
        lat, lon = 51.5074, -0.1278  # London

        retflag, tret = rise_trans(jd_start, MOON, CALC_MTRANSIT, [lon, lat, 0])
        jd_transit = tret[0]

        assert jd_transit > jd_start
        assert jd_transit < jd_start + 1.5  # Moon transit within ~1 day
        assert retflag == 0


class TestRiseTransPlanets:
    """Tests for planet rise/set calculations."""

    def test_mars_rise(self):
        """Test Mars rise calculation."""
        jd_start = julday(2024, 1, 15, 0)
        lat, lon = 40.7128, -74.0060  # New York

        retflag, tret = rise_trans(jd_start, MARS, CALC_RISE, [lon, lat, 0])
        jd_rise = tret[0]

        assert jd_rise > jd_start
        assert jd_rise < jd_start + 2
        assert retflag == 0

    def test_jupiter_transit(self):
        """Test Jupiter transit calculation."""
        jd_start = julday(2024, 1, 15, 0)
        lat, lon = 40.7128, -74.0060  # New York

        retflag, tret = rise_trans(jd_start, JUPITER, CALC_MTRANSIT, [lon, lat, 0])
        jd_transit = tret[0]

        assert jd_transit > jd_start
        assert jd_transit < jd_start + 1.5
        assert retflag == 0


class TestRiseTransCircumpolar:
    """Tests for circumpolar objects."""

    def test_sun_circumpolar_arctic_summer(self):
        """Test that Sun doesn't set in Arctic summer (midnight sun)."""
        # Svalbard in late June - Sun should be circumpolar
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 78.0, 16.0  # Svalbard

        retflag, tret = rise_trans(jd_start, SUN, CALC_SET, [lon, lat, 0])

        # Should return -2 for circumpolar (never sets)
        assert retflag == -2
        assert tret[0] == 0.0

    def test_sun_circumpolar_arctic_winter(self):
        """Test that Sun doesn't rise in Arctic winter (polar night)."""
        # Svalbard in late December - Sun should never rise
        jd_start = julday(2024, 12, 21, 0)
        lat, lon = 78.0, 16.0  # Svalbard

        retflag, tret = rise_trans(jd_start, SUN, CALC_RISE, [lon, lat, 0])

        # Should return -2 for circumpolar (never rises)
        assert retflag == -2
        assert tret[0] == 0.0

    def test_transit_still_works_for_circumpolar(self):
        """Test that transit calculations work for circumpolar objects."""
        # Even circumpolar Sun crosses meridian
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 78.0, 16.0  # Svalbard

        retflag, tret = rise_trans(jd_start, SUN, CALC_MTRANSIT, [lon, lat, 0])
        jd_transit = tret[0]

        # Transit should still be calculable
        assert jd_transit > jd_start
        assert jd_transit < jd_start + 1
        assert retflag == 0


class TestRiseTransFlags:
    """Tests for various flag combinations."""

    def test_disc_center_flag(self):
        """Test BIT_DISC_CENTER flag."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        # With upper limb (default)
        _, tret = rise_trans(jd_start, SUN, CALC_RISE, [lon, lat, 0])
        jd_rise_limb = tret[0]

        # With disc center
        _, tret = rise_trans(jd_start, SUN, CALC_RISE | BIT_DISC_CENTER, [lon, lat, 0])
        jd_rise_center = tret[0]

        # Center rise should be later than upper limb rise
        assert jd_rise_center > jd_rise_limb
        # Difference should be about 1-2 minutes (Sun's semi-diameter crossing)
        diff_minutes = (jd_rise_center - jd_rise_limb) * 24 * 60
        assert 0.5 < diff_minutes < 3.0

    def test_no_refraction_flag(self):
        """Test BIT_NO_REFRACTION flag."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        # With refraction (default)
        _, tret = rise_trans(jd_start, SUN, CALC_RISE, [lon, lat, 0])
        jd_rise_refr = tret[0]

        # Without refraction
        _, tret = rise_trans(
            jd_start, SUN, CALC_RISE | BIT_NO_REFRACTION, [lon, lat, 0]
        )
        jd_rise_no_refr = tret[0]

        # Without refraction, rise should be later (Sun appears lower)
        assert jd_rise_no_refr > jd_rise_refr
        # Difference should be about 2-4 minutes (34' refraction)
        diff_minutes = (jd_rise_no_refr - jd_rise_refr) * 24 * 60
        assert 1.0 < diff_minutes < 6.0

    def test_civil_twilight(self):
        """Test civil twilight calculation."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        # Civil twilight begins (Sun at -6 degrees)
        _, tret = rise_trans(
            jd_start, SUN, CALC_RISE | BIT_CIVIL_TWILIGHT, [lon, lat, 0]
        )
        jd_twilight = tret[0]

        # Regular sunrise
        _, tret = rise_trans(jd_start, SUN, CALC_RISE, [lon, lat, 0])
        jd_rise = tret[0]

        # Civil twilight should be before sunrise
        assert jd_twilight < jd_rise
        # Difference should be about 25-35 minutes
        diff_minutes = (jd_rise - jd_twilight) * 24 * 60
        assert 15 < diff_minutes < 60

    def test_nautical_twilight(self):
        """Test nautical twilight calculation."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        # Nautical twilight begins (Sun at -12 degrees)
        _, tret = rise_trans(
            jd_start, SUN, CALC_RISE | BIT_NAUTIC_TWILIGHT, [lon, lat, 0]
        )
        jd_nautical = tret[0]

        # Civil twilight
        _, tret = rise_trans(
            jd_start, SUN, CALC_RISE | BIT_CIVIL_TWILIGHT, [lon, lat, 0]
        )
        jd_civil = tret[0]

        # Nautical twilight should be before civil twilight
        assert jd_nautical < jd_civil

    def test_astronomical_twilight(self):
        """Test astronomical twilight calculation."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        # Astronomical twilight begins (Sun at -18 degrees)
        _, tret = rise_trans(
            jd_start, SUN, CALC_RISE | BIT_ASTRO_TWILIGHT, [lon, lat, 0]
        )
        jd_astro = tret[0]

        # Nautical twilight
        _, tret = rise_trans(
            jd_start, SUN, CALC_RISE | BIT_NAUTIC_TWILIGHT, [lon, lat, 0]
        )
        jd_nautical = tret[0]

        # Astronomical twilight should be before nautical twilight
        assert jd_astro < jd_nautical


class TestRiseTransLocations:
    """Tests for different geographic locations."""

    def test_equator(self):
        """Test rise/set at equator."""
        jd_start = julday(2024, 3, 20, 0)  # Near equinox
        lat, lon = 0.0, 0.0  # Equator, prime meridian

        _, tret = rise_trans(jd_start, SUN, CALC_RISE, [lon, lat, 0])
        jd_rise = tret[0]
        _, tret = rise_trans(jd_start, SUN, CALC_SET, [lon, lat, 0])
        jd_set = tret[0]

        # At equator during equinox, day and night are nearly equal
        day_length = (jd_set - jd_rise) * 24
        assert 11.5 < day_length < 12.5

    def test_southern_hemisphere(self):
        """Test rise/set in southern hemisphere."""
        jd_start = julday(2024, 12, 21, 0)  # Southern summer
        lat, lon = -33.8688, 151.2093  # Sydney

        _, tret = rise_trans(jd_start, SUN, CALC_RISE, [lon, lat, 0])
        jd_rise = tret[0]
        _, tret = rise_trans(jd_start, SUN, CALC_SET, [lon, lat, 0])
        jd_set = tret[0]

        # Should have valid times
        assert jd_rise > jd_start
        assert jd_set > jd_start

        # Find next sunrise/sunset that form a coherent pair
        # (starting from midnight UTC, for Sydney the set might be before or after rise)
        if jd_set < jd_rise:
            # Sunset happened before sunrise, get next sunset after sunrise
            _, tret = rise_trans(jd_rise + 0.01, SUN, CALC_SET, [lon, lat, 0])
            jd_set2 = tret[0]
            day_length = (jd_set2 - jd_rise) * 24
        else:
            day_length = (jd_set - jd_rise) * 24

        # Long summer day in Sydney
        assert day_length > 13  # Sydney has ~14h days in December

    def test_high_latitude_not_circumpolar(self):
        """Test high latitude location that's not circumpolar."""
        jd_start = julday(2024, 3, 20, 0)  # Equinox
        lat, lon = 65.0, 25.0  # Northern Finland

        flag_rise, tret = rise_trans(jd_start, SUN, CALC_RISE, [lon, lat, 0])
        jd_rise = tret[0]
        flag_set, tret = rise_trans(jd_start, SUN, CALC_SET, [lon, lat, 0])
        jd_set = tret[0]

        # At equinox, even high latitudes have sunrise/sunset
        assert flag_rise != -2
        assert flag_set != -2
        assert jd_rise > jd_start
        assert jd_set > jd_rise


class TestRiseTransErrors:
    """Tests for error handling."""

    def test_invalid_planet_raises_error(self):
        """An unknown body id raises the typed error on every backend."""
        from libephemeris.exceptions import Error

        jd_start = julday(2024, 6, 21, 0)

        with pytest.raises(Error):
            rise_trans(jd_start, 9999, CALC_RISE, [12.5, 41.9, 0])

    def test_bare_rsmi_defaults_to_rise(self):
        """rsmi without an event bit defaults to a rise search.

        Measured reference behavior (2.10.3.2): rsmi=0 and rsmi=16 (a bare
        modifier bit) both return the same event and instant as rsmi=1.
        """
        jd_start = julday(2024, 6, 21, 0)
        base = rise_trans(jd_start, SUN, 1, [12.5, 41.9, 0])
        for rsmi in (0, 16):
            got = rise_trans(jd_start, SUN, rsmi, [12.5, 41.9, 0])
            assert got[0] == base[0]
            assert abs(got[1][0] - base[1][0]) < 1e-9


class TestRiseTransSequential:
    """Tests for sequential event calculations."""

    def test_multiple_sunrises(self):
        """Test finding multiple consecutive sunrises."""
        jd = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964  # Rome
        sunrises = []

        for _ in range(3):
            _, tret = rise_trans(jd, SUN, CALC_RISE, [lon, lat, 0])
            jd_rise = tret[0]
            sunrises.append(jd_rise)
            jd = jd_rise + 0.5  # Start from noon after sunrise

        # Each sunrise should be about 1 day apart
        for i in range(1, len(sunrises)):
            diff = sunrises[i] - sunrises[i - 1]
            assert 0.99 < diff < 1.01  # About 1 day

    def test_sunrise_before_sunset_same_day(self):
        """Test that sunrise comes before sunset on the same day."""
        jd_start = julday(2024, 6, 21, 0)  # Midnight
        lat, lon = 41.9028, 12.4964  # Rome

        _, tret = rise_trans(jd_start, SUN, CALC_RISE, [lon, lat, 0])
        jd_rise = tret[0]
        _, tret = rise_trans(jd_start, SUN, CALC_SET, [lon, lat, 0])
        jd_set = tret[0]

        # Sunrise should be before sunset (starting from midnight)
        assert jd_rise < jd_set

        # Both should be on the same day
        y1, m1, d1, _ = revjul(jd_rise)
        y2, m2, d2, _ = revjul(jd_set)
        assert (y1, m1, d1) == (y2, m2, d2)

    def test_transit_between_rise_and_set(self):
        """Test that transit occurs between rise and set."""
        jd_start = julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        _, tret = rise_trans(jd_start, SUN, CALC_RISE, [lon, lat, 0])
        jd_rise = tret[0]
        _, tret = rise_trans(jd_start, SUN, CALC_MTRANSIT, [lon, lat, 0])
        jd_transit = tret[0]
        _, tret = rise_trans(jd_start, SUN, CALC_SET, [lon, lat, 0])
        jd_set = tret[0]

        # Transit should be between rise and set
        assert jd_rise < jd_transit < jd_set


class TestRiseTransStarTwilight:
    """Twilight bits apply to fixed stars (Sun-style depression crossing).

    Measured reference behavior honors CIVIL/NAUTIC/ASTRO for the Sun and
    for fixed stars — the event becomes a geometric center crossing of
    -6/-12/-18 degrees — while planets and the Moon keep the ordinary
    horizon crossing.
    """

    def test_star_twilight_equals_true_horizon_depression(self):
        from libephemeris import rise_trans_true_hor

        jd_start = julday(2019, 12, 21, 0)
        geopos = [36.8219, -1.2921, 1795.0]  # Nairobi
        for bits, depth in (
            (BIT_CIVIL_TWILIGHT, -6.0),
            (BIT_NAUTIC_TWILIGHT, -12.0),
            (BIT_ASTRO_TWILIGHT, -18.0),
        ):
            rf_tw, tret_tw = rise_trans(jd_start, "Aldebaran", CALC_RISE | bits, geopos)
            rf_hor, tret_hor = rise_trans_true_hor(
                jd_start,
                "Aldebaran",
                CALC_RISE | BIT_NO_REFRACTION | BIT_DISC_CENTER,
                geopos,
                horhgt=depth,
            )
            assert rf_tw == rf_hor == 0
            assert tret_tw[0] == pytest.approx(tret_hor[0], abs=1e-6)

    def test_star_twilight_ordering(self):
        jd_start = julday(2019, 12, 21, 0)
        geopos = [36.8219, -1.2921, 1795.0]
        times = []
        for bits in (BIT_ASTRO_TWILIGHT, BIT_NAUTIC_TWILIGHT, BIT_CIVIL_TWILIGHT):
            _, tret = rise_trans(jd_start, "Aldebaran", CALC_RISE | bits, geopos)
            times.append(tret[0])
        _, tret = rise_trans(
            jd_start, "Aldebaran", CALC_RISE | BIT_NO_REFRACTION, geopos
        )
        # Deeper depression classes rise earlier; the plain horizon last.
        assert times[0] < times[1] < times[2] < tret[0]

    def test_planet_twilight_bits_are_ignored(self):
        jd_start = julday(2019, 12, 21, 0)
        geopos = [36.8219, -1.2921, 1795.0]
        _, tret_tw = rise_trans(jd_start, MARS, CALC_RISE | BIT_CIVIL_TWILIGHT, geopos)
        _, tret_plain = rise_trans(jd_start, MARS, CALC_RISE, geopos)
        assert tret_tw[0] == pytest.approx(tret_plain[0], abs=1e-9)


class TestRiseTransGeoctrNoEclLat:
    """BIT_GEOCTR_NO_ECL_LAT (128): geocentric place, ecliptic latitude zeroed.

    The event uses the body's geocentric apparent place projected onto the
    ecliptic (latitude set to zero), so a body far from the ecliptic - most
    dramatically the Moon - rises/sets tens of minutes away from its ordinary
    time. These are self-contained property/semantic checks (no reference
    oracle): identity of BIT_HINDU_RISING, the minutes-scale Moon shift, the
    transit/twilight scope, and a direct check that the returned event really
    lands where the latitude-zeroed geocentric direction meets the horizon.
    """

    ROME = [12.4964, 41.9028, 0.0]

    def test_hindu_rising_is_no_ecl_lat_disc_center_no_refraction(self):
        """896 is exactly 128 | 256 | 512 and behaves identically."""
        assert BIT_HINDU_RISING == BIT_GEOCTR_NO_ECL_LAT | BIT_DISC_CENTER | 512
        jd_start = julday(2005, 3, 10, 0)
        for body in (SUN, MOON, MARS):
            for event in (CALC_RISE, CALC_SET):
                r1, t1 = rise_trans(jd_start, body, event | BIT_HINDU_RISING, self.ROME)
                r2, t2 = rise_trans(
                    jd_start,
                    body,
                    event | BIT_GEOCTR_NO_ECL_LAT | BIT_DISC_CENTER | BIT_NO_REFRACTION,
                    self.ROME,
                )
                assert r1 == r2 == 0
                assert t1[0] == pytest.approx(t2[0], abs=1e-9)

    def test_moon_rise_shifts_minutes_with_no_ecl_lat(self):
        """The Moon's high ecliptic latitude makes 128 move its rise by minutes."""
        jd_start = julday(2005, 3, 10, 0)
        base = CALC_RISE | BIT_DISC_CENTER | BIT_NO_REFRACTION
        _, t_plain = rise_trans(jd_start, MOON, base, self.ROME)
        _, t_nolat = rise_trans(jd_start, MOON, base | BIT_GEOCTR_NO_ECL_LAT, self.ROME)
        assert t_plain[0] > 0 and t_nolat[0] > 0
        shift_minutes = abs(t_nolat[0] - t_plain[0]) * 24 * 60
        # 2005-03-10 the Moon sits well off the ecliptic: the shift is ~25 min.
        assert shift_minutes > 5.0

    def test_sun_effect_is_small_but_nonzero(self):
        """The Sun rides the ecliptic, so 128 barely moves it (only parallax)."""
        jd_start = julday(2000, 6, 21, 0)
        base = CALC_RISE | BIT_DISC_CENTER | BIT_NO_REFRACTION
        _, t_plain = rise_trans(jd_start, SUN, base, self.ROME)
        _, t_nolat = rise_trans(jd_start, SUN, base | BIT_GEOCTR_NO_ECL_LAT, self.ROME)
        shift_seconds = abs(t_nolat[0] - t_plain[0]) * 86400.0
        assert shift_seconds < 5.0

    def test_no_ecl_lat_ignored_for_transits(self):
        """Measured reference behavior: 128 does not affect meridian transits."""
        jd_start = julday(2005, 3, 10, 0)
        for body in (SUN, MOON, MARS):
            for event in (CALC_MTRANSIT, CALC_ITRANSIT):
                _, t0 = rise_trans(jd_start, body, event, self.ROME)
                _, t1 = rise_trans(
                    jd_start, body, event | BIT_GEOCTR_NO_ECL_LAT, self.ROME
                )
                assert t1[0] == pytest.approx(t0[0], abs=1e-9)

    def test_no_ecl_lat_applies_to_twilight(self):
        """Measured reference behavior: 128 DOES shift the Sun's twilight."""
        jd_start = julday(2000, 6, 21, 0)
        _, t0 = rise_trans(jd_start, SUN, CALC_RISE | BIT_CIVIL_TWILIGHT, self.ROME)
        _, t1 = rise_trans(
            jd_start,
            SUN,
            CALC_RISE | BIT_CIVIL_TWILIGHT | BIT_GEOCTR_NO_ECL_LAT,
            self.ROME,
        )
        assert t0[0] != pytest.approx(t1[0], abs=1e-9)

    def test_event_lands_on_latitude_zeroed_geocentric_direction(self):
        """At the 896 event the latitude-zeroed geocentric direction is on the
        horizon: altitude ~ 0 (disc center, no refraction, flat horizon)."""
        from libephemeris import calc_ut, azalt
        from libephemeris.utils import ECL2HOR

        jd_start = julday(2005, 3, 10, 0)
        for body in (MOON, MARS):
            for event in (CALC_RISE, CALC_SET):
                r, tret = rise_trans(
                    jd_start, body, event | BIT_HINDU_RISING, self.ROME
                )
                assert r == 0
                t = tret[0]
                pos, _ = calc_ut(t, body)  # geocentric apparent ecliptic-of-date
                # azalt geopos order is (lon, lat, alt); latitude zeroed.
                _az, alt_true, _app = azalt(
                    t, ECL2HOR, tuple(self.ROME), 0.0, 0.0, (pos[0], 0.0, pos[2])
                )
                assert alt_true == pytest.approx(0.0, abs=0.01)

    def test_high_lat_circumpolar_decision_uses_projected_declination(self):
        """The circumpolar pre-check must use the SAME latitude-zeroed
        declination as the solve when BIT_GEOCTR_NO_ECL_LAT is set.

        At 64.13N in 1988 the Moon's TRUE declination (~+28.3 deg) makes it
        circumpolar (never sets -> plain rise res=-2), but its ecliptic-
        latitude-zeroed projected place (dec ~+23.4 deg) does rise and set.
        Judging circumpolarity on the true declination while solving on the
        projected place spuriously reported res=-2 for the projected event;
        the projected event exists and must be returned. Backend agnostic
        (leb and skyfield agree to the position precision).
        """
        from libephemeris import calc_ut, azalt
        from libephemeris.utils import ECL2HOR

        reykjavik = [-21.9, 64.13, 60.0]
        jd_start = 2447327.0  # 1988: Moon near maximum ecliptic latitude
        # Plain (true) place: the Moon is circumpolar here.
        r_plain, _ = rise_trans(jd_start, MOON, CALC_RISE, reykjavik)
        assert r_plain == -2
        # Projected place (128 alone and the full 896 Hindu combo): real event.
        for mask in (BIT_GEOCTR_NO_ECL_LAT, BIT_HINDU_RISING):
            r, tret = rise_trans(jd_start, MOON, CALC_RISE | mask, reykjavik)
            assert r == 0
            assert tret[0] > jd_start
        # For 896 (disc center, no refraction) the latitude-zeroed geocentric
        # direction is on the flat horizon at the returned time.
        _r, tret = rise_trans(jd_start, MOON, CALC_RISE | BIT_HINDU_RISING, reykjavik)
        pos, _ = calc_ut(tret[0], MOON)
        _az, alt_true, _app = azalt(
            tret[0], ECL2HOR, tuple(reykjavik), 0.0, 0.0, (pos[0], 0.0, pos[2])
        )
        assert alt_true == pytest.approx(0.0, abs=0.01)


class TestRiseTransExtendedBodies:
    """rise_trans/rise_trans_true_hor for bodies outside the classical
    _PLANET_MAP (lunar nodes and apsides, Chiron, the main-belt asteroids,
    numbered minor planets). These are placeable by calc_ut on either backend,
    so the reference computes their rise, set and transit; the engine routes
    them through the shared calc_ut pipeline. The properties below are backend
    agnostic (the suite runs under both Skyfield and sealed LEB) and reference
    free.
    """

    ROME = [12.5, 41.9028, 0.0]
    # Mean/True Node, Mean/Oscu Apogee, Chiron, Ceres, Vesta, Eros (10433).
    BODIES = [10, 11, 12, 13, 15, 17, 20, 10000 + 433]

    def test_all_events_resolve_without_raising(self):
        """Rise, set and both transits resolve (found or circumpolar) for every
        out-of-map body, on whichever backend the suite selected."""
        jd_start = julday(2015, 3, 20, 0)
        for body in self.BODIES:
            for event in (CALC_RISE, CALC_SET, CALC_MTRANSIT, CALC_ITRANSIT):
                res, tret = rise_trans(jd_start, body, event, self.ROME)
                assert res in (0, -2)
                if res == 0:
                    # Event within the next couple of days of the start.
                    assert 0.0 < tret[0] - jd_start < 2.5

    def test_rise_altitude_is_topocentric_horizon(self):
        """At the rise instant the body's TOPOCENTRIC true altitude sits on the
        horizon (no refraction => 0). This fails if the engine used the
        geocentric place: the ~3" parallax of a near asteroid would leave a
        residual altitude and shift the time by ~0.2 s.
        """
        from libephemeris import azalt, calc_ut, set_topo, FLG_TOPOCTR
        from libephemeris.state import get_topo
        from libephemeris.utils import ECL2HOR

        jd_start = julday(2015, 3, 20, 0)
        lon, lat, alt = self.ROME
        saved = get_topo()
        try:
            for body in self.BODIES:
                res, tret = rise_trans(
                    jd_start, body, CALC_RISE | BIT_NO_REFRACTION, self.ROME
                )
                if res != 0:
                    continue
                set_topo(lon, lat, alt)
                pos, _ = calc_ut(tret[0], body, FLG_TOPOCTR)
                _az, alt_true, _app = azalt(
                    tret[0],
                    ECL2HOR,
                    (lon, lat, alt),
                    0.0,
                    0.0,
                    (pos[0], pos[1], pos[2]),
                )
                assert alt_true == pytest.approx(0.0, abs=0.01)
        finally:
            # Restore the previous observer so global state does not leak.
            if saved is None:
                import libephemeris.state as _st

                _st._TOPO = None
            else:
                set_topo(
                    saved.longitude.degrees,
                    saved.latitude.degrees,
                    saved.elevation.m,
                )

    def test_transit_places_body_on_the_meridian(self):
        """At MTRANSIT the body's hour angle is ~0; at ITRANSIT it is ~180."""
        from libephemeris import calc_ut, sidtime, FLG_EQUATORIAL

        jd_start = julday(2015, 3, 20, 0)
        lon = self.ROME[0]
        for body in self.BODIES:
            for event, target_ha in ((CALC_MTRANSIT, 0.0), (CALC_ITRANSIT, 180.0)):
                res, tret = rise_trans(jd_start, body, event, self.ROME)
                assert res == 0
                eq, _ = calc_ut(tret[0], body, FLG_EQUATORIAL)
                lst_deg = (sidtime(tret[0]) * 15.0 + lon) % 360.0
                ha = ((lst_deg - eq[0]) - target_ha + 180.0) % 360.0 - 180.0
                assert abs(ha) < 0.02

    def test_true_hor_matches_zero_horizon_rise_trans(self):
        """rise_trans is the zero-horizon specialization of
        rise_trans_true_hor for these bodies too."""
        from libephemeris import rise_trans_true_hor

        jd_start = julday(2015, 3, 20, 0)
        for body in (15, 17, 10):
            r1, t1 = rise_trans(jd_start, body, CALC_RISE, self.ROME)
            r2, t2 = rise_trans_true_hor(
                jd_start, body, CALC_RISE, self.ROME, 0.0, 0.0, 0.0
            )
            assert r1 == r2
            if r1 == 0:
                assert t1[0] == pytest.approx(t2[0], abs=1e-6)

    def test_unplaceable_body_raises_typed_error(self):
        """A body no backend can place (a planetary moon with no registered
        SPK) raises the unified illegal-body error: both an Error and a
        ValueError carrying 'illegal planet number', identically on every
        backend."""
        from libephemeris.exceptions import Error

        jd_start = julday(2015, 3, 20, 0)
        with pytest.raises(Error, match="illegal planet number"):
            rise_trans(jd_start, 9999, CALC_RISE, self.ROME)
        with pytest.raises(ValueError, match="illegal planet number"):
            rise_trans(jd_start, 9999, CALC_RISE, self.ROME)


class TestRiseTransBackendEquivalence:
    """rise_trans must be backend-independent: the sealed LEB backend and the
    Skyfield backend evaluate the horizon geometry through the same
    calc_ut(FLG_TOPOCTR)/fixstar_ut + azalt chain and the same meridian-transit
    RA/Dec, so their event times must agree to the position precision. This is
    a pure equivalence property - no reference values are used.

    Regression guard for the earlier split, where the LEB path used
    fast_calc._topo_ecliptic + azalt while the Skyfield path used a native
    Skyfield topocentric altaz frame. Those chains diverged by ~0.3-0.5" at the
    horizon threshold even when the two backends agreed on the position to
    ~0.0003", which - amplified by the loose bisection exit - drove rise/set
    times apart by up to ~0.11 s (e.g. Moon 2081, Venus 2071).
    """

    # A grid of bodies spanning the classical planets (SUN..SATURN, with
    # SATURN exercising the planet-centre-vs-barycentre offset the Skyfield
    # backend applies) plus the mean lunar node, an out-of-_PLANET_MAP point
    # routed through the same calc_ut path. All are analytical or DE440-core
    # bodies whose data is identical on both backends; SPK-fitted minor bodies
    # (Chiron/Ceres/...) are deliberately excluded because the test suite runs
    # with strict precision and SPK auto-download disabled (see conftest), so
    # the Skyfield path would use the Keplerian fallback while the LEB path
    # uses its fitted Chebyshev coefficients - a data-layer difference, not a
    # rise/set backend-equivalence one.
    _BODIES = [SUN, MOON, 2, 3, MARS, JUPITER, 6, 10]  # 2=Mer 3=Ven 6=Sat 10=node
    _EVENTS = [CALC_RISE, CALC_SET, CALC_MTRANSIT, CALC_ITRANSIT]
    _GEOS = [
        [12.4964, 41.9028, 0.0],  # Rome
        [-58.38, -34.60, 25.0],  # Buenos Aires (southern)
        [103.82, 1.35, 15.0],  # Singapore (equatorial)
    ]
    _DATES = [
        julday(1950, 3, 20, 0.0),
        julday(2001, 6, 15, 0.0),
        julday(2071, 10, 3, 0.0),
    ]

    def _run_all_modes(self, jd, body, event, geo):
        """rise_trans in leb then skyfield, restoring the mode. Returns
        {mode: (retflag, jd_event)} or None if the leb backend is unavailable."""
        from libephemeris import set_calc_mode, get_calc_mode
        from libephemeris.state import get_leb_reader

        saved = get_calc_mode()
        out = {}
        try:
            for mode in ("leb", "skyfield"):
                set_calc_mode(mode)
                if mode == "leb" and get_leb_reader() is None:
                    return None  # no LEB data provisioned in this environment
                rf, tret = rise_trans(jd, body, event, geo)
                out[mode] = (rf, tret[0])
        finally:
            set_calc_mode(saved)
        return out

    def test_leb_equals_skyfield_on_grid(self):
        """leb and skyfield agree on the return flag everywhere and, when an
        event is found, on the time to well under the ~0.11 s pre-fix split."""
        # Tolerance: 0.02 s comfortably passes the residual floor (the
        # Saturn/Pluto planet-centre offset is ~0.0045 s of geocentric
        # longitude, everything else <0.001 s) while still catching the
        # 0.055-0.11 s regression this test guards against.
        tol_days = 0.02 / 86400.0
        checked = 0
        worst = 0.0
        for jd in self._DATES:
            for body in self._BODIES:
                for event in self._EVENTS:
                    for geo in self._GEOS:
                        res = self._run_all_modes(jd, body, event, geo)
                        if res is None:
                            pytest.skip("LEB backend not provisioned")
                        (rf_leb, t_leb) = res["leb"]
                        (rf_sky, t_sky) = res["skyfield"]
                        assert rf_leb == rf_sky, (
                            f"return-flag split body={body} event={event} "
                            f"geo={geo} jd={jd}: leb={rf_leb} sky={rf_sky}"
                        )
                        if rf_leb == 0:
                            checked += 1
                            d = abs(t_leb - t_sky)
                            worst = max(worst, d)
                            assert d < tol_days, (
                                f"backend split body={body} event={event} "
                                f"geo={geo} jd={jd}: "
                                f"{(t_leb - t_sky) * 86400:+.5f} s"
                            )
        assert checked > 0, "no events found on the grid"
