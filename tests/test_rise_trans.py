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
        """Test that invalid planet ID raises ValueError."""
        jd_start = julday(2024, 6, 21, 0)

        with pytest.raises(ValueError, match="illegal planet number"):
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
