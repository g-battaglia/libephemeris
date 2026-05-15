"""
Tests for coordinate transformation flags.

This module verifies that all coordinate transformation flags work correctly:
- FLG_EQUATORIAL: Equatorial coordinates (RA/Dec) instead of ecliptic
- FLG_J2000: J2000.0 reference frame (aequinox of J2000)
- FLG_ICRS: International Celestial Reference System
- FLG_SIDEREAL: Sidereal zodiac (requires ayanamsha mode to be set)
- Ecliptic of date (default): True ecliptic and equinox of date

Each test compares libephemeris results with pyswisseph to ensure 1:1 compatibility.
"""

import pytest
import math
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
    FLG_SWIEPH,
    FLG_SPEED,
    FLG_EQUATORIAL,
    FLG_J2000,
    FLG_ICRS,
    FLG_SIDEREAL,
    SIDM_LAHIRI,
    SIDM_FAGAN_BRADLEY,
)


# Test data: various Julian Days spanning different epochs
TEST_DATES = [
    2451545.0,  # J2000.0 (2000-01-01 12:00 TT)
    2440587.5,  # Unix epoch (1970-01-01 00:00)
    2460000.0,  # Recent date (2023)
    2415020.0,  # J1900.0
    2433282.4235,  # B1950.0
]

# Planets available in basic DE440 ephemeris (barycenters only for outer planets)
# Note: Mars center (499) may not be available in all DE440 variants
# The fallback mechanism uses barycenter (4) when planet center is unavailable
ALL_PLANETS = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    # Skip Mars - planet center (499) not available in standard DE440
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
    (URANUS, "Uranus"),
    (NEPTUNE, "Neptune"),
    (PLUTO, "Pluto"),
]


def angle_diff(a1: float, a2: float) -> float:
    """Calculate the smallest difference between two angles (handles wrap-around)."""
    diff = abs(a1 - a2)
    if diff > 180:
        diff = 360 - diff
    return diff


class TestEclipticOfDate:
    """Tests for default ecliptic coordinates (equinox of date)."""

    @pytest.mark.parametrize("jd", TEST_DATES)
    @pytest.mark.parametrize("planet_id,planet_name", ALL_PLANETS)
    def test_ecliptic_of_date_vs_swisseph(self, jd, planet_id, planet_name):
        """Test ecliptic coordinates (default) match pyswisseph."""
        flags = FLG_SWIEPH | FLG_SPEED

        res_swe, _ = swe.calc_ut(jd, planet_id, flags)
        res_lib, _ = ephem.calc_ut(jd, planet_id, flags)

        # Compare longitude
        lon_diff = angle_diff(res_swe[0], res_lib[0])
        assert lon_diff < 0.001, (
            f"{planet_name} at JD {jd}: longitude diff {lon_diff} >= 0.001 deg "
            f"(swe={res_swe[0]:.6f}, lib={res_lib[0]:.6f})"
        )

        # Compare latitude
        lat_diff = abs(res_swe[1] - res_lib[1])
        assert lat_diff < 0.001, (
            f"{planet_name} at JD {jd}: latitude diff {lat_diff} >= 0.001 deg "
            f"(swe={res_swe[1]:.6f}, lib={res_lib[1]:.6f})"
        )

        # Compare distance (relax tolerance for Pluto at distant dates)
        dist_diff = abs(res_swe[2] - res_lib[2])
        assert dist_diff < 0.0002, (
            f"{planet_name} at JD {jd}: distance diff {dist_diff} >= 0.0002 AU "
            f"(swe={res_swe[2]:.8f}, lib={res_lib[2]:.8f})"
        )

    def test_ecliptic_sun_velocity(self, standard_jd):
        """Test that Sun's ecliptic velocity is ~1 deg/day."""
        flags = FLG_SWIEPH | FLG_SPEED

        res_lib, _ = ephem.calc_ut(standard_jd, SUN, flags)

        # Sun moves approximately 1 degree per day
        assert 0.9 < res_lib[3] < 1.1, (
            f"Sun's longitude velocity should be ~1 deg/day, got {res_lib[3]}"
        )


class TestEquatorialCoordinates:
    """Tests for FLG_EQUATORIAL (RA/Dec instead of ecliptic lon/lat)."""

    @pytest.mark.parametrize("jd", TEST_DATES)
    @pytest.mark.parametrize("planet_id,planet_name", ALL_PLANETS)
    def test_equatorial_coordinates_vs_swisseph(self, jd, planet_id, planet_name):
        """Test equatorial (RA/Dec) coordinates match pyswisseph."""
        flags = FLG_SWIEPH | FLG_EQUATORIAL | FLG_SPEED

        res_swe, _ = swe.calc_ut(jd, planet_id, flags)
        res_lib, _ = ephem.calc_ut(jd, planet_id, flags)

        # Compare Right Ascension (stored in position[0])
        ra_diff = angle_diff(res_swe[0], res_lib[0])
        assert ra_diff < 0.001, (
            f"{planet_name} at JD {jd}: RA diff {ra_diff} >= 0.001 deg "
            f"(swe={res_swe[0]:.6f}, lib={res_lib[0]:.6f})"
        )

        # Compare Declination (stored in position[1])
        dec_diff = abs(res_swe[1] - res_lib[1])
        assert dec_diff < 0.001, (
            f"{planet_name} at JD {jd}: Dec diff {dec_diff} >= 0.001 deg "
            f"(swe={res_swe[1]:.6f}, lib={res_lib[1]:.6f})"
        )

        # Compare distance (should be same as ecliptic, relax for Pluto at distant dates)
        dist_diff = abs(res_swe[2] - res_lib[2])
        assert dist_diff < 0.0002, (
            f"{planet_name} at JD {jd}: distance diff {dist_diff} >= 0.0002 AU"
        )

    def test_equatorial_ra_range(self, standard_jd):
        """Test that RA values are in valid range [0, 360)."""
        flags = FLG_SWIEPH | FLG_EQUATORIAL

        for planet_id, planet_name in ALL_PLANETS:
            res_lib, _ = ephem.calc_ut(standard_jd, planet_id, flags)
            ra = res_lib[0]
            assert 0 <= ra < 360, f"{planet_name} RA {ra} out of range [0, 360)"

    def test_equatorial_dec_range(self, standard_jd):
        """Test that Dec values are in valid range [-90, 90]."""
        flags = FLG_SWIEPH | FLG_EQUATORIAL

        for planet_id, planet_name in ALL_PLANETS:
            res_lib, _ = ephem.calc_ut(standard_jd, planet_id, flags)
            dec = res_lib[1]
            assert -90 <= dec <= 90, f"{planet_name} Dec {dec} out of range [-90, 90]"

    def test_equatorial_differs_from_ecliptic(self, standard_jd):
        """Test that equatorial and ecliptic coordinates are different."""
        ecliptic_flags = FLG_SWIEPH
        equatorial_flags = FLG_SWIEPH | FLG_EQUATORIAL

        # Jupiter should have noticeably different coords in the two systems
        res_ecl, _ = ephem.calc_ut(standard_jd, JUPITER, ecliptic_flags)
        res_equ, _ = ephem.calc_ut(standard_jd, JUPITER, equatorial_flags)

        # The coordinates should differ (unless planet is at a special point)
        # Using a relaxed check since at certain positions they could be close
        lon_diff = angle_diff(res_ecl[0], res_equ[0])
        lat_diff = abs(res_ecl[1] - res_equ[1])

        # At least one coordinate should differ noticeably
        assert lon_diff > 0.01 or lat_diff > 0.01, (
            "Equatorial and ecliptic coordinates should differ"
        )


class TestJ2000Coordinates:
    """Tests for FLG_J2000 (J2000.0 reference frame)."""

    @pytest.mark.parametrize("jd", TEST_DATES)
    @pytest.mark.parametrize("planet_id,planet_name", ALL_PLANETS)
    def test_j2000_ecliptic_vs_swisseph(self, jd, planet_id, planet_name):
        """Test J2000 ecliptic coordinates match pyswisseph."""
        flags = FLG_SWIEPH | FLG_J2000 | FLG_SPEED

        res_swe, _ = swe.calc_ut(jd, planet_id, flags)
        res_lib, _ = ephem.calc_ut(jd, planet_id, flags)

        # Compare longitude
        lon_diff = angle_diff(res_swe[0], res_lib[0])
        assert lon_diff < 0.001, (
            f"{planet_name} at JD {jd}: J2000 ecliptic lon diff {lon_diff} >= 0.001 deg "
            f"(swe={res_swe[0]:.6f}, lib={res_lib[0]:.6f})"
        )

        # Compare latitude
        lat_diff = abs(res_swe[1] - res_lib[1])
        assert lat_diff < 0.001, (
            f"{planet_name} at JD {jd}: J2000 ecliptic lat diff {lat_diff} >= 0.001 deg"
        )

    @pytest.mark.parametrize("jd", TEST_DATES)
    @pytest.mark.parametrize("planet_id,planet_name", ALL_PLANETS)
    def test_j2000_equatorial_vs_swisseph(self, jd, planet_id, planet_name):
        """Test J2000 equatorial coordinates (ICRS-like) match pyswisseph."""
        flags = FLG_SWIEPH | FLG_J2000 | FLG_EQUATORIAL | FLG_SPEED

        res_swe, _ = swe.calc_ut(jd, planet_id, flags)
        res_lib, _ = ephem.calc_ut(jd, planet_id, flags)

        # Compare RA
        ra_diff = angle_diff(res_swe[0], res_lib[0])
        assert ra_diff < 0.001, (
            f"{planet_name} at JD {jd}: J2000 RA diff {ra_diff} >= 0.001 deg "
            f"(swe={res_swe[0]:.6f}, lib={res_lib[0]:.6f})"
        )

        # Compare Dec
        dec_diff = abs(res_swe[1] - res_lib[1])
        assert dec_diff < 0.001, (
            f"{planet_name} at JD {jd}: J2000 Dec diff {dec_diff} >= 0.001 deg"
        )

    def test_j2000_differs_from_date(self, standard_jd):
        """Test that J2000 and date coordinates differ due to precession."""
        # Use a date far from J2000 to see precession effect
        jd_2100 = 2488070.0  # ~2100

        date_flags = FLG_SWIEPH
        j2000_flags = FLG_SWIEPH | FLG_J2000

        res_date, _ = ephem.calc_ut(jd_2100, SUN, date_flags)
        res_j2000, _ = ephem.calc_ut(jd_2100, SUN, j2000_flags)

        # Due to precession, coordinates should differ by ~1.4 deg/century
        # Over 100 years, expect ~1.4 degree difference
        lon_diff = angle_diff(res_date[0], res_j2000[0])
        assert lon_diff > 0.5, (
            f"J2000 and date coordinates should differ due to precession, "
            f"got diff={lon_diff:.4f} deg"
        )

    def test_j2000_at_epoch_equals_date(self):
        """Test that at J2000.0 epoch, J2000 coords approximately equal date coords."""
        jd_j2000 = 2451545.0  # J2000.0 epoch

        date_flags = FLG_SWIEPH
        j2000_flags = FLG_SWIEPH | FLG_J2000

        res_date, _ = ephem.calc_ut(jd_j2000, SUN, date_flags)
        res_j2000, _ = ephem.calc_ut(jd_j2000, SUN, j2000_flags)

        # At J2000.0, the difference should be minimal (only nutation)
        lon_diff = angle_diff(res_date[0], res_j2000[0])
        assert lon_diff < 0.02, (
            f"At J2000.0, J2000 and date coords should be very close, "
            f"got diff={lon_diff:.4f} deg"
        )


class TestICRSCoordinates:
    """Tests for FLG_ICRS (International Celestial Reference System)."""

    @pytest.mark.parametrize("jd", TEST_DATES)
    @pytest.mark.parametrize("planet_id,planet_name", ALL_PLANETS[:5])  # Test subset
    def test_icrs_vs_swisseph(self, jd, planet_id, planet_name):
        """Test ICRS coordinates match pyswisseph."""
        flags = FLG_SWIEPH | FLG_ICRS | FLG_EQUATORIAL

        res_swe, _ = swe.calc_ut(jd, planet_id, flags)
        res_lib, _ = ephem.calc_ut(jd, planet_id, flags)

        # Compare RA
        ra_diff = angle_diff(res_swe[0], res_lib[0])
        assert ra_diff < 0.001, (
            f"{planet_name} at JD {jd}: ICRS RA diff {ra_diff} >= 0.001 deg "
            f"(swe={res_swe[0]:.6f}, lib={res_lib[0]:.6f})"
        )

        # Compare Dec
        dec_diff = abs(res_swe[1] - res_lib[1])
        assert dec_diff < 0.001, (
            f"{planet_name} at JD {jd}: ICRS Dec diff {dec_diff} >= 0.001 deg"
        )

    def test_icrs_is_fixed_reference_frame(self):
        """Test that ICRS is a fixed (non-rotating) reference frame."""
        # ICRS should give the same frame orientation regardless of date
        # (unlike equinox of date which precesses)
        jd1 = 2451545.0  # J2000.0
        jd2 = 2488070.0  # ~2100

        icrs_flags = FLG_SWIEPH | FLG_ICRS | FLG_EQUATORIAL

        # For a star-like distant object, ICRS coords should barely change
        # Using Neptune as a "slow" object
        res1, _ = ephem.calc_ut(jd1, NEPTUNE, icrs_flags)
        res2, _ = ephem.calc_ut(jd2, NEPTUNE, icrs_flags)

        # Neptune moves slowly; in ICRS, the frame doesn't precess
        # The position change is due to Neptune's orbital motion, not precession
        # This is more of a sanity check that ICRS is working
        assert res1[0] != res2[0], (
            "Neptune should have different ICRS RA at different dates"
        )


class TestSiderealCoordinates:
    """Tests for FLG_SIDEREAL (sidereal zodiac)."""

    @pytest.mark.parametrize("jd", TEST_DATES)
    @pytest.mark.parametrize("planet_id,planet_name", ALL_PLANETS)
    def test_sidereal_lahiri_vs_swisseph(self, jd, planet_id, planet_name):
        """Test sidereal (Lahiri) coordinates match pyswisseph."""
        # Set Lahiri ayanamsha
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.set_sid_mode(SIDM_LAHIRI)

        flags = FLG_SWIEPH | FLG_SIDEREAL | FLG_SPEED

        res_swe, _ = swe.calc_ut(jd, planet_id, flags)
        res_lib, _ = ephem.calc_ut(jd, planet_id, flags)

        # Compare sidereal longitude
        lon_diff = angle_diff(res_swe[0], res_lib[0])
        assert lon_diff < 0.01, (
            f"{planet_name} at JD {jd}: sidereal lon diff {lon_diff} >= 0.01 deg "
            f"(swe={res_swe[0]:.6f}, lib={res_lib[0]:.6f})"
        )

        # Compare latitude (should be same as tropical)
        lat_diff = abs(res_swe[1] - res_lib[1])
        assert lat_diff < 0.001, (
            f"{planet_name} at JD {jd}: sidereal lat diff {lat_diff} >= 0.001 deg"
        )

    @pytest.mark.parametrize("jd", TEST_DATES)
    def test_sidereal_fagan_bradley_vs_swisseph(self, jd):
        """Test sidereal (Fagan-Bradley) coordinates match pyswisseph."""
        swe.set_sid_mode(swe.SIDM_FAGAN_BRADLEY)
        ephem.set_sid_mode(SIDM_FAGAN_BRADLEY)

        flags = FLG_SWIEPH | FLG_SIDEREAL | FLG_SPEED

        res_swe, _ = swe.calc_ut(jd, SUN, flags)
        res_lib, _ = ephem.calc_ut(jd, SUN, flags)

        lon_diff = angle_diff(res_swe[0], res_lib[0])
        assert lon_diff < 0.01, (
            f"Sun at JD {jd}: Fagan-Bradley sidereal lon diff {lon_diff} >= 0.01 deg"
        )

    def test_sidereal_differs_from_tropical(self, standard_jd):
        """Test that sidereal coordinates differ from tropical by ~24 degrees."""
        ephem.set_sid_mode(SIDM_LAHIRI)

        tropical_flags = FLG_SWIEPH
        sidereal_flags = FLG_SWIEPH | FLG_SIDEREAL

        res_trop, _ = ephem.calc_ut(standard_jd, SUN, tropical_flags)
        res_sid, _ = ephem.calc_ut(standard_jd, SUN, sidereal_flags)

        # Lahiri ayanamsha at J2000 is ~23.85 degrees
        lon_diff = angle_diff(res_trop[0], res_sid[0])
        assert 22 < lon_diff < 26, (
            f"Sidereal should differ from tropical by ~24 deg, got {lon_diff:.2f} deg"
        )

    def test_sidereal_velocity_matches_tropical(self, standard_jd):
        """Test that sidereal velocity approximately matches tropical velocity."""
        ephem.set_sid_mode(SIDM_LAHIRI)

        tropical_flags = FLG_SWIEPH | FLG_SPEED
        sidereal_flags = FLG_SWIEPH | FLG_SIDEREAL | FLG_SPEED

        res_trop, _ = ephem.calc_ut(standard_jd, JUPITER, tropical_flags)
        res_sid, _ = ephem.calc_ut(standard_jd, JUPITER, sidereal_flags)

        # Velocity should be nearly the same (ayanamsha rate is ~50"/century)
        vel_diff = abs(res_trop[3] - res_sid[3])
        assert vel_diff < 0.001, (
            f"Sidereal and tropical velocity should match closely, diff={vel_diff}"
        )


class TestCombinedFlags:
    """Tests for combinations of coordinate transformation flags."""

    def test_j2000_equatorial(self, standard_jd):
        """Test J2000 + EQUATORIAL flags work together."""
        flags = FLG_SWIEPH | FLG_J2000 | FLG_EQUATORIAL | FLG_SPEED

        res_swe, _ = swe.calc_ut(standard_jd, SUN, flags)
        res_lib, _ = ephem.calc_ut(standard_jd, SUN, flags)

        ra_diff = angle_diff(res_swe[0], res_lib[0])
        dec_diff = abs(res_swe[1] - res_lib[1])

        assert ra_diff < 0.001, f"J2000+EQ RA diff {ra_diff} >= 0.001"
        assert dec_diff < 0.001, f"J2000+EQ Dec diff {dec_diff} >= 0.001"

    def test_sidereal_with_speed(self, standard_jd):
        """Test SIDEREAL + SPEED flags work together."""
        ephem.set_sid_mode(SIDM_LAHIRI)
        swe.set_sid_mode(swe.SIDM_LAHIRI)

        flags = FLG_SWIEPH | FLG_SIDEREAL | FLG_SPEED

        res_swe, _ = swe.calc_ut(standard_jd, JUPITER, flags)
        res_lib, _ = ephem.calc_ut(standard_jd, JUPITER, flags)

        # Check position
        lon_diff = angle_diff(res_swe[0], res_lib[0])
        assert lon_diff < 0.01, f"Sidereal+Speed lon diff {lon_diff} >= 0.01"

        # Check velocity
        vel_diff = abs(res_swe[3] - res_lib[3])
        assert vel_diff < 0.01, f"Sidereal+Speed velocity diff {vel_diff} >= 0.01"

    def test_equatorial_with_speed(self, standard_jd):
        """Test EQUATORIAL + SPEED flags work together."""
        flags = FLG_SWIEPH | FLG_EQUATORIAL | FLG_SPEED

        res_swe, _ = swe.calc_ut(standard_jd, MOON, flags)
        res_lib, _ = ephem.calc_ut(standard_jd, MOON, flags)

        # Check position
        ra_diff = angle_diff(res_swe[0], res_lib[0])
        assert ra_diff < 0.001, f"EQ+Speed RA diff {ra_diff} >= 0.001"

        # Check velocity (Moon RA motion is fast)
        vel_diff = abs(res_swe[3] - res_lib[3])
        assert vel_diff < 0.1, f"EQ+Speed velocity diff {vel_diff} >= 0.1"


class TestEdgeCases:
    """Edge case tests for coordinate transformations."""

    def test_pole_coordinates(self):
        """Test that coordinates near poles don't cause errors."""
        jd = 2451545.0
        flags = FLG_SWIEPH | FLG_EQUATORIAL

        # Calculate position of Pluto (can have high declination)
        res_lib, _ = ephem.calc_ut(jd, PLUTO, flags)

        # Should return valid coordinates
        assert -90 <= res_lib[1] <= 90, f"Declination {res_lib[1]} out of range"
        assert 0 <= res_lib[0] < 360, f"RA {res_lib[0]} out of range"

    def test_all_flags_return_tuples(self, standard_jd):
        """Test that all flag combinations return proper 6-element tuples."""
        flag_combos = [
            FLG_SWIEPH,
            FLG_SWIEPH | FLG_EQUATORIAL,
            FLG_SWIEPH | FLG_J2000,
            FLG_SWIEPH | FLG_J2000 | FLG_EQUATORIAL,
            FLG_SWIEPH | FLG_SPEED,
            FLG_SWIEPH | FLG_EQUATORIAL | FLG_SPEED,
        ]

        for flags in flag_combos:
            res, retflag = ephem.calc_ut(standard_jd, SUN, flags)
            assert len(res) == 6, f"Result should have 6 elements for flags {flags}"
            assert all(isinstance(x, float) for x in res), (
                "All elements should be floats"
            )

    def test_distant_dates(self):
        """Test coordinate transformations at distant dates."""
        # Test at distant past and future safely within DE440 range (1549-2650)
        dates = [
            2378497.0,  # 1800-01-01 (safely within DE440)
            2524594.0,  # 2200-01-01 (safely within DE440)
        ]

        for jd in dates:
            for flags in [FLG_SWIEPH, FLG_SWIEPH | FLG_J2000]:
                res_swe, _ = swe.calc_ut(jd, SUN, flags)
                res_lib, _ = ephem.calc_ut(jd, SUN, flags)

                lon_diff = angle_diff(res_swe[0], res_lib[0])
                assert lon_diff < 0.01, (
                    f"Sun at distant date JD {jd}: lon diff {lon_diff} >= 0.01"
                )


class TestCalcFunction:
    """Tests for calc (TT input) vs calc_ut (UT input)."""

    def test_calc_vs_calc_ut(self, standard_jd):
        """Test that calc and calc_ut give consistent results."""
        flags = FLG_SWIEPH | FLG_EQUATORIAL

        # calc takes TT, calc_ut takes UT
        # At J2000.0, the difference is small (~64 seconds)
        res_ut, _ = ephem.calc_ut(standard_jd, SUN, flags)
        res_tt, _ = ephem.calc(standard_jd, SUN, flags)

        # Should be very close (difference due to Delta T is small in position)
        lon_diff = angle_diff(res_ut[0], res_tt[0])
        assert lon_diff < 0.01, f"calc vs calc_ut diff {lon_diff} >= 0.01"

    def test_calc_with_j2000_flag(self, standard_jd):
        """Test calc with J2000 flag."""
        flags = FLG_SWIEPH | FLG_J2000

        res_swe, _ = swe.calc(standard_jd, JUPITER, flags)
        res_lib, _ = ephem.calc(standard_jd, JUPITER, flags)

        lon_diff = angle_diff(res_swe[0], res_lib[0])
        assert lon_diff < 0.001, f"calc J2000 lon diff {lon_diff} >= 0.001"
