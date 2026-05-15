"""Tests for multi-body consistency: positions should be physically consistent."""

from __future__ import annotations

import math
import pytest
import libephemeris as swe
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
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    OSCU_APOG,
    CHIRON,
    EARTH,
    FLG_SWIEPH,
    FLG_SPEED,
    FLG_HELCTR,
    FLG_EQUATORIAL,
)

JD_J2000 = 2451545.0


@pytest.mark.unit
class TestDistanceOrdering:
    """Verify planetary distances follow expected ordering."""

    def test_inner_planets_closer_than_outer(self):
        """Inner planets should generally be closer to Earth than outer planets."""
        flags = FLG_SWIEPH | FLG_SPEED
        bodies = [MERCURY, VENUS, MARS, JUPITER, SATURN]
        distances = []
        for body in bodies:
            result, _ = swe.calc_ut(JD_J2000, body, flags)
            distances.append(result[2])
        # Jupiter should be farther than Mars (generally true)
        mars_idx = bodies.index(MARS)
        jupiter_idx = bodies.index(JUPITER)
        saturn_idx = bodies.index(SATURN)
        assert distances[jupiter_idx] > distances[mars_idx]
        assert distances[saturn_idx] > distances[jupiter_idx]

    def test_outer_planet_distance_ordering(self):
        """Outer planets should be in approximate distance order."""
        flags = FLG_SWIEPH | FLG_SPEED
        bodies = [JUPITER, SATURN, URANUS, NEPTUNE]
        distances = []
        for body in bodies:
            result, _ = swe.calc_ut(JD_J2000, body, flags)
            distances.append(result[2])
        for i in range(1, len(distances)):
            assert distances[i] > distances[i - 1], (
                f"Distance ordering violated at index {i}"
            )

    def test_moon_closest(self):
        """Moon should be the closest body to Earth."""
        flags = FLG_SWIEPH | FLG_SPEED
        moon_result, _ = swe.calc_ut(JD_J2000, MOON, flags)
        moon_dist = moon_result[2]
        # Moon ~0.0026 AU
        assert moon_dist < 0.01, f"Moon distance {moon_dist} AU too large"
        # Compare to Mercury
        merc_result, _ = swe.calc_ut(JD_J2000, MERCURY, flags)
        assert merc_result[2] > moon_dist

    def test_sun_distance_1au(self):
        """Sun distance should be ~1 AU."""
        result, _ = swe.calc_ut(JD_J2000, SUN, FLG_SWIEPH | FLG_SPEED)
        sun_dist = result[2]
        assert 0.98 < sun_dist < 1.02, f"Sun distance {sun_dist} AU not ~1"


@pytest.mark.unit
class TestHeliocentricConsistency:
    """Heliocentric vs geocentric consistency."""

    def test_sun_helio_distance(self):
        """Sun in heliocentric mode should return distance 0 (Sun from itself)."""
        flags = FLG_SWIEPH | FLG_SPEED | FLG_HELCTR
        result, _ = swe.calc_ut(JD_J2000, SUN, flags)
        # Sun heliocentric = Sun from Sun = zero distance
        assert result[2] == pytest.approx(0.0, abs=1e-6)

    def test_earth_helio_valid(self):
        """Earth heliocentric returns a valid position."""
        flags = FLG_SWIEPH | FLG_SPEED | FLG_HELCTR
        result, _ = swe.calc_ut(JD_J2000, EARTH, flags)
        lon, lat, dist = result[0], result[1], result[2]
        assert 0 <= lon < 360
        assert -90 <= lat <= 90
        assert 0.98 < dist < 1.02  # ~1 AU from Sun

    def test_earth_helio_opposite_sun_geo(self):
        """Earth heliocentric longitude should be ~180° from geocentric Sun."""
        sun_geo, _ = swe.calc_ut(JD_J2000, SUN, FLG_SWIEPH | FLG_SPEED)
        earth_helio, _ = swe.calc_ut(
            JD_J2000, EARTH, FLG_SWIEPH | FLG_SPEED | FLG_HELCTR
        )
        diff = abs(earth_helio[0] - sun_geo[0])
        if diff > 180:
            diff = 360 - diff
        assert diff == pytest.approx(180.0, abs=1.0)

    @pytest.mark.parametrize(
        "body,name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_helio_distance_reasonable(self, body, name):
        """Heliocentric distances should be in known ranges."""
        flags = FLG_SWIEPH | FLG_SPEED | FLG_HELCTR
        result, _ = swe.calc_ut(JD_J2000, body, flags)
        dist = result[2]
        expected_ranges = {
            MERCURY: (0.3, 0.5),
            VENUS: (0.7, 0.73),
            MARS: (1.3, 1.7),
            JUPITER: (4.9, 5.5),
            SATURN: (9.0, 10.1),
        }
        lo, hi = expected_ranges[body]
        assert lo < dist < hi, f"{name} helio dist {dist} not in [{lo}, {hi}]"


@pytest.mark.unit
class TestSpeedConsistency:
    """Verify speed values are physically consistent."""

    def test_moon_fastest(self):
        """Moon should have the fastest longitude speed."""
        flags = FLG_SWIEPH | FLG_SPEED
        moon, _ = swe.calc_ut(JD_J2000, MOON, flags)
        sun, _ = swe.calc_ut(JD_J2000, SUN, flags)
        assert abs(moon[3]) > abs(sun[3]), "Moon should be faster than Sun"

    def test_outer_planets_slower(self):
        """Outer planets should have slower apparent speed than inner."""
        flags = FLG_SWIEPH | FLG_SPEED
        mars, _ = swe.calc_ut(JD_J2000, MARS, flags)
        jupiter, _ = swe.calc_ut(JD_J2000, JUPITER, flags)
        saturn, _ = swe.calc_ut(JD_J2000, SATURN, flags)
        # Mean speeds: Mars ~0.5°/d, Jupiter ~0.08°/d, Saturn ~0.03°/d
        # Use absolute values since bodies can be retrograde
        # Check over many dates to avoid retrograde coincidences
        mars_speeds = []
        jup_speeds = []
        sat_speeds = []
        for offset in range(0, 360, 30):
            jd = JD_J2000 + offset
            m, _ = swe.calc_ut(jd, MARS, flags)
            j, _ = swe.calc_ut(jd, JUPITER, flags)
            s, _ = swe.calc_ut(jd, SATURN, flags)
            mars_speeds.append(abs(m[3]))
            jup_speeds.append(abs(j[3]))
            sat_speeds.append(abs(s[3]))
        avg_mars = sum(mars_speeds) / len(mars_speeds)
        avg_jup = sum(jup_speeds) / len(jup_speeds)
        avg_sat = sum(sat_speeds) / len(sat_speeds)
        assert avg_mars > avg_jup > avg_sat

    def test_sun_speed_approximately_1_deg_per_day(self):
        """Sun's apparent speed should be ~1°/day."""
        flags = FLG_SWIEPH | FLG_SPEED
        result, _ = swe.calc_ut(JD_J2000, SUN, flags)
        sun_speed = result[3]
        assert 0.9 < abs(sun_speed) < 1.1, f"Sun speed {sun_speed}°/d not ~1"

    def test_moon_speed_approximately_13_deg_per_day(self):
        """Moon's apparent speed should be ~13°/day."""
        flags = FLG_SWIEPH | FLG_SPEED
        result, _ = swe.calc_ut(JD_J2000, MOON, flags)
        moon_speed = result[3]
        assert 11 < abs(moon_speed) < 15, f"Moon speed {moon_speed}°/d not ~13"

    def test_mean_node_always_retrograde(self):
        """Mean node should always have negative speed (retrograde)."""
        flags = FLG_SWIEPH | FLG_SPEED
        for offset in range(0, 3650, 365):
            jd = JD_J2000 + offset
            result, _ = swe.calc_ut(jd, MEAN_NODE, flags)
            assert result[3] < 0, (
                f"Mean node speed {result[3]} not retrograde at JD {jd}"
            )


@pytest.mark.unit
class TestEquatorialConsistency:
    """Verify equatorial coordinate transformations are consistent."""

    def test_equatorial_declination_range(self):
        """Declination should be in [-90, 90]."""
        flags = FLG_SWIEPH | FLG_SPEED | FLG_EQUATORIAL
        for body in [SUN, MOON, MARS, JUPITER]:
            result, _ = swe.calc_ut(JD_J2000, body, flags)
            dec = result[1]
            assert -90 <= dec <= 90, f"Dec {dec} out of range for body {body}"

    def test_equatorial_ra_range(self):
        """Right ascension should be in [0, 360)."""
        flags = FLG_SWIEPH | FLG_SPEED | FLG_EQUATORIAL
        for body in [SUN, MOON, MARS, JUPITER]:
            result, _ = swe.calc_ut(JD_J2000, body, flags)
            ra = result[0]
            assert 0 <= ra < 360, f"RA {ra} out of range for body {body}"

    def test_sun_max_declination(self):
        """Sun's declination should not exceed ~23.44° (obliquity)."""
        flags = FLG_SWIEPH | FLG_SPEED | FLG_EQUATORIAL
        max_dec = 0
        for offset in range(0, 365):
            jd = JD_J2000 + offset
            result, _ = swe.calc_ut(jd, SUN, flags)
            max_dec = max(max_dec, abs(result[1]))
        assert 23 < max_dec < 24, f"Sun max dec {max_dec}° not near obliquity"

    def test_equatorial_distance_same_as_ecliptic(self):
        """Distance should be the same regardless of coordinate system."""
        for body in [SUN, MOON, MARS]:
            ecl, _ = swe.calc_ut(JD_J2000, body, FLG_SWIEPH | FLG_SPEED)
            equ, _ = swe.calc_ut(
                JD_J2000, body, FLG_SWIEPH | FLG_SPEED | FLG_EQUATORIAL
            )
            assert ecl[2] == pytest.approx(equ[2], rel=1e-8), (
                f"Distance mismatch for body {body}"
            )


@pytest.mark.unit
class TestNodeConsistency:
    """Verify node positions are consistent."""

    def test_true_node_near_mean_node(self):
        """True node should be within ~2° of mean node."""
        flags = FLG_SWIEPH | FLG_SPEED
        mean, _ = swe.calc_ut(JD_J2000, MEAN_NODE, flags)
        true, _ = swe.calc_ut(JD_J2000, TRUE_NODE, flags)
        diff = abs(mean[0] - true[0])
        if diff > 180:
            diff = 360 - diff
        assert diff < 3.0, f"True-mean node difference {diff}° too large"

    def test_node_latitude_near_zero(self):
        """Node latitude should be ~0 (on the ecliptic by definition)."""
        flags = FLG_SWIEPH | FLG_SPEED
        mean, _ = swe.calc_ut(JD_J2000, MEAN_NODE, flags)
        assert abs(mean[1]) < 0.01, f"Mean node lat {mean[1]}° not near zero"

    def test_mean_node_regression_rate(self):
        """Mean node regresses ~19.35° per year."""
        flags = FLG_SWIEPH | FLG_SPEED
        r1, _ = swe.calc_ut(JD_J2000, MEAN_NODE, flags)
        r2, _ = swe.calc_ut(JD_J2000 + 365.25, MEAN_NODE, flags)
        # Node moves ~19.35° per year retrograde
        diff = r1[0] - r2[0]  # Should be positive (retrograde)
        if diff < 0:
            diff += 360
        assert 18 < diff < 21, f"Node regression rate {diff}°/yr not ~19.35"
