"""
Tests for orbit_max_min_true_distance function.

This function calculates the minimum and maximum geocentric distances
during a planet's orbit, useful for determining perigee and apogee.
"""

import pytest
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
    EARTH,
)


class TestOrbitMaxMinTrueDistanceBasic:
    """Test basic functionality of orbit_max_min_true_distance."""

    @pytest.mark.unit
    def test_returns_tuple_of_three_floats(self):
        """orbit_max_min_true_distance should return a tuple of three floats."""
        jd = 2451545.0  # J2000
        result = ephem.orbit_max_min_true_distance(jd, MARS, 0)

        assert isinstance(result, tuple)
        assert len(result) == 3
        assert isinstance(result[0], float)
        assert isinstance(result[1], float)
        assert isinstance(result[2], float)

    @pytest.mark.unit
    def test_min_less_than_max(self):
        """Minimum distance should be less than maximum distance."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(jd, MARS, 0)

        assert min_dist < max_dist

    @pytest.mark.unit
    def test_sun_returns_earth_orbit_distances(self):
        """Sun max/min should reflect Earth's orbital eccentricity (~0.983-1.017 AU)."""
        jd = 2451545.0
        result = ephem.orbit_max_min_true_distance(jd, SUN, 0)

        max_dist, min_dist, true_dist = result
        # Earth's aphelion (~1.017 AU) and perihelion (~0.983 AU)
        assert 1.01 < max_dist < 1.02, f"Sun max dist {max_dist} should be ~1.017 AU"
        assert 0.98 < min_dist < 0.99, f"Sun min dist {min_dist} should be ~0.983 AU"
        assert min_dist < max_dist
        # true_dist is the current Earth-Sun distance (~1 AU)
        assert 0.98 < true_dist < 1.02

    @pytest.mark.unit
    def test_earth_min_zero_max_orbit_diameter(self):
        """Earth is the observer: distance to its own orbit is 0 at minimum
        (coincident) and the orbit diameter (~2 AU) at maximum."""
        jd = 2451545.0
        max_dist, min_dist, _ = ephem.orbit_max_min_true_distance(jd, EARTH, 0)

        assert min_dist == 0.0  # coincident point on the same orbit
        assert 1.99 < max_dist < 2.01  # ~2 * semi-major axis of the EMB orbit


class TestOrbitMaxMinTrueDistanceMoon:
    """Test orbit_max_min_true_distance for Moon."""

    @pytest.mark.unit
    def test_moon_min_is_perigee(self):
        """Moon minimum distance should be approximately perigee (~0.0024 AU)."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(jd, MOON, 0)

        # Moon perigee is about 356,500 km = 0.00238 AU
        assert 0.0020 < min_dist < 0.0026

    @pytest.mark.unit
    def test_moon_max_is_apogee(self):
        """Moon maximum distance should be approximately apogee (~0.0027 AU)."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(jd, MOON, 0)

        # Moon apogee is about 406,700 km = 0.00272 AU
        assert 0.0025 < max_dist < 0.0030


class TestOrbitMaxMinTrueDistanceInnerPlanets:
    """Test orbit_max_min_true_distance for inner planets."""

    @pytest.mark.unit
    def test_mercury_distance_range(self):
        """Mercury min/max distances should be reasonable."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(
            jd, MERCURY, 0
        )

        # Mercury: min ~0.55 AU (inferior conjunction), max ~1.45 AU (superior conjunction)
        assert 0.4 < min_dist < 0.8
        assert 1.2 < max_dist < 1.6

    @pytest.mark.unit
    def test_venus_distance_range(self):
        """Venus min/max distances should be reasonable."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(jd, VENUS, 0)

        # Venus: min ~0.26 AU (inferior conjunction), max ~1.74 AU (superior conjunction)
        assert 0.2 < min_dist < 0.4
        assert 1.6 < max_dist < 1.9


class TestOrbitMaxMinTrueDistanceOuterPlanets:
    """Test orbit_max_min_true_distance for outer planets."""

    @pytest.mark.unit
    def test_mars_distance_range(self):
        """Mars min/max distances should be reasonable."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(jd, MARS, 0)

        # Mars: min ~0.37 AU (opposition at perihelion), max ~2.68 AU (conjunction at aphelion)
        assert 0.3 < min_dist < 0.6
        assert 2.4 < max_dist < 2.9

    @pytest.mark.unit
    def test_jupiter_distance_range(self):
        """Jupiter min/max distances should be reasonable."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(
            jd, JUPITER, 0
        )

        # Jupiter: min ~4.0 AU (opposition), max ~6.5 AU (conjunction)
        assert 3.8 < min_dist < 4.5
        assert 6.0 < max_dist < 6.8

    @pytest.mark.unit
    def test_saturn_distance_range(self):
        """Saturn min/max distances should be reasonable."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(jd, SATURN, 0)

        # Saturn: min ~8.0 AU, max ~11.1 AU
        assert 7.5 < min_dist < 9.0
        assert 10.0 < max_dist < 12.0

    @pytest.mark.unit
    def test_uranus_distance_range(self):
        """Uranus min/max distances should be reasonable."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(jd, URANUS, 0)

        # Uranus: min ~17.3 AU, max ~21.1 AU
        assert 16.5 < min_dist < 19.0
        assert 20.0 < max_dist < 22.0

    @pytest.mark.unit
    def test_neptune_distance_range(self):
        """Neptune min/max distances should be reasonable."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(
            jd, NEPTUNE, 0
        )

        # Neptune: min ~28.8 AU, max ~31.3 AU
        assert 28.0 < min_dist < 30.5
        assert 30.0 < max_dist < 32.5

    @pytest.mark.unit
    def test_pluto_distance_range(self):
        """Pluto min/max distances should be reasonable."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(jd, PLUTO, 0)

        # Pluto: min ~28.8 AU (at perihelion opposition), max ~50.3 AU (at aphelion conjunction)
        assert 28.0 < min_dist < 32.0
        assert 48.0 < max_dist < 52.0


class TestOrbitMaxMinTrueDistanceInvalidObjects:
    """Objects without an orbit raise the reference-format Error."""

    @pytest.mark.unit
    @pytest.mark.parametrize("ipl", [10, 11, 12, 13, 21, 22, -1])
    def test_points_without_orbit_raise(self, ipl):
        from libephemeris.exceptions import Error

        with pytest.raises(Error) as exc:
            ephem.orbit_max_min_true_distance(2451545.0, ipl, 0)
        assert str(exc.value) == (
            f"orbit_max_min_true_distance: body id {ipl} has no orbital elements"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("ipl", [23, 30, 39])
    def test_undefined_block_raises_illegal(self, ipl):
        from libephemeris.exceptions import Error

        with pytest.raises(Error) as exc:
            ephem.orbit_max_min_true_distance(2451545.0, ipl, 0)
        assert str(exc.value) == (
            f"orbit_max_min_true_distance: body id {ipl} is not a defined body"
        )

    @pytest.mark.unit
    def test_sun_is_valid_here(self):
        """Unlike get_orbital_elements, the Sun is a valid target for
        orbit_max_min_true_distance (Earth-orbit range)."""
        max_dist, min_dist, _ = ephem.orbit_max_min_true_distance(2451545.0, SUN, 0)
        assert 1.01 < max_dist < 1.02
        assert 0.98 < min_dist < 0.99


class TestOrbitMaxMinHelctrIgnored:
    """FLG_HELCTR is ignored for the Sun and Moon (hel == geo)."""

    @pytest.mark.unit
    @pytest.mark.parametrize("ipl", [SUN, MOON])
    def test_helctr_matches_geo(self, ipl):
        geo = ephem.orbit_max_min_true_distance(2451545.0, ipl, 0)
        hel = ephem.orbit_max_min_true_distance(2451545.0, ipl, ephem.FLG_HELCTR)
        for g, h in zip(geo, hel):
            assert g == pytest.approx(h, abs=1e-9)


class TestOrbitMaxMinGeocentricFictitious:
    """White Moon (56) / Waldemath (58): geocentric circular constructions."""

    @pytest.mark.unit
    @pytest.mark.parametrize("ipl", [56, 58])
    def test_max_min_true_are_the_model_radius(self, ipl):
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(
            2451545.0, ipl, 0
        )
        # Near-circular geocentric model: the three distances collapse onto the
        # small Earth-relative radius (~0.05 / ~0.007 AU), not the reference's
        # eccentric heliocentric range.
        assert max_dist < 0.1
        assert min_dist == pytest.approx(max_dist, rel=0.05)
        assert min_dist - 1e-9 <= true_dist <= max_dist + 1e-9


class TestOrbitMaxMinTrueDistanceAliases:
    """Test that function aliases work."""

    @pytest.mark.unit
    def test_swe_orbit_max_min_true_distance_alias(self):
        """orbit_max_min_true_distance should be available."""
        result = ephem.orbit_max_min_true_distance(2451545.0, MARS, 0)
        assert len(result) == 3

    @pytest.mark.unit
    def test_orbit_max_min_true_distance_alias(self):
        """orbit_max_min_true_distance should be available."""
        result = ephem.orbit_max_min_true_distance(2451545.0, MARS, 0)
        assert len(result) == 3


class TestOrbitMaxMinTrueDistanceAllPlanets:
    """Test orbit_max_min_true_distance for all planets."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
            (URANUS, "Uranus"),
            (NEPTUNE, "Neptune"),
            (PLUTO, "Pluto"),
        ],
    )
    def test_all_planets_return_positive_distances(self, planet_id, planet_name):
        """All planets should return positive min and max distances."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(
            jd, planet_id, 0
        )

        assert min_dist > 0, f"{planet_name} min distance should be positive"
        assert max_dist > 0, f"{planet_name} max distance should be positive"
        assert min_dist < max_dist, f"{planet_name} min should be less than max"
        assert true_dist > 0, f"{planet_name} true distance should be positive"
        assert min_dist <= true_dist <= max_dist, (
            f"{planet_name} true distance should be between min and max"
        )


class TestOrbitMaxMinTrueDistanceTimeScale:
    """Verify that orbit_max_min_true_distance uses TT (not UT1)."""

    @pytest.mark.unit
    def test_tt_time_scale_consistency(self):
        """The true distance from orbit_max_min_true_distance (TT input)
        should match calc (TT input), not calc_ut (UT input)."""
        jd_tt = 2451545.0  # J2000.0 in TT
        _, _, true_dist = ephem.orbit_max_min_true_distance(jd_tt, MARS, 0)
        # calc takes TT — should match closely
        pos_tt, _ = ephem.calc(jd_tt, MARS, ephem.FLG_SPEED)
        assert abs(true_dist - pos_tt[2]) < 0.001, (
            f"orbit_max_min true_dist ({true_dist:.6f}) should match "
            f"calc distance ({pos_tt[2]:.6f}) for same TT input"
        )


class TestCalcPctrTimeScale:
    """Verify that calc_pctr uses TT (not UT1)."""

    @pytest.mark.unit
    def test_pctr_uses_tt_time_scale(self):
        """calc_pctr should use TT internally, consistent with calc."""
        jd_tt = 2451545.0
        # calc_pctr: position of Moon as seen from Mars
        pos_pctr, _ = ephem.calc_pctr(jd_tt, MOON, MARS, ephem.FLG_SPEED)
        # The result should be a valid ecliptic longitude (not NaN/zero)
        assert 0.0 <= pos_pctr[0] < 360.0
        # Distance should be positive and reasonable (Mars-Moon is 0.5-3 AU)
        assert 0.3 < pos_pctr[2] < 4.0

    @pytest.mark.unit
    def test_pctr_sun_from_mars(self):
        """Sun position as seen from Mars should give ~180° different from
        geocentric Mars position (Mars is roughly opposite the Sun from Earth)."""
        jd_tt = 2451545.0
        pos_pctr, _ = ephem.calc_pctr(jd_tt, SUN, MARS, ephem.FLG_SPEED)
        assert 0.0 <= pos_pctr[0] < 360.0
        # Sun-Mars distance should be roughly 1.5 AU
        assert 1.0 < pos_pctr[2] < 2.0


class TestTrueDistanceSlotIsGeometric:
    """The 3rd slot is the GEOMETRIC (light-time-free) distance.

    Measured reference behavior returns the FLG_TRUEPOS distance in the
    true-distance slot; a missing flag used to leak the apparent,
    light-time-corrected value (up to ~1.6e-4 AU short on Mercury, where
    the radial rate is fastest).
    """

    @pytest.mark.unit
    @pytest.mark.parametrize("ipl", [2, 4, 9])
    def test_slot_matches_truepos_distance(self, ipl):
        jd = 2458849.5
        td = ephem.orbit_max_min_true_distance(jd, ipl, ephem.FLG_SWIEPH)[2]
        true = ephem.calc(jd, ipl, ephem.FLG_SWIEPH | ephem.FLG_TRUEPOS)[0][2]
        apparent = ephem.calc(jd, ipl, ephem.FLG_SWIEPH)[0][2]
        assert td == pytest.approx(true, abs=1e-12)
        # Sanity: the two reductions genuinely differ at this instant, so
        # the assertion above discriminates.
        assert abs(true - apparent) > 1e-9
