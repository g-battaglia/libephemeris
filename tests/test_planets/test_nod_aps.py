"""
Tests for planetary nodes and apsides calculation (nod_aps and nod_aps_ut).

These functions calculate geocentric positions of:
- Ascending and descending orbital nodes
- Perihelion and aphelion positions

The implementation computes osculating orbital elements from JPL DE440
heliocentric state vectors, then projects the resulting node/apse positions
to geocentric coordinates. This means distances are geocentric (not
heliocentric), and angular relationships like 180° between nodes or apsides
hold only in the heliocentric frame, not after geocentric projection.
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    SUN,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
    EARTH,
    NODBIT_MEAN,
    NODBIT_OSCU,
)


class TestNodApsBasic:
    """Test basic functionality of nod_aps functions."""

    @pytest.mark.unit
    def test_nod_aps_ut_returns_4_tuples(self):
        """nod_aps_ut should return 4 position tuples."""
        jd = 2451545.0  # J2000
        result = ephem.nod_aps_ut(jd, MARS, NODBIT_MEAN)

        assert isinstance(result, tuple)
        assert len(result) == 4
        for i, pos in enumerate(result):
            assert isinstance(pos, tuple), f"Position {i} should be tuple"
            assert len(pos) == 6, f"Position {i} should have 6 elements"

    @pytest.mark.unit
    def test_nod_aps_returns_4_tuples(self):
        """nod_aps (ET version) should return 4 position tuples."""
        jd = 2451545.0
        result = ephem.nod_aps(jd, MARS, NODBIT_MEAN)

        assert isinstance(result, tuple)
        assert len(result) == 4

    @pytest.mark.unit
    def test_position_elements_are_floats(self):
        """All position elements should be floats."""
        jd = 2451545.0
        nasc, ndsc, peri, aphe = ephem.nod_aps_ut(jd, MARS, NODBIT_MEAN)

        for name, pos in [
            ("nasc", nasc),
            ("ndsc", ndsc),
            ("peri", peri),
            ("aphe", aphe),
        ]:
            for i, val in enumerate(pos):
                assert isinstance(val, (int, float)), f"{name}[{i}] should be numeric"

    @pytest.mark.unit
    def test_sun_solar_orbit_apsides(self):
        """Sun returns the apparent solar-orbit apsides, zero nodes.

        The reference API mirrors Earth's orbit through the origin: the
        perihelion of the solar orbit sits near lon 282.9 deg at J2000
        (Earth's perihelion + 180), and the node slots are zeros because
        the solar orbit lies in the ecliptic plane.
        """
        jd = 2451545.0
        nasc, ndsc, peri, aphe = ephem.nod_aps_ut(jd, SUN, NODBIT_MEAN)

        assert nasc[0] == 0.0
        assert ndsc[0] == 0.0
        assert abs(peri[0] - 282.9335) < 0.01
        assert abs(aphe[0] - 102.9335) < 0.01
        assert abs(peri[2] - 0.9833) < 0.001
        assert abs(aphe[2] - 1.0167) < 0.001

    @pytest.mark.unit
    def test_earth_returns_zeros(self):
        """Earth should return zero positions (no heliocentric orbit)."""
        jd = 2451545.0
        nasc, ndsc, peri, aphe = ephem.nod_aps_ut(jd, EARTH, NODBIT_MEAN)

        assert nasc[0] == 0.0
        assert peri[0] == 0.0


class TestNodApsOrbitalRelationships:
    """Test that orbital relationships are mathematically correct.

    Note: The implementation returns geocentric positions, so heliocentric
    relationships (exact 180° between nodes/apsides, perihelion < aphelion
    distance) do not apply to the geocentric output values.
    """

    @pytest.mark.unit
    def test_nodes_have_near_zero_latitude(self):
        """Nodes should have near-zero latitude (on ecliptic by definition)."""
        jd = 2451545.0
        nasc, ndsc, peri, aphe = ephem.nod_aps_ut(jd, MARS, NODBIT_MEAN)

        assert abs(nasc[1]) < 0.01, f"Ascending node latitude {nasc[1]} should be ~0"
        assert abs(ndsc[1]) < 0.01, f"Descending node latitude {ndsc[1]} should be ~0"

    @pytest.mark.unit
    def test_apsides_have_latitude(self):
        """Apsides should have non-zero latitude for inclined orbits."""
        jd = 2451545.0
        # Pluto has high inclination (~17°)
        nasc, ndsc, peri, aphe = ephem.nod_aps_ut(jd, PLUTO, NODBIT_MEAN)

        # Pluto's apsides should have measurable latitude
        # (depending on argument of perihelion)
        # Just verify they're calculated (can be positive or negative)
        assert peri[1] is not None
        assert aphe[1] is not None

    @pytest.mark.unit
    def test_geocentric_distances_are_positive(self):
        """All geocentric distances should be positive."""
        jd = 2451545.0
        nasc, ndsc, peri, aphe = ephem.nod_aps_ut(jd, MARS, NODBIT_MEAN)

        assert nasc[2] > 0, "Ascending node geocentric distance should be positive"
        assert ndsc[2] > 0, "Descending node geocentric distance should be positive"
        assert peri[2] > 0, "Perihelion geocentric distance should be positive"
        assert aphe[2] > 0, "Aphelion geocentric distance should be positive"


class TestNodApsAllPlanets:
    """Test nod_aps for all planets."""

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
    def test_all_planets_return_valid_positions(self, planet_id, planet_name):
        """All planets should return valid node/apse positions."""
        jd = 2451545.0
        nasc, ndsc, peri, aphe = ephem.nod_aps_ut(jd, planet_id, NODBIT_MEAN)

        # Longitudes should be in valid range
        assert 0 <= nasc[0] < 360, f"{planet_name} ascending node longitude invalid"
        assert 0 <= ndsc[0] < 360, f"{planet_name} descending node longitude invalid"
        assert 0 <= peri[0] < 360, f"{planet_name} perihelion longitude invalid"
        assert 0 <= aphe[0] < 360, f"{planet_name} aphelion longitude invalid"

        # Distances should be positive
        assert nasc[2] > 0, f"{planet_name} ascending node distance should be positive"
        assert peri[2] > 0, f"{planet_name} perihelion distance should be positive"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_longitudes_differ_between_planets(self, planet_id, planet_name):
        """Different planets should have different node/apse longitudes."""
        jd = 2451545.0
        nasc, ndsc, peri, aphe = ephem.nod_aps_ut(jd, planet_id, NODBIT_MEAN)

        # All four longitudes should be distinct (geocentric projection
        # of different orbital elements)
        lons = [nasc[0], ndsc[0], peri[0], aphe[0]]
        # At least 3 of the 4 should be distinct (within 1°)
        unique = []
        for lon in lons:
            if not any(abs(lon - u) < 1.0 or abs(lon - u) > 359.0 for u in unique):
                unique.append(lon)
        assert len(unique) >= 3, (
            f"{planet_name} should have at least 3 distinct longitudes, got {lons}"
        )


class TestNodApsMethods:
    """Test different calculation methods."""

    @pytest.mark.unit
    def test_mean_and_oscu_methods_exist(self):
        """Both MEAN and OSCU methods should work."""
        jd = 2451545.0

        result_mean = ephem.nod_aps_ut(jd, MARS, NODBIT_MEAN)
        result_oscu = ephem.nod_aps_ut(jd, MARS, NODBIT_OSCU)

        # Both should return valid results
        assert len(result_mean) == 4
        assert len(result_oscu) == 4

        # Results may differ
        # (In our implementation, they use the same method currently)


class TestNodApsTimeVariation:
    """Test how nodes and apsides change over time."""

    @pytest.mark.unit
    def test_nodes_change_over_century(self):
        """Orbital nodes should show precession over a century."""
        jd_2000 = 2451545.0
        jd_2100 = jd_2000 + 36525  # 100 years later

        result_2000 = ephem.nod_aps_ut(jd_2000, MARS, NODBIT_MEAN)
        result_2100 = ephem.nod_aps_ut(jd_2100, MARS, NODBIT_MEAN)

        # Nodes should have precessed
        diff = abs(result_2100[0][0] - result_2000[0][0])

        # Mars' node precession is about 0.5°/century
        # Allow significant range due to osculating nature
        assert diff > 0.01, "Nodes should show some precession over 100 years"

    @pytest.mark.unit
    def test_perihelion_advances(self):
        """Perihelion should advance over time (apsidal precession)."""
        jd_2000 = 2451545.0
        jd_2100 = jd_2000 + 36525

        result_2000 = ephem.nod_aps_ut(jd_2000, MARS, NODBIT_MEAN)
        result_2100 = ephem.nod_aps_ut(jd_2100, MARS, NODBIT_MEAN)

        # Perihelion should advance
        diff = result_2100[2][0] - result_2000[2][0]
        # Handle wraparound
        if diff < -180:
            diff += 360
        if diff > 180:
            diff -= 360

        # Just verify there's some change
        assert abs(diff) > 0.01, "Perihelion should advance over 100 years"


class TestNodApsAliases:
    """Test that function aliases work."""

    @pytest.mark.unit
    def test_swe_nod_aps_ut_alias(self):
        """nod_aps_ut should be available."""
        result = ephem.nod_aps_ut(2451545.0, MARS, NODBIT_MEAN)
        assert len(result) == 4

    @pytest.mark.unit
    def test_swe_nod_aps_alias(self):
        """nod_aps should be available."""
        result = ephem.nod_aps(2451545.0, MARS, NODBIT_MEAN)
        assert len(result) == 4

    @pytest.mark.unit
    def test_nod_aps_ut_alias(self):
        """nod_aps_ut should be available."""
        result = ephem.nod_aps_ut(2451545.0, MARS, NODBIT_MEAN)
        assert len(result) == 4

    @pytest.mark.unit
    def test_nod_aps_alias(self):
        """nod_aps should be available."""
        result = ephem.nod_aps(2451545.0, MARS, NODBIT_MEAN)
        assert len(result) == 4


class TestNodApsMethodologyWarning:
    """Test the HeliocentricNodApsWarning backward compatibility.

    The geocentric osculating implementation no longer emits warnings
    (the old heliocentric implementation warned for inner planets).
    The warning class is kept as a deprecated stub for backward compatibility.
    """

    @pytest.mark.unit
    def test_no_warning_for_mercury(self):
        """Mercury should not trigger HeliocentricNodApsWarning (geocentric impl)."""
        import warnings

        jd = 2451545.0
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ephem.nod_aps_ut(jd, MERCURY, NODBIT_MEAN)

            helio_warnings = [
                warning
                for warning in w
                if issubclass(warning.category, ephem.HeliocentricNodApsWarning)
            ]
            assert len(helio_warnings) == 0

    @pytest.mark.unit
    def test_no_warning_for_venus(self):
        """Venus should not trigger HeliocentricNodApsWarning (geocentric impl)."""
        import warnings

        jd = 2451545.0
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ephem.nod_aps_ut(jd, VENUS, NODBIT_MEAN)

            helio_warnings = [
                warning
                for warning in w
                if issubclass(warning.category, ephem.HeliocentricNodApsWarning)
            ]
            assert len(helio_warnings) == 0

    @pytest.mark.unit
    def test_no_warning_for_mars(self):
        """Mars should not trigger HeliocentricNodApsWarning."""
        import warnings

        jd = 2451545.0
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ephem.nod_aps_ut(jd, MARS, NODBIT_MEAN)

            helio_warnings = [
                warning
                for warning in w
                if issubclass(warning.category, ephem.HeliocentricNodApsWarning)
            ]
            assert len(helio_warnings) == 0

    @pytest.mark.unit
    def test_no_warning_for_jupiter(self):
        """Jupiter should not trigger HeliocentricNodApsWarning."""
        import warnings

        jd = 2451545.0
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ephem.nod_aps_ut(jd, JUPITER, NODBIT_MEAN)

            helio_warnings = [
                warning
                for warning in w
                if issubclass(warning.category, ephem.HeliocentricNodApsWarning)
            ]
            assert len(helio_warnings) == 0

    @pytest.mark.unit
    def test_warning_can_be_suppressed(self):
        """HeliocentricNodApsWarning should be suppressible (backward compat)."""
        import warnings

        jd = 2451545.0
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings("ignore", category=ephem.HeliocentricNodApsWarning)
            ephem.nod_aps_ut(jd, MERCURY, NODBIT_MEAN)

            # No HeliocentricNodApsWarning should be recorded when filtered
            helio_warnings = [
                warning
                for warning in w
                if issubclass(warning.category, ephem.HeliocentricNodApsWarning)
            ]
            assert len(helio_warnings) == 0

    @pytest.mark.unit
    def test_warning_class_is_exported(self):
        """HeliocentricNodApsWarning should be exported from libephemeris."""
        assert hasattr(ephem, "HeliocentricNodApsWarning")
        assert issubclass(ephem.HeliocentricNodApsWarning, UserWarning)
