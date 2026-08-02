"""
Tests for orbital elements calculation (get_orbital_elements).

This function calculates Keplerian orbital elements:
- Semi-major axis (a)
- Eccentricity (e)
- Inclination (i)
- Longitude of ascending node (Omega)
- Argument of perihelion (omega)
- Mean anomaly (M)
- And other derived parameters

The function returns a 50-element tuple matching the reference ephemeris's format.
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
)


class TestOrbitalElementsBasic:
    """Test basic functionality of get_orbital_elements."""

    @pytest.mark.unit
    def test_returns_50_elements(self):
        """get_orbital_elements should return a tuple of 50 elements."""
        jd = 2451545.0  # J2000
        elements = ephem.get_orbital_elements(jd, MARS, 0)

        assert isinstance(elements, tuple)
        assert len(elements) == 50
        assert isinstance(elements[0], float)

    @pytest.mark.unit
    def test_returns_50_elements_ut(self):
        """get_orbital_elements_ut should return a tuple of 50 elements."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements_ut(jd, MARS, 0)

        assert isinstance(elements, tuple)
        assert len(elements) == 50

    @pytest.mark.unit
    def test_all_elements_are_floats(self):
        """All elements should be floats."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, JUPITER, 0)

        for i, val in enumerate(elements):
            assert isinstance(val, (int, float)), f"Element {i} should be numeric"

    @pytest.mark.unit
    def test_sun_raises_not_valid(self):
        """Sun has no heliocentric orbit: get_orbital_elements raises, matching
        the measured reference (``object 0 not valid``) rather than returning a
        silent all-zeros tuple."""
        from libephemeris.exceptions import Error

        jd = 2451545.0
        with pytest.raises(Error) as exc:
            ephem.get_orbital_elements(jd, SUN, 0)
        assert (
            str(exc.value) == "get_orbital_elements: error in get_orbital_elements(): "
            "object 0 not valid"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("ipl", [10, 11, 12, 13, 21, 22, -1, -2])
    def test_points_without_orbit_raise_not_valid(self, ipl):
        """Lunar nodes/apogees (10-13), interpolated apsides (21, 22) and
        negatives raise ``object N not valid`` for both time scales."""
        from libephemeris.exceptions import Error

        jd = 2451545.0
        expected = (
            f"get_orbital_elements: error in get_orbital_elements(): "
            f"object {ipl} not valid"
        )
        with pytest.raises(Error) as exc:
            ephem.get_orbital_elements(jd, ipl, 0)
        assert str(exc.value) == expected
        expected_ut = expected.replace(
            "get_orbital_elements:", "get_orbital_elements_ut:", 1
        )
        with pytest.raises(Error) as exc_ut:
            ephem.get_orbital_elements_ut(jd, ipl, 0)
        assert str(exc_ut.value) == expected_ut

    @pytest.mark.unit
    @pytest.mark.parametrize("ipl", [23, 30, 39])
    def test_undefined_block_raises_illegal_planet(self, ipl):
        """Ids 23-39 are an undefined block: ``illegal planet number N``."""
        from libephemeris.exceptions import Error

        jd = 2451545.0
        with pytest.raises(Error) as exc:
            ephem.get_orbital_elements(jd, ipl, 0)
        assert str(exc.value) == f"get_orbital_elements: illegal planet number {ipl}."

    @pytest.mark.unit
    @pytest.mark.parametrize("ipl", [59, 69, 99, 200, 9001, 9999])
    def test_sourceless_block_raises_not_valid(self, ipl):
        """Ids between the last fictitious body (58) and AST_OFFSET have no
        element source: the measured reference raises for the whole block,
        and the library's own position pipeline rejects them too. A missing
        guard used to leak a zero-initialized 50-tuple here (a "real orbit"
        at 0 AU)."""
        from libephemeris.exceptions import Error

        jd = 2451545.0
        with pytest.raises(Error) as exc:
            ephem.get_orbital_elements(jd, ipl, 0)
        assert f"object {ipl} not valid" in str(exc.value)

    @pytest.mark.unit
    def test_asteroid_zero_raises_typed_error(self):
        """AST_OFFSET itself ("asteroid 0") propagates the position
        pipeline's typed error (a subclass of Error) instead of returning
        zeros; real numbered asteroids keep their elements."""
        from libephemeris.exceptions import Error

        jd = 2451545.0
        with pytest.raises(Error):
            ephem.get_orbital_elements(jd, 10000, 0)
        elems = ephem.get_orbital_elements(jd, 10001, 0)
        assert elems[0] > 1.0  # Ceres semi-major axis ~2.77 AU

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "ipl,a_expected,e_expected",
        [
            (40, 40.99837, 0.00460),  # Cupido
            (44, 70.29949, 0.0),  # Apollon
            (47, 83.66907, 0.0),  # Poseidon
        ],
    )
    def test_uranian_elements_from_published_neely_rows(
        self, ipl, a_expected, e_expected
    ):
        """The Uranians (40-47) are propagated from Neely's published element
        rows (Hamburg school); the pinned values are this library's own
        regression output for those published elements."""
        jd = 2451545.0
        el = ephem.get_orbital_elements(jd, ipl, 0)
        assert el[0] == pytest.approx(a_expected, abs=1e-3)
        assert el[1] == pytest.approx(e_expected, abs=1e-3)

    @pytest.mark.unit
    def test_white_moon_is_geocentric_model(self):
        """White Moon (56) is a geocentric circular construction here: its
        osculating semi-major axis is Earth-relative (~0.053 AU), an intentional
        divergence from the reference's ~0.909 AU heliocentric element set."""
        jd = 2451545.0
        el = ephem.get_orbital_elements(jd, 56, 0)
        assert 0.04 < el[0] < 0.07  # geocentric radius, not the reference's 0.909
        assert el[1] < 0.01  # near-circular


class TestOrbitalElementsValues:
    """Test that orbital elements match known values."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name,expected_a,tolerance",
        [
            (MERCURY, "Mercury", 0.387, 0.01),
            (VENUS, "Venus", 0.723, 0.01),
            (MARS, "Mars", 1.524, 0.01),
            (JUPITER, "Jupiter", 5.203, 0.05),
            (SATURN, "Saturn", 9.537, 0.1),
            (URANUS, "Uranus", 19.19, 0.2),
            (NEPTUNE, "Neptune", 30.07, 0.3),
            (PLUTO, "Pluto", 39.48, 0.5),
        ],
    )
    def test_semi_major_axis(self, planet_id, planet_name, expected_a, tolerance):
        """Semi-major axis should match known values."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, planet_id, 0)

        a = elements[0]
        assert abs(a - expected_a) < tolerance, (
            f"{planet_name} semi-major axis {a:.4f} differs from expected {expected_a}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name,expected_e,tolerance",
        [
            (MERCURY, "Mercury", 0.2056, 0.01),
            (VENUS, "Venus", 0.0068, 0.005),
            (MARS, "Mars", 0.0934, 0.01),
            (JUPITER, "Jupiter", 0.0484, 0.01),
            (SATURN, "Saturn", 0.0539, 0.01),
            (PLUTO, "Pluto", 0.2488, 0.02),
        ],
    )
    def test_eccentricity(self, planet_id, planet_name, expected_e, tolerance):
        """Eccentricity should match known values."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, planet_id, 0)

        e = elements[1]
        assert abs(e - expected_e) < tolerance, (
            f"{planet_name} eccentricity {e:.4f} differs from expected {expected_e}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name,expected_i,tolerance",
        [
            (MERCURY, "Mercury", 7.0, 0.5),
            (VENUS, "Venus", 3.4, 0.5),
            (MARS, "Mars", 1.85, 0.2),
            (JUPITER, "Jupiter", 1.3, 0.2),
            (SATURN, "Saturn", 2.5, 0.3),
            (PLUTO, "Pluto", 17.1, 1.0),
        ],
    )
    def test_inclination(self, planet_id, planet_name, expected_i, tolerance):
        """Inclination should match known values."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, planet_id, 0)

        i = elements[2]
        assert abs(i - expected_i) < tolerance, (
            f"{planet_name} inclination {i:.4f}° differs from expected {expected_i}°"
        )


class TestOrbitalElementsRelationships:
    """Test mathematical relationships between elements."""

    @pytest.mark.unit
    def test_perihelion_from_a_and_e(self):
        """Perihelion distance should equal a*(1-e)."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, MARS, 0)

        a = elements[0]
        e = elements[1]
        q = elements[15]  # Perihelion distance

        expected_q = a * (1 - e)
        assert abs(q - expected_q) < 0.001, (
            f"Perihelion {q:.4f} should equal a*(1-e) = {expected_q:.4f}"
        )

    @pytest.mark.unit
    def test_aphelion_from_a_and_e(self):
        """Aphelion distance should equal a*(1+e)."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, MARS, 0)

        a = elements[0]
        e = elements[1]
        Q = elements[16]  # Aphelion distance

        expected_Q = a * (1 + e)
        assert abs(Q - expected_Q) < 0.001, (
            f"Aphelion {Q:.4f} should equal a*(1+e) = {expected_Q:.4f}"
        )

    @pytest.mark.unit
    def test_perihelion_less_than_aphelion(self):
        """Perihelion should always be less than aphelion."""
        jd = 2451545.0

        for planet_id in [MERCURY, MARS, JUPITER, SATURN]:
            elements = ephem.get_orbital_elements(jd, planet_id, 0)
            q = elements[15]
            Q = elements[16]
            assert q < Q, f"Perihelion {q} should be less than aphelion {Q}"

    @pytest.mark.unit
    def test_longitude_of_perihelion(self):
        """Longitude of perihelion should equal Omega + omega."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, MARS, 0)

        Omega = elements[3]  # Longitude of ascending node
        omega = elements[4]  # Argument of perihelion
        varpi = elements[5]  # Longitude of perihelion

        expected_varpi = (Omega + omega) % 360.0
        diff = abs(varpi - expected_varpi)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, (
            f"varpi {varpi:.4f} should equal Omega+omega = {expected_varpi:.4f}"
        )

    @pytest.mark.unit
    def test_mean_longitude(self):
        """Mean longitude should equal Omega + omega + M."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, JUPITER, 0)

        Omega = elements[3]
        omega = elements[4]
        M = elements[6]  # Mean anomaly
        L = elements[9]  # Mean longitude

        expected_L = (Omega + omega + M) % 360.0
        diff = abs(L - expected_L)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"L {L:.4f} should equal Omega+omega+M = {expected_L:.4f}"


class TestOrbitalElementsAnglesRange:
    """Test that angles are in valid ranges."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id",
        [MERCURY, VENUS, MARS, JUPITER, SATURN],
    )
    def test_angles_in_valid_range(self, planet_id):
        """All angular elements should be in 0-360 range."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, planet_id, 0)

        # Indices of angular elements
        angle_indices = [2, 3, 4, 5, 6, 7, 8, 9]  # i, Omega, omega, M, nu, E, L, varpi
        angle_names = ["i", "Omega", "omega", "M", "nu", "E", "L", "varpi"]

        for idx, name in zip(angle_indices, angle_names):
            val = elements[idx]
            assert 0 <= val < 360, f"{name} ({val}) should be in [0, 360)"


class TestOrbitalElementsOrbitalPeriod:
    """Test orbital period calculations."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name,expected_period_years,tolerance",
        [
            (MERCURY, "Mercury", 0.241, 0.01),
            (VENUS, "Venus", 0.615, 0.02),
            (MARS, "Mars", 1.881, 0.05),
            (JUPITER, "Jupiter", 11.86, 0.2),
            (SATURN, "Saturn", 29.46, 0.5),
        ],
    )
    def test_orbital_period(
        self, planet_id, planet_name, expected_period_years, tolerance
    ):
        """Orbital period should match known values."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, planet_id, 0)

        P = elements[10]  # Sidereal orbital period in tropical years
        assert abs(P - expected_period_years) < tolerance, (
            f"{planet_name} period {P:.3f} years differs from "
            f"expected {expected_period_years}"
        )

    @pytest.mark.unit
    def test_mean_motion_from_period(self):
        """Mean daily motion should be consistent with orbital period."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, MARS, 0)

        n = elements[11]  # Mean daily motion (degrees/day)
        P = elements[10]  # Sidereal orbital period (years)

        # n should equal 360 / (P * 365.24219)
        expected_n = 360.0 / (P * 365.24219)
        assert abs(n - expected_n) < 0.001, (
            f"Mean motion {n:.6f} deg/day differs from expected {expected_n:.6f}"
        )


class TestOrbitalElementsMoon:
    """Test orbital elements for the Moon (geocentric orbit)."""

    @pytest.mark.unit
    def test_moon_has_geocentric_elements(self):
        """Moon should have valid geocentric orbital elements."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, MOON, 0)

        a = elements[0]  # Semi-major axis
        e = elements[1]  # Eccentricity
        i = elements[2]  # Inclination

        # Moon's semi-major axis is about 384,400 km = 0.00257 AU
        assert 0.002 < a < 0.003, f"Moon semi-major axis {a} should be ~0.00257 AU"

        # Moon's eccentricity is about 0.0549
        assert 0.04 < e < 0.07, f"Moon eccentricity {e} should be ~0.0549"

        # Moon's inclination to ecliptic is about 5.145°
        assert 4.0 < i < 6.0, f"Moon inclination {i} should be ~5.145°"

    @pytest.mark.unit
    def test_moon_orbital_period(self):
        """Moon's orbital period should be about 27.3 days."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, MOON, 0)

        P_years = elements[10]  # Sidereal orbital period
        P_days = P_years * 365.24219

        # Sidereal month is about 27.32 days
        assert 25 < P_days < 30, (
            f"Moon orbital period {P_days:.2f} days should be ~27.3"
        )


class TestOrbitalElementsCurrentDistance:
    """Test current heliocentric distance."""

    @pytest.mark.unit
    def test_current_distance_positive(self):
        """Current heliocentric distance should be positive."""
        jd = 2451545.0

        for planet_id in [MARS, JUPITER, SATURN]:
            elements = ephem.get_orbital_elements(jd, planet_id, 0)
            r = elements[16]  # Current distance
            assert r > 0, f"Current distance should be positive, got {r}"

    @pytest.mark.unit
    def test_current_distance_between_q_and_Q(self):
        """Current distance should be between perihelion and aphelion."""
        jd = 2451545.0

        for planet_id in [MARS, JUPITER, SATURN]:
            elements = ephem.get_orbital_elements(jd, planet_id, 0)
            q = elements[11]
            Q = elements[12]
            r = elements[16]

            assert q <= r <= Q, (
                f"Distance {r} should be between perihelion {q} and aphelion {Q}"
            )


class TestOrbitalElementsAliases:
    """Test that function aliases work correctly."""

    @pytest.mark.unit
    def test_swe_get_orbital_elements_alias(self):
        """get_orbital_elements should be available."""
        result = ephem.get_orbital_elements(2451545.0, MARS, 0)
        assert len(result) == 50
        assert isinstance(result[0], float)

    @pytest.mark.unit
    def test_swe_get_orbital_elements_ut_alias(self):
        """get_orbital_elements_ut should be available."""
        result = ephem.get_orbital_elements_ut(2451545.0, MARS, 0)
        assert len(result) == 50
        assert isinstance(result[0], float)

    @pytest.mark.unit
    def test_get_orbital_elements_alias(self):
        """get_orbital_elements should be available."""
        result = ephem.get_orbital_elements(2451545.0, MARS, 0)
        assert len(result) == 50
        assert isinstance(result[0], float)

    @pytest.mark.unit
    def test_get_orbital_elements_ut_alias(self):
        """get_orbital_elements_ut should be available."""
        result = ephem.get_orbital_elements_ut(2451545.0, MARS, 0)
        assert len(result) == 50
        assert isinstance(result[0], float)


class TestOrbitalElementsAstronomicalAlmanac:
    """FLG_ORBEL_AA (Astronomical Almanac central mass) osculating elements.

    With FLG_ORBEL_AA the two-body fit is reduced about the Sun plus every major
    planet at or interior to the body's orbit (Earth = Earth+Moon barycentre).
    The semi-major axis shrinks, the perihelion-referred angles shift solidly,
    and the mean longitude and the orbit-plane elements (i, node) stay put.
    These properties are backend-independent and need no external reference.
    """

    JD = 2451545.0  # J2000

    @pytest.mark.unit
    def test_flag_is_exported(self):
        assert ephem.FLG_ORBEL_AA == 32768

    @pytest.mark.unit
    def test_mercury_aa_equals_default(self):
        """Mercury has no interior planet, so its Almanac central mass equals
        the default Sun+Mercury mass: every element is unchanged."""
        default = ephem.get_orbital_elements(self.JD, MERCURY, 0)
        almanac = ephem.get_orbital_elements(self.JD, MERCURY, ephem.FLG_ORBEL_AA)
        for i in range(17):
            assert almanac[i] == pytest.approx(default[i], abs=1e-12, rel=1e-12)

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id", [MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO]
    )
    def test_outer_planet_semimajor_axis_shrinks(self, planet_id):
        """Adding interior planetary mass to the central body raises the
        effective GM, so the osculating semi-major axis shrinks."""
        a_default = ephem.get_orbital_elements(self.JD, planet_id, 0)[0]
        a_almanac = ephem.get_orbital_elements(self.JD, planet_id, ephem.FLG_ORBEL_AA)[
            0
        ]
        assert a_almanac < a_default

    @pytest.mark.unit
    def test_saturn_jump_dominated_by_jupiter_mass(self):
        """The Almanac shrink jumps sharply from Jupiter to Saturn because
        Jupiter's own (large) mass enters the interior sum at Saturn: Saturn's
        fractional shrink is far larger than the terrestrial-only shrink at
        Mars or the still-interior-only shrink at Jupiter."""

        def frac_shrink(pid):
            a_d = ephem.get_orbital_elements(self.JD, pid, 0)[0]
            a_a = ephem.get_orbital_elements(self.JD, pid, ephem.FLG_ORBEL_AA)[0]
            return (a_d - a_a) / a_d

        mars = frac_shrink(MARS)
        jupiter = frac_shrink(JUPITER)
        saturn = frac_shrink(SATURN)
        # Jupiter's mass ratio (~1/1047) dwarfs the terrestrial interior sum
        # (~6e-6), so the Saturn jump is orders of magnitude larger.
        assert saturn > 20.0 * mars
        assert saturn > 20.0 * jupiter

    @pytest.mark.unit
    @pytest.mark.parametrize("planet_id", [MARS, JUPITER, SATURN, NEPTUNE, PLUTO])
    def test_mean_longitude_invariant(self, planet_id):
        """The perihelion-referred angles shift solidly (omega up, M down), so
        the mean longitude L = node + omega + M is left essentially invariant."""
        default = ephem.get_orbital_elements(self.JD, planet_id, 0)
        almanac = ephem.get_orbital_elements(self.JD, planet_id, ephem.FLG_ORBEL_AA)
        d_L = abs(((almanac[9] - default[9]) + 180.0) % 360.0 - 180.0)
        assert d_L < 0.01

    @pytest.mark.unit
    @pytest.mark.parametrize("planet_id", [MARS, JUPITER, SATURN, NEPTUNE, PLUTO])
    def test_orbit_plane_elements_unchanged(self, planet_id):
        """Inclination and node depend only on the angular-momentum direction,
        not on the central GM, so they are untouched by FLG_ORBEL_AA."""
        default = ephem.get_orbital_elements(self.JD, planet_id, 0)
        almanac = ephem.get_orbital_elements(self.JD, planet_id, ephem.FLG_ORBEL_AA)
        assert almanac[2] == pytest.approx(default[2], abs=1e-9)  # inclination
        d_node = abs(((almanac[3] - default[3]) + 180.0) % 360.0 - 180.0)
        assert d_node < 1e-6

    @pytest.mark.unit
    def test_perihelion_angle_shifts_solidly(self):
        """omega (argument of perihelion) and M (mean anomaly) shift by equal
        and opposite amounts under the Almanac central mass."""
        default = ephem.get_orbital_elements(self.JD, NEPTUNE, 0)
        almanac = ephem.get_orbital_elements(self.JD, NEPTUNE, ephem.FLG_ORBEL_AA)
        d_omega = ((almanac[4] - default[4]) + 180.0) % 360.0 - 180.0
        d_M = ((almanac[6] - default[6]) + 180.0) % 360.0 - 180.0
        assert abs(d_omega) > 1.0  # a real, several-degree shift
        assert d_omega == pytest.approx(-d_M, rel=0.02)

    @pytest.mark.unit
    def test_asteroid_semimajor_axis_shrinks(self):
        """A main-belt asteroid (interior planets: Mercury..Mars) also shrinks
        under the Almanac central mass, via the minor-body element path."""
        from libephemeris.constants import CERES

        a_default = ephem.get_orbital_elements(self.JD, CERES, 0)[0]
        a_almanac = ephem.get_orbital_elements(self.JD, CERES, ephem.FLG_ORBEL_AA)[0]
        assert a_almanac < a_default


class TestPeriodSlotConventions:
    """Slots [12]/[13] follow the published period definitions.

    The tropical-period slot is 360/(n + p) with p the IAU general
    precession in longitude, expressed in tropical years (Explanatory
    Supplement, 3rd ed., glossary); the synodic slot uses the textbook
    sidereal-period relation 1/P_syn = 1/P_E - 1/P (e.g. Murray &
    Dermott 1999) with Earth's sidereal year. The constant ~3.9e-5 unit factor
    an external implementation adds to slot [12] has no published basis
    and is a documented intentional divergence.
    """

    @pytest.mark.unit
    def test_tropical_slot_is_published_definition(self):
        el = ephem.get_orbital_elements(2451545.0, 5, 0)
        assert el[12] == pytest.approx(11.86743, abs=2e-4)
        assert el[13] == pytest.approx(398.851, abs=0.02)

    @pytest.mark.unit
    def test_deg_midp_stays_in_range(self):
        import libephemeris as le

        assert le.deg_midp(0.1, 359.9) == 0.0
        assert le.deg_midp(0.6, 359.4) == 0.0
        assert 0.0 <= le.deg_midp(2.0, 358.0) < 360.0
