"""
Tests for the orbital elements file parser.

This module tests the parse_orbital_elements() function and related utilities
for parsing hypothetical/fictitious body orbital elements files, as well
as the bundled fictitious_orbits.csv dataset included with libephemeris.
"""

import math
import tempfile
from pathlib import Path

import pytest

from libephemeris.hypothetical import (
    # Canonical names (SE-independent)
    TPolynomial,
    get_bundled_fictitious_orbits_path,
    load_bundled_fictitious_orbits,
    _parse_t_polynomial,
    _parse_epoch_or_equinox,
    parse_orbital_elements,
    OrbitalElements,
    get_orbital_body_by_name,
    calc_orbital_position,
    get_bundled_seorbel_path,
    load_bundled_seorbel,
    _parse_orbital_elements_line,
)


class TestTPolynomial:
    """Tests for the TPolynomial dataclass."""

    def test_constant_only(self):
        """Test polynomial with only constant term."""
        poly = TPolynomial(constant=42.5, linear=0.0)
        assert poly.evaluate(0.0) == 42.5
        assert poly.evaluate(1.0) == 42.5
        assert poly.evaluate(-1.0) == 42.5

    def test_with_linear_term(self):
        """Test polynomial with linear term."""
        poly = TPolynomial(constant=100.0, linear=50.0)
        assert poly.evaluate(0.0) == 100.0
        assert poly.evaluate(1.0) == 150.0
        assert poly.evaluate(2.0) == 200.0
        assert poly.evaluate(-1.0) == 50.0

    def test_negative_linear(self):
        """Test polynomial with negative linear term."""
        poly = TPolynomial(constant=100.0, linear=-10.0)
        assert poly.evaluate(0.0) == 100.0
        assert poly.evaluate(1.0) == 90.0
        assert poly.evaluate(5.0) == 50.0


class TestParseEpochOrEquinox:
    """Tests for _parse_epoch_or_equinox helper."""

    def test_j1900(self):
        """Test parsing J1900 epoch."""
        jd, is_jdate = _parse_epoch_or_equinox("J1900")
        assert jd == 2415020.0
        assert is_jdate is False

    def test_j2000(self):
        """Test parsing J2000 epoch."""
        jd, is_jdate = _parse_epoch_or_equinox("J2000")
        assert jd == 2451545.0
        assert is_jdate is False

    def test_b1950(self):
        """Test parsing B1950 epoch."""
        jd, is_jdate = _parse_epoch_or_equinox("B1950")
        assert jd == pytest.approx(2433282.42345905, abs=0.0001)
        assert is_jdate is False

    def test_jdate(self):
        """Test parsing JDATE equinox."""
        jd, is_jdate = _parse_epoch_or_equinox("JDATE")
        assert jd is None
        assert is_jdate is True

    def test_numeric_jd(self):
        """Test parsing numeric Julian Day."""
        jd, is_jdate = _parse_epoch_or_equinox("2451545.0")
        assert jd == 2451545.0
        assert is_jdate is False

    def test_case_insensitive(self):
        """Test that parsing is case-insensitive."""
        jd1, _ = _parse_epoch_or_equinox("j2000")
        jd2, _ = _parse_epoch_or_equinox("J2000")
        assert jd1 == jd2


class TestParseTPolynomial:
    """Tests for _parse_t_polynomial helper."""

    def test_simple_number(self):
        """Test parsing a simple number."""
        poly = _parse_t_polynomial("163.7409")
        assert poly.constant == pytest.approx(163.7409)
        assert poly.linear == 0.0

    def test_negative_number(self):
        """Test parsing a negative number."""
        poly = _parse_t_polynomial("-45.5")
        assert poly.constant == pytest.approx(-45.5)
        assert poly.linear == 0.0

    def test_polynomial_plus(self):
        """Test parsing polynomial with + operator."""
        poly = _parse_t_polynomial("252.8987988 + 707550.7341 * T")
        assert poly.constant == pytest.approx(252.8987988)
        assert poly.linear == pytest.approx(707550.7341)

    def test_polynomial_minus(self):
        """Test parsing polynomial with - operator."""
        poly = _parse_t_polynomial("47.787931 - 1670.056 * T")
        assert poly.constant == pytest.approx(47.787931)
        assert poly.linear == pytest.approx(-1670.056)

    def test_polynomial_compact(self):
        """Test parsing polynomial without spaces."""
        poly = _parse_t_polynomial("322.212069+1670.056*T")
        assert poly.constant == pytest.approx(322.212069)
        assert poly.linear == pytest.approx(1670.056)

    def test_polynomial_compact_minus(self):
        """Test parsing compact polynomial with minus."""
        poly = _parse_t_polynomial("47.787931-1670.056*T")
        assert poly.constant == pytest.approx(47.787931)
        assert poly.linear == pytest.approx(-1670.056)

    def test_zero_constant(self):
        """Test parsing with zero constant."""
        poly = _parse_t_polynomial("0.0")
        assert poly.constant == 0.0
        assert poly.linear == 0.0


class TestParseOrbitalElements:
    """Tests for parse_orbital_elements function."""

    def test_parse_real_seorbel_file(self):
        """Test parsing the actual Swiss Ephemeris seorbel.txt file."""
        seorbel_path = (
            Path(__file__).parent.parent / "ephe" / "seorbel.txt"
        )
        if not seorbel_path.exists():
            pytest.skip("seorbel.txt not found in reference ephemeris dir")

        elements = parse_orbital_elements(seorbel_path)

        # Should have parsed multiple elements
        assert len(elements) > 0

        # Check that Cupido is first (it's line #1 in the file)
        cupido = get_orbital_body_by_name(elements, "Cupido")
        assert cupido is not None
        assert cupido.semi_axis == pytest.approx(40.99837)
        assert cupido.eccentricity.constant == pytest.approx(0.00460)
        assert cupido.epoch_jd == 2415020.0  # J1900

    def test_parse_simple_file(self):
        """Test parsing a simple custom seorbel file."""
        content = """\
# Custom test file
J2000, J2000, 0.0, 100.0, 0.1, 45.0, 30.0, 5.0, TestPlanet
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_orbital_elements(filepath)
            assert len(elements) == 1

            elem = elements[0]
            assert elem.name == "TestPlanet"
            assert elem.epoch_jd == 2451545.0  # J2000
            assert elem.semi_axis == 100.0
            assert elem.eccentricity.constant == 0.1
            assert elem.arg_perihelion.constant == 45.0
            assert elem.asc_node.constant == 30.0
            assert elem.inclination.constant == 5.0
            assert elem.mean_anomaly.constant == 0.0
        finally:
            Path(filepath).unlink()

    def test_parse_geocentric_body(self):
        """Test parsing a geocentric body (with ', geo' suffix)."""
        content = """\
J2000, JDATE, 248.8833, 0.0029833, 0.0, 0.0, 0.0, 0.0, Waldemath, geo
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_orbital_elements(filepath)
            assert len(elements) == 1

            elem = elements[0]
            assert elem.name == "Waldemath"
            assert elem.is_geocentric is True
            assert elem.equinox_is_jdate is True
        finally:
            Path(filepath).unlink()

    def test_parse_t_polynomial_elements(self):
        """Test parsing elements with T-polynomial expressions."""
        # From the actual seorbel.txt Vulcan entry
        content = """\
J1900,JDATE, 252.8987988 + 707550.7341 * T, 0.13744, 0.019, 322.212069+1670.056*T, 47.787931-1670.056*T, 7.5, Vulcan
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_orbital_elements(filepath)
            assert len(elements) == 1

            elem = elements[0]
            assert elem.name == "Vulcan"
            assert elem.epoch_jd == 2415020.0  # J1900
            assert elem.equinox_is_jdate is True

            # Check T-polynomial parsing
            assert elem.mean_anomaly.constant == pytest.approx(252.8987988)
            assert elem.mean_anomaly.linear == pytest.approx(707550.7341)

            assert elem.arg_perihelion.constant == pytest.approx(322.212069)
            assert elem.arg_perihelion.linear == pytest.approx(1670.056)

            assert elem.asc_node.constant == pytest.approx(47.787931)
            assert elem.asc_node.linear == pytest.approx(-1670.056)

            assert elem.semi_axis == pytest.approx(0.13744)
            assert elem.eccentricity.constant == pytest.approx(0.019)
            assert elem.inclination.constant == pytest.approx(7.5)
        finally:
            Path(filepath).unlink()

    def test_parse_inline_comments(self):
        """Test that inline comments are handled correctly."""
        content = """\
J1900, J1900, 163.7409, 40.99837, 0.00460, 171.4333, 129.8325, 1.0833, Cupido   # 1
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_orbital_elements(filepath)
            assert len(elements) == 1
            assert elements[0].name == "Cupido"
        finally:
            Path(filepath).unlink()

    def test_skip_comment_lines(self):
        """Test that comment-only lines are skipped."""
        content = """\
# This is a comment
    # Indented comment
J2000, J2000, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0, Planet1
# Another comment
J2000, J2000, 90.0, 60.0, 0.0, 0.0, 0.0, 0.0, Planet2
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_orbital_elements(filepath)
            assert len(elements) == 2
            assert elements[0].name == "Planet1"
            assert elements[1].name == "Planet2"
        finally:
            Path(filepath).unlink()

    def test_skip_empty_lines(self):
        """Test that empty lines are skipped."""
        content = """\

J2000, J2000, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0, Planet1

J2000, J2000, 90.0, 60.0, 0.0, 0.0, 0.0, 0.0, Planet2

"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_orbital_elements(filepath)
            assert len(elements) == 2
        finally:
            Path(filepath).unlink()

    def test_file_not_found(self):
        """Test that FileNotFoundError is raised for missing file."""
        with pytest.raises(FileNotFoundError):
            parse_orbital_elements("/nonexistent/path/seorbel.txt")

    def test_parse_numeric_epoch(self):
        """Test parsing with numeric Julian Day epoch."""
        content = """\
2368547.66, 2431456.5, 0.0, 77.775, 0.3, 0.7, 0, 0, Isis-Transpluto
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_orbital_elements(filepath)
            assert len(elements) == 1

            elem = elements[0]
            assert elem.name == "Isis-Transpluto"
            assert elem.epoch_jd == pytest.approx(2368547.66)
            assert elem.equinox_jd == pytest.approx(2431456.5)
            assert elem.semi_axis == pytest.approx(77.775)
            assert elem.eccentricity.constant == pytest.approx(0.3)
        finally:
            Path(filepath).unlink()

    def test_preserves_order(self):
        """Test that elements are returned in file order."""
        content = """\
J2000, J2000, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0, First
J2000, J2000, 0.0, 60.0, 0.0, 0.0, 0.0, 0.0, Second
J2000, J2000, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, Third
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_orbital_elements(filepath)
            assert len(elements) == 3
            assert elements[0].name == "First"
            assert elements[1].name == "Second"
            assert elements[2].name == "Third"
        finally:
            Path(filepath).unlink()


class TestGetOrbitalBodyByName:
    """Tests for get_orbital_body_by_name function."""

    def test_find_existing_body(self):
        """Test finding an existing body by name."""
        elements = [
            OrbitalElements(
                name="TestBody",
                epoch_jd=2451545.0,
                equinox_jd=2451545.0,
                equinox_is_jdate=False,
                mean_anomaly=TPolynomial(0.0, 0.0),
                semi_axis=50.0,
                eccentricity=TPolynomial(0.1, 0.0),
                arg_perihelion=TPolynomial(0.0, 0.0),
                asc_node=TPolynomial(0.0, 0.0),
                inclination=TPolynomial(0.0, 0.0),
            )
        ]
        result = get_orbital_body_by_name(elements, "TestBody")
        assert result is not None
        assert result.name == "TestBody"

    def test_case_insensitive_search(self):
        """Test that name search is case-insensitive."""
        elements = [
            OrbitalElements(
                name="TestBody",
                epoch_jd=2451545.0,
                equinox_jd=2451545.0,
                equinox_is_jdate=False,
                mean_anomaly=TPolynomial(0.0, 0.0),
                semi_axis=50.0,
                eccentricity=TPolynomial(0.1, 0.0),
                arg_perihelion=TPolynomial(0.0, 0.0),
                asc_node=TPolynomial(0.0, 0.0),
                inclination=TPolynomial(0.0, 0.0),
            )
        ]
        result = get_orbital_body_by_name(elements, "testbody")
        assert result is not None
        assert result.name == "TestBody"

    def test_not_found_returns_none(self):
        """Test that searching for non-existent body returns None."""
        elements = [
            OrbitalElements(
                name="TestBody",
                epoch_jd=2451545.0,
                equinox_jd=2451545.0,
                equinox_is_jdate=False,
                mean_anomaly=TPolynomial(0.0, 0.0),
                semi_axis=50.0,
                eccentricity=TPolynomial(0.1, 0.0),
                arg_perihelion=TPolynomial(0.0, 0.0),
                asc_node=TPolynomial(0.0, 0.0),
                inclination=TPolynomial(0.0, 0.0),
            )
        ]
        result = get_orbital_body_by_name(elements, "NonExistent")
        assert result is None


class TestCalcOrbitalPosition:
    """Tests for calc_orbital_position function."""

    def test_circular_orbit(self):
        """Test position calculation for circular orbit."""
        elem = OrbitalElements(
            name="CircularPlanet",
            epoch_jd=2451545.0,  # J2000.0
            equinox_jd=2451545.0,
            equinox_is_jdate=False,
            mean_anomaly=TPolynomial(0.0, 0.0),  # At 0 deg at epoch
            semi_axis=1.0,  # 1 AU
            eccentricity=TPolynomial(0.0, 0.0),  # Circular
            arg_perihelion=TPolynomial(0.0, 0.0),
            asc_node=TPolynomial(0.0, 0.0),
            inclination=TPolynomial(0.0, 0.0),
        )

        pos = calc_orbital_position(elem, 2451545.0)  # At epoch
        lon, lat, dist, dlon, dlat, ddist = pos

        # At epoch with M=0, omega=0, Omega=0, longitude should be ~0
        assert lon == pytest.approx(0.0, abs=0.1)
        assert lat == pytest.approx(0.0, abs=0.01)
        assert dist == pytest.approx(1.0, abs=0.001)

        # Mean motion should be ~360/(365.25) deg/day = ~0.986 deg/day
        assert dlon == pytest.approx(360.0 / 365.25, rel=0.01)
        assert dlat == pytest.approx(0.0, abs=0.001)
        assert ddist == pytest.approx(0.0, abs=0.001)

    def test_elliptic_orbit(self):
        """Test position calculation for elliptic orbit."""
        elem = OrbitalElements(
            name="EllipticPlanet",
            epoch_jd=2451545.0,  # J2000.0
            equinox_jd=2451545.0,
            equinox_is_jdate=False,
            mean_anomaly=TPolynomial(0.0, 0.0),  # At perihelion
            semi_axis=2.0,  # 2 AU
            eccentricity=TPolynomial(0.5, 0.0),  # High eccentricity
            arg_perihelion=TPolynomial(0.0, 0.0),
            asc_node=TPolynomial(0.0, 0.0),
            inclination=TPolynomial(0.0, 0.0),
        )

        pos = calc_orbital_position(elem, 2451545.0)  # At epoch
        lon, lat, dist, dlon, dlat, ddist = pos

        # At perihelion, distance = a(1-e) = 2(1-0.5) = 1 AU
        assert dist == pytest.approx(1.0, abs=0.001)

    def test_inclined_orbit(self):
        """Test position calculation for inclined orbit."""
        elem = OrbitalElements(
            name="InclinedPlanet",
            epoch_jd=2451545.0,  # J2000.0
            equinox_jd=2451545.0,
            equinox_is_jdate=False,
            mean_anomaly=TPolynomial(90.0, 0.0),  # 90 deg from perihelion
            semi_axis=1.0,
            eccentricity=TPolynomial(0.0, 0.0),  # Circular
            arg_perihelion=TPolynomial(0.0, 0.0),
            asc_node=TPolynomial(0.0, 0.0),
            inclination=TPolynomial(30.0, 0.0),  # 30 deg inclination
        )

        pos = calc_orbital_position(elem, 2451545.0)  # At epoch
        lon, lat, dist, dlon, dlat, ddist = pos

        # With i=30 deg and body at 90 deg from node, latitude should be significant
        # At u = 90 deg (ascending node + 90), latitude = arcsin(sin(i) * sin(u))
        # = arcsin(sin(30) * sin(90)) = arcsin(0.5 * 1) = 30 deg
        assert lat == pytest.approx(30.0, abs=1.0)

    def test_time_dependent_elements(self):
        """Test position calculation with time-dependent elements."""
        # Create an element with a high mean motion polynomial
        elem = OrbitalElements(
            name="FastMover",
            epoch_jd=2451545.0,  # J2000.0
            equinox_jd=2451545.0,
            equinox_is_jdate=False,
            # 100 + 36000*T means 36000 deg/century motion = ~0.986 deg/day secular
            mean_anomaly=TPolynomial(100.0, 36000.0),
            semi_axis=1.0,
            eccentricity=TPolynomial(0.0, 0.0),
            arg_perihelion=TPolynomial(0.0, 0.0),
            asc_node=TPolynomial(0.0, 0.0),
            inclination=TPolynomial(0.0, 0.0),
        )

        # At epoch
        pos1 = calc_orbital_position(elem, 2451545.0)
        # 1 century later
        pos2 = calc_orbital_position(elem, 2451545.0 + 36525.0)

        # The longitude should have advanced by 36000 deg (modulo 360)
        # 36000 % 360 = 0, so longitudes should be the same
        expected_advance = 36000.0 % 360.0
        actual_advance = (pos2[0] - pos1[0]) % 360.0
        assert actual_advance == pytest.approx(expected_advance, abs=1.0)

    def test_position_matches_cupido(self):
        """Test that calculated position is reasonable for Cupido."""
        seorbel_path = (
            Path(__file__).parent.parent / "ephe" / "seorbel.txt"
        )
        if not seorbel_path.exists():
            pytest.skip("seorbel.txt not found in reference ephemeris dir")

        elements = parse_orbital_elements(seorbel_path)
        cupido = get_orbital_body_by_name(elements, "Cupido")
        assert cupido is not None

        # Calculate position at J2000.0
        pos = calc_orbital_position(cupido, 2451545.0)
        lon, lat, dist, dlon, dlat, ddist = pos

        # Check reasonable values
        assert 0 <= lon < 360
        assert -90 <= lat <= 90
        assert dist > 0
        assert dist == pytest.approx(cupido.semi_axis, rel=0.1)

        # Mean motion should be very slow (Cupido has ~262 year period)
        # n = 360 / (a^1.5 * 365.25) deg/day
        expected_n = 360.0 / (cupido.semi_axis**1.5 * 365.25)
        assert dlon == pytest.approx(expected_n, rel=0.1)


class TestOrbitalElementsDataclass:
    """Tests for the OrbitalElements dataclass."""

    def test_get_mean_motion(self):
        """Test mean motion calculation from semi-major axis."""
        elem = OrbitalElements(
            name="Test",
            epoch_jd=2451545.0,
            equinox_jd=2451545.0,
            equinox_is_jdate=False,
            mean_anomaly=TPolynomial(0.0, 0.0),
            semi_axis=1.0,  # 1 AU = 1 year period
            eccentricity=TPolynomial(0.0, 0.0),
            arg_perihelion=TPolynomial(0.0, 0.0),
            asc_node=TPolynomial(0.0, 0.0),
            inclination=TPolynomial(0.0, 0.0),
        )

        n = elem.get_mean_motion()
        # For a = 1 AU, period = 1 year, n = 360/365.25 deg/day
        expected = 360.0 / 365.25
        assert n == pytest.approx(expected, rel=0.001)

    def test_line_number_tracking(self):
        """Test that line numbers are tracked correctly."""
        content = """\
# Comment
J2000, J2000, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0, Planet1
# Another comment
J2000, J2000, 0.0, 60.0, 0.0, 0.0, 0.0, 0.0, Planet2
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_orbital_elements(filepath)
            assert elements[0].line_number == 2  # After the first comment
            assert elements[1].line_number == 4  # After another comment
        finally:
            Path(filepath).unlink()


class TestBundledFictitiousOrbits:
    """Tests for the bundled fictitious_orbits.csv dataset."""

    def test_get_bundled_fictitious_orbits_path_exists(self):
        """Test that the bundled fictitious_orbits.csv path exists."""
        path = get_bundled_fictitious_orbits_path()
        assert path.exists(), "Bundled fictitious_orbits.csv should exist"
        assert path.is_file(), "Bundled fictitious_orbits.csv should be a file"
        assert path.name == "fictitious_orbits.csv", (
            "File should be named fictitious_orbits.csv"
        )

    def test_get_bundled_fictitious_orbits_path_in_package(self):
        """Test that the bundled file is inside the libephemeris package."""
        path = get_bundled_fictitious_orbits_path()
        assert "libephemeris" in str(path), "Path should contain 'libephemeris'"

    def test_load_bundled_fictitious_orbits_returns_list(self):
        """Test that load_bundled_fictitious_orbits returns a list of elements."""
        elements = load_bundled_fictitious_orbits()
        assert isinstance(elements, list), "Should return a list"
        assert len(elements) > 0, "Should have at least one element"
        for elem in elements:
            assert isinstance(elem, OrbitalElements), (
                "Each element should be OrbitalElements"
            )

    def test_load_bundled_fictitious_orbits_contains_expected_bodies(self):
        """Test that the dataset contains the expected hypothetical bodies."""
        elements = load_bundled_fictitious_orbits()
        names = [elem.name for elem in elements]

        # Hamburg School Uranian planets
        expected_bodies = [
            "Cupido",
            "Hades",
            "Zeus",
            "Kronos",
            "Apollon",
            "Admetos",
            "Vulkanus",
            "Poseidon",
        ]
        for body in expected_bodies:
            assert any(body.lower() in name.lower() for name in names), (
                f"Expected body '{body}' not found in bundled dataset"
            )

    def test_load_bundled_fictitious_orbits_cupido_elements(self):
        """Test that Cupido has correct orbital elements from bundled dataset."""
        elements = load_bundled_fictitious_orbits()
        cupido = get_orbital_body_by_name(elements, "Cupido")

        assert cupido is not None, "Cupido should be found"
        assert abs(cupido.semi_axis - 40.99837) < 0.001, (
            "Cupido semi-axis should match Witte/Neely value"
        )
        assert abs(cupido.eccentricity.constant - 0.00460) < 0.0001, (
            "Cupido eccentricity should match Witte/Neely value"
        )
        assert cupido.epoch_jd == pytest.approx(2415020.0), (
            "Cupido epoch should be J1900"
        )

    def test_load_bundled_fictitious_orbits_transpluto_elements(self):
        """Test that Isis-Transpluto has correct orbital elements (Strubell 1952)."""
        elements = load_bundled_fictitious_orbits()
        transpluto = get_orbital_body_by_name(elements, "Isis-Transpluto")

        assert transpluto is not None, "Isis-Transpluto should be found"
        assert abs(transpluto.semi_axis - 77.775) < 0.001, (
            "Transpluto semi-axis should match Strubell 1952 value"
        )
        assert abs(transpluto.eccentricity.constant - 0.3) < 0.01, (
            "Transpluto eccentricity should match Strubell 1952 value"
        )

    def test_load_bundled_fictitious_orbits_geocentric_bodies(self):
        """Test that geocentric bodies are correctly identified."""
        elements = load_bundled_fictitious_orbits()

        # Waldemath: geocentric (orbits Earth)
        waldemath = get_orbital_body_by_name(elements, "Waldemath")
        assert waldemath is not None, "Waldemath should be found"
        assert waldemath.is_geocentric, "Waldemath should be geocentric"

        # Selena: geocentric
        selena = get_orbital_body_by_name(elements, "Selena")
        assert selena is not None, "Selena should be found"
        assert selena.is_geocentric, "Selena should be geocentric"

    def test_load_bundled_fictitious_orbits_can_calculate_position(self):
        """Test that positions can be calculated from the bundled elements."""
        elements = load_bundled_fictitious_orbits()
        cupido = get_orbital_body_by_name(elements, "Cupido")

        assert cupido is not None, "Cupido should be found"

        pos = calc_orbital_position(cupido, 2451545.0)

        assert 0 <= pos[0] < 360, "Longitude should be in range [0, 360)"
        assert -90 <= pos[1] <= 90, "Latitude should be in range [-90, 90]"
        assert pos[2] > 0, "Distance should be positive"

    def test_load_bundled_fictitious_orbits_vulcan_tpoly(self):
        """Test that Vulcan's T-polynomial elements are parsed correctly."""
        elements = load_bundled_fictitious_orbits()
        vulcan = get_orbital_body_by_name(elements, "Vulcan")

        assert vulcan is not None, "Vulcan should be found"
        assert vulcan.equinox_is_jdate is True
        assert vulcan.mean_anomaly.constant == pytest.approx(252.8987988)
        assert vulcan.mean_anomaly.linear == pytest.approx(707550.7341)
        assert vulcan.arg_perihelion.linear == pytest.approx(1670.056)
        assert vulcan.asc_node.linear == pytest.approx(-1670.056)

    # --- Backward-compatibility shims ---

    def test_get_bundled_seorbel_path_delegates_to_new_function(self):
        """Test that get_bundled_seorbel_path() is a shim to the new API."""
        new_path = get_bundled_fictitious_orbits_path()
        compat_path = get_bundled_seorbel_path()
        assert compat_path == new_path, (
            "Compat shim should return the same path as the new function"
        )

    def test_load_bundled_seorbel_delegates_to_new_function(self):
        """Test that load_bundled_seorbel() is a shim to the new API."""
        new_elements = load_bundled_fictitious_orbits()
        compat_elements = load_bundled_seorbel()
        assert len(compat_elements) == len(new_elements), (
            "Compat shim should return same number of elements"
        )
        for a, b in zip(compat_elements, new_elements):
            assert a.name == b.name, "Element names should match"


class TestParseTPolynomialEdgeCases:
    """Additional edge case tests for _parse_t_polynomial helper.

    These tests cover various spacing patterns, format variations,
    and edge cases that might occur in seorbel.txt files.
    """

    # Tests for various spacing patterns
    def test_polynomial_extra_spaces(self):
        """Test parsing polynomial with extra spaces around operators."""
        poly = _parse_t_polynomial("252.8987988  +  707550.7341  *  T")
        assert poly.constant == pytest.approx(252.8987988)
        assert poly.linear == pytest.approx(707550.7341)

    def test_polynomial_leading_trailing_spaces(self):
        """Test parsing polynomial with leading and trailing spaces."""
        poly = _parse_t_polynomial("   100.5 + 50.25 * T   ")
        assert poly.constant == pytest.approx(100.5)
        assert poly.linear == pytest.approx(50.25)

    def test_polynomial_no_space_before_minus(self):
        """Test parsing polynomial without space before minus."""
        poly = _parse_t_polynomial("100.0- 50.0 * T")
        assert poly.constant == pytest.approx(100.0)
        assert poly.linear == pytest.approx(-50.0)

    def test_polynomial_no_space_after_minus(self):
        """Test parsing polynomial without space after minus."""
        poly = _parse_t_polynomial("100.0 -50.0 * T")
        assert poly.constant == pytest.approx(100.0)
        assert poly.linear == pytest.approx(-50.0)

    def test_polynomial_mixed_spacing_patterns(self):
        """Test various mixed spacing patterns."""
        test_cases = [
            ("100.0+50.0*T", 100.0, 50.0),
            ("100.0 +50.0*T", 100.0, 50.0),
            ("100.0+ 50.0*T", 100.0, 50.0),
            ("100.0 + 50.0*T", 100.0, 50.0),
            ("100.0+50.0 *T", 100.0, 50.0),
            ("100.0+50.0* T", 100.0, 50.0),
            ("100.0+50.0 * T", 100.0, 50.0),
        ]
        for expr, expected_const, expected_linear in test_cases:
            poly = _parse_t_polynomial(expr)
            assert poly.constant == pytest.approx(expected_const), f"Failed for: {expr}"
            assert poly.linear == pytest.approx(expected_linear), f"Failed for: {expr}"

    def test_polynomial_negative_constant(self):
        """Test polynomial with negative constant."""
        poly = _parse_t_polynomial("-100.5 + 50.25 * T")
        assert poly.constant == pytest.approx(-100.5)
        assert poly.linear == pytest.approx(50.25)

    def test_polynomial_both_negative(self):
        """Test polynomial with negative constant and negative rate."""
        poly = _parse_t_polynomial("-100.5 - 50.25 * T")
        assert poly.constant == pytest.approx(-100.5)
        assert poly.linear == pytest.approx(-50.25)

    def test_polynomial_lowercase_t(self):
        """Test that lowercase 't' is also accepted."""
        poly = _parse_t_polynomial("100.0 + 50.0 * t")
        assert poly.constant == pytest.approx(100.0)
        assert poly.linear == pytest.approx(50.0)

    def test_polynomial_integer_values(self):
        """Test parsing polynomial with integer values (no decimal point)."""
        poly = _parse_t_polynomial("100 + 50 * T")
        assert poly.constant == pytest.approx(100.0)
        assert poly.linear == pytest.approx(50.0)

    def test_polynomial_zero_linear(self):
        """Test parsing polynomial with zero linear term."""
        poly = _parse_t_polynomial("100.5 + 0.0 * T")
        assert poly.constant == pytest.approx(100.5)
        assert poly.linear == pytest.approx(0.0)

    def test_polynomial_very_large_values(self):
        """Test parsing polynomial with very large values."""
        poly = _parse_t_polynomial("999999.999 + 123456789.0 * T")
        assert poly.constant == pytest.approx(999999.999)
        assert poly.linear == pytest.approx(123456789.0)

    def test_polynomial_very_small_values(self):
        """Test parsing polynomial with very small (but non-zero) values."""
        poly = _parse_t_polynomial("0.000001 + 0.000000001 * T")
        assert poly.constant == pytest.approx(0.000001)
        assert poly.linear == pytest.approx(0.000000001)

    def test_polynomial_scientific_notation(self):
        """Test that scientific notation in values is handled correctly."""
        # Note: This may or may not be supported depending on implementation
        # Testing to document behavior
        try:
            poly = _parse_t_polynomial("1.5e2")
            assert poly.constant == pytest.approx(150.0)
            assert poly.linear == 0.0
        except ValueError:
            pytest.skip("Scientific notation not supported")


class TestParseTPolynomialErrors:
    """Tests for error handling in _parse_t_polynomial."""

    def test_invalid_expression_raises_error(self):
        """Test that completely invalid expressions raise ValueError."""
        with pytest.raises(ValueError):
            _parse_t_polynomial("not a number")

    def test_empty_string_raises_error(self):
        """Test that empty string raises ValueError."""
        with pytest.raises(ValueError):
            _parse_t_polynomial("")

    def test_whitespace_only_raises_error(self):
        """Test that whitespace-only string raises ValueError."""
        with pytest.raises(ValueError):
            _parse_t_polynomial("   ")

    def test_duplicate_operator_handled_gracefully(self):
        """Test that duplicate operators are handled gracefully.

        The parser is tolerant and handles edge cases like '++' by
        treating it as a single '+'. This documents expected behavior.
        """
        poly = _parse_t_polynomial("100 ++ 50 * T")
        assert poly.constant == pytest.approx(100.0)
        assert poly.linear == pytest.approx(50.0)

    def test_bare_t_term_handled_gracefully(self):
        """Test that a bare T term (with implied coefficient) is handled.

        The parser treats '+ * T' as having coefficient 1.0.
        This documents expected behavior.
        """
        poly = _parse_t_polynomial("+ * T")
        assert poly.constant == pytest.approx(0.0)
        assert poly.linear == pytest.approx(1.0)


class TestParseOrbitalElementsLine:
    """Direct tests for _parse_orbital_elements_line function."""

    def test_parse_simple_line(self):
        """Test parsing a simple valid line."""
        line = "J2000, J2000, 0.0, 100.0, 0.1, 45.0, 30.0, 5.0, TestPlanet"
        elem = _parse_orbital_elements_line(line, 1)

        assert elem is not None
        assert elem.name == "TestPlanet"
        assert elem.epoch_jd == 2451545.0  # J2000
        assert elem.equinox_jd == 2451545.0
        assert elem.semi_axis == 100.0
        assert elem.eccentricity.constant == pytest.approx(0.1)
        assert elem.arg_perihelion.constant == pytest.approx(45.0)
        assert elem.asc_node.constant == pytest.approx(30.0)
        assert elem.inclination.constant == pytest.approx(5.0)
        assert elem.mean_anomaly.constant == pytest.approx(0.0)

    def test_parse_line_with_t_polynomial(self):
        """Test parsing a line with T-polynomial expressions."""
        line = "J1900, JDATE, 252.8987988 + 707550.7341 * T, 0.13744, 0.019, 322.212069+1670.056*T, 47.787931-1670.056*T, 7.5, Vulcan"
        elem = _parse_orbital_elements_line(line, 1)

        assert elem is not None
        assert elem.name == "Vulcan"
        assert elem.epoch_jd == 2415020.0  # J1900
        assert elem.equinox_is_jdate is True
        assert elem.mean_anomaly.constant == pytest.approx(252.8987988)
        assert elem.mean_anomaly.linear == pytest.approx(707550.7341)
        assert elem.arg_perihelion.constant == pytest.approx(322.212069)
        assert elem.arg_perihelion.linear == pytest.approx(1670.056)
        assert elem.asc_node.constant == pytest.approx(47.787931)
        assert elem.asc_node.linear == pytest.approx(-1670.056)

    def test_parse_geocentric_body(self):
        """Test parsing a geocentric body line."""
        line = "J2000, JDATE, 248.8833, 0.0029833, 0.0, 0.0, 0.0, 0.0, Waldemath, geo"
        elem = _parse_orbital_elements_line(line, 1)

        assert elem is not None
        assert elem.name == "Waldemath"
        assert elem.is_geocentric is True

    def test_parse_line_with_inline_comment(self):
        """Test parsing a line with inline comment."""
        line = "J1900, J1900, 163.7409, 40.99837, 0.00460, 171.4333, 129.8325, 1.0833, Cupido   # 1"
        elem = _parse_orbital_elements_line(line, 1)

        assert elem is not None
        assert elem.name == "Cupido"
        assert elem.semi_axis == pytest.approx(40.99837)

    def test_parse_line_numeric_epoch(self):
        """Test parsing a line with numeric Julian Day epoch."""
        line = "2368547.66, 2431456.5, 0.0, 77.775, 0.3, 0.7, 0, 0, Isis-Transpluto"
        elem = _parse_orbital_elements_line(line, 1)

        assert elem is not None
        assert elem.name == "Isis-Transpluto"
        assert elem.epoch_jd == pytest.approx(2368547.66)
        assert elem.equinox_jd == pytest.approx(2431456.5)

    def test_line_number_tracking(self):
        """Test that line number is tracked correctly."""
        line = "J2000, J2000, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0, TestBody"
        elem = _parse_orbital_elements_line(line, 42)

        assert elem is not None
        assert elem.line_number == 42

    def test_parse_line_b1950_epoch(self):
        """Test parsing a line with B1950 epoch."""
        line = "B1950, J2000, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0, TestBody"
        elem = _parse_orbital_elements_line(line, 1)

        assert elem is not None
        assert elem.epoch_jd == pytest.approx(2433282.42345905, abs=0.0001)

    def test_parse_line_with_spaces_in_name(self):
        """Test parsing a body with spaces in the name."""
        line = "J2000, J2000, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0, Selena/White Moon, geo"
        elem = _parse_orbital_elements_line(line, 1)

        assert elem is not None
        assert elem.name == "Selena/White Moon"
        assert elem.is_geocentric is True


class TestParseOrbitalElementsLineErrors:
    """Tests for error handling in _parse_orbital_elements_line."""

    def test_too_few_fields_raises_error(self):
        """Test that lines with too few fields raise ValueError."""
        line = "J2000, J2000, 0.0, 100.0, 0.1, 45.0, 30.0"  # Only 7 fields
        with pytest.raises(ValueError, match="Expected at least 9"):
            _parse_orbital_elements_line(line, 1)

    def test_invalid_epoch_raises_error(self):
        """Test that invalid epoch raises ValueError."""
        line = "INVALID_EPOCH, J2000, 0.0, 100.0, 0.1, 45.0, 30.0, 5.0, TestBody"
        with pytest.raises(ValueError, match="Cannot parse epoch"):
            _parse_orbital_elements_line(line, 1)

    def test_jdate_as_epoch_raises_error(self):
        """Test that JDATE as epoch raises ValueError (epoch cannot be JDATE)."""
        line = "JDATE, J2000, 0.0, 100.0, 0.1, 45.0, 30.0, 5.0, TestBody"
        with pytest.raises(ValueError, match="Epoch cannot be JDATE"):
            _parse_orbital_elements_line(line, 1)

    def test_invalid_semi_axis_raises_error(self):
        """Test that invalid semi-major axis raises ValueError."""
        line = "J2000, J2000, 0.0, not_a_number, 0.1, 45.0, 30.0, 5.0, TestBody"
        with pytest.raises(ValueError, match="Cannot parse semi-major axis"):
            _parse_orbital_elements_line(line, 1)

    def test_invalid_mean_anomaly_raises_error(self):
        """Test that invalid mean anomaly expression raises ValueError."""
        line = "J2000, J2000, invalid_expr, 100.0, 0.1, 45.0, 30.0, 5.0, TestBody"
        with pytest.raises(ValueError):
            _parse_orbital_elements_line(line, 1)


class TestParseOrbitalElementsRobustness:
    """Integration tests for parse_orbital_elements robustness with various formats."""

    def test_parse_all_documented_formats(self):
        """Test parsing a file with all documented seorbel.txt formats."""
        content = """\
# Test file with various documented formats

# Standard format with J2000 epoch
J2000, J2000, 0.0, 100.0, 0.1, 45.0, 30.0, 5.0, StandardBody

# J1900 epoch (common for Uranian planets)
J1900, J1900, 163.7409, 40.99837, 0.00460, 171.4333, 129.8325, 1.0833, CupidoStyle

# B1950 epoch
B1950, J2000, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0, B1950Body

# JDATE equinox (date of the day)
J2000, JDATE, 0.0, 60.0, 0.0, 0.0, 0.0, 0.0, JdateEquinox

# Numeric epochs
2451545.0, 2451545.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, NumericEpoch

# T-polynomial mean anomaly with various spacing
J1900, JDATE, 252.8987988 + 707550.7341 * T, 0.13744, 0.019, 0.0, 0.0, 7.5, TPolySpaced
J1900, JDATE, 100.0+50.0*T, 0.13744, 0.019, 0.0, 0.0, 7.5, TPolyCompact
J1900, JDATE, 100.0 - 50.0 * T, 0.13744, 0.019, 0.0, 0.0, 7.5, TPolyMinus

# Multiple T-polynomials in one line
J1900, JDATE, 252.0+707550.0*T, 0.13744, 0.019, 322.0+1670.0*T, 47.0-1670.0*T, 7.5, MultiPoly

# Geocentric body
J2000, JDATE, 248.8833, 0.0029833, 0.0, 0.0, 0.0, 0.0, GeoBody, geo

# Inline comment
J2000, J2000, 0.0, 80.0, 0.0, 0.0, 0.0, 0.0, CommentBody # This is a comment

"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_orbital_elements(filepath)

            # Should have parsed all bodies
            assert len(elements) == 11

            # Verify key bodies were parsed correctly
            names = [e.name for e in elements]
            assert "StandardBody" in names
            assert "CupidoStyle" in names
            assert "B1950Body" in names
            assert "JdateEquinox" in names
            assert "NumericEpoch" in names
            assert "TPolySpaced" in names
            assert "TPolyCompact" in names
            assert "TPolyMinus" in names
            assert "MultiPoly" in names
            assert "GeoBody" in names
            assert "CommentBody" in names

            # Verify T-polynomial parsing
            tpoly_spaced = get_orbital_body_by_name(elements, "TPolySpaced")
            assert tpoly_spaced is not None, "TPolySpaced not found"
            assert tpoly_spaced.mean_anomaly.constant == pytest.approx(252.8987988)
            assert tpoly_spaced.mean_anomaly.linear == pytest.approx(707550.7341)

            tpoly_compact = get_orbital_body_by_name(elements, "TPolyCompact")
            assert tpoly_compact is not None, "TPolyCompact not found"
            assert tpoly_compact.mean_anomaly.constant == pytest.approx(100.0)
            assert tpoly_compact.mean_anomaly.linear == pytest.approx(50.0)

            tpoly_minus = get_orbital_body_by_name(elements, "TPolyMinus")
            assert tpoly_minus is not None, "TPolyMinus not found"
            assert tpoly_minus.mean_anomaly.constant == pytest.approx(100.0)
            assert tpoly_minus.mean_anomaly.linear == pytest.approx(-50.0)

            multi_poly = get_orbital_body_by_name(elements, "MultiPoly")
            assert multi_poly is not None, "MultiPoly not found"
            assert multi_poly.arg_perihelion.linear == pytest.approx(1670.0)
            assert multi_poly.asc_node.linear == pytest.approx(-1670.0)

            # Verify geocentric body
            geo_body = get_orbital_body_by_name(elements, "GeoBody")
            assert geo_body is not None, "GeoBody not found"
            assert geo_body.is_geocentric is True

        finally:
            Path(filepath).unlink()

    def test_parse_whitespace_variations(self):
        """Test parsing handles various whitespace patterns."""
        content = """\
# Different whitespace patterns
J2000,J2000,0.0,100.0,0.1,45.0,30.0,5.0,NoSpaces
J2000 , J2000 , 0.0 , 100.0 , 0.1 , 45.0 , 30.0 , 5.0 , ExtraSpaces
J2000,  J2000,  0.0,  100.0,  0.1,  45.0,  30.0,  5.0,  TwoSpaces
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_orbital_elements(filepath)
            assert len(elements) == 3

            for elem in elements:
                assert elem.semi_axis == pytest.approx(100.0)
                assert elem.eccentricity.constant == pytest.approx(0.1)
        finally:
            Path(filepath).unlink()

    def test_parse_malformed_file_with_error_context(self):
        """Test that parsing errors include helpful context."""
        content = """\
J2000, J2000, 0.0, 100.0, 0.1, 45.0, 30.0, 5.0, ValidBody
J2000, J2000, invalid_value, 100.0, 0.1, 45.0, 30.0, 5.0, InvalidBody
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            with pytest.raises(ValueError, match="line 2"):
                parse_orbital_elements(filepath)
        finally:
            Path(filepath).unlink()

    def test_empty_file(self):
        """Test parsing an empty file returns empty list."""
        content = ""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_orbital_elements(filepath)
            assert elements == []
        finally:
            Path(filepath).unlink()

    def test_comments_only_file(self):
        """Test parsing a file with only comments returns empty list."""
        content = """\
# This is a comment
# Another comment
   # Indented comment
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_orbital_elements(filepath)
            assert elements == []
        finally:
            Path(filepath).unlink()
