"""
API compatibility and return type tests.

Verifies that all public API functions return the expected types
and formats, ensuring compatibility with the reference ephemeris.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    MARS,
    JUPITER,
    EARTH,
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    OSCU_APOG,
    INTP_APOG,
    INTP_PERG,
    CHIRON,
    FLG_SPEED,
    SIDM_LAHIRI,
    GREG_CAL,
)


class TestCalcUtReturnFormat:
    """Test calc_ut return format."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MERCURY, "Mercury"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (MEAN_NODE, "MeanNode"),
            (TRUE_NODE, "TrueNode"),
            (MEAN_APOG, "MeanApog"),
            (OSCU_APOG, "OscuApog"),
            (INTP_APOG, "IntpApog"),
            (INTP_PERG, "IntpPerg"),
            (CHIRON, "Chiron"),
            (EARTH, "Earth"),
        ],
    )
    def test_calc_ut_returns_tuple_pair(self, body_id: int, name: str):
        """calc_ut returns (result_tuple, retflag)."""
        jd = 2451545.0
        retval = swe.calc_ut(jd, body_id, FLG_SPEED)
        assert isinstance(retval, tuple), f"{name}: not a tuple"
        assert len(retval) == 2, f"{name}: expected 2 elements"
        result, retflag = retval
        assert isinstance(result, tuple), f"{name}: result not a tuple"
        assert len(result) == 6, f"{name}: result has {len(result)} elements"
        assert isinstance(retflag, int), f"{name}: retflag not int"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
        ],
    )
    def test_calc_ut_result_elements_are_float(self, body_id: int, name: str):
        """All result elements should be native Python floats."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
        for i, val in enumerate(result):
            assert isinstance(val, float), (
                f"{name}: result[{i}] is {type(val).__name__}, not float"
            )

    @pytest.mark.unit
    def test_earth_geocentric_is_zeros(self):
        """Earth geocentric should return all zeros."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, EARTH, 0)
        for i, val in enumerate(result):
            assert val == 0.0, f"Earth geocentric result[{i}] = {val}"


class TestCalcReturnFormat:
    """Test calc (TT) return format."""

    @pytest.mark.unit
    def test_calc_returns_tuple_pair(self):
        """calc returns (result_tuple, retflag)."""
        jd_tt = 2451545.0
        retval = swe.calc(jd_tt, SUN, FLG_SPEED)
        assert isinstance(retval, tuple)
        assert len(retval) == 2
        result, retflag = retval
        assert len(result) == 6
        assert isinstance(retflag, int)


class TestJuldayRevjul:
    """Test Julian Day conversion functions."""

    @pytest.mark.unit
    def test_julday_returns_float(self):
        """julday returns a float."""
        jd = swe.julday(2000, 1, 1, 12.0)
        assert isinstance(jd, float)
        assert jd == pytest.approx(2451545.0, abs=0.5)

    @pytest.mark.unit
    def test_revjul_returns_4_values(self):
        """revjul returns (year, month, day, hour)."""
        result = swe.revjul(2451545.0)
        assert len(result) == 4
        y, m, d, h = result
        assert y == 2000
        assert m == 1
        assert d == 1
        assert abs(h - 12.0) < 0.001

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "year,month,day,hour",
        [
            (2000, 1, 1, 12.0),
            (1999, 12, 31, 0.0),
            (2024, 7, 15, 18.5),
            (1900, 1, 1, 0.0),
            (1582, 10, 15, 12.0),
        ],
    )
    def test_julday_revjul_roundtrip(
        self, year: int, month: int, day: int, hour: float
    ):
        """julday -> revjul should recover the original date."""
        jd = swe.julday(year, month, day, hour)
        y, m, d, h = swe.revjul(jd)
        assert y == year
        assert m == month
        assert d == day
        assert abs(h - hour) < 0.001


class TestDeltatFormat:
    """Test Delta T function return format."""

    @pytest.mark.unit
    def test_deltat_returns_float(self):
        """deltat returns a float (days)."""
        dt = swe.deltat(2451545.0)
        assert isinstance(dt, float)
        assert math.isfinite(dt)


class TestSidtimeFormat:
    """Test sidereal time function."""

    @pytest.mark.unit
    def test_sidtime_returns_float(self):
        """sidtime returns a float (hours)."""
        st = swe.sidtime(2451545.0)
        assert isinstance(st, float)
        assert 0 <= st < 24, f"Sidereal time {st} not in [0, 24)"

    @pytest.mark.unit
    @pytest.mark.parametrize("jd", [2451545.0, 2440587.5, 2460000.0])
    def test_sidtime_in_range(self, jd: float):
        """Sidereal time always in [0, 24) hours."""
        st = swe.sidtime(jd)
        assert 0 <= st < 24


class TestHousesFormat:
    """Test house calculation return format."""

    @pytest.mark.unit
    def test_houses_returns_cusps_and_ascmc(self):
        """houses returns (cusps, ascmc)."""
        jd = 2451545.0
        result = swe.houses(jd, 41.9, 12.5, ord("P"))
        assert len(result) == 2
        cusps, ascmc = result
        assert len(cusps) >= 12
        assert len(ascmc) >= 2

    @pytest.mark.unit
    def test_houses_cusps_are_float(self):
        """All cusp values should be floats."""
        jd = 2451545.0
        cusps, ascmc = swe.houses(jd, 41.9, 12.5, ord("P"))
        for i, val in enumerate(cusps[:12]):
            assert isinstance(val, float), f"Cusp {i + 1} is {type(val).__name__}"


class TestUtcConversion:
    """Test UTC/JD conversion functions."""

    @pytest.mark.unit
    def test_utc_to_jd_returns_tuple(self):
        """utc_to_jd returns a tuple of JD values."""
        result = swe.utc_to_jd(2000, 1, 1, 12, 0, 0.0, GREG_CAL)
        assert isinstance(result, tuple)
        assert len(result) >= 2

    @pytest.mark.unit
    def test_jdut1_to_utc_returns_tuple(self):
        """jdut1_to_utc returns UTC components."""
        result = swe.jdut1_to_utc(2451545.0, GREG_CAL)
        assert isinstance(result, tuple)
        assert len(result) >= 5  # year, month, day, hour, min, sec


class TestGetPlanetName:
    """Test planet name function."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,expected_fragment",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MERCURY, "Merc"),
            (MARS, "Mars"),
            (JUPITER, "Jup"),
            (MEAN_NODE, "Node"),
        ],
    )
    def test_get_planet_name(self, body_id: int, expected_fragment: str):
        """get_planet_name returns a string containing the planet name."""
        name = swe.get_planet_name(body_id)
        assert isinstance(name, str)
        assert len(name) > 0
        assert expected_fragment.lower() in name.lower(), (
            f"Body {body_id}: got '{name}', expected to contain '{expected_fragment}'"
        )


class TestAyanamshaFormat:
    """Test ayanamsha function return formats."""

    @pytest.mark.unit
    def test_get_ayanamsa_ut_returns_float(self):
        """get_ayanamsa_ut returns a float."""
        swe.set_sid_mode(SIDM_LAHIRI)
        ayan = swe.get_ayanamsa_ut(2451545.0)
        assert isinstance(ayan, float)
        assert math.isfinite(ayan)

    @pytest.mark.unit
    def test_get_ayanamsa_ex_ut_returns_tuple(self):
        """get_ayanamsa_ex_ut returns (retflag, ayanamsa)."""
        swe.set_sid_mode(SIDM_LAHIRI)
        result = swe.get_ayanamsa_ex_ut(2451545.0, 0)
        assert isinstance(result, tuple)
        assert len(result) == 2
        retflag, ayan = result
        assert isinstance(retflag, int)
        assert isinstance(ayan, float)


class TestOrbitalElementsFormat:
    """Test orbital elements return format."""

    @pytest.mark.unit
    def test_get_orbital_elements_returns_50(self):
        """get_orbital_elements returns 50-element tuple."""
        result = swe.get_orbital_elements(2451545.0, MARS, 0)
        assert isinstance(result, tuple)
        assert len(result) == 50

    @pytest.mark.unit
    def test_get_orbital_elements_ut_returns_50(self):
        """get_orbital_elements_ut returns 50-element tuple."""
        result = swe.get_orbital_elements_ut(2451545.0, MARS, 0)
        assert isinstance(result, tuple)
        assert len(result) == 50


class TestCotransFormat:
    """Test coordinate transform return format."""

    @pytest.mark.unit
    def test_cotrans_returns_3_tuple(self):
        """cotrans returns a 3-element tuple."""
        result = swe.cotrans((100.0, 20.0, 1.0), 23.44)
        assert isinstance(result, tuple)
        assert len(result) == 3
        for val in result:
            assert isinstance(val, float)
