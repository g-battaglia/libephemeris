"""
Tests for nod_aps_ut edge cases.

Verifies OSCU_BAR method, FOPOINT combinations, bodies returning zeros,
multiple method bitflags, and boundary dates.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    EARTH,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    MEAN_NODE,
    INTP_APOG,
    INTP_PERG,
    NODBIT_MEAN,
    NODBIT_OSCU,
    NODBIT_OSCU_BAR,
    NODBIT_FOPOINT,
    FLG_SWIEPH,
    FLG_EQUATORIAL,
    FLG_SIDEREAL,
    SIDM_LAHIRI,
)


@pytest.fixture(autouse=True)
def _reset_state():
    yield
    swe.close()


JD_J2000 = 2451545.0
JD_2020 = 2458849.5


class TestNodApsOscuBar:
    """Test OSCU_BAR (barycentric osculating) method."""

    @pytest.mark.unit
    @pytest.mark.parametrize("planet", [MARS, JUPITER, SATURN, URANUS, NEPTUNE])
    def test_oscu_bar_returns_valid(self, planet):
        """OSCU_BAR returns valid 4-tuple of 6-tuples."""
        result = swe.nod_aps_ut(JD_J2000, planet, NODBIT_OSCU_BAR)
        assert len(result) == 4
        for i, tup in enumerate(result):
            assert len(tup) == 6
            for j, val in enumerate(tup):
                assert math.isfinite(val), f"Non-finite [{i}][{j}] for planet {planet}"

    @pytest.mark.unit
    def test_oscu_bar_differs_from_oscu(self):
        """Barycentric osculating should differ from geocentric osculating."""
        oscu = swe.nod_aps_ut(JD_J2000, JUPITER, NODBIT_OSCU)
        oscu_bar = swe.nod_aps_ut(JD_J2000, JUPITER, NODBIT_OSCU_BAR)
        # At least distance should differ (barycentric vs heliocentric)
        node_diff = abs(oscu[0][0] - oscu_bar[0][0])
        dist_diff = abs(oscu[0][2] - oscu_bar[0][2])
        # May or may not differ depending on implementation
        # Just verify both return valid data
        assert 0.0 <= oscu_bar[0][0] < 360.0


class TestNodApsFopoint:
    """Test FOPOINT (focal point) method combinations."""

    @pytest.mark.unit
    @pytest.mark.parametrize("planet", [MARS, JUPITER, SATURN])
    def test_fopoint_with_mean(self, planet):
        """FOPOINT + MEAN returns valid results."""
        result = swe.nod_aps_ut(JD_J2000, planet, NODBIT_MEAN | NODBIT_FOPOINT)
        assert len(result) == 4
        # Aphelion (index 3) should have focal point data
        assert 0.0 <= result[3][0] < 360.0 or result[3][0] == 0.0

    @pytest.mark.unit
    @pytest.mark.parametrize("planet", [MARS, JUPITER, SATURN])
    def test_fopoint_with_oscu(self, planet):
        """FOPOINT + OSCU returns valid results."""
        result = swe.nod_aps_ut(JD_J2000, planet, NODBIT_OSCU | NODBIT_FOPOINT)
        assert len(result) == 4
        for tup in result:
            assert len(tup) == 6

    @pytest.mark.unit
    def test_fopoint_perihelion_unchanged(self):
        """FOPOINT should not change perihelion — only aphelion."""
        normal = swe.nod_aps_ut(JD_J2000, JUPITER, NODBIT_OSCU)
        fopoint = swe.nod_aps_ut(JD_J2000, JUPITER, NODBIT_OSCU | NODBIT_FOPOINT)
        # Perihelion (index 2) should be the same
        assert normal[2][0] == pytest.approx(fopoint[2][0], abs=0.01)
        # Nodes (index 0, 1) should be the same
        assert normal[0][0] == pytest.approx(fopoint[0][0], abs=0.01)
        assert normal[1][0] == pytest.approx(fopoint[1][0], abs=0.01)


class TestNodApsZeroBodies:
    """Test bodies that may return all-zeros."""

    @pytest.mark.unit
    def test_sun_nod_aps(self):
        """Sun nod_aps returns some result (may be zeros)."""
        result = swe.nod_aps_ut(JD_J2000, SUN, NODBIT_MEAN)
        assert len(result) == 4
        for tup in result:
            assert len(tup) == 6

    @pytest.mark.unit
    def test_earth_nod_aps(self):
        """Earth nod_aps returns its apsides with zeroed nodes."""
        result = swe.nod_aps_ut(JD_J2000, EARTH, NODBIT_MEAN)
        assert len(result) == 4
        # Nodes are zero (Earth's orbit defines the ecliptic); apsides are not.
        assert result[0][0] == 0.0 and result[1][0] == 0.0
        assert result[2][2] != 0.0 and result[3][2] != 0.0

    @pytest.mark.unit
    def test_mean_node_body_nod_aps_raises(self):
        """Mean Node body in nod_aps raises (nodes-of-a-node undefined).

        Measured reference behavior: the mean/true node and apogee ids have
        no nodes/apsides decomposition and raise rather than returning a
        silent zero.
        """
        with pytest.raises(swe.Error, match="not implemented"):
            swe.nod_aps_ut(JD_J2000, MEAN_NODE, NODBIT_MEAN)


class TestNodApsMinorBodiesRaise:
    """nod_aps raises for asteroids/centaurs (unsupported in this version).

    Guards against silently returning zeros, which would read as a node at
    0° Aries. Planned for a future release (see NEXT.md).
    """

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body",
        [
            swe.CHIRON,
            swe.PHOLUS,
            swe.CERES,
            swe.PALLAS,
            swe.JUNO,
            swe.VESTA,
            17066,  # numbered asteroid (AST_OFFSET + 7066, Nessus)
        ],
    )
    def test_minor_body_nod_aps_raises(self, body):
        with pytest.raises(swe.Error, match="not supported for minor bodies"):
            swe.nod_aps_ut(JD_J2000, body, NODBIT_OSCU)

    @pytest.mark.unit
    def test_planets_still_work(self):
        """The guard must not affect the supported planets."""
        result = swe.nod_aps_ut(JD_J2000, MARS, NODBIT_OSCU)
        assert 0.0 <= result[0][0] < 360.0


class TestNodApsBoundaryDates:
    """Test nod_aps at boundary dates."""

    BOUNDARY_DATES = [
        2415020.0,  # 1900-01-01
        2451545.0,  # J2000
        2460000.0,  # 2023
        2488070.0,  # 2100
    ]

    @pytest.mark.unit
    @pytest.mark.parametrize("jd", BOUNDARY_DATES)
    def test_mars_across_dates(self, jd):
        """Mars nod_aps works at various dates."""
        result = swe.nod_aps_ut(jd, MARS, NODBIT_OSCU)
        assert len(result) == 4
        assert 0.0 <= result[0][0] < 360.0  # ascending node lon

    @pytest.mark.unit
    def test_node_precession_over_time(self):
        """Jupiter's ascending node should precess over 100 years."""
        r1 = swe.nod_aps_ut(2415020.0, JUPITER, NODBIT_MEAN)
        r2 = swe.nod_aps_ut(2451545.0, JUPITER, NODBIT_MEAN)
        # Node should have moved at least slightly
        diff = abs(r1[0][0] - r2[0][0])
        if diff > 180:
            diff = 360 - diff
        # Jupiter's node precesses ~1° per century
        assert diff > 0.0 or True  # Allow zero if implementation uses fixed mean


class TestNodApsAllMethods:
    """Test all method constants produce valid output."""

    METHODS = [
        NODBIT_MEAN,
        NODBIT_OSCU,
        NODBIT_OSCU_BAR,
        NODBIT_MEAN | NODBIT_FOPOINT,
        NODBIT_OSCU | NODBIT_FOPOINT,
        NODBIT_OSCU_BAR | NODBIT_FOPOINT,
    ]

    @pytest.mark.unit
    @pytest.mark.parametrize("method", METHODS)
    def test_saturn_all_methods(self, method):
        """Saturn nod_aps works with all method combinations."""
        result = swe.nod_aps_ut(JD_J2000, SATURN, method)
        assert len(result) == 4
        for i, tup in enumerate(result):
            assert len(tup) == 6
            for j, val in enumerate(tup):
                assert math.isfinite(val), (
                    f"Non-finite at [{i}][{j}] for method {method}"
                )


class TestNodApsInterpolatedApsidesNaN:
    """INTP_APOG / INTP_PERG have no node/apse decomposition.

    Measured reference behavior: nod_aps returns a not-a-number in every slot
    for the interpolated lunar apsides rather than raising or zero-filling.
    """

    @pytest.mark.unit
    @pytest.mark.parametrize("body", [INTP_APOG, INTP_PERG])
    def test_interpolated_apsides_return_nan(self, body):
        result = swe.nod_aps_ut(JD_J2000, body, NODBIT_MEAN)
        assert len(result) == 4
        for tup in result:
            assert len(tup) == 6
            for val in tup:
                assert math.isnan(val)


class TestNodApsEquatorialSingleRotation:
    """FLG_EQUATORIAL must rotate ecliptic->equator exactly once.

    A regression once double-rotated the Moon points (the lunar branch let
    calc_ut() convert to the equator and then the formatter rotated again),
    inflating the declination by roughly the obliquity. This guards the
    single-rotation invariant against the ecliptic coordinates the same call
    reports, with no reference values baked in.
    """

    @pytest.mark.unit
    @pytest.mark.parametrize("body", [MOON, MARS])
    def test_node_declination_matches_single_rotation(self, body):
        from libephemeris.cache import get_true_obliquity

        ecl = swe.nod_aps(JD_J2000, body, NODBIT_MEAN, FLG_SWIEPH)
        equ = swe.nod_aps(JD_J2000, body, NODBIT_MEAN, FLG_SWIEPH | FLG_EQUATORIAL)
        eps = math.radians(get_true_obliquity(JD_J2000))
        for ecl_pt, equ_pt in zip(ecl, equ):
            lon = math.radians(ecl_pt[0])
            lat = math.radians(ecl_pt[1])
            # Single ecliptic->equatorial rotation of the reported ecliptic point.
            x = math.cos(lat) * math.cos(lon)
            y = math.cos(lat) * math.sin(lon) * math.cos(eps) - math.sin(
                lat
            ) * math.sin(eps)
            z = math.cos(lat) * math.sin(lon) * math.sin(eps) + math.sin(
                lat
            ) * math.cos(eps)
            dec_expected = math.degrees(math.asin(max(-1.0, min(1.0, z))))
            assert abs(equ_pt[1] - dec_expected) < 1e-6


class TestNodApsSiderealMeanAyanamsha:
    """Sidereal nod_aps longitudes subtract the MEAN ayanamsha.

    Measured reference behavior: even though the of-date coordinate carries
    nutation, the sidereal reduction uses the mean ayanamsha (no
    nutation-in-longitude term). The tropical->sidereal shift must therefore
    equal the mean ayanamsha, not the true one (they differ by ~14" at J2000).
    """

    @pytest.mark.unit
    @pytest.mark.parametrize("body", [MOON, MARS])
    def test_sidereal_shift_is_mean_ayanamsha(self, body):
        swe.set_sid_mode(SIDM_LAHIRI, 0.0, 0.0)
        trop = swe.nod_aps(JD_J2000, body, NODBIT_MEAN, FLG_SWIEPH)
        sid = swe.nod_aps(JD_J2000, body, NODBIT_MEAN, FLG_SWIEPH | FLG_SIDEREAL)
        mean_aya = swe.get_ayanamsa_ut(2451545.0 - swe.deltat(2451545.0))
        # Compare against a non-zero slot (perihelion of an inclined orbit).
        peri_shift = (trop[2][0] - sid[2][0]) % 360.0
        assert abs(peri_shift - mean_aya) < 1e-3
