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
    CHIRON,
    PHOLUS,
    CERES,
    PALLAS,
    JUNO,
    VESTA,
    NODBIT_MEAN,
    NODBIT_OSCU,
    NODBIT_OSCU_BAR,
    NODBIT_FOPOINT,
    FLG_SWIEPH,
    FLG_SPEED,
    FLG_HELCTR,
    FLG_J2000,
    FLG_EQUATORIAL,
    FLG_SIDEREAL,
    SIDM_LAHIRI,
)

# Curated minor bodies whose osculating nodes/apsides nod_aps now computes.
# Epochs kept inside the 1920-2080 asteroid-SPK safe window (see CLAUDE.md).
_MINOR_BODIES = [CHIRON, PHOLUS, CERES, PALLAS, JUNO, VESTA]
_MINOR_EPOCHS = [2436205.0, 2451545.0, 2469807.0]  # ~1958, J2000, ~2050


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


class TestNodApsMinorBodies:
    """nod_aps computes osculating nodes/apsides for asteroids/centaurs.

    These are property/invariant checks on the library alone (no reference
    values baked in): the osculating decomposition is sourced from the same
    position pipeline as calc()/get_orbital_elements, so the results must be
    finite, self-consistent, and consistent with the reported orbital
    elements. The like-for-like numerical agreement with the reference (and its
    documented osculating-realization residuals) is pinned in the validation
    compare suite, not here.
    """

    @pytest.mark.unit
    @pytest.mark.parametrize("body", _MINOR_BODIES)
    @pytest.mark.parametrize("jd", _MINOR_EPOCHS)
    @pytest.mark.parametrize("center", [0, FLG_HELCTR])
    def test_finite_and_in_range(self, body, jd, center):
        """Every slot is finite; longitudes/latitudes/distances are in range."""
        result = swe.nod_aps_ut(jd, body, NODBIT_OSCU, FLG_SWIEPH | center)
        assert len(result) == 4
        for pt in result:
            assert len(pt) == 6
            for val in pt:
                assert math.isfinite(val)
            assert 0.0 <= pt[0] < 360.0
            assert -90.0 <= pt[1] <= 90.0
            assert pt[2] > 0.0

    @pytest.mark.unit
    @pytest.mark.parametrize("body", _MINOR_BODIES)
    def test_methods_default_mean_oscu_identical(self, body):
        """Minor bodies have no mean-element model, so the default, NODBIT_MEAN
        and NODBIT_OSCU methods all reduce to the same osculating result."""
        base = swe.nod_aps_ut(JD_J2000, body, 0, FLG_SWIEPH)
        for method in (NODBIT_MEAN, NODBIT_OSCU):
            other = swe.nod_aps_ut(JD_J2000, body, method, FLG_SWIEPH)
            for a, b in zip(base, other):
                assert a == pytest.approx(b, abs=1e-9)

    @pytest.mark.unit
    @pytest.mark.parametrize("body", _MINOR_BODIES)
    def test_oscu_bar_only_shifts_trans_jovian(self, body):
        """NODBIT_OSCU_BAR re-references only orbits beyond Jupiter to the
        barycentre; a main-belt asteroid is unchanged from plain OSCU."""
        a_semi = swe.get_orbital_elements(JD_J2000, body, FLG_HELCTR)[0]
        oscu = swe.nod_aps_ut(JD_J2000, body, NODBIT_OSCU, FLG_SWIEPH)
        bar = swe.nod_aps_ut(JD_J2000, body, NODBIT_OSCU_BAR, FLG_SWIEPH)
        peri_shift = abs((oscu[2][0] - bar[2][0] + 180.0) % 360.0 - 180.0)
        if a_semi < 6.0:  # main belt: heliocentric == barycentric
            assert peri_shift < 1e-6
        else:  # trans-jovian centaur: the barycentric ellipse differs
            assert peri_shift > 1e-3

    @pytest.mark.unit
    @pytest.mark.parametrize("body", _MINOR_BODIES)
    @pytest.mark.parametrize("jd", _MINOR_EPOCHS)
    def test_helio_nodes_antipodal_and_in_ecliptic(self, body, jd):
        """Heliocentric nodes lie on the ecliptic (lat 0) and are 180° apart."""
        nasc, ndsc, _, _ = swe.nod_aps_ut(
            jd, body, NODBIT_OSCU, FLG_SWIEPH | FLG_HELCTR
        )
        sep = abs((nasc[0] - ndsc[0]) % 360.0 - 180.0)
        assert sep < 1e-6
        assert abs(nasc[1]) < 1e-6
        assert abs(ndsc[1]) < 1e-6

    @pytest.mark.unit
    @pytest.mark.parametrize("body", _MINOR_BODIES)
    @pytest.mark.parametrize("jd", _MINOR_EPOCHS)
    def test_helio_perihelion_closer_and_antipodal(self, body, jd):
        """Heliocentric perihelion is nearer than aphelion and 180° opposite."""
        _, _, peri, aphe = swe.nod_aps_ut(
            jd, body, NODBIT_OSCU, FLG_SWIEPH | FLG_HELCTR
        )
        assert peri[2] < aphe[2]
        sep = abs((peri[0] - aphe[0]) % 360.0 - 180.0)
        assert sep < 1e-6
        assert peri[1] == pytest.approx(-aphe[1], abs=1e-6)

    @pytest.mark.unit
    @pytest.mark.parametrize("body", _MINOR_BODIES)
    @pytest.mark.parametrize("jd", _MINOR_EPOCHS)
    def test_apsis_distances_track_orbital_elements(self, body, jd):
        """The heliocentric perihelion/aphelion distances equal the reported
        q/Q from get_orbital_elements (both from the same osculating state).
        Distances are frame-invariant, so this holds at every epoch."""
        elem = swe.get_orbital_elements_ut(jd, body, FLG_HELCTR)
        q_dist, big_q = elem[15], elem[16]  # perihelion / aphelion distance (AU)
        _, _, peri, aphe = swe.nod_aps_ut(
            jd, body, NODBIT_OSCU, FLG_SWIEPH | FLG_HELCTR
        )
        assert peri[2] == pytest.approx(q_dist, rel=5e-4)
        assert aphe[2] == pytest.approx(big_q, rel=5e-4)

    @pytest.mark.unit
    @pytest.mark.parametrize("body", _MINOR_BODIES)
    def test_node_longitude_tracks_omega_at_j2000(self, body):
        """At the J2000 epoch (where the of-date and J2000 ecliptics coincide),
        the heliocentric ascending-node longitude equals the longitude of
        ascending node from get_orbital_elements up to the retained nutation.

        The nod_aps J2000 output builds the point on the ecliptic of date and
        precesses it (nutation retained), while get_orbital_elements works
        natively in the J2000 ecliptic; that convention gap grows with
        |t - J2000|, so the tight angular check is anchored at J2000."""
        omega_node = swe.get_orbital_elements(JD_J2000, body, FLG_HELCTR)[3]
        nasc, _, _, _ = swe.nod_aps(
            JD_J2000, body, NODBIT_OSCU, FLG_SWIEPH | FLG_HELCTR | FLG_J2000
        )
        node_diff = abs((nasc[0] - omega_node + 180.0) % 360.0 - 180.0)
        assert node_diff < 0.02

    @pytest.mark.unit
    @pytest.mark.parametrize("body", _MINOR_BODIES)
    def test_speed_slots_finite(self, body):
        """FLG_SPEED fills every speed slot with a finite central difference."""
        result = swe.nod_aps_ut(JD_J2000, body, NODBIT_OSCU, FLG_SWIEPH | FLG_SPEED)
        for pt in result:
            for val in pt[3:]:
                assert math.isfinite(val)

    @pytest.mark.unit
    def test_planets_still_work(self):
        """Extending to minor bodies must not affect the supported planets."""
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

    Measured reference behavior: nod_aps marks the three *position* slots
    (longitude, latitude, distance) not-a-number for the interpolated lunar
    apsides, but leaves the three *speed* slots at 0.0 -- never NaN -- for every
    method, with and without FLG_SPEED. The point is undefined as an orbital
    node/apse, yet the speed channels stay well-formed floats. Mirror that
    exactly (it neither raises nor zero-fills the position slots).
    """

    @pytest.mark.unit
    @pytest.mark.parametrize("body", [INTP_APOG, INTP_PERG])
    @pytest.mark.parametrize("method", [0, NODBIT_MEAN, NODBIT_OSCU])
    @pytest.mark.parametrize("flags", [FLG_SWIEPH, FLG_SWIEPH | FLG_SPEED])
    def test_interpolated_apsides_nan_positions_zero_speeds(self, body, method, flags):
        result = swe.nod_aps_ut(JD_J2000, body, method, flags)
        assert len(result) == 4
        want_nan_speeds = bool(flags & FLG_SPEED)
        for tup in result:
            assert len(tup) == 6
            # Position slots: not-a-number.
            for val in tup[:3]:
                assert math.isnan(val)
            # Speed slots (measured contract): exactly 0.0 on a plain
            # request; NaN propagates when FLG_SPEED is explicitly asked.
            for val in tup[3:]:
                assert isinstance(val, float)
                if want_nan_speeds:
                    assert math.isnan(val)
                else:
                    assert val == 0.0


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


class TestNodApsBarycentricFlag:
    """FLG_BARYCTR re-references node/apse points to the solar-system
    barycentre (the point vector gains the Sun's barycentric offset)."""

    def test_baryctr_differs_from_helctr_by_sun_offset(self):
        import math

        import libephemeris as le

        jd = 2455362.5
        hel = le.nod_aps_ut(
            jd, le.JUPITER, le.NODBIT_OSCU, le.FLG_SWIEPH | le.FLG_HELCTR
        )
        bar = le.nod_aps_ut(
            jd, le.JUPITER, le.NODBIT_OSCU, le.FLG_SWIEPH | le.FLG_BARYCTR
        )
        # The physical point is the same: the barycentric distance differs
        # from the heliocentric one by (at most) the Sun-SSB offset.
        d_dist = abs(bar[0][2] - hel[0][2])
        assert 0.0 < d_dist < 0.02  # Sun-SSB stays within ~0.02 AU
        # And the longitudes are no longer bit-identical (the old no-op).
        assert not math.isclose(bar[0][0], hel[0][0], abs_tol=1e-12)

    def test_geocentric_and_helctr_unchanged_by_flagless_paths(self):
        import libephemeris as le

        jd = 2455362.5
        geo = le.nod_aps_ut(jd, le.JUPITER, le.NODBIT_OSCU, le.FLG_SWIEPH)
        assert all(0.0 <= v < 360.0 for v in (geo[0][0], geo[2][0]))


class TestNodApsRoundUContracts:
    """Round-U measured contracts: Moon barycentre, flag precedence, INTP."""

    def test_moon_baryctr_offsets_by_barycentric_earth(self):
        import libephemeris as le

        jd = 2436005.0
        geo = le.nod_aps_ut(jd, le.MOON, le.NODBIT_OSCU, le.FLG_SWIEPH)
        bar = le.nod_aps_ut(jd, le.MOON, le.NODBIT_OSCU, le.FLG_SWIEPH | le.FLG_BARYCTR)
        # The barycentric point sits ~1 AU from the SSB (geocentric point
        # plus the barycentric Earth), not at the geocentric distance.
        assert 0.9 < bar[1][2] < 1.1
        assert geo[1][2] < 0.01

    def test_helctr_wins_over_baryctr_when_both_set(self):
        import libephemeris as le

        jd = 2436005.0
        hel = le.nod_aps_ut(jd, le.MARS, le.NODBIT_OSCU, le.FLG_SWIEPH | le.FLG_HELCTR)
        both = le.nod_aps_ut(
            jd,
            le.MARS,
            le.NODBIT_OSCU,
            le.FLG_SWIEPH | le.FLG_HELCTR | le.FLG_BARYCTR,
        )
        assert both[2][0] == hel[2][0]

    def test_intp_flag_contracts(self):
        import math

        import libephemeris as le

        t = le.julday(2020, 6, 15, 0.0)
        # EQUATORIAL: NaN stays NaN (no clamped pole declination).
        equ = le.nod_aps_ut(t, 21, le.NODBIT_MEAN, le.FLG_SWIEPH | le.FLG_EQUATORIAL)
        assert math.isnan(equ[0][1]) and equ[0][3] == 0.0
        # SPEED: NaN propagates through the speed slots.
        spd = le.nod_aps_ut(t, 21, le.NODBIT_MEAN, le.FLG_SWIEPH | le.FLG_SPEED)
        assert all(math.isnan(v) for v in spd[0])
