# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Guards on the shared helpers every house construction is built out of.

None of the eleven is exported and nothing outside the module calls one, so
the behaviour they carry is observable only through the house entry points.
Three groups of checks live here.

*The statements that need no reference values.* The degree trigonometry is its
own inverse, the frame rotation is undone by the opposite tilt, the Ascendant
always follows the midheaven, feeding the reported Ascendant and midheaven
back to :func:`house_pos` returns the cusps they open, and the Placidus answer
is at least seven exactly while the body is above the horizon.

*The three configurations no recorded chart reaches.* A geographic pole with
the ARMC at exactly half a turn puts an equinox on the horizon, so neither
intersection of the ecliptic with the horizon is rising and the eastern test
is a tie; a body exactly on the horizon is the same tie for the above-horizon
indicator. Both are decided here -- the boundary is closed, the candidate is
kept, on the horizon counts as above -- and pinned so the choice stops being
unwitnessed. The third is the obliquity of ninety degrees, where the ecliptic
passes through the celestial poles and the midheaven does not exist.

*A large and a small obliquity.* The recorded house-position family carries
two hundred obliquities and every one of them is within a degree and a half of
23.4, so the obliquity axis of the rotation, the anchoring and the midheaven is
barely exercised by it. The same identities are asserted here at five and at
sixty degrees.
"""

from __future__ import annotations

import importlib
import math

import pytest

import libephemeris as ephem

houses_mod = importlib.import_module("libephemeris.houses")

_sin_deg = houses_mod._sin_deg
_cos_deg = houses_mod._cos_deg
_tan_deg = houses_mod._tan_deg
_atan_deg = houses_mod._atan_deg
_asin_deg = houses_mod._asin_deg
_acos_deg = houses_mod._acos_deg
_rotate_frame = houses_mod._rotate_frame
_ascendant_on_eastern_horizon = houses_mod._ascendant_on_eastern_horizon
_above_horizon_test = houses_mod._above_horizon_test
_asc_diff_saturated = houses_mod._asc_diff_saturated
_armc_to_mc = houses_mod._armc_to_mc
_calc_ascendant = houses_mod._calc_ascendant

#: The obliquity the recorded families carry, and one far below and one far
#: above it. Nothing in the recorded house-position family leaves 22.7-24.2.
OBLIQUITIES = (5.0, 23.4392911, 60.0)

#: Latitudes spanning both hemispheres, both polar circles and both poles.
LATITUDES = (-90.0, -80.0, -66.56, -45.0, -0.001, 0.0, 23.44, 66.56, 80.0, 90.0)


class TestDegreeTrigonometry:
    """The six one-argument wrappers."""

    @pytest.mark.unit
    def test_direct_functions_are_the_radian_ones(self):
        for x in (-720.0, -180.0, -90.0, -30.0, 0.0, 1e-9, 45.0, 90.0, 359.5):
            assert _sin_deg(x) == math.sin(math.radians(x))
            assert _cos_deg(x) == math.cos(math.radians(x))
            assert _tan_deg(x) == math.tan(math.radians(x))

    @pytest.mark.unit
    def test_inverse_functions_return_degrees_in_their_principal_range(self):
        for r in (-1.0, -0.5, 0.0, 0.25, 1.0):
            assert -90.0 <= _asin_deg(r) <= 90.0
            assert 0.0 <= _acos_deg(r) <= 180.0
            assert _asin_deg(r) == math.degrees(math.asin(r))
            assert _acos_deg(r) == math.degrees(math.acos(r))
        for r in (-1e9, -1.0, 0.0, 3.0, 1e9):
            assert -90.0 < _atan_deg(r) < 90.0
            assert _atan_deg(r) == math.degrees(math.atan(r))

    @pytest.mark.unit
    def test_the_inverse_functions_saturate_at_both_ends(self):
        """An argument carried past the domain answers at the nearer end.

        The ratios these helpers are fed are quotients of quantities equal in
        exact arithmetic, so an argument a rounding step outside ``[-1, 1]``
        means "at the boundary" and not "no answer". The two ends are treated
        alike.
        """
        one_up = math.nextafter(1.0, 2.0)
        one_down = math.nextafter(-1.0, -2.0)
        assert _asin_deg(one_up) == _asin_deg(1.0) == 90.0
        assert _asin_deg(one_down) == _asin_deg(-1.0) == -90.0
        assert _acos_deg(one_up) == _acos_deg(1.0) == 0.0
        assert _acos_deg(one_down) == _acos_deg(-1.0) == 180.0
        for far in (2.0, 1e6, math.inf):
            assert _asin_deg(far) == 90.0
            assert _asin_deg(-far) == -90.0
            assert _acos_deg(far) == 0.0
            assert _acos_deg(-far) == 180.0

    @pytest.mark.unit
    def test_the_arccosine_of_an_argument_one_step_above_unity(self):
        """The 30 recorded APC cells at the equator reach exactly this.

        Without the domain rule they raise a ``ValueError`` out of the public
        API instead of answering.
        """
        assert _acos_deg(1.0 + 2.220446049250313e-16) == 0.0
        assert ephem.house_pos(
            131.21184972525845,
            0.0,
            23.454715835385162,
            (312.21199529610425, 7.343736916801047e-05),
            b"Y",
        ) == pytest.approx(4.115542837916699, abs=1e-12)


class TestRotateFrame:
    """The rotation of a direction into a frame tilted about longitude zero."""

    @pytest.mark.unit
    def test_a_zero_tilt_is_the_identity(self):
        for lon, lat in ((0.0, 0.0), (137.5, -42.0), (359.9, 89.0)):
            out_lon, out_lat = _rotate_frame(lon, lat, 0.0)
            assert out_lon == pytest.approx(lon % 360.0, abs=1e-12)
            assert out_lat == pytest.approx(lat, abs=1e-12)

    @pytest.mark.unit
    @pytest.mark.parametrize("eps", OBLIQUITIES)
    def test_the_opposite_tilt_undoes_it(self, eps):
        for lon, lat in ((0.0, 0.0), (73.25, 12.5), (211.0, -35.0), (300.0, 5.0)):
            mid = _rotate_frame(lon, lat, eps)
            back = _rotate_frame(mid[0], mid[1], -eps)
            assert back[0] == pytest.approx(lon % 360.0, abs=1e-10)
            assert back[1] == pytest.approx(lat, abs=1e-10)

    @pytest.mark.unit
    @pytest.mark.parametrize("eps", OBLIQUITIES)
    def test_the_positive_tilt_carries_the_equator_onto_the_ecliptic(self, eps):
        """A tilt of ``+eps`` takes an equatorial place to an ecliptic one.

        The solstitial colure is the check that fixes the sign: the point at
        right ascension 90 and declination ``eps`` is ecliptic longitude 90 and
        latitude 0.
        """
        lon, lat = _rotate_frame(90.0, eps, eps)
        assert lon == pytest.approx(90.0, abs=1e-10)
        assert lat == pytest.approx(0.0, abs=1e-10)

    @pytest.mark.unit
    def test_the_longitude_is_normalised_and_the_latitude_bounded(self):
        for lon in (0.0, 95.0, 180.0, 275.0, 359.99):
            for lat in (-90.0, -45.0, 0.0, 45.0, 90.0):
                out_lon, out_lat = _rotate_frame(lon, lat, 37.5)
                assert 0.0 <= out_lon < 360.0
                assert -90.0 <= out_lat <= 90.0


class TestAscendantOnEasternHorizon:
    """Keeping the rising intersection rather than the setting one."""

    @pytest.mark.unit
    @pytest.mark.parametrize("eps", OBLIQUITIES)
    def test_it_never_fires_equatorward_of_the_polar_circle(self, eps):
        """The culminating degree is above the horizon there, always."""
        limit = 90.0 - eps
        for share in (-0.99, -0.6, 0.0, 0.25, 0.99):
            lat = share * limit
            for k in range(0, 360, 5):
                assert _ascendant_on_eastern_horizon(123.0, float(k), eps, lat) == 123.0

    @pytest.mark.unit
    @pytest.mark.parametrize("eps", OBLIQUITIES)
    def test_it_fires_only_poleward_of_it_and_by_half_a_turn(self, eps):
        fired = 0
        for lat in (-89.0, -70.0, 70.0, 89.0):
            if abs(lat) < 90.0 - eps:
                continue
            for k in range(0, 360, 5):
                out = _ascendant_on_eastern_horizon(123.0, float(k), eps, lat)
                assert out in (123.0, 303.0)
                fired += out == 303.0
        assert fired > 0

    @pytest.mark.unit
    @pytest.mark.parametrize("eps", OBLIQUITIES)
    @pytest.mark.parametrize("lat", LATITUDES)
    def test_the_ascendant_follows_the_midheaven(self, eps, lat):
        """``0 < (Asc - MC) mod 360 < 180`` -- the invariant the unit exists for.

        The two degenerate frames are excluded: at a geographic pole with an
        equinox on the meridian the Ascendant *is* the midheaven, and at zero
        obliquity there is nothing to anchor.
        """
        for k in range(120):
            armc = k * 3.0
            _, ascmc = ephem.houses_armc(armc, lat, eps, ord("W"))
            asc, mc = ascmc[0], ascmc[1]
            if abs(lat) == 90.0 and armc % 180.0 == 0.0:
                # An equinox on the meridian of a chart at a geographic pole:
                # the culminating degree is the rising one and the difference
                # the invariant measures has collapsed to nothing.
                continue
            assert 0.0 < (asc - mc) % 360.0 < 180.0

    @pytest.mark.unit
    @pytest.mark.parametrize("eps", OBLIQUITIES)
    def test_the_two_ascendant_paths_agree(self, eps):
        """The cusp side and the position side reach the same rising degree."""
        for lat in LATITUDES:
            for k in range(0, 360, 3):
                armc = float(k)
                if abs(lat) == 90.0 and (armc % 180.0 == 0.0 or eps == 0.0):
                    # The two configurations with no rising point at all: an
                    # equinox on the horizon, and the ecliptic lying along it.
                    continue
                cusp_asc = ephem.houses_armc(armc, lat, eps, ord("W"))[1][0]
                candidate = _calc_ascendant((armc + 90.0) % 360.0, eps, lat, lat)
                pos_asc = _ascendant_on_eastern_horizon(candidate, armc, eps, lat)
                gap = abs((cusp_asc - pos_asc + 180.0) % 360.0 - 180.0)
                assert gap < 1e-9

    @pytest.mark.unit
    def test_the_recorded_six_house_shift_at_high_latitude(self):
        """One frame, one body, three latitudes: the anchoring is the frame's.

        Without it the two polar answers land exactly six houses from the
        recorded ones.
        """
        body = (285.38339935017507, 0.00014133248984594553)
        eps = 23.45654434517177
        assert ephem.house_pos(
            105.21808180875615, 0.0, eps, body, b"E"
        ) == pytest.approx(3.9622100040763675, abs=1e-12)
        assert ephem.house_pos(
            285.21808180875615, 80.0, eps, body, b"E"
        ) == pytest.approx(10.873743298344568, abs=1e-12)
        assert ephem.house_pos(
            105.21808180875615, -80.0, eps, body, b"E"
        ) == pytest.approx(4.87374329834457, abs=1e-12)

    @pytest.mark.unit
    def test_the_boundary_is_closed_and_keeps_the_candidate(self):
        """The culminating degree exactly on the horizon counts as above it.

        ``cos(lat) + sin(lat) tan(eps) sin(armc)`` is exactly zero for this
        frame in binary64: the culminating degree is on the horizon and neither
        intersection of the ecliptic with the horizon is rising. The closed
        boundary keeps the candidate, which is the branch the cusp-side
        Ascendant takes, so the two paths answer alike rather than half a turn
        apart.
        """
        armc, eps, geolat = 270.0, 83.8, 6.2
        quantity = _cos_deg(geolat) + _sin_deg(geolat) * _tan_deg(eps) * _sin_deg(armc)
        assert quantity == 0.0
        assert _ascendant_on_eastern_horizon(41.0, armc, eps, geolat) == 41.0

    @pytest.mark.unit
    def test_the_tie_at_a_geographic_pole_with_the_armc_at_half_a_turn(self):
        """No recorded call reaches it: not one has an integer ARMC.

        The ARMC at exactly half a turn puts an equinox on the meridian and, at
        a geographic pole, on the horizon: the ecliptic cuts the horizon at the
        two equinoxes and neither is rising. The choice is the closed boundary
        of the test above -- the candidate is kept -- which gives 180 at the
        north pole and 0 at the south, the one-sided limit from just inside the
        pole and the branch the cusp path returns, so the two paths answer
        alike instead of half a turn apart.

        The tie is never reached exactly: ``cos(90)`` and ``sin(180)`` are two
        binary64 residues, and the test compares them rather than two zeros. It
        keeps the candidate for every obliquity a chart can carry, and takes
        the other branch at the south pole only past an obliquity of 26.565
        degrees, where ``tan(eps)`` passes one half. That is pinned too, so the
        two paths are never seen to part company without a witness.
        """
        # A zero obliquity is excluded: the ecliptic then lies along the
        # horizon of a polar chart and there is no rising point at all.
        for eps in (5.0, 23.4392911, 24.5, 26.0):
            for lat, expected in ((90.0, 180.0), (-90.0, 0.0)):
                candidate = _calc_ascendant(270.0, eps, lat, lat)
                assert (
                    _ascendant_on_eastern_horizon(candidate, 180.0, eps, lat)
                    == expected
                )
                cusp_asc = ephem.houses_armc(180.0, lat, eps, ord("W"))[1][0]
                assert abs((cusp_asc - expected + 180.0) % 360.0 - 180.0) < 1e-9

        # Past the half of the tangent the south-pole cell reads the residues
        # the other way round; the north-pole cell never does.
        assert _ascendant_on_eastern_horizon(0.0, 180.0, 60.0, -90.0) == 180.0
        assert _ascendant_on_eastern_horizon(180.0, 180.0, 60.0, 90.0) == 180.0

    @pytest.mark.unit
    def test_the_ascendant_at_a_pole_follows_the_sign_of_the_armc(self):
        """The ecliptic meets the horizon at the equinoxes there.

        The Ascendant can only be 0 or 180, and it is 180 while the sine of the
        ARMC is positive: the one-sided limit from just inside the pole.
        """
        eps = 23.4392911
        for lat in (90.0, -90.0):
            for armc in (30.0, 90.0, 150.0):
                assert ephem.houses_armc(armc, lat, eps, ord("W"))[1][0] == 180.0
            for armc in (210.0, 270.0, 330.0):
                assert ephem.houses_armc(armc, lat, eps, ord("W"))[1][0] == 0.0


class TestAboveHorizonTest:
    """The divided sign test, and where it stops being an altitude."""

    @pytest.mark.unit
    def test_its_sign_is_the_altitudes_away_from_the_poles(self):
        for lat in (-80.0, -45.0, -0.5, 0.0, 23.44, 66.0, 80.0):
            for dec in (-70.0, -23.44, 0.0, 15.0, 70.0):
                for md in (0.0, 37.0, 89.0, 91.0, 179.0, 180.0, 270.0):
                    altitude = _asin_deg(
                        _sin_deg(lat) * _sin_deg(dec)
                        + _cos_deg(lat) * _cos_deg(dec) * _cos_deg(md)
                    )
                    indicator = _above_horizon_test(dec, lat, md)
                    if abs(altitude) < 1e-9:
                        continue
                    assert (indicator >= 0.0) == (altitude >= 0.0)

    @pytest.mark.unit
    def test_it_is_bounded_for_a_circumpolar_body(self):
        """No arccosine is taken, so a body that never sets still answers."""
        for md in (0.0, 90.0, 180.0):
            value = _above_horizon_test(80.0, 80.0, md)
            assert math.isfinite(value)
            assert value > 0.0

    @pytest.mark.unit
    def test_the_body_exactly_on_the_horizon_counts_as_above(self):
        """The unwitnessed tie: no recorded call places a body on the horizon.

        The interval is closed on the horizon, which is what makes the two
        branches of a semi-arc construction meet without a gap at the rising
        and the setting points. The Placidus answer shows it: a body exactly at
        the setting point reads seven, the cusp it opens, and one exactly at
        the rising point reads one.
        """
        assert _above_horizon_test(83.0, 7.0, 180.0) == 0.0
        eps = 23.4392911
        assert ephem.house_pos(270.0, 0.0, eps, (0.0, 0.0), b"P") == 1.0
        assert ephem.house_pos(90.0, 0.0, eps, (0.0, 0.0), b"P") == 7.0

    @pytest.mark.unit
    @pytest.mark.parametrize("eps", OBLIQUITIES)
    def test_the_placidus_answer_carries_the_sign(self, eps):
        """At least seven above the horizon, below seven beneath it."""
        for lat in (-80.0, -45.0, 0.0, 45.0, 80.0):
            for k in range(0, 360, 9):
                for j in range(0, 360, 17):
                    armc, lon = float(k), float(j)
                    position = ephem.house_pos(armc, lat, eps, (lon, 0.0), b"P")
                    ra, dec = houses_mod._house_pos_equatorial(lon, 0.0, eps)
                    md = houses_mod._meridian_distance(ra, armc)
                    indicator = _above_horizon_test(dec, lat, md)
                    if abs(indicator) < 1e-12:
                        continue
                    assert (indicator >= 0.0) == (position >= 7.0)


class TestAscensionalDifference:
    """The arcsine of the product of the two tangents, and nothing else."""

    @pytest.mark.unit
    def test_it_is_the_published_relation(self):
        for lat in (-60.0, -23.44, 0.0, 12.0, 55.0):
            for dec in (-20.0, 0.0, 10.0, 23.44):
                assert _asc_diff_saturated(dec, lat) == pytest.approx(
                    math.degrees(
                        math.asin(
                            math.tan(math.radians(lat)) * math.tan(math.radians(dec))
                        )
                    ),
                    abs=1e-12,
                )

    @pytest.mark.unit
    def test_it_saturates_symmetrically(self):
        """A circumpolar body answers ``+90``, one that never rises ``-90``.

        The semi-diurnal arc is then a whole 180 degrees or nothing at all,
        which is the right geometry, and the two ends are alike.
        """
        assert _asc_diff_saturated(40.0, 80.0) == 90.0
        assert _asc_diff_saturated(-40.0, 80.0) == -90.0
        assert _asc_diff_saturated(40.0, -80.0) == -90.0
        assert _asc_diff_saturated(-40.0, -80.0) == 90.0

    @pytest.mark.unit
    def test_the_polar_value(self):
        """At a geographic pole every body off the equator is a right angle."""
        for dec in (-60.0, -1.0, 1.0, 60.0):
            assert abs(_asc_diff_saturated(dec, 90.0)) == 90.0
            assert abs(_asc_diff_saturated(dec, -90.0)) == 90.0
        assert _asc_diff_saturated(0.0, 90.0) == 0.0
        assert _asc_diff_saturated(0.0, -90.0) == 0.0

    @pytest.mark.unit
    def test_a_body_on_the_equator_or_an_observer_on_it(self):
        for lat in (-70.0, 0.0, 70.0):
            assert _asc_diff_saturated(0.0, lat) == 0.0
        for dec in (-70.0, 0.0, 70.0):
            assert _asc_diff_saturated(dec, 0.0) == 0.0


class TestArmcToMc:
    """The degree of the ecliptic whose right ascension is the ARMC."""

    @pytest.mark.unit
    def test_at_zero_obliquity_the_midheaven_is_the_armc(self):
        """Exactly, including on the seam of the circle."""
        for k in range(0, 360, 5):
            assert _armc_to_mc(float(k), 0.0) == pytest.approx(float(k), abs=1e-12)
        assert _armc_to_mc(330.0, 0.0) == 330.0
        assert _armc_to_mc(0.0, 0.0) == 0.0
        assert ephem.houses_armc(330.0, 0.0, 0.0, ord("O"))[1][1] == 330.0

    @pytest.mark.unit
    @pytest.mark.parametrize("eps", OBLIQUITIES)
    def test_it_stays_in_the_quadrant_of_the_armc(self, eps):
        """Below a right angle the two never part company on a cardinal point."""
        for k in range(1, 360):
            if k % 90 == 0:
                continue
            armc = float(k)
            assert int(_armc_to_mc(armc, eps) // 90.0) == int(armc // 90.0)

    @pytest.mark.unit
    @pytest.mark.parametrize("eps", OBLIQUITIES)
    def test_the_cardinal_points_are_shared(self, eps):
        for armc, expected in (
            (0.0, 0.0),
            (90.0, 90.0),
            (180.0, 180.0),
            (270.0, 270.0),
        ):
            assert _armc_to_mc(armc, eps) == pytest.approx(expected, abs=1e-9)

    @pytest.mark.unit
    def test_it_never_sees_the_latitude(self):
        for lat in LATITUDES:
            assert (
                ephem.houses_armc(75.0, lat, 23.4392911, ord("W"))[1][1]
                == ephem.houses_armc(75.0, 0.0, 23.4392911, ord("W"))[1][1]
            )

    @pytest.mark.unit
    def test_the_recorded_midheavens(self):
        for armc, lat, eps, expected in (
            (30.0, -45.0, 23.4392911, 32.18125916628427),
            (120.0, -45.0, 23.4392911, 117.91054978805977),
            (225.0, 45.0, 22.0, 227.16381750390525),
            (315.0, -45.0, 24.5, 312.3009805063376),
            (89.999, 45.0, 23.4392911, 89.99908251793784),
            (269.999, 45.0, 23.4392911, 269.99908251793784),
            (359.999, 45.0, 23.4392911, 359.9989100604347),
        ):
            got = ephem.houses_armc(armc, lat, eps, ord("W"))[1][1]
            assert got == pytest.approx(expected, abs=1e-12)

    @pytest.mark.unit
    @pytest.mark.parametrize("eps", OBLIQUITIES)
    def test_feeding_the_midheaven_back_returns_the_tenth_cusp(self, eps):
        """Exactly, on the system whose houses begin a quarter turn before it."""
        for lat in LATITUDES:
            for k in range(0, 360, 9):
                armc = float(k)
                mc = ephem.houses_armc(armc, lat, eps, ord("D"))[1][1]
                assert ephem.house_pos(armc, lat, eps, (mc, 0.0), b"D") == 10.0

    @pytest.mark.unit
    def test_the_degenerate_obliquity_answers_a_celestial_pole(self):
        """At ninety degrees the ecliptic passes through the celestial poles.

        The meridian then meets it only at the poles themselves, so the
        midheaven does not exist: the two candidates are ecliptic longitude 90,
        the north pole, and 270, the south. The two-component relation the unit
        evaluates everywhere else resolves to the one on the side of the
        ascending half of the ecliptic -- 90 while the sine of the ARMC is
        positive, 270 while it is negative -- and nothing recorded there is a
        value of any rule.
        """
        for armc in (30.0, 105.0, 150.0):
            assert _armc_to_mc(armc, 90.0) == pytest.approx(90.0, abs=1e-9)
        for armc in (210.0, 285.0, 330.0):
            assert _armc_to_mc(armc, 90.0) == pytest.approx(270.0, abs=1e-9)


class TestRoundTripsThroughThePublicApi:
    """What the reported angles answer when fed back to :func:`house_pos`."""

    @pytest.mark.unit
    @pytest.mark.parametrize("eps", OBLIQUITIES)
    def test_the_ascendant_opens_the_first_house(self, eps):
        for lat in LATITUDES:
            for k in range(0, 360, 9):
                armc = float(k)
                if abs(lat) == 90.0 and (armc % 180.0 == 0.0 or eps == 0.0):
                    continue
                asc = ephem.houses_armc(armc, lat, eps, ord("A"))[1][0]
                for selector, expected in (
                    (b"A", 1.0),
                    (b"V", 1.5),
                    (b"W", 1.0 + (asc % 30.0) / 30.0),
                ):
                    if selector == b"W" and min(asc % 30.0, 30.0 - asc % 30.0) < 1e-9:
                        # The two Ascendant paths differ in the last place; on a
                        # sign boundary that step is a whole house of the
                        # whole-sign scale.
                        continue
                    got = ephem.house_pos(armc, lat, eps, (asc, 0.0), selector)
                    if got > 12.5:
                        # The same last-place step, with the body a hair before
                        # the first cusp, at the far end of the half-open
                        # interval.
                        got -= 12.0
                    assert got == pytest.approx(expected, abs=1e-9)
