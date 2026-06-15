# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage completion tests for ``libephemeris.hypothetical``.

These tests exercise the orbital-elements parsers, the Keplerian propagation
helpers, the equinox-precession path, the velocity wrap-around branches, and
the dispatch table in :mod:`libephemeris.hypothetical`. They target the
specific lines/branches that the broad suite leaves uncovered, using
deterministic Julian Days and synthetic orbital elements so that no network
access or ephemeris file is required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from libephemeris import hypothetical as H
from libephemeris.hypothetical import (
    CUPIDO,
    HARRINGTON,
    ISIS,
    NIBIRU,
    PROSERPINA,
    PLANET_X_LOWELL,
    PLANET_X_PICKERING,
    VULCAN,
    WALDEMATH,
    WHITE_MOON,
    OrbitalElements,
    TPolynomial,
    UranianKeplerianElements,
    _calc_keplerian_position,
    _calc_keplerian_position_raw,
    _calc_orbital_position_raw,
    _calc_transpluto_raw,
    _calc_uranian_planet_raw,
    _calc_vulcan_raw,
    _keplerian_to_ecliptic_j2000,
    _parse_fictitious_orbits_csv,
    _parse_orbital_elements_line,
    _parse_t_polynomial,
    _solve_kepler_equation,
    calc_fictitious_position,
    calc_hypothetical_position,
    calc_orbital_position,
    calc_planet_x_lowell,
    calc_planet_x_pickering,
    calc_proserpina,
    calc_transpluto,
    calc_transpluto_position,
    calc_uranian_planet,
    calc_vulcan,
    calc_waldemath_position,
    calc_white_moon_position,
    get_bundled_fictitious_orbits_path,
    get_hypothetical_name,
    is_hypothetical_body,
    list_hypothetical_bodies,
    parse_orbital_elements,
)

J2000 = 2451545.0


# ---------------------------------------------------------------------------
# _parse_t_polynomial: pattern-2 minus and fallback-split branches
# ---------------------------------------------------------------------------


class TestParseTPolynomial:
    """Exercise the regex pattern-2 and the fallback string-split path."""

    def test_pattern2_rate_times_t_minus_constant(self):
        """``rate * T - const`` -> negative constant (lines 1490-1496)."""
        poly = _parse_t_polynomial("707550.7341 * T - 252.0")
        assert poly.linear == pytest.approx(707550.7341)
        assert poly.constant == pytest.approx(-252.0)

    def test_pattern2_rate_times_t_plus_constant(self):
        """``rate * T + const`` -> positive constant."""
        poly = _parse_t_polynomial("3.0 * T + 9.0")
        assert poly.linear == pytest.approx(3.0)
        assert poly.constant == pytest.approx(9.0)

    def test_fallback_split_constant_minus_t(self):
        """``const - T`` falls through regex into split path (lines 1517, 1522)."""
        poly = _parse_t_polynomial("5.0 - T")
        assert poly.constant == pytest.approx(5.0)
        assert poly.linear == pytest.approx(-1.0)

    def test_fallback_split_bare_minus_t(self):
        """``-T`` (no constant) goes through split with negative sign."""
        poly = _parse_t_polynomial("-T")
        assert poly.constant == pytest.approx(0.0)
        assert poly.linear == pytest.approx(-1.0)

    def test_fallback_split_constant_plus_coeff_t(self):
        """``const - coeff T`` (no '*') hits split with explicit coefficient."""
        poly = _parse_t_polynomial("5.0 - 2.0 T")
        assert poly.constant == pytest.approx(5.0)
        assert poly.linear == pytest.approx(-2.0)

    def test_simple_number(self):
        """Plain numbers skip the regex entirely."""
        poly = _parse_t_polynomial("42.5")
        assert poly.constant == pytest.approx(42.5)
        assert poly.linear == 0.0


# ---------------------------------------------------------------------------
# parse_orbital_elements: the "parsed is falsy" loop-back arc (1625->1608)
# ---------------------------------------------------------------------------


class TestParseOrbitalElementsLoopBack:
    """Cover the branch where the per-line parser yields a falsy result."""

    def test_parsed_none_skips_append(self, tmp_path, monkeypatch):
        """When the line parser returns None the loop continues (1625->1608)."""
        f = tmp_path / "orbits.txt"
        f.write_text(
            "J2000, J2000, 0.0, 100.0, 0.1, 45.0, 30.0, 5.0, BodyA\n"
            "J2000, J2000, 0.0, 100.0, 0.1, 45.0, 30.0, 5.0, BodyB\n",
            encoding="utf-8",
        )
        monkeypatch.setattr(H, "_parse_orbital_elements_line", lambda line, n: None)
        result = parse_orbital_elements(f)
        assert result == []

    def test_geocentric_space_geo_suffix(self):
        """``Foobar geo`` (space, no comma) -> geocentric (lines 1912-1913)."""
        line = "J2000, J2000, 0.0, 100.0, 0.1, 45.0, 30.0, 5.0, Foobar geo"
        elem = _parse_orbital_elements_line(line, 1)
        assert elem is not None
        assert elem.is_geocentric is True
        assert elem.name == "Foobar"


# ---------------------------------------------------------------------------
# CSV parser error paths
# ---------------------------------------------------------------------------


class TestFictitiousOrbitsCsvErrors:
    """Drive the CSV parser's inline-comment and validation error branches."""

    def test_bundled_path_exists(self):
        """The packaged dataset path resolves to an existing file."""
        path = get_bundled_fictitious_orbits_path()
        assert path.exists()

    def test_bundled_path_missing_raises(self, monkeypatch):
        """A missing bundled dataset raises FileNotFoundError (line 1657)."""
        real_exists = Path.exists

        def fake_exists(self):
            if self.name == "fictitious_orbits.csv":
                return False
            return real_exists(self)

        monkeypatch.setattr(Path, "exists", fake_exists)
        with pytest.raises(FileNotFoundError):
            get_bundled_fictitious_orbits_path()

    def test_csv_missing_file_raises(self, tmp_path):
        """A nonexistent CSV path raises FileNotFoundError (line 1691)."""
        with pytest.raises(FileNotFoundError):
            _parse_fictitious_orbits_csv(tmp_path / "nope.csv")

    def test_csv_inline_comment_and_blank(self, tmp_path):
        """Inline ``#`` comments are stripped; blank/comment lines skipped."""
        f = tmp_path / "orbits.csv"
        f.write_text(
            "# header comment\n"
            "\n"
            "Body, J2000, J2000, 0.0, 100.0, 0.1, 45.0, 30.0, 5.0, 0  # trailing\n",
            encoding="utf-8",
        )
        elems = _parse_fictitious_orbits_csv(f)
        assert len(elems) == 1
        assert elems[0].name == "Body"
        assert elems[0].is_geocentric is False

    def test_csv_too_few_columns(self, tmp_path):
        """Fewer than 10 columns raises a ValueError (lines 1707-1711)."""
        f = tmp_path / "short.csv"
        f.write_text("Body, J2000, J2000, 0.0\n", encoding="utf-8")
        with pytest.raises(ValueError, match="expected at least 10"):
            _parse_fictitious_orbits_csv(f)

    def test_csv_epoch_jdate_rejected(self, tmp_path):
        """A JDATE epoch is rejected with a wrapped ValueError (1727, 1748-1749)."""
        f = tmp_path / "jdate.csv"
        f.write_text(
            "Body, JDATE, J2000, 0.0, 100.0, 0.1, 45.0, 30.0, 5.0, 0\n",
            encoding="utf-8",
        )
        with pytest.raises(ValueError, match="Error parsing line"):
            _parse_fictitious_orbits_csv(f)

    def test_csv_bad_semi_axis(self, tmp_path):
        """A non-numeric semi-major axis raises a ValueError (1734-1735)."""
        f = tmp_path / "badaxis.csv"
        f.write_text(
            "Body, J2000, J2000, 0.0, notnum, 0.1, 45.0, 30.0, 5.0, 0\n",
            encoding="utf-8",
        )
        with pytest.raises(ValueError, match="semi-major axis"):
            _parse_fictitious_orbits_csv(f)

    def test_csv_bad_geocentric_flag(self, tmp_path):
        """A non-integer geocentric flag raises a ValueError (1743-1747)."""
        f = tmp_path / "badgeo.csv"
        f.write_text(
            "Body, J2000, J2000, 0.0, 100.0, 0.1, 45.0, 30.0, 5.0, maybe\n",
            encoding="utf-8",
        )
        with pytest.raises(ValueError, match="geocentric flag"):
            _parse_fictitious_orbits_csv(f)


# ---------------------------------------------------------------------------
# calc_orbital_position: velocity wrap branches (line 2111 and 2113)
# ---------------------------------------------------------------------------


class TestCalcOrbitalPositionWrap:
    """Exercise both longitude wrap-around branches in calc_orbital_position."""

    @staticmethod
    def _fast_elem(semi_axis):
        return OrbitalElements(
            name="Fast",
            epoch_jd=J2000,
            equinox_jd=J2000,
            equinox_is_jdate=False,
            mean_anomaly=TPolynomial(0.0, 0.0),
            semi_axis=semi_axis,
            eccentricity=TPolynomial(0.0, 0.0),
            arg_perihelion=TPolynomial(0.0, 0.0),
            asc_node=TPolynomial(0.0, 0.0),
            inclination=TPolynomial(0.0, 0.0),
        )

    def test_dlon_positive_wrap(self):
        """Fast prograde body whose central diff exceeds +90 (line 2111)."""
        elem = self._fast_elem(0.013)
        result = calc_orbital_position(elem, J2000)
        # After wrap correction dlon must be within (-180, 180).
        assert -180.0 < result[3] < 180.0

    def test_dlon_negative_wrap(self):
        """Fast body whose central diff drops below -90 (line 2113)."""
        elem = self._fast_elem(0.01)
        result = calc_orbital_position(elem, 2451552.77)
        assert -180.0 < result[3] < 180.0

    def test_raw_helper_positive_and_negative(self):
        """Confirm the raw helper produces the expected wrap inputs."""
        elem = self._fast_elem(0.013)
        prev = _calc_orbital_position_raw(elem, J2000 - 1.0)[0]
        nxt = _calc_orbital_position_raw(elem, J2000 + 1.0)[0]
        assert (nxt - prev) / 2.0 > 90.0


# ---------------------------------------------------------------------------
# is_hypothetical_body / get_hypothetical_name / list_hypothetical_bodies
# ---------------------------------------------------------------------------


class TestSimpleLookups:
    def test_is_hypothetical_body_true_and_false(self):
        """Range check helper (line 2214)."""
        assert is_hypothetical_body(CUPIDO) is True
        assert is_hypothetical_body(0) is False

    def test_get_hypothetical_name_unknown_default(self):
        """Unknown id falls back to the ``Unknown (..)`` default (line 2232)."""
        assert get_hypothetical_name(CUPIDO) == "Cupido"
        assert get_hypothetical_name(99999) == "Unknown (99999)"

    def test_list_hypothetical_bodies_copy(self):
        """The listing returns a copy of the name map (line 3809)."""
        bodies = list_hypothetical_bodies()
        assert bodies[CUPIDO] == "Cupido"
        bodies[CUPIDO] = "mutated"
        assert list_hypothetical_bodies()[CUPIDO] == "Cupido"


# ---------------------------------------------------------------------------
# calc_fictitious_position: unknown body + equinox precession paths
# ---------------------------------------------------------------------------


class TestFictitiousPosition:
    def test_unknown_body_raises(self):
        """Missing elements raise UnknownBodyError (lines 2335-2337)."""
        from libephemeris.exceptions import UnknownBodyError

        with pytest.raises(UnknownBodyError):
            calc_fictitious_position(99999, J2000)

    def test_harrington_uses_j2000_equinox_branch(self):
        """Harrington's equinox is J2000 -> else branch (lines 2916-2918)."""
        lon, lat, dist = (calc_fictitious_position(HARRINGTON, J2000))[:3]
        assert 0.0 <= lon < 360.0
        assert dist > 0.0

    def test_nibiru_uses_precession_branch(self):
        """Nibiru's equinox differs from J2000 -> precession path executes."""
        lon, lat, dist = (calc_fictitious_position(NIBIRU, J2000))[:3]
        assert 0.0 <= lon < 360.0
        assert dist > 0.0


# ---------------------------------------------------------------------------
# _keplerian_to_ecliptic_j2000 zero-distance latitude branch (line 2927)
# ---------------------------------------------------------------------------


class TestKeplerianZeroDistance:
    def test_zero_semi_axis_latitude_else(self):
        """a=0 gives r=0 so latitude falls into the else branch (line 2927)."""
        elem = UranianKeplerianElements(
            name="Zero",
            epoch=J2000,
            equinox_jd=J2000,
            a=0.0,
            e=0.0,
            i=0.0,
            omega=0.0,
            Omega=0.0,
            M0=0.0,
            n=0.0,
        )
        lon, lat, r = _keplerian_to_ecliptic_j2000(elem, J2000)
        assert r == 0.0
        assert lat == 0.0


# ---------------------------------------------------------------------------
# calc_uranian_planet + raw helper (lines 2989-2995, 3014-3015)
# ---------------------------------------------------------------------------


class TestUranianPlanet:
    def test_calc_uranian_planet_runs(self):
        """A full Uranian calculation exercises the velocity block."""
        lon, lat, dist, dlon, dlat, ddist = calc_uranian_planet(CUPIDO, J2000)
        assert 0.0 <= lon < 360.0
        assert dist > 0.0
        assert -180.0 < dlon < 180.0

    def test_raw_helper(self):
        """The raw helper returns a 3-tuple (lines 3014-3015)."""
        out = _calc_uranian_planet_raw(CUPIDO, J2000)
        assert len(out) == 3
        assert 0.0 <= out[0] < 360.0


# ---------------------------------------------------------------------------
# calc_transpluto negative wrap (line 2743)
# ---------------------------------------------------------------------------


class TestTranspluto:
    def test_calc_transpluto_negative_wrap(self):
        """At the 0/360 crossing the central diff wraps negative (line 2743)."""
        # jd located where Transpluto crosses the longitude origin.
        result = calc_transpluto(2368022.0)
        assert -180.0 < result[3] < 180.0

    def test_raw_helper(self):
        """The raw Transpluto helper returns a plausible longitude."""
        lon, lat, dist = _calc_transpluto_raw(J2000)
        assert 0.0 <= lon < 360.0
        assert dist > 0.0

    def test_calc_transpluto_position_delegates(self):
        """calc_transpluto_position routes through the Keplerian path (3229)."""
        result = calc_transpluto_position(J2000)
        assert 0.0 <= result[0] < 360.0


# ---------------------------------------------------------------------------
# _solve_kepler_equation non-convergence return (line 3258)
# ---------------------------------------------------------------------------


class TestSolveKepler:
    def test_non_convergence_returns_last_estimate(self):
        """tol=0 forces the loop to exhaust and return E (line 3258)."""
        e = _solve_kepler_equation(0.5, 0.5, tol=0.0)
        assert isinstance(e, float)

    def test_high_eccentricity_initial_guess(self):
        """e>=0.8 selects the pi initial guess and still converges."""
        e = _solve_kepler_equation(0.5, 0.9)
        assert isinstance(e, float)


# ---------------------------------------------------------------------------
# _calc_keplerian_position and its raw helper (3276-3340, 3347-3392)
# ---------------------------------------------------------------------------


class TestKeplerianPosition:
    def test_calc_keplerian_position_isis(self):
        """The eccentric ISIS orbit drives the full Keplerian block."""
        out = _calc_keplerian_position(ISIS, J2000)
        assert len(out) == 6
        assert 0.0 <= out[0] < 360.0
        assert out[2] > 0.0

    def test_calc_keplerian_position_negative_wrap(self):
        """At the ISIS 0/360 crossing the per-second diff wraps (line 3335)."""
        # jd located precisely at the ISIS longitude origin crossing so the
        # 1-second central difference straddles 360 -> 0.
        out = _calc_keplerian_position(ISIS, 2117498.2575203027)
        assert -180.0 < out[3] < 180.0

    def test_calc_keplerian_position_unknown_raises(self):
        """An id with no Keplerian elements raises ValueError (line 3277)."""
        with pytest.raises(ValueError, match="Keplerian elements"):
            _calc_keplerian_position(99999, J2000)

    def test_calc_keplerian_position_raw_isis(self):
        """The raw helper returns a 3-tuple for ISIS (lines 3347-3392)."""
        out = _calc_keplerian_position_raw(ISIS, J2000)
        assert len(out) == 3
        assert out[2] > 0.0

    def test_calc_keplerian_position_raw_unknown_raises(self):
        """The raw helper also validates the id (line 3348)."""
        with pytest.raises(ValueError, match="Keplerian elements"):
            _calc_keplerian_position_raw(99999, J2000)


# ---------------------------------------------------------------------------
# calc_white_moon_position, calc_waldemath_position, calc_proserpina
# ---------------------------------------------------------------------------


class TestSimpleCircularBodies:
    def test_white_moon(self):
        """White Moon uses the published circular-orbit formula (3420-3423)."""
        lon, lat, dist, dlon, dlat, ddist = calc_white_moon_position(J2000)
        assert 0.0 <= lon < 360.0
        assert dist == pytest.approx(0.05280098949)
        assert lat == 0.0
        assert dlat == 0.0
        assert ddist == 0.0

    def test_white_moon_ignores_flag(self):
        """The legacy ``use_true_lilith`` flag is ignored."""
        a = calc_white_moon_position(J2000, use_true_lilith=False)
        b = calc_white_moon_position(J2000, use_true_lilith=True)
        assert a == b

    def test_waldemath_position(self):
        """Waldemath delegates to calc_waldemath (line 3461)."""
        lon, lat, dist, dlon, dlat, ddist = calc_waldemath_position(J2000)
        assert 0.0 <= lon < 360.0
        assert dist > 0.0
        assert lat == 0.0

    def test_proserpina(self):
        """Proserpina propagates a circular orbit (lines 3504-3526)."""
        lon, lat, dist, dlon, dlat, ddist = calc_proserpina(J2000)
        assert 0.0 <= lon < 360.0
        assert dist > 0.0
        assert lat == 0.0
        assert dlat == 0.0
        assert ddist == 0.0


# ---------------------------------------------------------------------------
# calc_vulcan velocity wrap (line 3071) + raw helper
# ---------------------------------------------------------------------------


class TestVulcan:
    def test_calc_vulcan_negative_wrap(self):
        """Fast prograde Vulcan wraps negative at the origin (line 3071)."""
        result = calc_vulcan(2449553.0)
        assert -180.0 < result[3] < 180.0

    def test_calc_vulcan_basic(self):
        """A baseline Vulcan calculation succeeds."""
        result = calc_vulcan(J2000)
        assert 0.0 <= result[0] < 360.0
        assert result[2] > 0.0

    def test_raw_helper(self):
        """The raw Vulcan helper returns a 3-tuple."""
        out = _calc_vulcan_raw(J2000)
        assert len(out) == 3
        assert out[2] > 0.0


# ---------------------------------------------------------------------------
# Planet X Lowell / Pickering delegators and legacy raw helpers
# ---------------------------------------------------------------------------


class TestPlanetX:
    def test_lowell_delegates(self):
        """calc_planet_x_lowell delegates to the fictitious path (line 3572)."""
        out = calc_planet_x_lowell(J2000)
        assert 0.0 <= out[0] < 360.0
        assert out[2] > 0.0

    def test_lowell_raw_legacy(self):
        """The retained legacy Lowell raw helper still computes (3579-3621)."""
        out = H._calc_planet_x_lowell_raw(J2000)
        assert len(out) == 3
        assert out[2] > 0.0

    def test_pickering_delegates(self):
        """calc_planet_x_pickering delegates to the fictitious path (3667)."""
        out = calc_planet_x_pickering(J2000)
        assert 0.0 <= out[0] < 360.0
        assert out[2] > 0.0

    def test_pickering_raw_legacy(self):
        """The retained legacy Pickering raw helper computes (3674-3716)."""
        out = H._calc_planet_x_pickering_raw(J2000)
        assert len(out) == 3
        assert out[2] > 0.0


# ---------------------------------------------------------------------------
# calc_hypothetical_position dispatch table (lines 3749-3793)
# ---------------------------------------------------------------------------


class TestHypotheticalDispatch:
    def test_non_hypothetical_raises(self):
        """A non-hypothetical id raises ValueError (lines 3749-3750)."""
        with pytest.raises(ValueError, match="not a hypothetical body"):
            calc_hypothetical_position(0, J2000)

    def test_uranian_dispatch(self):
        """Uranian ids route to calc_uranian_planet (lines 3754-3755)."""
        out = calc_hypothetical_position(CUPIDO, J2000)
        assert 0.0 <= out[0] < 360.0

    def test_fictitious_dispatch(self):
        """Fictitious ids route to calc_fictitious_position (3757-3758)."""
        out = calc_hypothetical_position(NIBIRU, J2000)
        assert 0.0 <= out[0] < 360.0

    def test_isis_dispatch(self):
        """ISIS routes to calc_transpluto_position (lines 3761-3762)."""
        out = calc_hypothetical_position(ISIS, J2000)
        assert 0.0 <= out[0] < 360.0

    def test_vulcan_dispatch(self):
        """VULCAN routes to calc_vulcan (lines 3765-3766)."""
        out = calc_hypothetical_position(VULCAN, J2000)
        assert 0.0 <= out[0] < 360.0

    def test_white_moon_dispatch(self):
        """WHITE_MOON routes to calc_white_moon_position (3769-3770)."""
        out = calc_hypothetical_position(WHITE_MOON, J2000)
        assert 0.0 <= out[0] < 360.0

    def test_waldemath_dispatch(self):
        """WALDEMATH routes to calc_waldemath_position (3773-3774)."""
        out = calc_hypothetical_position(WALDEMATH, J2000)
        assert 0.0 <= out[0] < 360.0

    def test_proserpina_dispatch(self):
        """PROSERPINA routes to calc_proserpina (lines 3777-3778)."""
        out = calc_hypothetical_position(PROSERPINA, J2000)
        assert 0.0 <= out[0] < 360.0

    def test_lowell_dispatch(self):
        """PLANET_X_LOWELL routes to calc_planet_x_lowell (3781-3782)."""
        out = calc_hypothetical_position(PLANET_X_LOWELL, J2000)
        assert 0.0 <= out[0] < 360.0

    def test_pickering_dispatch(self):
        """PLANET_X_PICKERING routes to calc_planet_x_pickering (3785-3786)."""
        out = calc_hypothetical_position(PLANET_X_PICKERING, J2000)
        assert 0.0 <= out[0] < 360.0

    def test_unmapped_hypothetical_returns_zeros(self):
        """A hypothetical id with no handler returns zeros (line 3793)."""
        # 69 = FICT_OFFSET + 29: inside the hypothetical range but unmapped.
        out = calc_hypothetical_position(69, J2000)
        assert out == (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def _unwrapped_central_dlon(body_id, jd, h=1.0):
    """Independent truth: central finite-difference of the UNWRAPPED longitude.

    calc_uranian_planet returns the longitude normalised to [0, 360), so a
    correct daily speed is the central difference of the longitude after
    unwrapping the 0/360 jump.
    """
    p0 = calc_uranian_planet(body_id, jd - h)[0]
    p1 = calc_uranian_planet(body_id, jd + h)[0]
    d = (p1 - p0 + 180.0) % 360.0 - 180.0
    return d / (2.0 * h)


class TestUranianLongitudeSpeed:
    """Longitude speed of calc_uranian_planet across the 0/360 boundary.

    Regression for the wrap-around defect (audit round 11): the central
    difference was wrapped with a +/-180 threshold *after* dividing by
    2*dt_step, so a boundary crossing (raw jump ~360 -> halved to ~180) was not
    corrected and the reported speed jumped to ~-180 deg/day instead of the
    body's mean motion.
    """

    def test_speed_matches_unwrapped_finite_difference(self):
        from libephemeris.hypothetical import URANIAN_KEPLERIAN_ELEMENTS

        for body_id in URANIAN_KEPLERIAN_ELEMENTS:
            for jd in (J2000, J2000 - 36525.0, J2000 + 36525.0):
                reported = calc_uranian_planet(body_id, jd)[3]
                expected = _unwrapped_central_dlon(body_id, jd)
                assert reported == pytest.approx(expected, abs=1e-9), (
                    f"body {body_id} jd {jd}: dlon {reported} != {expected}"
                )

    def test_speed_is_mean_motion_across_0_360_crossing(self):
        # CUPIDO crosses the 0/360 boundary near JD 2386563.5; its longitude
        # speed there is the mean motion (~0.0038 deg/day), not the ~180 deg/day
        # artefact the old threshold produced. |dlon| must stay tiny.
        for jd in (2386562.0, 2386563.0, 2386564.0, 2386565.0):
            dlon = calc_uranian_planet(CUPIDO, jd)[3]
            assert abs(dlon) < 0.01, f"jd {jd}: dlon {dlon} (wrap artefact)"
            assert dlon == pytest.approx(
                _unwrapped_central_dlon(CUPIDO, jd), abs=1e-9
            )
