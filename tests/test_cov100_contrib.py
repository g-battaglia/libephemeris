# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for ``libephemeris.contrib``.

These tests exercise the contrib helper functions deterministically,
monkeypatching the ephemeris seams (``_calc_ut``, ``fixstar_ut``,
``houses``) so no network access or DE440 data is required. The goal is
to drive every conditional branch in the contrib module.
"""

from __future__ import annotations

import math

import pytest

import libephemeris.contrib as contrib


# ---------------------------------------------------------------------------
# house_system_name (lines 401-403)
# ---------------------------------------------------------------------------


def test_house_system_name_known_code():
    """A known letter code resolves to its human-readable name."""
    assert contrib.house_system_name("P") == "Placidus"


def test_house_system_name_unknown_code():
    """An unknown but str code returns ``Unknown``."""
    assert contrib.house_system_name("z") == "Unknown"


def test_house_system_name_non_str_raises():
    """A non-str argument raises a TypeError naming the bad type."""
    with pytest.raises(TypeError) as exc:
        contrib.house_system_name(ord("P"))  # type: ignore[arg-type]
    assert "expects a str code" in str(exc.value)
    assert "int" in str(exc.value)


# ---------------------------------------------------------------------------
# ochchabala (lines 458, 464)
# ---------------------------------------------------------------------------


def test_ochchabala_unknown_planet_returns_zero():
    """A planet absent from the exaltation table scores 0.0."""
    assert contrib.ochchabala(contrib.RAHU, 123.0) == 0.0


def test_ochchabala_at_exaltation_is_sixty():
    """Sun exalted at 10 Aries scores the full 60 strength units."""
    ex_long = contrib.ARIES * 30.0 + 10.0
    assert contrib.ochchabala(contrib.SURYA, ex_long) == pytest.approx(60.0)


def test_ochchabala_far_side_takes_short_arc():
    """A separation > 180 wraps to the short arc (line 464)."""
    # Sun exalted at 10 deg. Place it at 350 deg => raw diff 340 -> 20.
    val = contrib.ochchabala(contrib.SURYA, 350.0)
    expected = 60.0 * (1.0 - 20.0 / 180.0)
    assert val == pytest.approx(expected)


# ---------------------------------------------------------------------------
# residential_strength (lines 475-476)
# ---------------------------------------------------------------------------


def test_residential_strength_start_of_sign():
    """A planet at the start of a sign has full (1.0) residential strength."""
    assert contrib.residential_strength(0.0) == pytest.approx(1.0)


def test_residential_strength_mid_sign():
    """At 15 deg into a sign the strength is halfway (0.5)."""
    assert contrib.residential_strength(45.0) == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# naisargika_relation (lines 498-499)
# ---------------------------------------------------------------------------


def test_naisargika_relation_known_pair():
    """A pair present in the table returns its tabulated relation."""
    assert contrib.naisargika_relation(contrib.SURYA, contrib.CHANDRA) == 2


def test_naisargika_relation_unknown_pair_is_neutral():
    """A pair missing from the table defaults to neutral (1)."""
    assert contrib.naisargika_relation(999, 998) == 1


# ---------------------------------------------------------------------------
# tatkalika_relation (line 510)
# ---------------------------------------------------------------------------


def test_tatkalika_relation_same_house_neutral():
    """Same house separation is neutral (1)."""
    assert contrib.tatkalika_relation(1, 1) == 1


def test_tatkalika_relation_friend_house():
    """A friendly house separation returns 2."""
    # separation of 2 houses -> friend
    assert contrib.tatkalika_relation(1, 3) == 2


def test_tatkalika_relation_enemy_house():
    """An unfriendly house separation returns 0."""
    # separation of 6 houses -> enemy
    assert contrib.tatkalika_relation(1, 7) == 0


# ---------------------------------------------------------------------------
# raman_houses (lines 520-523)
# ---------------------------------------------------------------------------


def test_raman_houses_wraps_houses(monkeypatch):
    """raman_houses calls the canonical houses with Sripati code 'S'."""
    captured = {}

    def fake_houses(jd, geolat, geolon, hsys):
        captured["args"] = (jd, geolat, geolon, hsys)
        cusps = tuple(float(i * 30) for i in range(12))
        ascmc = (0.0,) * 10
        return cusps, ascmc

    monkeypatch.setattr("libephemeris.houses", fake_houses, raising=True)

    result = contrib.raman_houses(2451545.0, 12.5, 41.9)

    assert captured["args"] == (2451545.0, 41.9, 12.5, ord("S"))
    assert result == tuple(float(i * 30) for i in range(12))


# ---------------------------------------------------------------------------
# saturn_4_stars (lines 533-553, including the exception path 545-552)
# ---------------------------------------------------------------------------


def test_saturn_4_stars_success(monkeypatch):
    """All four stars resolve to their fixstar_ut longitudes."""
    lookup = {
        "Aldebaran": 69.0,
        "Regulus": 149.0,
        "Antares": 249.0,
        "Fomalhaut": 333.0,
    }

    def fake_fixstar_ut(name, jd):
        return ((lookup[name], 0.0, 1.0, 0.0, 0.0, 0.0), name, 0)

    monkeypatch.setattr("libephemeris.fixstar_ut", fake_fixstar_ut, raising=True)

    out = contrib.saturn_4_stars(2451545.0)
    assert out == (69.0, 149.0, 249.0, 333.0)


def test_saturn_4_stars_failure_returns_nan(monkeypatch):
    """A fixstar_ut failure yields NaN for that star and logs (545-552)."""

    def boom(name, jd):
        raise RuntimeError("no star data")

    monkeypatch.setattr("libephemeris.fixstar_ut", boom, raising=True)

    out = contrib.saturn_4_stars(2451545.0)
    assert len(out) == 4
    assert all(math.isnan(v) for v in out)


# ---------------------------------------------------------------------------
# jdnow (lines 563-566)
# ---------------------------------------------------------------------------


def test_jdnow_returns_float():
    """jdnow returns a finite Julian Day around the modern era."""
    jd = contrib.jdnow()
    assert isinstance(jd, float)
    # Some time after J2000 and well before far future.
    assert 2451545.0 < jd < 2600000.0


# ---------------------------------------------------------------------------
# years_diff (lines 606-608)
# ---------------------------------------------------------------------------


def test_years_diff_one_year():
    """Exactly one calendar year apart is ~1.0 decimal years."""
    val = contrib.years_diff("2000-01-01", "2001-01-01")
    assert val == pytest.approx(366.0 / 365.25, abs=1e-9)


def test_years_diff_space_separator():
    """A space-separated datetime parses (the .replace branch)."""
    val = contrib.years_diff("2000-01-01 00:00:00", "2000-01-02 00:00:00")
    assert val == pytest.approx(1.0 / 365.25, abs=1e-9)


# ---------------------------------------------------------------------------
# geoc2d (lines 619-629)
# ---------------------------------------------------------------------------


def test_geoc2d_north():
    """'12N34' parses to positive decimal degrees."""
    assert contrib.geoc2d("12N34") == pytest.approx(12 + 34 / 60.0)


def test_geoc2d_south_negative():
    """'12S30' parses to a negative value (S branch)."""
    assert contrib.geoc2d("12S30") == pytest.approx(-(12 + 30 / 60.0))


def test_geoc2d_west_negative():
    """'5W15' parses to a negative value (W branch)."""
    assert contrib.geoc2d("5W15") == pytest.approx(-(5 + 15 / 60.0))


def test_geoc2d_empty_parts_default_to_zero():
    """Missing degree/minute parts default to 0.0."""
    # "N" with no degree and no minute -> 0.0
    assert contrib.geoc2d("N") == pytest.approx(0.0)


def test_geoc2d_plain_decimal_fallback():
    """A bare decimal with no hemisphere letter falls through to float()."""
    assert contrib.geoc2d("41.9") == pytest.approx(41.9)


# ---------------------------------------------------------------------------
# match_aspect3 / match_aspect4 (lines 708-712, 724-726)
# ---------------------------------------------------------------------------


def test_match_aspect3_positive_signed():
    """A small positive separation keeps the signed delta positive."""
    is_match, delta, signed = contrib.match_aspect3(0, 0.0, 0.0, 1, 10.0, 1.0)
    assert is_match == 0  # 10 deg separation, conjunction with 1 deg orb
    assert delta == pytest.approx(10.0)
    assert signed == pytest.approx(10.0)


def test_match_aspect3_signed_wraps_negative():
    """A separation past 180 wraps the signed delta negative (710-711)."""
    is_match, delta, signed = contrib.match_aspect3(0, 0.0, 0.0, 1, 200.0, 1.0)
    assert signed == pytest.approx(-160.0)


def test_match_aspect4_returns_abs_separation():
    """match_aspect4 appends the absolute angular separation."""
    is_match, delta, signed, sep = contrib.match_aspect4(0, 0.0, 0.0, 1, 200.0, 1.0)
    # angular separation between 0 and 200 is 160
    assert sep == pytest.approx(160.0)
    assert signed == pytest.approx(-160.0)


# ---------------------------------------------------------------------------
# _walk_until_aspect (lines 743-750) via next_aspect family
# ---------------------------------------------------------------------------


def _fixed_calc_ut(lon, speed=1.0):
    """Build a fake _calc_ut returning a constant longitude/speed."""

    def _impl(jd, body, flags):
        return (lon, 0.0, 1.0, speed, 0.0, 0.0), 0

    return _impl


def test_walk_until_aspect_found(monkeypatch):
    """next_aspect returns the first JD whose orb is within tolerance."""
    # Body sits permanently at 120 deg => trine (120) with target 0 deg.
    monkeypatch.setattr(contrib, "_calc_ut", _fixed_calc_ut(120.0))
    jd, delta = contrib.next_aspect(0, 120.0, 0.0, 2451545.0, orb=1.0)
    assert jd == pytest.approx(2451545.0)  # matches on first sample
    assert delta == pytest.approx(0.0)


def test_walk_until_aspect_not_found(monkeypatch):
    """An aspect never reached returns (nan, inf) after the window (750)."""
    # Body fixed at 10 deg; ask for an exact 0 orb conjunction it never makes.
    monkeypatch.setattr(contrib, "_calc_ut", _fixed_calc_ut(10.0))
    jd, delta = contrib.next_aspect(0, 0.0, 0.0, 2451545.0, orb=0.0)
    assert math.isnan(jd)
    assert math.isinf(delta)


# ---------------------------------------------------------------------------
# next_aspect / next_aspect2 / next_aspect_cusp / next_aspect_cusp2
# (lines 775, 787, 798, 808)
# ---------------------------------------------------------------------------


def test_next_aspect2_default_orb(monkeypatch):
    """next_aspect2 delegates to next_aspect with the default 1 deg orb."""
    monkeypatch.setattr(contrib, "_calc_ut", _fixed_calc_ut(60.0))
    jd, delta = contrib.next_aspect2(0, 60.0, 0.0, 2451545.0)
    assert jd == pytest.approx(2451545.0)
    assert delta == pytest.approx(0.0)


def test_next_aspect_cusp(monkeypatch):
    """next_aspect_cusp resolves an aspect to a fixed cusp longitude."""
    monkeypatch.setattr(contrib, "_calc_ut", _fixed_calc_ut(90.0))
    jd, delta = contrib.next_aspect_cusp(0, 90.0, 0.0, 2451545.0, orb=1.0)
    assert jd == pytest.approx(2451545.0)
    assert delta == pytest.approx(0.0)


def test_next_aspect_cusp2_default_orb(monkeypatch):
    """next_aspect_cusp2 uses the default orb of 1 deg."""
    monkeypatch.setattr(contrib, "_calc_ut", _fixed_calc_ut(90.0))
    jd, delta = contrib.next_aspect_cusp2(0, 90.0, 0.0, 2451545.0)
    assert jd == pytest.approx(2451545.0)
    assert delta == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# next_aspect_with / next_aspect_with2 (lines 827-836, 846)
# ---------------------------------------------------------------------------


def test_next_aspect_with_found(monkeypatch):
    """Two bodies forming the aspect on the first sample return that JD."""

    def fake_calc(jd, body, flags):
        # body 0 at 0 deg, body 1 at 90 deg -> square (90)
        lon = 0.0 if body == 0 else 90.0
        return (lon, 0.0, 1.0, 1.0, 0.0, 0.0), 0

    monkeypatch.setattr(contrib, "_calc_ut", fake_calc)
    jd, delta = contrib.next_aspect_with(0, 1, 90.0, 2451545.0, orb=1.0)
    assert jd == pytest.approx(2451545.0)
    assert delta == pytest.approx(0.0)


def test_next_aspect_with_not_found(monkeypatch):
    """An unreachable two-body aspect returns (nan, inf) (line 836)."""

    def fake_calc(jd, body, flags):
        lon = 0.0 if body == 0 else 10.0
        return (lon, 0.0, 1.0, 1.0, 0.0, 0.0), 0

    monkeypatch.setattr(contrib, "_calc_ut", fake_calc)
    jd, delta = contrib.next_aspect_with(0, 1, 0.0, 2451545.0, orb=0.0)
    assert math.isnan(jd)
    assert math.isinf(delta)


def test_next_aspect_with2_default_orb(monkeypatch):
    """next_aspect_with2 delegates with the default orb of 1 deg."""

    def fake_calc(jd, body, flags):
        lon = 0.0 if body == 0 else 180.0
        return (lon, 0.0, 1.0, 1.0, 0.0, 0.0), 0

    monkeypatch.setattr(contrib, "_calc_ut", fake_calc)
    jd, delta = contrib.next_aspect_with2(0, 1, 180.0, 2451545.0)
    assert jd == pytest.approx(2451545.0)
    assert delta == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# next_retro (lines 871-881)
# ---------------------------------------------------------------------------


def test_next_retro_detects_station(monkeypatch):
    """next_retro returns the JD where speed first goes negative."""
    start = 2451545.0

    def fake_calc(jd, body, flags):
        # Speed positive for the first step, negative afterward.
        speed = 1.0 if jd < start + 0.5 else -1.0
        return (0.0, 0.0, 1.0, speed, 0.0, 0.0), 0

    monkeypatch.setattr(contrib, "_calc_ut", fake_calc)
    jd = contrib.next_retro(0, start)
    # First sample speed >= 0, second sample speed < 0 -> detected at start+0.5
    assert jd == pytest.approx(start + 0.5)


def test_next_retro_not_found(monkeypatch):
    """next_retro returns NaN if no station occurs within max_days (881)."""

    monkeypatch.setattr(contrib, "_calc_ut", _fixed_calc_ut(0.0, speed=1.0))
    jd = contrib.next_retro(0, 2451545.0, max_days=2.0)
    assert math.isnan(jd)
