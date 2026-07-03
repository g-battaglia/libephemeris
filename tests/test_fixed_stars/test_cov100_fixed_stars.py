# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for :mod:`libephemeris.fixed_stars`.

These tests exercise the resolver tiers, the LEB and Skyfield calculation
paths, the post-calculation flag transformations, the velocity wraparound
handling, the batch path, and the magnitude lookups. They drive specific
branches that the broader test suite does not reach, using deterministic
inputs and targeted monkeypatching of internal seams (fake LEB readers,
stubbed calculation helpers, injected catalog/alias entries).

All tests run fully offline: the DE440 kernels are loaded from disk and no
network access is required.
"""

from __future__ import annotations

import pytest

import libephemeris.fixed_stars as fs
import libephemeris.state as state
from libephemeris.constants import (
    FLG_EQUATORIAL,
    FLG_J2000,
    FLG_NONUT,
    FLG_SPEED,
    FLG_SPEED3,
    FLG_TOPOCTR,
    FLG_XYZ,
    REGULUS,
)
from libephemeris.exceptions import Error

JD = 2451545.0  # J2000.0 TT


# Realistic Earth barycentric ICRS position/velocity near J2000 (AU, AU/day).
# Values need only be finite for the LEB pipeline math; they are not used to
# assert precise positions here.
_EARTH_POS = (-0.1771, 0.8876, 0.3848)
_EARTH_VEL = (-0.01720, -0.00292, -0.00127)


class _FakeLEBReader:
    """Minimal LEB reader stub for the star LEB pipeline.

    Returns Earth state for the EARTH body and raises ``KeyError`` for every
    other body so the gravitational-deflection deflector lookups are skipped.
    Deliberately omits ``eval_nutation`` so ``_get_leb_frame_data`` raises and
    the erfa frame-data fallback branch is taken.
    """

    def eval_body(self, body_id, jd):
        if body_id == 14:  # EARTH
            return (_EARTH_POS, _EARTH_VEL)
        raise KeyError(body_id)


# ---------------------------------------------------------------------------
# propagate_proper_motion
# ---------------------------------------------------------------------------


def test_propagate_proper_motion_basic():
    """Propagation returns a normalized RA and clamped declination."""
    ra, dec = fs.propagate_proper_motion(10.0, 5.0, 0.5, 0.3, JD, JD + 36525.0)
    assert 0.0 <= ra < 360.0
    assert -90.0 <= dec <= 90.0


def test_propagate_proper_motion_clamps_and_pole():
    """Declination clamps at +/-90 and RA stays unchanged at the pole."""
    # Force dec_target above +90 (clamp high branch) and cos(dec)~0 at pole.
    ra_hi, dec_hi = fs.propagate_proper_motion(
        100.0, 89.999, 1.0, 36000.0, JD, JD + 36525.0
    )
    assert dec_hi == 90.0
    # Start exactly at the pole so cos(dec_epoch) ~ 0 -> RA unchanged.
    ra_pole, _ = fs.propagate_proper_motion(123.0, 90.0, 5.0, -36000.0, JD, JD + 36525.0)
    assert ra_pole == 123.0
    # Force dec_target below -90 (clamp low branch).
    _, dec_lo = fs.propagate_proper_motion(
        50.0, -89.999, 0.0, -36000.0, JD, JD + 36525.0
    )
    assert dec_lo == -90.0


# ---------------------------------------------------------------------------
# _fuzzy_match_star
# ---------------------------------------------------------------------------


def test_fuzzy_match_too_short():
    """Inputs that normalize too short return None without matching."""
    assert fs._fuzzy_match_star("ab") is None


def test_fuzzy_match_duplicate_catalog_entry(monkeypatch):
    """A duplicate catalog entry (same id) is not appended twice (1941->1939)."""
    entry = next(e for e in fs.STAR_CATALOG if e.name == "Betelgeuse")
    dup = fs.StarCatalogEntry(
        id=entry.id,
        name=entry.name,
        nomenclature=entry.nomenclature,
        hip_number=entry.hip_number,
        data=entry.data,
        magnitude=entry.magnitude,
    )
    monkeypatch.setattr(fs, "STAR_CATALOG", [entry, dup])
    monkeypatch.setattr(fs, "STAR_ALIASES", {})
    assert fs._fuzzy_match_star(entry.name) == entry.id


def test_fuzzy_match_ambiguous(monkeypatch):
    """Two aliases normalizing alike to different ids => None (1956)."""
    monkeypatch.setattr(fs, "STAR_CATALOG", [])
    monkeypatch.setattr(fs, "STAR_ALIASES", {"QWXYZ": 111, "QWXYZ ": 222})
    assert fs._normalize_phonetic("QWXYZ") == fs._normalize_phonetic("QWXYZ ")
    assert fs._fuzzy_match_star("QWXYZ") is None


# ---------------------------------------------------------------------------
# resolve_star_name
# ---------------------------------------------------------------------------


def test_resolve_star_name_empty():
    """Empty name resolves to None."""
    assert fs.resolve_star_name("") is None


def test_resolve_star_name_comma_prefix_alias():
    """Comma-prefix matches an alias prefix (2018-2020)."""
    # "DOG STAR" is a Sirius alias; no catalog name starts with "DOG".
    assert fs.resolve_star_name(",DOG STAR") is not None


def test_resolve_star_name_comma_prefix_nomenclature():
    """Comma-prefix matches a nomenclature prefix (2023-2025)."""
    # "ALAND" matches Alpheratz's nomenclature "alAnd" but no name/alias.
    assert fs.resolve_star_name(",alAnd") is not None


def test_resolve_star_name_comma_prefix_none():
    """Comma-prefix with no match returns None (2027) and empty prefix None."""
    assert fs.resolve_star_name(",zzzzzqqq") is None
    assert fs.resolve_star_name(",   ") is None


def test_resolve_star_name_exact_nomenclature():
    """Exact nomenclature match resolves (2047)."""
    assert fs.resolve_star_name("alAnd") is not None


def test_resolve_star_name_bayer_designation():
    """Bayer designation parsing resolves (2052-2054).

    "Theta Octantis" is not an exact name/alias/nomenclature, so resolution
    reaches the Bayer-parser tier rather than short-circuiting earlier.
    """
    assert fs.resolve_star_name("Theta Octantis") is not None


def test_resolve_star_name_flamsteed_miss():
    """Parsed Flamsteed not present in aliases falls through (2060->2066)."""
    # "99 Leonis" parses to "99 LEO" which is not an alias and resolves None.
    assert fs.resolve_star_name("99 Leonis") is None


def test_resolve_star_name_prefix_catalog_name():
    """Prefix match against a catalog name resolves (2072)."""
    # No alias starts with "VINDEMI" but the catalog name "Vindemiatrix" does.
    assert fs.resolve_star_name("VINDEMI") is not None


def test_resolve_star_name_fuzzy():
    """Fuzzy phonetic fallback resolves a misspelling (2077).

    "Reegulus" is not an exact alias nor a prefix of any alias/name, so it
    falls through to the phonetic fuzzy matcher, which maps it to Regulus.
    """
    assert fs.resolve_star_name("Reegulus") == REGULUS


# ---------------------------------------------------------------------------
# _calc_star_position_from_observer / velocity guards
# ---------------------------------------------------------------------------


def test_calc_from_observer_unknown_star():
    """Unknown star id raises ValueError (2107)."""
    with pytest.raises(ValueError):
        fs._calc_star_position_from_observer(987654321, None)


def test_calc_velocity_unknown_star():
    """Unknown star id raises ValueError (2430)."""
    with pytest.raises(ValueError):
        fs.calc_fixed_star_velocity(987654321, JD)


def test_calc_velocity_wraparound_positive(monkeypatch):
    """speed_lon > 180 wraps down by 360 (2453)."""

    def fake(star_id, jd, *a, **k):
        if abs(jd - JD) < 1e-9:
            return (0.0, 0.0, 1e6)
        if jd < JD:
            return (1.0, 0.0, 1e6)  # prev
        return (359.0, 0.0, 1e6)  # next -> 358 -> -2

    monkeypatch.setattr(fs, "calc_fixed_star_position", fake)
    result = fs.calc_fixed_star_velocity(REGULUS, JD)
    assert result[3] == pytest.approx(-2.0)


def test_calc_velocity_wraparound_negative(monkeypatch):
    """speed_lon < -180 wraps up by 360 (2455)."""

    def fake(star_id, jd, *a, **k):
        if abs(jd - JD) < 1e-9:
            return (0.0, 0.0, 1e6)
        if jd < JD:
            return (359.0, 0.0, 1e6)  # prev
        return (1.0, 0.0, 1e6)  # next -> -358 -> +2

    monkeypatch.setattr(fs, "calc_fixed_star_position", fake)
    result = fs.calc_fixed_star_velocity(REGULUS, JD)
    assert result[3] == pytest.approx(2.0)


# ---------------------------------------------------------------------------
# _calc_star_position_leb
# ---------------------------------------------------------------------------


def test_calc_star_position_leb_no_reader(monkeypatch):
    """Missing LEB reader raises KeyError (2188)."""
    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    with pytest.raises(KeyError):
        fs._calc_star_position_leb(REGULUS, JD)


def test_calc_star_position_leb_unknown_star(monkeypatch):
    """Unknown star id with reader present raises ValueError (2190-2191)."""
    monkeypatch.setattr(state, "get_leb_reader", lambda: object())
    with pytest.raises(ValueError):
        fs._calc_star_position_leb(987654321, JD)


def test_calc_star_position_leb_erfa_fallback(monkeypatch):
    """Reader without nutation triggers the erfa frame fallback (2260-2266)."""
    monkeypatch.setattr(state, "get_leb_reader", lambda: _FakeLEBReader())
    lon, lat, dist = fs._calc_star_position_leb(REGULUS, JD)
    assert 0.0 <= lon < 360.0
    assert dist > 0.0


def test_calc_star_position_leb_noaberr_and_j2000(monkeypatch):
    """LEB path with noaberr and J2000 frame both succeed."""
    monkeypatch.setattr(state, "get_leb_reader", lambda: _FakeLEBReader())
    lon_na, _, _ = fs._calc_star_position_leb(REGULUS, JD, noaberr=True)
    assert 0.0 <= lon_na < 360.0
    lon_j, _, _ = fs._calc_star_position_leb(REGULUS, JD, j2000_frame=True)
    assert 0.0 <= lon_j < 360.0


# ---------------------------------------------------------------------------
# _calc_star_position_skyfield topocentric branch
# ---------------------------------------------------------------------------


def test_calc_star_position_skyfield_topo():
    """Topocentric observer path runs (2318-2321)."""
    lon, lat, dist = fs._calc_star_position_skyfield(
        REGULUS, JD, topo=(12.5, 41.9, 0.0)
    )
    assert 0.0 <= lon < 360.0
    assert dist > 0.0


# ---------------------------------------------------------------------------
# calc_fixed_star_position LEB dispatch
# ---------------------------------------------------------------------------


def test_calc_position_leb_keyerror_falls_back(monkeypatch):
    """KeyError from LEB path falls back to Skyfield (2380)."""
    monkeypatch.setattr(state, "get_leb_reader", lambda: object())

    def raise_key(*a, **k):
        raise KeyError("no earth")

    monkeypatch.setattr(fs, "_calc_star_position_leb", raise_key)
    lon, lat, dist = fs.calc_fixed_star_position(REGULUS, JD)
    assert 0.0 <= lon < 360.0


def test_calc_position_leb_outside_range_falls_back(monkeypatch):
    """ValueError 'outside range' logs and falls back to Skyfield (2384-2386)."""
    monkeypatch.setattr(state, "get_leb_reader", lambda: object())

    def raise_oor(*a, **k):
        raise ValueError("jd is outside range of LEB file")

    monkeypatch.setattr(fs, "_calc_star_position_leb", raise_oor)
    lon, lat, dist = fs.calc_fixed_star_position(REGULUS, JD)
    assert 0.0 <= lon < 360.0


def test_calc_position_leb_valueerror_falls_back(monkeypatch):
    """Any LEB ValueError (out-of-range OR corruption) falls back to Skyfield,
    mirroring the planet path (calc_ut catches ValueError from the LEB
    reader and falls through)."""
    monkeypatch.setattr(state, "get_leb_reader", lambda: object())

    def raise_other(*a, **k):
        raise ValueError("Corrupted LEB2 data: simulated")

    monkeypatch.setattr(fs, "_calc_star_position_leb", raise_other)
    lon, lat, dist = fs.calc_fixed_star_position(REGULUS, JD)
    assert 0.0 <= lon < 360.0  # Skyfield fallback produced a real position


# ---------------------------------------------------------------------------
# _resolve_star_se
# ---------------------------------------------------------------------------


def test_resolve_star_se_empty():
    """Empty search yields the 'star name empty' error (2520)."""
    star_id, error, name = fs._resolve_star_se("   ")
    assert star_id == -1
    assert error == "star name empty"


def test_resolve_star_se_comma_empty_key():
    """Trailing comma with empty nomenclature key errors (2529)."""
    star_id, error, name = fs._resolve_star_se("Regulus,")
    assert star_id == -1
    assert error is not None


def test_resolve_star_se_digit_then_nondigit():
    """Sequential number stops at the first non-digit char (2538->2543)."""
    star_id, error, name = fs._resolve_star_se("5x")
    assert star_id != -1
    assert name is not None


def test_resolve_star_se_sequential_out_of_range():
    """Out-of-range sequential number errors (2546)."""
    star_id, error, name = fs._resolve_star_se("999999")
    assert star_id == -1
    assert "sequential" in error


def test_resolve_star_se_wildcard_invalid():
    """A '%' that is not a sole trailing char is invalid (2555-2556)."""
    star_id, error, _ = fs._resolve_star_se("Reg%ul")
    assert star_id == -1
    assert "invalid search string" in error
    star_id2, error2, _ = fs._resolve_star_se("Reg%%")
    assert star_id2 == -1
    assert "invalid search string" in error2


def test_resolve_star_se_wildcard_prefix():
    """Trailing '%' wildcard matches a prefix (2557-2560)."""
    star_id, error, name = fs._resolve_star_se("Reg%")
    assert star_id != -1
    assert name is not None


def test_resolve_star_se_wildcard_no_match():
    """Trailing '%' wildcard with no match errors (2561)."""
    star_id, error, _ = fs._resolve_star_se("Zzzzq%")
    assert star_id == -1
    assert "did not match" in error


def test_resolve_star_se_alias_id_without_entry(monkeypatch):
    """An alias whose id has no catalog entry exhausts the inner loop (2571->2569)."""
    aliases = dict(fs.STAR_ALIASES)
    aliases["zzzfakestar"] = 999999999
    monkeypatch.setattr(fs, "STAR_ALIASES", aliases)
    star_id, error, _ = fs._resolve_star_se("zzzfakestar")
    assert star_id == -1


# ---------------------------------------------------------------------------
# _preprocess_flags / _fixstar_topo
# ---------------------------------------------------------------------------


def test_preprocess_flags_speed3():
    """FLG_SPEED3 is converted so speed is computed (2615)."""
    result = fs.fixstar("Regulus", JD, FLG_SPEED3)
    assert result[0][3] != 0.0


def test_fixstar_topo_set():
    """A configured topo position is returned as a tuple (2646-2650)."""
    import libephemeris as ephem

    ephem.set_topo(12.5, 41.9, 0.0)
    lon, lat, alt = fs._fixstar_topo()
    assert lon == pytest.approx(12.5)
    assert lat == pytest.approx(41.9)


def test_fixstar_topo_unset_raises(monkeypatch):
    """Missing topo position raises Error (2641-2645)."""
    monkeypatch.setattr(state, "_TOPO", None)
    with pytest.raises(Error):
        fs.fixstar("Regulus", JD, FLG_TOPOCTR)


# ---------------------------------------------------------------------------
# _apply_fixstar_flags
# ---------------------------------------------------------------------------


def test_apply_flags_j2000_not_native():
    """J2000 precession path runs when not natively J2000 (2693-2695)."""
    res = (150.0, 0.4, 1e6, 0.001, 0.0, 0.0)
    out = fs._apply_fixstar_flags(res, JD, FLG_J2000, j2000_native=False)
    assert len(out) == 6


def test_apply_flags_nonut():
    """NONUT subtracts nutation in longitude (2701-2704)."""
    result = fs.fixstar("Regulus", JD, FLG_NONUT)
    assert 0.0 <= result[0][0] < 360.0


def test_apply_flags_nonut_equatorial():
    """NONUT + EQUATORIAL uses mean obliquity (2713)."""
    result = fs.fixstar("Regulus", JD, FLG_NONUT | FLG_EQUATORIAL)
    assert 0.0 <= result[0][0] < 360.0


def test_apply_flags_xyz():
    """XYZ output produces Cartesian coordinates and velocities (2750-2776)."""
    result = fs.fixstar("Regulus", JD, FLG_XYZ | FLG_SPEED)
    x, y, z, vx, vy, vz = result[0]
    assert abs(x) > 0.0 or abs(y) > 0.0 or abs(z) > 0.0


# ---------------------------------------------------------------------------
# fixstar_ut / fixstar exception handlers
# ---------------------------------------------------------------------------


def test_fixstar_ut_valueerror_wrapped(monkeypatch):
    """A ValueError from the calc path is wrapped as Error (2859-2860)."""

    def raise_val(*a, **k):
        raise ValueError("boom")

    monkeypatch.setattr(fs, "calc_fixed_star_position", raise_val)
    with pytest.raises(Error):
        fs.fixstar_ut("Regulus", JD, 0)


def test_fixstar_ut_error_passthrough(monkeypatch):
    """An Error from the calc path propagates unchanged (2857-2858)."""

    def raise_err(*a, **k):
        raise Error("passthrough")

    monkeypatch.setattr(fs, "calc_fixed_star_position", raise_err)
    with pytest.raises(Error, match="passthrough"):
        fs.fixstar_ut("Regulus", JD, 0)


def test_fixstar_valueerror_wrapped(monkeypatch):
    """fixstar wraps a ValueError as Error (3075-3076)."""

    def raise_val(*a, **k):
        raise ValueError("boom")

    monkeypatch.setattr(fs, "calc_fixed_star_position", raise_val)
    with pytest.raises(Error):
        fs.fixstar("Regulus", JD, 0)


def test_fixstar_error_passthrough(monkeypatch):
    """fixstar passes an Error through (3073-3074)."""

    def raise_err(*a, **k):
        raise Error("passthrough")

    monkeypatch.setattr(fs, "calc_fixed_star_position", raise_err)
    with pytest.raises(Error, match="passthrough"):
        fs.fixstar("Regulus", JD, 0)


def test_fixstar_speed_path():
    """fixstar with FLG_SPEED computes a velocity (3060-3063)."""
    result = fs.fixstar("Regulus", JD, FLG_SPEED)
    assert result[0][3] != 0.0


# ---------------------------------------------------------------------------
# batch_fixstars_ut
# ---------------------------------------------------------------------------


def test_batch_basic():
    """Skyfield batch path returns positions for resolved stars (2969-2994)."""
    results = fs.batch_fixstars_ut(["Regulus", "Spica"], JD, 0)
    assert results[0][1] == "Regulus,alLeo"
    assert results[1][1] == "Spica,alVir"


def test_batch_speed():
    """Skyfield batch path with FLG_SPEED computes velocity (2963-2967, 2984+)."""
    results = fs.batch_fixstars_ut(["Regulus"], JD, FLG_SPEED)
    assert results[0][0][3] != 0.0


def test_batch_all_unresolved_skip():
    """A batch of only unresolved stars (skip) returns early (resolved empty)."""
    results = fs.batch_fixstars_ut(["ZzzNope"], JD, 0, skip_errors=True)
    assert results == (None,)


def test_batch_skyfield_wraparound_positive(monkeypatch):
    """Skyfield batch speed_lon > 180 wraps down (2986)."""

    def fake(star_id, observer, *a, **k):
        # observer identity distinguishes prev/cur/next via closure counters.
        return fake.queue.pop(0)

    fake.queue = [
        (0.0, 0.0, 1e6),  # current
        (1.0, 0.0, 1e6),  # prev
        (359.0, 0.0, 1e6),  # next -> 358 -> -2
    ]
    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    monkeypatch.setattr(fs, "_calc_star_position_from_observer", fake)
    results = fs.batch_fixstars_ut(["Regulus"], JD, FLG_SPEED)
    assert results[0][0][3] == pytest.approx(-2.0)


def test_batch_skyfield_wraparound_negative(monkeypatch):
    """Skyfield batch speed_lon < -180 wraps up (2988)."""

    def fake(star_id, observer, *a, **k):
        return fake.queue.pop(0)

    fake.queue = [
        (0.0, 0.0, 1e6),  # current
        (359.0, 0.0, 1e6),  # prev
        (1.0, 0.0, 1e6),  # next -> -358 -> +2
    ]
    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    monkeypatch.setattr(fs, "_calc_star_position_from_observer", fake)
    results = fs.batch_fixstars_ut(["Regulus"], JD, FLG_SPEED)
    assert results[0][0][3] == pytest.approx(2.0)


def test_batch_skyfield_error_skip(monkeypatch):
    """A ValueError during a star calc is skipped when skip_errors (3003-3005)."""

    def raise_val(*a, **k):
        raise ValueError("boom")

    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    monkeypatch.setattr(fs, "_calc_star_position_from_observer", raise_val)
    results = fs.batch_fixstars_ut(["Regulus"], JD, 0, skip_errors=True)
    assert results == (None,)


def test_batch_skyfield_error_raise(monkeypatch):
    """A ValueError during a star calc is wrapped as Error without skip (3006)."""

    def raise_val(*a, **k):
        raise ValueError("boom")

    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    monkeypatch.setattr(fs, "_calc_star_position_from_observer", raise_val)
    with pytest.raises(Error):
        fs.batch_fixstars_ut(["Regulus"], JD, 0)


def test_batch_skyfield_error_passthrough(monkeypatch):
    """An Error during a star calc propagates without skip (2999, 3002)."""

    def raise_err(*a, **k):
        raise Error("calc failed")

    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    monkeypatch.setattr(fs, "_calc_star_position_from_observer", raise_err)
    with pytest.raises(Error, match="calc failed"):
        fs.batch_fixstars_ut(["Regulus"], JD, 0)


def test_batch_skyfield_error_skip_continue(monkeypatch):
    """An Error during a star calc is skipped when skip_errors (3001)."""

    def raise_err(*a, **k):
        raise Error("calc failed")

    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    monkeypatch.setattr(fs, "_calc_star_position_from_observer", raise_err)
    results = fs.batch_fixstars_ut(["Regulus"], JD, 0, skip_errors=True)
    assert results == (None,)


def test_batch_leb_keyerror_falls_back(monkeypatch):
    """KeyError on the LEB batch path falls back to Skyfield (2948-2949)."""
    monkeypatch.setattr(state, "get_leb_reader", lambda: object())

    def raise_key(*a, **k):
        raise KeyError("no earth")

    monkeypatch.setattr(fs, "_calc_star_position_leb", raise_key)
    results = fs.batch_fixstars_ut(["Regulus"], JD, 0)
    assert results[0][1] == "Regulus,alLeo"


def test_batch_leb_outside_range_falls_back(monkeypatch):
    """ValueError 'outside range' on the LEB batch path falls back (2950-2952)."""
    monkeypatch.setattr(state, "get_leb_reader", lambda: object())

    def raise_oor(*a, **k):
        raise ValueError("jd outside range of LEB")

    monkeypatch.setattr(fs, "_calc_star_position_leb", raise_oor)
    results = fs.batch_fixstars_ut(["Regulus"], JD, 0)
    assert results[0][1] == "Regulus,alLeo"


def test_batch_leb_valueerror_falls_back(monkeypatch):
    """Any LEB ValueError on the batch path falls back to Skyfield
    (out-of-range AND corrupted/truncated data alike)."""
    monkeypatch.setattr(state, "get_leb_reader", lambda: object())

    def raise_other(*a, **k):
        raise ValueError("Corrupted LEB2 data: simulated")

    monkeypatch.setattr(fs, "_calc_star_position_leb", raise_other)
    results = fs.batch_fixstars_ut(["Regulus"], JD, 0)
    assert results[0] is not None
    assert 0.0 <= results[0][0][0] < 360.0


def test_batch_leb_speed_wraparound_positive(monkeypatch):
    """LEB batch speed_lon > 180 wraps down (2932)."""

    def fake(star_id, jd, *a, **k):
        if abs(jd - JD) < 1e-9:
            return (0.0, 0.0, 1e6)
        if jd < JD:
            return (1.0, 0.0, 1e6)  # prev
        return (359.0, 0.0, 1e6)  # next -> 358 -> -2

    monkeypatch.setattr(state, "get_leb_reader", lambda: object())
    monkeypatch.setattr(fs, "_calc_star_position_leb", fake)
    results = fs.batch_fixstars_ut(["Regulus"], JD, FLG_SPEED)
    assert results[0][0][3] == pytest.approx(-2.0)


def test_batch_leb_speed_wraparound_negative(monkeypatch):
    """LEB batch speed_lon < -180 wraps up (2934)."""

    def fake(star_id, jd, *a, **k):
        if abs(jd - JD) < 1e-9:
            return (0.0, 0.0, 1e6)
        if jd < JD:
            return (359.0, 0.0, 1e6)  # prev
        return (1.0, 0.0, 1e6)  # next -> -358 -> +2

    monkeypatch.setattr(state, "get_leb_reader", lambda: object())
    monkeypatch.setattr(fs, "_calc_star_position_leb", fake)
    results = fs.batch_fixstars_ut(["Regulus"], JD, FLG_SPEED)
    assert results[0][0][3] == pytest.approx(2.0)


# ---------------------------------------------------------------------------
# _resolve_star2 (used by fixstar2*)
# ---------------------------------------------------------------------------


def test_resolve_star2_name_hip_fix():
    """A curated name->HIP correction resolves (3179-3181)."""
    entry, error = fs._resolve_star2("Alaraph")
    assert entry is not None
    assert error is None


def test_resolve_star2_bayer():
    """Bayer designation resolves via nomenclature (3188).

    "Theta Octantis" is not an exact name/alias/nomenclature/name-fix, so the
    Bayer-parser tier handles it.
    """
    entry, error = fs._resolve_star2("Theta Octantis")
    assert entry is not None
    assert entry.nomenclature == "thOct"


def test_resolve_star2_partial_nomenclature_single():
    """A unique partial nomenclature match resolves (3215, 3218)."""
    entry, error = fs._resolve_star2("zePh")
    assert entry is not None
    assert entry.nomenclature == "zePhe"


def test_resolve_star2_partial_nomenclature_ambiguous():
    """An ambiguous partial nomenclature match errors (3220-3221)."""
    entry, error = fs._resolve_star2("ALE")
    assert entry is None
    assert "Ambiguous" in error


def test_resolve_star2_alias_id_without_entry(monkeypatch):
    """An alias id without a catalog entry falls through (3169->3177)."""
    aliases = dict(fs.STAR_ALIASES)
    aliases["FAKEALIASXYZ"] = 99999999
    monkeypatch.setattr(fs, "STAR_ALIASES", aliases)
    entry, error = fs._resolve_star2("FAKEALIASXYZ")
    assert entry is None


def test_resolve_star2_flamsteed_id_without_entry(monkeypatch):
    """A Flamsteed alias id without a catalog entry falls through (3196->3201)."""
    aliases = dict(fs.STAR_ALIASES)
    aliases["99 LEO"] = 88888888
    monkeypatch.setattr(fs, "STAR_ALIASES", aliases)
    entry, error = fs._resolve_star2("99 Leonis")
    assert entry is None


# ---------------------------------------------------------------------------
# fixstar2_ut / fixstar2 speed + exception handlers
# ---------------------------------------------------------------------------


def test_fixstar2_ut_speed():
    """fixstar2_ut FLG_SPEED computes a velocity (3292-3295)."""
    result = fs.fixstar2_ut("Regulus", JD, FLG_SPEED)
    assert result[0][3] != 0.0


def test_fixstar2_ut_valueerror_wrapped(monkeypatch):
    """fixstar2_ut wraps a ValueError as Error (3308-3309)."""

    def raise_val(*a, **k):
        raise ValueError("boom")

    monkeypatch.setattr(fs, "calc_fixed_star_position", raise_val)
    with pytest.raises(Error):
        fs.fixstar2_ut("Regulus", JD, 0)


def test_fixstar2_ut_error_passthrough(monkeypatch):
    """fixstar2_ut passes an Error through (3306-3307)."""

    def raise_err(*a, **k):
        raise Error("passthrough")

    monkeypatch.setattr(fs, "calc_fixed_star_position", raise_err)
    with pytest.raises(Error, match="passthrough"):
        fs.fixstar2_ut("Regulus", JD, 0)


def test_fixstar2_speed():
    """fixstar2 FLG_SPEED computes a velocity (3365-3368)."""
    result = fs.fixstar2("Regulus", JD, FLG_SPEED)
    assert result[0][3] != 0.0


def test_fixstar2_valueerror_wrapped(monkeypatch):
    """fixstar2 wraps a ValueError as Error (3381-3382)."""

    def raise_val(*a, **k):
        raise ValueError("boom")

    monkeypatch.setattr(fs, "calc_fixed_star_position", raise_val)
    with pytest.raises(Error):
        fs.fixstar2("Regulus", JD, 0)


def test_fixstar2_error_passthrough(monkeypatch):
    """fixstar2 passes an Error through (3379-3380)."""

    def raise_err(*a, **k):
        raise Error("passthrough")

    monkeypatch.setattr(fs, "calc_fixed_star_position", raise_err)
    with pytest.raises(Error, match="passthrough"):
        fs.fixstar2("Regulus", JD, 0)


# ---------------------------------------------------------------------------
# fixstar_mag
# ---------------------------------------------------------------------------


def test_fixstar_mag_basic():
    """Magnitude lookup returns the catalog magnitude and full name."""
    mag, name = fs.fixstar_mag("Regulus")
    assert name == "Regulus,alLeo"
    assert isinstance(mag, float)


def test_fixstar_mag_missing_magnitude(monkeypatch):
    """A resolved id with no magnitude raises Error (4547)."""
    monkeypatch.setattr(fs, "_resolve_star_id", lambda s: (77777, None, "Fake"))
    with pytest.raises(Error, match="Magnitude not available"):
        fs.fixstar_mag("whatever")


def test_fixstar_mag_fallback_name(monkeypatch):
    """A magnitude id without a catalog entry uses the canonical name (4556-4557)."""
    mags = dict(fs._STAR_MAGNITUDES)
    mags[77777] = 3.14
    monkeypatch.setattr(fs, "_STAR_MAGNITUDES", mags)
    monkeypatch.setattr(fs, "_resolve_star_id", lambda s: (77777, None, "FakeCanon"))
    mag, name = fs.fixstar_mag("whatever")
    assert mag == pytest.approx(3.14)
    assert name == "FakeCanon"
