"""Independent invariants for heliocentric/barycentric fixed-star places.

A Sun- or SSB-centered fixed-star result is an astrometric place: annual
aberration and gravitational deflection are not observer corrections for
those centers.  Tests below use geometry, finite differences, and cross-entry
point consistency; no numeric result was recorded from another ephemeris.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as lib

_CENTER_FLAG = {"HELCTR": lib.FLG_HELCTR, "BARYCTR": lib.FLG_BARYCTR}


@pytest.mark.parametrize(
    "star,jd,center",
    [
        ("Regulus", lib.J2000, "HELCTR"),
        ("Regulus", lib.J2000, "BARYCTR"),
        ("Regulus", lib.julday(2020, 6, 15, 0.0), "HELCTR"),
        ("Sirius", lib.J2000, "HELCTR"),
        ("Sirius", lib.julday(2020, 6, 15, 0.0), "BARYCTR"),
        ("Aldebaran", lib.J2000, "HELCTR"),
        ("Aldebaran", lib.J2000, "BARYCTR"),
        ("Vega", lib.julday(2020, 6, 15, 0.0), "HELCTR"),
    ],
)
def test_centered_star_positions_are_valid(star, jd, center) -> None:
    pos, name, _ret = lib.fixstar_ut(star, jd, lib.FLG_SWIEPH | _CENTER_FLAG[center])
    assert name.startswith(star)
    assert all(math.isfinite(value) for value in pos)
    assert 0.0 <= pos[0] < 360.0
    assert -90.0 <= pos[1] <= 90.0
    assert pos[2] > 0.0


@pytest.mark.parametrize("center", ["HELCTR", "BARYCTR"])
def test_center_retflag_implies_astrometric_corrections(center) -> None:
    flags = lib.FLG_SWIEPH | _CENTER_FLAG[center]
    _pos, _name, ret = lib.fixstar_ut("Regulus", lib.J2000, flags)
    assert ret & lib.FLG_NOGDEFL
    assert ret & lib.FLG_NOABERR


def test_center_priority_topoctr_over_bary_over_helctr() -> None:
    lib.set_topo(12.5, 41.9, 100.0)
    try:
        p_all, _n, r_all = lib.fixstar_ut(
            "Regulus",
            lib.J2000,
            lib.FLG_SWIEPH | lib.FLG_TOPOCTR | lib.FLG_BARYCTR | lib.FLG_HELCTR,
        )
        p_topo, _n2, r_topo = lib.fixstar_ut(
            "Regulus", lib.J2000, lib.FLG_SWIEPH | lib.FLG_TOPOCTR
        )
    finally:
        lib.set_topo(0.0, 0.0, 0.0)
    assert r_all == r_topo
    assert p_all == pytest.approx(p_topo, rel=1e-12, abs=1e-9)

    p_bh, _n3, r_bh = lib.fixstar_ut(
        "Regulus", lib.J2000, lib.FLG_SWIEPH | lib.FLG_BARYCTR | lib.FLG_HELCTR
    )
    p_b, _n4, r_b = lib.fixstar_ut(
        "Regulus", lib.J2000, lib.FLG_SWIEPH | lib.FLG_BARYCTR
    )
    assert r_bh == r_b
    assert p_bh == pytest.approx(p_b, rel=1e-12, abs=1e-9)


@pytest.mark.parametrize("center", ["HELCTR", "BARYCTR"])
def test_geocentric_difference_is_stellar_parallax_scale(center) -> None:
    p_geo, _n, _r = lib.fixstar_ut(
        "Sirius", lib.J2000, lib.FLG_SWIEPH | lib.FLG_TRUEPOS
    )
    p_center, _n2, _r2 = lib.fixstar_ut(
        "Sirius", lib.J2000, lib.FLG_SWIEPH | _CENTER_FLAG[center]
    )
    dlon_arcsec = abs((p_center[0] - p_geo[0] + 180.0) % 360.0 - 180.0) * 3600.0
    # A one-AU baseline applied to a stellar parallax is sub-arcsecond; the
    # generous bound only guards against accidentally ignoring/changing center.
    assert 0.0 < dlon_arcsec < 2.0


@pytest.mark.parametrize("center", ["HELCTR", "BARYCTR"])
def test_batch_and_fixstar2_match_single(center) -> None:
    flags = lib.FLG_SWIEPH | _CENTER_FLAG[center]
    stars = ("Regulus", "Sirius", "Aldebaran")
    batch = lib.batch_fixstars_ut(stars, lib.J2000, flags)
    for star, brow in zip(stars, batch):
        assert brow is not None
        single = lib.fixstar_ut(star, lib.J2000, flags)
        modern = lib.fixstar2_ut(star, lib.J2000, flags)
        assert brow[0] == pytest.approx(single[0], rel=1e-12, abs=1e-9)
        assert brow[2] == single[2]
        assert modern[0] == pytest.approx(single[0], rel=1e-12, abs=1e-9)


def _wrapped_delta(after: float, before: float) -> float:
    return (after - before + 180.0) % 360.0 - 180.0


class TestTopocentricStarSpeeds:
    """Topocentric speeds are derivatives including the diurnal observer term."""

    def test_speed_matches_position_finite_difference(self):
        flags = lib.FLG_SWIEPH | lib.FLG_TOPOCTR
        h = 0.0025
        lib.set_topo(12.5, 41.9, 100.0)
        try:
            before = lib.fixstar2_ut("Regulus", lib.J2000 - h, flags)[0]
            after = lib.fixstar2_ut("Regulus", lib.J2000 + h, flags)[0]
            position = lib.fixstar2_ut("Regulus", lib.J2000, flags | lib.FLG_SPEED)[0]
        finally:
            lib.set_topo(0.0, 0.0, 0.0)

        dlon = _wrapped_delta(after[0], before[0]) / (2.0 * h)
        dlat = (after[1] - before[1]) / (2.0 * h)
        assert position[3] == pytest.approx(dlon, abs=2e-7)
        assert position[4] == pytest.approx(dlat, abs=2e-7)

    def test_topocentric_speed_differs_from_geocentric(self):
        geo = lib.fixstar2_ut("Regulus", lib.J2000, lib.FLG_SWIEPH | lib.FLG_SPEED)[0]
        lib.set_topo(12.5, 41.9, 100.0)
        try:
            topo = lib.fixstar2_ut(
                "Regulus",
                lib.J2000,
                lib.FLG_SWIEPH | lib.FLG_TOPOCTR | lib.FLG_SPEED,
            )[0]
        finally:
            lib.set_topo(0.0, 0.0, 0.0)
        assert abs(topo[3] - geo[3]) > 1.0e-5
        assert abs(topo[4] - geo[4]) > 1.0e-5


class TestEtRetflagVerbatim:
    """TT entries echo requests; UT entries add documented implied bits."""

    def test_et_entries_echo_flags_verbatim(self):
        flags_to_check = (
            0,
            lib.FLG_J2000,
            lib.FLG_SPEED,
            lib.FLG_HELCTR,
            lib.FLG_TRUEPOS,
            lib.FLG_SWIEPH,
            lib.FLG_JPLEPH,
        )
        for flags in flags_to_check:
            assert lib.fixstar("Regulus", lib.J2000, flags)[2] == flags
            assert lib.fixstar2("Regulus", lib.J2000, flags)[2] == flags

    def test_et_sidereal_fixed_epoch_echoes_verbatim(self):
        lib.set_sid_mode(lib.SIDM_J2000, 0.0, 0.0)
        try:
            ret = lib.fixstar("Regulus", lib.J2000, lib.FLG_SIDEREAL)[2]
            assert ret == lib.FLG_SIDEREAL
        finally:
            lib.set_sid_mode(lib.SIDM_FAGAN_BRADLEY, 0.0, 0.0)

    def test_ut_entries_add_selection_and_implied_bits(self):
        base = lib.fixstar_ut("Regulus", lib.J2000, 0)[2]
        assert base & lib.FLG_SWIEPH

        j2000 = lib.fixstar_ut("Regulus", lib.J2000, lib.FLG_J2000)[2]
        assert j2000 & lib.FLG_SWIEPH
        assert j2000 & lib.FLG_J2000
        assert j2000 & lib.FLG_NONUT

        truepos = lib.fixstar2_ut("Regulus", lib.J2000, lib.FLG_TRUEPOS)[2]
        assert truepos & lib.FLG_SWIEPH
        assert truepos & lib.FLG_TRUEPOS
        assert truepos & lib.FLG_NOABERR
        assert truepos & lib.FLG_NOGDEFL
