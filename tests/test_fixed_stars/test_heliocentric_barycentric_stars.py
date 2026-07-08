"""Heliocentric / barycentric fixed-star places (parity).

A heliocentric (``FLG_HELCTR``) or barycentric (``FLG_BARYCTR``) fixed-star
place is the star's astrometric position observed from the Sun / solar-system
barycentre, with NO annual aberration and NO gravitational deflection — the
reference API implies (and echoes) ``FLG_NOGDEFL | FLG_NOABERR`` for these
centers. Before this was implemented the star path ignored the center flags
and always returned the geocentric place.

The pinned longitudes/latitudes below are the ecliptic-of-date places under
that convention (they agree with the reference to <0.005"). The other tests
assert the structural properties: the implied retflag bits, the center
priority (TOPOCTR > BARYCTR > HELCTR), and cross-entry-point / cross-backend
consistency — so they hold on both the LEB and Skyfield backends the gates run.
"""

from __future__ import annotations

import pytest

import libephemeris as lib

# Ecliptic-of-date (lon, lat) in degrees for the Sun-centred (HELCTR) and
# SSB-centred (BARYCTR) astrometric place (no aberration, no deflection).
# Parity values: they match the reference apparent-place to <0.005".
_GOLDEN = {
    ("Regulus", 2451545.0, "HELCTR"): (149.8252659, 0.4648488),
    ("Regulus", 2451545.0, "BARYCTR"): (149.8252660, 0.4648488),
    ("Regulus", 2459015.5, "HELCTR"): (150.1086936, 0.4655100),
    ("Sirius", 2451545.0, "HELCTR"): (104.0777982, -39.6052375),
    ("Sirius", 2459015.5, "BARYCTR"): (104.3587180, -39.6099439),
    ("Aldebaran", 2451545.0, "HELCTR"): (69.7853208, -5.4673200),
    ("Aldebaran", 2451545.0, "BARYCTR"): (69.7853209, -5.4673200),
    ("Vega", 2459015.5, "HELCTR"): (285.5983810, 61.7318100),
}

_CENTER_FLAG = {"HELCTR": lib.FLG_HELCTR, "BARYCTR": lib.FLG_BARYCTR}


@pytest.mark.parametrize(("key", "expected"), list(_GOLDEN.items()))
def test_helctr_barctr_star_positions(key, expected) -> None:
    star, jd, center = key
    flags = lib.FLG_SWIEPH | _CENTER_FLAG[center]
    pos, _name, _ret = lib.fixstar_ut(star, jd, flags)
    # Separation tolerance: the LEB/Skyfield star pipeline floor is ~0.005".
    assert pos[0] == pytest.approx(expected[0], abs=2e-6)
    assert pos[1] == pytest.approx(expected[1], abs=2e-6)


@pytest.mark.parametrize("center", ["HELCTR", "BARYCTR"])
def test_helctr_barctr_retflag_implies_nogdefl_noaberr(center) -> None:
    # The UT entry points echo the reference's implied NOGDEFL|NOABERR bits
    # for heliocentric/barycentric output.
    flags = lib.FLG_SWIEPH | _CENTER_FLAG[center]
    _pos, _name, ret = lib.fixstar_ut("Regulus", 2451545.0, flags)
    assert ret & lib.FLG_NOGDEFL
    assert ret & lib.FLG_NOABERR


def test_center_priority_topoctr_over_bary_over_helctr() -> None:
    jd = 2451545.0
    lib.set_topo(12.5, 41.9, 100.0)
    # TOPOCTR wins over both BARYCTR and HELCTR (position + retflag).
    p_all, _n, r_all = lib.fixstar_ut(
        "Regulus", jd, lib.FLG_SWIEPH | lib.FLG_TOPOCTR | lib.FLG_BARYCTR | lib.FLG_HELCTR
    )
    p_topo, _n2, r_topo = lib.fixstar_ut("Regulus", jd, lib.FLG_SWIEPH | lib.FLG_TOPOCTR)
    assert r_all == r_topo
    assert p_all[0] == pytest.approx(p_topo[0], abs=1e-9)
    # BARYCTR wins over HELCTR.
    p_bh, _n3, r_bh = lib.fixstar_ut(
        "Regulus", jd, lib.FLG_SWIEPH | lib.FLG_BARYCTR | lib.FLG_HELCTR
    )
    p_b, _n4, r_b = lib.fixstar_ut("Regulus", jd, lib.FLG_SWIEPH | lib.FLG_BARYCTR)
    assert r_bh == r_b
    assert p_bh[0] == pytest.approx(p_b[0], abs=1e-9)


@pytest.mark.parametrize("center", ["HELCTR", "BARYCTR"])
def test_helctr_geocentric_differ_only_by_observer(center) -> None:
    # A heliocentric/barycentric place differs from the geocentric TRUEPOS
    # (also no aberration/deflection) only by the observer parallax, which for
    # these distant stars is well under an arcminute — but it is non-zero, so
    # the center flag is genuinely honoured (was previously ignored: identical).
    jd = 2451545.0
    p_geo, _n, _r = lib.fixstar_ut("Sirius", jd, lib.FLG_SWIEPH | lib.FLG_TRUEPOS)
    p_c, _n2, _r2 = lib.fixstar_ut("Sirius", jd, lib.FLG_SWIEPH | _CENTER_FLAG[center])
    dlon = abs((p_c[0] - p_geo[0] + 180.0) % 360.0 - 180.0) * 3600.0
    assert 0.0 < dlon < 60.0


@pytest.mark.parametrize("center", ["HELCTR", "BARYCTR"])
def test_batch_and_fixstar2_match_single(center) -> None:
    jd = 2451545.0
    flags = lib.FLG_SWIEPH | _CENTER_FLAG[center]
    stars = ("Regulus", "Sirius", "Aldebaran")
    batch = lib.batch_fixstars_ut(stars, jd, flags)
    for star, brow in zip(stars, batch):
        assert brow is not None
        single = lib.fixstar_ut(star, jd, flags)
        assert brow[0][0] == pytest.approx(single[0][0], rel=1e-12, abs=1e-9)
        assert brow[2] == single[2]
        # fixstar2_ut resolves the same star and honours the same center.
        two = lib.fixstar2_ut(star, jd, flags)
        assert two[0][0] == pytest.approx(single[0][0], rel=1e-12, abs=1e-9)
