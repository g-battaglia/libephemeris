"""Heliocentric/barycentric barycentre reports and osculating-element GM parity.

FIX B — heliocentric planet places:
  * the light-time is retarded over the BARYCENTRIC distance |planet - SSB|/c
    (not the heliocentric |planet - Sun|/c), which removes a ~0.5" longitude
    error for the inner planets (and the matching speed bias);
  * the outer planets are reported as their system BARYCENTRE (no
    centre-of-body offset), which removes a ~0.05" error (Pluto).
  The pinned longitudes are parity values (they match the reference to
  <0.001"). A regression of either half shifts them well beyond the tolerance.

FIX B2 — barycentric (FLG_BARYCTR) giant places:
  * the giants are likewise reported as their raw system BARYCENTRE (no COB)
    under FLG_BARYCTR, which removes up to ~0.09" (Pluto) / ~0.05" (Saturn).
    The BARYCTR light-time already runs on the barycentric distance, so only
    the COB is skipped. The pinned longitudes are parity values (<0.001").

FIX C — osculating orbital elements:
  * the two-body mu is G(M_sun + M_planet_system) (not GM_sun alone), which
    tightens the giants' semi-major axis / aphelion from ~1e-2 AU to <=~2.5e-5;
  * Earth's heliocentric orbit is derived from the Earth-Moon barycentre.
"""

from __future__ import annotations

import pytest

import libephemeris as lib

_J2000 = 2451545.0

# FIX B: heliocentric ecliptic (lon, lat) in degrees at J2000 — the
# barycentre place with barycentric light-time (parity to <0.001").
_HELCTR_LONLAT = {
    "MERCURY": (253.7715869, -3.0219453),
    "JUPITER": (36.2881236, -1.1746124),
    "SATURN": (45.7164531, -2.3032261),
    "PLUTO": (250.5412181, 11.1617071),
}

# FIX B2: barycentric ecliptic (lon, lat) in degrees at J2000 — the giants
# reported as their raw system barycentre (no COB), parity to <0.001".
_BARYCTR_LONLAT = {
    "JUPITER": (36.3109046, -1.1739836),
    "SATURN": (45.7361844, -2.3036915),
    "URANUS": (316.3935418, -0.6843648),
    "NEPTUNE": (303.9097217, 0.2424320),
    "PLUTO": (250.5300157, 11.1602509),
}

# FIX C: osculating (a, q, Q) in AU at J2000 (parity values).
_ELEMENTS = {
    "JUPITER": (5.204266630, 4.950429163, 5.458104098),
    "SATURN": (9.582017171, 9.048074649, 10.115959693),
    "NEPTUNE": (30.103639649, 29.766023027, 30.441256271),
    "EARTH": (0.999996427, 0.983294125, 1.016698730),
}


def _ipl(name: str) -> int:
    return getattr(lib, name)


@pytest.mark.parametrize(("name", "lonlat"), list(_HELCTR_LONLAT.items()))
def test_helctr_position_is_barycentre_with_barycentric_lighttime(name, lonlat) -> None:
    pos, _ret = lib.calc(_J2000, _ipl(name), lib.FLG_SWIEPH | lib.FLG_HELCTR)
    # 5e-6 deg = 0.018": below it the Mercury barycentric-LT shift (~0.28")
    # and the outer-planet COB shift (~0.05" Pluto/Saturn, ~0.025" Jupiter)
    # would all fail, while the LEB/Skyfield floor (~0.003") stays inside.
    assert pos[0] == pytest.approx(lonlat[0], abs=5e-6)
    assert pos[1] == pytest.approx(lonlat[1], abs=5e-6)


@pytest.mark.parametrize("name", ["MERCURY", "JUPITER", "PLUTO"])
def test_helctr_speed_is_derivative_of_position(name) -> None:
    # The reported HELCTR longitude speed must be the exact derivative of the
    # reported HELCTR longitude (the barycentric light-time rate fix). The
    # reference's own speed is self-inconsistent here, so we compare against a
    # finite difference of the library's position, not the reference.
    flags = lib.FLG_SWIEPH | lib.FLG_HELCTR | lib.FLG_SPEED
    pos, _ret = lib.calc(_J2000, _ipl(name), flags)
    h = 1e-3
    p_lo, _ = lib.calc(_J2000 - h, _ipl(name), flags)
    p_hi, _ = lib.calc(_J2000 + h, _ipl(name), flags)
    fd = (p_hi[0] - p_lo[0]) / (2.0 * h)
    # 0.02"/day: finite-difference noise floor; a wrong light-time rate would
    # leave a ~0.5"/day Mercury bias.
    assert pos[3] * 3600.0 == pytest.approx(fd * 3600.0, abs=0.02)


@pytest.mark.parametrize(("name", "lonlat"), list(_BARYCTR_LONLAT.items()))
def test_baryctr_position_is_barycentre(name, lonlat) -> None:
    # The giants' FLG_BARYCTR place is their raw system barycentre (no COB).
    pos, _ret = lib.calc(_J2000, _ipl(name), lib.FLG_SWIEPH | lib.FLG_BARYCTR)
    # 5e-6 deg = 0.018": below it the skipped-COB shift (~0.09" Pluto,
    # ~0.05" Saturn, ~0.02" Jupiter) fails, while the LEB/Skyfield floor
    # (~0.003") stays inside.
    assert pos[0] == pytest.approx(lonlat[0], abs=5e-6)
    assert pos[1] == pytest.approx(lonlat[1], abs=5e-6)


def test_baryctr_truepos_matches_apparent_position() -> None:
    # The COB is skipped on both the TRUEPOS (no light-time) and apparent
    # BARYCTR paths, so Pluto's barycentre place differs only by the small
    # light-time shift — a stray COB on either path would open a ~0.09" gap.
    tp, _ = lib.calc(_J2000, lib.PLUTO, lib.FLG_SWIEPH | lib.FLG_BARYCTR | lib.FLG_TRUEPOS)
    ap, _ = lib.calc(_J2000, lib.PLUTO, lib.FLG_SWIEPH | lib.FLG_BARYCTR)
    # Light-time over ~32 AU moves Pluto's barycentric longitude by ~0.02 deg;
    # a COB left on only one path would add a further ~0.09" (2.5e-5 deg).
    assert abs(tp[0] - ap[0]) < 0.03


@pytest.mark.parametrize("name", ["JUPITER", "SATURN", "PLUTO"])
def test_baryctr_speed_is_derivative_of_position(name) -> None:
    # The reported BARYCTR longitude speed must be the exact derivative of the
    # reported (barycentre) longitude: skipping the position COB must also skip
    # the COB velocity rate so the two stay mutually consistent.
    flags = lib.FLG_SWIEPH | lib.FLG_BARYCTR | lib.FLG_SPEED
    pos, _ret = lib.calc(_J2000, _ipl(name), flags)
    h = 1e-3
    p_lo, _ = lib.calc(_J2000 - h, _ipl(name), flags)
    p_hi, _ = lib.calc(_J2000 + h, _ipl(name), flags)
    fd = (p_hi[0] - p_lo[0]) / (2.0 * h)
    # 0.02"/day: finite-difference noise floor; a stray COB velocity rate would
    # leave a ~0.07"/day Pluto bias.
    assert pos[3] * 3600.0 == pytest.approx(fd * 3600.0, abs=0.02)


def test_sun_heliocentric_stays_origin() -> None:
    # The Sun is excluded from the barycentric-LT set; its heliocentric place
    # must remain the origin (a regression would light-time-retard it away).
    pos, _ret = lib.calc(_J2000, lib.SUN, lib.FLG_SWIEPH | lib.FLG_HELCTR)
    assert pos[2] == pytest.approx(0.0, abs=1e-9)


@pytest.mark.parametrize(("name", "aqQ"), list(_ELEMENTS.items()))
def test_orbital_elements_gm_and_emb(name, aqQ) -> None:
    el = lib.get_orbital_elements(_J2000, _ipl(name), lib.FLG_SWIEPH | lib.FLG_HELCTR)
    a, q, Q = el[0], el[15], el[16]
    # Giants: GM = G(M_sun+M_planet) holds a/Q to ~2.5e-5; a GM_sun-only
    # regression biases them by ~1e-3..1e-2 AU. Earth: EMB target holds a/Q to
    # ~5e-7; the Earth-centre orbit wobbles by ~2..9e-4 AU.
    tol = 3e-5 if name != "EARTH" else 5e-6
    assert a == pytest.approx(aqQ[0], abs=tol)
    assert q == pytest.approx(aqQ[1], abs=tol)
    assert Q == pytest.approx(aqQ[2], abs=tol)
