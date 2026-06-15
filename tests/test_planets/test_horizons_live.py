"""Live (network) numerical validation of the Horizons backend.

The other Horizons tests use full mocks (project policy, no network in CI). Those
exercise the plumbing but never check that the *numbers* are right — exactly the
gap that let a different backend (ASSIST) ship a wrong-direction bug. These tests
opt into the network (``@pytest.mark.network``, deselected by ``-m "not network"``
and skipped when offline) and assert that live Horizons output agrees numerically
with the independent Skyfield/DE440 backend for real bodies.
"""
from __future__ import annotations


import pytest

import libephemeris as L
from libephemeris.constants import (
    MARS,
    JUPITER,
    SATURN,
    MOON,
    CERES,
    FLG_SWIEPH,
    FLG_EQUATORIAL,
    FLG_HELCTR,
    FLG_J2000,
    FLG_TRUEPOS,
)

J2000 = 2451545.0


def _ang_diff(a: float, b: float) -> float:
    d = abs(a - b) % 360.0
    return min(d, 360.0 - d)


@pytest.fixture
def horizons_vs_skyfield():
    """Return a comparator (jd_tt, body, flags) -> (horizons, skyfield) tuples,
    or skip the test if Horizons can't be reached."""
    def _both(jd_tt, body, flags):
        L.set_calc_mode("skyfield")
        sky = L.calc(jd_tt, body, flags)[0]
        L.set_calc_mode("horizons")
        try:
            hor = L.calc(jd_tt, body, flags)[0]
        except (ConnectionError, OSError) as exc:  # offline / JPL unreachable
            pytest.skip(f"Horizons not reachable: {exc}")
        return hor, sky

    try:
        yield _both
    finally:
        L.set_calc_mode("auto")


@pytest.mark.network
def test_horizons_apparent_ecliptic_matches_skyfield(horizons_vs_skyfield):
    """Geocentric apparent ecliptic longitude/latitude from Horizons must match
    the DE440/Skyfield backend to < 1" (both are full apparent-place pipelines)."""
    for body in (MARS, JUPITER, SATURN, MOON):
        hor, sky = horizons_vs_skyfield(J2000, body, FLG_SWIEPH)
        dlon = _ang_diff(hor[0], sky[0]) * 3600.0
        dlat = abs(hor[1] - sky[1]) * 3600.0
        ddist = abs(hor[2] - sky[2])
        assert dlon < 1.0, f"body {body}: Horizons-Skyfield dlon={dlon:.3f}\""
        assert dlat < 1.0, f"body {body}: Horizons-Skyfield dlat={dlat:.3f}\""
        assert ddist < 1e-5, f"body {body}: dist diff {ddist:.2e} AU"


@pytest.mark.network
def test_horizons_equatorial_matches_skyfield(horizons_vs_skyfield):
    """RA/Dec (FLG_EQUATORIAL) from Horizons must match Skyfield to < 1"."""
    for body in (MARS, JUPITER):
        hor, sky = horizons_vs_skyfield(J2000, body, FLG_SWIEPH | FLG_EQUATORIAL)
        assert _ang_diff(hor[0], sky[0]) * 3600.0 < 1.0
        assert abs(hor[1] - sky[1]) * 3600.0 < 1.0


@pytest.mark.network
def test_horizons_asteroid_heliocentric_matches_skyfield(horizons_vs_skyfield):
    """A small body (Ceres) heliocentric geometric position from Horizons must
    match the Skyfield/SPK backend to < 1" — the regression guard for the kind of
    wrong-frame / wrong-direction backend bug that motivated this file."""
    flags = FLG_SWIEPH | FLG_HELCTR | FLG_J2000 | FLG_TRUEPOS
    hor, sky = horizons_vs_skyfield(J2000, CERES, flags)
    assert _ang_diff(hor[0], sky[0]) * 3600.0 < 1.0
    assert abs(hor[2] - sky[2]) < 1e-4
    # Heliocentric distance must be physically sane for Ceres (~2.6-3.0 AU),
    # a coarse direction/sanity invariant independent of Skyfield.
    assert 2.4 < hor[2] < 3.1


@pytest.mark.network
def test_horizons_returns_native_floats(horizons_vs_skyfield):
    """Horizons output must be native Python floats (no numpy leak)."""
    hor, _ = horizons_vs_skyfield(J2000, MARS, FLG_SWIEPH)
    for v in hor:
        assert type(v) is float
