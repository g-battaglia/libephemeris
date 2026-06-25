"""Regression guards for calc_pctr / nod_aps flag handling (issue #34 + audit).

These cover the gaps that the audit in ``scripts/validate_nonut_nodaps.py``
surfaced and fixed, in BOTH the LEB and Skyfield backends:

1. ``calc_pctr(X, Earth, flags) == calc(X, flags)`` — the planet-centric path
   must reduce exactly to the geocentric apparent position. This guards the
   missing-stellar-aberration bug (calc_pctr was ~20.5" off the geocentric
   apparent Moon) AND the FLG_NONUT handling (the equatorial/ecliptic mean-frame
   split now mirrors ``_calc_body``).

2. ``nod_aps`` honors FLG_NONUT without changing the (nutation-invariant) node
   longitude, matching the reference ephemeris.

(The sidereal+equatorial *speed* was investigated in the same audit and left
unchanged: the reference ephemeris computes that speed in the tropical/true-equator frame --
the same as the plain equatorial speed -- even though the sidereal position uses
the mean equator. libephemeris already matched that, so there is nothing new to
guard there.)
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import (
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    EARTH,
    FLG_EQUATORIAL,
    FLG_J2000,
    FLG_NONUT,
    NODBIT_OSCU,
)

from .conftest import CompareHelper, lon_error_arcsec

_PCTR_BODIES = [(MOON, "Moon"), (MERCURY, "Mercury"), (VENUS, "Venus"), (MARS, "Mars")]
_PCTR_FLAGS = [
    (0, "ecliptic"),
    (FLG_NONUT, "ecliptic+NONUT"),
    (FLG_EQUATORIAL, "equatorial"),
    (FLG_EQUATORIAL | FLG_NONUT, "equatorial+NONUT"),
    (FLG_J2000, "J2000"),
]


@pytest.mark.leb_compare
@pytest.mark.parametrize("body,bname", _PCTR_BODIES)
@pytest.mark.parametrize("flags,fname", _PCTR_FLAGS)
def test_calc_pctr_earth_reduces_to_geocentric(
    compare: CompareHelper, body: int, bname: str, flags: int, fname: str
):
    """calc_pctr(X, Earth, flags) must equal calc(X, flags) in each backend.

    Catches (a) the missing stellar aberration in the planet-centric path and
    (b) FLG_NONUT not being honored there -- both would break this identity.
    """
    jds = [ephem.julday(y, 6, 15, 12.0) for y in (1700, 1850, 1950, 2000, 2100, 2300)]
    for mode_name, run in (("skyfield", compare.skyfield), ("leb", compare.leb)):
        for jd in jds:
            pctr, _ = run(ephem.calc, jd, body, flags)  # warm any state first
            pctr, _ = run(ephem.calc_pctr, jd, body, EARTH, flags)
            geo, _ = run(ephem.calc, jd, body, flags)
            err = lon_error_arcsec(pctr[0], geo[0])
            err_lat = abs(pctr[1] - geo[1]) * 3600.0
            assert err < 0.1 and err_lat < 0.1, (
                f"{mode_name} calc_pctr({bname},Earth,{fname}) != calc: "
                f'Δlon={err:.4f}" Δlat={err_lat:.4f}" @JD{jd:.0f}'
            )


@pytest.mark.leb_compare
@pytest.mark.parametrize("body,bname", [(MARS, "Mars"), (JUPITER, "Jupiter")])
def test_nod_aps_nonut_node_invariant(compare: CompareHelper, body: int, bname: str):
    """nod_aps honors FLG_NONUT and the node longitude stays nutation-invariant.

    The osculating node longitude does not carry nutation (verified vs
    the reference ephemeris), so NONUT must run cleanly and leave the node unchanged.
    """
    jds = [ephem.julday(y, 1, 1, 12.0) for y in (1900, 2000, 2100)]
    for mode_name, run in (("skyfield", compare.skyfield), ("leb", compare.leb)):
        for jd in jds:
            base = run(ephem.nod_aps, jd, body, 0, NODBIT_OSCU)
            nonut = run(ephem.nod_aps, jd, body, FLG_NONUT, NODBIT_OSCU)
            for idx in (0, 1):  # ascending, descending node
                err = lon_error_arcsec(base[idx][0], nonut[idx][0])
                assert err < 0.01, (
                    f"{mode_name} nod_aps({bname}) node[{idx}] moved under NONUT: "
                    f'Δ={err:.4f}" @JD{jd:.0f}'
                )
