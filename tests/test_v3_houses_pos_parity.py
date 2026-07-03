"""house_pos() parity at high latitudes for previously-broken house systems.

These tests pin ``house_pos`` for the house systems whose placement was wrong
before the v3 review fixes:

* Equal-based systems ``E`` / ``A`` / ``W`` / ``V`` used an ascendant without
  the eastern-horizon hemisphere correction, so beyond the polar circle the
  result was 180 deg (six houses) off.
* Topocentric ``T`` reused the Placidus semi-arc placement instead of the
  Polich-Page position-circle method.
* The generic fallback used ecliptic-longitude cusp interpolation for
  ``H`` / ``F`` / ``Y`` / ``I`` / ``i`` / ``J`` / ``S`` / ``U``, dropping the
  body's ecliptic latitude and (for ``H``) the azimuth-based placement.

The expected values below were produced by the reference ephemeris API and are
hardcoded so the test does not depend on it at run time. Inputs are chosen at
high latitude in both hemispheres and away from exact house cusps (the
reference snaps a body sitting exactly on a cusp with zero latitude to the
integer cusp number, a measure-zero tie-break the closed-form placement -- like
the pre-existing Placidus/Koch/Campanus branches -- does not reproduce).
"""

from __future__ import annotations

import importlib

import pytest

# ``libephemeris.houses`` the module is shadowed in the package namespace by the
# re-exported ``houses()`` function, so import the module object explicitly.
H = importlib.import_module("libephemeris.houses")


# (armc, geolat, eps, lon, lat_body)
CASES = [
    (123.0, 66.0, 23.4393, 47.3, 2.5),
    (271.0, 71.0, 23.4393, 251.9, -4.5),
    (34.0, -63.0, 23.4393, 133.6, 4.5),
    (317.0, 68.5, 23.4368, 299.3, -3.0),
    (200.0, -72.0, 23.4393, 88.1, 0.0),
    (95.0, 74.0, 23.437, 331.7, 5.0),
]

# Ground-truth house positions from the reference ephemeris API (hardcoded).
EXPECTED = {
    "E": [
        7.9721676638,
        9.5365316027,
        3.0497092808,
        7.0533076763,
        4.9841248455,
        5.9843637673,
    ],
    "A": [
        7.9721676638,
        9.5365316027,
        3.0497092808,
        7.0533076763,
        4.9841248455,
        5.9843637673,
    ],
    "W": [
        8.5766666759,
        10.3966666759,
        3.4533333426,
        7.9766666759,
        5.9366666759,
        6.0566666759,
    ],
    "V": [
        8.4721676638,
        10.0365316027,
        3.5497092808,
        7.5533076763,
        5.4841248455,
        6.4843637673,
    ],
    "T": [
        8.2724705338,
        6.0382533669,
        2.2953538895,
        6.7335947752,
        5.1949177682,
        5.5270591974,
    ],
    "H": [
        7.0632183236,
        9.3613830185,
        12.9362441805,
        9.5448839052,
        7.4519223097,
        5.9812613942,
    ],
    "F": [
        7.9120787544,
        9.4523642355,
        3.2277053808,
        7.078290528,
        4.9067434921,
        5.9975955476,
    ],
    "Y": [
        8.6355714554,
        5.1794996562,
        2.3506957109,
        4.8575115367,
        4.6300080994,
        5.4261273727,
    ],
    "I": [
        8.5016546343,
        5.1514604727,
        2.5148310268,
        6.146154829,
        4.9488693128,
        5.4064545708,
    ],
    "i": [
        8.5016546343,
        5.1514604727,
        2.5148310268,
        6.146154829,
        4.9488693128,
        5.4064545708,
    ],
    "J": [
        7.9157142377,
        5.9927105303,
        2.0415926809,
        6.5956903901,
        5.8327082027,
        6.3427522129,
    ],
    "S": [
        8.3523437436,
        9.9001637764,
        2.7791706797,
        7.7849967281,
        6.0707264955,
        6.4562969293,
    ],
    "U": [
        8.6269599269,
        9.765813081,
        3.0858483526,
        7.7394760922,
        4.6072820966,
        4.8249208222,
    ],
}

TOL = 1e-6


def _circular_house_diff(a: float, b: float) -> float:
    """Absolute house-position difference on the wrap-around circle [1, 13)."""
    d = abs(a - b)
    return min(d, 12.0 - d)


@pytest.mark.parametrize("hsys", sorted(EXPECTED))
def test_house_pos_matches_reference_high_latitude(hsys):
    """Every fixed system matches the reference oracle at high latitude."""
    expected = EXPECTED[hsys]
    for (armc, geolat, eps, lon, lat_body), exp in zip(CASES, expected):
        got = H.house_pos(armc, geolat, eps, [lon, lat_body], hsys.encode())
        assert _circular_house_diff(got, exp) < TOL, (
            f"hsys={hsys} armc={armc} geolat={geolat} eps={eps} "
            f"lon={lon} lat={lat_body}: got {got}, expected {exp}"
        )


@pytest.mark.parametrize("hsys", ["E", "A", "W", "V"])
def test_equal_systems_polar_hemisphere_correction(hsys):
    """Equal systems apply the eastern-horizon fix (no six-house polar error).

    Before the fix these returned a value ~6 houses away from the oracle
    beyond the polar circle (armc=330, geolat=75).
    """
    got = H.house_pos(330.0, 75.0, 23.4393, [300.0, 0.0], hsys.encode())
    # Oracle values for this configuration (hardcoded).
    expected = {"E": 6.339098, "A": 6.339098, "W": 7.0, "V": 6.839098}[hsys]
    assert _circular_house_diff(got, expected) < 1e-4


def test_topocentric_uses_polich_page_not_placidus():
    """Topocentric 'T' now differs from Placidus and matches the oracle.

    At high latitude the Polich-Page placement is not identical to the
    Placidus semi-arc placement; the old code used the latter for 'T'.
    """
    armc, geolat, eps, lon, lat_body = 271.0, 71.0, 23.4393, 251.9, -4.5
    t = H.house_pos(armc, geolat, eps, [lon, lat_body], b"T")
    p = H.house_pos(armc, geolat, eps, [lon, lat_body], b"P")
    assert _circular_house_diff(t, p) > 0.1  # genuinely different systems
    assert _circular_house_diff(t, 6.0382533669) < TOL  # matches the oracle


def test_horizon_uses_azimuth_depends_on_latitude():
    """Horizontal 'H' is azimuth-based, so it depends on ecliptic latitude."""
    armc, geolat, eps, lon = 123.0, 66.0, 23.4393, 47.3
    h_north = H.house_pos(armc, geolat, eps, [lon, 5.0], b"H")
    h_zero = H.house_pos(armc, geolat, eps, [lon, 0.0], b"H")
    h_south = H.house_pos(armc, geolat, eps, [lon, -5.0], b"H")
    # A pure ecliptic-longitude interpolation would be identical for all three.
    assert h_north != h_zero
    assert h_south != h_zero


def test_sunshine_and_apc_depend_on_latitude():
    """Sunshine 'I'/'i' and APC 'Y' use the body's equatorial position."""
    armc, geolat, eps, lon = 271.0, 71.0, 23.4393, 251.9
    for hsys in ("I", "i", "Y"):
        top = H.house_pos(armc, geolat, eps, [lon, 6.0], hsys.encode())
        bot = H.house_pos(armc, geolat, eps, [lon, -6.0], hsys.encode())
        assert top != bot, f"{hsys} should depend on ecliptic latitude"
