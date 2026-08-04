"""Native-float contract + global-state isolation guards.

Two cross-cutting invariants that no single feature test covers, and whose
violations are silent (a numpy-typed return serialises fine; a leaked global
setting only shows as a wrong number in an unrelated later call):

  * every public position component is a native Python ``float`` (a numpy-float
    leak was found before in the fixstar path);
  * the global setters (set_topo / set_sid_mode / set_calc_mode) scope their
    effect — a geocentric call is unaffected by a topo change, a tropical call by
    a sidereal mode, and switching backend/ayanamsha and back is reproducible.
"""

from __future__ import annotations

import numpy as np
import pytest

import libephemeris as L
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    PLUTO,
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    FLG_SWIEPH,
    FLG_SPEED,
    FLG_EQUATORIAL,
    FLG_XYZ,
    FLG_HELCTR,
    FLG_J2000,
    FLG_SIDEREAL,
    FLG_TOPOCTR,
)

JD = 2451545.0
CERES = 10000 + 1
CHIRON = 10000 + 2060

_BODIES = [SUN, MOON, MARS, PLUTO, MEAN_NODE, TRUE_NODE, MEAN_APOG, CERES, CHIRON]
_FLAGS = [
    FLG_SWIEPH,
    FLG_SWIEPH | FLG_SPEED,
    FLG_SWIEPH | FLG_EQUATORIAL,
    FLG_SWIEPH | FLG_XYZ,
    FLG_SWIEPH | FLG_HELCTR,
    FLG_SWIEPH | FLG_J2000,
]


@pytest.fixture(autouse=True)
def _skyfield_mode():
    L.set_calc_mode("skyfield")
    yield
    L.set_calc_mode("auto")
    L.set_topo(0.0, 0.0, 0.0)


@pytest.mark.parametrize("body", _BODIES)
@pytest.mark.parametrize("flags", _FLAGS)
def test_calc_returns_native_floats(body, flags):
    """Every component of calc/calc_ut must be a native float (no numpy leak)."""
    for fn in (L.calc, L.calc_ut):
        try:
            comps, _ = fn(JD, body, flags)
        except Exception:  # body/flag combo unsupported — not a contract case
            continue
        for v in comps:
            assert type(v) is float, f"{fn.__name__}({body},{flags}): {type(v)}"
            assert not isinstance(v, np.generic)


@pytest.mark.parametrize("body", _BODIES)
def test_pheno_returns_native_floats(body):
    """Every pheno/pheno_ut slot must be a native float.

    Guards the Moon apparent-diameter (slot 3) numpy ``float64`` leak on the
    skyfield path, where ``moon_dist.au`` flowed into ``_calc_apparent_diameter``
    without a ``float()`` cast.
    """
    for fn in (L.pheno, L.pheno_ut):
        try:
            result = fn(JD, body, FLG_SWIEPH)
        except Exception:  # body unsupported by pheno — not a contract case
            continue
        for i, v in enumerate(result):
            assert type(v) is float, f"{fn.__name__}({body})[{i}]: {type(v)}"
            assert not isinstance(v, np.generic)


def test_geocentric_unaffected_by_topo():
    L.set_topo(0.0, 0.0, 0.0)
    geo = L.calc_ut(JD, MOON, FLG_SWIEPH)[0][0]
    L.set_topo(100.0, 60.0, 3000.0)  # change observer
    geo2 = L.calc_ut(JD, MOON, FLG_SWIEPH)[0][0]  # still geocentric (no TOPOCTR)
    assert abs(geo - geo2) < 1e-9


def test_topocentric_actually_differs():
    """Sanity: with FLG_TOPOCTR the observer DOES matter (guards the above)."""
    L.set_topo(0.0, 0.0, 0.0)
    a = L.calc_ut(JD, MOON, FLG_SWIEPH | FLG_TOPOCTR)[0][0]
    L.set_topo(0.0, 80.0, 0.0)
    b = L.calc_ut(JD, MOON, FLG_SWIEPH | FLG_TOPOCTR)[0][0]
    assert abs(a - b) * 3600.0 > 1.0  # parallax shifts the Moon measurably


def test_sidereal_mode_isolation():
    L.set_sid_mode(L.SIDM_J2000)
    j2000 = L.calc_ut(JD, MOON, FLG_SWIEPH | FLG_SIDEREAL)[0][0]
    L.set_sid_mode(L.SIDM_J1900)
    j1900 = L.calc_ut(JD, MOON, FLG_SWIEPH | FLG_SIDEREAL)[0][0]
    L.set_sid_mode(L.SIDM_J2000)
    j2000_again = L.calc_ut(JD, MOON, FLG_SWIEPH | FLG_SIDEREAL)[0][0]
    assert abs(j2000 - j2000_again) < 1e-9, (
        "stale ayanamsha state after switching modes"
    )
    assert abs((j2000 - j1900 + 180) % 360 - 180) * 3600.0 > 100.0


def test_tropical_not_polluted_by_sidereal_mode():
    geo = L.calc_ut(JD, MOON, FLG_SWIEPH)[0][0]
    L.set_sid_mode(1, 0.0, 0.0)
    trop = L.calc_ut(JD, MOON, FLG_SWIEPH)[0][0]  # no SIDEREAL flag
    assert abs(geo - trop) < 1e-9, "sidereal mode leaked into a tropical call"


def test_calc_mode_switch_reproducible():
    sky = L.calc_ut(JD, MARS, FLG_SWIEPH)[0][0]
    L.set_calc_mode("leb")
    try:
        L.calc_ut(JD, MARS, FLG_SWIEPH)
    except Exception:
        pass
    L.set_calc_mode("skyfield")
    sky2 = L.calc_ut(JD, MARS, FLG_SWIEPH)[0][0]
    assert abs(sky - sky2) < 1e-9
