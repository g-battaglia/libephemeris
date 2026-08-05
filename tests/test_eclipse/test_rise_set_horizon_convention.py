"""What a rise/set event actually IS, measured with refraction left ON.

The neighbouring `test_rise_set_independent.py` adjudicates rise/set against an
erfa observed-altitude reference, but it deliberately passes
``BIT_NO_REFRACTION`` "so lib and erfa share the same altitude condition — pure
geometry". That makes it blind by construction to everything this module cares
about: the refraction model, the disc-limb offset, the topocentric parallax and
the observer-altitude convention are exactly the terms it switches off. A
17-second temperature-default change or a dropped parallax would sail through it.

So this module asks a different question. It does not compare our event time to
someone else's event time — a comparison that is ill-conditioned near the poles,
where the Sun grazes the horizon and a milliarcsecond of altitude becomes minutes
of clock. It takes the instant we return and asks an INDEPENDENT computation
where the Sun actually was. That is well conditioned everywhere, it is stated in
arcminutes rather than seconds, and the answer must be the same number over the
whole globe:

    at every sunrise and sunset we return, the Sun's true topocentric UPPER LIMB
    sits on the refracted horizon.

Independence comes from erfa (IAU SOFA) for the sidereal rotation, fed the
library's own topocentric apparent RA/Dec — the same split the Moon case in
`test_rise_set_independent.py` already uses, and the one that keeps the ~8.8"
solar diurnal parallax in the measurement instead of quietly dropping it.

Network-free; needs pyerfa.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as L
from libephemeris.constants import (
    SUN,
    FLG_SWIEPH,
    FLG_EQUATORIAL,
    FLG_TOPOCTR,
    CALC_RISE,
    CALC_SET,
    BIT_DISC_CENTER,
    BIT_NO_REFRACTION,
)

erfa = pytest.importorskip("erfa")

_AU_KM = 149597870.7
_RSUN_KM = 696000.0

#: Standard atmosphere. Passed EXPLICITLY everywhere below, because the library
#: default is a literal 0 C (see `test_the_zero_default_is_zero_celsius`).
_STD_PRESSURE = 1013.25
_STD_TEMPERATURE = 15.0

#: Where the Sun's true upper limb sits at the event under the standard
#: atmosphere: minus the refraction at the apparent horizon. Measured across the
#: grid below, the spread is 0.0069' (0.42"), so this is an invariant and not an
#: average. The tolerance is ~4x that spread — tight enough that dropping the
#: solar parallax (0.147') or switching to a fixed 16' semidiameter would fail it,
#: loose enough that ordinary ephemeris noise cannot.
_LIMB_ON_HORIZON_ARCMIN = -33.590
_LIMB_TOL_ARCMIN = 0.03

#: Same quantity with the library default atmosphere (0 C): colder air refracts
#: more, so the Sun has further to climb.
_LIMB_AT_ZERO_C_ARCMIN = -36.736

#: Equator to just inside the polar circle, both hemispheres, all four seasons.
#: The point is that ONE number has to hold across all of it.
_GRID = [
    ("rome_equinox", 12.5, 41.9, 2024, 3, 20),
    ("rome_solstice", 12.5, 41.9, 2024, 6, 21),
    ("equator", 0.0, 0.0, 2024, 9, 22),
    ("sydney", 151.2, -33.9, 2024, 12, 21),
    ("reykjavik", -21.9, 64.1, 2024, 4, 15),
    ("nairobi", 36.8, -1.3, 2025, 1, 15),
    ("ushuaia", -68.3, -54.8, 2024, 9, 22),
    ("singapore", 103.8, 1.35, 2025, 8, 5),
]


def _true_altitude_arcmin(jd_ut: float, lon: float, lat: float) -> float:
    """Sun's true topocentric altitude, rotated by erfa rather than by us.

    The position is the library's own topocentric apparent RA/Dec (validated
    elsewhere); the hour angle comes from erfa's GST. Asking for FLG_TOPOCTR is
    what keeps the diurnal parallax in — the term an ICRS-fed ``erfa.atco13``
    with ``px=0`` would silently discard, moving the answer by 8.8".
    """
    L.set_topo(lon, lat, 0.0)
    radec, _ = L.calc_ut(jd_ut, SUN, FLG_SWIEPH | FLG_EQUATORIAL | FLG_TOPOCTR)
    ra, dec = math.radians(radec[0]), math.radians(radec[1])
    gst = erfa.gst06a(2400000.5, jd_ut - 2400000.5, 2400000.5, jd_ut - 2400000.5)
    hour_angle = gst + math.radians(lon) - ra
    phi = math.radians(lat)
    altitude = math.asin(
        math.sin(dec) * math.sin(phi)
        + math.cos(dec) * math.cos(phi) * math.cos(hour_angle)
    )
    return math.degrees(altitude) * 60.0


def _semidiameter_arcmin(jd_ut: float, lon: float, lat: float) -> float:
    """Topocentric solar semidiameter, from the actual distance of the day."""
    L.set_topo(lon, lat, 0.0)
    position, _ = L.calc_ut(jd_ut, SUN, FLG_SWIEPH | FLG_TOPOCTR)
    return math.degrees(math.asin(_RSUN_KM / (position[2] * _AU_KM))) * 60.0


def _event(jd0, lon, lat, rsmi, atpress=_STD_PRESSURE, attemp=_STD_TEMPERATURE, alt=0.0):
    retflag, tret = L.rise_trans(jd0, SUN, rsmi, (lon, lat, alt), atpress, attemp)
    assert retflag == 0, f"expected an event, got retflag={retflag}"
    return tret[0]


@pytest.mark.parametrize("name,lon,lat,y,m,d", _GRID, ids=lambda v: v if isinstance(v, str) else "")
@pytest.mark.parametrize("rsmi,label", [(CALC_RISE, "rise"), (CALC_SET, "set")])
def test_the_upper_limb_sits_on_the_refracted_horizon(name, lon, lat, y, m, d, rsmi, label):
    """The one invariant: independent of latitude, season and hemisphere."""
    jd0 = L.julday(y, m, d, 0.0, L.GREG_CAL)
    jd = _event(jd0, lon, lat, rsmi)
    limb = _true_altitude_arcmin(jd, lon, lat) + _semidiameter_arcmin(jd, lon, lat)
    assert abs(limb - _LIMB_ON_HORIZON_ARCMIN) < _LIMB_TOL_ARCMIN, (
        f"{name} {label}: upper limb at {limb:.4f}', expected "
        f"{_LIMB_ON_HORIZON_ARCMIN:.4f}' +/- {_LIMB_TOL_ARCMIN}'. The horizon "
        f"convention moved — check the refraction model, the semidiameter term "
        f"or the topocentric parallax before regenerating anything."
    )


@pytest.mark.parametrize("name,lon,lat,y,m,d", _GRID, ids=lambda v: v if isinstance(v, str) else "")
def test_disc_center_moves_the_target_from_the_limb_to_the_centre(name, lon, lat, y, m, d):
    """BIT_DISC_CENTER drops the semidiameter and nothing else.

    Guards the flag against being wired to the wrong sign, which would look
    correct in an ordering test and be a full solar diameter out here.
    """
    jd0 = L.julday(y, m, d, 0.0, L.GREG_CAL)
    jd = _event(jd0, lon, lat, CALC_RISE | BIT_DISC_CENTER)
    centre = _true_altitude_arcmin(jd, lon, lat)
    assert abs(centre - _LIMB_ON_HORIZON_ARCMIN) < _LIMB_TOL_ARCMIN, (
        f"{name}: with BIT_DISC_CENTER the CENTRE should sit where the limb "
        f"otherwise does ({_LIMB_ON_HORIZON_ARCMIN:.4f}'), got {centre:.4f}'"
    )


@pytest.mark.parametrize("name,lon,lat,y,m,d", _GRID, ids=lambda v: v if isinstance(v, str) else "")
def test_the_geometric_horizon_is_exactly_zero(name, lon, lat, y, m, d):
    """Disc centre with refraction off must land on 0.000', not merely near it.

    This is the tightest statement the suite can make, and it is the one that
    would catch the search itself drifting: with both physical terms removed the
    only thing left is the root-finder and the position pipeline.
    """
    jd0 = L.julday(y, m, d, 0.0, L.GREG_CAL)
    jd = _event(jd0, lon, lat, CALC_RISE | BIT_DISC_CENTER | BIT_NO_REFRACTION)
    centre = _true_altitude_arcmin(jd, lon, lat)
    assert abs(centre) < 0.002, f"{name}: geometric horizon crossing at {centre:.5f}', expected 0"


def test_the_zero_default_is_zero_celsius_not_the_standard_atmosphere():
    """attemp=0 means 0 C, as in the reference — not an implicit 15 C.

    Worth its own test because it is the cheapest mistake to make in a caller and
    the hardest to notice: the event still looks perfectly reasonable, it is just
    17 seconds off a standard-atmosphere convention.
    """
    jd0 = L.julday(2024, 3, 20, 0.0, L.GREG_CAL)
    lon, lat = 12.5, 41.9

    default_jd = _event(jd0, lon, lat, CALC_RISE, atpress=0.0, attemp=0.0)
    standard_jd = _event(jd0, lon, lat, CALC_RISE)

    limb = _true_altitude_arcmin(default_jd, lon, lat) + _semidiameter_arcmin(default_jd, lon, lat)
    assert abs(limb - _LIMB_AT_ZERO_C_ARCMIN) < _LIMB_TOL_ARCMIN, (
        f"default-atmosphere limb at {limb:.4f}', expected {_LIMB_AT_ZERO_C_ARCMIN:.4f}' (0 C)"
    )

    seconds = (standard_jd - default_jd) * 86400.0
    assert 15.0 < seconds < 19.0, (
        f"0 C -> 15 C should move sunrise by about 17 s, got {seconds:.1f} s"
    )


def test_observer_elevation_does_not_lower_the_horizon():
    """geopos[2] never depresses the horizon, whatever it appears to do.

    Three statements, because the naive one is wrong in an instructive way:

    1. With an EXPLICIT pressure, elevation changes the event by nothing at all
       (only the parallax sees it, and that is far below a second of time).
    2. With the default atpress=0 the elevation is what the barometric estimate
       reads, so a higher observer gets thinner air, less refraction and a LATER
       sunrise — the opposite of the physical intuition, and the reason this
       looks like a bug when it is a convention. rise_trans keeps the level sea
       horizon at any elevation, which is exactly what published rise/set tables
       do.
    3. Lowering the horizon is opt-in through horhgt=-100, and it is worth
       minutes, not seconds.
    """
    jd0 = L.julday(2024, 3, 20, 0.0, L.GREG_CAL)
    lon, lat = 12.5, 41.9

    # (1) explicit pressure: elevation is inert
    fixed_sea = _event(jd0, lon, lat, CALC_RISE)
    fixed_high = _event(jd0, lon, lat, CALC_RISE, alt=1000.0)
    assert abs(fixed_high - fixed_sea) * 86400.0 < 0.5, (
        "with an explicit atpress the observer altitude must not move the event"
    )

    # (2) estimated pressure: elevation thins the air, so the event moves LATER
    est_sea = _event(jd0, lon, lat, CALC_RISE, atpress=0.0, attemp=0.0)
    est_high = _event(jd0, lon, lat, CALC_RISE, atpress=0.0, attemp=0.0, alt=1000.0)
    later_seconds = (est_high - est_sea) * 86400.0
    assert 18.0 < later_seconds < 30.0, (
        f"1000 m with an estimated pressure should push sunrise ~24 s LATER, got {later_seconds:+.1f} s"
    )

    # (3) the dip is the only thing that actually lowers the horizon
    retflag, tret = L.rise_trans_true_hor(
        jd0, SUN, CALC_RISE, (lon, lat, 1000.0), _STD_PRESSURE, _STD_TEMPERATURE, -100.0
    )
    assert retflag == 0
    dip_seconds = (tret[0] - fixed_high) * 86400.0
    assert -400.0 < dip_seconds < -320.0, (
        f"horhgt=-100 at 1000 m should pull sunrise ~360 s EARLIER, got {dip_seconds:+.1f} s"
    )
