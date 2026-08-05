"""What a rise/set event actually IS, measured with refraction left ON.

The neighbouring `test_rise_set_independent.py` adjudicates rise/set against an
erfa observed-altitude reference, but it deliberately passes
``BIT_NO_REFRACTION`` "so lib and erfa share the same altitude condition". That
makes it blind by construction to everything this module cares about: the
refraction model, the disc-limb offset and the topocentric parallax are exactly
the terms it switches off.

So this module asks a different question. It does not compare our event time to
someone else's event time — a comparison that is ill-conditioned near the poles,
where the Sun grazes the horizon and a milliarcsecond of altitude becomes minutes
of clock. It takes the instant we return and asks where the Sun actually was.
That is well conditioned everywhere, it is stated in arcminutes rather than
seconds, and the answer must be the same number over the whole globe:

    at every sunrise and sunset we return, the Sun's true UPPER LIMB sits on the
    refracted horizon.

WHY THE ORACLE IS SKYFIELD AND NOT ERFA
---------------------------------------
The first version of this module measured the altitude from the library's own
``calc_ut(FLG_TOPOCTR)`` and rotated it with erfa, and claimed that asking for
``FLG_TOPOCTR`` "is what keeps the diurnal parallax in". That claim was false in
the only way that matters. The event solves ``alt_true + SD + refraction = 0``
using the library's altitude; measuring ``alt_true + SD`` from the *same* source
cancels any error common to both, so the invariant held identically whether or
not the parallax was there. Demonstrated by mutation: stripping ``FLG_TOPOCTR``
inside ``calc_ut`` — i.e. deleting the diurnal parallax globally — moved sunrise
by 0.77 s and every assertion still passed.

Skyfield reads its own JPL kernel through its own frame pipeline, so the
position is genuinely a second opinion and the parallax becomes an independent
term. Under the same mutation the measured limb now moves to −33.734′, which is
0.141′ from target — seven times the tolerance below. The refraction and
semidiameter terms were already independent (they are the test's own arithmetic,
not the library's) and mutating either still fails, as before.

The cost is a coverage window: `de421.bsp` spans 1899-07-29 to 2053, and outside
it Skyfield RAISES rather than quietly drifting. That is the right failure — a
maintainer adding a 1750 or 2145 case gets "outside the oracle's coverage"
instead of a misleading message about refraction. The library itself reaches far
wider; this module deliberately does not claim to check those epochs.

Needs `skyfield` and `skyfield-data`, both hard install dependencies.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as L
from libephemeris.constants import (
    SUN,
    FLG_SWIEPH,
    CALC_RISE,
    CALC_SET,
    BIT_DISC_CENTER,
    BIT_NO_REFRACTION,
)

skyfield_api = pytest.importorskip("skyfield.api")
skyfield_data = pytest.importorskip("skyfield_data")

_AU_KM = 149597870.7

#: Deliberately the library's own solar radius, so what the limb assertion
#: checks is the semidiameter TERM (computed here from Skyfield's distance),
#: not the radius constant — which is common mode: a change of the library's
#: value by up to ~800 km (the IAU 2015 nominal 695700 included) moves the limb
#: by <= 0.007' and passes inside the 0.02' tolerance unseen.
_RSUN_KM = 696000.0

#: Standard atmosphere. Passed EXPLICITLY everywhere below, because the library
#: default is a literal 0 C (see `test_the_zero_default_is_zero_celsius`).
_STD_PRESSURE = 1013.25
_STD_TEMPERATURE = 15.0

#: Where the Sun's true upper limb sits at the event: minus the refraction at the
#: apparent horizon. DERIVED, not measured — it is exactly
#: ``_sinclair_refraction_deg(0.0, 1013.25, 15.0)`` = 33.59332', which is what
#: every event away from sea level returns to five decimal places.
_LIMB_ON_HORIZON_ARCMIN = -33.5933

#: The grid below is at ``geopos[2] = 0``, where the measured values spread over
#: 0.0069' (0.41") rather than landing on the constant exactly. That spread is
#: not ephemeris noise: at sea level ``calc_dip`` returns exactly 0, so
#: ``_rise_true_to_apparent``'s clamp puts a 0.56-degree step in the objective
#: function right at the root, the solver classifies it as the dip discontinuity
#: and returns a bisection midpoint instead of refining. Move the observer to
#: 50 m and the same grid spreads 0.0014" — 300x tighter. Sea level is kept
#: because it is the convention rise/set tables use, so the tolerance has to
#: cover the artifact: 0.02' is ~3x it, and still 7x below the 0.141' a dropped
#: parallax costs.
_LIMB_TOL_ARCMIN = 0.02

#: Same quantity with the library default atmosphere (0 C): colder air refracts
#: more, so the Sun has further to climb.
_LIMB_AT_ZERO_C_ARCMIN = -36.7395

#: Equator to just inside the polar circle, both hemispheres, all four seasons,
#: all inside the oracle's coverage window (see the module docstring).
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

_EVENTS = [("rise", CALC_RISE), ("set", CALC_SET)]


@pytest.fixture(scope="module")
def oracle():
    """Skyfield, loaded once: the second opinion everything is checked against."""
    loader = skyfield_api.Loader(skyfield_data.get_skyfield_data_path(), verbose=False)
    ephemeris = loader("de421.bsp")
    return {
        "ts": loader.timescale(),
        "earth": ephemeris["earth"],
        "sun": ephemeris["sun"],
        "wgs84": skyfield_api.wgs84,
    }


# NOTE ON GLOBAL STATE. Nothing in this module touches the global observer: the
# oracle is Skyfield, which takes the site as an argument. An earlier version
# added an autouse fixture that called ``L.set_topo(0, 0, 0)`` on teardown,
# claiming to "assert" the invariant. It asserted nothing, it duplicated
# `tests/conftest.py`'s `reset_ephemeris_state` (which saves and restores
# `_TOPO` at a wider scope, so its teardown runs later anyway), and it moved the
# state in the wrong direction: `_TOPO = None` makes `calc_ut(FLG_TOPOCTR)`
# raise, while `Topos(0, 0, 0)` silently returns a Null-Island position. Four
# other tests depend on that raise.


def _true_altitude_and_semidiameter(oracle, jd_ut: float, lon: float, lat: float):
    """Sun's true topocentric altitude and semidiameter, in arcminutes.

    Both come from Skyfield — its kernel, its frames, its topocentric geometry —
    so nothing on this side of the comparison is computed by the code under test.
    """
    time = oracle["ts"].ut1_jd(jd_ut)
    observer = oracle["earth"] + oracle["wgs84"].latlon(lat, lon, elevation_m=0)
    apparent = observer.at(time).observe(oracle["sun"]).apparent()
    altitude, _azimuth, distance = apparent.altaz()
    semidiameter = math.degrees(math.asin(_RSUN_KM / (distance.au * _AU_KM)))
    return altitude.degrees * 60.0, semidiameter * 60.0


def _altitude_rate_arcsec_per_second(oracle, jd_ut: float, lon: float, lat: float) -> float:
    """Signed vertical rate at `jd_ut`. Positive rising, negative setting."""
    before, _ = _true_altitude_and_semidiameter(oracle, jd_ut - 30.0 / 86400.0, lon, lat)
    after, _ = _true_altitude_and_semidiameter(oracle, jd_ut + 30.0 / 86400.0, lon, lat)
    return (after - before) * 60.0 / 60.0


def _event(jd0, lon, lat, rsmi, atpress=_STD_PRESSURE, attemp=_STD_TEMPERATURE, alt=0.0):
    retflag, tret = L.rise_trans(jd0, SUN, rsmi, (lon, lat, alt), atpress, attemp, FLG_SWIEPH)
    assert retflag == 0, f"expected an event, got retflag={retflag}"
    return tret[0]


def _ids(value):
    return value if isinstance(value, str) else None


@pytest.mark.parametrize("name,lon,lat,y,m,d", _GRID, ids=_ids)
@pytest.mark.parametrize("label,rsmi", _EVENTS, ids=_ids)
def test_the_upper_limb_sits_on_the_refracted_horizon(oracle, name, lon, lat, y, m, d, label, rsmi):
    """The one invariant: independent of latitude, season and hemisphere."""
    jd0 = L.julday(y, m, d, 0.0, L.GREG_CAL)
    jd = _event(jd0, lon, lat, rsmi)
    altitude, semidiameter = _true_altitude_and_semidiameter(oracle, jd, lon, lat)
    limb = altitude + semidiameter
    assert abs(limb - _LIMB_ON_HORIZON_ARCMIN) < _LIMB_TOL_ARCMIN, (
        f"{name} {label}: upper limb at {limb:.5f}', expected "
        f"{_LIMB_ON_HORIZON_ARCMIN:.4f}' +/- {_LIMB_TOL_ARCMIN}'. The horizon "
        f"convention moved — check the refraction model, the semidiameter term "
        f"and the topocentric parallax, in that order of likelihood."
    )

    # The altitude invariant is SYMMETRIC in rise and set, so on its own the two
    # parametrised directions are one assertion run twice. A mutant that returns
    # the rising crossing for CALC_SET passed all 51 tests. The sign of the
    # vertical rate is the one thing that tells them apart, so it is asserted.
    rate = _altitude_rate_arcsec_per_second(oracle, jd, lon, lat)
    expected_sign = 1.0 if rsmi == CALC_RISE else -1.0
    assert rate * expected_sign > 0.0, (
        f"{name} {label}: the Sun's altitude is moving {rate:+.2f}\"/s here, which is the "
        f"wrong direction for a {label} — the event type is not selecting what it claims"
    )


@pytest.mark.parametrize("name,lon,lat,y,m,d", _GRID, ids=_ids)
@pytest.mark.parametrize("label,rsmi", _EVENTS, ids=_ids)
def test_disc_center_moves_the_target_from_the_limb_to_the_centre(oracle, name, lon, lat, y, m, d, label, rsmi):
    """BIT_DISC_CENTER drops the semidiameter and nothing else.

    Both directions, because a sign error that only affects setting would look
    correct in every rise-only assertion above.
    """
    jd0 = L.julday(y, m, d, 0.0, L.GREG_CAL)
    jd = _event(jd0, lon, lat, rsmi | BIT_DISC_CENTER)
    centre, _semidiameter = _true_altitude_and_semidiameter(oracle, jd, lon, lat)
    assert abs(centre - _LIMB_ON_HORIZON_ARCMIN) < _LIMB_TOL_ARCMIN, (
        f"{name} {label}: with BIT_DISC_CENTER the CENTRE should sit where the limb "
        f"otherwise does ({_LIMB_ON_HORIZON_ARCMIN:.4f}'), got {centre:.5f}'"
    )


@pytest.mark.parametrize("name,lon,lat,y,m,d", _GRID, ids=_ids)
@pytest.mark.parametrize("label,rsmi", _EVENTS, ids=_ids)
def test_the_geometric_horizon_is_exactly_zero(oracle, name, lon, lat, y, m, d, label, rsmi):
    """Disc centre with refraction off must land on 0.0000', not merely near it.

    With both physical terms removed the only things left are the root-finder and
    the two position pipelines, so this is the tightest statement the module can
    make. 2e-4' is ~10x the worst observed residual (1.9e-5') and ~1 ms of time;
    it is deliberately not looser, because loosening it is how a search
    regression hides. Disabling the solver's Newton refinement fails 15 of the 16
    cases at this tolerance.
    """
    jd0 = L.julday(y, m, d, 0.0, L.GREG_CAL)
    jd = _event(jd0, lon, lat, rsmi | BIT_DISC_CENTER | BIT_NO_REFRACTION)
    centre, _semidiameter = _true_altitude_and_semidiameter(oracle, jd, lon, lat)
    assert abs(centre) < 2e-4, f"{name} {label}: geometric horizon crossing at {centre:.6f}', expected 0"


def test_the_zero_default_is_zero_celsius_not_the_standard_atmosphere(oracle):
    """attemp=0 means 0 C, as in the reference — not an implicit 15 C.

    Worth its own test because it is the cheapest mistake to make in a caller and
    the hardest to notice: the event still looks perfectly reasonable, it is just
    17 seconds off a standard-atmosphere convention.
    """
    jd0 = L.julday(2024, 3, 20, 0.0, L.GREG_CAL)
    lon, lat = 12.5, 41.9

    default_jd = _event(jd0, lon, lat, CALC_RISE, atpress=0.0, attemp=0.0)
    standard_jd = _event(jd0, lon, lat, CALC_RISE)

    altitude, semidiameter = _true_altitude_and_semidiameter(oracle, default_jd, lon, lat)
    limb = altitude + semidiameter
    assert abs(limb - _LIMB_AT_ZERO_C_ARCMIN) < _LIMB_TOL_ARCMIN, (
        f"default-atmosphere limb at {limb:.5f}', expected {_LIMB_AT_ZERO_C_ARCMIN:.4f}' (0 C)"
    )

    seconds = (standard_jd - default_jd) * 86400.0
    assert 15.0 < seconds < 19.0, (
        f"0 C -> 15 C should move sunrise by about 17 s, got {seconds:.1f} s"
    )


def test_observer_elevation_does_not_lower_the_horizon():
    """geopos[2] never depresses the horizon, whatever it appears to do.

    Three statements, because the naive one is wrong in an instructive way:

    1. With an EXPLICIT pressure, elevation changes the event by nothing that
       matters (13 ms, and not for the reason one would guess — see the note in
       `rise_trans`).
    2. With the default atpress=0 the elevation is what the barometric estimate
       reads, so a higher observer gets thinner air, less refraction and a LATER
       sunrise — the opposite of the physical intuition, and the reason this
       looks like a bug when it is a convention.
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

    # (3) the dip is the only thing that actually lowers the horizon. Measured
    # under the SAME standard atmosphere as `fixed_high`, because the effect is
    # atmosphere-dependent: the default 0 C gives -362.3 s here, not -360.0 s.
    retflag, tret = L.rise_trans_true_hor(
        jd0, SUN, CALC_RISE, (lon, lat, 1000.0), _STD_PRESSURE, _STD_TEMPERATURE, -100.0
    )
    assert retflag == 0
    dip_seconds = (tret[0] - fixed_high) * 86400.0
    assert -380.0 < dip_seconds < -340.0, (
        f"horhgt=-100 at 1000 m should pull sunrise ~360 s EARLIER, got {dip_seconds:+.1f} s"
    )


def test_a_shallow_horhgt_at_sea_level_does_nothing():
    """The dead band beside horhgt=0, pinned because the docs now warn about it.

    `_rise_true_to_apparent` clamps to the unrefracted altitude whenever the
    apparent altitude falls below the dip, and at sea level the dip is exactly 0.
    So the whole 0.56 degrees of refraction between horhgt=0 and horhgt=-0.56 is
    a range in which asking for a lower horizon changes nothing, and past its
    edge the response is a KINK rather than a step — 322 s per degree, so -0.57
    moves the event 3.3 s and -0.60 moves it 12.9 s. A caller who asks for a
    horizon depressed by 30 arcminutes gets the undepressed sunrise.

    Pre-existing behaviour and not obviously wrong to fix — but it is exactly
    the trap `horhgt` documentation should stop, so it is written down as a fact
    rather than left to be rediscovered.
    """
    jd0 = L.julday(2024, 3, 20, 0.0, L.GREG_CAL)
    geopos = (12.5, 41.9, 0.0)

    def rise(horhgt):
        retflag, tret = L.rise_trans_true_hor(
            jd0, SUN, CALC_RISE, geopos, _STD_PRESSURE, _STD_TEMPERATURE, horhgt
        )
        assert retflag == 0
        return tret[0]

    flat = rise(0.0)
    for inside_band in (-0.10, -0.30, -0.50):
        assert abs(rise(inside_band) - flat) * 86400.0 < 1.0, (
            f"horhgt={inside_band} moved the event; the dead band this test documents is gone, "
            f"which is an improvement — update the docstring in eclipse.py to match"
        )

    beyond = (rise(-0.60) - flat) * 86400.0
    assert beyond < -5.0, (
        f"past the refraction the horizon should finally move; got {beyond:+.1f} s"
    )
