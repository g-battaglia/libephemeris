"""
Tests for the SIDBIT projection-flag handling.

The implemented projection flags SIDBIT_ECL_T0, SIDBIT_SSY_PLANE and
SIDBIT_USER_UT are applied by the module-level setter (they are retained and
consumed by the sidereal path). The remaining defined projection flags
(SIDBIT_ECL_DATE, SIDBIT_NO_PREC_OFFSET, SIDBIT_PREC_ORIG) are accepted
silently for compatibility -- no warning, matching the reference -- but
currently reduce to the base ayanamsha value. Only a genuinely unknown high bit
still strips-and-warns. EphemerisContext strips every projection flag (it does
not carry them).
"""

from __future__ import annotations

import warnings

import pytest

import libephemeris as ephem
from libephemeris import state


@pytest.fixture(autouse=True)
def _reset_state():
    yield
    ephem.set_sid_mode(ephem.SIDM_LAHIRI)
    ephem.close()


JD = 2451545.0


@pytest.mark.unit
def test_unknown_high_bit_stripped_and_warned():
    """A genuinely unknown high bit (not a defined SIDBIT) warns and stores only
    the base mode."""
    unknown_bit = 0x4000  # 16384: above every defined SIDBIT flag
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ephem.set_sid_mode(ephem.SIDM_FAGAN_BRADLEY | unknown_bit)
    assert state.get_sid_mode() == ephem.SIDM_FAGAN_BRADLEY
    assert any(issubclass(w.category, UserWarning) for w in caught)


@pytest.mark.unit
@pytest.mark.parametrize(
    "bit",
    [
        ephem.SIDBIT_NO_PREC_OFFSET,
        ephem.SIDBIT_PREC_ORIG,
    ],
)
def test_accepted_sidbit_no_warn_and_reduces_to_base(bit):
    """SIDBIT_NO_PREC_OFFSET / PREC_ORIG are accepted silently (no warning,
    matching the reference) and currently reduce to the base ayanamsha value,
    consistently regardless of the active backend."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ephem.set_sid_mode(ephem.SIDM_FAGAN_BRADLEY | bit)
    assert state.get_sid_mode() == ephem.SIDM_FAGAN_BRADLEY
    assert not caught
    with_flag = ephem.get_ayanamsa_ut(JD)

    ephem.set_sid_mode(ephem.SIDM_FAGAN_BRADLEY)
    base = ephem.get_ayanamsa_ut(JD)
    assert with_flag == pytest.approx(base, abs=1e-12)


@pytest.mark.unit
def test_ecl_date_applies_the_geometric_delta():
    """SIDBIT_ECL_DATE shifts the ayanamsha by _ecl_date_ayanamsha_delta.

    The shift is applied for defining-pair modes and is identical on the
    LEB fast path (which delegates) and the direct path.
    """
    from libephemeris.ayanamsha_definitions import AYANAMSHA_DEFINING
    from libephemeris.planets import _ecl_date_ayanamsha_delta
    from libephemeris.state import get_timescale

    mode = 21  # Suryasiddhanta: defining epoch ~499 CE, delta well resolved
    ephem.set_sid_mode(mode)
    base = ephem.get_ayanamsa_ut(JD)
    ephem.set_sid_mode(mode | ephem.SIDBIT_ECL_DATE)
    with_flag = ephem.get_ayanamsa_ut(JD)

    dv, t0 = AYANAMSHA_DEFINING[mode]
    jd_tt = float(get_timescale().ut1_jd(JD).tt)
    expected = _ecl_date_ayanamsha_delta(dv, t0, jd_tt)
    measured = (with_flag - base + 180.0) % 360.0 - 180.0
    assert measured == pytest.approx(expected, abs=1e-12)
    assert abs(expected) > 1e-5  # the delta is genuinely non-zero here


@pytest.mark.unit
@pytest.mark.parametrize(
    "bit",
    [ephem.SIDBIT_ECL_T0, ephem.SIDBIT_SSY_PLANE, ephem.SIDBIT_USER_UT],
)
def test_implemented_sidbit_no_warn_and_retained(bit):
    """Implemented projection flags do not warn; the base mode is stored and
    the flag is retained for the sidereal path to apply."""
    from libephemeris.planets import _get_sidereal_bits

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ephem.set_sid_mode(ephem.SIDM_FAGAN_BRADLEY | bit)
    assert state.get_sid_mode() == ephem.SIDM_FAGAN_BRADLEY
    assert not caught
    assert _get_sidereal_bits() & bit


@pytest.mark.unit
def test_ecl_date_delta_geometry():
    """The staged SIDBIT_ECL_DATE realization is geometrically self-consistent.

    ``_ecl_date_ayanamsha_delta`` returns the amount to add to the Method-B
    ayanamsha so the sidereal zero is referred to the mean ecliptic of date. It
    must vanish at the defining epoch (the two ecliptics coincide there), grow
    with the epoch separation, and the corrected value must equal the tropical
    longitude of the sidereal zero measured on the mean ecliptic of date. This
    checks the geometry directly (no reference ephemeris involved).
    """
    import math

    from libephemeris.planets import _ecl_date_ayanamsha_delta
    from libephemeris.precession_vondrak import (
        _ltp_ecliptic_frame,
        method_b_accumulated_precession,
    )

    dv, t0 = 24.0, 2451545.0  # arbitrary defining value (deg) and epoch (JD TT)

    # (1) Vanishes at the defining epoch.
    assert abs(_ecl_date_ayanamsha_delta(dv, t0, t0)) < 1e-9

    # (2) Grows with epoch separation.
    d_near = abs(_ecl_date_ayanamsha_delta(dv, t0, t0 + 50 * 365.25))
    d_far = abs(_ecl_date_ayanamsha_delta(dv, t0, t0 + 500 * 365.25))
    assert 0.0 < d_near < d_far

    # (3) Corrected value == tropical longitude of the sidereal zero on the
    #     mean ecliptic of date.
    tt = t0 + 150 * 365.25
    method_b = dv + method_b_accumulated_precession(tt, t0)
    corrected = (method_b + _ecl_date_ayanamsha_delta(dv, t0, tt)) % 360.0
    e0 = _ltp_ecliptic_frame(t0)
    lon0 = math.radians(dv)
    cx, sx = math.cos(lon0), math.sin(lon0)
    z = (
        e0[0][0] * cx + e0[1][0] * sx,
        e0[0][1] * cx + e0[1][1] * sx,
        e0[0][2] * cx + e0[1][2] * sx,
    )
    ed = _ltp_ecliptic_frame(tt)
    wx = ed[0][0] * z[0] + ed[0][1] * z[1] + ed[0][2] * z[2]
    wy = ed[1][0] * z[0] + ed[1][1] * z[1] + ed[1][2] * z[2]
    lon_date = math.degrees(math.atan2(wy, wx)) % 360.0
    assert corrected == pytest.approx(lon_date, abs=1e-9)


@pytest.mark.unit
def test_composite_uses_base_mode():
    """The composite value must compute its independently defined base mode."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ephem.set_sid_mode(ephem.SIDM_J1900 | ephem.SIDBIT_SSY_PLANE)
    a_composite = ephem.get_ayanamsa_ut(JD)

    ephem.set_sid_mode(ephem.SIDM_J1900)
    a_base = ephem.get_ayanamsa_ut(JD)

    ephem.set_sid_mode(ephem.SIDM_J2000)
    a_other = ephem.get_ayanamsa_ut(JD)

    assert abs(a_composite - a_base) < 1e-9
    assert abs(a_composite - a_other) > 0.5


@pytest.mark.unit
def test_plain_mode_unaffected():
    """A plain mode (no SIDBIT) is stored and warned-free."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ephem.set_sid_mode(ephem.SIDM_LAHIRI)
    assert state.get_sid_mode() == ephem.SIDM_LAHIRI
    assert not caught


@pytest.mark.unit
def test_user_mode_preserved():
    """SIDM_USER (255) must survive masking untouched."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ephem.set_sid_mode(255, JD, 25.0)
    assert state.get_sid_mode() == 255
    assert not caught


# =============================================================================
# EphemerisContext parity with the module setter
# =============================================================================
#
# EphemerisContext.set_sid_mode() must apply the identical SIDBIT strip-and-warn
# guard as state.set_sid_mode(); otherwise ctx.get_sid_mode() would report the
# raw composite (e.g. 257 instead of 1) even though the two stores must stay in
# lockstep (known-differences 10.3).


@pytest.mark.unit
def test_context_strips_sidbit_and_warns():
    """Context set_sid_mode strips SIDBIT flags, warns, and stores the base mode."""
    from libephemeris.context import EphemerisContext

    ctx = EphemerisContext()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ctx.set_sid_mode(ephem.SIDM_FAGAN_BRADLEY | ephem.SIDBIT_ECL_T0)
    assert ctx.get_sid_mode() == ephem.SIDM_FAGAN_BRADLEY
    assert any(issubclass(w.category, UserWarning) for w in caught)


@pytest.mark.unit
def test_context_matches_module_stripped_mode():
    """ctx.get_sid_mode() must equal state.get_sid_mode() for a composite value."""
    from libephemeris.context import EphemerisContext

    composite = ephem.SIDM_LAHIRI | ephem.SIDBIT_ECL_T0  # 1 | 256 = 257

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ephem.set_sid_mode(composite)
        ctx = EphemerisContext()
        ctx.set_sid_mode(composite)

    assert ctx.get_sid_mode() == state.get_sid_mode() == ephem.SIDM_LAHIRI
    # Full tuple parity too (t0/ayan_t0 carried through unchanged).
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ephem.set_sid_mode(composite, JD, 23.5)
        ctx2 = EphemerisContext()
        ctx2.set_sid_mode(composite, JD, 23.5)
    assert ctx2.get_sid_mode(full=True) == state.get_sid_mode(full=True)


@pytest.mark.unit
def test_context_plain_mode_unaffected():
    """A plain mode is stored verbatim in the context and emits no warning."""
    from libephemeris.context import EphemerisContext

    ctx = EphemerisContext()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ctx.set_sid_mode(ephem.SIDM_LAHIRI)
    assert ctx.get_sid_mode() == ephem.SIDM_LAHIRI
    assert not caught


@pytest.mark.unit
def test_context_user_mode_preserved():
    """SIDM_USER (255) survives masking untouched in the context."""
    from libephemeris.context import EphemerisContext

    ctx = EphemerisContext()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ctx.set_sid_mode(255, JD, 25.0)
    assert ctx.get_sid_mode() == 255
    assert not caught
