"""
Tests for the SIDBIT projection-flag guard (this version).

The SIDBIT_* flags (>= 256: ECL_T0, SSY_PLANE, USER_UT, ECL_DATE,
NO_PREC_OFFSET, PREC_ORIG) select alternative ecliptic/equinox projections
that are not implemented in this version. set_sid_mode() must strip them,
warn, and keep the base ayanamsha mode — instead of silently falling back
to Lahiri for the composite value. Planned for a future release (NEXT.md).
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
def test_sidbit_stripped_and_warned():
    """A composite mode warns and stores only the base mode."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ephem.set_sid_mode(ephem.SIDM_FAGAN_BRADLEY | ephem.SIDBIT_ECL_T0)
    assert state.get_sid_mode() == ephem.SIDM_FAGAN_BRADLEY
    assert any(issubclass(w.category, UserWarning) for w in caught)


@pytest.mark.unit
def test_composite_uses_base_mode_not_lahiri():
    """The composite value must compute the base ayanamsha, not Lahiri."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ephem.set_sid_mode(ephem.SIDM_FAGAN_BRADLEY | ephem.SIDBIT_SSY_PLANE)
    a_composite = ephem.get_ayanamsa_ut(JD)

    ephem.set_sid_mode(ephem.SIDM_FAGAN_BRADLEY)
    a_base = ephem.get_ayanamsa_ut(JD)

    ephem.set_sid_mode(ephem.SIDM_LAHIRI)
    a_lahiri = ephem.get_ayanamsa_ut(JD)

    assert abs(a_composite - a_base) < 1e-9
    assert abs(a_composite - a_lahiri) > 0.5  # clearly not the old Lahiri fallback


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
