"""Sealed-mode errors for a custom LEB file that omits a core body.

In calculation mode ``leb`` with a partial file installed via
``set_leb_file()``, a core body absent from the active reader must raise the
declared :class:`UnknownBodyError` from every public entry point — the
module-level ``calc_ut``/``calc``, the ``EphemerisContext`` entry points, and
the vector consumers — never an undeclared ``KeyError``.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import MARS, MOON, NODBIT_OSCU, SUN
from libephemeris.exceptions import UnknownBodyError
from libephemeris.time_utils import julday


@pytest.fixture
def sealed_partial_leb(test_leb_file_minimal):
    """Force sealed LEB mode onto the Sun+Earth-only fixture file."""
    from libephemeris import state

    prev_mode = state.get_calc_mode()
    ephem.set_calc_mode("leb")
    ephem.set_leb_file(test_leb_file_minimal)

    yield

    ephem.set_calc_mode(prev_mode)
    ephem.set_leb_file("")


@pytest.fixture
def jd_mid(test_leb_file_minimal):
    """A JD inside the minimal fixture's stored 2024 coverage."""
    return julday(2024, 6, 15, 12.0)


class TestSealedMissingBodyErrors:
    """Missing core body in mode=leb -> declared error, not KeyError."""

    def test_module_calc_ut_raises_unknown_body(self, sealed_partial_leb, jd_mid):
        with pytest.raises(UnknownBodyError) as exc_info:
            ephem.calc_ut(jd_mid, MARS, 0)
        assert exc_info.value.body_id == MARS
        assert not isinstance(exc_info.value, KeyError)

    def test_module_calc_raises_unknown_body(self, sealed_partial_leb, jd_mid):
        with pytest.raises(UnknownBodyError) as exc_info:
            ephem.calc(jd_mid, MOON, 0)
        assert exc_info.value.body_id == MOON
        assert not isinstance(exc_info.value, KeyError)

    def test_context_calc_raises_unknown_body(self, sealed_partial_leb, jd_mid):
        from libephemeris import EphemerisContext

        ctx = EphemerisContext()
        with pytest.raises(UnknownBodyError) as exc_info:
            ctx.calc(jd_mid, MARS, 0)
        assert exc_info.value.body_id == MARS
        assert not isinstance(exc_info.value, KeyError)

    def test_context_calc_ut_with_context_leb_raises_unknown_body(
        self, sealed_partial_leb, jd_mid, test_leb_file_minimal
    ):
        from libephemeris import EphemerisContext

        ctx = EphemerisContext()
        ctx.set_leb_file(test_leb_file_minimal)
        try:
            with pytest.raises(UnknownBodyError) as exc_info:
                ctx.calc_ut(jd_mid, MARS, 0)
        finally:
            ctx.set_leb_file(None)
        assert exc_info.value.body_id == MARS
        assert not isinstance(exc_info.value, KeyError)

    def test_nod_aps_vector_path_raises_unknown_body(self, sealed_partial_leb, jd_mid):
        with pytest.raises(UnknownBodyError) as exc_info:
            ephem.nod_aps(jd_mid, MARS, NODBIT_OSCU, 0)
        assert exc_info.value.body_id == MARS
        assert not isinstance(exc_info.value, KeyError)

    def test_stored_body_still_served(self, sealed_partial_leb, jd_mid):
        """Bodies present in the partial file keep working in sealed mode."""
        position, _flags = ephem.calc_ut(jd_mid, SUN, 0)
        assert 0.0 <= position[0] < 360.0
        assert position[2] > 0.9

    def test_context_range_miss_without_global_raises_range_error(
        self, test_leb_file_minimal
    ):
        """A context-file date-range miss reports the context file's bounds.

        With no global LEB configured, classifying the miss against the
        global configuration would misreport a plain range miss on a stored
        body as a missing-body condition.
        """
        from unittest.mock import patch

        from libephemeris import EphemerisContext, state
        from libephemeris.exceptions import EphemerisRangeError

        prev_mode = state.get_calc_mode()
        ephem.set_calc_mode("leb")
        ephem.set_leb_file(None)
        ctx = EphemerisContext()
        ctx.set_leb_file(test_leb_file_minimal)
        jd_out = julday(1900, 1, 1, 0.0)
        try:
            with (
                patch("libephemeris.state._discover_leb_file", return_value=None),
                patch(
                    "libephemeris.state._discover_reviewed_leb_tier_cores",
                    return_value={},
                ),
            ):
                with pytest.raises(EphemerisRangeError) as exc_info:
                    ctx.calc_ut(jd_out, SUN, 0)
            err = exc_info.value
            assert err.start_jd is not None and err.end_jd is not None
            assert not (err.start_jd <= jd_out <= err.end_jd)
            assert err.body_id == SUN
        finally:
            ctx.set_leb_file(None)
            ephem.set_calc_mode(prev_mode)
            ephem.set_leb_file("")
