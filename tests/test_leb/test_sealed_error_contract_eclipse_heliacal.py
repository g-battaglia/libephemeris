"""Sealed-mode error contract for eclipse/heliacal LEB dispatch.

In sealed ``leb`` mode a probe outside the active LEB core must surface as
the documented :class:`EphemerisRangeError` (the type ``calc_ut``/``calc``
already raise), never as the reader's internal ValueError/KeyError. In
``auto`` mode the directly-callable pythonic lunar functions (re-exported
from ``libephemeris``) must keep the pre-existing behaviour: an LEB range
miss retries transparently on the Skyfield path and returns a result.
"""

from __future__ import annotations

import pytest

import libephemeris
from libephemeris.constants import CALC_RISE, SUN
from libephemeris.exceptions import EphemerisRangeError, UnknownBodyError
from libephemeris.time_utils import julday

GEOPOS_ROME = (12.5, 41.9, 0.0)
ATMO = (1013.25, 15.0, 50.0, 0.25)
OBSERVER = (25.0, 1.0, 1.0, 1.0, 0.0, 0.0)

# Outside the 2023-2028 test LEB core, inside the DE440s kernel (1849-2150).
JD_OUT = julday(2018, 7, 1, 0.0)
# Maximum of the 2018-07-27 total lunar eclipse (visible from Rome).
JD_LUN_2018 = julday(2018, 7, 27, 20.37)


@pytest.fixture
def sealed_leb(test_leb_file):
    """Sealed leb mode backed by the range-limited (2023-2028) test core."""
    from libephemeris import state

    prev_mode = state.get_calc_mode()
    libephemeris.set_calc_mode("leb")
    libephemeris.set_leb_file(test_leb_file)
    yield
    libephemeris.set_calc_mode(prev_mode)
    libephemeris.set_leb_file("")


@pytest.fixture
def auto_leb(test_leb_file):
    """Auto mode with the range-limited test core as the active LEB file."""
    from libephemeris import state

    prev_mode = state.get_calc_mode()
    libephemeris.set_calc_mode("auto")
    libephemeris.set_leb_file(test_leb_file)
    yield
    libephemeris.set_calc_mode(prev_mode)
    libephemeris.set_leb_file("")


class TestSealedEclipseErrorContract:
    def test_rise_trans_out_of_range_raises_typed_error(self, sealed_leb):
        from libephemeris.eclipse import rise_trans

        with pytest.raises(EphemerisRangeError) as exc_info:
            rise_trans(JD_OUT, SUN, CALC_RISE, GEOPOS_ROME)
        err = exc_info.value
        assert err.start_jd is not None and err.end_jd is not None
        assert err.start_jd < err.end_jd
        assert not (err.start_jd <= JD_OUT <= err.end_jd)

    def test_lun_eclipse_how_out_of_range_raises_typed_error(self, sealed_leb):
        from libephemeris.eclipse import lun_eclipse_how

        with pytest.raises(EphemerisRangeError):
            lun_eclipse_how(JD_LUN_2018, list(GEOPOS_ROME))

    def test_lun_eclipse_how_pythonic_direct_raises_typed_error(self, sealed_leb):
        from libephemeris.eclipse import _lun_eclipse_how_pythonic

        with pytest.raises(EphemerisRangeError):
            _lun_eclipse_how_pythonic(JD_LUN_2018, 41.9, 12.5)


class TestSealedHeliacalErrorContract:
    def test_vis_limit_mag_out_of_range_raises_typed_error(self, sealed_leb):
        from libephemeris.heliacal import vis_limit_mag

        with pytest.raises(EphemerisRangeError):
            vis_limit_mag(JD_OUT, GEOPOS_ROME, ATMO, OBSERVER, "Venus", 0)


class TestSealedMissTranslationHelper:
    """Unit contract of the shared eclipse/heliacal translation helper."""

    def _sealed(self):
        from libephemeris import state

        prev_mode = state.get_calc_mode()
        libephemeris.set_calc_mode("leb")
        return prev_mode

    def test_range_miss_preserves_bounds(self):
        from libephemeris.eclipse import _raise_if_sealed_leb_miss

        cause = ValueError(
            "JD 2458300.5 outside range [2459945.5, 2461771.5] for body 0"
        )
        prev_mode = self._sealed()
        try:
            with pytest.raises(EphemerisRangeError) as exc_info:
                _raise_if_sealed_leb_miss(cause)
        finally:
            libephemeris.set_calc_mode(prev_mode)
        err = exc_info.value
        assert err.requested_jd == 2458300.5
        assert err.start_jd == 2459945.5
        assert err.end_jd == 2461771.5
        assert err.body_id == 0
        assert err.__cause__ is cause

    def test_body_map_miss_raises_unknown_body_error(self):
        # Same type planets._raise_leb_range_miss uses for an unstored
        # core body: absent body -> UnknownBodyError, not a range error.
        from libephemeris.eclipse import _raise_if_sealed_leb_miss

        prev_mode = self._sealed()
        try:
            with pytest.raises(UnknownBodyError) as exc_info:
                _raise_if_sealed_leb_miss(KeyError("Body 3 not in LEB file"))
        finally:
            libephemeris.set_calc_mode(prev_mode)
        err = exc_info.value
        assert err.body_id == 3
        assert "not stored" in err.message

    def test_unrelated_error_reraised_verbatim(self):
        from libephemeris.eclipse import _raise_if_sealed_leb_miss

        cause = ValueError("geopos must have at least 3 elements")
        prev_mode = self._sealed()
        try:
            with pytest.raises(ValueError) as exc_info:
                _raise_if_sealed_leb_miss(cause)
        finally:
            libephemeris.set_calc_mode(prev_mode)
        assert exc_info.value is cause

    def test_auto_mode_returns_without_raising(self):
        from libephemeris import state
        from libephemeris.eclipse import _raise_if_sealed_leb_miss

        prev_mode = state.get_calc_mode()
        libephemeris.set_calc_mode("auto")
        try:
            _raise_if_sealed_leb_miss(
                ValueError("JD 1.0 outside range [2.0, 3.0] for body 0")
            )
        finally:
            libephemeris.set_calc_mode(prev_mode)


class TestDirectPythonicAutoFallback:
    """Direct calls to the re-exported lunar pythonic API keep the
    pre-branch auto-mode semantics: an out-of-core date computes via the
    Skyfield retry instead of surfacing the reader's ValueError."""

    def test_lun_eclipse_how_pythonic_direct_out_of_core(self, auto_leb):
        from libephemeris.eclipse import _lun_eclipse_how_pythonic

        rc, attr = _lun_eclipse_how_pythonic(JD_LUN_2018, 41.9, 12.5)
        assert rc != 0
        assert isinstance(attr[0], float)
        # Total lunar eclipse near maximum: umbral magnitude well above 1.
        assert attr[0] > 1.0

    def test_lun_eclipse_when_loc_pythonic_direct_out_of_core(self, auto_leb):
        from libephemeris.eclipse import _lun_eclipse_when_loc_pythonic

        retflag, tret, attr = _lun_eclipse_when_loc_pythonic(
            julday(2018, 7, 20, 0.0), 41.9, 12.5
        )
        assert retflag > 0
        # Finds the 2018-07-27 eclipse (maximum ~20:22 UT).
        assert abs(tret[0] - JD_LUN_2018) < 0.1


class TestSealedBodyMissTypedErrors:
    """Sealed body-map miss -> UnknownBodyError, never RuntimeError or zeros."""

    @pytest.fixture
    def sealed_partial_leb(self, test_leb_file_minimal):
        """Sealed leb mode backed by the Sun+Earth-only fixture file."""
        from libephemeris import state

        prev_mode = state.get_calc_mode()
        libephemeris.set_calc_mode("leb")
        libephemeris.set_leb_file(test_leb_file_minimal)
        yield
        libephemeris.set_calc_mode(prev_mode)
        libephemeris.set_leb_file("")

    def test_rise_trans_missing_body_raises_unknown_body(self, sealed_partial_leb):
        from libephemeris.constants import MOON
        from libephemeris.eclipse import rise_trans

        with pytest.raises(UnknownBodyError):
            rise_trans(julday(2024, 6, 1, 0.0), MOON, CALC_RISE, GEOPOS_ROME)

    def test_lun_eclipse_how_pythonic_missing_moon_raises_unknown_body(
        self, sealed_partial_leb
    ):
        """No fabricated "no eclipse" zeros for a body the file cannot serve."""
        from libephemeris.eclipse import _lun_eclipse_how_pythonic

        with pytest.raises(UnknownBodyError):
            _lun_eclipse_how_pythonic(julday(2024, 6, 1, 0.0), 41.9, 12.5)

    def test_body_miss_regex_matches_tiered_reader_message(self):
        from libephemeris.eclipse import _LEB_BODY_MISS_RE

        for message in (
            "Body 4 not in LEB file",
            "Body 4 not in any LEB file",
            "Body 4 not in any installed LEB tier",
        ):
            match = _LEB_BODY_MISS_RE.search(message)
            assert match is not None, message
            assert match.group("body") == "4"
