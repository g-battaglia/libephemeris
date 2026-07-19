"""
Regression tests for the sealed-``leb`` typed-error contract at the two
eclipse coverage probes that previously leaked a raw ``RuntimeError``.

Both ``_calc_gamma`` and ``_calculate_local_eclipse_phases_impl`` probe the
active LEB reader and, on a coverage/body miss, used to swallow the
``KeyError``/``ValueError`` and fall through to ``state.get_planets()``. In
sealed ``leb`` mode ``get_planets()`` raises

    RuntimeError("JPL/SPICE ephemeris access is disabled in calculation
                  mode 'leb'")

so a coverage-edge eclipse aborted with an untyped ``RuntimeError`` instead of
the documented ``EphemerisRangeError`` / ``UnknownBodyError``. The sibling site
in ``_rise_trans_true_hor_impl`` was already guarded with
``_raise_if_sealed_leb_miss`` / ``_reraise_if_leb_range_error``; these tests pin
the same contract for the two eclipse probes.

The Skyfield/JPL branch is stubbed with a marker exception, so the tests need
no ephemeris data and never open a kernel: reaching the marker *is* the
assertion that the auto-mode fallback still runs.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from libephemeris.eclipse import (
    _calc_gamma,
    _calculate_local_eclipse_phases,
)
from libephemeris.exceptions import EphemerisRangeError, UnknownBodyError


# The exact message sealed mode produces when the kernel path is reached.
_SEALED_RUNTIME_ERROR = RuntimeError(
    "JPL/SPICE ephemeris access is disabled in calculation mode 'leb'"
)

_RANGE_MISS = "JD 2460310.5000 outside range [2396758.5, 2469807.5] for body 10"
_BODY_MISS = "Body 15 not in any LEB file"

_JD_2024 = 2460310.5
_GEOPOS = (12.5, 41.9, 0.0)


class _ReachedSkyfield(Exception):
    """Marker: the non-LEB (Skyfield/JPL) fallback branch was entered."""


def _sealed_get_planets(*args, **kwargs):
    """Stand in for state.get_planets() under a sealed leb runtime."""
    raise _SEALED_RUNTIME_ERROR


def _marker_get_planets(*args, **kwargs):
    """Stand in for state.get_planets() when the fallback is legitimate."""
    raise _ReachedSkyfield()


class TestCalcGammaSealedContract:
    """_calc_gamma must not degrade a sealed-mode miss into RuntimeError."""

    def test_range_miss_raises_typed_error_not_runtime_error(self):
        reader = MagicMock(name="LEBReader")

        with (
            patch("libephemeris.eclipse.get_calc_mode", return_value="leb"),
            patch("libephemeris.state.get_leb_reader", return_value=reader),
            patch(
                "libephemeris.fast_calc._apparent_icrs_cartesian",
                side_effect=ValueError(_RANGE_MISS),
            ),
            patch("libephemeris.state.get_planets", side_effect=_sealed_get_planets),
        ):
            with pytest.raises(EphemerisRangeError) as excinfo:
                _calc_gamma(_JD_2024)

        assert excinfo.value.start_jd == 2396758.5
        assert excinfo.value.end_jd == 2469807.5
        assert excinfo.value.body_id == 10

    def test_body_miss_raises_unknown_body_not_runtime_error(self):
        reader = MagicMock(name="LEBReader")

        with (
            patch("libephemeris.eclipse.get_calc_mode", return_value="leb"),
            patch("libephemeris.state.get_leb_reader", return_value=reader),
            patch(
                "libephemeris.fast_calc._apparent_icrs_cartesian",
                side_effect=KeyError(_BODY_MISS),
            ),
            patch("libephemeris.state.get_planets", side_effect=_sealed_get_planets),
        ):
            with pytest.raises(UnknownBodyError) as excinfo:
                _calc_gamma(_JD_2024)

        assert excinfo.value.body_id == 15

    def test_auto_mode_still_falls_back_to_skyfield(self):
        """Auto mode keeps the inline non-LEB fallback (unchanged behaviour)."""
        reader = MagicMock(name="LEBReader")

        with (
            patch("libephemeris.eclipse.get_calc_mode", return_value="auto"),
            patch("libephemeris.state.get_leb_reader", return_value=reader),
            patch(
                "libephemeris.fast_calc._apparent_icrs_cartesian",
                side_effect=ValueError(_RANGE_MISS),
            ),
            patch("libephemeris.state.get_planets", side_effect=_marker_get_planets),
        ):
            # Reaching the marker proves the fallback branch ran rather than
            # the range miss propagating out of auto mode.
            with pytest.raises(_ReachedSkyfield):
                _calc_gamma(_JD_2024)


class TestLocalEclipsePhasesSealedContract:
    """The local-phase probe must honour the same sealed contract."""

    def test_range_miss_raises_typed_error_not_runtime_error(self):
        reader = MagicMock(name="LEBReader")

        with (
            patch("libephemeris.eclipse.get_calc_mode", return_value="leb"),
            patch(
                "libephemeris.eclipse._get_leb_reader_safe",
                return_value=reader,
            ),
            patch(
                "libephemeris.fast_calc._topo_ecliptic",
                side_effect=ValueError(_RANGE_MISS),
            ),
            patch("libephemeris.state.get_planets", side_effect=_sealed_get_planets),
        ):
            with pytest.raises(EphemerisRangeError) as excinfo:
                _calculate_local_eclipse_phases(_JD_2024, *_GEOPOS)

        assert excinfo.value.start_jd == 2396758.5
        assert excinfo.value.end_jd == 2469807.5

    def test_body_miss_raises_unknown_body_not_runtime_error(self):
        reader = MagicMock(name="LEBReader")

        with (
            patch("libephemeris.eclipse.get_calc_mode", return_value="leb"),
            patch(
                "libephemeris.eclipse._get_leb_reader_safe",
                return_value=reader,
            ),
            patch(
                "libephemeris.fast_calc._topo_ecliptic",
                side_effect=KeyError(_BODY_MISS),
            ),
            patch("libephemeris.state.get_planets", side_effect=_sealed_get_planets),
        ):
            with pytest.raises(UnknownBodyError) as excinfo:
                _calculate_local_eclipse_phases(_JD_2024, *_GEOPOS)

        assert excinfo.value.body_id == 15

    def test_auto_mode_still_retries_on_the_skyfield_path(self):
        """Auto mode still reaches the wrapper's reader=None retry."""
        reader = MagicMock(name="LEBReader")

        with (
            patch("libephemeris.eclipse.get_calc_mode", return_value="auto"),
            patch(
                "libephemeris.eclipse._get_leb_reader_safe",
                return_value=reader,
            ),
            patch(
                "libephemeris.fast_calc._topo_ecliptic",
                side_effect=ValueError(_RANGE_MISS),
            ),
            patch("libephemeris.state.get_planets", side_effect=_marker_get_planets),
        ):
            with pytest.raises(_ReachedSkyfield):
                _calculate_local_eclipse_phases(_JD_2024, *_GEOPOS)
