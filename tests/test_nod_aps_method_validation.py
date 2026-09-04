# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Regression tests for issue #58: an unknown nod_aps method must be refused.

``nod_aps()`` and ``nod_aps_ut()`` did not validate ``method``. The
method-bit precedence rule resolves a mask with no model bit to
``NODBIT_MEAN``, so every unrecognised value — 8, -1, anything with a stray
bit — reached the calculation and came back as a mean-element answer. The
caller got a plausible, well-formed, wrong result with no signal that the
model they asked for had been ignored.

Zero is not an unknown value: it is the documented absence of a model bit
and keeps resolving to ``NODBIT_MEAN``.
"""

from __future__ import annotations

import pytest

import libephemeris as eph
from libephemeris.constants import (
    NODBIT_FOPOINT,
    NODBIT_MEAN,
    NODBIT_OSCU,
    NODBIT_OSCU_BAR,
)
from libephemeris.exceptions import Error, InputValidationError

_JD = 2451545.0

_VALID_METHODS = [
    0,
    NODBIT_MEAN,
    NODBIT_OSCU,
    NODBIT_OSCU_BAR,
    NODBIT_FOPOINT,
    NODBIT_MEAN | NODBIT_FOPOINT,
    NODBIT_OSCU | NODBIT_FOPOINT,
    NODBIT_OSCU_BAR | NODBIT_FOPOINT,
    NODBIT_MEAN | NODBIT_OSCU | NODBIT_OSCU_BAR | NODBIT_FOPOINT,
]

# Values that carry a bit outside {1, 2, 4, 256}, or are negative.
_UNKNOWN_METHODS = [8, 16, 32, 64, 128, 255, 264, 512, 1024, -1, -2]


@pytest.mark.unit
@pytest.mark.parametrize("method", _UNKNOWN_METHODS)
@pytest.mark.parametrize("func", [eph.nod_aps_ut, eph.nod_aps], ids=["ut", "et"])
def test_unknown_method_is_refused(func, method):
    """The reported case: a wrong constant used to return the mean result."""
    with pytest.raises(InputValidationError):
        func(_JD, eph.MARS, method)


@pytest.mark.unit
@pytest.mark.parametrize("method", _VALID_METHODS)
@pytest.mark.parametrize("func", [eph.nod_aps_ut, eph.nod_aps], ids=["ut", "et"])
def test_documented_methods_are_accepted(func, method):
    """Narrowing the input must not reject any documented combination."""
    result = func(_JD, eph.MARS, method)
    assert len(result) == 4
    assert all(len(component) == 6 for component in result)


@pytest.mark.unit
@pytest.mark.parametrize("method", [3.5, "mean", None, [1], object()])
@pytest.mark.parametrize("func", [eph.nod_aps_ut, eph.nod_aps], ids=["ut", "et"])
def test_non_integer_method_is_refused_with_the_library_error(func, method):
    """A non-integer must not leak TypeError from a bitwise operation."""
    with pytest.raises(InputValidationError):
        func(_JD, eph.MARS, method)


@pytest.mark.unit
def test_zero_still_means_mean():
    """Zero is the absence of a model bit, not an unknown one."""
    assert eph.nod_aps_ut(_JD, eph.MARS, 0) == eph.nod_aps_ut(
        _JD, eph.MARS, NODBIT_MEAN
    )


@pytest.mark.unit
def test_error_is_catchable_as_the_library_base_error():
    """The typed error stays inside the documented hierarchy."""
    with pytest.raises(Error):
        eph.nod_aps_ut(_JD, eph.MARS, 8)
    assert issubclass(InputValidationError, Error)


@pytest.mark.unit
def test_message_names_the_valid_flags():
    """The message has to say what a valid method looks like."""
    with pytest.raises(InputValidationError) as excinfo:
        eph.nod_aps_ut(_JD, eph.MARS, 8)
    message = str(excinfo.value)
    for name in ("NODBIT_MEAN", "NODBIT_OSCU", "NODBIT_OSCU_BAR", "NODBIT_FOPOINT"):
        assert name in message


@pytest.mark.unit
def test_boolean_is_accepted_as_the_integer_it_is():
    """bool is an int subclass; True is NODBIT_MEAN and must not be refused."""
    assert eph.nod_aps_ut(_JD, eph.MARS, True) == eph.nod_aps_ut(
        _JD, eph.MARS, NODBIT_MEAN
    )
