# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Fictitious-body dispatch must not warn in sealed leb mode.

The runtime analytical model is the declared source for bodies 40-58; routing
to it is by-design dispatch, not a sealed-mode degradation. Kerykeion pins
mode="leb", so a WARNING here floods every downstream consumer that requests
Uranian points (kerykeion issue #240).
"""

from __future__ import annotations

import logging
import math

import pytest

import libephemeris as ephem
from libephemeris import state
from libephemeris.constants import FLG_SPEED, FLG_SWIEPH
from libephemeris.exceptions import FictitiousRuntimeDispatch
from libephemeris.logging_config import LOGGER_NAME
from libephemeris.time_utils import julday

FLAGS = FLG_SWIEPH | FLG_SPEED
JD = julday(2024, 3, 15, 7.25)  # inside the test LEB file range (2023-2028)

# The supported fictitious set: the eight Hamburg bodies plus White Moon.
SUPPORTED_FICTITIOUS = tuple(range(40, 48)) + (56,)


@pytest.fixture
def sealed_leb(test_leb_file):
    """Sealed leb mode backed by the small session test file."""
    prev_mode = state.get_calc_mode()
    ephem.set_calc_mode("leb")
    ephem.set_leb_file(test_leb_file)
    yield
    ephem.set_calc_mode(prev_mode)
    ephem.set_leb_file(None)


@pytest.fixture
def libephemeris_records():
    """Collect libephemeris logger records directly (propagate is False)."""
    logger = logging.getLogger(LOGGER_NAME)
    records: list[logging.LogRecord] = []

    class _Collector(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            records.append(record)

    handler = _Collector(level=logging.DEBUG)
    prev_level = logger.level
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)
    yield records
    logger.removeHandler(handler)
    logger.setLevel(prev_level)


@pytest.mark.parametrize("body_id", SUPPORTED_FICTITIOUS)
def test_no_warning_for_fictitious_dispatch_in_sealed_mode(
    sealed_leb, libephemeris_records, body_id
) -> None:
    """calc_ut on a fictitious body emits nothing at WARNING or above."""
    pos, _ = ephem.calc_ut(JD, body_id, FLAGS)
    assert all(math.isfinite(v) for v in pos)
    warnings = [r for r in libephemeris_records if r.levelno >= logging.WARNING]
    assert warnings == [], [r.getMessage() for r in warnings]


def test_dispatch_is_still_logged_at_debug(sealed_leb, libephemeris_records) -> None:
    """The by-design dispatch stays observable, one tier down."""
    ephem.calc_ut(JD, 40, FLAGS)
    debug_msgs = [
        r.getMessage()
        for r in libephemeris_records
        if r.levelno == logging.DEBUG and "runtime-model dispatch" in r.getMessage()
    ]
    assert debug_msgs, "expected the typed dispatch to be logged at DEBUG"


def test_dispatch_signal_is_an_ordinary_keyerror() -> None:
    """Callers catching (KeyError, ValueError) must keep working unchanged."""
    assert issubclass(FictitiousRuntimeDispatch, KeyError)
    err = FictitiousRuntimeDispatch(
        "Fictitious body 40 is calculated from its runtime analytical model"
    )
    assert "runtime analytical model" in str(err)
