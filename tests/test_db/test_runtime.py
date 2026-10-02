# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Public API parity and source-purity tests, independent of PostgreSQL."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

import pytest

import libephemeris as ephemeris
from libephemeris import fast_calc, state
from libephemeris.db import backend
from libephemeris.db.reader import DBReader
from libephemeris.exceptions import ConfigurationError, DBDataError, NetworkSealedError
from libephemeris.db.policy import require_database
from libephemeris.net import require_network

from .conftest import BUNDLED_FILE, DATASET_ID


@pytest.mark.parametrize("body", [0, 1, 2, 4, 5, 9, 10, 11, 12, 14, 21, 40])
@pytest.mark.parametrize(
    "flags",
    [
        ephemeris.FLG_SPEED,
        ephemeris.FLG_TRUEPOS,
        ephemeris.FLG_SPEED | ephemeris.FLG_EQUATORIAL,
        ephemeris.FLG_SPEED | ephemeris.FLG_SIDEREAL,
    ],
)
def test_public_db_matches_file_backend(db_runtime, body, flags):
    """Compare both modes ephemerally using identical project-native inputs.

    Args:
        db_runtime: Instrumented native transport double.
        body: Public body identifier.
        flags: Requested representation and calculation options.
    """
    ephemeris.set_leb_file(str(BUNDLED_FILE))
    ephemeris.set_calc_mode("leb")
    expected_ut = ephemeris.calc_ut(2451545.0, body, flags)
    expected_tt = ephemeris.calc(2451545.0, body, flags)
    ephemeris.set_calc_mode("db")
    actual_ut = ephemeris.calc_ut(2451545.0, body, flags)
    actual_tt = ephemeris.calc(2451545.0, body, flags)
    assert actual_ut == expected_ut
    assert actual_tt == expected_tt


def test_inputs_are_read_again_and_released(db_runtime):
    """Repeated requests use DB reads without retaining records or frames.

    Args:
        db_runtime: Instrumented native transport double.
    """
    ephemeris.calc_ut(2451545.0, ephemeris.MARS, ephemeris.FLG_SPEED)
    assert db_runtime.metadata_calls == 1
    first_count = db_runtime.segment_calls
    assert first_count <= 3
    ephemeris.calc_ut(2451545.0, ephemeris.MARS, ephemeris.FLG_SPEED)
    assert db_runtime.metadata_calls == 2
    assert db_runtime.segment_calls == 2 * first_count
    assert getattr(backend._operation_state, "reader", None) is None
    assert not fast_calc._leb_frame_cache
    assert not isinstance(getattr(fast_calc._active_local, "reader", None), DBReader)
    assert state.get_leb_reader() is None
    assert state.get_horizons_client() is None


def test_runtime_cannot_open_files_or_other_sources(db_runtime, monkeypatch):
    """Fail if file discovery or JPL/Horizons resolution is attempted.

    Args:
        db_runtime: Instrumented native transport double.
        monkeypatch: Pytest patch helper.
    """

    def forbidden(*args, **kwargs):
        """Reject access to another source.

        Args:
            *args: Unused source arguments.
            **kwargs: Unused source options.

        Raises:
            AssertionError: A non-DB source was requested.
        """
        raise AssertionError("Non-DB source accessed")

    monkeypatch.setattr(state, "_get_leb_reader_locked", forbidden)
    monkeypatch.setattr(state, "get_planets", forbidden)
    monkeypatch.setattr(state, "_get_or_create_horizons_client", forbidden)
    trace = ephemeris.start_tracing()
    try:
        ephemeris.calc_ut(2451545.0, ephemeris.MARS, ephemeris.FLG_SPEED)
        assert ephemeris.get_trace_results()[ephemeris.MARS] == "DB"
        context = ephemeris.EphemerisContext()
        context.set_leb_file("/nonexistent/context.leb")
        assert context.get_leb_reader() is None
        assert context.calc_ut(
            2451545.0, ephemeris.MARS, ephemeris.FLG_SPEED
        ) == ephemeris.calc_ut(2451545.0, ephemeris.MARS, ephemeris.FLG_SPEED)
        assert ephemeris.get_current_file_data(0) == ("", 0.0, 0.0, 0)
        coverage = ephemeris.get_body_coverage(ephemeris.MARS)
        assert coverage.source == "DB"
        assert coverage.data_file is None
    finally:
        trace.var.reset(trace)


def test_failure_releases_operation_and_never_falls_back(db_runtime, monkeypatch):
    """Preserve DB failures and release ownership through exceptional exits.

    Args:
        db_runtime: Instrumented transport double.
        monkeypatch: Pytest patch helper.
    """

    def fail(dataset_id, keys):
        """Simulate incomplete persisted inputs.

        Args:
            dataset_id: Requested version.
            keys: Requested segments.

        Raises:
            DBDataError: Simulated corruption.
        """
        raise DBDataError("Missing persisted segment")

    monkeypatch.setattr(db_runtime, "fetch_segments", fail)
    with pytest.raises(DBDataError, match="Missing"):
        ephemeris.calc_ut(2451545.0, ephemeris.MARS, ephemeris.FLG_SPEED)
    assert getattr(backend._operation_state, "reader", None) is None
    assert not isinstance(getattr(fast_calc._active_local, "reader", None), DBReader)


def test_context_observer_and_thread_inputs_are_isolated(db_runtime):
    """Concurrent contexts share only the transport, not operation inputs.

    Args:
        db_runtime: Native transport double.
    """

    def calculate(longitude):
        """Run one topocentric calculation with explicit context settings.

        Args:
            longitude: Observer longitude in degrees.

        Returns:
            The calculation result and post-operation reader slot.
        """
        context = ephemeris.EphemerisContext()
        context.set_topo(longitude, 45.0, 0.0)
        result = context.calc_ut(
            2451545.0, ephemeris.MOON, ephemeris.FLG_SPEED | ephemeris.FLG_TOPOCTR
        )
        return result, getattr(backend._operation_state, "reader", None)

    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(calculate, [0.0, 90.0]))
    assert results[0][0] != results[1][0]
    assert results[0][1] is None
    assert results[1][1] is None


def test_db_network_policy_distinguishes_database_and_http(db_runtime):
    """Automatic DB policy permits only explicit database transport.

    Args:
        db_runtime: Configured DB runtime.
    """
    require_database()
    with pytest.raises(NetworkSealedError):
        require_network("Unexpected download")
    ephemeris.set_network_policy("sealed")
    with pytest.raises(NetworkSealedError):
        require_database()


def test_configuration_uuid_and_secrets_are_validated(monkeypatch):
    """Configuration failures do not quote DSN secrets.

    Args:
        monkeypatch: Pytest patch helper.
    """
    secret = "postgresql://user:secret_password@host/database"
    with pytest.raises(ConfigurationError) as failure:
        backend.set_db_config(secret, "not-a-uuid")
    assert "secret_password" not in str(failure.value)
    backend.set_db_config(None)
    monkeypatch.setenv("LIBEPHEMERIS_DB_URL", secret)
    monkeypatch.setenv("LIBEPHEMERIS_DB_DATASET", DATASET_ID)
    assert backend.get_db_config() == (secret, DATASET_ID)
