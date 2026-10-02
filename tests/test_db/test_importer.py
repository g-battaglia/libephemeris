# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Offline coefficient extraction and opt-in PostgreSQL round-trip tests."""

from __future__ import annotations

import json
import os
import struct
from pathlib import Path
from uuid import uuid4

import pytest

import libephemeris as ephemeris
from libephemeris.db.contract import Series, segment_byte_count
from libephemeris.db.importer import (
    Artifact,
    connect_provisioner,
    import_artifacts,
    provision_schema,
)
from libephemeris.db.reader import DBReader
from libephemeris.db.store import PostgresStore
from libephemeris.exceptions import DBDataError
from libephemeris.leb_reader import open_leb

from .artifacts import make_artifacts


def test_all_native_formats_extract_identical_coefficients(tmp_path):
    """Cover bounded LEB1, monolithic LEB2 and chunked LEB2 decoding.

    Args:
        tmp_path: Disposable artifact directory.
    """
    expected = None
    for path in make_artifacts(tmp_path):
        with open_leb(str(path)) as native:
            artifact = Artifact(path, native)
            series = list(artifact.body_series())
            segments = list(artifact.body_segments(series[0]))
            assert len(segments) == 2
            assert [index for index, payload in segments] == [0, 1]
            assert all(
                len(payload) == segment_byte_count(series[0])
                for index, payload in segments
            )
            if expected is None:
                expected = segments
            else:
                assert segments == expected
            auxiliary = dict(artifact.auxiliary_sections())
            assert {0, 2, 3, 4}.issubset(auxiliary)
            assert artifact.manifest_record()["sha256"]


def test_portable_fixture_documents_exact_byte_layout():
    """Check the language-neutral fixture against mathematical identities."""
    path = Path(__file__).parent / "fixtures/contract_v1.json"
    fixture = json.loads(path.read_text())
    series = Series(**fixture["series"])
    payload = bytes.fromhex(fixture["payload_hex"])
    assert len(payload) == segment_byte_count(series)
    assert struct.unpack("<9d", payload) == (
        1.0,
        2.0,
        3.0,
        -4.0,
        0.5,
        0.0,
        7.0,
        0.0,
        0.0,
    )
    from libephemeris.leb_reader import _clenshaw_with_derivative

    positions = []
    velocities = []
    coefficients = struct.unpack("<9d", payload)
    for component in range(3):
        offset = component * 3
        value, derivative = _clenshaw_with_derivative(
            coefficients[offset : offset + 3], 0.0
        )
        positions.append(value)
        velocities.append(derivative)
    assert positions == fixture["evaluation"]["position"]
    assert velocities == fixture["evaluation"]["velocity"]


@pytest.mark.db_integration
def test_postgres_import_roundtrip_and_immutable_publication(tmp_path):
    """Verify real SQL, import idempotence and read-only transport on opt-in DB.

    Args:
        tmp_path: Disposable source directory.
    """
    dsn = os.environ.get("LIBEPHEMERIS_TEST_DB_URL")
    if not dsn:
        pytest.skip("Set LIBEPHEMERIS_TEST_DB_URL to a disposable PostgreSQL database")
    pytest.importorskip("psycopg")
    previous_policy = ephemeris.get_configured_network_policy()
    ephemeris.set_network_policy("allow")
    try:
        with connect_provisioner(dsn) as connection:
            with connection.transaction():
                provision_schema(connection)
            for path in make_artifacts(tmp_path):
                dataset_id = str(uuid4())
                import_artifacts(connection, dataset_id, [path], "base")
                import_artifacts(connection, dataset_id, [path], "base")
                store = PostgresStore(dsn)
                reader = DBReader(store, dataset_id)
                try:
                    with open_leb(str(path)) as native:
                        for jd in [10.0, 11.0, 12.0, 13.0, 14.0]:
                            assert reader.eval_body(0, jd) == native.eval_body(0, jd)
                            assert reader.eval_nutation(jd) == native.eval_nutation(jd)
                            assert reader.delta_t(jd) == native.delta_t(jd)
                        assert reader.get_star(42) == native.get_star(42)
                    import psycopg

                    with pytest.raises(psycopg.Error, match="immutable"):
                        connection.execute(
                            "DELETE FROM libephemeris.segments WHERE dataset_id=%s",
                            (dataset_id,),
                        )
                finally:
                    reader.close()
                    store.close()
                with pytest.raises(DBDataError):
                    import_artifacts(connection, dataset_id, [path], "medium")
    finally:
        ephemeris.set_network_policy(previous_policy)
