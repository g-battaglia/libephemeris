# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Regression tests for DB lifetime, source errors and persistent-result leaks."""

from __future__ import annotations

import math

import pytest

import libephemeris as ephemeris
from libephemeris import eclipse, fast_calc, fixed_stars, planets
from libephemeris.db import backend
from libephemeris.db.contract import NUTATION_ID
from libephemeris.db.reader import DBReader
from libephemeris.exceptions import DBDataError, UnknownBodyError

from .conftest import DATASET_ID


@pytest.mark.parametrize(
    "function, arguments",
    [
        (ephemeris.rise_trans, (2451545.0, ephemeris.SUN, 1, (0.0, 45.0, 0.0))),
        (
            ephemeris.rise_trans_true_hor,
            (2451545.0, ephemeris.SUN, 1, (0.0, 45.0, 0.0)),
        ),
        (ephemeris.sol_eclipse_max_time, (2460409.0,)),
        (
            ephemeris.get_orbital_elements,
            (2451545.0, ephemeris.MARS, ephemeris.FLG_HELCTR),
        ),
        (
            ephemeris.get_orbital_elements_ut,
            (2451545.0, ephemeris.MARS, ephemeris.FLG_HELCTR),
        ),
        (ephemeris.orbit_max_min_true_distance, (2451545.0, ephemeris.MARS, 0)),
        (ephemeris.calc_true_lunar_node, (2451545.0,)),
        (ephemeris.calc_true_lilith, (2451545.0,)),
        (fixed_stars.calc_fixed_star_position, (ephemeris.SPICA_STAR, 2451545.0)),
        (fixed_stars.calc_fixed_star_velocity, (ephemeris.SPICA_STAR, 2451545.0)),
        (eclipse._calc_gamma, (2460409.0,)),
        (eclipse._calc_local_eclipse_max_time, (2460409.0, 45.0, 0.0, 0.0, 0.125)),
    ],
)
def test_other_public_calculations_release_all_inputs(
    db_runtime, monkeypatch, function, arguments
):
    """Verify one owner and immediate cleanup beyond calc/calc_ut entry points.

    Args:
        db_runtime: Native coefficient transport double.
        monkeypatch: Pytest restoration helper.
        function: Public entry point under test.
        arguments: Arguments exercising its stored-state calculations.
    """
    readers = []

    def create_reader(store, dataset_id):
        """Observe allocated input owners without changing their behavior.

        Args:
            store: Operation transport.
            dataset_id: Explicit dataset identity.

        Returns:
            Newly allocated reader, retained only by this test for inspection.
        """
        reader = DBReader(store, dataset_id)
        readers.append(reader)
        return reader

    monkeypatch.setattr(backend, "DBReader", create_reader)
    function(*arguments)
    assert len(readers) == 1
    assert readers[0]._closed
    assert readers[0]._segments == {}
    assert readers[0]._series == {}
    assert getattr(backend._operation_state, "reader", None) is None
    assert not isinstance(getattr(fast_calc._active_local, "reader", None), DBReader)


@pytest.mark.parametrize("flags", [0, ephemeris.FLG_NOGDEFL])
def test_star_ayanamsha_never_uses_global_result_cache(db_runtime, monkeypatch, flags):
    """Cover both exact apparent anchors and interpolated deflection-free nodes.

    Args:
        db_runtime: Native coefficient transport double.
        monkeypatch: Pytest restoration helper.
        flags: Anchor flags selecting exact or interpolation evaluation.
    """

    def forbidden(*args, **kwargs):
        """Reject access to the persistent fixed-star result LRU.

        Args:
            *args: Unused cache key fields.
            **kwargs: Unused cache options.

        Raises:
            AssertionError: A DB calculation tried to use a global result cache.
        """
        raise AssertionError("DB-derived anchor entered global result cache")

    monkeypatch.setattr(planets, "_star_position_ecliptic_cached", forbidden)
    ephemeris.set_sid_mode(ephemeris.SIDM_TRUE_CITRA)
    first = ephemeris.get_ayanamsa_ex_ut(2451545.0, flags)
    first_reads = db_runtime.segment_calls
    second = ephemeris.get_ayanamsa_ex_ut(2451545.0, flags)
    assert first == second
    assert first_reads > 0
    assert db_runtime.segment_calls > first_reads
    assert db_runtime.metadata_calls == 2


def test_absent_optional_sections_keep_existing_local_models(db_runtime, monkeypatch):
    """Missing optional inputs differ from corrupt declared coefficient channels.

    Args:
        db_runtime: Complete native coefficient transport double.
        monkeypatch: Pytest restoration helper.
    """
    original_metadata = db_runtime.metadata

    def metadata_without_nutation(dataset_id):
        """Describe a valid dataset that never declared stored nutation.

        Args:
            dataset_id: Requested dataset.

        Returns:
            Body metadata without the optional nutation series.
        """
        metadata = original_metadata(dataset_id)
        del metadata[NUTATION_ID]
        return metadata

    def no_delta_t(dataset_id, jd):
        """Describe an absent optional table, not a transport error.

        Args:
            dataset_id: Requested dataset.
            jd: Requested epoch.

        Returns:
            No persisted samples.
        """
        return []

    monkeypatch.setattr(db_runtime, "metadata", metadata_without_nutation)
    monkeypatch.setattr(db_runtime, "delta_t_points", no_delta_t)
    position, name, flags = ephemeris.fixstar_ut("Spica", 2451545.0)
    assert all(math.isfinite(value) for value in position)
    assert math.isfinite(ephemeris.get_ayanamsa_ut(2451545.0))
    assert getattr(backend._operation_state, "reader", None) is None


def test_missing_body_in_search_is_a_public_source_error(db_runtime, monkeypatch):
    """Translate DB missing-channel signals without a file/Horizons fallback.

    Args:
        db_runtime: Complete native coefficient transport double.
        monkeypatch: Pytest restoration helper.
    """
    original_metadata = db_runtime.metadata

    def metadata_without_sun(dataset_id):
        """Omit a required target from the source inventory.

        Args:
            dataset_id: Requested dataset.

        Returns:
            Valid metadata for all channels except the Sun.
        """
        metadata = original_metadata(dataset_id)
        del metadata[ephemeris.SUN]
        return metadata

    monkeypatch.setattr(db_runtime, "metadata", metadata_without_sun)
    with pytest.raises(UnknownBodyError, match="DB dataset"):
        ephemeris.rise_trans(2451545.0, ephemeris.SUN, 1, (0.0, 45.0, 0.0))
    assert getattr(backend._operation_state, "reader", None) is None
    assert not isinstance(getattr(fast_calc._active_local, "reader", None), DBReader)


def test_closed_reader_cannot_resume_auxiliary_reads(artifact_store):
    """Ensure cleanup is final even for methods that do not load coefficients.

    Args:
        artifact_store: Native coefficient transport double.
    """
    reader = DBReader(artifact_store, DATASET_ID)
    reader.close()
    with pytest.raises(DBDataError, match="closed"):
        reader.delta_t(2451545.0)
    with pytest.raises(DBDataError, match="closed"):
        reader.get_star(42)
    with pytest.raises(DBDataError, match="closed"):
        reader.eval_nutation(2451545.0)
    with pytest.raises(DBDataError, match="closed"):
        reader.eval_body(ephemeris.SUN, 2451545.0)
