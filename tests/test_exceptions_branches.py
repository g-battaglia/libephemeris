# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for libephemeris.exceptions helpers."""

from __future__ import annotations

import pytest

from libephemeris import state
from libephemeris.exceptions import (
    EphemerisRangeError,
    SPKNotFoundError,
    validate_jd_range,
)


class TestSPKNotFoundFromFilepath:
    def test_message_without_body_id_uses_placeholders(self):
        err = SPKNotFoundError.from_filepath("/tmp/missing.bsp")
        msg = str(err)
        assert "missing.bsp" in msg
        # else-branches when body_id / body_name are absent:
        assert "<body_id_or_name>" in msg
        assert "<body_id>" in msg

    def test_message_with_body_id_and_name_uses_concrete_values(self):
        err = SPKNotFoundError.from_filepath(
            "/tmp/missing.bsp", body_name="ceres", body_id="2000001"
        )
        msg = str(err)
        assert "2000001" in msg
        assert "CERES" in msg


class TestValidateJdRange:
    @pytest.fixture(autouse=True)
    def _no_real_ephemeris(self, monkeypatch):
        # Avoid loading any real ephemeris while exercising the range logic.
        monkeypatch.setattr(state, "get_planets", lambda *a, **k: None)

    def test_zero_range_skips_validation(self, monkeypatch):
        # start == end == 0 -> the "no range info" early return (no raise).
        monkeypatch.setattr(
            state, "get_current_file_data", lambda _i: ("any.bsp", 0.0, 0.0, 440)
        )
        assert validate_jd_range(3_000_000.0) is None

    def test_out_of_range_uses_denum_filename(self, monkeypatch):
        # path empty but denum set -> ephemeris_file = "de{denum}.bsp".
        monkeypatch.setattr(
            state,
            "get_current_file_data",
            lambda _i: ("", 2451545.0, 2451546.0, 440),
        )
        with pytest.raises(EphemerisRangeError) as exc:
            validate_jd_range(3_000_000.0)
        assert "de440.bsp" in str(exc.value)

    def test_out_of_range_without_file_info(self, monkeypatch):
        # No path and no denum -> the error message omits the file line.
        monkeypatch.setattr(
            state,
            "get_current_file_data",
            lambda _i: (None, 2451545.0, 2451546.0, 0),
        )
        with pytest.raises(EphemerisRangeError) as exc:
            validate_jd_range(3_000_000.0)
        assert "Ephemeris file:" not in str(exc.value)
