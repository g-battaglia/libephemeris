# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Unit tests for libephemeris.angles (astrological angle calculations)."""

from __future__ import annotations

import pytest

from libephemeris import angles
from libephemeris.angles import calc_angles, get_angle_value
from libephemeris.constants import (
    ANTIVERTEX,
    ASCENDANT,
    DESCENDANT,
    IC,
    MC,
    VERTEX,
    _MC_ANGLE_ID,
    _VERTEX_ANGLE_ID,
)

_JD = 2451545.0  # J2000.0
_LAT = 41.9  # Rome
_LON = 12.5

_EXPECTED_KEYS = {
    "Ascendant",
    "MC",
    "ARMC",
    "Vertex",
    "Equatorial_Ascendant",
    "Descendant",
    "IC",
    "AntiVertex",
}


class TestCalcAngles:
    """Exercise the real houses pipeline (no planetary ephemeris needed)."""

    def test_returns_all_keys_in_range(self):
        result = calc_angles(_JD, _LAT, _LON)
        assert set(result) == _EXPECTED_KEYS
        for val in result.values():
            assert isinstance(val, float)
            assert 0.0 <= val < 360.0

    def test_derived_points_are_opposite(self):
        r = calc_angles(_JD, _LAT, _LON)
        assert r["Descendant"] == pytest.approx((r["Ascendant"] + 180.0) % 360.0)
        assert r["IC"] == pytest.approx((r["MC"] + 180.0) % 360.0)
        assert r["AntiVertex"] == pytest.approx((r["Vertex"] + 180.0) % 360.0)

    def test_short_ascmc_defaults_equatorial_ascendant(self, monkeypatch):
        # Cover the `ascmc[4] if len(ascmc) > 4 else 0.0` false branch: an
        # ascmc tuple with <=4 entries must yield Equatorial_Ascendant == 0.0.
        def fake(*_args, **_kwargs):
            return ((), (10.0, 20.0, 30.0, 40.0), True, None)

        monkeypatch.setattr(angles, "houses_with_fallback", fake)
        r = calc_angles(_JD, _LAT, _LON)
        assert r["Equatorial_Ascendant"] == 0.0
        assert r["Ascendant"] == 10.0
        assert r["MC"] == 20.0
        assert r["Vertex"] == 40.0
        assert r["AntiVertex"] == (40.0 + 180.0) % 360.0


class TestGetAngleValue:
    @pytest.fixture(autouse=True)
    def _patch_houses(self, monkeypatch):
        # Deterministic ascmc so each mapping resolves to a unique value
        # (index 4 present -> exercises the len(ascmc) > 4 true branch too).
        def fake(*_args, **_kwargs):
            return ((), (11.0, 22.0, 33.0, 44.0, 55.0), True, None)

        monkeypatch.setattr(angles, "houses_with_fallback", fake)

    def test_each_recognized_angle(self):
        assert get_angle_value(ASCENDANT, _JD, _LAT, _LON) == 11.0
        assert get_angle_value(_MC_ANGLE_ID, _JD, _LAT, _LON) == 22.0
        assert get_angle_value(DESCENDANT, _JD, _LAT, _LON) == (11.0 + 180.0) % 360.0
        assert get_angle_value(IC, _JD, _LAT, _LON) == (22.0 + 180.0) % 360.0
        assert get_angle_value(_VERTEX_ANGLE_ID, _JD, _LAT, _LON) == 44.0
        assert get_angle_value(ANTIVERTEX, _JD, _LAT, _LON) == (44.0 + 180.0) % 360.0
        # Public aliases MC (=1) and VERTEX (=3) resolve to the same angles.
        assert get_angle_value(MC, _JD, _LAT, _LON) == 22.0
        assert get_angle_value(VERTEX, _JD, _LAT, _LON) == 44.0

    def test_unrecognized_angle_returns_zero(self):
        assert get_angle_value(999_999, _JD, _LAT, _LON) == 0.0
