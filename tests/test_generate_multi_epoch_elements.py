# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Tests for the audited multi-epoch element generator conventions."""

from __future__ import annotations

import sys
from types import SimpleNamespace

import pytest

from scripts import generate_multi_epoch_elements as generator


def test_reviewed_default_grid_has_all_81_ten_year_nodes() -> None:
    """The generator defaults must reproduce the shipped 1650--2450 grid."""
    epochs = generator.generate_epochs()

    assert generator.DEFAULT_SPACING == 10
    assert len(epochs) == 81
    assert epochs[0] == generator._year_to_jd(1650)
    assert epochs[-1] == generator._year_to_jd(2450)
    assert all(
        right - left == pytest.approx(10 * 365.25)
        for left, right in zip(epochs, epochs[1:])
    )


@pytest.mark.parametrize(
    ("start", "end", "spacing", "message"),
    [
        (2000, 1990, 10, "end_year"),
        (2000, 2010, 0, "spacing"),
        (2000, 2010, -1, "spacing"),
    ],
)
def test_epoch_grid_rejects_ambiguous_or_nonprogressing_ranges(
    start: int, end: int, spacing: int, message: str
) -> None:
    """Invalid grids fail explicitly instead of hanging or emitting no data."""
    with pytest.raises(ValueError, match=message):
        generator.generate_epochs(start, end, spacing)


class _FakeKernel:
    """Minimal type-21 kernel double exposing only audited metadata and state."""

    def __init__(self, segments: list[SimpleNamespace]) -> None:
        self.segments = segments
        self.compute_called = False
        self.closed = False

    def compute_type21(self, center: int, target: int, jd: float):
        """Return an axis-aligned state while recording that validation passed."""
        self.compute_called = True
        assert center == 10
        assert target == 12345
        assert jd == 2451545.0
        return (
            (generator.AU_KM, 0.0, 0.0),
            (generator.AU_KM / 86400.0, 0.0, 0.0),
        )

    def close(self) -> None:
        """Record the generator's unconditional resource cleanup."""
        self.closed = True


def _install_fake_spktype21(monkeypatch, kernel: _FakeKernel) -> None:
    """Install a deterministic local module double for the function-local import."""

    class _FakeSPKType21:
        @staticmethod
        def open(path: str) -> _FakeKernel:
            assert path == "model.bsp"
            return kernel

    monkeypatch.setitem(
        sys.modules,
        "spktype21",
        SimpleNamespace(SPKType21=_FakeSPKType21),
    )


@pytest.mark.parametrize(
    "segments",
    [
        [],
        [SimpleNamespace(center=0, target=12345, frame=1)],
        [SimpleNamespace(center=10, target=12345, frame=2)],
        [
            SimpleNamespace(center=10, target=12345, frame=1),
            SimpleNamespace(center=10, target=54321, frame=1),
        ],
    ],
)
def test_spk_state_rejects_non_heliocentric_or_non_j2000_metadata(
    monkeypatch, segments: list[SimpleNamespace]
) -> None:
    """Center, target, and frame metadata are validated before reading states."""
    kernel = _FakeKernel(segments)
    _install_fake_spktype21(monkeypatch, kernel)

    result = generator._get_spk_state_vector("model.bsp", 2451545.0)

    assert result is None
    assert not kernel.compute_called
    assert kernel.closed


def test_spk_state_accepts_only_declared_solar_j2000_segments(monkeypatch) -> None:
    """A valid type-21 state is converted from km and km/s to au and au/day."""
    kernel = _FakeKernel(
        [
            SimpleNamespace(center=10, target=12345, frame=1),
            SimpleNamespace(center=10, target=12345, frame=1),
        ]
    )
    _install_fake_spktype21(monkeypatch, kernel)

    result = generator._get_spk_state_vector("model.bsp", 2451545.0)

    assert result == pytest.approx((1.0, 0.0, 0.0, 1.0, 0.0, 0.0))
    assert kernel.compute_called
    assert kernel.closed
