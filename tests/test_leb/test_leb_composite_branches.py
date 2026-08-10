# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for CompositeLEBReader error / companion paths."""

from __future__ import annotations

import pytest

from libephemeris.leb_composite import CompositeLEBReader


class FakeReader:
    """Minimal LEB reader stand-in for composite dispatch tests."""

    def __init__(
        self,
        bodies=None,
        path="fake.leb",
        jd_range=(0.0, 1.0),
        nutation=False,
        delta_t=None,
        stars=None,
        cool_raises=False,
        close_raises=False,
    ):
        self._bodies = bodies if bodies is not None else {0: object()}
        self.path = path
        self.jd_range = jd_range
        self._delta_t_jds = delta_t if delta_t is not None else []
        self._stars = stars if stars is not None else {}
        self._nutation = nutation
        self._cool_raises = cool_raises
        self._close_raises = close_raises

    def has_nutation(self):
        return self._nutation

    def eval_body(self, body_id, jd):
        return ((1.0, 2.0, 3.0), (0.1, 0.2, 0.3))

    def eval_nutation(self, jd):
        return (0.0, 0.0)

    def delta_t(self, jd):
        return 69.0

    def get_star(self, star_id):
        return "star"

    def warm(self, a, b):
        pass

    def cool(self):
        if self._cool_raises:
            raise OSError("cool failed")

    def close(self):
        if self._close_raises:
            raise ValueError("close failed")


def test_empty_readers_raises():
    with pytest.raises(ValueError, match="at least one reader"):
        CompositeLEBReader([])


def test_path_property_uses_first_reader():
    comp = CompositeLEBReader([FakeReader(path="primary.leb")])
    assert comp.path == "primary.leb"


def test_cool_swallows_reader_errors():
    comp = CompositeLEBReader([FakeReader(cool_raises=True)])
    comp.cool()  # must not raise


def test_aux_data_absent_raises_clean_errors():
    comp = CompositeLEBReader([FakeReader(nutation=False, delta_t=[], stars={})])
    assert comp.has_nutation() is False
    with pytest.raises(ValueError, match="No nutation data"):
        comp.eval_nutation(2451545.0)
    with pytest.raises(ValueError, match="No Delta-T data"):
        comp.delta_t(2451545.0)
    with pytest.raises(KeyError):
        comp.get_star(1)


def test_close_swallows_reader_errors():
    comp = CompositeLEBReader([FakeReader(close_raises=True)])
    comp.close()  # must not raise
    assert comp._readers == []


def test_aux_data_present_success_paths():
    reader = FakeReader(
        bodies={0: object()},
        nutation=True,
        delta_t=[2451545.0],
        stars={1: "vega"},
    )
    comp = CompositeLEBReader([reader])
    assert comp.has_nutation() is True
    assert comp.eval_nutation(2451545.0) == (0.0, 0.0)
    assert comp.delta_t(2451545.0) == 69.0
    assert comp.get_star(1) == "star"
    # body dispatch
    assert comp.has_body(0) is True
    pos, vel = comp.eval_body(0, 2451545.0)
    assert pos == (1.0, 2.0, 3.0)
    with pytest.raises(KeyError):
        comp.eval_body(999, 2451545.0)
    # warm fans out to all readers without error
    comp.warm(0.0, 1.0)
    # widest jd range across readers
    assert comp.jd_range == (0.0, 1.0)


def test_from_file_single_non_group_file(tmp_path, monkeypatch):
    import libephemeris.leb_reader as lr

    primary = tmp_path / "ephemeris.leb"
    primary.write_bytes(b"x")

    monkeypatch.setattr(lr, "open_leb", lambda path: FakeReader(path=path))
    # name has no group suffix -> companion discovery is skipped entirely.
    comp = CompositeLEBReader.from_file_with_companions(str(primary))
    assert len(comp._readers) == 1


def test_from_file_with_companions_filters_and_skips_bad(tmp_path, monkeypatch):
    import libephemeris.leb_reader as lr

    (tmp_path / "base_core.leb2").write_bytes(b"x")
    (tmp_path / "base_asteroids.leb2").write_bytes(b"x")  # opens but fails
    (tmp_path / "base_weird.leb2").write_bytes(b"x")  # non-group suffix, filtered
    primary = tmp_path / "base_core.leb2"

    def fake_open(path):
        if path.endswith("base_core.leb2"):
            return FakeReader(bodies={0: object()}, path=path)
        if path.endswith("base_asteroids.leb2"):
            raise ValueError("corrupt companion")
        return FakeReader(path=path)

    monkeypatch.setattr(lr, "open_leb", fake_open)

    comp = CompositeLEBReader.from_file_with_companions(str(primary))
    # asteroids failed to open, weird filtered by suffix -> only primary remains.
    assert len(comp._readers) == 1


def test_reviewed_core_attaches_every_same_prefix_sibling(tmp_path):
    """Discovery: a reviewed core attaches all its same-prefix siblings.

    Exercises the production wiring (_has_pinned_sibling_companions +
    pinned_only=True) end to end. Since runtime content verification was
    removed, presence is the only requirement: every same-prefix group file
    on disk attaches. Ported from the retired uranians-companion suite to
    the asteroids/apogee groups.
    """
    import os
    import shutil

    import libephemeris as ephem
    from libephemeris import state

    bundled_core = os.path.join(
        os.path.dirname(state.__file__), "data", "leb2", "base_core.leb2"
    )
    core = tmp_path / "base_core.leb2"
    shutil.copy(bundled_core, core)
    # Same-prefix manifest-named siblings; bytes are irrelevant (a valid LEB2
    # payload is enough), so reuse the core bytes for both.
    shutil.copy(bundled_core, tmp_path / "base_asteroids.leb2")
    shutil.copy(bundled_core, tmp_path / "base_apogee.leb2")

    prev_mode = state.get_calc_mode()
    prev_leb_file = state._LEB_FILE
    ephem.set_calc_mode("auto")
    ephem.set_leb_file(str(core))
    try:
        reader = state.get_leb_reader()
        assert len(reader._readers) == 3
    finally:
        ephem.set_leb_file(prev_leb_file)
        ephem.set_calc_mode(prev_mode)


def test_legacy_uranians_companions_are_never_attached(tmp_path):
    """A leftover {tier}_uranians file is ignored by every discovery path.

    Pins the 3.1.0 contract: fictitious provenance is invariantly
    "Analytical" — a stale companion must not resurface source="LEB"
    coverage for bodies 40-47 (attaching it would, via get_body_coverage).
    """
    import os
    import shutil

    import libephemeris as ephem
    from libephemeris import state

    bundled_core = os.path.join(
        os.path.dirname(state.__file__), "data", "leb2", "base_core.leb2"
    )
    core = tmp_path / "base_core.leb2"
    shutil.copy(bundled_core, core)
    # Any valid LEB2 payload works: the skip must key on the filename.
    shutil.copy(bundled_core, tmp_path / "base_uranians.leb2")

    # from_file_with_companions (both trust modes)
    for pinned_only in (False, True):
        comp = CompositeLEBReader.from_file_with_companions(
            str(core), pinned_only=pinned_only
        )
        try:
            assert len(comp._readers) == 1
        finally:
            comp.close()

    # from_directory
    comp = CompositeLEBReader.from_directory(str(tmp_path))
    try:
        assert len(comp._readers) == 1
    finally:
        comp.close()

    # End to end: activated reader reports no coverage for fictitious ids.
    prev_mode = state.get_calc_mode()
    prev_leb_file = state._LEB_FILE
    ephem.set_calc_mode("leb")
    ephem.set_leb_file(str(core))
    try:
        for body_id in (40, 47, 56):
            assert ephem.get_body_coverage(body_id) is None
    finally:
        ephem.set_leb_file(prev_leb_file)
        ephem.set_calc_mode(prev_mode)
