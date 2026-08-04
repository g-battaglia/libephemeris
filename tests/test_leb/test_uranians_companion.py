# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Regenerated Hamburg-body (uranians) LEB companion: generation, trust, dispatch.

Covers the full chain introduced with the regenerated ``uranians`` group:

- LEB1 partial generation samples only the runtime Neely propagation;
- LEB2 conversion authenticates the exact {40..47} inventory;
- companion attach and calculation sourcing require the manifest SHA-256 pin
  (arbitrary same-named artifacts stay unused);
- the public file-backed output is equivalent to the pure runtime model;
- fail-closed behavior for IDs 48-58 and the fast_calc guard are unchanged.
"""

from __future__ import annotations

import math
import os
import shutil

import pytest

import libephemeris as ephem
from libephemeris import download, state
from libephemeris.constants import (
    FLG_EQUATORIAL,
    FLG_HELCTR,
    FLG_SPEED,
    FLG_SWIEPH,
    FLG_XYZ,
)
from libephemeris.exceptions import UnknownBodyError
from libephemeris.hypothetical import _calc_uranian_planet_raw
from libephemeris.leb_composite import CompositeLEBReader
from libephemeris.leb_format import (
    BODY_PARAMS,
    COORD_HELIO_ECL,
    SECTION_BODY_INDEX,
    SECTION_COMPRESSED_CHEBYSHEV,
)
from libephemeris.leb_groups import LEB2_GROUPS
from libephemeris.leb_reader import open_leb
from libephemeris.time_utils import julday

from scripts.generate_leb import (
    COMPANION_BODY_GROUPS,
    assemble_leb,
    generate_body_helio,
)
from scripts.generate_leb2 import (
    LEB2_GROUPS as LEB2_GROUP_BODIES,
    _ALLOWED_GROUP_INVENTORIES,
    convert_leb1_to_leb2,
    resolve_group_source,
)

JD_START = julday(2023, 1, 1, 0.0)
JD_END = julday(2028, 1, 1, 0.0)
SAMPLE_JDS = (julday(2024, 3, 15, 7.25), julday(2026, 11, 2, 18.5))


@pytest.fixture(scope="module")
def uranians_leb1(tmp_path_factory):
    """Standalone uranians LEB1 partial over a small test range."""
    path = tmp_path_factory.mktemp("uranians") / "ephemeris_test_uranians.leb"
    assemble_leb(
        output=str(path),
        jd_start=JD_START,
        jd_end=JD_END,
        bodies=list(range(40, 48)),
        workers=1,
        verbose=False,
        # Real regeneration includes shared sections in the LEB1 partial; the
        # LEB2 companion converter must strip those duplicates.
        skip_aux=False,
    )
    return str(path)


@pytest.fixture(scope="module")
def uranians_leb2(uranians_leb1, tmp_path_factory):
    """LEB2 companion converted from the standalone partial."""
    path = tmp_path_factory.mktemp("uranians2") / "base_uranians.leb2"
    convert_leb1_to_leb2(
        input_path=uranians_leb1,
        output_path=str(path),
        group="uranians",
        verbose=False,
    )
    return str(path)


def _clear_trust_caches() -> None:
    state._FICTITIOUS_SOURCE_CACHE.clear()
    state._FICTITIOUS_TRUST_BY_READER.clear()


# ---------------------------------------------------------------------------
# Registry pins
# ---------------------------------------------------------------------------


def test_uranians_registered_consistently() -> None:
    """Format params, generator group, and LEB2 inventory stay in sync."""
    assert COMPANION_BODY_GROUPS["uranians"] == list(range(40, 48))
    assert LEB2_GROUPS[-1] == "uranians"
    assert LEB2_GROUP_BODIES["uranians"] == sorted(range(40, 48))
    assert _ALLOWED_GROUP_INVENTORIES["uranians"] == (frozenset(range(40, 48)),)
    for bid in range(40, 48):
        assert BODY_PARAMS[bid] == (3652.5, 8, COORD_HELIO_ECL, 3)


def test_uranians_compression_target_is_tight() -> None:
    """Angular channels keep the 1e-12 native-component target."""
    from libephemeris.leb_compression import BODY_TARGET_AU

    for bid in range(40, 48):
        assert BODY_TARGET_AU[bid] == 1e-12


def test_companion_groups_stay_out_of_main_generation() -> None:
    """Full generation and merged mains never include fictitious IDs."""
    from scripts.generate_leb import BODY_GROUPS

    assert "uranians" not in BODY_GROUPS
    for bodies in BODY_GROUPS.values():
        assert not any(40 <= bid <= 58 for bid in bodies)


def test_resolve_group_source_uses_sibling_partial() -> None:
    """convert-all reads companion-only groups from the standalone partial."""
    main = "data/leb/ephemeris_base.leb"
    assert resolve_group_source(main, "base", "core") == main
    assert (
        resolve_group_source(main, "base", "uranians")
        == "data/leb/ephemeris_base_uranians.leb"
    )


# ---------------------------------------------------------------------------
# Generation
# ---------------------------------------------------------------------------


def test_generate_body_helio_rejects_non_hamburg_ids() -> None:
    """Only IDs 40-47 have an independently sourced heliocentric model."""
    for bid in (39, 48, 50, 56, 58):
        with pytest.raises(ValueError, match="retired"):
            generate_body_helio(bid, JD_START, JD_END, 3652.5, 8)


def test_leb1_partial_matches_runtime_model(uranians_leb1) -> None:
    """Fitted coefficients reproduce the Neely propagation off the fit nodes."""
    reader = open_leb(uranians_leb1)
    try:
        for jd in SAMPLE_JDS:
            for bid in range(40, 48):
                (lon, lat, dist), _vel = reader.eval_body(bid, jd)
                ref_lon, ref_lat, ref_dist = _calc_uranian_planet_raw(bid, jd)
                assert lon == pytest.approx(ref_lon, abs=1e-9)
                assert lat == pytest.approx(ref_lat, abs=1e-9)
                assert dist == pytest.approx(ref_dist, abs=1e-10)
    finally:
        reader.close()


def test_leb2_conversion_keeps_precision_and_inventory(uranians_leb2) -> None:
    """The compressed companion stays within the 1e-12-target error class."""
    reader = open_leb(uranians_leb2)
    try:
        assert sorted(reader._bodies) == list(range(40, 48))
        for jd in SAMPLE_JDS:
            for bid in (40, 43, 47):
                (lon, lat, dist), _vel = reader.eval_body(bid, jd)
                ref_lon, ref_lat, ref_dist = _calc_uranian_planet_raw(bid, jd)
                assert lon == pytest.approx(ref_lon, abs=1e-9)
                assert lat == pytest.approx(ref_lat, abs=1e-9)
                assert dist == pytest.approx(ref_dist, abs=1e-10)
    finally:
        reader.close()


def test_leb2_companion_omits_shared_core_sections(uranians_leb2) -> None:
    """Named companions contain body channels, never duplicated core tables."""
    reader = open_leb(uranians_leb2)
    try:
        assert set(reader._sections) == {
            SECTION_BODY_INDEX,
            SECTION_COMPRESSED_CHEBYSHEV,
        }
    finally:
        reader.close()


def test_leb2_conversion_rejects_partial_inventory(tmp_path) -> None:
    """A source missing one Hamburg body cannot become a named-group artifact."""
    source = tmp_path / "ephemeris_test_uranians.leb"
    assemble_leb(
        output=str(source),
        jd_start=JD_START,
        jd_end=JD_END,
        bodies=list(range(40, 47)),  # Poseidon missing
        workers=1,
        verbose=False,
        skip_aux=True,
    )
    with pytest.raises(ValueError, match="invalid 'uranians' body inventory"):
        convert_leb1_to_leb2(
            input_path=str(source),
            output_path=str(tmp_path / "base_uranians.leb2"),
            group="uranians",
            verbose=False,
        )


# ---------------------------------------------------------------------------
# Trust gate
# ---------------------------------------------------------------------------


def test_companion_attaches_without_a_content_check(
    uranians_leb2, test_leb_file, tmp_path
) -> None:
    """A manifest-named companion attaches on presence alone.

    Runtime content verification was removed in 3.0.0rc15: integrity is bought
    once, when the installer writes the file (see
    ``state._matches_pinned_data_file``). An artifact this test never pinned
    therefore attaches and serves its fictitious channels — the deliberate
    inverse of the previous trust gate.
    """
    primary = tmp_path / "base_core.leb"
    shutil.copy(test_leb_file, primary)
    shutil.copy(uranians_leb2, tmp_path / "base_uranians.leb2")

    _clear_trust_caches()
    composite = CompositeLEBReader.from_file_with_companions(str(primary))
    try:
        assert composite.has_body(40)
        assert state.leb_fictitious_source_trusted(composite, 40)
    finally:
        composite.close()


def test_companion_channels_are_trusted_without_leaking_sideways(
    uranians_leb2, test_leb_file, tmp_path
) -> None:
    """The companion's own channels are trusted; the primary's are not."""
    primary = tmp_path / "base_core.leb"
    shutil.copy(test_leb_file, primary)
    companion = tmp_path / "base_uranians.leb2"
    shutil.copy(uranians_leb2, companion)

    _clear_trust_caches()
    composite = CompositeLEBReader.from_file_with_companions(str(primary))
    try:
        assert composite.has_body(40)
        for bid in range(40, 48):
            assert state.leb_fictitious_source_trusted(composite, bid)
        # Non-fictitious bodies resolve through the primary: trust does not
        # leak sideways from the pinned companion.
        assert not state.leb_fictitious_source_trusted(composite, 48)
    finally:
        composite.close()
        _clear_trust_caches()


def test_trust_no_longer_depends_on_content(uranians_leb2) -> None:
    """A direct reader over fictitious channels is trusted on presence alone."""
    _clear_trust_caches()
    reader = open_leb(uranians_leb2)
    try:
        assert state.leb_fictitious_source_trusted(reader, 40)
    finally:
        reader.close()
        _clear_trust_caches()


# ---------------------------------------------------------------------------
# Public dispatch
# ---------------------------------------------------------------------------


@pytest.fixture
def trusted_composite(uranians_leb2, test_leb_file, monkeypatch):
    """Activate a composite (aux bodies + uranians) as the LEB reader."""
    _clear_trust_caches()
    composite = CompositeLEBReader([open_leb(test_leb_file), open_leb(uranians_leb2)])
    with state._INIT_LOCK:
        state._LEB_READER = composite
        state._LEB_FILE = "injected-uranians-test"
    ephem.set_calc_mode("auto")
    yield composite
    with state._INIT_LOCK:
        state._LEB_READER = None
        state._LEB_FILE = None
    composite.close()
    _clear_trust_caches()


def test_filebacked_output_matches_runtime(trusted_composite) -> None:
    """File-backed sourcing reproduces the pure runtime public output."""
    flag_sets = (
        FLG_SPEED,
        FLG_HELCTR | FLG_SPEED,
        FLG_EQUATORIAL | FLG_SPEED,
        FLG_XYZ | FLG_SPEED,
    )
    for jd in SAMPLE_JDS:
        for bid in (40, 43, 47):
            for flags in flag_sets:
                filebacked, ret_fb = ephem.calc_ut(jd, bid, FLG_SWIEPH | flags)
                ephem.set_calc_mode("skyfield")
                runtime, ret_rt = ephem.calc_ut(jd, bid, FLG_SWIEPH | flags)
                ephem.set_calc_mode("auto")
                assert ret_fb == ret_rt
                for got, want in zip(filebacked, runtime):
                    assert got == pytest.approx(want, abs=1e-9)


def test_out_of_range_dates_fall_back_to_model(trusted_composite) -> None:
    """Dates outside the companion range still produce runtime-model output."""
    jd = julday(1900, 6, 1, 12.0)  # far outside the 2023-2028 test range
    filebacked, _ = ephem.calc_ut(jd, 42, FLG_SWIEPH | FLG_SPEED)
    ephem.set_calc_mode("skyfield")
    runtime, _ = ephem.calc_ut(jd, 42, FLG_SWIEPH | FLG_SPEED)
    for got, want in zip(filebacked, runtime):
        assert got == pytest.approx(want, abs=1e-9)


def test_unverified_fictitious_ids_still_fail_closed(trusted_composite) -> None:
    """Nibiru and Waldemath keep raising even with trusted uranians channels."""
    for bid in (49,):
        with pytest.raises(UnknownBodyError):
            ephem.calc_ut(SAMPLE_JDS[0], bid, FLG_SWIEPH | FLG_SPEED)


def test_fast_calc_guard_still_rejects_fictitious(trusted_composite) -> None:
    """The persisted-channel guard is untouched: fast_calc never serves 40-58."""
    from libephemeris import fast_calc

    with pytest.raises(KeyError, match="runtime analytical model"):
        fast_calc.fast_calc_ut(trusted_composite, SAMPLE_JDS[0], 40, FLG_SWIEPH)


def _activate(composite) -> None:
    with state._INIT_LOCK:
        state._LEB_READER = composite
        state._LEB_FILE = "injected-uranians-test"
    ephem.set_calc_mode("auto")


def _deactivate(composite) -> None:
    with state._INIT_LOCK:
        state._LEB_READER = None
        state._LEB_FILE = None
    composite.close()
    _clear_trust_caches()


def _count_body_evals(composite, body_id, monkeypatch):
    counter = {"n": 0}
    real = composite.eval_body

    def spy(bid, jd):
        if bid == body_id:
            counter["n"] += 1
        return real(bid, jd)

    monkeypatch.setattr(composite, "eval_body", spy)
    return counter


def test_reader_is_consulted_for_the_body(
    uranians_leb2, test_leb_file, monkeypatch
) -> None:
    """The body channel is read from the reader, not from the runtime model."""
    _clear_trust_caches()
    composite = CompositeLEBReader([open_leb(test_leb_file), open_leb(uranians_leb2)])
    calls = _count_body_evals(composite, 40, monkeypatch)
    _activate(composite)
    try:
        ephem.calc_ut(SAMPLE_JDS[0], 40, FLG_SWIEPH | FLG_SPEED)
        assert calls["n"] > 0, "trusted companion must be consulted for the body"
    finally:
        _deactivate(composite)


def test_reader_absent_from_the_manifest_still_serves_the_body(
    uranians_leb2, test_leb_file, monkeypatch
) -> None:
    """A file the manifest does not even list serves its fictitious channels.

    The pre-3.0.0rc15 gate refused to read 40-47 from an artifact whose bytes
    did not match a pin, falling back to the runtime model. That gate is gone
    by design: presence is now the whole requirement, and this test pins the
    new contract in the strongest available form — the manifest entry itself
    is deleted, and the channel is still consulted.
    """
    monkeypatch.delitem(download.DATA_FILES, "base_uranians.leb2", raising=False)
    _clear_trust_caches()
    composite = CompositeLEBReader([open_leb(test_leb_file), open_leb(uranians_leb2)])
    calls = _count_body_evals(composite, 40, monkeypatch)
    _activate(composite)
    try:
        filebacked, _ = ephem.calc_ut(SAMPLE_JDS[0], 40, FLG_SWIEPH | FLG_SPEED)
        assert calls["n"] > 0, "the channel must serve the body"
        assert all(math.isfinite(v) for v in filebacked)
    finally:
        _deactivate(composite)


def test_stencil_mixes_reader_and_model_at_range_edges(trusted_composite) -> None:
    """A speed stencil straddling the companion edge falls back per sample.

    At jd = JD_END - 0.5 the jd+1 sample is beyond the companion range, and at
    jd = JD_START + epsilon the light-time-retarded Earth eval steps below
    jd_start; both must degrade to the runtime model per sample (not raise),
    staying within the accepted LEB delta of the pure-model reference.
    """
    for jd in (JD_END - 0.5, JD_START + 1e-3):
        filebacked, _ = ephem.calc_ut(jd, 40, FLG_SWIEPH | FLG_SPEED)
        assert all(math.isfinite(v) for v in filebacked)
        ephem.set_calc_mode("skyfield")
        runtime, _ = ephem.calc_ut(jd, 40, FLG_SWIEPH | FLG_SPEED)
        ephem.set_calc_mode("auto")
        for got, want in zip(filebacked, runtime):
            assert got == pytest.approx(want, abs=1e-9)


def test_reviewed_core_attaches_every_same_prefix_sibling(
    uranians_leb2, tmp_path
) -> None:
    """Discovery: a reviewed core attaches all its same-prefix siblings.

    Exercises the production wiring (_has_pinned_sibling_companions +
    pinned_only=True) that the injected-reader fixtures bypass. Since runtime
    content verification was removed, ``pinned_only`` no longer narrows the
    set: every same-prefix group file present on disk attaches, including one
    whose bytes belong to a different group entirely.
    """
    bundled_core = os.path.join(
        os.path.dirname(state.__file__), "data", "leb2", "base_core.leb2"
    )
    core = tmp_path / "base_core.leb2"
    shutil.copy(bundled_core, core)
    uran = tmp_path / "base_uranians.leb2"
    shutil.copy(uranians_leb2, uran)
    # Same-prefix companion carrying the wrong group's bytes: it attaches too.
    shutil.copy(uranians_leb2, tmp_path / "base_asteroids.leb2")

    _clear_trust_caches()
    # Discovery is defined for LEB-permitting modes only: under a forced
    # skyfield/horizons backend get_leb_reader() returns None by contract,
    # so select 'auto' explicitly (as trusted_composite/_activate do; the
    # autouse reset_ephemeris_state fixture restores the previous mode).
    ephem.set_calc_mode("auto")
    ephem.set_leb_file(str(core))
    reader = state.get_leb_reader()
    try:
        assert reader.has_body(40), "uranians sibling must attach"
        assert state.leb_fictitious_source_trusted(reader, 40)
        # core + uranians + asteroids: presence is the only requirement now.
        assert len(reader._readers) == 3
    finally:
        ephem.set_leb_file(None)
        _clear_trust_caches()
