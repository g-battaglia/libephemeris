# SPDX-License-Identifier: AGPL-3.0
"""Actual native readers preserve distinct body domains and refuse source gaps."""

import math

import numpy as np
import pytest

from tools.leb_generation.generate import generate_leb1, generate_leb2
from tools.leb_generation.read import read_body_index
from tools.leb_generation.source import check_input_hashes, continuous_overlap
from tools.leb_generation.verify import verify_source_file
from tools.tests.test_leb_generation import native_reader
from tools.tests.test_leb_generation_processes import source  # noqa: F401 -- shared native-data fixture
from tools.tests.test_leb_source_generation import LinearSource


class PartialSource(LinearSource):
    """Two analytic bodies with different available intervals."""

    def coverage(self, body_id):
        return [(2451546.0, 2451548.0)] if body_id == 14 else [(2451545.0, 2451549.0)]

    def available_domain(self, body_id, start, end):
        return continuous_overlap(self.coverage(body_id), start, end)

    def require_domain(self, body_id, start, end):
        if not any(low <= start < end <= high for low, high in self.coverage(body_id)):
            raise ValueError("synthetic body domain exceeded")


@pytest.mark.parametrize("generate", [generate_leb1, generate_leb2])
@pytest.mark.parametrize("inventory", [[0, 14], [14]])
def test_distinct_body_domains_are_stored_and_enforced_by_native_reader(tmp_path, generate, inventory):
    import libephemeris as engine

    path = tmp_path / "partial.leb"
    checks = generate(PartialSource(), path, inventory, 2451545.0, 2451549.0,
                      skip_aux=True, available_domains=True)
    assert checks.positions == checks.cartesian_rates == 9 * len(inventory)
    header, entries = read_body_index(path)
    assert (header.jd_start, header.jd_end) == (2451545.0, 2451549.0)
    expected = {0: (2451545.0, 2451549.0), 14: (2451546.0, 2451548.0)}
    assert [(row.body_id, row.jd_start, row.jd_end) for row in entries] == [
        (body, *expected[body]) for body in inventory]
    reader = native_reader(path)
    try:
        for row in entries:
            for epoch in (row.jd_start, (row.jd_start + row.jd_end) / 2.0, row.jd_end):
                actual = np.array(reader.eval_body(row.body_id, epoch)).ravel()
                np.testing.assert_allclose(actual, PartialSource().states(row.body_id, [epoch])[0],
                                           rtol=1e-10, atol=0.0)
        for epoch in (math.nextafter(2451546.0, -math.inf), math.nextafter(2451548.0, math.inf)):
            with pytest.raises(engine.Error):
                reader.eval_body(14, epoch)
            if 0 in inventory:
                assert reader.eval_body(0, epoch)
    finally:
        reader.close()
    assert verify_source_file(path, PartialSource(), 17) == {
        "positions": 17 * len(inventory), "cartesian_rates": 17 * len(inventory)}


def test_missing_body_overlap_does_not_publish_a_smaller_inventory(tmp_path):
    class AbsentBodySource(PartialSource):
        def coverage(self, body_id):
            return [] if body_id == 14 else super().coverage(body_id)

    path = tmp_path / "existing.leb"
    path.write_bytes(b"preserved")
    with pytest.raises(ValueError, match="no continuous source overlap"):
        generate_leb1(AbsentBodySource(), path, [0, 14], 2451545.0, 2451549.0,
                      skip_aux=True, available_domains=True, overwrite=True)
    assert path.read_bytes() == b"preserved"
    assert list(tmp_path.iterdir()) == [path]


@pytest.mark.parametrize("available", [False, True])
def test_default_strict_domain_and_explicit_gapped_overlap_preserve_output(tmp_path, available):
    class GappedSource(PartialSource):
        def coverage(self, body_id):
            return [(2451545.0, 2451546.0), (2451548.0, 2451549.0)]

    path = tmp_path / "existing.leb"
    path.write_bytes(b"preserved")
    with pytest.raises(ValueError, match="source gap" if available else "domain exceeded"):
        generate_leb1(GappedSource(), path, [0], 2451545.0, 2451549.0,
                      skip_aux=True, available_domains=available, overwrite=True)
    assert path.read_bytes() == b"preserved"
    assert list(tmp_path.iterdir()) == [path]


def test_continuous_overlap_keeps_endpoints_and_never_selects_only_one_island():
    assert continuous_overlap([(1.0, 3.0), (3.0, 7.0)], 0.0, 8.0) == (1.0, 7.0)
    assert continuous_overlap([(-math.inf, math.inf)], 1.0, 7.0) == (1.0, 7.0)
    with pytest.raises(ValueError, match="no continuous source overlap"):
        continuous_overlap([(1.0, 3.0)], 3.0, 5.0)
    with pytest.raises(ValueError, match="source gap"):
        continuous_overlap([(1.0, 3.0), (4.0, 7.0)], 2.0, 6.0)
    for start, end in ((math.nan, 3.0), (1.0, math.inf), (3.0, 3.0), (4.0, 1.0)):
        with pytest.raises(ValueError, match="finite with positive width"):
            continuous_overlap([(1.0, 7.0)], start, end)




@pytest.mark.parametrize("workers,single", [(1, False), (1, True), (2, False)])
def test_real_main_and_minor_spk_endpoints_keep_distinct_domains(tmp_path, source, workers, single):  # noqa: F811 -- imported pytest fixture
    seed = source.kernel.parent / "spk/2060_190001_210001.bsp"
    if not seed.is_file():
        pytest.skip("per-body generation needs the existing explicit Chiron image")
    check_input_hashes([(seed, "9cb35103833f3437c89eb02f333595e3dffe0819c9db35dd390e946dc8f50cc7")])
    source.add_spk_source(15, seed, 20002060)
    low, high = source.coverage(15)[0]
    start, end = low - 2.0, low + 2.0
    assert end < high
    path = tmp_path / "native.leb2"
    generate_leb2(source, path, [0, 15], start, end, skip_aux=True,
                  available_domains=True, workers=workers, single=single)
    header, entries = read_body_index(path)
    assert (header.jd_start, header.jd_end) == (start, end)
    assert [(row.body_id, row.jd_start, row.jd_end) for row in entries] == [
        (0, start, end), (15, low, end)]
    assert verify_source_file(path, source, 17) == {"positions": 34, "cartesian_rates": 34}


def test_tt_source_states_match_public_geometric_icrs_at_the_same_tt(source):  # noqa: F811 -- imported pytest fixture
    import libephemeris as engine

    # The Moon magnifies a TT/TDB labeling mistake that absolute barycentric
    # fitting errors can conceal. TRUEPOS/ICRS isolate one exact geometric ray.
    epochs = [2415108.6333333333, 2451545.125, 2460000.75]
    moon = source.states(1, epochs)
    earth = source.states(14, epochs)
    flags = (engine.FLG_SWIEPH | engine.FLG_TRUEPOS | engine.FLG_ICRS
             | engine.FLG_J2000 | engine.FLG_EQUATORIAL | engine.FLG_XYZ | engine.FLG_SPEED)
    expected = np.array([engine.calc(epoch, engine.MOON, flags)[0] for epoch in epochs])
    # Independent public reduction and the raw TT source must describe the
    # same ICRS subtraction, with only binary64 composition noise remaining.
    np.testing.assert_allclose(moon[:, :3] - earth[:, :3], expected[:, :3], rtol=0.0, atol=1e-14)
    np.testing.assert_allclose(moon[:, 3:] - earth[:, 3:], expected[:, 3:], rtol=0.0, atol=1e-14)


@pytest.mark.parametrize("mode", [["--merge", "/unused"], ["--verify-only", "/unused"]])
def test_available_domain_option_is_not_silently_ignored_by_other_modes(mode, capsys):
    from tools.leb_generation.cli_leb1 import main

    with pytest.raises(SystemExit) as result:
        main([*mode, "--available-domains", "--output", "/unused-output"])
    assert result.value.code == 1
    assert "applies only to generation" in capsys.readouterr().err
