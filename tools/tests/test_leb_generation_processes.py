# SPDX-License-Identifier: AGPL-3.0
"""Real spawned native fitting, coefficient preservation and failed publication."""

import os
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from tools.leb_generation.generate import generate_leb1, generate_leb2
from tools.leb_generation.read import read_leb1
from tools.leb_generation.source import NativeSource, check_input_hashes
from tools.leb_generation.verify import verify_source_file


@pytest.fixture
def source():
    directory = Path(os.environ.get("LIBEPHEMERIS_TEST_DATA_DIR", Path.home() / ".libephemeris"))
    path = directory / "de440.bsp"
    if not path.is_file():
        pytest.skip("isolated native fitting requires the existing local DE440 image")
    check_input_hashes([(path, "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5")])
    return NativeSource(path, "medium")


@pytest.mark.parametrize("workers,single", [(2, False), (1, True)])
def test_spawned_body_fits_preserve_coefficients_and_requested_order(tmp_path, source, workers, single):
    original, isolated = tmp_path / "original.leb", tmp_path / "isolated.leb"
    bodies = [14, 0, 3, 21, 22]  # Mixed Cartesian/angular domains, more tasks than slots.
    expected = generate_leb1(source, original, bodies, 2451545.0, 2451549.0, skip_aux=True)
    actual = generate_leb1(source, isolated, bodies, 2451545.0, 2451549.0,
                           skip_aux=True, workers=workers, single=single)
    assert actual == expected
    parts = [read_leb1(path)[1] for path in (original, isolated)]
    for rows in parts:
        assert [body.entry.body_id for body in rows] == bodies
    for left, right in zip(*parts, strict=True):
        assert left.entry == right.entry
        for (li, lv), (ri, rv) in zip(left.chunks(1), right.chunks(1), strict=True):
            assert li == ri
            np.testing.assert_array_equal(lv, rv)
    assert verify_source_file(isolated, source, 17) == {"positions": 85, "cartesian_rates": 51}


@pytest.mark.parametrize("nbody", [False, True])
def test_isolated_minor_sources_retain_explicit_seed_and_force_images(tmp_path, source, nbody):
    directory = source.kernel.parent
    seed = directory / "spk/2060_190001_210001.bsp"
    planets = directory / "assist/linux_p1550p2650.440"
    asteroids = directory / "assist/sb441-n16.bsp"
    identities = {seed: "9cb35103833f3437c89eb02f333595e3dffe0819c9db35dd390e946dc8f50cc7"}
    if nbody:
        identities.update({planets: "29915576d0a6555766b99485ac3056ee415e86df4fce282611c31afb329ad062",
                           asteroids: "919d612ce3c72a78fc7158f9120156542d0f21e6b8b052e4c1339c759747fd90"})
    if not all(path.is_file() for path in identities):
        pytest.skip("isolated minor fitting requires the explicit local seed/force images")
    check_input_hashes(identities.items())
    source.add_spk_source(15, seed, 20002060)
    if nbody:
        source.add_nbody_sources([15], planets, asteroids)
    selection = source.isolated_selection(15)
    assert dict(selection.hashes) == {source.kernel: source.sha256, **identities}
    original, isolated = tmp_path / "original.leb", tmp_path / "isolated.leb"
    expected = generate_leb1(source, original, [15], 2451545.0, 2451549.0, skip_aux=True)
    actual = generate_leb1(source, isolated, [15], 2451545.0, 2451549.0, skip_aux=True, single=True)
    assert actual == expected
    left, right = (read_leb1(path)[1][0] for path in (original, isolated))
    np.testing.assert_array_equal(next(left.chunks(1))[1], next(right.chunks(1))[1])
    assert verify_source_file(isolated, source, 17) == {"positions": 17, "cartesian_rates": 17}


def test_spawned_fits_publish_final_source_checked_leb2(tmp_path, source):
    output = tmp_path / "isolated.leb2"
    counts = generate_leb2(source, output, [0, 14], 2451545.0, 2451549.0,
                           skip_aux=True, workers=2)
    assert counts.positions == counts.cartesian_rates == 18
    assert verify_source_file(output, source, 17) == {"positions": 34, "cartesian_rates": 34}


def test_changed_child_source_identity_preserves_existing_destination(tmp_path, source, monkeypatch):
    output = tmp_path / "preserved.leb"
    output.write_bytes(b"existing destination")
    selection = source.isolated_selection(0)
    mismatched = replace(selection, hashes=((source.kernel, "0" * 64),))
    monkeypatch.setattr(source, "isolated_selection", lambda body: mismatched)
    with pytest.raises(ValueError, match="source identity changed"):
        generate_leb1(source, output, [0], 2451545.0, 2451549.0,
                      skip_aux=True, single=True, overwrite=True)
    assert output.read_bytes() == b"existing destination"
    assert list(tmp_path.iterdir()) == [output]


@pytest.mark.parametrize("workers,single", [(0, False), (-1, False), (2, True)])
def test_invalid_process_limits_refuse_before_writing(tmp_path, workers, single):
    with pytest.raises(ValueError, match="workers"):
        generate_leb1(None, tmp_path / "unused.leb", [0], 2451545.0, 2451549.0,
                      workers=workers, single=single)
    assert list(tmp_path.iterdir()) == []
