# SPDX-License-Identifier: AGPL-3.0
"""Actual local-source discovery, declared tier policy and ambiguity refusals."""

import argparse
from dataclasses import replace
from pathlib import Path

import pytest

from tools.leb_generation.discovery import (
    LocalSpk,
    automatic_nbody_ids,
    scan_local_spks,
    select_local_spk,
    target_aliases,
)
from tools.leb_generation.options import resolve_bodies
from tools.leb_generation.read import read_body_index
from tools.leb_generation.source import check_input_hashes
from tools.leb_generation.verify import verify_source_file
from tools.tests.test_leb_generation_processes import source  # noqa: F401 -- shared native-data fixture


def candidate(name, low, high, *, seed=True, digest=None):
    return LocalSpk(Path(name), 20002060, ((low, high),), seed, digest or name)


def test_selection_uses_complete_native_coverage_instead_of_filename_or_size():
    narrow = candidate("apparently-wide-1600-2500.bsp", 1900.0, 2100.0)
    wide = candidate("unhelpful-name.bsp", 1600.0, 2500.0)
    assert select_local_spk(15, [narrow, wide], (1950.0, 2050.0)) == wide
    assert (
        select_local_spk(15, [narrow], (1850.0, 2150.0), available_domains=True)
        == narrow
    )
    with pytest.raises(ValueError, match="no local SPK"):
        select_local_spk(15, [narrow], (1850.0, 2150.0))
    with pytest.raises(ValueError, match="no local SPK"):
        select_local_spk(
            15,
            [replace(wide, windows_tt=((1600.0, 1900.0), (2000.0, 2500.0)))],
            (1800.0, 2200.0),
            available_domains=True,
        )


def test_ambiguous_solutions_and_incomparable_domains_require_explicit_source():
    first = candidate("first.bsp", 1800.0, 2200.0)
    second = candidate("second.bsp", 1800.0, 2200.0)
    with pytest.raises(ValueError, match="ambiguous.*--spk"):
        select_local_spk(15, [first, second], (1950.0, 2050.0))
    with pytest.raises(ValueError, match="ambiguous"):
        select_local_spk(
            15, [first, candidate("crossing.bsp", 1900.0, 2300.0)], (1950.0, 2050.0)
        )
    identical = replace(first, path=Path("copy.bsp"))
    assert select_local_spk(15, [first, identical], (1950.0, 2050.0)) == first


def test_nbody_seed_selection_requires_one_actual_heliocentric_j2000_realization():
    valid = candidate("seed.bsp", 1900.0, 2100.0)
    invalid = replace(
        valid, path=Path("not-a-seed.bsp"), seed_at_j2000=False, sha256="different"
    )
    assert select_local_spk(15, [invalid, valid], (1500.0, 2700.0), nbody=True) == valid
    with pytest.raises(ValueError, match="no local SPK.*J2000 seed"):
        select_local_spk(15, [invalid], (1500.0, 2700.0), nbody=True)
    with pytest.raises(ValueError, match="ambiguous"):
        select_local_spk(
            15,
            [valid, candidate("another-seed.bsp", 1600.0, 2500.0)],
            (1500.0, 2700.0),
            nbody=True,
        )


def test_native_aliases_and_extended_policy_preserve_requested_inventory():
    import libephemeris as engine
    from libephemeris.exotic_bodies import (
        EXOTIC_ASSIST_PERTURBER_IDS,
        EXOTIC_EXTENDED_IDS,
        EXOTIC_IDS,
    )

    assert target_aliases(15) == {2002060, 20002060}
    assert target_aliases(17) == {2000001, 20000001}
    with pytest.raises(ValueError, match="no native SPK identity"):
        target_aliases(0)
    args = argparse.Namespace(bodies=None, group="all", tier="extended")
    requested = resolve_bodies(args, engine)
    assert set(requested) & set(EXOTIC_IDS) == set(EXOTIC_EXTENDED_IDS)
    assert {0, 1, 14, 15, 17, 21, 22} <= set(requested)
    assert automatic_nbody_ids(requested, "extended") == set(EXOTIC_EXTENDED_IDS) - set(
        EXOTIC_ASSIST_PERTURBER_IDS
    )
    assert not automatic_nbody_ids(requested, "medium")
    # Explicit body requests are never silently removed by automatic assembly.
    args.bodies, args.group = [str(body) for body in EXOTIC_IDS], None
    assert set(resolve_bodies(args, engine)) == set(EXOTIC_IDS)


@pytest.fixture
def cache(tmp_path, source):  # noqa: F811 -- imported pytest fixture
    import shutil

    seed = source.kernel.parent / "spk/2060_190001_210001.bsp"
    if not seed.is_file():
        pytest.skip("automatic source discovery requires the existing Chiron image")
    check_input_hashes(
        [(seed, "9cb35103833f3437c89eb02f333595e3dffe0819c9db35dd390e946dc8f50cc7")]
    )
    root = tmp_path / "data"
    (root / "spk").mkdir(parents=True)
    (root / "de440.bsp").symlink_to(source.kernel)
    # A deliberately uninformative name also makes any accidental attempted
    # replacement test safe: this copy is owned by the test, never shared data.
    shutil.copyfile(seed, root / "spk/unnamed.bsp")
    return root


def test_native_scan_deduces_target_and_true_clock_domain(cache, source):  # noqa: F811 -- imported pytest fixture
    rows = scan_local_spks(source, [15, 17], [cache / "spk", cache / "spk"])
    assert rows[17] == [] and len(rows[15]) == 1
    row = rows[15][0]
    assert row.naif_id == 20002060 and row.seed_at_j2000
    assert (
        row.sha256 == "9cb35103833f3437c89eb02f333595e3dffe0819c9db35dd390e946dc8f50cc7"
    )
    assert row.windows_tt == (source.clock.tdb_interval_to_tt(2415020.5, 2488069.5),)


@pytest.mark.parametrize("compressed", [False, True])
def test_cli_automatically_generates_and_verifies_from_native_metadata(
    cache, tmp_path, source, compressed
):  # noqa: F811 -- imported pytest fixture
    from tools.leb_generation.cli_leb1 import main as leb1
    from tools.leb_generation.cli_leb2 import main as leb2

    output = tmp_path / ("auto.leb2" if compressed else "auto.leb")
    common = ["--tier", "medium", "--data-dir", str(cache), "--auto-sources"]
    fit = [
        "--start-jd",
        "2451545",
        "--end-jd",
        "2451549",
        "--bodies",
        "SUN",
        "CHIRON",
        "--skip-aux",
        "--workers",
        "2",
        "-o",
        str(output),
    ]
    assert (
        leb2(["generate", *common, *fit]) if compressed else leb1([*common, *fit])
    ) == 0
    assert [entry.body_id for entry in read_body_index(output)[1]] == [0, 15]
    source.add_spk_source(15, cache / "spk/unnamed.bsp", 20002060)
    assert verify_source_file(output, source, 17) == {
        "positions": 34,
        "cartesian_rates": 34,
    }
    assert (
        leb2(["verify", str(output), *common])
        if compressed
        else leb1(["--verify-only", str(output), *common])
    ) == 0


def test_automatic_selection_does_not_overwrite_a_selected_source(
    cache, source, capsys
):  # noqa: F811 -- imported pytest fixture
    from tools.leb_generation.cli_leb1 import main

    path = cache / "spk/unnamed.bsp"
    original = path.read_bytes()
    with pytest.raises(SystemExit) as result:
        main(
            [
                "--tier",
                "medium",
                "--data-dir",
                str(cache),
                "--auto-sources",
                "--bodies",
                "CHIRON",
                "--start-jd",
                "2451545",
                "--end-jd",
                "2451549",
                "--skip-aux",
                "--force",
                "-o",
                str(path),
            ]
        )
    assert result.value.code == 1
    assert "output must differ" in capsys.readouterr().err
    assert path.read_bytes() == original


def test_explicit_source_bypasses_ambiguous_discovery_for_that_body(
    cache, tmp_path, capsys
):
    from tools.leb_generation.cli_leb1 import main

    unrelated = cache / "spk/invalid.bsp"
    unrelated.write_bytes(
        b"invalid metadata must not be used for an already explicit body"
    )
    output = tmp_path / "explicit.leb"
    assert (
        main(
            [
                "--tier",
                "medium",
                "--data-dir",
                str(cache),
                "--auto-sources",
                "--bodies",
                "CHIRON",
                "--spk",
                "15",
                str(cache / "spk/unnamed.bsp"),
                "20002060",
                "--start-jd",
                "2451545",
                "--end-jd",
                "2451549",
                "--skip-aux",
                "-o",
                str(output),
            ]
        )
        == 0
    )
    assert output.is_file()
    assert "Minor-body source 15" in capsys.readouterr().out


def test_absent_body_or_corrupt_metadata_cannot_publish_a_partial_inventory(
    cache, tmp_path
):
    from tools.leb_generation.cli_leb1 import main

    output = tmp_path / "existing.leb"
    output.write_bytes(b"preserved")
    common = [
        "--tier",
        "medium",
        "--data-dir",
        str(cache),
        "--auto-sources",
        "--start-jd",
        "2451545",
        "--end-jd",
        "2451549",
        "--skip-aux",
        "--force",
        "-o",
        str(output),
    ]
    with pytest.raises(SystemExit) as result:
        main([*common, "--bodies", "SUN", "CERES"])
    assert result.value.code == 1 and output.read_bytes() == b"preserved"
    (cache / "spk/corrupt.bsp").write_bytes(b"DAF/SPK " + b"invalid SPK metadata")
    with pytest.raises(SystemExit) as result:
        main([*common, "--bodies", "SUN", "CHIRON"])
    assert result.value.code == 1 and output.read_bytes() == b"preserved"


def test_non_spk_force_or_cache_files_are_reported_without_becoming_sources(
    cache, tmp_path, capsys
):
    from tools.leb_generation.cli_leb1 import main

    force = cache / "spk/force.bsp"
    force.write_bytes(b"\x00\x05\x16\x07\x00\x02\x00\x00not a DAF SPK")
    output = tmp_path / "real-source.leb"
    assert (
        main(
            [
                "--tier",
                "medium",
                "--data-dir",
                str(cache),
                "--auto-sources",
                "--bodies",
                "CHIRON",
                "--start-jd",
                "2451545",
                "--end-jd",
                "2451549",
                "--skip-aux",
                "-o",
                str(output),
            ]
        )
        == 0
    )
    assert f"Ignored non-SPK cache file: {force}" in capsys.readouterr().out
    assert force.read_bytes().startswith(b"\x00\x05\x16\x07")


def test_automatic_nbody_uses_local_force_files_and_explicit_seed(
    cache, tmp_path, source
):  # noqa: F811 -- imported pytest fixture
    from tools.leb_generation.cli_leb1 import main

    forces = source.kernel.parent / "assist"
    identities = {
        forces
        / "linux_p1550p2650.440": "29915576d0a6555766b99485ac3056ee415e86df4fce282611c31afb329ad062",
        forces
        / "sb441-n16.bsp": "919d612ce3c72a78fc7158f9120156542d0f21e6b8b052e4c1339c759747fd90",
    }
    if not all(path.is_file() for path in identities):
        pytest.skip("automatic N-body assembly needs the existing force files")
    check_input_hashes(identities.items())
    (cache / "assist").mkdir()
    for path in identities:
        (cache / "assist" / path.name).symlink_to(path)
    output = tmp_path / "nbody.leb"
    common = [
        "--tier",
        "medium",
        "--data-dir",
        str(cache),
        "--auto-sources",
        "--nbody",
        "CHIRON",
    ]
    assert (
        main(
            [
                *common,
                "--start-jd",
                "2451545",
                "--end-jd",
                "2451549",
                "--bodies",
                "CHIRON",
                "--skip-aux",
                "--single",
                "-o",
                str(output),
            ]
        )
        == 0
    )
    source.add_spk_source(15, cache / "spk/unnamed.bsp", 20002060)
    source.add_nbody_sources([15], *identities.keys())
    assert verify_source_file(output, source, 17) == {
        "positions": 17,
        "cartesian_rates": 17,
    }
    assert main(["--verify-only", str(output), *common]) == 0


@pytest.mark.parametrize(
    "arguments", [["--spk-dir", "/unused"], ["--merge", "/unused", "--auto-sources"]]
)
def test_unused_automatic_source_selections_are_explicit_errors(arguments, capsys):
    from tools.leb_generation.cli_leb1 import main

    with pytest.raises(SystemExit) as result:
        main([*arguments, "--output", "/unused-output"])
    assert result.value.code == 1
    assert "--auto-sources" in capsys.readouterr().err
