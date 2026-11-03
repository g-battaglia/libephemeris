"""Lossless offline coefficient serialization, without astronomical evaluation."""
import ast
import hashlib
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys

import pytest

from tools.generators import gen_iau2000_tables as generator


ROOT = Path(__file__).resolve().parents[2]
RUST = ROOT / "core/src/tables/iau2000_tables.rs"


def arrays(source):
    """Read literal Rust arrays, retaining declared types and dimensions."""
    result = []
    for name, declaration, literal in re.findall(
        r"pub const (IAU2000_\w+): (.*?) = (\[.*?\n\]);", source, re.S
    ):
        result.append((name, declaration, ast.literal_eval(literal.replace("_i32", ""))))
    return result


def words(values):
    result = b""
    for value in values:
        if isinstance(value, list):
            result += words(value)
        else:
            result += struct.pack("<d" if type(value) is float else "<q", value)
    return result


def test_all_declared_arrays_preserve_every_shipped_word_and_shape():
    raw = generator.DEFAULT_INPUT.read_bytes()
    assert len(raw) == 15728
    assert hashlib.sha256(raw).hexdigest() == generator.EXPECTED_SHA256
    actual = arrays(generator.generate(raw))
    expected = arrays(RUST.read_text())
    assert len(actual) == len(expected) == 10
    for (name, declaration, values), (old_name, old_declaration, old_values) in zip(
        actual, expected, strict=True
    ):
        assert (name, declaration) == (old_name, old_declaration)
        assert words(values) == words(old_values), name
        assert [len(row) if isinstance(row, list) else None for row in values] == [
            len(row) if isinstance(row, list) else None for row in old_values
        ]


def test_copied_generator_and_data_run_without_site_packages_git_or_checkout(tmp_path):
    script = tmp_path / "gen_iau2000_tables.py"
    shutil.copyfile(generator.__file__, script)
    data = tmp_path / "data"
    data.mkdir()
    shutil.copyfile(generator.DEFAULT_INPUT, data / generator.DEFAULT_INPUT.name)
    environment = {**os.environ, "PATH": "", "PYTHONDONTWRITEBYTECODE": "1"}
    first, second = tmp_path / "one.rs", tmp_path / "two.rs"
    for options, output in (([], first), (["--npz", str(data / generator.DEFAULT_INPUT.name)], second)):
        result = subprocess.run(
            [sys.executable, "-I", "-S", str(script), *options, "--out", str(output)],
            env=environment, cwd=tmp_path, capture_output=True, text=True, check=False,
        )
        assert result.returncode == 0, result.stderr
    assert first.read_bytes() == second.read_bytes() == generator.generate(
        generator.DEFAULT_INPUT.read_bytes()
    ).encode()


@pytest.mark.parametrize("payload", [b"", b"not an archive", b"PK\x03\x04", None])
@pytest.mark.parametrize("optimized", [False, True])
def test_wrong_identity_refuses_before_output_even_with_assertions_disabled(tmp_path, payload, optimized):
    source, output = tmp_path / "input.npz", tmp_path / "output.rs"
    source.write_bytes(payload if payload is not None else generator.DEFAULT_INPUT.read_bytes() + b"\n")
    output.write_text("keep previous coefficients")
    result = subprocess.run(
        [sys.executable, *(["-O"] if optimized else []), "-B", generator.__file__,
         "--npz", str(source), "--out", str(output)],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 2
    assert "input SHA-256" in result.stderr
    assert output.read_text() == "keep previous coefficients"
    assert sorted(path.name for path in tmp_path.iterdir()) == ["input.npz", "output.rs"]


def test_missing_input_and_output_parent_are_explicit_errors(tmp_path, capsys):
    assert generator.main(["--npz", str(tmp_path / "missing"), "--out", str(tmp_path / "out.rs")]) == 2
    assert generator.main(["--out", str(tmp_path / "missing/out.rs")]) == 2
    assert "IAU2000 generator:" in capsys.readouterr().err
    assert list(tmp_path.iterdir()) == []


def test_input_output_aliases_never_destroy_the_dataset(tmp_path):
    source = tmp_path / "input.npz"
    shutil.copyfile(generator.DEFAULT_INPUT, source)
    original = source.read_bytes()
    hardlink, symlink = tmp_path / "hardlink", tmp_path / "symlink"
    hardlink.hardlink_to(source)
    symlink.symlink_to(source)
    for output in (source, hardlink, symlink):
        assert generator.main(["--npz", str(source), "--out", str(output)]) == 2
        assert source.read_bytes() == original


def test_publication_failure_preserves_previous_output_and_removes_temporary(tmp_path, monkeypatch):
    output = tmp_path / "out.rs"
    output.write_text("prior table")
    def refuse(*args):
        raise OSError("controlled publication failure")
    monkeypatch.setattr(generator.os, "replace", refuse)
    assert generator.main(["--out", str(output)]) == 2
    assert output.read_text() == "prior table"
    assert list(tmp_path.iterdir()) == [output]
