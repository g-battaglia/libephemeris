"""The installed console entry executes the same native application as Rust."""
import sys

import pytest

from libephemeris._core import cli_main


def test_cli_version_returns_success(monkeypatch, capfd):
    monkeypatch.setattr(sys, "argv", ["libephemeris", "--version"])
    assert cli_main() == 0
    output = capfd.readouterr()
    assert "libephemeris" in output.out.lower()
    assert "4.0.0b4" in output.out
    assert output.err == ""


def test_cli_help_returns_success(monkeypatch, capfd):
    monkeypatch.setattr(sys, "argv", ["libephemeris", "--help"])
    assert cli_main() == 0
    output = capfd.readouterr()
    assert "leb" in output.out.lower()
    assert "download" in output.out.lower()
    assert output.err == ""


@pytest.mark.parametrize("arguments", [[], ["unknown-command"]])
def test_cli_usage_errors_return_two(monkeypatch, capfd, arguments):
    monkeypatch.setattr(sys, "argv", ["libephemeris", *arguments])
    assert cli_main() == 2
    output = capfd.readouterr()
    assert "libephemeris" in output.err.lower()
