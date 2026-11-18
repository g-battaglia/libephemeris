"""The installed console entry executes the same native application as Rust."""
import sys

import pytest

from libephemeris._core import cli_main


def test_cli_version_returns_success(monkeypatch, capfd):
    monkeypatch.setattr(sys, "argv", ["libephemeris", "--version"])
    assert cli_main() == 0
    output = capfd.readouterr()
    # Exact output rejects a stale beta extension with the same version prefix.
    assert output.out.strip() == "libephemeris, version 4.1.0"
    assert output.err == ""


def test_cli_help_returns_success(monkeypatch, capfd):
    monkeypatch.setattr(sys, "argv", ["libephemeris", "--help"])
    assert cli_main() == 0
    output = capfd.readouterr()
    assert "leb" in output.out.lower()
    assert "download" in output.out.lower()
    assert output.err == ""


def test_verbose_status_reports_opt_in_cache_budget(monkeypatch, capfd, tmp_path):
    monkeypatch.setenv("LIBEPHEMERIS_LEB_CACHE_MB", "32")
    monkeypatch.setattr(sys, "argv", ["libephemeris", "status", "-vv", "--data-dir", str(tmp_path)])
    assert cli_main() == 0
    assert "LEB process cache: 0 of 33554432 coefficient bytes" in capfd.readouterr().out


@pytest.mark.parametrize("arguments", [[], ["unknown-command"]])
def test_cli_usage_errors_return_two(monkeypatch, capfd, arguments):
    monkeypatch.setattr(sys, "argv", ["libephemeris", *arguments])
    assert cli_main() == 2
    output = capfd.readouterr()
    assert "libephemeris" in output.err.lower()
