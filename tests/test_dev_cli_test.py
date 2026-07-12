"""Regression tests for backend selection in the developer test CLI."""

from __future__ import annotations

from unittest.mock import patch

from libephemeris.dev_cli import cmd_test


def test_skyfield_runner_forces_skyfield_backend() -> None:
    """The command group must not inherit an auto-selected LEB backend."""
    with patch.object(cmd_test, "_pytest") as run_pytest:
        cmd_test._skyfield_pytest(["tests/example.py", "-q"])

    run_pytest.assert_called_once_with(
        ["tests/example.py", "-q", "--calc-mode", "skyfield"]
    )
