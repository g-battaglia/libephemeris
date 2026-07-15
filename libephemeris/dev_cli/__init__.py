# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Dev CLI for libephemeris — development tools and task runner.

Provides the ``leph`` command with hierarchical subcommands for testing,
LEB generation, diagnostics, data downloads, release management, and more.

Shell completion + shell function (works with uv, no venv activation needed):

    eval "$(uv run ./leph completion zsh)"      # zsh — test it now
    uv run ./leph completion zsh >> ~/.zshrc     # zsh — make permanent
    uv run ./leph completion bash >> ~/.bashrc   # bash
    uv run ./leph completion fish > ~/.config/fish/conf.d/leph.fish  # fish
"""

from __future__ import annotations

import os
from pathlib import Path

import click

from .cmd_code import code_group
from .cmd_completion import completion_group
from .cmd_diag import diag_group
from .cmd_download import download_group
from .cmd_generate import generate_group
from .cmd_leb import leb_group
from .cmd_leb2 import leb2_group
from .cmd_manual import manual_group
from .cmd_release import release_group
from .cmd_status import status
from .cmd_test import test_group


@click.group(
    name="leph",
    help="""\
libephemeris developer CLI — over 100 commands for testing, data generation,
code quality, LEB binary ephemeris management, and release workflows.

\b
Quick start:
  leph test skyfield essential     # Fast sanity check (~900 tests, ~20s)
  leph test leb-backend unit-fast  # Recommended daily driver (~1 min)
  leph code lint                   # Ruff linter with auto-fix

\b
Subgroups at a glance:
  test        7 test suites: skyfield, leb-backend, lunar,
              leb-format, leb2-format, horizons, coverage
  code        Ruff linter/formatter, mypy type checker
  leb / leb2  Generate, verify, and compress LEB binary ephemeris files
  download    Fetch DE/SPK/IERS data, ASSIST files, and generation source kernels
  status      Show current data/configuration status for the local setup
  generate    Planet-center SPKs and Keplerian elements
  diag        Print body positions per tier, verify data integrity
  release     Retired LEB publication commands (fail closed)
  manual      Build user manuals (EPUB/PDF, Italian/English)
  completion  Generate shell completion scripts (zsh, bash, fish)

\b
TAB completion setup (works with uv, no venv activation needed):
  eval "$(uv run ./leph completion zsh)"      # try it now
  uv run ./leph completion zsh >> ~/.zshrc    # make permanent

Use -h on any subcommand for details:  leph test -h, leph test skyfield -h, etc.
""",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 120},
)
@click.version_option(package_name="libephemeris", prog_name="leph")
def cli() -> None:
    """Root CLI group."""


# Register all subgroups
cli.add_command(code_group)
cli.add_command(test_group)
cli.add_command(leb_group)
cli.add_command(leb2_group)
cli.add_command(download_group)
cli.add_command(diag_group)
cli.add_command(generate_group)
cli.add_command(release_group)
cli.add_command(manual_group)
cli.add_command(completion_group)
cli.add_command(status)


def main() -> None:
    """Run the repo-local ``leph`` developer CLI."""
    # Every command targets files in this checkout. ``uv run --project``
    # selects the environment but deliberately preserves the caller's cwd, so
    # anchor relative tests/scripts/data paths here for shell functions invoked
    # from another directory.
    os.chdir(Path(__file__).resolve().parents[2])
    # A stable program name is required for Click's ``_LEPH_COMPLETE``
    # protocol when the CLI is launched through ``python -m``.
    cli(prog_name="leph")
