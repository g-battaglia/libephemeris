# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Parseable retirement stubs for legacy LEB1 release commands."""

from __future__ import annotations

import click


RETIREMENT_MESSAGE = "LEB1 release upload is retired pending clean-room regeneration."


def _retired(version: str) -> None:
    """Reject a legacy release request without inspecting local assets."""
    del version
    raise click.ClickException(RETIREMENT_MESSAGE)


@click.group(
    "release",
    short_help="Retired LEB1 release commands; uploads are disabled.",
    help=(
        "Legacy LEB1 release commands remain parseable so automation fails "
        "safely. Upload is retired pending clean-room regeneration and review."
    ),
)
def release_group() -> None:
    """Retired release commands."""


def _retired_command(name: str) -> click.Command:
    """Create a legacy release subcommand that always fails closed."""

    @click.command(name, short_help="Retired pending clean-room regeneration.")
    @click.argument("version")
    def command(version: str) -> None:
        """Reject the legacy release request."""
        _retired(version)

    return command


for _command_name in (
    "leb",
    "leb-base",
    "leb-medium",
    "leb-extended",
    "leb-dry-run",
):
    release_group.add_command(_retired_command(_command_name))
