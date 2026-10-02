# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Explicit provisioning CLI: python -m libephemeris.db schema|import.

Provenance:
    Project-authored command orchestration. No provisioning occurs implicitly
    during library import or runtime calculations. This module defines no
    astronomical model or coefficient.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from .importer import connect_provisioner, import_artifacts, provision_schema


def main() -> None:
    """Parse an explicit offline provisioning action and execute it.

    Raises:
        SystemExit: Arguments are invalid or the provisioning action fails.
    """
    parser = argparse.ArgumentParser(description="Provision PostgreSQL ephemeris data")
    parser.add_argument("--dsn", default=os.environ.get("LIBEPHEMERIS_DB_URL"))
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("schema", help="Explicitly install the versioned schema")
    importer = commands.add_parser("import", help="Publish an immutable LEB dataset")
    importer.add_argument("--dataset", required=True)
    importer.add_argument(
        "--tier", choices=("base", "medium", "extended"), required=True
    )
    importer.add_argument("paths", nargs="+", type=Path)
    arguments = parser.parse_args()
    if not arguments.dsn:
        parser.error("Set LIBEPHEMERIS_DB_URL or supply --dsn")

    # Provisioning is an explicit network action, unlike a calculation.
    from ..net import get_configured_network_policy, set_network_policy

    previous_policy = get_configured_network_policy()
    set_network_policy("allow")
    try:
        with connect_provisioner(arguments.dsn) as connection:
            if arguments.command == "schema":
                with connection.transaction():
                    provision_schema(connection)
            else:
                import_artifacts(
                    connection,
                    arguments.dataset,
                    arguments.paths,
                    arguments.tier,
                )
    except Exception as error:
        # Only library/input errors are safe to display. Driver diagnostics
        # can quote connection secrets or query parameters.
        from ..exceptions import Error

        message = (
            str(error)
            if isinstance(error, (Error, ValueError))
            else "PostgreSQL provisioning failed"
        )
        parser.exit(1, message + "\n")
    finally:
        set_network_policy(previous_policy)


if __name__ == "__main__":
    main()
