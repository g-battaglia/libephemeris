# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Provision content-addressed PostgreSQL LEB2 storage."""

from __future__ import annotations

import argparse
import sys
from importlib.resources import files

import psycopg
from libephemeris import CoefficientSourceError

from .config import admin_dsn
from .importer import upload_files


def _schema(dsn: str | None) -> None:
    with psycopg.connect(admin_dsn(dsn)) as conn:
        conn.execute(files("libephemeris_postgres").joinpath("schema.sql").read_text())


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="libephemeris-postgres")
    sub = parser.add_subparsers(dest="command", required=True)
    for command in ("schema", "upload", "list"):
        child = sub.add_parser(command)
        child.add_argument("--dsn")
        if command == "upload":
            child.add_argument("files", nargs="+")
    args = parser.parse_args(argv)
    try:
        if args.command == "schema":
            _schema(args.dsn)
        elif args.command == "upload":
            upload_files(args.files, dsn=args.dsn)
        else:
            with psycopg.connect(admin_dsn(args.dsn)) as conn:
                for row in conn.execute(
                    "SELECT name,sha256,size,complete FROM libephemeris.files ORDER BY name"
                ):
                    print(*row)
        return 0
    except Exception as exc:
        print(
            str(exc)
            if isinstance(exc, CoefficientSourceError)
            else "PostgreSQL provisioning failed",
            file=sys.stderr,
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
