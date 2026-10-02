# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Command line interface for PostgreSQL coefficient datasets."""

from __future__ import annotations

import argparse
import sys

from .config import admin_dsn
from .importer import import_files
from .verify import verify_files


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m libephemeris_postgres")
    sub = parser.add_subparsers(dest="command", required=True)
    schema = sub.add_parser("schema")
    schema.add_argument("--dsn")
    imp = sub.add_parser("import")
    imp.add_argument("--tier", required=True, choices=("base", "medium", "extended"))
    imp.add_argument("--dsn")
    imp.add_argument("--dataset")
    imp.add_argument("--resume")
    imp.add_argument("files", nargs="+")
    verify = sub.add_parser("verify")
    verify.add_argument("--dataset", required=True)
    verify.add_argument("--dsn")
    verify.add_argument("files", nargs="+")
    listing = sub.add_parser("list")
    listing.add_argument("--dsn")
    return parser


def _schema(dsn: str | None) -> None:
    import psycopg
    from importlib.resources import files

    sql = files("libephemeris_postgres").joinpath("schema.sql").read_text()
    try:
        with psycopg.connect(admin_dsn(dsn), autocommit=True) as conn:
            with conn.cursor() as cur:
                cur.execute(sql)
    except Exception:
        raise RuntimeError("PostgreSQL provisioning failed") from None


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "schema":
            _schema(args.dsn)
        elif args.command == "import":
            dataset = import_files(
                args.files,
                tier=args.tier,
                dsn=args.dsn,
                dataset_id=args.dataset,
                resume=args.resume,
            )
            print(dataset)
        elif args.command == "verify":
            verify_files(args.files, dataset_id=args.dataset, dsn=args.dsn)
        elif args.command == "list":
            import psycopg
            from . import _artifact_reviewed

            with psycopg.connect(admin_dsn(args.dsn)) as conn:
                for row in conn.execute(
                    "SELECT d.dataset_id, d.tier, d.complete, a.name, a.sha256 "
                    "FROM libephemeris.datasets d LEFT JOIN libephemeris.artifacts a "
                    "USING(dataset_id) ORDER BY d.created_at, a.artifact_no"
                ):
                    print(*row, "reviewed=" + str(_artifact_reviewed(row[3], row[4])))
        return 0
    except Exception as exc:
        from libephemeris import CoefficientSourceError

        # Only the explicit source-error contract is safe to expose. Driver
        # diagnostics (including wrapped connection errors) can quote secrets.
        message = (
            str(exc)
            if isinstance(exc, CoefficientSourceError)
            else "PostgreSQL provisioning failed"
        )
        print(message, file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
