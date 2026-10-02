# SPDX-License-Identifier: AGPL-3.0-only
from __future__ import annotations

import argparse
import sys

from .importer import upload_files


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="libephemeris-postgres")
    parser.add_argument("command", choices=["upload"])
    parser.add_argument("files", nargs="+")
    args = parser.parse_args(argv)
    try:
        upload_files(args.files)
        return 0
    except Exception:
        print("PostgreSQL provisioning failed", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
