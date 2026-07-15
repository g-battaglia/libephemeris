#!/usr/bin/env python3
"""Retired entry point for legacy LEB1 release uploads.

The command-line shape is retained so existing automation fails safely and
clearly. Publishing prebuilt LEB1 assets remains disabled until those files
have been regenerated from independently sourced inputs and reviewed under the
current clean-room policy.

Provenance:
    Project-authored fail-closed compatibility entry point. It performs no
    upload, generation, or astronomical calculation. The retained command name
    is workflow metadata only; any future release implementation must use the
    registered generators, source attestations, hashes, and release review.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence


RETIREMENT_MESSAGE = "LEB1 release upload is retired pending clean-room regeneration."


def _build_parser() -> argparse.ArgumentParser:
    """Build the legacy parser without exposing any upload implementation."""
    parser = argparse.ArgumentParser(
        description=RETIREMENT_MESSAGE,
        epilog=(
            "Generate files locally for development with `leph leb generate`, "
            "but do not publish them until their provenance review is complete."
        ),
    )
    parser.add_argument(
        "--version",
        required=True,
        help="Version string retained for command-line compatibility.",
    )
    parser.add_argument(
        "--tier",
        choices=["base", "medium", "extended", "all"],
        default="all",
        help="Legacy tier selector; publishing remains disabled.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Parse the legacy option; no assets are inspected or uploaded.",
    )
    parser.add_argument(
        "--tag",
        default="data-v1",
        help="Legacy release tag; publishing remains disabled.",
    )
    parser.add_argument(
        "--update-hashes",
        action="store_true",
        help="Legacy option; the retired command never changes download metadata.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Parse the legacy command and reject every release attempt."""
    _build_parser().parse_args(argv)
    print(f"ERROR: {RETIREMENT_MESSAGE}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
