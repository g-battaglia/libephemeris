#!/usr/bin/env python3
"""
Download Swiss Ephemeris reference data files for comparison tests.

These files are used by pyswisseph as the oracle in compare_scripts/tests/.
Without them, pyswisseph falls back to Moshier analytical ephemeris which
does not support asteroids, barycentric coordinates, or fixed stars.

Files are downloaded to data/reference/ (gitignored).

Usage:
    python scripts/download_reference_data.py          # Download all
    python scripts/download_reference_data.py --force  # Re-download
"""

from __future__ import annotations

import argparse
import os
import sys
import urllib.request

BASE_URL = "https://raw.githubusercontent.com/aloistr/swisseph/master/ephe"

FILES = [
    ("sepl_18.se1", "Planets 1800-2400 CE"),
    ("semo_18.se1", "Moon 1800-2400 CE"),
    ("seas_18.se1", "Asteroids 1800-2400 CE"),
    ("sepl_12.se1", "Planets 1200-1800 CE (for DE440 start ~1550)"),
    ("semo_12.se1", "Moon 1200-1800 CE"),
    ("seas_12.se1", "Asteroids 1200-1800 CE"),
    ("sepl_24.se1", "Planets 2400-3000 CE (for DE440 end ~2650)"),
    ("semo_24.se1", "Moon 2400-3000 CE"),
    ("seas_24.se1", "Asteroids 2400-3000 CE"),
    ("sefstars.txt", "Fixed star catalog (True Citra, True Revati, etc.)"),
]

DEST_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "reference")


def download_file(filename: str, description: str, force: bool = False) -> bool:
    dest = os.path.join(DEST_DIR, filename)
    if os.path.exists(dest) and not force:
        print(f"  {filename}: already exists (skip)")
        return True

    url = f"{BASE_URL}/{filename}"
    print(f"  {filename}: downloading ({description})...", end="", flush=True)
    try:
        urllib.request.urlretrieve(url, dest)
        size_kb = os.path.getsize(dest) / 1024
        print(f" OK ({size_kb:.0f} KB)")
        return True
    except Exception as e:
        print(f" FAILED: {e}")
        if os.path.exists(dest):
            os.remove(dest)
        return False


def main() -> int:
    parser = argparse.ArgumentParser(description="Download reference data for comparison tests.")
    parser.add_argument("--force", action="store_true", help="Re-download even if files exist")
    args = parser.parse_args()

    os.makedirs(DEST_DIR, exist_ok=True)

    print(f"Downloading reference data to {DEST_DIR}/\n")

    ok = 0
    for filename, description in FILES:
        if download_file(filename, description, force=args.force):
            ok += 1

    print(f"\nDone: {ok}/{len(FILES)} files ready.")
    return 0 if ok == len(FILES) else 1


if __name__ == "__main__":
    sys.exit(main())
