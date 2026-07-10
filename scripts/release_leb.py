#!/usr/bin/env python3
"""
Upload LEB binary ephemeris files to GitHub Releases.

This script calculates SHA256 hashes, uploads LEB files to a GitHub release
using the gh CLI, and optionally updates download.py with the new hashes.

Requirements:
    - gh CLI installed and authenticated (gh auth login)

Usage:
    # Upload all available LEB files
    python scripts/release_leb.py --version 0.22.0

    # Upload specific tier
    python scripts/release_leb.py --version 0.22.0 --tier medium

    # Dry run (show what would be uploaded)
    python scripts/release_leb.py --version 0.22.0 --dry-run

    # Upload and auto-update download.py with new SHA256 hashes
    python scripts/release_leb.py --version 0.22.0 --update-hashes

Files are uploaded to:
    https://github.com/g-battaglia/libephemeris/releases/download/data-v1/
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import subprocess
import sys
from pathlib import Path


TIER_FILES = {
    "base": "ephemeris_base.leb",
    "medium": "ephemeris_medium.leb",
    "extended": "ephemeris_extended.leb",
}

DEFAULT_TIERS = ["base", "medium", "extended"]

RELEASE_TAG = "data-v1"
REPO = "g-battaglia/libephemeris"


def get_data_dir() -> Path:
    """Get the libephemeris data directory.

    Checks LIBEPHEMERIS_DATA_DIR env var, falls back to ~/.libephemeris.
    """
    env_dir = os.environ.get("LIBEPHEMERIS_DATA_DIR")
    if env_dir:
        return Path(env_dir)
    return Path.home() / ".libephemeris"


def find_leb_file(filename: str) -> Path | None:
    """Find a LEB file in common locations.

    Search order:
        1. data/leb/ in the repository (generation output directory)
        2. Data directory (~/.libephemeris/leb/)

    Args:
        filename: LEB filename to find.

    Returns:
        Path to the file, or None if not found.
    """
    # 1. Repository data/leb/ directory (where generate_leb.py writes output)
    repo_root = Path(__file__).parent.parent
    repo_path = repo_root / "data" / "leb" / filename
    if repo_path.exists():
        return repo_path

    # 2. User data directory (~/.libephemeris/leb/)
    data_path = get_data_dir() / "leb" / filename
    if data_path.exists():
        return data_path

    return None


def calculate_sha256(filepath: Path) -> str:
    """Calculate SHA256 hash of a file."""
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256_hash.update(chunk)
    return sha256_hash.hexdigest()


def get_file_size_mb(filepath: Path) -> float:
    """Get file size in MB."""
    return filepath.stat().st_size / (1024 * 1024)


def check_gh_auth() -> bool:
    """Check if gh CLI is authenticated."""
    try:
        result = subprocess.run(
            ["gh", "auth", "status"],
            capture_output=True,
            text=True,
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


def check_release_exists(tag: str) -> bool:
    """Check if a release with the given tag exists."""
    try:
        result = subprocess.run(
            ["gh", "release", "view", tag, "--repo", REPO],
            capture_output=True,
            text=True,
        )
        return result.returncode == 0
    except Exception:
        return False


def upload_file(filepath: Path, tag: str, dry_run: bool = False) -> bool:
    """Upload a file to a GitHub release."""
    if dry_run:
        print(f"  [DRY RUN] Would upload: {filepath.name}")
        return True

    try:
        result = subprocess.run(
            [
                "gh",
                "release",
                "upload",
                tag,
                str(filepath),
                "--repo",
                REPO,
                "--clobber",
            ],
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            print(f"  [OK] Uploaded: {filepath.name}")
            return True
        else:
            print(f"  [FAIL] Upload failed: {result.stderr}")
            return False
    except Exception as e:
        print(f"  [FAIL] Exception: {e}")
        return False


def update_download_py(
    updates: list[tuple[str, str, float]],
    dry_run: bool = False,
) -> bool:
    """Update SHA256 hashes and sizes in download.py.

    Args:
        updates: List of (filename, sha256, size_mb) tuples.
        dry_run: If True, show changes without writing.

    Returns:
        True if all updates succeeded.
    """
    download_py = Path(__file__).parent.parent / "libephemeris" / "download.py"
    if not download_py.exists():
        print(f"  [FAIL] {download_py} not found")
        return False

    content = download_py.read_text()
    original = content

    for filename, sha256, size_mb in updates:
        # Update sha256 value (handles both None and quoted string)
        sha256_pattern = re.compile(
            rf'("{filename}":\s*\{{[\s\S]*?"sha256":\s*)'
            rf'(None|"[a-f0-9]+")'
        )
        match = sha256_pattern.search(content)
        if match:
            content = (
                content[: match.start(2)] + f'"{sha256}"' + content[match.end(2) :]
            )
        else:
            print(f"  [WARN] Could not find sha256 entry for {filename}")

        # Update size_mb value (handles both None and float)
        size_pattern = re.compile(
            rf'("{filename}":\s*\{{[\s\S]*?"size_mb":\s*)'
            rf"(None|[\d.]+)"
        )
        match = size_pattern.search(content)
        if match:
            content = (
                content[: match.start(2)] + f"{size_mb:.1f}" + content[match.end(2) :]
            )
        else:
            print(f"  [WARN] Could not find size_mb entry for {filename}")

    if content == original:
        print("  [INFO] No changes needed in download.py")
        return True

    if dry_run:
        print("  [DRY RUN] Would update download.py with new hashes")
        return True

    download_py.write_text(content)
    print(f"  [OK] Updated {download_py}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Upload LEB binary ephemeris files to GitHub Releases",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
    python scripts/release_leb.py --version 0.22.0
    python scripts/release_leb.py --version 0.22.0 --tier medium
    python scripts/release_leb.py --version 0.22.0 --dry-run
    python scripts/release_leb.py --version 0.22.0 --update-hashes
        """,
    )

    parser.add_argument(
        "--version",
        required=True,
        help="Version string for logging purposes (e.g., 0.22.0)",
    )
    parser.add_argument(
        "--tier",
        choices=["base", "medium", "extended", "all"],
        default="all",
        help="Which tier(s) to upload (default: all)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be uploaded without actually uploading",
    )
    parser.add_argument(
        "--tag",
        default=RELEASE_TAG,
        help=f"Release tag (default: {RELEASE_TAG})",
    )
    parser.add_argument(
        "--update-hashes",
        action="store_true",
        help="Auto-update SHA256 hashes and sizes in libephemeris/download.py",
    )

    args = parser.parse_args()

    # Check authentication
    if not args.dry_run and not check_gh_auth():
        print("ERROR: gh CLI is not authenticated.")
        print("Run: gh auth login")
        return 1

    # Determine files to upload
    if args.tier == "all":
        tiers = DEFAULT_TIERS
    else:
        tiers = [args.tier]

    files_to_upload: list[tuple[str, Path]] = []
    for tier in tiers:
        filename = TIER_FILES[tier]
        filepath = find_leb_file(filename)
        if filepath is not None:
            files_to_upload.append((tier, filepath))
        else:
            search_paths = [
                Path(__file__).parent.parent / "data" / "leb",
                get_data_dir() / "leb",
            ]
            locations = ", ".join(str(p) for p in search_paths)
            print(f"[SKIP] {filename} not found in: {locations}")

    if not files_to_upload:
        print("ERROR: No LEB files found to upload.")
        print()
        print("Generate them first with:")
        print("  poe leb:generate:base")
        print("  poe leb:generate:medium")
        print("  poe leb:generate:extended")
        return 1

    # Compute hashes and sizes
    file_details: list[tuple[str, Path, str, float]] = []
    for tier, filepath in files_to_upload:
        sha256 = calculate_sha256(filepath)
        size_mb = get_file_size_mb(filepath)
        file_details.append((tier, filepath, sha256, size_mb))

    # Print summary
    print("=" * 70)
    print("LEB BINARY EPHEMERIS RELEASE UPLOAD")
    print("=" * 70)
    print(f"Version: {args.version}")
    print(f"Release tag: {args.tag}")
    print(f"Repository: {REPO}")
    print(f"Dry run: {args.dry_run}")
    print()
    print("Files to upload:")
    print()

    for tier, filepath, sha256, size_mb in file_details:
        print(f"  {filepath.name}")
        print(f"    Path:   {filepath}")
        print(f"    Size:   {size_mb:.1f} MB")
        print(f"    SHA256: {sha256}")
        print()

    # Check release exists
    if not args.dry_run:
        if not check_release_exists(args.tag):
            print(f"WARNING: Release '{args.tag}' does not exist.")
            print(f"Create it first with: gh release create {args.tag} --repo {REPO}")
            return 1

    # Upload files
    print("Uploading...")
    print()

    success_count = 0
    for tier, filepath, sha256, size_mb in file_details:
        if upload_file(filepath, args.tag, args.dry_run):
            success_count += 1

    # Summary
    print()
    print("=" * 70)
    if args.dry_run:
        print(
            f"DRY RUN COMPLETE: {success_count}/{len(file_details)} files would be uploaded"
        )
    else:
        print(f"UPLOAD COMPLETE: {success_count}/{len(file_details)} files uploaded")
    print("=" * 70)

    # Print URLs
    if not args.dry_run:
        print()
        print("Download URLs:")
        base_url = f"https://github.com/{REPO}/releases/download/{args.tag}"
        for tier, filepath, sha256, size_mb in file_details:
            print(f"  {base_url}/{filepath.name}")

    # Update download.py
    if args.update_hashes:
        print()
        print("Updating download.py...")
        hash_updates = [
            (filepath.name, sha256, size_mb)
            for tier, filepath, sha256, size_mb in file_details
        ]
        update_download_py(hash_updates, dry_run=args.dry_run)
    else:
        # Always print the manual instructions
        print()
        print("To update download.py with SHA256 hashes, either:")
        print("  1. Re-run with --update-hashes to auto-update")
        print("  2. Manually update DATA_FILES in libephemeris/download.py:")
        print()
        for tier, filepath, sha256, size_mb in file_details:
            print(
                f'    "{filepath.name}": {{"sha256": "{sha256}", "size_mb": {size_mb:.1f}, ...}},'
            )

    return 0


if __name__ == "__main__":
    sys.exit(main())
