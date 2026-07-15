#!/usr/bin/env python3
"""
Download SPK files for common celestial bodies.

Provenance:
    Project-authored selective downloader calling ``libephemeris.spk``, whose
    data source is NASA/JPL Horizons and whose format is the public NAIF SPK
    standard. Body groups, chunk sizes, retries, and cache paths are workflow
    choices; this script creates no astronomical coefficient or derived orbit.

This script pre-downloads SPK (SPICE kernel) files from JPL Horizons for common
minor bodies, allowing users to populate the cache ahead of time rather than
downloading on demand.

For full initialization (all bodies, all chunks, 1550-2650), use:
    libephemeris init

This script is useful for selective downloads (specific bodies, custom ranges).

Usage:
    python scripts/download_spk.py                    # Download all common bodies
    python scripts/download_spk.py --bodies chiron ceres  # Download specific bodies
    python scripts/download_spk.py --all              # Download all bodies
    python scripts/download_spk.py --list             # List available bodies
    python scripts/download_spk.py --cache-dir /path  # Use custom cache directory
    python scripts/download_spk.py --chunk 20         # Use 20-year chunks

"""

import argparse
import os
import sys
import time
from typing import Optional

# Ensure libephemeris can be imported from the project root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# Body definitions: name -> (horizons_id, ipl, naif_id)
# Populated lazily from SPK_BODY_NAME_MAP
AVAILABLE_BODIES: dict[str, tuple[str, int, int]] = {}


def _init_bodies() -> None:
    """Initialize body definitions from libephemeris constants.

    Reads all bodies from SPK_BODY_NAME_MAP in constants.py, ensuring
    the body list is always synchronized with the library.
    """
    global AVAILABLE_BODIES

    try:
        from libephemeris.constants import SPK_BODY_NAME_MAP
        from libephemeris.spk import _get_body_name

        AVAILABLE_BODIES = {}
        for ipl, (horizons_id, naif_id) in SPK_BODY_NAME_MAP.items():
            name = _get_body_name(ipl)
            if name is not None:
                AVAILABLE_BODIES[name.lower()] = (horizons_id, ipl, naif_id)
            else:
                # Fallback to a generic name based on ipl
                AVAILABLE_BODIES[f"body_{ipl}"] = (horizons_id, ipl, naif_id)

    except ImportError as e:
        print(f"Error importing libephemeris constants: {e}", file=sys.stderr)
        sys.exit(1)


# Common body groups for convenience
COMMON_BODIES = [
    "chiron",
    "pholus",
    "ceres",
    "pallas",
    "juno",
    "vesta",
    "eris",
    "sedna",
]


def download_spk_for_body(
    body_name: str,
    start_date: str,
    end_date: str,
    cache_dir: Optional[str] = None,
    force: bool = False,
    verbose: bool = True,
    chunk_years: Optional[int] = None,
) -> tuple[int, int, int]:
    """
    Download SPK file(s) for a specific body.

    Args:
        body_name: Name of the body (e.g., 'chiron', 'ceres')
        start_date: Start date in YYYY-MM-DD format
        end_date: End date in YYYY-MM-DD format
        cache_dir: Optional custom cache directory
        force: If True, re-download even if file exists
        verbose: If True, print progress messages
        chunk_years: If set, split the range into chunks of this many years.
            None means download as a single file.

    Returns:
        Tuple of (success_count, skipped_count, failed_count)
    """
    from libephemeris import spk_auto

    body_name_lower = body_name.lower()
    if body_name_lower not in AVAILABLE_BODIES:
        if verbose:
            print(f"Unknown body: {body_name}", file=sys.stderr)
        return (0, 0, 1)

    horizons_id, ipl, naif_id = AVAILABLE_BODIES[body_name_lower]

    # Determine cache directory
    if cache_dir is None:
        cache_dir = spk_auto.ensure_cache_dir()
    else:
        os.makedirs(cache_dir, exist_ok=True)

    success = 0
    skipped = 0
    failed = 0

    if chunk_years is not None:
        # Parse start/end years and generate chunks
        start_year = int(start_date.split("-")[0])
        end_year = int(end_date.split("-")[0])

        for chunk_start in range(start_year, end_year, chunk_years):
            chunk_end = min(chunk_start + chunk_years, end_year)
            chunk_start_date = f"{chunk_start}-01-01"
            chunk_end_date = f"{chunk_end}-01-01"

            s, sk, f = _download_single_chunk(
                body_name=body_name_lower,
                horizons_id=horizons_id,
                ipl=ipl,
                naif_id=naif_id,
                start_date=chunk_start_date,
                end_date=chunk_end_date,
                cache_dir=cache_dir,
                force=force,
                verbose=verbose,
            )
            success += s
            skipped += sk
            failed += f

            # Rate limiting between Horizons requests
            if s > 0:
                time.sleep(1.5)
    else:
        s, sk, f = _download_single_chunk(
            body_name=body_name_lower,
            horizons_id=horizons_id,
            ipl=ipl,
            naif_id=naif_id,
            start_date=start_date,
            end_date=end_date,
            cache_dir=cache_dir,
            force=force,
            verbose=verbose,
        )
        success += s
        skipped += sk
        failed += f

    return (success, skipped, failed)


def _download_single_chunk(
    body_name: str,
    horizons_id: str,
    ipl: int,
    naif_id: int,
    start_date: str,
    end_date: str,
    cache_dir: str,
    force: bool,
    verbose: bool,
) -> tuple[int, int, int]:
    """Download a single SPK chunk for a body."""
    from libephemeris import spk_auto
    from libephemeris.spk_auto import _iso_to_jd

    jd_start = _iso_to_jd(start_date)
    jd_end = _iso_to_jd(end_date)

    # Generate filename
    filename = spk_auto._generate_spk_cache_filename(horizons_id, jd_start, jd_end)
    output_path = os.path.join(cache_dir, filename)

    # Check if already cached
    if os.path.exists(output_path) and not force:
        if verbose:
            print(f"    {body_name} [{start_date} to {end_date}]: cached")
        return (0, 1, 0)

    # Download
    if verbose:
        print(
            f"    {body_name} [{start_date} to {end_date}]: downloading...",
            end="",
            flush=True,
        )

    try:
        spk_path = spk_auto.download_spk_from_horizons(
            body_id=horizons_id,
            jd_start=jd_start,
            jd_end=jd_end,
            output_path=output_path,
        )
        if verbose:
            size_mb = os.path.getsize(spk_path) / (1024 * 1024)
            print(f" OK ({size_mb:.2f} MB)")
        return (1, 0, 0)
    except Exception as e:
        if verbose:
            print(f" FAILED: {e}")
        return (0, 0, 1)


def list_available_bodies() -> None:
    """Print list of available bodies."""
    print("\nAvailable bodies for SPK download:")
    print("=" * 60)

    # Group by category using known classifications
    centaurs = ["chiron", "pholus", "nessus", "asbolus", "chariklo"]
    main_belt = ["ceres", "pallas", "juno", "vesta", "hygiea"]
    near_earth = ["eros", "apophis"]
    tnos = [
        "eris",
        "sedna",
        "haumea",
        "makemake",
        "ixion",
        "orcus",
        "quaoar",
        "varuna",
        "gonggong",
    ]

    def _print_group(title: str, names: list[str]) -> None:
        print(f"\n{title}:")
        for name in names:
            if name in AVAILABLE_BODIES:
                horizons_id, ipl, naif_id = AVAILABLE_BODIES[name]
                print(f"  {name:12s} - Horizons ID: {horizons_id:8s} (NAIF: {naif_id})")

    _print_group("Centaurs", centaurs)
    _print_group("Main Belt Asteroids", main_belt)
    _print_group("Near-Earth Asteroids", near_earth)
    _print_group("Trans-Neptunian Objects (TNOs)", tnos)

    # List any bodies not in the above categories
    all_categorized = set(centaurs + main_belt + near_earth + tnos)
    uncategorized = [n for n in sorted(AVAILABLE_BODIES) if n not in all_categorized]
    if uncategorized:
        _print_group("Other", uncategorized)

    print("\nPreset groups:")
    print("  --common     Chiron, Pholus, Ceres, Pallas, Juno, Vesta, Eris, Sedna")
    print("  --all        All available bodies")


def list_cache_contents(cache_dir: Optional[str] = None) -> None:
    """List contents of the SPK cache."""
    from libephemeris import spk_auto

    cached = spk_auto.list_cached_spk(cache_dir)

    if not cached:
        print("\nSPK cache is empty.")
        return

    print(f"\nCached SPK files ({len(cached)} files):")
    print("=" * 80)

    total_size = 0.0
    for entry in cached:
        filename = entry["filename"]
        size_mb = float(entry["size_mb"])
        date_start = entry.get("date_start", "Unknown")
        date_end = entry.get("date_end", "Unknown")
        total_size += size_mb

        print(f"  {filename}")
        print(f"    Size: {size_mb:.2f} MB | Coverage: {date_start} to {date_end}")

    print(f"\nTotal cache size: {total_size:.2f} MB")


def main() -> int:
    """Main entry point."""
    # Initialize body definitions
    _init_bodies()

    parser = argparse.ArgumentParser(
        description="Download SPK files for common celestial bodies.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  python scripts/download_spk.py                        # Download common bodies
  python scripts/download_spk.py --all                  # Download all bodies
  python scripts/download_spk.py --bodies chiron ceres  # Download specific bodies
  python scripts/download_spk.py --list                 # List available bodies
  python scripts/download_spk.py --list-cache           # Show cache contents
  python scripts/download_spk.py --start 2020-01-01 --end 2050-01-01
  python scripts/download_spk.py --chunk 20             # 20-year chunks
  python scripts/download_spk.py --all --chunk 20 --start 1550-01-01 --end 2650-01-01

For full initialization, prefer:  libephemeris init
        """,
    )

    parser.add_argument(
        "--bodies",
        nargs="+",
        help="Specific bodies to download (e.g., chiron ceres eris)",
    )
    parser.add_argument(
        "--common",
        action="store_true",
        help="Download common bodies (Chiron, Pholus, Ceres, Pallas, Juno, Vesta, Eris, Sedna)",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        dest="download_all",
        help="Download all available bodies",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        dest="list_bodies",
        help="List available bodies and exit",
    )
    parser.add_argument(
        "--list-cache",
        action="store_true",
        help="List cached SPK files and exit",
    )
    parser.add_argument(
        "--start",
        default="2000-01-01",
        help="Start date for SPK coverage (default: 2000-01-01)",
    )
    parser.add_argument(
        "--end",
        default="2100-01-01",
        help="End date for SPK coverage (default: 2100-01-01)",
    )
    parser.add_argument(
        "--chunk",
        type=int,
        default=None,
        metavar="YEARS",
        help="Split downloads into chunks of N years (e.g., 20)",
    )
    parser.add_argument(
        "--cache-dir",
        help="Custom cache directory (default: ~/.libephemeris/spk/)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download even if file already exists",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress progress output",
    )

    args = parser.parse_args()

    # Handle list options
    if args.list_bodies:
        list_available_bodies()
        return 0

    if args.list_cache:
        list_cache_contents(args.cache_dir)
        return 0

    # Determine which bodies to download
    bodies_to_download: list[str] = []

    if args.bodies:
        bodies_to_download = [b.lower() for b in args.bodies]
        # Validate bodies
        for body in bodies_to_download:
            if body not in AVAILABLE_BODIES:
                print(f"Error: Unknown body '{body}'", file=sys.stderr)
                print("Use --list to see available bodies", file=sys.stderr)
                return 1
    elif args.download_all:
        bodies_to_download = sorted(AVAILABLE_BODIES.keys())
    else:
        # Default to common bodies
        bodies_to_download = [b for b in COMMON_BODIES if b in AVAILABLE_BODIES]

    if not bodies_to_download:
        print("No bodies to download.", file=sys.stderr)
        return 1

    verbose = not args.quiet

    if verbose:
        print(f"\nDownloading SPK files for {len(bodies_to_download)} bodies...")
        print(f"Date range: {args.start} to {args.end}")
        if args.chunk:
            print(f"Chunk size: {args.chunk} years")
        if args.cache_dir:
            print(f"Cache directory: {args.cache_dir}")
        print()

    # Download each body
    total_success = 0
    total_skipped = 0
    total_failed = 0

    for body_name in bodies_to_download:
        if verbose:
            print(f"  {body_name}:")

        success, skipped, failed = download_spk_for_body(
            body_name=body_name,
            start_date=args.start,
            end_date=args.end,
            cache_dir=args.cache_dir,
            force=args.force,
            verbose=verbose,
            chunk_years=args.chunk,
        )
        total_success += success
        total_skipped += skipped
        total_failed += failed

    if verbose:
        print()
        parts = []
        if total_success:
            parts.append(f"{total_success} downloaded")
        if total_skipped:
            parts.append(f"{total_skipped} cached")
        if total_failed:
            parts.append(f"{total_failed} failed")
        print(f"Completed: {', '.join(parts)}")

    return 0 if total_failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
