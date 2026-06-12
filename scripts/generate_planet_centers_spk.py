#!/usr/bin/env python3
"""
Generate tier-specific planet_centers SPK files for libephemeris.

This script downloads satellite SPK files from JPL NAIF and extracts planet
center segments (NAIF IDs 599, 699, 799, 899, 999) to create compact SPK files
for each precision tier.

USAGE
-----
    # Generate for a specific tier
    python scripts/generate_planet_centers_spk.py --tier base
    python scripts/generate_planet_centers_spk.py --tier medium
    python scripts/generate_planet_centers_spk.py --tier extended

    # Generate all tiers
    python scripts/generate_planet_centers_spk.py --all

    # Force re-download of source files
    python scripts/generate_planet_centers_spk.py --tier medium --force

OUTPUT
------
    planet_centers_base.bsp      (~15-20 MB, 1850-2150)
    planet_centers_medium.bsp    (~40-50 MB, 1550-2650)
    planet_centers_extended.bsp  (~80-100 MB, -12000 to +17000 partial)

REQUIREMENTS
------------
    - spiceypy >= 6.0.0 (pip install spiceypy)
    - Internet connection (to download source SPK files)
    - Disk space: ~500 MB (base), ~4 GB (medium), ~6.5 GB (extended)

TIER DETAILS
------------
    base:     Uses old compact satellite SPKs (jup204, sat319, etc.)
              Coverage: 1850-2150 for all planets

    medium:   Uses standard satellite SPKs (jup365, sat441, etc.)
              Coverage: 1550-2650 (partial for Pluto: 1800-2200)

    extended: Uses "xl" extended-range satellite SPKs where available
              Coverage: -12000 to +17000 for Uranus/Neptune
                        -502 to +4500 for Saturn
                        1600-2200 for Jupiter (max available)
                        1800-2200 for Pluto (max available)

REFERENCES
----------
    - JPL NAIF Generic Kernels: https://naif.jpl.nasa.gov/naif/data_generic.html
    - SPK File Format: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html

Author: libephemeris team
"""

from __future__ import annotations

import argparse
import os
import ssl
import struct
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np


def _is_valid_bsp(filepath: str) -> bool:
    """Check if a BSP file can be opened and its segments enumerated by jplephem."""
    try:
        from jplephem.spk import SPK

        with SPK.open(filepath) as kernel:
            for segment in kernel.segments:
                _ = segment.center, segment.target
        return True
    except Exception:
        return False


# =============================================================================
# CONFIGURATION: TIER-SPECIFIC SOURCE FILES
# =============================================================================

NAIF_BASE = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites"

# Source SPK files for each tier
# Format: (url, center_id, target_id) or list of such tuples for multi-file merge
TIER_SOURCES: Dict[
    str, Dict[str, Union[Tuple[str, int, int], List[Tuple[str, int, int]]]]
] = {
    "base": {
        # Old compact versions - good for 1850-2150 range
        "jupiter": (f"{NAIF_BASE}/a_old_versions/jup204.bsp", 5, 599),
        "saturn": (f"{NAIF_BASE}/a_old_versions/sat319.bsp", 6, 699),
        "uranus": (f"{NAIF_BASE}/a_old_versions/ura083.bsp", 7, 799),
        "neptune": (f"{NAIF_BASE}/a_old_versions/nep050.bsp", 8, 899),
        "pluto": (f"{NAIF_BASE}/a_old_versions/plu017.bsp", 9, 999),
    },
    "medium": {
        # Standard versions - good for 1550-2650 range
        "jupiter": (f"{NAIF_BASE}/jup365.bsp", 5, 599),  # 1600-2200
        "saturn": (f"{NAIF_BASE}/sat441.bsp", 6, 699),  # 1750-2250
        "uranus": (
            f"{NAIF_BASE}/ura184_part-3.bsp",
            7,
            799,
        ),  # 1600-2400 (main moons segment)
        "neptune": (f"{NAIF_BASE}/nep105.bsp", 8, 899),  # 1600-2400
        "pluto": (f"{NAIF_BASE}/plu060.bsp", 9, 999),  # 1800-2200
    },
    "extended": {
        # Extended-range "xl" versions where available
        "jupiter": (f"{NAIF_BASE}/jup365.bsp", 5, 599),  # 1600-2200 (max available)
        "saturn": [  # Merge 2 files: -502 to +4500
            (f"{NAIF_BASE}/sat441xl_part-1.bsp", 6, 699),  # -502 to 2014
            (f"{NAIF_BASE}/sat441xl_part-2.bsp", 6, 699),  # 2014 to 4500
        ],
        "uranus": (f"{NAIF_BASE}/ura111xl-799.bsp", 7, 799),  # -12001 to +17000
        "neptune": (f"{NAIF_BASE}/nep097xl-899.bsp", 8, 899),  # -12001 to +17000
        "pluto": (f"{NAIF_BASE}/plu060.bsp", 9, 999),  # 1800-2200 (max available)
    },
}

# Output filenames for each tier
TIER_OUTPUT = {
    "base": "planet_centers_base.bsp",
    "medium": "planet_centers_medium.bsp",
    "extended": "planet_centers_extended.bsp",
}

# Expected coverage descriptions
TIER_COVERAGE = {
    "base": "1850-2150 (full coverage)",
    "medium": "1550-2650 (partial: Pluto 1800-2200)",
    "extended": "Partial: Jupiter 1600-2200, Saturn -502 to +4500, Uranus/Neptune full, Pluto 1800-2200",
}

# Leap seconds kernel URL (required by SPICE for time conversions)
LEAPSECONDS_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"


# =============================================================================
# SSL CONTEXT
# =============================================================================


def _get_ssl_context(verify: bool = True) -> ssl.SSLContext:
    """Create an SSL context for HTTPS downloads."""
    ctx = ssl.create_default_context()
    if not verify:
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
    return ctx


# =============================================================================
# DOWNLOAD FUNCTIONS
# =============================================================================


def download_file(
    url: str,
    dest_path: str,
    ssl_ctx: ssl.SSLContext,
    fallback_ctx: Optional[ssl.SSLContext] = None,
) -> str:
    """Download a file from URL with progress reporting.

    Uses atomic writes (tempfile + os.replace) to prevent partial files
    from being left at the destination path if the download is interrupted.
    """
    import urllib.request

    filename = os.path.basename(url)
    print(f"  Downloading {filename}...")

    # Get file size
    total_size = 0
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=30, context=ssl_ctx) as response:
            total_size = int(response.headers.get("Content-Length", 0))
    except Exception:
        if fallback_ctx:
            try:
                req = urllib.request.Request(url, method="HEAD")
                with urllib.request.urlopen(
                    req, timeout=30, context=fallback_ctx
                ) as response:
                    total_size = int(response.headers.get("Content-Length", 0))
            except Exception:
                pass

    # Download to temp file, then atomically move to dest
    dest_dir = os.path.dirname(dest_path)
    temp_fd, temp_path = tempfile.mkstemp(dir=dest_dir, suffix=".download")

    try:
        downloaded = 0
        chunk_size = 1024 * 1024  # 1 MB chunks
        ctx = ssl_ctx

        try:
            with urllib.request.urlopen(url, timeout=300, context=ctx) as response:
                with os.fdopen(temp_fd, "wb") as f:
                    while True:
                        chunk = response.read(chunk_size)
                        if not chunk:
                            break
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total_size > 0:
                            pct = downloaded / total_size * 100
                            print(
                                f"\r    {pct:5.1f}% ({downloaded / 1024 / 1024:.1f} MB)",
                                end="",
                            )
                        else:
                            print(f"\r    {downloaded / 1024 / 1024:.1f} MB", end="")
        except Exception:
            if fallback_ctx is None:
                raise
            print("    Warning: SSL verification failed; using permissive mode")
            # temp_fd was already closed by os.fdopen above on error, reopen
            temp_fd2, temp_path2 = tempfile.mkstemp(dir=dest_dir, suffix=".download")
            # Clean up first temp file
            if os.path.exists(temp_path):
                os.unlink(temp_path)
            temp_path = temp_path2
            downloaded = 0
            with urllib.request.urlopen(
                url, timeout=300, context=fallback_ctx
            ) as response:
                with os.fdopen(temp_fd2, "wb") as f:
                    while True:
                        chunk = response.read(chunk_size)
                        if not chunk:
                            break
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total_size > 0:
                            pct = downloaded / total_size * 100
                            print(
                                f"\r    {pct:5.1f}% ({downloaded / 1024 / 1024:.1f} MB)",
                                end="",
                            )

        print(f"\r    100.0% ({downloaded / 1024 / 1024:.1f} MB) - Done")
        os.replace(temp_path, dest_path)
        return dest_path

    except Exception:
        if os.path.exists(temp_path):
            os.unlink(temp_path)
        raise


def download_source_files(
    tier: str,
    cache_dir: str,
    ssl_ctx: ssl.SSLContext,
    fallback_ctx: ssl.SSLContext,
    force: bool = False,
) -> Dict[str, List[str]]:
    """Download all source SPK files for a tier.

    Returns:
        Dictionary mapping planet name to list of local file paths.
    """

    sources = TIER_SOURCES[tier]
    local_files: Dict[str, List[str]] = {}

    print(f"\n=== Downloading Source Files for '{tier}' tier ===\n")

    for planet, source_info in sources.items():
        # Handle single file or list of files
        if isinstance(source_info, list):
            file_list = source_info
        else:
            file_list = [source_info]

        local_files[planet] = []

        for url, _, _ in file_list:
            filename = os.path.basename(url)
            local_path = os.path.join(cache_dir, filename)

            if os.path.exists(local_path) and not force:
                if _is_valid_bsp(local_path):
                    print(f"  {filename}: Already cached")
                    local_files[planet].append(local_path)
                else:
                    print(f"  {filename}: Cached file is corrupted, re-downloading")
                    try:
                        os.remove(local_path)
                    except OSError:
                        pass
                    download_file(url, local_path, ssl_ctx, fallback_ctx)
                    local_files[planet].append(local_path)
            else:
                download_file(url, local_path, ssl_ctx, fallback_ctx)
                local_files[planet].append(local_path)

    # Download leap seconds kernel
    leapseconds_path = os.path.join(cache_dir, "naif0012.tls")
    if not os.path.exists(leapseconds_path):
        print("\n  Downloading leap seconds kernel...")
        download_file(LEAPSECONDS_URL, leapseconds_path, ssl_ctx, fallback_ctx)

    return local_files


# =============================================================================
# SPK EXTRACTION FUNCTIONS
# =============================================================================


def pack_spk_descriptor(dc: np.ndarray, ic: np.ndarray) -> np.ndarray:
    """Pack double and integer components into a 5-element SPK descriptor."""
    descr = np.zeros(5, dtype=np.float64)
    descr[0] = dc[0]
    descr[1] = dc[1]
    for i in range(3):
        i1 = int(ic[2 * i])
        i2 = int(ic[2 * i + 1])
        packed = struct.pack("<ii", i1, i2)
        descr[2 + i] = struct.unpack("<d", packed)[0]
    return descr


def extract_segment(
    source_file: str,
    output_handle: int,
    target_id: int,
    planet_name: str,
) -> Optional[Tuple[float, float]]:
    """Extract a planet center segment from source SPK and write to output.

    Returns:
        Tuple of (start_et, end_et) or None if failed.
    """
    import spiceypy as spice

    print(f"  Extracting {planet_name} center (NAIF {target_id})...")

    try:
        source_handle = spice.spklef(source_file)

        try:
            # Get coverage for the target body
            cover = spice.spkcov(source_file, target_id)
            n_windows = spice.wncard(cover)

            if n_windows == 0:
                print(
                    f"    WARNING: No coverage for {target_id} in {os.path.basename(source_file)}"
                )
                return None

            # Get the coverage window
            start_et, end_et = spice.wnfetd(cover, 0)
            print(f"    Coverage: ET {start_et:.0f} to {end_et:.0f}")

            # Convert to date for display
            start_utc = spice.et2utc(start_et, "ISOC", 0)
            end_utc = spice.et2utc(end_et, "ISOC", 0)
            print(f"    Date range: {start_utc} to {end_utc}")

            # Find the segment
            spice.dafbfs(source_handle)
            found_segment = False

            while True:
                found = spice.daffna()
                if not found:
                    break

                raw_descr = spice.dafgs()
                dc, ic = spice.dafus(raw_descr, 2, 6)

                if int(ic[0]) == target_id:
                    found_segment = True
                    seg_type = int(ic[3])
                    ident = spice.dafgn()

                    print(f"    Segment type: {seg_type}, ID: {ident}")

                    # Pack descriptor and copy segment
                    descr5 = pack_spk_descriptor(dc, ic)
                    spice.spksub(
                        source_handle,
                        descr5,
                        ident,
                        dc[0],
                        dc[1],
                        output_handle,
                    )
                    print("    Written to output")
                    break

            if not found_segment:
                print(f"    ERROR: Segment for {target_id} not found")
                return None

            return (start_et, end_et)

        finally:
            spice.spkuef(source_handle)

    except Exception as e:
        print(f"    ERROR: {e}")
        return None


def generate_tier_spk(
    tier: str,
    source_files: Dict[str, List[str]],
    output_path: str,
    cache_dir: str,
) -> str:
    """Generate planet_centers SPK for a specific tier."""
    import spiceypy as spice

    print("\n=== Extracting Planet Center Segments ===\n")

    # Load leap seconds kernel
    leapseconds_path = os.path.join(cache_dir, "naif0012.tls")
    spice.furnsh(leapseconds_path)

    # Remove existing output file (SPICE spkopn cannot overwrite)
    if os.path.exists(output_path):
        os.remove(output_path)

    # Create output SPK
    ifname = f"Planet Centers SPK for libephemeris ({tier} tier)"
    output_handle = spice.spkopn(output_path, ifname, 0)

    results = {}
    sources = TIER_SOURCES[tier]

    for planet, source_info in sources.items():
        # Get center/target IDs
        if isinstance(source_info, list):
            _, center_id, target_id = source_info[0]
        else:
            _, center_id, target_id = source_info

        # Extract from each source file (handle merged files)
        for source_file in source_files[planet]:
            result = extract_segment(
                source_file,
                output_handle,
                target_id,
                planet,
            )
            if result:
                if planet not in results:
                    results[planet] = result
                else:
                    # Extend range for merged files
                    old_start, old_end = results[planet]
                    new_start, new_end = result
                    results[planet] = (min(old_start, new_start), max(old_end, new_end))

    spice.spkcls(output_handle)
    spice.kclear()

    return output_path


def verify_spk(output_path: str, leapseconds_path: str) -> None:
    """Verify the generated SPK file."""
    import spiceypy as spice

    print("\n=== Verification ===\n")

    # Load leap seconds kernel for time conversions
    spice.furnsh(leapseconds_path)

    file_size = os.path.getsize(output_path)
    print(f"  Output file: {output_path}")
    print(f"  File size: {file_size / 1024 / 1024:.2f} MB")

    # List bodies
    bodies = spice.spkobj(output_path)
    print("\n  Bodies in file:")
    for body_id in sorted(bodies):
        cover = spice.spkcov(output_path, body_id)
        n_windows = spice.wncard(cover)
        if n_windows > 0:
            start_et, end_et = spice.wnfetd(cover, 0)
            start_utc = spice.et2utc(start_et, "ISOC", 0)
            end_utc = spice.et2utc(end_et, "ISOC", 0)
            print(f"    {body_id}: {start_utc} to {end_utc}")

    spice.kclear()


# =============================================================================
# MAIN
# =============================================================================


def generate_for_tier(
    tier: str,
    output_dir: Path,
    cache_dir: str,
    force: bool = False,
) -> str:
    """Generate planet_centers SPK for a single tier."""
    ssl_ctx = _get_ssl_context(verify=True)
    fallback_ctx = _get_ssl_context(verify=False)

    # Download source files
    source_files = download_source_files(tier, cache_dir, ssl_ctx, fallback_ctx, force)

    # Generate output
    output_path = str(output_dir / TIER_OUTPUT[tier])
    result = generate_tier_spk(tier, source_files, output_path, cache_dir)

    # Verify (leap seconds kernel is in cache_dir)
    leapseconds_path = os.path.join(cache_dir, "naif0012.tls")
    verify_spk(result, leapseconds_path)

    return result


def download_sources_for_tier(tier: str, cache_dir: str, force: bool = False) -> None:
    """Download only the source kernels needed to generate one tier later."""
    ssl_ctx = _get_ssl_context(verify=True)
    fallback_ctx = _get_ssl_context(verify=False)
    download_source_files(tier, cache_dir, ssl_ctx, fallback_ctx, force)


def main():
    parser = argparse.ArgumentParser(
        description="Generate tier-specific planet_centers SPK files for libephemeris",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python scripts/generate_planet_centers_spk.py --tier base
    python scripts/generate_planet_centers_spk.py --tier medium
    python scripts/generate_planet_centers_spk.py --tier extended
    python scripts/generate_planet_centers_spk.py --all
        """,
    )

    parser.add_argument(
        "--tier",
        choices=["base", "medium", "extended"],
        help="Generate SPK for a specific tier",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Generate SPK for all tiers",
    )
    parser.add_argument(
        "--download-only",
        action="store_true",
        help="Only download and cache source kernels, do not generate outputs",
    )
    parser.add_argument(
        "--force",
        "-f",
        action="store_true",
        help="Force re-download of source files",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Output directory (default: ~/.libephemeris)",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        help="Cache directory for source files (default: ~/.libephemeris/tmp/planet_centers)",
    )

    args = parser.parse_args()

    if not args.tier and not args.all:
        parser.error("Specify --tier or --all")

    print("=" * 70)
    print("PLANET CENTERS SPK GENERATOR")
    print("=" * 70)
    print()

    tiers = ["base", "medium", "extended"] if args.all else [args.tier]

    if not args.download_only:
        # Check for spiceypy only when generation is requested.
        try:
            import spiceypy

            print(f"Using spiceypy version {spiceypy.__version__}")
        except ImportError:
            print("ERROR: spiceypy is required but not installed.")
            print("Install with: pip install spiceypy>=6.0.0")
            return 1

        default_output = Path.home() / ".libephemeris"
        output_dir = args.output_dir or default_output
        print(f"Output directory: {output_dir}")
        print()

    for tier in tiers:
        print(f"\n{'=' * 70}")
        print(f"TIER: {tier}")
        print(f"Expected coverage: {TIER_COVERAGE[tier]}")
        if args.download_only:
            print("Mode: source download only")
        else:
            print(f"Output file: {TIER_OUTPUT[tier]}")
        print("=" * 70)

        if args.cache_dir:
            cache_dir = str(args.cache_dir)
        else:
            cache_dir = os.path.join(
                os.path.expanduser("~"), ".libephemeris", "tmp", "planet_centers"
            )
        os.makedirs(cache_dir, exist_ok=True)
        print(f"Cache directory: {cache_dir}")

        if args.download_only:
            download_sources_for_tier(tier, cache_dir, args.force)
        else:
            generate_for_tier(tier, output_dir, cache_dir, args.force)

    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)

    return 0


if __name__ == "__main__":
    sys.exit(main())
