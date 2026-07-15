#!/usr/bin/env python3
"""
Generate tier-specific planet_centers SPK files for libephemeris.

Provenance:
    Source segments are downloaded from NASA/JPL NAIF's public satellite SPK
    archive and selected by documented NAIF center IDs. Project code copies the
    required segments into compact tier artifacts without fitting a correction
    model. Source URLs, kernel identities, coverage intersection, hashes, SPK
    structure, and verification procedure are documented in
    ``docs/methodology/planet-centers-spk.md``.

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
    planet_centers_medium.bsp    (~191 MB; coverage varies by planet)
    planet_centers_extended.bsp  (~80-100 MB, -12000 to +17000 partial)

REQUIREMENTS
------------
    - spiceypy >= 6.0.0 (pip install spiceypy)
    - Internet connection (to download source SPK files)
    - Disk space: ~500 MB (base), ~6.5 GB (medium/extended)

TIER DETAILS
------------
    base:     Uses old compact satellite SPKs (jup204, sat319, etc.)
              Coverage: 1850-2150 for all planets

    medium:   Uses the widest JPL planet-center SPKs that cover the tier.
              Coverage is full for Saturn/Uranus/Neptune and partial for
              Jupiter (1600-2200) and Pluto (1800-2200).

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
import hashlib
import os
import ssl
import struct
import sys
import tempfile
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple, Union

import numpy as np
import certifi


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
        # Jupiter and Pluto have no JPL planet-center kernel spanning the
        # complete DE440 interval. Use JPL's widest published center products
        # and report the bounded gaps explicitly. The XL kernels cover the
        # complete medium interval for Saturn, Uranus, and Neptune.
        "jupiter": (f"{NAIF_BASE}/jup365.bsp", 5, 599),  # 1600-2200
        "saturn": [
            (f"{NAIF_BASE}/sat441xl_part-1.bsp", 6, 699),
            (f"{NAIF_BASE}/sat441xl_part-2.bsp", 6, 699),
        ],  # -502 to +4500
        "uranus": (f"{NAIF_BASE}/ura111xl-799.bsp", 7, 799),
        "neptune": (f"{NAIF_BASE}/nep097xl-899.bsp", 8, 899),
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
    "medium": ("Saturn/Uranus/Neptune 1550-2650; Jupiter 1600-2200; Pluto 1800-2200"),
    "extended": "Partial: Jupiter 1600-2200, Saturn -502 to +4500, Uranus/Neptune full, Pluto 1800-2200",
}

# Output clipping windows. XL source kernels extend far beyond the base and
# medium DE tiers; copying those complete descriptors would bloat the release
# artifact and can make SPICE's human-readable UTC formatter reject extreme
# years. Segment states are clipped with spksub() to the relevant tier window.
# Extended intentionally retains each source's full independently published
# coverage.
TIER_CLIP_DATES: Dict[str, Optional[Tuple[str, str]]] = {
    "base": ("1850 JAN 01 00:00:00 TDB", "2150 JAN 01 00:00:00 TDB"),
    "medium": ("1550 JAN 01 00:00:00 TDB", "2650 JAN 01 00:00:00 TDB"),
    "extended": None,
}

# Required continuous coverage for reviewed release artifacts.  Medium is
# intentionally asymmetric because no published Jupiter- or Pluto-center SPK
# spans the complete DE440 interval.  Verification converts these TDB dates to
# ET with the same loaded NAIF leap-seconds kernel used by the generator.
TIER_REQUIRED_COVERAGE_DATES: Dict[str, Dict[int, Tuple[str, str, str]]] = {
    "medium": {
        599: (
            "Jupiter",
            "1600 JAN 01 00:00:00 TDB",
            "2200 JAN 01 00:00:00 TDB",
        ),
        699: (
            "Saturn",
            "1550 JAN 01 00:00:00 TDB",
            "2650 JAN 01 00:00:00 TDB",
        ),
        799: (
            "Uranus",
            "1550 JAN 01 00:00:00 TDB",
            "2650 JAN 01 00:00:00 TDB",
        ),
        899: (
            "Neptune",
            "1550 JAN 01 00:00:00 TDB",
            "2650 JAN 01 00:00:00 TDB",
        ),
        999: (
            "Pluto",
            "1800 JAN 01 00:00:00 TDB",
            "2200 JAN 01 00:00:00 TDB",
        ),
    },
}

# SPK descriptor endpoints and str2et() results are double-precision ET
# seconds.  One second is generous for numeric round-off while remaining far
# too narrow to hide a scientifically meaningful truncation.
_COVERAGE_ENDPOINT_TOLERANCE_ET = 1.0

# Leap seconds kernel URL (required by SPICE for time conversions)
LEAPSECONDS_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"

# Exact official JPL/NAIF inputs used for the reviewed medium artifact. The
# extended tier uses the same source products. Pinning the source bytes makes
# a future in-place archive update fail visibly instead of silently producing
# a different release artifact.
REVIEWED_SOURCE_SHA256 = {
    "jup365.bsp": "dbf016c01ba4d022154838000cf3f06962cf958ddc503a366f7fe8f81495c5cb",
    "sat441xl_part-1.bsp": "45e6f6687cf9a4d7c94c206c13650d35583740019df939c19f2ee0f0f7d2750d",
    "sat441xl_part-2.bsp": "d021819387e6045e2cd2819a7907ca62cc8cb3091ca925dfe894d2e0f4f7af1e",
    "ura111xl-799.bsp": "fa55665eb129613574a5d029bf2cc7b521d3f0d2edc7de5790add3cf950714bf",
    "nep097xl-899.bsp": "95e815cb8584e530d3c83f35e39af4bf67d0874f01b65538cbe9dbf279818df8",
    "plu060.bsp": "dfbb102491a26ed41ae08ca3f8963f22f0219df1d8f265ab87b9ad825a826fc6",
    "naif0012.tls": "678e32bdb5a744117a467cd9601cd6b373f0e9bc9bbde1371d5eee39600a039b",
}


# =============================================================================
# SSL CONTEXT
# =============================================================================


def _get_ssl_context() -> ssl.SSLContext:
    """Create a certificate-verifying SSL context for source downloads."""
    return ssl.create_default_context(cafile=certifi.where())


# =============================================================================
# DOWNLOAD FUNCTIONS
# =============================================================================


def download_file(
    url: str,
    dest_path: str,
    ssl_ctx: ssl.SSLContext,
    expected_sha256: Optional[str] = None,
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
        pass

    # Download to temp file, then atomically move to dest
    dest_dir = os.path.dirname(dest_path)
    temp_fd, temp_path = tempfile.mkstemp(dir=dest_dir, suffix=".download")

    try:
        downloaded = 0
        digest = hashlib.sha256()
        chunk_size = 1024 * 1024  # 1 MB chunks
        ctx = ssl_ctx

        with urllib.request.urlopen(url, timeout=300, context=ctx) as response:
            with os.fdopen(temp_fd, "wb") as f:
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    digest.update(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        pct = downloaded / total_size * 100
                        print(
                            f"\r    {pct:5.1f}% ({downloaded / 1024 / 1024:.1f} MB)",
                            end="",
                        )
                    else:
                        print(f"\r    {downloaded / 1024 / 1024:.1f} MB", end="")

        actual_sha256 = digest.hexdigest()
        if expected_sha256 is not None and actual_sha256 != expected_sha256:
            raise ValueError(
                f"Source hash mismatch for {filename}: expected "
                f"{expected_sha256}, got {actual_sha256}"
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
            expected_sha256 = REVIEWED_SOURCE_SHA256.get(filename)

            if os.path.exists(local_path) and not force:
                hash_matches = expected_sha256 is None or (
                    _sha256_file(local_path) == expected_sha256
                )
                if _is_valid_bsp(local_path) and hash_matches:
                    print(f"  {filename}: Already cached")
                    local_files[planet].append(local_path)
                else:
                    print(
                        f"  {filename}: Cached file is corrupt or has the wrong "
                        "source hash, re-downloading"
                    )
                    try:
                        os.remove(local_path)
                    except OSError:
                        pass
                    download_file(url, local_path, ssl_ctx, expected_sha256)
                    local_files[planet].append(local_path)
            else:
                download_file(url, local_path, ssl_ctx, expected_sha256)
                local_files[planet].append(local_path)

    # Download leap seconds kernel
    leapseconds_path = os.path.join(cache_dir, "naif0012.tls")
    leapseconds_sha256 = REVIEWED_SOURCE_SHA256["naif0012.tls"]
    leapseconds_valid = os.path.exists(leapseconds_path) and (
        _sha256_file(leapseconds_path) == leapseconds_sha256
    )
    if not leapseconds_valid:
        print("\n  Downloading leap seconds kernel...")
        download_file(
            LEAPSECONDS_URL,
            leapseconds_path,
            ssl_ctx,
            leapseconds_sha256,
        )

    return local_files


def _sha256_file(path: str) -> str:
    """Return the SHA-256 of a cached source without loading it into memory."""
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def _format_et(spice, et: float) -> str:
    """Format an ET for logs without making extraction depend on UTC limits."""
    try:
        return str(spice.et2utc(et, "ISOC", 0))
    except Exception:
        # CSPICE et2utc() supports a narrower civil-year range than several
        # official XL satellite kernels. ET is the native SPK coordinate and
        # remains unambiguous for diagnostic output.
        return f"ET {float(et):.3f}"


def extract_segment(
    source_file: str,
    output_handle: int,
    target_id: int,
    planet_name: str,
    clip_start_et: Optional[float] = None,
    clip_end_et: Optional[float] = None,
    center_id: Optional[int] = None,
) -> Optional[Tuple[float, float]]:
    """Extract every planet-center segment from a source SPK.

    JPL kernels commonly split a target across several consecutive DAF
    segments. Copying only the first descriptor silently truncated coverage
    (for example at 1997 for Jupiter and 2014 for Saturn). Every matching
    descriptor must be preserved; Skyfield/SPICE select the applicable one at
    evaluation time.

    Returns:
        Tuple of (start_et, end_et) or None if failed.
    """
    import spiceypy as spice

    print(f"  Extracting {planet_name} center (NAIF {target_id})...")

    try:
        source_handle = spice.spklef(source_file)

        try:
            # Report the source's declared target coverage. The authoritative
            # copied interval below is accumulated from the descriptors that
            # are actually written.
            cover = spice.spkcov(source_file, target_id)
            n_windows = spice.wncard(cover)

            if n_windows == 0:
                print(
                    f"    WARNING: No coverage for {target_id} in {os.path.basename(source_file)}"
                )
                return None

            for window_index in range(n_windows):
                start_et, end_et = spice.wnfetd(cover, window_index)
                start_utc = _format_et(spice, start_et)
                end_utc = _format_et(spice, end_et)
                print(
                    f"    Coverage window {window_index + 1}/{n_windows}: "
                    f"{start_utc} to {end_utc}"
                )

            # Find the segment
            spice.dafbfs(source_handle)
            copied_count = 0
            copied_start: Optional[float] = None
            copied_end: Optional[float] = None

            while True:
                found = spice.daffna()
                if not found:
                    break

                raw_descr = spice.dafgs()
                dc, ic = spice.dafus(raw_descr, 2, 6)

                descriptor_matches = int(ic[0]) == target_id and (
                    center_id is None or int(ic[1]) == center_id
                )
                if descriptor_matches:
                    seg_type = int(ic[3])
                    ident = spice.dafgn()

                    segment_start = float(dc[0])
                    segment_end = float(dc[1])
                    copy_start = (
                        segment_start
                        if clip_start_et is None
                        else max(segment_start, clip_start_et)
                    )
                    copy_end = (
                        segment_end
                        if clip_end_et is None
                        else min(segment_end, clip_end_et)
                    )
                    if copy_start >= copy_end:
                        continue

                    print(
                        f"    Copying segment {copied_count + 1}: "
                        f"type {seg_type}, ID: {ident}, "
                        f"{_format_et(spice, copy_start)} to "
                        f"{_format_et(spice, copy_end)}"
                    )

                    # Pack descriptor and copy segment
                    descr5 = pack_spk_descriptor(dc, ic)
                    spice.spksub(
                        source_handle,
                        descr5,
                        ident,
                        copy_start,
                        copy_end,
                        output_handle,
                    )
                    copied_count += 1
                    copied_start = (
                        copy_start
                        if copied_start is None
                        else min(copied_start, copy_start)
                    )
                    copied_end = (
                        copy_end if copied_end is None else max(copied_end, copy_end)
                    )

            if copied_count == 0 or copied_start is None or copied_end is None:
                print(f"    ERROR: Segment for {target_id} not found")
                return None

            print(f"    Written {copied_count} segment(s) to output")
            return (copied_start, copied_end)

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

    # Create output SPK. Any source extraction failure invalidates the whole
    # artifact; a structurally valid but partial kernel is worse than no file
    # because its filename would overstate the reviewed tier coverage.
    ifname = f"Planet Centers SPK for libephemeris ({tier} tier)"
    output_handle = spice.spkopn(output_path, ifname, 0)
    output_open = True
    try:
        results = {}
        failed_sources = []
        sources = TIER_SOURCES[tier]
        clip_dates = TIER_CLIP_DATES[tier]
        if clip_dates is None:
            clip_start_et = clip_end_et = None
        else:
            clip_start_et = float(spice.str2et(clip_dates[0]))
            clip_end_et = float(spice.str2et(clip_dates[1]))

        for planet, source_info in sources.items():
            # Get center/target IDs
            if isinstance(source_info, list):
                _, center_id, target_id = source_info[0]
            else:
                _, center_id, target_id = source_info

            # Extract from each source file (handle merged files). Every
            # configured part is required, including both Saturn XL halves.
            for source_file in source_files[planet]:
                result = extract_segment(
                    source_file,
                    output_handle,
                    target_id,
                    planet,
                    clip_start_et,
                    clip_end_et,
                    center_id,
                )
                if result is None:
                    failed_sources.append(os.path.basename(source_file))
                    continue
                if planet not in results:
                    results[planet] = result
                else:
                    # Extend range for merged files
                    old_start, old_end = results[planet]
                    new_start, new_end = result
                    results[planet] = (
                        min(old_start, new_start),
                        max(old_end, new_end),
                    )

        missing_planets = sorted(set(sources) - set(results))
        if failed_sources or missing_planets:
            details = []
            if failed_sources:
                details.append("failed sources: " + ", ".join(failed_sources))
            if missing_planets:
                details.append("missing planets: " + ", ".join(missing_planets))
            raise RuntimeError(
                "Incomplete planet-center artifact (" + "; ".join(details) + ")"
            )

        spice.spkcls(output_handle)
        output_open = False
        return output_path
    except BaseException:
        if output_open:
            try:
                spice.spkcls(output_handle)
            except Exception:
                pass
        try:
            os.remove(output_path)
        except FileNotFoundError:
            pass
        raise
    finally:
        spice.kclear()


def _validate_required_coverage(
    tier: str,
    coverage_windows: Dict[int, List[Tuple[float, float]]],
    to_et: Callable[[str], float],
) -> None:
    """Fail when a tier artifact does not span every required interval."""
    requirements = TIER_REQUIRED_COVERAGE_DATES.get(tier)
    if requirements is None:
        return

    problems = []
    tolerance = _COVERAGE_ENDPOINT_TOLERANCE_ET
    for target_id, (planet, start_text, end_text) in requirements.items():
        required_start = float(to_et(start_text))
        required_end = float(to_et(end_text))
        windows = sorted(coverage_windows.get(target_id, []))

        # Walk the union rather than comparing only min/max endpoints: two
        # distant descriptors with a gap must not masquerade as full coverage.
        covered_until: Optional[float] = None
        for start_et, end_et in windows:
            if end_et < required_start - tolerance:
                continue
            if covered_until is None:
                if start_et > required_start + tolerance:
                    break
                covered_until = end_et
            elif start_et <= covered_until + tolerance:
                covered_until = max(covered_until, end_et)
            else:
                break

            if covered_until >= required_end - tolerance:
                break

        if covered_until is not None and covered_until >= required_end - tolerance:
            continue

        if windows:
            actual_start = min(window[0] for window in windows)
            actual_end = max(window[1] for window in windows)
            actual = f"actual ET envelope {actual_start:.3f} to {actual_end:.3f}"
        else:
            actual = "no coverage windows"
        problems.append(
            f"{planet} ({target_id}) requires ET {required_start:.3f} to "
            f"{required_end:.3f}; {actual}"
        )

    if problems:
        raise RuntimeError(
            f"Generated SPK fails required {tier} coverage: " + "; ".join(problems)
        )


def verify_spk(output_path: str, leapseconds_path: str, tier: str) -> None:
    """Verify the generated SPK inventory and tier coverage."""
    import spiceypy as spice

    print("\n=== Verification ===\n")

    # Load leap seconds kernel for time conversions
    spice.furnsh(leapseconds_path)

    try:
        file_size = os.path.getsize(output_path)
        print(f"  Output file: {output_path}")
        print(f"  File size: {file_size / 1024 / 1024:.2f} MB")

        # Verify exact center/target inventory independently through jplephem.
        from jplephem.spk import SPK

        expected_pairs = {(5, 599), (6, 699), (7, 799), (8, 899), (9, 999)}
        with SPK.open(output_path) as kernel:
            actual_pairs = {
                (segment.center, segment.target) for segment in kernel.segments
            }
        if actual_pairs != expected_pairs:
            raise RuntimeError(
                "Generated SPK has unexpected center/target inventory: "
                f"expected {sorted(expected_pairs)}, got {sorted(actual_pairs)}"
            )

        # List and retain every disjoint coverage window for validation.
        bodies = spice.spkobj(output_path)
        if set(bodies) != {target for _, target in expected_pairs}:
            raise RuntimeError(
                f"Generated SPK has unexpected target IDs: {sorted(bodies)}"
            )
        coverage_windows: Dict[int, List[Tuple[float, float]]] = {}
        print("\n  Bodies in file:")
        for body_id in sorted(bodies):
            cover = spice.spkcov(output_path, body_id)
            n_windows = spice.wncard(cover)
            windows = []
            for window_index in range(n_windows):
                start_et, end_et = spice.wnfetd(cover, window_index)
                start_et = float(start_et)
                end_et = float(end_et)
                windows.append((start_et, end_et))
                start_utc = _format_et(spice, start_et)
                end_utc = _format_et(spice, end_et)
                print(
                    f"    {body_id} window {window_index + 1}/{n_windows}: "
                    f"{start_utc} to {end_utc}"
                )
            coverage_windows[int(body_id)] = windows

        _validate_required_coverage(tier, coverage_windows, spice.str2et)
    finally:
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
    ssl_ctx = _get_ssl_context()

    output_dir.mkdir(parents=True, exist_ok=True)

    # Download source files
    source_files = download_source_files(tier, cache_dir, ssl_ctx, force)

    # Generate output
    output_path = str(output_dir / TIER_OUTPUT[tier])
    result = generate_tier_spk(tier, source_files, output_path, cache_dir)

    # Verify (leap seconds kernel is in cache_dir)
    leapseconds_path = os.path.join(cache_dir, "naif0012.tls")
    verify_spk(result, leapseconds_path, tier)

    return result


def download_sources_for_tier(tier: str, cache_dir: str, force: bool = False) -> None:
    """Download only the source kernels needed to generate one tier later."""
    ssl_ctx = _get_ssl_context()
    download_source_files(tier, cache_dir, ssl_ctx, force)


def main():
    """Generate or download inputs for tier-specific planet-center SPKs.

    The command makes the tier, source-kernel cache, output directory, and
    verification phase explicit.  Barycentric states originate in the public
    NAIF/JPL kernels named by this module; segment assembly and tier partition
    are reproducible project packaging decisions.
    """
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
