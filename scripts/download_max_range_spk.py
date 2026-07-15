#!/usr/bin/env python3
"""
Download SPK files with the maximum date range supported by JPL Horizons.

Provenance:
    Project-authored batch wrapper over the public NASA/JPL Horizons SPK
    service and the registered body-ID map. The 1600-2500 interval is a measured
    service-availability boundary, not an orbital model or accuracy claim.
    Downloaded bytes retain their JPL/NAIF provenance and are never fitted here.

JPL Horizons supports SPK generation for minor bodies in the range
1600-01-01 to 2500-01-01 (verified empirically 2026-02-20).

This script downloads a single SPK file per body covering the full
1600-2500 range and saves it to the repository root, where
discover_local_spks() will find and register it.

Usage:
    python scripts/download_max_range_spk.py              # All 21 bodies
    python scripts/download_max_range_spk.py --body chiron eris  # Specific bodies
    python scripts/download_max_range_spk.py --dry-run    # Show what would be done
    python scripts/download_max_range_spk.py --output-dir /path  # Custom output dir
"""

from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# JPL Horizons maximum SPK date range for minor bodies
MAX_START = "1600-01-01"
MAX_END = "2500-01-01"


def _build_expected_skipped() -> dict[str, str]:
    """Derive name-keyed skip dict from SPK_AUTO_DOWNLOAD_BLOCKED."""
    from libephemeris.constants import SPK_AUTO_DOWNLOAD_BLOCKED, SPK_BODY_NAME_MAP
    from libephemeris.spk import _get_body_name

    result: dict[str, str] = {}
    for ipl, reason in SPK_AUTO_DOWNLOAD_BLOCKED.items():
        if ipl in SPK_BODY_NAME_MAP:
            horizons_id = SPK_BODY_NAME_MAP[ipl][0]
            key = horizons_id.rstrip(";").lower()
        else:
            name = _get_body_name(ipl) or str(ipl)
            key = name.lower()
        result[key] = reason
    return result


def _get_all_bodies() -> list[tuple[str, str, int, int]]:
    """Return (name, horizons_id, ipl, naif_id) for all SPK bodies."""
    from libephemeris.constants import SPK_BODY_NAME_MAP
    from libephemeris.spk import _get_body_name

    bodies = []
    for ipl, (horizons_id, naif_id) in SPK_BODY_NAME_MAP.items():
        name = _get_body_name(ipl)
        if name is None:
            cleaned = horizons_id.rstrip(";")
            name = cleaned if not cleaned.isdigit() else f"body_{ipl}"
        bodies.append((name, horizons_id, ipl, naif_id))
    return sorted(bodies, key=lambda b: b[0])


def _make_filename(horizons_id: str, start: str, end: str) -> str:
    """Generate filename like '2060_160001_250001.bsp'."""
    from libephemeris.spk_auto import _iso_to_jd

    jd_start = int(_iso_to_jd(start))
    jd_end = int(_iso_to_jd(end))

    # Sanitize horizons_id for filename (e.g. "Ceres;" -> "ceres_")
    safe_id = "".join(
        c if c.isalnum() or c in "-_" else "_" for c in horizons_id.lower()
    )
    safe_id = safe_id.rstrip("_")

    return f"{safe_id}_{jd_start}_{jd_end}.bsp"


def download_body(
    name: str,
    horizons_id: str,
    ipl: int,
    naif_id: int,
    output_dir: str,
    progress: str = "",
    force: bool = False,
) -> bool:
    """Download max-range SPK for a single body. Returns True on success."""
    from libephemeris.spk_auto import download_spk_from_horizons, _iso_to_jd

    filename = _make_filename(horizons_id, MAX_START, MAX_END)
    output_path = os.path.join(output_dir, filename)

    if os.path.exists(output_path) and not force:
        size_mb = os.path.getsize(output_path) / (1024 * 1024)
        print(f"  {progress}{name:14s} CACHED  {filename} ({size_mb:.1f} MB)")
        return True

    print(
        f"  {progress}{name:14s} downloading {MAX_START} to {MAX_END} ...",
        end="",
        flush=True,
    )

    try:
        jd_start = _iso_to_jd(MAX_START)
        jd_end = _iso_to_jd(MAX_END)

        download_spk_from_horizons(
            body_id=horizons_id,
            jd_start=jd_start,
            jd_end=jd_end,
            output_path=output_path,
        )

        size_mb = os.path.getsize(output_path) / (1024 * 1024)
        print(f" OK ({size_mb:.1f} MB)")
        return True
    except Exception as e:
        print(f" FAILED: {e}")
        return False


def main() -> int:
    """Download or inventory maximum-range public Horizons small-body SPKs.

    Parses optional body filters, creates the requested cache directory, and
    requests the fixed 1600–2500 service range with a configurable courtesy
    delay. Existing files are preserved unless ``--force`` is supplied;
    ``--dry-run`` performs no network writes.

    Returns:
        Zero when every selected body succeeds (or is intentionally skipped),
        otherwise one after reporting the failed bodies.
    """
    parser = argparse.ArgumentParser(
        description="Download max-range SPK files (1600-2500) for all minor bodies.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--body",
        nargs="+",
        metavar="NAME",
        help="Download only specific bodies (by name, case-insensitive)",
    )
    parser.add_argument(
        "--output-dir",
        default=os.path.join(os.path.expanduser("~"), ".libephemeris", "spk"),
        help="Output directory (default: ~/.libephemeris/spk)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download even if file already exists",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be downloaded without downloading",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=2.0,
        help="Delay in seconds between Horizons requests (default: 2.0)",
    )
    args = parser.parse_args()

    EXPECTED_SKIPPED_BODIES = _build_expected_skipped()
    all_bodies = _get_all_bodies()

    # Filter bodies if requested
    if args.body:
        requested = {b.lower() for b in args.body}
        bodies = [b for b in all_bodies if b[0].lower() in requested]
        unknown = requested - {b[0].lower() for b in bodies}
        if unknown:
            print(f"Unknown bodies: {', '.join(unknown)}", file=sys.stderr)
            print(f"Available: {', '.join(b[0] for b in all_bodies)}", file=sys.stderr)
            return 1
    else:
        bodies = all_bodies

    print(f"SPK max-range download: {MAX_START} to {MAX_END}")
    print(f"Output directory: {args.output_dir}")
    print(f"Bodies: {len(bodies)}")
    print()

    if args.dry_run:
        for name, horizons_id, ipl, naif_id in bodies:
            filename = _make_filename(horizons_id, MAX_START, MAX_END)
            output_path = os.path.join(args.output_dir, filename)
            skip_key = horizons_id.rstrip(";").lower()
            if skip_key in EXPECTED_SKIPPED_BODIES:
                exists = "SKIPPED"
            else:
                exists = "EXISTS" if os.path.exists(output_path) else "TO DOWNLOAD"
            print(f"  {name:14s} {filename:40s} [{exists}]")
        return 0

    success = 0
    failed = 0
    skipped = 0

    total = len(bodies)
    for i, (name, horizons_id, ipl, naif_id) in enumerate(bodies):
        progress = f"[{i + 1:2d}/{total}] "
        reason = EXPECTED_SKIPPED_BODIES.get(horizons_id.rstrip(";").lower())
        if reason is not None:
            print(f"  {progress}{name:14s} SKIPPED  {reason}")
            skipped += 1
            continue
        ok = download_body(
            name=name,
            horizons_id=horizons_id,
            ipl=ipl,
            naif_id=naif_id,
            output_dir=args.output_dir,
            progress=progress,
            force=args.force,
        )
        if ok:
            success += 1
        else:
            failed += 1

        # Rate limit between Horizons requests (skip for cached files)
        if i < len(bodies) - 1:
            time.sleep(args.delay)

    print()
    print(f"Done: {success} success, {failed} failed, {skipped} skipped")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
