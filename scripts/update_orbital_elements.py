#!/usr/bin/env python3
"""
Update orbital elements from JPL Small-Body Database (SBDB).

Provenance:
    Candidate elements and epochs come directly from NASA/JPL's public SBDB API.
    Project code maps documented fields and reports or applies a textual update;
    every applied change requires review of units, frame, epoch, uncertainty,
    and targeted tests. Compatibility output and residual fitting are not
    inputs to this maintenance path.

This maintenance script fetches current orbital elements from the JPL SBDB API
for all bodies in MINOR_BODY_ELEMENTS, compares with current values, and either
updates minor_bodies.py directly or generates a report of needed changes.

RECOMMENDED: Run this script quarterly to keep orbital elements current.

Usage:
    python scripts/update_orbital_elements.py                  # Report mode (default)
    python scripts/update_orbital_elements.py --report         # Generate change report
    python scripts/update_orbital_elements.py --update         # Update minor_bodies.py
    python scripts/update_orbital_elements.py --body chiron    # Check single body
    python scripts/update_orbital_elements.py --threshold 0.01 # Custom diff threshold

Requirements:
    pip install requests

Quarterly Maintenance Schedule:
    Run this script at least once per quarter (every 3 months) to ensure
    orbital elements remain accurate. The orbital elements drift over time
    due to perturbations not captured in the Keplerian model. Recommended
    schedule:
        - January 1st
        - April 1st
        - July 1st
        - October 1st

    Or set up a cron job:
        0 0 1 1,4,7,10 * cd /path/to/libephemeris && python scripts/update_orbital_elements.py --report

API Reference:
    JPL SBDB API: https://ssd-api.jpl.nasa.gov/doc/sbdb.html
"""

import argparse
import json
import os
import re
import sys
from dataclasses import dataclass
from datetime import datetime
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from libephemeris.minor_bodies import OrbitalElements  # noqa: F401 (availability probe)

# Ensure libephemeris can be imported from the project root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Try to import requests
try:
    import requests
except ImportError:
    requests = None  # type: ignore


# JPL SBDB API endpoint
SBDB_API_URL = "https://ssd-api.jpl.nasa.gov/sbdb.api"

# Mapping from body name to JPL SBDB query identifiers
# These are the SPK-ID or designation that SBDB recognizes
BODY_SBDB_IDS: dict[str, str] = {
    "Chiron": "2060",  # Numbered asteroid 2060 Chiron
    "Pholus": "5145",  # Numbered asteroid 5145 Pholus
    "Ceres": "1",  # 1 Ceres (dwarf planet)
    "Pallas": "2",  # 2 Pallas
    "Juno": "3",  # 3 Juno
    "Vesta": "4",  # 4 Vesta
    "Eris": "136199",  # 136199 Eris (dwarf planet)
    "Sedna": "90377",  # 90377 Sedna
    "Haumea": "136108",  # 136108 Haumea (dwarf planet)
    "Makemake": "136472",  # 136472 Makemake (dwarf planet)
    "Ixion": "28978",  # 28978 Ixion
    "Orcus": "90482",  # 90482 Orcus
    "Quaoar": "50000",  # 50000 Quaoar
    "Varuna": "20000",  # 20000 Varuna (classical KBO)
    "Nessus": "7066",  # 7066 Nessus (centaur)
    "Asbolus": "8405",  # 8405 Asbolus (centaur)
    "Chariklo": "10199",  # 10199 Chariklo (centaur)
    "Gonggong": "225088",  # 225088 Gonggong (TNO, dwarf planet candidate)
    "Apophis": "99942",  # 99942 Apophis (Near-Earth asteroid)
    "Hygiea": "10",  # 10 Hygiea (fourth largest asteroid, dwarf planet candidate)
    "Interamnia": "704",  # 704 Interamnia (fifth largest asteroid)
    "Davida": "511",  # 511 Davida (seventh largest asteroid)
    "Europa": "52",  # 52 Europa (main belt asteroid, not Jupiter's moon)
    "Sylvia": "87",  # 87 Sylvia (triple asteroid system with moons Romulus and Remus)
    "Psyche": "16",  # 16 Psyche (metallic M-type asteroid, NASA Psyche mission target)
    "Eros": "433",  # 433 Eros (near-Earth asteroid, NEAR Shoemaker mission target)
    "Amor": "1221",  # 1221 Amor (prototype of the Amor near-Earth asteroid class)
    "Icarus": "1566",  # 1566 Icarus (Apollo asteroid, highly eccentric orbit)
    "Toro": "1685",  # 1685 Toro (Apollo asteroid, near-Earth asteroid)
    "Sappho": "80",  # 80 Sappho (main belt asteroid)
    "Pandora": "55",  # 55 Pandora (main belt asteroid, not Saturn moon)
    "Lilith": "1181",  # 1181 Lilith (main belt asteroid)
    "Hidalgo": "944",  # 944 Hidalgo (Centaur-class asteroid)
    "Toutatis": "4179",  # 4179 Toutatis (Apollo PHA)
    "Itokawa": "25143",  # 25143 Itokawa (Hayabusa sample return mission target)
    "Bennu": "101955",  # 101955 Bennu (OSIRIS-REx sample return mission target)
    "Ryugu": "162173",  # 162173 Ryugu (Hayabusa2 sample return mission target)
}


@dataclass
class FetchedElements:
    """Orbital elements fetched from JPL SBDB."""

    name: str
    epoch: float  # Julian Day TDB
    a: float  # Semi-major axis (AU)
    e: float  # Eccentricity
    i: float  # Inclination (degrees)
    omega: float  # Argument of perihelion (degrees)
    Omega: float  # Longitude of ascending node (degrees)
    M0: float  # Mean anomaly at epoch (degrees)
    n: float  # Mean motion (degrees/day)


def check_requests() -> bool:
    """Check if requests library is available."""
    return requests is not None


def fetch_orbital_elements(
    body_name: str, verbose: bool = False
) -> Optional[FetchedElements]:
    """
    Fetch orbital elements for a body from JPL SBDB API.

    Args:
        body_name: Name of the body (e.g., 'Chiron', 'Ceres')
        verbose: If True, print debug information

    Returns:
        FetchedElements if successful, None otherwise
    """
    if body_name not in BODY_SBDB_IDS:
        if verbose:
            print(f"  Unknown body: {body_name}", file=sys.stderr)
        return None

    sbdb_id = BODY_SBDB_IDS[body_name]

    # Build API request
    # Request orbital elements with full precision
    params = {
        "sstr": sbdb_id,  # Search string (SPK-ID or name)
        "phys-par": "false",  # Don't need physical parameters
        "full-prec": "true",  # Full precision for orbital elements
    }

    try:
        if verbose:
            print(f"  Fetching {body_name} (ID: {sbdb_id})...")

        response = requests.get(SBDB_API_URL, params=params, timeout=30)  # type: ignore[union-attr]
        response.raise_for_status()
        data = response.json()

        # Check for API errors
        if "error" in data:
            if verbose:
                print(f"  API error for {body_name}: {data['error']}", file=sys.stderr)
            return None

        # Extract orbital elements
        orbit = data.get("orbit", {})
        elements = orbit.get("elements", [])

        # Build a dict of element name -> value
        elem_dict: dict[str, float] = {}
        for elem in elements:
            elem_dict[elem["name"]] = float(elem["value"])

        # Get epoch (Julian Day TDB)
        epoch_jd = float(orbit.get("epoch", 0))

        # Map SBDB element names to our structure
        # SBDB uses: e, a, q, i, om, w, ma, n, tp, per, ...
        # Our structure: a, e, i, omega (w), Omega (om), M0 (ma), n
        fetched = FetchedElements(
            name=body_name,
            epoch=epoch_jd,
            a=elem_dict.get("a", 0.0),  # Semi-major axis
            e=elem_dict.get("e", 0.0),  # Eccentricity
            i=elem_dict.get("i", 0.0),  # Inclination
            omega=elem_dict.get("w", 0.0),  # Argument of perihelion (w in SBDB)
            Omega=elem_dict.get("om", 0.0),  # Longitude of ascending node (om in SBDB)
            M0=elem_dict.get("ma", 0.0),  # Mean anomaly (ma in SBDB)
            n=elem_dict.get("n", 0.0),  # Mean motion
        )

        if verbose:
            print(f"    Epoch: JD {epoch_jd:.1f}")

        return fetched

    except requests.exceptions.RequestException as e:  # type: ignore[union-attr]
        if verbose:
            print(f"  Network error for {body_name}: {e}", file=sys.stderr)
        return None
    except (KeyError, ValueError, json.JSONDecodeError) as e:
        if verbose:
            print(f"  Parse error for {body_name}: {e}", file=sys.stderr)
        return None


def compare_elements(
    current: "OrbitalElements", fetched: FetchedElements, threshold: float = 0.001
) -> dict[str, tuple[float, float, float]]:
    """
    Compare current and fetched orbital elements.

    Args:
        current: Current OrbitalElements from minor_bodies.py
        fetched: Newly fetched elements from SBDB
        threshold: Relative difference threshold for reporting (default 0.1%)

    Returns:
        Dict of element_name -> (current_value, fetched_value, percent_diff)
        Only includes elements that differ by more than threshold
    """
    differences: dict[str, tuple[float, float, float]] = {}

    # Compare each element
    comparisons = [
        ("epoch", current.epoch, fetched.epoch),
        ("a", current.a, fetched.a),
        ("e", current.e, fetched.e),
        ("i", current.i, fetched.i),
        ("omega", current.omega, fetched.omega),
        ("Omega", current.Omega, fetched.Omega),
        ("M0", current.M0, fetched.M0),
        ("n", current.n, fetched.n),
    ]

    for name, curr_val, new_val in comparisons:
        if curr_val == 0:
            # Avoid division by zero
            if new_val != 0:
                differences[name] = (curr_val, new_val, float("inf"))
        else:
            rel_diff = abs(new_val - curr_val) / abs(curr_val)
            if rel_diff > threshold:
                pct_diff = rel_diff * 100
                differences[name] = (curr_val, new_val, pct_diff)

    return differences


def generate_report(
    all_differences: dict[str, dict[str, tuple[float, float, float]]],
    fetched_epochs: dict[str, float],
) -> str:
    """
    Generate a human-readable report of orbital element differences.

    Args:
        all_differences: Dict of body_name -> element differences
        fetched_epochs: Dict of body_name -> new epoch JD

    Returns:
        Formatted report string
    """
    lines = []
    lines.append("=" * 80)
    lines.append("ORBITAL ELEMENTS UPDATE REPORT")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("=" * 80)
    lines.append("")

    bodies_needing_update = 0
    for body_name, diffs in all_differences.items():
        if diffs:
            bodies_needing_update += 1

    if bodies_needing_update == 0:
        lines.append("All orbital elements are up to date! No changes needed.")
    else:
        lines.append(f"Found {bodies_needing_update} bodies with outdated elements:")
        lines.append("")

        for body_name, diffs in sorted(all_differences.items()):
            if not diffs:
                continue

            lines.append(f"  {body_name}:")
            if body_name in fetched_epochs:
                lines.append(f"    New epoch: JD {fetched_epochs[body_name]:.1f}")

            for elem_name, (curr, new, pct) in sorted(diffs.items()):
                if pct == float("inf"):
                    pct_str = "INF"
                else:
                    pct_str = f"{pct:.4f}%"
                lines.append(
                    f"    {elem_name:8s}: {curr:20.14g} -> {new:20.14g} ({pct_str})"
                )
            lines.append("")

    lines.append("=" * 80)
    lines.append("")
    lines.append("To update minor_bodies.py, run:")
    lines.append("  python scripts/update_orbital_elements.py --update")
    lines.append("")

    return "\n".join(lines)


def generate_python_code(fetched: FetchedElements, body_const: str) -> str:
    """
    Generate Python code for updating a body's orbital elements.

    Args:
        fetched: Fetched orbital elements
        body_const: Constant name (e.g., 'CHIRON')

    Returns:
        Python code string for the OrbitalElements entry
    """
    lines = []
    lines.append(f"    {body_const}: OrbitalElements(")
    lines.append(f'        name="{fetched.name}",')
    lines.append(f"        epoch={fetched.epoch},")
    lines.append(f"        a={fetched.a},  # AU")
    lines.append(f"        e={fetched.e},")
    lines.append(f"        i={fetched.i},")
    lines.append(f"        omega={fetched.omega},  # argument of perihelion")
    lines.append(f"        Omega={fetched.Omega},  # longitude of ascending node")
    lines.append(f"        M0={fetched.M0},")
    lines.append(f"        n={fetched.n},")
    lines.append("    ),")
    return "\n".join(lines)


def update_minor_bodies_file(
    fetched_elements: dict[str, FetchedElements],
    dry_run: bool = True,
    verbose: bool = True,
) -> bool:
    """
    Update the minor_bodies.py file with new orbital elements.

    Args:
        fetched_elements: Dict of body_name -> FetchedElements
        dry_run: If True, only show what would be changed
        verbose: If True, print progress

    Returns:
        True if update was successful (or dry run completed)
    """
    # Path to minor_bodies.py
    minor_bodies_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "libephemeris",
        "minor_bodies.py",
    )

    if not os.path.exists(minor_bodies_path):
        print(f"Error: Cannot find {minor_bodies_path}", file=sys.stderr)
        return False

    # Read current file
    with open(minor_bodies_path, "r", encoding="utf-8") as f:
        content = f.read()

    # Track changes
    changes_made = 0

    # For each body, find and replace its orbital elements
    for body_name, fetched in fetched_elements.items():
        # Map body name to constant name
        const_name = f"SE_{body_name.upper()}"

        # Pattern to match the body's OrbitalElements entry
        # This matches from "SE_BODYNAME: OrbitalElements(" to the closing "),"
        pattern = rf"({const_name}: OrbitalElements\()([^)]+)(\),)"

        match = re.search(pattern, content, re.DOTALL)
        if not match:
            if verbose:
                print(f"  Warning: Could not find {const_name} in minor_bodies.py")
            continue

        # Generate new orbital elements code (just the inner part)
        new_inner = f"""
        name="{fetched.name}",
        epoch={fetched.epoch},
        a={fetched.a},  # AU
        e={fetched.e},
        i={fetched.i},
        omega={fetched.omega},  # argument of perihelion
        Omega={fetched.Omega},  # longitude of ascending node
        M0={fetched.M0},
        n={fetched.n},
    """

        new_entry = f"{match.group(1)}{new_inner}{match.group(3)}"
        content = content[: match.start()] + new_entry + content[match.end() :]
        changes_made += 1

        if verbose:
            print(f"  Updated {body_name}")

    if changes_made == 0:
        if verbose:
            print("No changes to make.")
        return True

    # Update the source comment with new date
    today = datetime.now().strftime("%Y-%b-%d")
    content = re.sub(
        r"# Source: NASA JPL Small-Body Database \(sbdb\.api\), retrieved \d{4}-[A-Za-z]{3}-\d{2}",
        f"# Source: NASA JPL Small-Body Database (sbdb.api), retrieved {today}",
        content,
    )

    if dry_run:
        print(f"\nDry run: Would update {changes_made} bodies in minor_bodies.py")
        print("Run with --update to apply changes.")
        return True

    # Write updated file
    with open(minor_bodies_path, "w", encoding="utf-8") as f:
        f.write(content)

    if verbose:
        print(f"\nUpdated {changes_made} bodies in {minor_bodies_path}")

    return True


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Update orbital elements from JPL Small-Body Database.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/update_orbital_elements.py                  # Report mode (safe)
  python scripts/update_orbital_elements.py --report         # Same as above
  python scripts/update_orbital_elements.py --update         # Update minor_bodies.py
  python scripts/update_orbital_elements.py --body chiron    # Check single body
  python scripts/update_orbital_elements.py --dry-run        # Show what would change

Quarterly Maintenance:
  This script should be run at least once per quarter to keep orbital
  elements current. Set up a reminder or cron job for:
    - January 1st
    - April 1st
    - July 1st
    - October 1st
        """,
    )

    parser.add_argument(
        "--report",
        action="store_true",
        default=True,
        help="Generate a report of needed changes (default)",
    )
    parser.add_argument(
        "--update",
        action="store_true",
        help="Update minor_bodies.py with new values",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be changed without modifying files",
    )
    parser.add_argument(
        "--body",
        type=str,
        help="Check a single body (e.g., 'chiron', 'ceres')",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.001,
        help="Relative difference threshold for reporting (default: 0.001 = 0.1%%)",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output results as JSON",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress progress output",
    )

    args = parser.parse_args()
    verbose = not args.quiet

    # Check requests library
    if not check_requests():
        print(
            "Error: requests library is required.",
            file=sys.stderr,
        )
        print("Install it with: pip install requests", file=sys.stderr)
        return 1

    # Import MINOR_BODY_ELEMENTS
    try:
        from libephemeris.minor_bodies import (  # noqa: F401 (probe)
            MINOR_BODY_ELEMENTS,
            OrbitalElements,
        )
        from libephemeris.constants import (
            CHIRON,
            PHOLUS,
            CERES,
            PALLAS,
            JUNO,
            VESTA,
            ERIS,
            SEDNA,
            HAUMEA,
            MAKEMAKE,
            IXION,
            ORCUS,
            QUAOAR,
        )
    except ImportError as e:
        print(f"Error importing libephemeris: {e}", file=sys.stderr)
        return 1

    # Map body IDs to names and constants
    body_id_to_name = {
        CHIRON: "Chiron",
        PHOLUS: "Pholus",
        CERES: "Ceres",
        PALLAS: "Pallas",
        JUNO: "Juno",
        VESTA: "Vesta",
        ERIS: "Eris",
        SEDNA: "Sedna",
        HAUMEA: "Haumea",
        MAKEMAKE: "Makemake",
        IXION: "Ixion",
        ORCUS: "Orcus",
        QUAOAR: "Quaoar",
    }

    # Determine which bodies to check
    if args.body:
        body_name = args.body.capitalize()
        if body_name not in BODY_SBDB_IDS:
            print(f"Error: Unknown body '{args.body}'", file=sys.stderr)
            print(f"Available: {', '.join(BODY_SBDB_IDS.keys())}", file=sys.stderr)
            return 1
        bodies_to_check = [body_name]
    else:
        bodies_to_check = list(BODY_SBDB_IDS.keys())

    if verbose:
        print(f"Fetching orbital elements for {len(bodies_to_check)} bodies...")
        print()

    # Fetch and compare all bodies
    all_differences: dict[str, dict[str, tuple[float, float, float]]] = {}
    fetched_elements: dict[str, FetchedElements] = {}
    fetched_epochs: dict[str, float] = {}

    for body_name in bodies_to_check:
        fetched = fetch_orbital_elements(body_name, verbose=verbose)
        if fetched is None:
            if verbose:
                print(f"  {body_name}: Failed to fetch", file=sys.stderr)
            continue

        fetched_elements[body_name] = fetched
        fetched_epochs[body_name] = fetched.epoch

        # Find corresponding current elements
        body_id = None
        for bid, bname in body_id_to_name.items():
            if bname == body_name:
                body_id = bid
                break

        if body_id is None or body_id not in MINOR_BODY_ELEMENTS:
            if verbose:
                print(f"  {body_name}: Not in MINOR_BODY_ELEMENTS", file=sys.stderr)
            continue

        current = MINOR_BODY_ELEMENTS[body_id]
        diffs = compare_elements(current, fetched, args.threshold)
        all_differences[body_name] = diffs

    # Output results
    if args.json:
        # JSON output
        output = {
            "timestamp": datetime.now().isoformat(),
            "bodies_checked": len(bodies_to_check),
            "bodies_with_changes": sum(1 for d in all_differences.values() if d),
            "differences": {
                name: {
                    elem: {"current": curr, "new": new, "percent_diff": pct}
                    for elem, (curr, new, pct) in diffs.items()
                }
                for name, diffs in all_differences.items()
            },
            "fetched_epochs": fetched_epochs,
        }
        print(json.dumps(output, indent=2))
    elif args.update:
        # Update mode
        if not fetched_elements:
            print("No elements fetched. Cannot update.", file=sys.stderr)
            return 1
        success = update_minor_bodies_file(
            fetched_elements, dry_run=args.dry_run, verbose=verbose
        )
        return 0 if success else 1
    else:
        # Report mode (default)
        report = generate_report(all_differences, fetched_epochs)
        print(report)

    return 0


if __name__ == "__main__":
    sys.exit(main())
