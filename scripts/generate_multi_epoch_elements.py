#!/usr/bin/env python3
"""
Generate multi-epoch orbital elements for all minor bodies.

L1/L2 from KEPLERIAN_TODO.md: Creates MINOR_BODY_ELEMENTS_MULTI entries
at 20-year intervals (1650-2450) for all 37 bodies with SPK support.

The script:
1. Downloads/finds wide-range SPK files for each body
2. Computes heliocentric state vectors at 20-year epochs
3. Converts state vectors to osculating Keplerian elements
4. Outputs Python code for MINOR_BODY_ELEMENTS_MULTI

Usage:
    python scripts/generate_multi_epoch_elements.py                # All bodies with SPK
    python scripts/generate_multi_epoch_elements.py --body chiron  # Single body
    python scripts/generate_multi_epoch_elements.py --dry-run      # Show epochs only
    python scripts/generate_multi_epoch_elements.py --output /path/to/output.py
    python scripts/generate_multi_epoch_elements.py --spacing 20   # 20-year intervals
    python scripts/generate_multi_epoch_elements.py --start 1650 --end 2450

Requirements:
    pip install spktype21

Algorithm:
    State vector (x,y,z,vx,vy,vz) from SPK → osculating Keplerian elements
    using the vis-viva equation and standard orbital mechanics.

    The SPK provides heliocentric ICRS (J2000 equatorial) positions and
    velocities. We convert to ecliptic J2000 and then extract (a, e, i,
    ω, Ω, M₀, n) in the same conventions as MINOR_BODY_ELEMENTS.

References:
    Curtis "Orbital Mechanics for Engineering Students" Ch. 4
    Bate, Mueller, White "Fundamentals of Astrodynamics" Ch. 2
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

# Ensure libephemeris can be imported from the project root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Sun's gravitational parameter in AU³/day²
# GM_sun = 1.32712440018e20 m³/s² → convert to AU³/day²
# 1 AU = 149597870700 m, 1 day = 86400 s
# k² = GM / (AU³/day²) = 0.01720209895² (Gaussian gravitational constant squared)
MU_SUN_AU3_DAY2 = 2.9591220828559093e-04  # GM_sun in AU³/day²

# Obliquity of J2000 ecliptic
OBLIQUITY_J2000_RAD = math.radians(23.4392911)
COS_EPS = math.cos(OBLIQUITY_J2000_RAD)
SIN_EPS = math.sin(OBLIQUITY_J2000_RAD)

# AU in km
AU_KM = 149597870.7

# Default epoch grid: 1650–2450 at 20-year spacing
# Julian dates for Jan 1 of each year (approximately)
DEFAULT_START_YEAR = 1650
DEFAULT_END_YEAR = 2450
DEFAULT_SPACING = 20


@dataclass
class KeplerianElements:
    """Osculating Keplerian elements at a given epoch."""

    name: str
    epoch: float  # Julian Day TDB
    a: float  # Semi-major axis (AU)
    e: float  # Eccentricity
    i: float  # Inclination (degrees)
    omega: float  # Argument of perihelion (degrees)
    Omega: float  # Longitude of ascending node (degrees)
    M0: float  # Mean anomaly at epoch (degrees)
    n: float  # Mean motion (degrees/day)


def _year_to_jd(year: float) -> float:
    """Convert a calendar year to approximate Julian Day.

    Uses the standard formula: JD of J2000.0 = 2451545.0 (2000 Jan 1.5 TDB).
    """
    return 2451545.0 + (year - 2000.0) * 365.25


def _jd_to_approx_year(jd: float) -> float:
    """Convert Julian Day to approximate calendar year."""
    return 2000.0 + (jd - 2451545.0) / 365.25


def _equatorial_to_ecliptic(
    x_eq: float, y_eq: float, z_eq: float
) -> tuple[float, float, float]:
    """Rotate from J2000 equatorial (ICRS) to ecliptic J2000."""
    x_ecl = x_eq
    y_ecl = y_eq * COS_EPS + z_eq * SIN_EPS
    z_ecl = -y_eq * SIN_EPS + z_eq * COS_EPS
    return x_ecl, y_ecl, z_ecl


def _state_to_keplerian(
    x: float,
    y: float,
    z: float,
    vx: float,
    vy: float,
    vz: float,
    mu: float,
    epoch: float,
    name: str,
) -> Optional[KeplerianElements]:
    """Convert Cartesian state vector to Keplerian orbital elements.

    All inputs in ecliptic J2000 frame, AU and AU/day units.

    Args:
        x, y, z: Position in AU (ecliptic J2000)
        vx, vy, vz: Velocity in AU/day (ecliptic J2000)
        mu: Gravitational parameter (AU³/day²)
        epoch: Julian Day TDB
        name: Body name for the output

    Returns:
        KeplerianElements or None if conversion fails (e.g. degenerate orbit)
    """
    # Position and velocity vectors
    r_vec = (x, y, z)
    v_vec = (vx, vy, vz)

    r = math.sqrt(x**2 + y**2 + z**2)
    v = math.sqrt(vx**2 + vy**2 + vz**2)

    if r < 1e-15:
        return None

    # Specific angular momentum h = r × v
    hx = y * vz - z * vy
    hy = z * vx - x * vz
    hz = x * vy - y * vx
    h = math.sqrt(hx**2 + hy**2 + hz**2)

    if h < 1e-20:
        return None

    # Node vector n = k × h (k = [0, 0, 1])
    nx = -hy
    ny = hx
    n_mag = math.sqrt(nx**2 + ny**2)

    # Eccentricity vector e = (1/μ)[(v² - μ/r)r - (r·v)v]
    rdotv = x * vx + y * vy + z * vz
    coeff1 = (v**2 - mu / r) / mu
    coeff2 = rdotv / mu

    ex = coeff1 * x - coeff2 * vx
    ey = coeff1 * y - coeff2 * vy
    ez = coeff1 * z - coeff2 * vz
    e = math.sqrt(ex**2 + ey**2 + ez**2)

    # Semi-major axis from vis-viva: v² = μ(2/r - 1/a)
    energy = v**2 / 2.0 - mu / r
    if abs(energy) < 1e-20:
        # Parabolic — skip
        return None

    a = -mu / (2.0 * energy)

    # Inclination
    i_rad = math.acos(max(-1.0, min(1.0, hz / h)))
    i_deg = math.degrees(i_rad)

    # Longitude of ascending node (Ω)
    if n_mag > 1e-15:
        Omega_rad = math.atan2(ny, nx)
        if Omega_rad < 0:
            Omega_rad += 2.0 * math.pi
    else:
        Omega_rad = 0.0
    Omega_deg = math.degrees(Omega_rad)

    # Argument of perihelion (ω)
    if n_mag > 1e-15 and e > 1e-10:
        # cos(ω) = (n · e) / (|n| |e|)
        ndote = nx * ex + ny * ey
        omega_rad = math.acos(max(-1.0, min(1.0, ndote / (n_mag * e))))
        if ez < 0:
            omega_rad = 2.0 * math.pi - omega_rad
    elif e > 1e-10:
        # Zero inclination, use longitude of perihelion
        omega_rad = math.atan2(ey, ex)
        if omega_rad < 0:
            omega_rad += 2.0 * math.pi
    else:
        omega_rad = 0.0
    omega_deg = math.degrees(omega_rad)

    # True anomaly (ν)
    if e > 1e-10:
        # cos(ν) = (e · r) / (|e| |r|)
        edotr = ex * x + ey * y + ez * z
        nu_rad = math.acos(max(-1.0, min(1.0, edotr / (e * r))))
        if rdotv < 0:
            nu_rad = 2.0 * math.pi - nu_rad
    else:
        # Circular orbit — use argument of latitude
        if n_mag > 1e-15:
            ndotr = nx * x + ny * y
            nu_rad = math.acos(max(-1.0, min(1.0, ndotr / (n_mag * r))))
            if z < 0:
                nu_rad = 2.0 * math.pi - nu_rad
            # Subtract ω since for circular orbits ω is arbitrary
            nu_rad = nu_rad - omega_rad
        else:
            nu_rad = math.atan2(y, x)
        if nu_rad < 0:
            nu_rad += 2.0 * math.pi

    # Eccentric anomaly and mean anomaly
    if e < 1.0:
        # Elliptic orbit
        # tan(E/2) = sqrt((1-e)/(1+e)) * tan(ν/2)
        E_rad = 2.0 * math.atan2(
            math.sqrt(1.0 - e) * math.sin(nu_rad / 2.0),
            math.sqrt(1.0 + e) * math.cos(nu_rad / 2.0),
        )
        # Mean anomaly: M = E - e·sin(E)
        M_rad = E_rad - e * math.sin(E_rad)
    else:
        # Hyperbolic orbit
        # tanh(H/2) = sqrt((e-1)/(e+1)) * tan(ν/2)
        sinh_nu_half = math.sin(nu_rad / 2.0)
        cosh_nu_half = math.cos(nu_rad / 2.0)
        if abs(cosh_nu_half) < 1e-15:
            return None
        tan_nu_half = sinh_nu_half / cosh_nu_half
        tanh_H_half = math.sqrt((e - 1.0) / (e + 1.0)) * tan_nu_half
        # H = 2 * atanh(tanh_H_half)
        if abs(tanh_H_half) >= 1.0:
            return None
        H = 2.0 * math.atanh(tanh_H_half)
        M_rad = e * math.sinh(H) - H

    M_deg = math.degrees(M_rad) % 360.0

    # Mean motion (degrees/day)
    if e < 1.0:
        # n = sqrt(μ/a³) in rad/day → degrees/day
        n_rad = math.sqrt(mu / abs(a) ** 3)
        n_deg = math.degrees(n_rad)
    else:
        # Hyperbolic: n = sqrt(μ/|a|³)
        n_rad = math.sqrt(mu / abs(a) ** 3)
        n_deg = math.degrees(n_rad)

    return KeplerianElements(
        name=name,
        epoch=epoch,
        a=a,
        e=e,
        i=i_deg,
        omega=omega_deg,
        Omega=Omega_deg,
        M0=M_deg,
        n=n_deg,
    )


def _get_spk_state_vector(
    spk_file: str, jd: float
) -> Optional[tuple[float, float, float, float, float, float]]:
    """Get heliocentric ecliptic J2000 state vector from SPK at a given JD.

    Returns (x, y, z, vx, vy, vz) in AU and AU/day, or None on failure.
    """
    try:
        from spktype21 import SPKType21
    except ImportError:
        print(
            "Error: spktype21 is required. Install with: pip install spktype21",
            file=sys.stderr,
        )
        return None

    try:
        kernel = SPKType21.open(spk_file)
        try:
            center_id = kernel.segments[0].center
            target_id = kernel.segments[0].target

            # compute_type21 returns (position_km, velocity_km_per_s)
            pos_km, vel_km_s = kernel.compute_type21(center_id, target_id, jd)

            # Convert km → AU, km/s → AU/day
            x_eq = pos_km[0] / AU_KM
            y_eq = pos_km[1] / AU_KM
            z_eq = pos_km[2] / AU_KM

            vx_eq = vel_km_s[0] / AU_KM * 86400.0
            vy_eq = vel_km_s[1] / AU_KM * 86400.0
            vz_eq = vel_km_s[2] / AU_KM * 86400.0

            # Rotate from equatorial ICRS to ecliptic J2000
            x_ecl, y_ecl, z_ecl = _equatorial_to_ecliptic(x_eq, y_eq, z_eq)
            vx_ecl, vy_ecl, vz_ecl = _equatorial_to_ecliptic(vx_eq, vy_eq, vz_eq)

            return x_ecl, y_ecl, z_ecl, vx_ecl, vy_ecl, vz_ecl
        finally:
            kernel.close()
    except Exception as exc:
        # SPK may not cover this epoch
        return None


def _find_spk_file(body_name: str, horizons_id: str) -> Optional[str]:
    """Find an existing wide-range SPK file for a body.

    Searches the project root, spk/ subdirectory, and the libephemeris
    SPK cache (~/.libephemeris/spk/) for BSP files matching the body's
    Horizons ID. Prefers the widest-range file available.
    """
    project_root = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    search_dirs = [project_root, project_root / "spk"]

    # Also search the libephemeris SPK cache directory
    cache_dir = Path.home() / ".libephemeris" / "spk"
    if cache_dir.exists():
        search_dirs.append(cache_dir)

    # Normalize horizons_id for filename matching
    safe_id = "".join(
        c if c.isalnum() or c in "-_" else "_" for c in horizons_id.lower()
    ).rstrip("_")

    # Collect all matching files across all directories, prefer widest range
    all_matches: list[Path] = []

    for search_dir in search_dirs:
        if not search_dir.exists():
            continue

        # Try ID-based files (e.g., 2060_160001_250001.bsp)
        for pattern in [f"{safe_id}_*.bsp", f"{horizons_id}_*.bsp"]:
            all_matches.extend(search_dir.glob(pattern))

        # Try name-based files (e.g., chiron_*.bsp, ceres_*.bsp)
        name_lower = body_name.lower()
        for pattern in [f"{name_lower}_*.bsp"]:
            all_matches.extend(search_dir.glob(pattern))

    if not all_matches:
        return None

    # Deduplicate and prefer widest range (largest file typically = widest range)
    unique = list({str(p): p for p in all_matches}.values())
    # Sort by file size descending — widest-range SPK files are largest
    unique.sort(key=lambda p: p.stat().st_size, reverse=True)
    return str(unique[0])


def _get_spk_jd_range(spk_file: str) -> Optional[tuple[float, float]]:
    """Get the JD coverage range of an SPK file.

    Scans ALL segments and returns the overall min/max JD range,
    since SPK files often contain many consecutive segments.
    """
    try:
        from spktype21 import SPKType21

        kernel = SPKType21.open(spk_file)
        try:
            if not kernel.segments:
                return None
            jd_min = min(seg.start_jd for seg in kernel.segments)
            jd_max = max(seg.end_jd for seg in kernel.segments)
            return jd_min, jd_max
        finally:
            kernel.close()
    except Exception:
        return None


def generate_epochs(
    start_year: int = DEFAULT_START_YEAR,
    end_year: int = DEFAULT_END_YEAR,
    spacing: int = DEFAULT_SPACING,
) -> list[float]:
    """Generate Julian Day epochs at regular year intervals.

    Uses the same epoch values as the existing MINOR_BODY_ELEMENTS_MULTI
    for backward compatibility where possible.
    """
    epochs = []
    year = start_year
    while year <= end_year:
        jd = _year_to_jd(year)
        epochs.append(jd)
        year += spacing
    return epochs


def generate_elements_for_body(
    body_name: str,
    spk_file: str,
    epochs: list[float],
    verbose: bool = True,
) -> list[KeplerianElements]:
    """Generate Keplerian elements at all epochs for a single body.

    Args:
        body_name: Human-readable body name
        spk_file: Path to SPK file
        epochs: List of Julian Day epochs
        verbose: Print progress

    Returns:
        List of KeplerianElements, one per valid epoch
    """
    # Get SPK coverage range
    jd_range = _get_spk_jd_range(spk_file)
    if jd_range is None:
        if verbose:
            print(f"  {body_name}: Cannot determine SPK range", file=sys.stderr)
        return []

    spk_start, spk_end = jd_range

    elements_list = []
    for jd in epochs:
        # Skip epochs outside SPK range (with small margin)
        if jd < spk_start + 1.0 or jd > spk_end - 1.0:
            if verbose:
                yr = _jd_to_approx_year(jd)
                print(f"    {yr:.0f}: SKIP (outside SPK range)")
            continue

        state = _get_spk_state_vector(spk_file, jd)
        if state is None:
            if verbose:
                yr = _jd_to_approx_year(jd)
                print(f"    {yr:.0f}: SKIP (SPK read error)")
            continue

        x, y, z, vx, vy, vz = state
        elems = _state_to_keplerian(x, y, z, vx, vy, vz, MU_SUN_AU3_DAY2, jd, body_name)

        if elems is None:
            if verbose:
                yr = _jd_to_approx_year(jd)
                print(f"    {yr:.0f}: SKIP (degenerate orbit)")
            continue

        elements_list.append(elems)
        if verbose:
            yr = _jd_to_approx_year(jd)
            print(f"    {yr:.0f}: a={elems.a:.6f} e={elems.e:.6f} i={elems.i:.4f}")

    return elements_list


def format_elements_python(
    body_const: str,
    elements_list: list[KeplerianElements],
) -> str:
    """Format a body's multi-epoch elements as Python code.

    Generates code for inclusion in MINOR_BODY_ELEMENTS_MULTI.
    """
    lines = [f"    {body_const}: ["]

    for elem in elements_list:
        yr = _jd_to_approx_year(elem.epoch)
        lines.append("        OrbitalElements(")
        lines.append(f'            name="{elem.name}",')
        lines.append(f"            epoch={elem.epoch},")
        lines.append(f"            a={elem.a:.15g},")
        lines.append(f"            e={elem.e:.16g},")
        lines.append(f"            i={elem.i:.13g},")
        lines.append(f"            omega={elem.omega:.13g},")
        lines.append(f"            Omega={elem.Omega:.13g},")
        lines.append(f"            M0={elem.M0:.13g},")
        lines.append(f"            n={elem.n:.17g},")
        lines.append(f"        ),  # ~{yr:.0f}")

    lines.append("    ],")
    return "\n".join(lines)


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Generate multi-epoch orbital elements from SPK files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/generate_multi_epoch_elements.py                # All bodies
  python scripts/generate_multi_epoch_elements.py --body chiron  # Single body
  python scripts/generate_multi_epoch_elements.py --spacing 20   # 20-year intervals
  python scripts/generate_multi_epoch_elements.py --dry-run      # Show plan only
  python scripts/generate_multi_epoch_elements.py --output out.py  # Write to file
        """,
    )

    parser.add_argument(
        "--body",
        type=str,
        nargs="*",
        help="Specific body name(s) to generate (e.g., 'chiron', 'ceres')",
    )
    parser.add_argument(
        "--start",
        type=int,
        default=DEFAULT_START_YEAR,
        help=f"Start year (default: {DEFAULT_START_YEAR})",
    )
    parser.add_argument(
        "--end",
        type=int,
        default=DEFAULT_END_YEAR,
        help=f"End year (default: {DEFAULT_END_YEAR})",
    )
    parser.add_argument(
        "--spacing",
        type=int,
        default=DEFAULT_SPACING,
        help=f"Year spacing between epochs (default: {DEFAULT_SPACING})",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output file path (default: stdout)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show epoch grid and available SPK files without generating",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress progress output",
    )

    args = parser.parse_args()
    verbose = not args.quiet

    # Import libephemeris constants
    try:
        from libephemeris.constants import SPK_BODY_NAME_MAP
        from libephemeris.minor_bodies import MINOR_BODY_ELEMENTS
        from libephemeris.spk import _get_body_name
    except ImportError as e:
        print(f"Error importing libephemeris: {e}", file=sys.stderr)
        return 1

    # Build body list: (se_const_name, body_name, horizons_id, body_id)
    all_bodies: list[tuple[str, str, str, int]] = []
    for body_id, elem in MINOR_BODY_ELEMENTS.items():
        body_name = elem.name
        if body_id in SPK_BODY_NAME_MAP:
            horizons_id = SPK_BODY_NAME_MAP[body_id][0]
        else:
            horizons_id = ""

        # Derive bare constant name from body name
        const_name = body_name.upper()
        # Handle special cases (collision avoidance with planet/moon names)
        name_to_const = {
            "Europa": "EUROPA_AST",
            "Pandora": "PANDORA_AST",
            "Lilith": "LILITH_AST",
        }
        if body_name in name_to_const:
            const_name = name_to_const[body_name]

        all_bodies.append((const_name, body_name, horizons_id, body_id))

    # Filter by --body argument if specified
    if args.body:
        requested = {b.lower() for b in args.body}
        all_bodies = [
            (c, n, h, bid) for c, n, h, bid in all_bodies if n.lower() in requested
        ]
        if not all_bodies:
            print(f"Error: No matching bodies found for: {args.body}", file=sys.stderr)
            available = sorted(elem.name for elem in MINOR_BODY_ELEMENTS.values())
            print(f"Available: {', '.join(available)}", file=sys.stderr)
            return 1

    # Sort alphabetically by body name
    all_bodies.sort(key=lambda x: x[1])

    # Generate epoch grid
    epochs = generate_epochs(args.start, args.end, args.spacing)

    if verbose:
        print("Multi-epoch element generation")
        print(f"  Epoch range: {args.start}–{args.end}")
        print(f"  Spacing: {args.spacing} years")
        print(f"  Epochs: {len(epochs)}")
        print(f"  Bodies: {len(all_bodies)}")
        print()

    if args.dry_run:
        print("Epoch grid (JD TDB):")
        for jd in epochs:
            yr = _jd_to_approx_year(jd)
            print(f"  {yr:.0f}: JD {jd:.1f}")
        print()

        print("Body SPK availability:")
        for const_name, body_name, horizons_id, body_id in all_bodies:
            spk = _find_spk_file(body_name, horizons_id)
            if spk:
                jd_range = _get_spk_jd_range(spk)
                if jd_range:
                    yr0 = _jd_to_approx_year(jd_range[0])
                    yr1 = _jd_to_approx_year(jd_range[1])
                    n_epochs = sum(
                        1 for jd in epochs if jd_range[0] + 1 < jd < jd_range[1] - 1
                    )
                    print(
                        f"  {body_name:14s} {const_name:20s} "
                        f"SPK: {yr0:.0f}–{yr1:.0f} "
                        f"({n_epochs}/{len(epochs)} epochs)"
                    )
                else:
                    print(f"  {body_name:14s} {const_name:20s} SPK: (range unknown)")
            else:
                print(f"  {body_name:14s} {const_name:20s} NO SPK")
        return 0

    # Generate elements for each body
    output_parts: list[str] = []
    success_count = 0
    skip_count = 0

    for const_name, body_name, horizons_id, body_id in all_bodies:
        if verbose:
            print(f"  {body_name} ({const_name}):")

        spk_file = _find_spk_file(body_name, horizons_id)
        if spk_file is None:
            if verbose:
                print("    SKIP: No SPK file found")
            skip_count += 1
            continue

        if verbose:
            print(f"    SPK: {os.path.basename(spk_file)}")

        elements_list = generate_elements_for_body(
            body_name, spk_file, epochs, verbose=verbose
        )

        if not elements_list:
            if verbose:
                print("    SKIP: No valid elements generated")
            skip_count += 1
            continue

        code = format_elements_python(const_name, elements_list)
        output_parts.append(code)
        success_count += 1

        if verbose:
            print(f"    Generated {len(elements_list)} epoch entries")
            print()

    # Assemble output
    header = f"""# =============================================================================
# MULTI-EPOCH ORBITAL ELEMENTS (L1/L2)
# =============================================================================
# Generated by scripts/generate_multi_epoch_elements.py
# Epoch range: {args.start}–{args.end}, spacing: {args.spacing} years
# Source: JPL SPK type 21 state vectors → osculating Keplerian elements
#
# Each body has elements at {args.spacing}-year intervals; the consumer
# (_get_epoch_elements_blend) picks the nearest epoch and cross-fades
# positions near epoch midpoints (Hermite interpolation of osculating
# elements was tested and rejected for overshoot).
#
# Bodies with SPK data: {success_count}/{len(all_bodies)}
# Bodies without SPK: {skip_count}/{len(all_bodies)}

MINOR_BODY_ELEMENTS_MULTI: dict[int, list[OrbitalElements]] = {{
"""

    footer = "}\n"

    full_output = header + "\n".join(output_parts) + "\n" + footer

    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(full_output)
        if verbose:
            print(f"\nOutput written to {args.output}")
    else:
        print(full_output)

    if verbose:
        print(f"\nSummary: {success_count} bodies generated, {skip_count} skipped")

    return 0


if __name__ == "__main__":
    sys.exit(main())
