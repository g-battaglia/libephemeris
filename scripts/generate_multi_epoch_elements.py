#!/usr/bin/env python3
"""
Generate multi-epoch orbital elements for all minor bodies.

Provenance:
    Input position/velocity records are public NASA/JPL Horizons SPK type-21
    states. The J2000-equatorial-to-J2000-ecliptic rotation and
    state-to-osculating-element equations are stated in this script and
    referenced to Curtis and Bate-Mueller-White. Epoch spacing, selection, and
    runtime cross-fade policy are project choices. The script emits the body-map
    fragment used by the reviewed ``multi_epoch_data.py`` artifact; it does not
    fit foreign output.

The reviewed artifact in ``libephemeris/multi_epoch_data.py`` uses the default
ten-year grid from 1650 through 2450 inclusive. The generator attempts every
body in ``MINOR_BODY_ELEMENTS`` and emits only bodies for which a compatible
heliocentric J2000 type-21 kernel covers at least one requested epoch. The
reviewed artifact contains 36 of 37 bodies; Bennu is intentionally represented
only by its separately documented curated epoch in ``minor_bodies.py``.

The script:
1. Downloads/finds wide-range SPK files for each body
2. Computes heliocentric state vectors at the requested epochs (ten years by
   default)
3. Converts state vectors to osculating Keplerian elements
4. Outputs Python code for MINOR_BODY_ELEMENTS_MULTI

Usage:
    python scripts/generate_multi_epoch_elements.py                # All bodies with SPK
    python scripts/generate_multi_epoch_elements.py --body chiron  # Single body
    python scripts/generate_multi_epoch_elements.py --dry-run      # Show epochs only
    python scripts/generate_multi_epoch_elements.py --output /path/to/output.py
    python scripts/generate_multi_epoch_elements.py --spacing 5    # Optional denser grid
    python scripts/generate_multi_epoch_elements.py --start 1650 --end 2450

Requirements:
    pip install spktype21

Algorithm:
    State vector (x,y,z,vx,vy,vz) from SPK → osculating Keplerian elements
    using the vis-viva equation, angular-momentum and eccentricity vectors,
    and the elliptic/hyperbolic anomaly identities written beside the code.

    The accepted SPK segments declare NAIF frame ID 1 (J2000 equatorial) and
    center ID 10 (Sun). We convert to ecliptic J2000 and then extract (a, e, i,
    ω, Ω, M₀, n) in the same conventions as MINOR_BODY_ELEMENTS.

References:
    * NASA/JPL Horizons System manual, "SPK File Generation".
    * NASA/JPL/NAIF, ``SPK Required Reading``, type 21.
    * IAU (1976), System of Astronomical Constants; IAU (2012), Resolution B2.
    * Curtis, *Orbital Mechanics for Engineering Students*, chapter 4.
    * Bate, Mueller & White, *Fundamentals of Astrodynamics*, chapter 2.
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

# Historical Gaussian solar mass parameter in AU³/day².
#
# The IAU 1976 System of Astronomical Constants defined the Gaussian
# gravitational constant as k = 0.01720209895 when length, mass, and time are
# measured in au, solar masses, and days. The value below is k², rounded to
# the precision used when the reviewed multi_epoch_data.py table was produced.
# IAU 2012 Resolution B2 subsequently made the au an exact SI length and
# removed k from the current system of constants. We retain this historical
# value *only* as a declared generation convention so regenerating the shipped
# table does not silently change every mean motion. It must not be described
# as a current IAU nominal or observational GM of the Sun.
MU_SUN_AU3_DAY2 = 2.9591220828559093e-04

# Mean obliquity of the J2000 ecliptic used by the table's frame convention.
# The IAU 1976 resolution gives epsilon_0 = 23°26'21.448" =
# 23.439291111... degrees. The stored 23.4392911-degree value is its explicit
# 0.0000001-degree rounding. Position and velocity receive the identical
# passive x-axis rotation, preserving their common reference frame.
OBLIQUITY_J2000_RAD = math.radians(23.4392911)
COS_EPS = math.cos(OBLIQUITY_J2000_RAD)
SIN_EPS = math.sin(OBLIQUITY_J2000_RAD)

# Astronomical unit in kilometres. IAU 2012 Resolution B2 defines
# 1 au = 149,597,870,700 m exactly, hence this conversion is exact.
AU_KM = 149597870.7

# Reviewed artifact grid: Julian-year labels 1650–2450 at ten-year spacing.
# These labels are mapped to uniform 365.25-day offsets from J2000.0; they are
# not civil-calendar January 1 instants. See _year_to_jd().
DEFAULT_START_YEAR = 1650
DEFAULT_END_YEAR = 2450
DEFAULT_SPACING = 10


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
    """Convert a Julian-year label to its uniform-grid Julian Day.

    Uses the IAU/JPL epoch definition J2000.0 = JD 2451545.0 TDB
    (2000 January 1, 12:00 TDB) and 365.25 Julian days per Julian year. Thus
    labels other than 2000 are not assertions about civil-calendar January 1.
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
    """Get a heliocentric ecliptic-J2000 state from a type-21 SPK.

    Type-21 records store a target state relative to the segment's declared
    center and reference frame. This generator accepts only center ID 10
    (the Sun) and frame ID 1 (J2000). Enforcing those metadata is essential:
    merely calling a relative-state reader does not make a barycentric or
    planet-centered segment heliocentric, and rotating an arbitrary SPICE
    frame by the J2000 obliquity would be scientifically invalid.

    The type-21 reader returns kilometres and kilometres per second. The
    function converts them to au and au/day, then applies the documented
    equatorial-J2000 to ecliptic-J2000 x-axis rotation to both vectors.

    Returns:
        ``(x, y, z, vx, vy, vz)`` in ecliptic-J2000 au and au/day, or
        ``None`` when the kernel is incompatible, unreadable, or out of range.
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
            if not kernel.segments:
                return None

            center_id = kernel.segments[0].center
            target_id = kernel.segments[0].target

            # Horizons small-body kernels used for the reviewed artifact are
            # Sun-centered J2000 segments. Reject mixed or differently framed
            # files rather than silently assigning the wrong physical meaning
            # to their coordinates.
            if any(
                segment.center != center_id
                or segment.target != target_id
                or segment.frame != 1
                for segment in kernel.segments
            ):
                return None
            if center_id != 10:
                return None

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
    except Exception:
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
    """Generate the declared uniform Julian-year epoch grid.

    ``start_year`` and ``end_year`` are human-readable Julian-year labels,
    not Gregorian calendar dates. Each integer year advances exactly 365.25
    TDB days from J2000.0, matching the epoch convention embedded in the
    reviewed :data:`MINOR_BODY_ELEMENTS_MULTI` artifact. With the defaults,
    both endpoints are included and 81 epochs are emitted.

    Args:
        start_year: First inclusive Julian-year label.
        end_year: Last inclusive Julian-year label.
        spacing: Positive interval between labels in Julian years.

    Returns:
        Increasing Julian Dates on the declared TDB-like generation grid.

    Raises:
        ValueError: If spacing is not positive or the end precedes the start.
    """
    if spacing <= 0:
        raise ValueError("spacing must be positive")
    if end_year < start_year:
        raise ValueError("end_year must be greater than or equal to start_year")

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
  python scripts/generate_multi_epoch_elements.py --spacing 5    # Optional denser grid
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
        from libephemeris.spk import _get_body_name  # noqa: F401 (availability probe)
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
# Provenance: public NASA/JPL Horizons SPK states; explicit state-to-element
# equations in this script; no compatibility-output fit. This command emits the
# body-map fragment incorporated into the reviewed multi_epoch_data.py module.
# Generation conventions: IAU-1976 k² solar mass parameter and J2000 mean
# obliquity, plus the IAU-2012 exact au. The script rejects center != Sun and
# frame != J2000 so "heliocentric ecliptic J2000" is enforced, not assumed.
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
