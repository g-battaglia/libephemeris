#!/usr/bin/env python3
"""
Derive hypothetical body orbital elements from primary published sources.

This script independently computes L0, n, M0 values for the Hamburg School
Uranian planets and other hypothetical bodies using standard Keplerian
mechanics applied to the published orbital elements in fictitious_orbits.csv.

Primary sources:
    - Witte/Sieggrun, 'Regelwerk fur Planetenbilder' (1928), refined by Neely (1988)
    - Strubell, 'Die Sterne' 3/1952, p. 70ff (Transpluto/Isis)
    - L.H. Weston / Theosophical tradition (Vulcan)
    - Georg Waldemath, 1898 (Waldemath hypothetical second moon)

Usage:
    python scripts/derive_hypothetical_elements.py

Output:
    Prints derived values with comparison to current calibrated values.
"""

from __future__ import annotations

import os
import sys

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# Constants
J1900_JD = 2415020.0  # Julian Day of J1900.0
J2000_JD = 2451545.0  # Julian Day of J2000.0
DAYS_PER_CENTURY = 36525.0


def mean_motion_from_kepler(a_au: float) -> float:
    """Compute mean motion in deg/day from semi-major axis via Kepler's 3rd law.

    For heliocentric orbits:
        Period = a^1.5 years (with a in AU)
        n = 360 / (Period * 365.25) deg/day

    Args:
        a_au: Semi-major axis in AU.

    Returns:
        Mean motion in degrees per day.
    """
    period_years = a_au**1.5
    period_days = period_years * 365.25
    return 360.0 / period_days


def mean_longitude_at_epoch(
    mean_anomaly_deg: float,
    arg_perihelion_deg: float,
    asc_node_deg: float,
) -> float:
    """Compute mean longitude from mean anomaly and orbital angles.

    L = M + omega + Omega

    For circular orbits (e=0, omega=0, Omega=0), L = M.

    Args:
        mean_anomaly_deg: Mean anomaly at epoch (degrees).
        arg_perihelion_deg: Argument of perihelion (degrees).
        asc_node_deg: Longitude of ascending node (degrees).

    Returns:
        Mean longitude in degrees (0-360).
    """
    return (mean_anomaly_deg + arg_perihelion_deg + asc_node_deg) % 360.0


def propagate_mean_longitude(
    L0: float, n_deg_day: float, epoch_jd: float, target_jd: float
) -> float:
    """Propagate mean longitude from one epoch to another.

    L(t) = L0 + n * (t - t0)

    Args:
        L0: Mean longitude at epoch (degrees).
        n_deg_day: Mean motion (degrees per day).
        epoch_jd: Epoch Julian Day.
        target_jd: Target Julian Day.

    Returns:
        Mean longitude at target epoch (degrees, 0-360).
    """
    dt = target_jd - epoch_jd
    return (L0 + n_deg_day * dt) % 360.0


def load_csv_elements(csv_path: str) -> list[dict]:
    """Load orbital elements from fictitious_orbits.csv.

    Args:
        csv_path: Path to the CSV file.

    Returns:
        List of dicts with parsed orbital elements.
    """
    elements = []
    with open(csv_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 10:
                continue

            name = parts[0]
            try:
                epoch_jd = float(parts[1])
            except ValueError:
                continue

            # Parse mean_anomaly (may contain T-polynomial)
            ma_str = parts[3].strip()
            if "+" in ma_str and "T" in ma_str:
                # T-polynomial: skip for now
                mean_anomaly = None
            else:
                try:
                    mean_anomaly = float(ma_str)
                except ValueError:
                    mean_anomaly = None

            try:
                semi_axis = float(parts[4])
                eccentricity = float(parts[5])
                arg_perihelion = float(parts[6])
                asc_node = float(parts[7])
                inclination = float(parts[8])
                geocentric = int(parts[9])
            except (ValueError, IndexError):
                continue

            source = parts[10] if len(parts) > 10 else ""

            elements.append(
                {
                    "name": name,
                    "epoch_jd": epoch_jd,
                    "mean_anomaly": mean_anomaly,
                    "semi_axis": semi_axis,
                    "eccentricity": eccentricity,
                    "arg_perihelion": arg_perihelion,
                    "asc_node": asc_node,
                    "inclination": inclination,
                    "geocentric": geocentric,
                    "source": source,
                }
            )

    return elements


def main() -> None:
    """Derive and compare hypothetical body orbital elements."""
    csv_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "libephemeris",
        "data",
        "fictitious_orbits.csv",
    )

    if not os.path.exists(csv_path):
        print(f"ERROR: CSV file not found at {csv_path}")
        sys.exit(1)

    elements = load_csv_elements(csv_path)

    # Hamburg School Uranian planets
    uranian_names = [
        "Cupido",
        "Hades",
        "Zeus",
        "Kronos",
        "Apollon",
        "Admetos",
        "Vulkanus",
        "Poseidon",
    ]

    print("=" * 80)
    print("DERIVATION OF HYPOTHETICAL BODY ELEMENTS FROM PRIMARY SOURCES")
    print("=" * 80)
    print()
    print("Source: fictitious_orbits.csv (compiled from Witte/Sieggrun 1928,")
    print("        refined by Neely 1988)")
    print()

    # Current calibrated values for comparison
    current_values = {
        "Cupido": {"L0": 105.301693, "n": 0.0037945179},
        "Hades": {"L0": 336.363662, "n": 0.0027875901},
        "Zeus": {"L0": 104.289095, "n": 0.0022203750},
        "Kronos": {"L0": 17.111353, "n": 0.0019351856},
        "Apollon": {"L0": 138.565328, "n": 0.0017177599},
        "Admetos": {"L0": 350.613913, "n": 0.0016016766},
        "Vulkanus": {"L0": 55.397715, "n": 0.0015069325},
        "Poseidon": {"L0": 166.140256, "n": 0.0013256078},
    }

    print(
        f"{'Body':<12} {'n_kepler':>12} {'n_current':>12} {'delta_n':>12} "
        f"{'L0_derived':>12} {'L0_current':>12} {'delta_L0':>10}"
    )
    print("-" * 80)

    for name in uranian_names:
        elem = next((e for e in elements if e["name"] == name), None)
        if elem is None:
            print(f"{name:<12} NOT FOUND in CSV")
            continue

        # Derive mean motion from Kepler's 3rd law
        n_kepler = mean_motion_from_kepler(elem["semi_axis"])

        # Derive mean longitude at J1900
        if elem["mean_anomaly"] is not None:
            L0_j1900 = mean_longitude_at_epoch(
                elem["mean_anomaly"],
                elem["arg_perihelion"],
                elem["asc_node"],
            )
        else:
            L0_j1900 = 0.0

        cur = current_values.get(name, {})
        n_cur = cur.get("n", 0.0)
        L0_cur = cur.get("L0", 0.0)

        delta_n = (n_kepler - n_cur) * 1e6  # micro-deg/day
        delta_L0 = ((L0_j1900 - L0_cur + 180) % 360) - 180

        print(
            f"{name:<12} {n_kepler:>12.10f} {n_cur:>12.10f} {delta_n:>+10.2f}μ "
            f"{L0_j1900:>12.6f} {L0_cur:>12.6f} {delta_L0:>+10.4f}°"
        )

    print()
    print("Note: 'n_kepler' is derived purely from Kepler's 3rd law (a^1.5).")
    print("      'delta_n' is in micro-degrees/day (1e-6 deg/day).")
    print("      'delta_L0' is derived_L0 - current_L0 in degrees.")
    print()
    print("The calibrated values in hypothetical.py were fitted for exact")
    print("compatibility with the reference implementation. The Keplerian")
    print("derivation provides independent verification from primary sources.")

    # Special case: Hades discrepancy
    print()
    print("=" * 80)
    print("HADES DISCREPANCY ANALYSIS")
    print("=" * 80)
    hades = next((e for e in elements if e["name"] == "Hades"), None)
    if hades:
        print(f"  CSV mean_anomaly at J1900: {hades['mean_anomaly']:.4f} deg")
        L0_hades = mean_longitude_at_epoch(
            hades["mean_anomaly"],
            hades["arg_perihelion"],
            hades["asc_node"],
        )
        print(f"  Derived mean longitude (M + omega + Omega): {L0_hades:.6f} deg")
        print()
        print("  HADES_KEPLERIAN_ELEMENTS:")
        print("    M0 = 26.850162 (mean anomaly semantics)")
        print("    n  = 0.00278759")
        print()
        print("  URANIAN_KEPLERIAN_ELEMENTS['Hades']:")
        print("    M0 = 336.363662 (mean longitude semantics)")
        print("    n  = 0.0027875901")
        print()
        print("  Resolution: The standalone uses mean anomaly (M), while the")
        print("  unified dict uses mean longitude (L = M + omega + Omega).")
        print("  These are different but self-consistent within each usage context.")
        print(
            f"  Difference in n: {abs(0.00278759 - 0.0027875901) * 1e6:.2f} micro-deg/day"
        )
        print("  (negligible over astrological timeframes)")

    # Transpluto
    print()
    print("=" * 80)
    print("TRANSPLUTO/ISIS DERIVATION")
    print("=" * 80)
    transpluto = next(
        (e for e in elements if "Transpluto" in e["name"] or "Isis" in e["name"]),
        None,
    )
    if transpluto:
        print("  Source: Strubell, 'Die Sterne' 3/1952, p. 70ff")
        print(f"  CSV epoch: JD {transpluto['epoch_jd']}")
        print(f"  a = {transpluto['semi_axis']} AU, e = {transpluto['eccentricity']}")
        n_tp = mean_motion_from_kepler(transpluto["semi_axis"])
        print(f"  Kepler n = {n_tp:.10f} deg/day")
        print("  Current n = 0.0011968259 deg/day")
        print(f"  Delta: {abs(n_tp - 0.0011968259) * 1e6:.2f} micro-deg/day")


if __name__ == "__main__":
    main()
