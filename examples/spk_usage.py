#!/usr/bin/env python3
"""
Example: Using SPK kernels for high-precision minor body calculations.

This script demonstrates how to:
1. Download SPK kernels from JPL Horizons
2. Register bodies for SPK-based calculations
3. Compare SPK precision vs Keplerian approximation
4. Use SPK with EphemerisContext for thread-safe calculations

Requirements:
    - Internet connection (for downloading SPK files)
    - libephemeris installed

Usage:
    python examples/spk_usage.py
"""

import os
import sys
import tempfile

# Add parent directory to path for development
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import libephemeris as eph


def example_basic_download():
    """Example 1: Basic SPK download and registration."""
    print("=" * 60)
    print("Example 1: Basic SPK Download and Registration")
    print("=" * 60)

    # Create temporary directory for SPK files
    with tempfile.TemporaryDirectory() as tmpdir:
        print(f"\nDownloading Chiron SPK to {tmpdir}...")

        try:
            # Download SPK for Chiron (asteroid 2060)
            # Date range: 2020-2030 (10 years is enough for demo)
            path = eph.download_spk(
                body="2060",  # Chiron's asteroid number
                start="2020-01-01",
                end="2030-01-01",
                path=tmpdir,
            )
            print(f"Downloaded: {os.path.basename(path)}")
            print(f"File size: {os.path.getsize(path) / 1024:.1f} KB")

            # Register the SPK for Chiron
            eph.register_spk_body(
                ipl=eph.CHIRON,
                spk_file=path,
                naif_id=eph.NAIF_CHIRON,  # 2002060
            )
            print(f"Registered CHIRON with NAIF ID {eph.NAIF_CHIRON}")

            # Check registration
            info = eph.get_spk_body_info(eph.CHIRON)
            print(f"Registration info: {info}")

            # Calculate position using SPK
            jd = 2459215.5  # 2021-01-01
            pos, _ = eph.calc_ut(jd, eph.CHIRON, eph.FLG_SPEED)
            print(f"\nChiron position (SPK) at JD {jd}:")
            print(f"  Longitude: {pos[0]:.6f} deg")
            print(f"  Latitude:  {pos[1]:.6f} deg")
            print(f"  Distance:  {pos[2]:.6f} AU")
            print(f"  Speed:     {pos[3]:.6f} deg/day")

            # Unregister to go back to Keplerian
            eph.unregister_spk_body(eph.CHIRON)
            print("\nUnregistered SPK, back to Keplerian model")

        except Exception as e:
            print(f"Error: {e}")
            print("(This example requires internet connection)")


def example_compare_precision():
    """Example 2: Compare SPK vs Keplerian precision."""
    print("\n" + "=" * 60)
    print("Example 2: Compare SPK vs Keplerian Precision")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            # First, calculate with Keplerian (default)
            jd = 2459215.5  # 2021-01-01
            pos_kep, _ = eph.calc_ut(jd, eph.CHIRON, 0)
            print(f"\nChiron at JD {jd}:")
            print(f"  Keplerian: {pos_kep[0]:.4f} deg")

            # Download and register SPK
            path = eph.download_and_register_spk(
                body="2060",
                ipl=eph.CHIRON,
                start="2020-01-01",
                end="2030-01-01",
                path=tmpdir,
            )

            # Calculate with SPK
            pos_spk, _ = eph.calc_ut(jd, eph.CHIRON, 0)
            print(f"  SPK:       {pos_spk[0]:.4f} deg")

            # Compare
            diff_arcsec = abs(pos_spk[0] - pos_kep[0]) * 3600
            print(f"\n  Difference: {diff_arcsec:.1f} arcseconds")
            print(f"              ({diff_arcsec / 60:.2f} arcminutes)")

            if diff_arcsec > 60:
                print("  -> SPK provides significantly higher precision!")
            else:
                print(
                    "  -> Results are close (Keplerian model is accurate for this date)"
                )

            # Cleanup
            eph.unregister_spk_body(eph.CHIRON)

        except Exception as e:
            print(f"Error: {e}")
            print("(This example requires internet connection)")


def example_multiple_bodies():
    """Example 3: Register multiple minor bodies."""
    print("\n" + "=" * 60)
    print("Example 3: Multiple Minor Bodies")
    print("=" * 60)

    bodies = [
        ("1", eph.CERES, eph.NAIF_CERES, "Ceres"),
        ("4", eph.VESTA, eph.NAIF_VESTA, "Vesta"),
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            print(f"\nDownloading SPK files for {len(bodies)} bodies...")

            for number, ipl, naif_id, name in bodies:
                path = eph.download_and_register_spk(
                    body=number,
                    ipl=ipl,
                    start="2020-01-01",
                    end="2030-01-01",
                    path=tmpdir,
                )
                print(f"  {name}: registered")

            # List all registered bodies
            print(f"\nRegistered bodies: {eph.list_spk_bodies()}")

            # Calculate positions
            jd = 2459215.5
            print(f"\nPositions at JD {jd}:")
            for number, ipl, naif_id, name in bodies:
                pos, _ = eph.calc_ut(jd, ipl, 0)
                print(f"  {name}: {pos[0]:.4f} deg")

            # Cleanup
            for _, ipl, _, _ in bodies:
                eph.unregister_spk_body(ipl)

        except Exception as e:
            print(f"Error: {e}")


def example_ephemeris_context():
    """Example 4: Using SPK with EphemerisContext (thread-safe)."""
    print("\n" + "=" * 60)
    print("Example 4: SPK with EphemerisContext (Thread-Safe)")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            # Download SPK file first
            path = eph.download_spk(
                body="2060",
                start="2020-01-01",
                end="2030-01-01",
                path=tmpdir,
            )

            # Create two contexts with different SPK registrations
            ctx1 = eph.EphemerisContext()
            ctx2 = eph.EphemerisContext()

            # Register SPK only in ctx1
            ctx1.register_spk_body(eph.CHIRON, path, eph.NAIF_CHIRON)

            print("\nContext 1: SPK registered")
            print(f"  list_spk_bodies(): {ctx1.list_spk_bodies()}")

            print("\nContext 2: No SPK (uses Keplerian)")
            print(f"  list_spk_bodies(): {ctx2.list_spk_bodies()}")

            # Calculate in both contexts
            jd = 2459215.5
            pos1, _ = ctx1.calc_ut(jd, eph.CHIRON, 0)
            pos2, _ = ctx2.calc_ut(jd, eph.CHIRON, 0)

            print(f"\nChiron at JD {jd}:")
            print(f"  Context 1 (SPK):      {pos1[0]:.4f} deg")
            print(f"  Context 2 (Keplerian): {pos2[0]:.4f} deg")

            diff = abs(pos1[0] - pos2[0]) * 3600
            print(f"  Difference: {diff:.1f} arcseconds")

        except Exception as e:
            print(f"Error: {e}")


def example_coverage_check():
    """Example 5: Check SPK coverage before calculation."""
    print("\n" + "=" * 60)
    print("Example 5: SPK Coverage Check")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            # Download SPK with limited date range
            path = eph.download_spk(
                body="2060",
                start="2020-01-01",
                end="2025-01-01",  # Only 5 years
                path=tmpdir,
            )

            # Check coverage
            coverage = eph.get_spk_coverage(path)
            if coverage:
                start_jd, end_jd = coverage
                print("\nSPK coverage:")
                print(f"  Start: JD {start_jd:.1f}")
                print(f"  End:   JD {end_jd:.1f}")

                # Register
                eph.register_spk_body(eph.CHIRON, path, eph.NAIF_CHIRON)

                # Calculate within coverage
                jd_ok = 2459215.5  # 2021-01-01 (within range)
                pos, _ = eph.calc_ut(jd_ok, eph.CHIRON, 0)
                print(f"\nCalculation at JD {jd_ok} (within coverage): OK")
                print(f"  Longitude: {pos[0]:.4f} deg")

                # Calculate outside coverage
                jd_bad = 2465000.0  # ~2029 (outside range)
                try:
                    pos, _ = eph.calc_ut(jd_bad, eph.CHIRON, 0)
                    print(f"\nCalculation at JD {jd_bad} (outside coverage):")
                    print("  Should have raised error!")
                except ValueError as e:
                    print(f"\nCalculation at JD {jd_bad} (outside coverage):")
                    print(f"  Correctly raised error: {e}")

                eph.unregister_spk_body(eph.CHIRON)

        except Exception as e:
            print(f"Error: {e}")


def main():
    """Run all examples."""
    print("\n" + "#" * 60)
    print("# libephemeris SPK Examples")
    print("#" * 60)

    # Check if we should skip network examples
    skip_network = os.environ.get("SKIP_NETWORK_EXAMPLES", "0") == "1"

    if skip_network:
        print("\nNetwork examples skipped (SKIP_NETWORK_EXAMPLES=1)")
        print("Showing offline functionality only...\n")

        # Show basic offline functionality
        print("Current SPK registrations:", eph.list_spk_bodies())
        print("NAIF_CHIRON constant:", eph.NAIF_CHIRON)

        # Keplerian fallback always works
        pos, _ = eph.calc_ut(2451545.0, eph.CHIRON, 0)
        print(f"Chiron (Keplerian): {pos[0]:.4f} deg")

    else:
        # Run examples that require network
        example_basic_download()
        example_compare_precision()
        example_multiple_bodies()
        example_ephemeris_context()
        example_coverage_check()

    print("\n" + "#" * 60)
    print("# Examples completed!")
    print("#" * 60 + "\n")


if __name__ == "__main__":
    main()
