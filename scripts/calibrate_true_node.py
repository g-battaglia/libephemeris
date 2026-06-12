#!/usr/bin/env python3
"""
True Lunar Node Calibration Script.

This script calibrates the empirical precession correction coefficient
in libephemeris to better match Swiss Ephemeris True Node calculations.

The goal is to find the optimal value for `precession_correction` in
lunar.py:1416 that minimizes the difference between libephemeris and
Swiss Ephemeris.

Usage:
    python scripts/calibrate_true_node.py
"""

import random
import numpy as np
from typing import Tuple, List
import swisseph as swe

# We need to temporarily modify libephemeris to test different coefficients
# So we'll compute the correction externally

# Import libephemeris components we need
from libephemeris.lunar import calc_true_lunar_node
from libephemeris import julday


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


def get_swe_true_node(jd: float) -> float:
    """Get True Node from Swiss Ephemeris."""
    pos, _ = swe.calc_ut(jd, swe.TRUE_NODE, 0)
    return pos[0]


def get_libephem_true_node(jd: float) -> float:
    """Get True Node from libephemeris."""
    lon, lat, dist = calc_true_lunar_node(jd)
    return lon


def generate_test_dates(n: int = 500, seed: int = 42) -> List[float]:
    """Generate random Julian Days in the valid range (1900-2100)."""
    random.seed(seed)
    dates = []

    for _ in range(n):
        year = random.randint(1900, 2100)
        month = random.randint(1, 12)
        day = random.randint(1, 28)
        hour = random.uniform(0, 24)
        jd = julday(year, month, day, hour)
        dates.append(jd)

    return sorted(dates)


def compute_errors(dates: List[float]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute errors between libephemeris and Swiss Ephemeris.

    Returns:
        Tuple of (errors_arcsec, T_values, jd_values)
    """
    errors = []
    T_values = []
    jd_values = []

    J2000 = 2451545.0

    for jd in dates:
        swe_node = get_swe_true_node(jd)
        lib_node = get_libephem_true_node(jd)

        # Signed difference (libephemeris - swisseph)
        diff = lib_node - swe_node
        if diff > 180:
            diff -= 360
        elif diff < -180:
            diff += 360

        # Convert to arcsec
        diff_arcsec = diff * 3600

        # T in Julian centuries from J2000
        T = (jd - J2000) / 36525.0

        errors.append(diff_arcsec)
        T_values.append(T)
        jd_values.append(jd)

    return np.array(errors), np.array(T_values), np.array(jd_values)


def fit_linear_correction(
    T_values: np.ndarray, errors: np.ndarray
) -> Tuple[float, float]:
    """
    Fit a linear model: error = offset + slope * T

    Returns:
        Tuple of (offset, slope) in arcsec
    """
    # Use numpy polyfit for linear regression
    slope, offset = np.polyfit(T_values, errors, 1)
    return offset, slope


def fit_quadratic_correction(
    T_values: np.ndarray, errors: np.ndarray
) -> Tuple[float, float, float]:
    """
    Fit a quadratic model: error = c0 + c1*T + c2*T^2

    Returns:
        Tuple of (c0, c1, c2) in arcsec
    """
    coeffs = np.polyfit(T_values, errors, 2)
    return coeffs[2], coeffs[1], coeffs[0]  # c0, c1, c2


def main():
    print("=" * 70)
    print("TRUE LUNAR NODE CALIBRATION")
    print("=" * 70)
    print()

    # Generate test dates
    print("Generating 500 random test dates (1900-2100)...")
    dates = generate_test_dates(500)
    print(f"  Date range: JD {dates[0]:.1f} to {dates[-1]:.1f}")
    print()

    # Compute current errors
    print("Computing errors (libephemeris vs Swiss Ephemeris)...")
    errors, T_values, jd_values = compute_errors(dates)
    print()

    # Current statistics
    print("CURRENT ERROR STATISTICS (before calibration):")
    print("-" * 50)
    print(
        f"  Mean error:   {np.mean(errors):+.2f} arcsec ({np.mean(errors) / 3600:+.6f}°)"
    )
    print(f"  Std dev:      {np.std(errors):.2f} arcsec")
    print(f"  Min error:    {np.min(errors):+.2f} arcsec")
    print(f"  Max error:    {np.max(errors):+.2f} arcsec")
    print(f"  RMS error:    {np.sqrt(np.mean(errors**2)):.2f} arcsec")
    print()

    # Fit linear correction
    print("FITTING LINEAR CORRECTION: error = offset + slope * T")
    print("-" * 50)
    offset, slope = fit_linear_correction(T_values, errors)
    print(f"  Offset (c0): {offset:+.2f} arcsec ({offset / 3600:+.8f}°)")
    print(f"  Slope (c1):  {slope:+.2f} arcsec/century ({slope / 3600:+.8f}°/century)")
    print()

    # Calculate residuals after linear correction
    linear_corrected = errors - (offset + slope * T_values)
    print("RESIDUALS AFTER LINEAR CORRECTION:")
    print("-" * 50)
    print(f"  Mean residual: {np.mean(linear_corrected):+.2f} arcsec")
    print(f"  Std dev:       {np.std(linear_corrected):.2f} arcsec")
    print(f"  RMS residual:  {np.sqrt(np.mean(linear_corrected**2)):.2f} arcsec")
    print()

    # Fit quadratic correction
    print("FITTING QUADRATIC CORRECTION: error = c0 + c1*T + c2*T^2")
    print("-" * 50)
    c0, c1, c2 = fit_quadratic_correction(T_values, errors)
    print(f"  c0 (offset):    {c0:+.2f} arcsec ({c0 / 3600:+.8f}°)")
    print(f"  c1 (linear):    {c1:+.2f} arcsec/century ({c1 / 3600:+.8f}°/century)")
    print(f"  c2 (quadratic): {c2:+.2f} arcsec/century² ({c2 / 3600:+.8f}°/century²)")
    print()

    # Calculate residuals after quadratic correction
    quad_corrected = errors - (c0 + c1 * T_values + c2 * T_values**2)
    print("RESIDUALS AFTER QUADRATIC CORRECTION:")
    print("-" * 50)
    print(f"  Mean residual: {np.mean(quad_corrected):+.2f} arcsec")
    print(f"  Std dev:       {np.std(quad_corrected):.2f} arcsec")
    print(f"  RMS residual:  {np.sqrt(np.mean(quad_corrected**2)):.2f} arcsec")
    print()

    # Current coefficient in lunar.py
    print("CURRENT IMPLEMENTATION (lunar.py:1416):")
    print("-" * 50)
    print("  precession_correction = 0.00006 * T")
    print(f"  This equals: {0.00006 * 3600:.2f} arcsec/century")
    print()

    # Recommended fix
    print("=" * 70)
    print("RECOMMENDED FIX")
    print("=" * 70)
    print()

    # Convert correction to degrees (what we need to ADD to libephemeris to match swe)
    # If error is positive (lib > swe), we need to subtract
    # So correction = -error

    correction_offset_deg = -offset / 3600
    correction_slope_deg = -slope / 3600

    print("To match Swiss Ephemeris, modify lunar.py line 1416:")
    print()
    print("  # BEFORE:")
    print("  precession_correction = 0.00006 * T")
    print()
    print("  # AFTER (linear correction):")
    print(
        f"  precession_correction = {correction_offset_deg:.8f} + {correction_slope_deg:.8f} * T"
    )
    print()

    # Or just adjust the coefficient if offset is small
    if abs(offset) < 10:  # If offset < 10 arcsec, just adjust coefficient
        current_coeff = 0.00006
        new_coeff = (
            current_coeff + correction_slope_deg / 100
        )  # per century -> per century
        # Actually, the current coefficient is in degrees, and T is in centuries
        # So: current gives 0.00006 deg/century = 0.216 arcsec/century
        # We need to add correction_slope_deg per century
        new_coeff = current_coeff - (slope / 3600)

        print("  # SIMPLIFIED (adjust coefficient only):")
        print(f"  precession_correction = {new_coeff:.8f} * T")
        print()

    print("=" * 70)
    print("EXPECTED IMPROVEMENT")
    print("=" * 70)
    print()
    print(f"  Current RMS error:  {np.sqrt(np.mean(errors**2)):.2f} arcsec")
    print(f"  After linear fix:   {np.sqrt(np.mean(linear_corrected**2)):.2f} arcsec")
    print(f"  After quadratic:    {np.sqrt(np.mean(quad_corrected**2)):.2f} arcsec")
    print()
    print(
        f"  Improvement factor: {np.sqrt(np.mean(errors**2)) / np.sqrt(np.mean(linear_corrected**2)):.1f}x"
    )
    print()

    return {
        "offset": offset,
        "slope": slope,
        "c0": c0,
        "c1": c1,
        "c2": c2,
        "current_rms": np.sqrt(np.mean(errors**2)),
        "linear_rms": np.sqrt(np.mean(linear_corrected**2)),
        "quad_rms": np.sqrt(np.mean(quad_corrected**2)),
    }


if __name__ == "__main__":
    results = main()
