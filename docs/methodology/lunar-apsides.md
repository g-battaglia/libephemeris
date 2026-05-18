# Lunar Apsides: Computational Methodology

LibEphemeris computes the interpolated lunar apsides (perigee and apogee) using a passage-interpolated harmonic fitting method anchored to JPL DE440/DE441 numerical integrations, differing fundamentally from the analytical term-selection approach used by Swiss Ephemeris.

## Background

The lunar apsides — perigee (closest approach) and apogee (farthest point) — are among the most computationally challenging quantities in positional astronomy. The Moon's orbit is strongly perturbed by the Sun, causing the instantaneous (osculating) perigee to oscillate by approximately ±30 degrees over a single anomalistic month.

For practical ephemeris use, this oscillation must be smoothed to produce an "interpolated" or "natural" apsidal position that reflects the genuine long-term motion of the apsidal line without the spurious short-period volatility inherent in the two-body approximation.

The choice of smoothing methodology constitutes the most significant computational difference between LibEphemeris and Swiss Ephemeris.

## Method

LibEphemeris constructs the interpolated apsides from the physical geometry of the JPL DE440/DE441 numerical integrations:

1. **Passage identification.** All perigee passages (local Earth-Moon distance minima) are identified from JPL state vectors over a 1000-year calibration span (1500–2500 CE). At each passage, the Moon's ecliptic longitude is an unambiguous physical measurement of the perigee direction. Over 12,000 passages are used.

2. **Spline interpolation.** A cubic spline is fitted through the passage longitudes (with angle unwrapping) to produce a smooth, continuous perigee longitude function at arbitrary times.

3. **Harmonic series calibration.** A 61-term trigonometric perturbation series, constructed from the standard Delaunay arguments (D, M, M', F), is fitted to the spline via least squares. Terms with amplitudes below 0.001 degrees are discarded.

4. **Residual correction.** A precomputed correction table (~15,000 entries) absorbs the remaining difference between the harmonic model and the JPL ground truth.

The result is a smooth apsidal curve anchored to the physical distance extrema of the Moon as computed by modern numerical integration.

## Precision and Validation

### Measured Discrepancy

The interpolated perigee (`INTP_PERG`) in LibEphemeris differs from Swiss Ephemeris by up to approximately 5 degrees. This is the largest single discrepancy between the two libraries.

The difference arises from two distinct smoothing philosophies applied to the same underlying phenomenon:

| Property                      | Swiss Ephemeris                 | LibEphemeris                            |
| ----------------------------- | ------------------------------- | --------------------------------------- |
| Ground truth                  | ELP2000-82B analytical theory   | JPL DE440/DE441 numerical integration   |
| Smoothing method              | Analytical term selection       | Physical passage interpolation          |
| Perigee oscillation amplitude | ~15 deg from mean               | ~25 deg from mean                       |
| Apogee oscillation amplitude  | ~5 deg from mean                | ~5 deg from mean                        |
| Date range                    | ~-5400 to +5400 CE              | 1550–2650 (DE440) / -13200 to +17191 (DE441) |

The interpolated apogee (`INTP_APOG`) shows a smaller discrepancy (~0.36 degrees maximum), as both approaches produce similar results for the apogee where perturbation amplitudes are smaller.

### Rationale

LibEphemeris adopts JPL numerical integrations as the primary reference for orbital geometry. The DE440/DE441 ephemerides incorporate lunar laser ranging data accurate to approximately 1 milliarcsecond and represent the current standard for planetary and lunar ephemeris computation (Park et al., 2021).

The ELP2000-82B theory, while a significant achievement of 20th-century celestial mechanics, is a truncated analytical approximation fitted to an earlier generation of observations. Where the analytical smoothing and the physical passage interpolation disagree, the JPL-grounded approach more closely represents the actual state of the Earth-Moon system.

This choice prioritizes physical accuracy over backward compatibility with the analytical framework.

## Comparison with Swiss Ephemeris

Swiss Ephemeris computes the interpolated apsides using the analytical method developed by S.L. Moshier, based on the ELP2000-82B lunar theory (Chapront-Touzé & Chapront, 1988).

This approach works within the analytical framework of the lunar theory itself: the thousands of trigonometric terms in ELP2000-82B are classified by their physical origin, and terms associated with the mean anomaly (the Moon's monthly orbital cycle) are excluded. The remaining terms define the smoothed apsidal position.

This produces a mathematically coherent result within its theoretical framework. However, the output is constrained by the truncation level and fitting epoch of the analytical theory (1988). The resulting curve can differ from the physical geometry of the Earth-Moon system — as represented by modern numerical integrations — by several degrees.

## References

1. Park, R.S. et al. (2021). "The JPL Planetary and Lunar Ephemerides DE440 and DE441." *Astronomical Journal*, 161(3), 105.
2. Chapront-Touzé, M. & Chapront, J. (1988). "ELP 2000-82B: A semi-analytical lunar ephemeris." *Astronomy & Astrophysics*, 190, 342-352.
3. Moshier, S.L. (1992). "Comparison of a 7000-year lunar ephemeris with analytical theory." *Astronomy & Astrophysics*, 262, 613-616.
4. Meeus, J. (1998). *Astronomical Algorithms*, 2nd edition. Willmann-Bell.

See also: [interpolated-perigee.md](interpolated-perigee.md) (calibration details),
[interpolated-apogee.md](interpolated-apogee.md) (apogee-specific methodology).
