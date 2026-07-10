# Lunar Apsides: Computational Methodology

LibEphemeris computes the interpolated lunar apsides (perigee and apogee) using a passage-interpolated harmonic fitting method: passage geometry from the JPL DE440/DE441 numerical integrations, with the fitted coefficients and residual table calibrated against reference-API output used as a black-box oracle (NOTICE.md, Calibration Data Disclosure) — rather than the classic analytical term-selection approach.

## Background

The lunar apsides — perigee (closest approach) and apogee (farthest point) — are among the most computationally challenging quantities in positional astronomy. The Moon's orbit is strongly perturbed by the Sun, causing the instantaneous (osculating) perigee to oscillate by approximately ±30 degrees over a single anomalistic month.

For practical ephemeris use, this oscillation must be smoothed to produce an "interpolated" or "natural" apsidal position that reflects the genuine long-term motion of the apsidal line without the spurious short-period volatility inherent in the two-body approximation.

The choice of smoothing methodology is the single most consequential design decision for the interpolated apsides, and the place where ephemeris engines differ most.

## Method

LibEphemeris constructs the interpolated apsides from the physical geometry of the JPL DE440/DE441 numerical integrations:

1. **Passage identification.** All perigee passages (local Earth-Moon distance minima) are identified from JPL state vectors over a 1000-year calibration span (1500–2500 CE). At each passage, the Moon's ecliptic longitude is an unambiguous physical measurement of the perigee direction. Over 12,000 passages are used.

2. **Spline interpolation.** A cubic spline is fitted through the passage longitudes (with angle unwrapping) to produce a smooth, continuous perigee longitude function at arbitrary times.

3. **Harmonic series calibration.** A 66-term trigonometric perturbation series, constructed from the standard Delaunay arguments (D, M, M', F), is fitted to the spline via least squares. Terms with amplitudes below 0.001 degrees are discarded.

4. **Residual correction.** A precomputed correction table (`PERIGEE_CORRECTIONS`, 201,249 entries at a 2-day step over 1549–2650) absorbs the remaining difference between the harmonic model and the calibration target (reference-API output; black-box oracle, disclosed in NOTICE.md).

The result is a smooth apsidal curve anchored to the physical distance extrema of the Moon as computed by modern numerical integration.

## Precision and Validation

### Smoothing characteristics

The passage-interpolated perigee follows the full ~25° oscillation amplitude of the
physical apsidal line; an analytical term-selection method that removes
mean-anomaly terms yields a smaller ~15° amplitude. For the apogee, where
perturbation amplitudes are smaller, the two philosophies converge.

| Property                      | LibEphemeris                            |
| ----------------------------- | --------------------------------------- |
| Passage geometry              | JPL DE440/DE441 numerical integration   |
| Calibration target            | Reference-API output (black-box oracle) |
| Smoothing method              | Physical passage interpolation          |
| Perigee oscillation amplitude | ~25 deg from mean                       |
| Apogee oscillation amplitude  | ~5 deg from mean                        |
| Date range                    | 1550–2650 (DE440) / -13200 to +17191 (DE441) |

Because the two approaches smooth the same phenomenon differently, the interpolated
perigee (`INTP_PERG`) is where libephemeris differs most from an analytical engine
(up to ~5°), while the interpolated apogee (`INTP_APOG`) differs by ~0.36° at most.
These head-to-head figures are catalogued in the
[Swiss Ephemeris Comparison](../comparison/known-differences.md).

### Rationale

LibEphemeris adopts JPL numerical integrations as the primary reference for orbital geometry. The DE440/DE441 ephemerides incorporate lunar laser ranging data accurate to approximately 1 milliarcsecond and represent the current standard for planetary and lunar ephemeris computation (Park et al., 2021).

The ELP2000-82B theory, while a significant achievement of 20th-century celestial mechanics, is a truncated analytical approximation fitted to an earlier generation of observations. Where the analytical smoothing and the physical passage interpolation disagree, the JPL-grounded approach more closely represents the actual state of the Earth-Moon system.

This choice prioritizes physical accuracy over backward compatibility with the analytical framework.

## References

1. Park, R.S. et al. (2021). "The JPL Planetary and Lunar Ephemerides DE440 and DE441." *Astronomical Journal*, 161(3), 105.
2. Chapront-Touzé, M. & Chapront, J. (1988). "ELP 2000-82B: A semi-analytical lunar ephemeris." *Astronomy & Astrophysics*, 190, 342-352.
3. Moshier, S.L. (1992). "Comparison of a 7000-year lunar ephemeris with analytical theory." *Astronomy & Astrophysics*, 262, 613-616.
4. Meeus, J. (1998). *Astronomical Algorithms*, 2nd edition. Willmann-Bell.

See also: [interpolated-perigee.md](interpolated-perigee.md) (calibration details),
[interpolated-apogee.md](interpolated-apogee.md) (apogee-specific methodology).
