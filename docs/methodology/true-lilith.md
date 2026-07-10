# True Lilith (Osculating Lunar Apogee) Calculation

LibEphemeris computes the osculating lunar apogee (True Lilith) directly from the Moon's geocentric state vectors obtained via JPL DE440/DE441: the eccentricity vector of the instantaneous two-body orbit points toward perigee, and the apogee is taken 180° opposite. No analytical perturbation series is applied — the perturbations are already contained in the numerical ephemeris state vectors.

## Table of Contents

- [Background](#background)
- [Method](#method)
  - [Eccentricity Vector Method](#1-eccentricity-vector-method)
  - [Orbital Elements Method](#2-orbital-elements-method)
- [Precision and Validation](#precision-and-validation)
  - [Measured Precision](#measured-precision)
  - [Lilith Method Selection Guide](#lilith-method-selection-guide)
- [The Osculating Apogee Concept](#the-osculating-apogee-concept)
- [The Interpolated Apogee Alternative](#the-interpolated-apogee-alternative)
- [References](#references)

## Background

True Lilith, also known as the osculating lunar apogee, is the apogee point of the Moon's instantaneous (osculating) Keplerian orbit. Unlike Mean Lilith (which moves smoothly along the mean lunar orbit), True Lilith oscillates significantly due to solar gravitational perturbations.

The standard method for computing the osculating apogee derives it directly from the Moon's **position and velocity vectors** (state vectors), from which the full set of osculating Keplerian orbital elements is obtained. This avoids the error introduced by analytical approximations and works directly with the numerical ephemeris. It is the standard osculating-apogee definition used by the `OSCU_APOG` body.

### Physical Context

The lunar orbit is far from a simple Keplerian ellipse because it is strongly perturbed by solar gravity (three-body problem). When the osculating ellipse is computed from instantaneous state vectors, the resulting apogee direction oscillates with the period of a synodic month as the Sun's tidal influence alternately stretches and compresses the apparent orbital ellipse.

Key characteristics of the osculating apogee:

1. **Two-Body Approximation**: The osculating apogee treats the Moon's orbit as a simple Keplerian ellipse (two-body problem: Earth and Moon only).

2. **Large Oscillations**: The osculating apogee oscillates **+/- 30 degrees** around the mean apogee. This is an artifact of the insufficient two-body model when applied to the strongly perturbed Earth-Moon-Sun system.

3. **Not a Real Motion**: The 30-degree oscillation does not represent actual motion of the lunar apsides. It is a mathematical consequence of forcing the three-body Earth-Moon-Sun system into a two-body framework.

4. **Instantaneous Validity**: The osculating apogee coincides with actual lunar maximum distance only twice per month — at conjunction and opposition with the Moon.

This is a well-understood artifact of the osculating-elements formalism, documented in standard celestial mechanics texts (Brouwer & Clemence 1961; Murray & Dermott 1999). The +/-30 degree monthly swing reflects the strength of solar perturbations, not real apsidal motion.

## Method

### 1. Eccentricity Vector Method

The implementation (`calc_true_lilith`) computes the apogee direction from the eccentricity vector:

```
Algorithm:
1. Get Moon's geocentric position (r) and velocity (v) from JPL DE ephemeris,
   directly in the true ecliptic frame of date (precession and nutation
   included by the frame)
2. Compute angular momentum: h = r x v
3. Compute eccentricity vector: e = (v x h)/mu - r/|r| (points toward perigee)
4. Apogee direction = -e (opposite to perigee)
5. Convert from Cartesian to spherical (ecliptic longitude, latitude);
   apogee distance = p/(1 - e) from the osculating conic
```

No analytical correction series (evection, variation, annual equation, etc.) is
applied on top of this: those perturbations are already present in the DE440/DE441
state vectors, and the osculating apogee is *defined* as the two-body apogee of
the instantaneous orbit. (Analytical perturbation series of that kind are used
elsewhere, for the mean-node/mean-apse ELP2000-style models — not here.)

### 2. Orbital Elements Method

`calc_true_lilith_orbital_elements` is retained for backward compatibility as an
alias of `calc_true_lilith`: deriving the apogee longitude from the classical
elements (Omega + omega + 180°) is mathematically equivalent to the eccentricity
vector method, so both entry points share the single implementation above.

## Precision and Validation

### Measured Precision

Precision is measured via 500-date random sampling for the osculating apogee
(`OSCU_APOG`), as measured differences:

| Method | Mean Difference | Max Difference | RMS Difference |
|--------|-----------------|----------------|----------------|
| **True Lilith** (LibEphemeris) | ~52 arcsec (~0.015 deg) | ~235 arcsec (~0.065 deg) | ~60 arcsec (~0.017 deg) |
| **Mean Lilith** | ~12 arcsec (~0.003 deg) | ~18 arcsec (~0.005 deg) | smooth / stable |
| Interpolated Apogee | ~6 arcsec RMS (within correction table, 1549–2651) | ~40 arcsec | see [interpolated-apogee.md](interpolated-apogee.md) |

For further precision benchmarks across all bodies, see [Precision Reference](../reference/precision.md).

### Lilith Method Selection Guide

| Use Case | Recommended Method | Precision (measured) |
|----------|-------------------|----------------------|
| Osculating lunar apogee (`OSCU_APOG`) | **True Lilith** | ~0.015 deg mean |
| Smooth, predictable motion (`MEAN_APOG`) | **Mean Lilith** | ~0.003 deg mean |
| Physical apogee passages (`INTP_APOG`) | Interpolated Apogee | ~6" RMS within the correction table (1549–2651) |

For applications requiring the osculating (instantaneous) lunar apogee, **True Lilith** (`calc_true_lilith`) is recommended. Its sub-arcminute precision makes it suitable for all practical astrological applications.

For smooth motion without 30-degree oscillations, **Mean Lilith** (`calc_mean_lilith`) provides predictable behaviour with ~0.003° mean difference from the standard mean-apogee definition.

## The Osculating Apogee Concept

The osculating apogee is derived from the Moon's instantaneous state vectors — the
standard approach. Small residual differences against other engines (~0.015° mean)
come from the ephemeris generation (DE440 vs older data), perturbation-model
details, and minor precession/nutation/obliquity choices; the measured head-to-head
is catalogued in the [Swiss Ephemeris Comparison](../comparison/index.md).

### The Osculating Apogee Paradox

The +/-30 degree monthly oscillation of the osculating apogee is an artifact of the two-body approximation (see Brouwer & Clemence 1961, Ch. XI). This has important implications:

1. **No "Correct" Answer**: There is no single "correct" osculating apogee because the concept itself is a simplification that does not accurately represent lunar motion.

2. **Different Implementations Valid**: Different engines' implementations are all valid mathematical representations of the osculating apogee concept, even though they differ slightly.

3. **True Lilith is Model-Dependent**: When solar perturbations are significant (as they always are for the Moon), the "osculating apogee" concept is inherently ambiguous.

The osculating apogee is a mathematical construct with limited physical meaning. The "true" apogee in terms of actual lunar distance extremes only aligns with the osculating apogee twice per month.

## The Interpolated Apogee Alternative

The "interpolated" or "natural" apogee (`INTP_APOG`) smooths out the 30-degree oscillations by interpolating between actual lunar apogee passages. This results in:

- Oscillation amplitude: **+/- 5 degrees** (vs. +/- 30 degrees for osculating)
- More physically meaningful: represents actual variation of apogee passages
- Less dependent on two-body approximation artifacts

For details on the interpolated apogee implementation, see [Interpolated Apogee Methodology](interpolated-apogee.md).

## References

1. Brouwer, D. & Clemence, G.M. "Methods of Celestial Mechanics" (1961),
   Academic Press, Chapter XI (osculating elements and perturbations)

2. Murray, C.D. & Dermott, S.F. "Solar System Dynamics" (1999),
   Cambridge University Press, Chapter 2

3. Chapront-Touze, M. & Chapront, J. "Lunar Tables and Programs from 4000 B.C.
   to A.D. 8000" (1991), Willmann-Bell

4. Meeus, J. "Astronomical Algorithms" (2nd ed., 1998), Willmann-Bell, Chapter 47

5. Brown, E.W. "An Introductory Treatise on the Lunar Theory" (1896)

6. Santoni, Francis. "Ephemerides de la lune noire vraie 1910-2010"
   (Editions St. Michel, 1993)

7. Koch, Dieter. "Was ist Lilith und welche Ephemeride ist richtig", Meridian 1/95

8. Vallado, D. "Fundamentals of Astrodynamics and Applications" (4th ed., 2013)
