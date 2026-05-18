# Computational Methodology

LibEphemeris is a pure-Python astronomical ephemeris library built on NASA JPL numerical integrations and modern IAU standards. This document describes the principal areas in which LibEphemeris's computational approach differs from that of Swiss Ephemeris and the rationale for each choice.

## Table of Contents

- [Background](#background)
- [Method](#method)
  - [Ephemeris Foundation and Analytical Fallback](#ephemeris-foundation-and-analytical-fallback)
  - [Outer Planet Body Centers](#outer-planet-body-centers)
  - [Lunar Apsides](#lunar-apsides)
  - [Delta T (TT - UT1)](#delta-t-tt---ut1)
  - [Minor Body Dynamics](#minor-body-dynamics)
- [Precision and Validation](#precision-and-validation)
- [Comparison with Swiss Ephemeris](#comparison-with-swiss-ephemeris)
- [References](#references)

## Background

These differences are methodological, not errors. Both libraries produce scientifically valid results. LibEphemeris prioritizes consistency with current JPL ephemerides and IAU models; Swiss Ephemeris prioritizes backward compatibility and computational speed.

All differences produce sub-arcsecond effects for the planets and sub-degree effects for the lunar apsides. LibEphemeris chooses to align with the most current JPL and IAU standards at the cost of strict numerical agreement with Swiss Ephemeris.

## Method

### Ephemeris Foundation and Analytical Fallback

Swiss Ephemeris supports three computation modes: (1) its own binary repackaging of JPL DE431, (2) a fallback to the Moshier semi-analytical ephemeris when binary files are absent (`FLG_MOSEPH`), and (3) direct JPL file reading.

The Moshier mode uses VSOP87 for the planets and ELP2000-82B for the Moon — analytical theories developed in the 1980s. These introduce errors on the order of 1 arcsecond for the inner planets and 10+ arcseconds for outer planets at dates far from the fitting epoch.

LibEphemeris uses exclusively JPL DE440 (Park et al., 2021) or DE441, the most recent numerical planetary ephemerides. These are the same models used by NASA for spacecraft navigation.

The `FLG_MOSEPH` flag is accepted for API compatibility but silently ignored. Every calculation — without exception — uses the full JPL numerical integration. There is no reduced-precision fallback.

| Property            | Swiss Ephemeris                    | LibEphemeris           |
| ------------------- | ---------------------------------- | ---------------------- |
| Primary ephemeris   | DE431 (2013)                       | DE440/DE441 (2021)     |
| Analytical fallback | Moshier (VSOP87 + ELP2000-82B)     | None (JPL always)      |
| Reference frame     | ICRF 2.0                           | ICRF 3.0               |

Users migrating from pyswisseph do not need to manage ephemeris files or worry about silent accuracy degradation. LibEphemeris produces consistent precision regardless of configuration.

### Outer Planet Body Centers

JPL DE ephemerides provide positions for Jupiter, Saturn, Uranus, Neptune, and Pluto as *system barycenters* — the center of mass of the planet and all its satellites. This is not the same as the physical center of the planet.

The angular offset between barycenter and body center can be significant:

| Planet  | Maximum offset | Primary contributor  |
| ------- | -------------- | -------------------- |
| Jupiter | ~0.6"          | Galilean satellites  |
| Saturn  | ~0.2"          | Titan                |
| Neptune | ~0.05"         | Triton               |
| Pluto   | ~0.3"          | Charon (binary)      |
| Uranus  | ~0.01"         | Major satellites     |

Swiss Ephemeris returns system barycenters by default. A separate flag (`FLG_CENTER_BODY`) and additional satellite ephemeris files are required to obtain planet body centers.

LibEphemeris corrects to the true planet body center *automatically* for every calculation, using a three-tier fallback:

1. **Tier 1 — SPK-based** (<0.001"): Bundled planet center segments extracted from JPL satellite ephemerides (jup204, sat319, ura083, nep050, plu017).
2. **Tier 2 — Analytical satellite models** (sub-arcsecond): Rigorous theories for major moons — Lieske E5 for Jupiter's Galilean satellites, TASS 1.7 (Vienne & Duriez, 1995) for Saturn's system, Keplerian models for Neptune's Triton and Pluto's Charon.
3. **Tier 3 — Raw barycenter**: Used only when both higher tiers are unavailable.

No user configuration is required. The correction is transparent. For Jupiter, this eliminates up to 0.6 arcseconds of systematic error without any user action.

### Lunar Apsides

The interpolated (or "natural") lunar perigee and apogee represent smoothed apsidal positions, removing the ~30 degree oscillations inherent in the osculating (instantaneous) elements. Computing this smoothed position requires a methodological choice about what constitutes the "true" apsidal motion.

Swiss Ephemeris uses the analytical approach developed by S.L. Moshier, based on ELP2000-82B lunar theory (Chapront-Touzé & Chapront, 1988). This method isolates specific trigonometric terms in the analytical series and removes those associated with short-period perturbations.

LibEphemeris grounds the interpolated apsides in the physical geometry of the JPL DE440/DE441 ephemeris:

1. Identify all perigee and apogee passages (local distance extrema) from JPL state vectors over a 1000-year span.
2. At each passage, the Moon's ecliptic longitude defines the true apsidal position unambiguously.
3. Cubic spline interpolation through these passage points produces a smooth, continuous apsidal longitude function.
4. A 61-term harmonic perturbation series is fitted to this function via least squares.
5. A residual correction table absorbs remaining model error.

The interpolated perigee (`INTP_PERG`) differs from Swiss Ephemeris by up to ~5 degrees. This is the largest single discrepancy between the two libraries and reflects the different smoothing philosophies applied to the same physical phenomenon.

| Aspect                | Swiss Ephemeris                  | LibEphemeris                         |
| --------------------- | -------------------------------- | ------------------------------------ |
| Ground truth          | Analytical lunar theory (1988)   | JPL numerical ephemeris (2021)       |
| Smoothing             | Term-based filtering             | Passage interpolation                |
| Perigee amplitude     | ~15 deg from mean                | ~25 deg from mean                    |
| Apogee amplitude      | ~5 deg from mean                 | ~5 deg from mean                     |

For detailed perigee calibration methodology, see [interpolated-perigee.md](interpolated-perigee.md).
For apogee-specific methodology, see [interpolated-apogee.md](interpolated-apogee.md).

### Delta T (TT - UT1)

The difference between Terrestrial Time and Universal Time (Delta T) is essential for converting between dynamical and civil timescales, particularly for historical dates.

Swiss Ephemeris uses the polynomial model of Espenak & Meeus (2006), which is fitted to telescopic observations from the 17th century onward and extrapolated for earlier dates.

LibEphemeris uses the Stephenson, Morrison & Hohenkerk (2016) model via Skyfield. This is a more recent model that incorporates re-analyzed pre-telescopic observations — Babylonian, Chinese, Arabic, and European eclipse records — and is considered the current standard for historical Delta T computation.

| Date range            | Method                                                   |
| --------------------- | -------------------------------------------------------- |
| 720 BCE — ~2016 CE    | Cubic spline interpolation (Stephenson et al. Table S15) |
| Outside spline range  | Parabolic: Delta T = -320 + 32.5 × u²                   |
| 1973 — present        | IERS observed values (optional)                          |

For modern dates (1900–2100), both models agree to within ~1 second. The difference is most pronounced for dates before ~1600 CE, where the Stephenson model's broader observational basis provides improved accuracy.

Historical eclipse timing, ancient chart calculations, and long-range date conversions benefit from the updated model. Contemporary calculations are unaffected.

### Minor Body Dynamics

Swiss Ephemeris uses proprietary binary files (`.se1`) containing pre-computed Chebyshev polynomial approximations for asteroid positions. These are fast to evaluate but limited to the date range covered by the file, and precision can degrade at the boundaries.

LibEphemeris supports three approaches for minor bodies, in order of precision:

1. **JPL SPK kernels**: Downloaded directly from NASA JPL Horizons via `astroquery`. These contain the same numerical integration results used by JPL and provide sub-arcsecond accuracy across their full date range.
2. **N-body integration** (optional, via `rebound` + `assist`): For dates outside SPK kernel coverage, LibEphemeris can perform real-time gravitational N-body integration, computing the gravitational influence of the major planets on the asteroid at each timestep.
3. **Keplerian propagation** (fallback): Two-body propagation with secular perturbation corrections from Laplace-Lagrange theory.

The first approach (SPK kernels) is recommended for production use. The N-body option requires the optional `[nbody]` extra.

Asteroid positions are derived from the same JPL data products used by professional astronomers, with automatic download capability. For work at extreme date ranges, real-time N-body integration provides a physically rigorous alternative to pre-computed polynomial tables.

## Precision and Validation

All differences produce sub-arcsecond effects for the planets and sub-degree effects for the lunar apsides. Both libraries are scientifically valid. LibEphemeris chooses to align with the most current JPL and IAU standards at the cost of strict numerical agreement with Swiss Ephemeris.

### Cross-validation against professional astronomy standards

LibEphemeris photometric models, body constants, and coordinate transformations are cross-validated against independent professional astronomy tools during development:

- **IAU ERFA/SOFA** (via pyerfa) — Official IAU standard routines for nutation, precession, obliquity, and frame transformations. Used as a runtime dependency for all coordinate calculations.
- **Astropy** — The professional astronomy library used by Hubble, JWST, and major observatories worldwide. Used as a development dependency for independent verification of magnitude formulas, body constants, and coordinate pipelines.
- **Astronomical Almanac** — Published photometric formulas (Mallama & Hilton 2018) used for all planetary magnitude calculations.
- **IAU 2015 Resolution B3** — Official nominal values for solar and planetary radii, used for apparent diameter calculations.

Where LibEphemeris and Swiss Ephemeris produce different results, the cross-validation determines which value is more astronomically accurate. In several cases — notably apparent diameters (IAU equatorial radii vs mean volumetric), Delta T model (2016 vs 2006), and ephemeris generation (DE440 vs DE431) — LibEphemeris uses the more current and accurate standard.

## Comparison with Swiss Ephemeris

| Area                     | Swiss Ephemeris                        | LibEphemeris                               |
| ------------------------ | -------------------------------------- | ------------------------------------------ |
| Ephemeris source         | DE431 (2013) + Moshier fallback        | DE440/DE441 (2021), no fallback            |
| Outer planet positions   | System barycenters (default)           | Planet body centers (default)              |
| Lunar apsides            | ELP2000-82B analytical filtering       | JPL passage-interpolated fitting           |
| Delta T model            | Espenak & Meeus (2006)                 | Stephenson, Morrison & Hohenkerk (2016)    |
| Minor body dynamics      | Pre-computed Chebyshev polynomials     | JPL SPK + optional N-body integration      |
| Nutation model           | IAU 2006/2000A                         | IAU 2006/2000A (via pyerfa)                |
| Magnitude formulas       | Internal implementation                | Mallama & Hilton 2018 (peer-reviewed)      |
| Body radii               | Mean volumetric                        | IAU 2015 equatorial (Astronomical Almanac) |
| Validation               | Internal test suite                    | Cross-validated against astropy, ERFA, IAU |

## References

1. Park, R.S. et al. (2021). "The JPL Planetary and Lunar Ephemerides DE440 and DE441." *Astronomical Journal*, 161(3), 105.
2. Chapront-Touzé, M. & Chapront, J. (1988). "ELP 2000-82B: A semi-analytical lunar ephemeris." *Astronomy & Astrophysics*, 190, 342-352.
3. Stephenson, F.R., Morrison, L.V. & Hohenkerk, C.Y. (2016). "Measurement of the Earth's rotation: 720 BC to AD 2015." *Proceedings of the Royal Society A*, 472, 20160404.
4. Espenak, F. & Meeus, J. (2006). "Five Millennium Canon of Solar Eclipses." NASA/TP-2006-214141.
5. Lieske, J.H. (1998). "Galilean Satellites of Jupiter. Theory E5." *Astronomy & Astrophysics Supplement*, 129, 205-217.
6. Vienne, A. & Duriez, L. (1995). "TASS1.6: Ephemerides of the major Saturnian satellites." *Astronomy & Astrophysics*, 297, 588-605.
7. Moshier, S.L. (1992). "Comparison of a 7000-year lunar ephemeris with analytical theory." *Astronomy & Astrophysics*, 262, 613-616.
8. IERS Conventions (2010). IERS Technical Note No. 36, ed. Petit, G. & Luzum, B.
