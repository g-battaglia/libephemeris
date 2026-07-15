# Compatibility Precision and Independent Validation

LibEphemeris validates astronomical accuracy against independent sources:
NASA JPL DE440/DE441 and Horizons, ERFA/SOFA, IERS, Hipparcos/Gaia/SIMBAD,
Astropy, and published defining equations. PySwissEphemeris may be called only
as an external reference API for compatibility comparison. Those observations
are ephemeral: no returned value is stored in this page, a fixture, a
coefficient, a table, or a generated artifact.

This page therefore reports model classes and non-reconstructive aggregate
behavior rather than per-date reference output. Exact regression tolerances
are documented next to their tests and justified from the relevant source's
accuracy, a propagated numerical error budget, or a declared project acceptance
policy. A finite sampled test is not presented as proof of a continuous maximum.

## Precision summary

- **Major planets and Moon:** positions originate in JPL numerical
  integrations. Differences from another ephemeris primarily reflect the
  selected JPL solution, time-scale realization, apparent-place reductions,
  and body-center convention. Modern-era angular differences are generally at
  arcsecond scale; they grow toward kernel endpoints, especially for the Moon.
- **Coordinate frames:** precession, nutation, obliquity, and frame rotations
  are checked against ERFA/SOFA. Long-term precession uses the published
  Vondrák model outside the useful range of the IAU 2006 polynomial.
- **Earth rotation:** modern ΔT uses IERS observations; historical and future
  values use a published reconstruction and extrapolation. Remote-epoch
  uncertainty is dominated by Earth rotation, not numerical solver error.
- **Fixed stars:** catalog astrometry comes from Hipparcos/Gaia-class sources
  and IAU/WGSN names. Positions and space motion are cross-checked with
  Astropy/SkyCoord and catalog services.
- **Minor bodies:** registered SPK states are checked against JPL Horizons.
  Keplerian fallback is explicitly lower precision and degrades away from its
  element epoch.
- **Houses and event solvers:** defining geometric conditions, round trips,
  monotonicity, bracketing, and numerical-derivative convergence are tested
  directly.
- **Refraction and visibility:** compact refraction uses published
  Sæmundsson/Bennett formulas; extended refraction uses an ICAO Standard
  Atmosphere ray tracer; heliacal visibility uses the published Schaefer
  VISLIMIT family and physical observer optics.
- **Eclipses and occultations:** apparent-disc contacts and Sun–Earth shadow
  cones are checked from independent geometry and published shadow conventions.

## Lunar nodes and apsides

| Body | Independent basis | Verification |
|---|---|---|
| True Node | JPL Moon state and angular momentum | Direct DE441 state-vector construction |
| Mean Node | ERFA/IERS Delaunay node argument | Standards values and smooth rates |
| Osculating Apogee | JPL eccentricity vector | Direct two-body extraction |
| Mean Apogee | ERFA/IERS Delaunay identities | Inclined-plane geometry and derivative checks |
| Interpolated Apogee/Perigee | DE440 apsis passages + Delaunay series + hash-pinned residual tables | Exactness at the DE440 passages, endpoint, continuity, and derivative checks |

Interpolated points are anchored to the actual DE440 apsis passages
(exact there by construction); ephemeral comparison with external
implementations between passages stays within ~0.05° in longitude and
`0.25°` in latitude, without retaining any sampled output.

## Independent triangulation

The main adjudication paths are:

1. Compare LEB output with the Skyfield/JPL source state used to generate it.
2. Compare frame and time transformations with ERFA/SOFA and IERS.
3. Compare minor-body states with JPL Horizons over their SPK coverage.
4. Compare fixed-star propagation with Astropy/SkyCoord and authoritative
   catalog astrometry.
5. Check derived points and event times against their defining geometry using
   independent numerical implementations, step halving, and metamorphic
   relations.

Compatibility calls can reveal API-shape, flag, and convention differences,
but they are not an accuracy authority and their output is discarded after the
comparison.

## Expected model-dependent differences

The largest predictable families are:

- lunar positions near the ends of a JPL kernel;
- deep-time ΔT and sidereal-time extrapolation;
- minor bodies when only a Keplerian fallback is available;
- star distances when catalog radial velocities or parallax policies differ;
- abstract mean/smoothed lunar points with different definitions;
- below-horizon refraction and heliacal visibility under different atmosphere
  or observer models;
- grazing eclipse/occultation classification under different limb and shadow
  conventions;
- historical ayanamshas under different published zero points or precession
  realizations.

See [Known differences](known-differences.md) for causes and
[Intentional divergences](intentional-divergences.md) for deliberately retained
semantic choices.
