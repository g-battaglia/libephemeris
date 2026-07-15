# Computational Methodology

LibEphemeris is a pure-Python astronomical ephemeris library built on NASA JPL
numerical integrations and modern IAU standards. This document describes the
principal modelling choices behind the library and the rationale for each.

## Table of Contents

- [Philosophy](#philosophy)
- [Method](#method)
  - [Ephemeris Foundation](#ephemeris-foundation)
  - [Outer Planet Body Centers](#outer-planet-body-centers)
  - [Lunar Apsides](#lunar-apsides)
  - [Delta T (TT - UT1)](#delta-t-tt---ut1)
  - [Minor Body Dynamics](#minor-body-dynamics)
- [Precision and Validation](#precision-and-validation)
- [References](#references)

## Philosophy

LibEphemeris aligns every component with the most current JPL and IAU models, so a
chart is accurate and internally self-consistent from antiquity to the far future.
The guiding principles are: use NASA JPL numerical ephemerides for planetary
and lunar state vectors, make every fallback explicit, use the official IAU/ERFA
routines for the precession–nutation–obliquity chain, and share a single time
argument (ΔT) and obliquity realization between body positions and house
angles. Model accuracy depends on body class and source; the documented strict
precision controls prevent an unnoticed minor-body downgrade.

No numerical model is admitted only because it agrees with another software
package. Each model, coefficient family, generated table, and public-data input
has an explicit provenance boundary; every project-selected threshold or
approximation is identified separately and justified by a mathematical
invariant, convergence argument, or measured error budget. The complete,
machine-checked inventory is
[Algorithm and Data Provenance](algorithm-provenance.md). It also classifies
infrastructure modules explicitly so that “no scientific content” cannot be
confused with missing documentation.

## Method

### Ephemeris Foundation

The main planetary/lunar sources are JPL **DE440** (Park et al., 2021) and
**DE441** in the ICRF. The local Skyfield path reads those kernels directly; LEB
evaluates precomputed Chebyshev representations; and the remote backend requests
state vectors from JPL Horizons. Dedicated analytical or catalog pipelines
serve nodes, apsides, hypothetical bodies, houses, and fixed stars. Minor bodies
may use SPK, optional N-body, or Keplerian sources. The `FLG_MOSEPH` flag is
accepted for API compatibility but does not select a Moshier ephemeris.

| Property | LibEphemeris |
| -------- | ------------ |
| Primary ephemeris | DE440/DE441 (2021) |
| Runtime backends | LEB, Horizons, local Skyfield |
| Analytical models | Derived bodies and minor-body fallback |
| Reference frame | ICRF 3.0 |

### Outer Planet Body Centers

JPL DE ephemerides provide positions for Jupiter, Saturn, Uranus, Neptune, and
Pluto as *system barycenters* — the center of mass of the planet and all its
satellites — which is not the same as the physical center of the planet. The
angular offset can be significant:

| Planet  | Typical offset (geocentric) | Primary contributor  |
| ------- | --------------------------- | -------------------- |
| Jupiter | ~0.02"                      | Galilean satellites  |
| Saturn  | ~0.03"                      | Titan                |
| Neptune | ~0.01"                      | Triton               |
| Pluto   | ~0.15"                      | Charon (binary)      |
| Uranus  | ~0.003"                     | Major satellites     |

(Measured values; see
[planet-centers-spk.md](planet-centers-spk.md) for the per-planet
kilometre offsets these correspond to.)

Outer planets stored as system barycenters use a physical planet-center segment
from JPL whenever the active segment covers the epoch. `FLG_HELCTR` and
`FLG_BARYCTR` light time is computed from the observer to that body center.

If the JPL center segment is unavailable or out of range, LibEphemeris uses the
system barycenter explicitly and records that backend choice. It does not
invent a center correction from an analytical satellite theory.

### Lunar Apsides

True/osculating lunar points are obtained from JPL DE440/DE441 state vectors:
the angular-momentum vector gives the node and the eccentricity vector gives
the instantaneous apsidal direction. Mean node and apogee use published lunar
mean-element polynomials and conventional orbital-plane geometry. Interpolated
apogee and perigee use separate Delaunay perturbation series and a versioned,
hash-pinned compatibility refinement.

See [Lunar Nodes and Apsides](lunar-apsides.md),
[Interpolated Perigee](interpolated-perigee.md), and
[Interpolated Apogee](interpolated-apogee.md).

### Delta T (TT - UT1)

The difference between Terrestrial Time and Universal Time (ΔT) converts between
dynamical and civil timescales and matters most for historical dates.
LibEphemeris uses the **Stephenson, Morrison & Hohenkerk (2016)** reconstruction
(via Skyfield) blended with observed IERS values — a multi-era model that
incorporates re-analyzed pre-telescopic eclipse records (Babylonian, Chinese,
Arabic, European) and is the current scientific standard for historical ΔT.
Contemporary calculations are unaffected; historical eclipse timing and ancient
charts benefit from the updated reconstruction. The model, its piecewise structure
and the selectable alternatives are documented in full in
[delta-t.md](delta-t.md).

### Minor Body Dynamics

LibEphemeris supports three approaches for minor bodies, in order of precision:

1. **JPL SPK kernels**: downloaded by the built-in direct HTTPS client from NASA
   JPL Horizons. These contain JPL numerical-integration results and provide
   sub-arcsecond accuracy across their coverage interval.
2. **N-body integration** (optional, via `rebound` + `assist`): for dates outside
   SPK coverage, real-time gravitational N-body integration including the major
   planets' influence at each timestep.
3. **Keplerian propagation** (fallback): two-body propagation with secular
   perturbation corrections from Laplace-Lagrange theory.

SPK kernels are recommended for production use; the N-body option requires the
`[nbody]` extra. Asteroid positions thus derive from the same JPL data products
used by professional astronomers, with automatic download capability.

## Precision and Validation

Planetary reductions and derived lunar points have separate, measured error
budgets documented in the precision reference. LibEphemeris's photometric models, body constants,
and coordinate transformations are cross-validated against independent professional
astronomy tools during development:

- **IAU ERFA/SOFA** (via pyerfa) — the official IAU standard routines for
  nutation, precession, obliquity, and frame transformations. Used as a runtime
  dependency for all coordinate calculations.
- **Astropy** — the professional astronomy library used by Hubble, JWST, and major
  observatories. Used as a development dependency for independent verification of
  magnitude formulas, body constants, and coordinate pipelines.
- **Astronomical Almanac** — published photometric formulas (Mallama & Hilton
  2018) used for all planetary magnitude calculations.
- **IAU 2015 Resolution B3** — official nominal values for solar and planetary
  radii, used for apparent-diameter calculations.

For API semantics, known differences, and the clean-room validation methodology,
see
[Swiss Ephemeris Comparison](../comparison/index.md).

## References

1. Park, R.S. et al. (2021). "The JPL Planetary and Lunar Ephemerides DE440 and DE441." *Astronomical Journal*, 161(3), 105.
2. Stephenson, F.R., Morrison, L.V. & Hohenkerk, C.Y. (2016). "Measurement of the Earth's rotation: 720 BC to AD 2015." *Proceedings of the Royal Society A*, 472, 20160404.
3. Lieske, J.H. (1998). "Galilean Satellites of Jupiter. Theory E5." *Astronomy & Astrophysics Supplement*, 129, 205-217.
4. Vienne, A. & Duriez, L. (1995). "TASS1.6: Ephemerides of the major Saturnian satellites." *Astronomy & Astrophysics*, 297, 588-605.
5. IERS Conventions (2010). IERS Technical Note No. 36, ed. Petit, G. & Luzum, B.

This overview lists the principal foundations. The exhaustive per-module and
per-generator bibliography is the
[provenance registry](algorithm-provenance.md), not this short reference list.
