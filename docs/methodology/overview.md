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
The guiding principles are: use NASA JPL numerical ephemerides directly (no
reduced-precision fallback), use the official IAU/ERFA routines for the
precession–nutation–obliquity chain, and share a single time argument (ΔT) and a
single obliquity realization between body positions and house angles. All
modelling choices produce sub-arcsecond effects for the planets and sub-degree
effects for the lunar apsides.

## Method

### Ephemeris Foundation

LibEphemeris uses exclusively JPL **DE440** (Park et al., 2021) or **DE441**, the
most recent numerical planetary ephemerides — the same data NASA uses for
spacecraft navigation. Every calculation, without exception, uses the full JPL
numerical integration in the ICRF 3.0 reference frame. There is no analytical
reduced-precision fallback: the `FLG_MOSEPH` flag is accepted for API
compatibility but ignored, so users never have to manage ephemeris files or worry
about silent accuracy degradation.

| Property | LibEphemeris |
| -------- | ------------ |
| Primary ephemeris | DE440/DE441 (2021) |
| Analytical fallback | None (JPL always) |
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

LibEphemeris corrects to the true planet body center **automatically** for every
calculation, with no user configuration, using a three-tier fallback:

1. **Tier 1 — SPK-based** (<0.001"): bundled planet-center segments extracted from
   JPL satellite ephemerides (jup204, sat319, ura083, nep050, plu017).
2. **Tier 2 — Analytical satellite models** (sub-arcsecond): Lieske E5 for
   Jupiter's Galilean satellites, TASS 1.7 (Vienne & Duriez, 1995) for Saturn's
   system, Keplerian models for Neptune's Triton and Pluto's Charon.
3. **Tier 3 — Raw barycenter**: only when both higher tiers are unavailable.

For Jupiter this transparently eliminates up to 0.6 arcseconds of systematic error.

### Lunar Apsides

The interpolated (or "natural") lunar perigee and apogee represent smoothed
apsidal positions, removing the ~30° oscillations inherent in the osculating
(instantaneous) elements. LibEphemeris derives the interpolated apsides from passage geometry on the
JPL DE440/DE441 ephemeris, with the fitted coefficients and the live
residual table calibrated against reference-API output (black-box oracle;
NOTICE.md, Calibration Data Disclosure):

1. Identify all perigee and apogee passages (local distance extrema) from JPL
   state vectors over a 1000-year span.
2. At each passage, the Moon's ecliptic longitude defines the true apsidal
   position unambiguously.
3. Cubic spline interpolation through these passage points produces a smooth,
   continuous apsidal longitude function.
4. A 66-term harmonic perturbation series is fitted to this function via least
   squares.
5. A residual correction table absorbs remaining model error.

For detailed perigee calibration methodology, see
[interpolated-perigee.md](interpolated-perigee.md); for apogee-specific
methodology, see [interpolated-apogee.md](interpolated-apogee.md).

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

1. **JPL SPK kernels**: downloaded directly from NASA JPL Horizons via
   `astroquery`. These contain the same numerical integration results used by JPL
   and provide sub-arcsecond accuracy across their full date range.
2. **N-body integration** (optional, via `rebound` + `assist`): for dates outside
   SPK coverage, real-time gravitational N-body integration including the major
   planets' influence at each timestep.
3. **Keplerian propagation** (fallback): two-body propagation with secular
   perturbation corrections from Laplace-Lagrange theory.

SPK kernels are recommended for production use; the N-body option requires the
`[nbody]` extra. Asteroid positions thus derive from the same JPL data products
used by professional astronomers, with automatic download capability.

## Precision and Validation

All modelling choices produce sub-arcsecond effects for the planets and sub-degree
effects for the lunar apsides. LibEphemeris's photometric models, body constants,
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

For a head-to-head comparison with Swiss Ephemeris — measured precision, known
differences, and validation methodology — see
[Swiss Ephemeris Comparison](../comparison/index.md).

## References

1. Park, R.S. et al. (2021). "The JPL Planetary and Lunar Ephemerides DE440 and DE441." *Astronomical Journal*, 161(3), 105.
2. Chapront-Touzé, M. & Chapront, J. (1988). "ELP 2000-82B: A semi-analytical lunar ephemeris." *Astronomy & Astrophysics*, 190, 342-352.
3. Stephenson, F.R., Morrison, L.V. & Hohenkerk, C.Y. (2016). "Measurement of the Earth's rotation: 720 BC to AD 2015." *Proceedings of the Royal Society A*, 472, 20160404.
4. Lieske, J.H. (1998). "Galilean Satellites of Jupiter. Theory E5." *Astronomy & Astrophysics Supplement*, 129, 205-217.
5. Vienne, A. & Duriez, L. (1995). "TASS1.6: Ephemerides of the major Saturnian satellites." *Astronomy & Astrophysics*, 297, 588-605.
6. Moshier, S.L. (1992). "Comparison of a 7000-year lunar ephemeris with analytical theory." *Astronomy & Astrophysics*, 262, 613-616.
7. IERS Conventions (2010). IERS Technical Note No. 36, ed. Petit, G. & Luzum, B.
