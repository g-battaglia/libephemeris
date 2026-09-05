# PyERFA Integration

LibEphemeris integrates PyERFA (the Python wrapper for the IAU SOFA/ERFA library) to provide IAU-standard nutation and precession models with milliarcsecond (mas) precision for coordinate transformations.

## Table of Contents

- [Background](#background)
- [Method](#method)
  - [Nutation Models](#nutation-models)
  - [Error Growth Over Time](#error-growth-over-time)
  - [Obliquity](#obliquity)
  - [Precession-Nutation Matrix](#precession-nutation-matrix)
  - [Long-Term Precession (Vondrák 2011)](#long-term-precession-vondrák-2011)
  - [Cached Nutation](#cached-nutation)
- [API Reference](#api-reference)
- [Precision and Validation](#precision-and-validation)
- [Runtime Role](#runtime-role)
- [Installation](#installation)
- [References](#references)

## Background

Precession and nutation are fundamental corrections required for transforming celestial coordinates between epochs. Precession describes the slow, secular drift of Earth's rotational axis over millennia, while nutation describes the shorter-period oscillations superimposed on this drift, driven primarily by the gravitational torques of the Moon and Sun on Earth's equatorial bulge.

The International Astronomical Union (IAU) has adopted progressively refined
models for these effects. The IAU 2000A nutation model uses ~1300 terms to
achieve ~0.2 mas precision. The IAU 2006 precession model improves upon the IAU
2000 precession by incorporating updated rate corrections derived from VLBI
observations. PyERFA provides reference implementations of these models and is
a required runtime dependency, giving every LibEphemeris installation the same
observatory-grade reduction chain.

## Method

### Nutation Models

PyERFA exposes three relevant nutation functions. LibEphemeris always uses
`nut06a()` in its runtime reduction; the others support validation and model
comparison:

| Model | Precision | Terms | LibEphemeris role |
|-------|-----------|-------|-------------------|
| IAU 2000A | ~0.2 mas | ~1300 | Validation/model comparison |
| IAU 2000B | ~1 mas | ~80 | Validation/model comparison |
| IAU 2006/2000A (`nut06a`) | ~0.01-0.05 mas | ~1300 | Required runtime model |

The nut06a model combines the IAU 2006 precession with IAU 2000A nutation, providing the best available precision for equinox-based coordinates. It applies small corrections to the IAU 2000A nutation to account for the updated precession rates.

### Error Growth Over Time

The differences between models grow with distance from J2000.0:

| Years from J2000 | IAU 2000B vs 2000A | nut06a vs nut00a |
|------------------|-------------------|------------------|
| 0 | ~0.2 mas | ~0.01 mas |
| 10 | ~0.5 mas | ~0.02 mas |
| 50 | ~2 mas | ~0.05 mas |
| 100 | ~5 mas | ~0.1 mas |

### Obliquity

PyERFA provides the IAU 2006 obliquity model with milliarcsecond accuracy:

```python
from libephemeris.erfa_nutation import get_erfa_obliquity_iau2006
obliquity = get_erfa_obliquity_iau2006(jd_tt)  # Returns radians
```

This helper is used for validation and model comparison. The mean obliquity
used in runtime coordinate reductions is the long-term Vondrák pole-angle
realization described under
[Long-Term Precession](#long-term-precession-vondrák-2011); near J2000 the two
agree at the sub-milliarcsecond level.

### Precession-Nutation Matrix

For complete coordinate transformations from J2000 to a date of interest:

```python
from libephemeris.erfa_nutation import get_erfa_pnm06a_matrix
# Returns the precession-nutation matrix for J2000 to date transformation
matrix = get_erfa_pnm06a_matrix(jd_tt)
```

### Long-Term Precession (Vondrák 2011)

The IAU 2006 precession is a polynomial only valid for a few centuries around
J2000; at remote epochs it diverges rapidly (for example the Sun's ecliptic
longitude is ~36" off at year -3000). For its apparent-place reduction
LibEphemeris therefore uses the **long-term precession model of Vondrák,
Capitaine & Wallace (2011)**, "New precession expressions, valid for long time
intervals" (A&A 534, A22), which is fitted to a numerical integration and stays
accurate over ±200,000 years while agreeing with IAU 2006 to sub-milliarcsecond
precision near J2000.

The precession matrix is evaluated through PyERFA's reference Vondrák
routines:

```python
import erfa
epj = 2000.0 + (jd_tt - 2451545.0) / 365.25
P = erfa.ltpb(epj)      # ICRS -> mean equator/equinox of date (frame bias included)
```

`libephemeris/precession_vondrak.py` wraps these into the matrix builders used by
every reduction path (the LEB fast path, the Skyfield reference path, and the
ecliptic-body / SPK / fixed-star paths). The of-date mean obliquity is the
angle between the Vondrák ecliptic pole (`erfa.ltpecl`) and equator pole
(`erfa.ltpequ`), evaluated in `sidereal_longterm.mean_obliquity_rad()`. Nutation remains IAU
2006/2000A and is layered on top of the Vondrák precession.

The remaining model floor at deep-BCE dates is the underlying ephemeris generation
(DE441), which the precession model does not affect.

#### Of-date mean obliquity — pole-angle realization

The of-date mean obliquity used in every reduction and reported by the
`ECL_NUT` pseudo-body is the **pole-angle realization**: the angle between the
Vondrák ecliptic-pole and equator-pole series
(`sidereal_longterm.mean_obliquity_rad()`) — the same two poles that build the
long-term precession matrix. Deriving the obliquity from those poles keeps the
precession and every equator↔ecliptic-of-date rotation in one self-consistent
frame, so a direction lying in the mean ecliptic of date (the Sun, by
definition) reduces to ~0 ecliptic latitude at every epoch. Vondrák 2011 also
publishes a direct `ε_A` polynomial-plus-periodic series; it is a separate fit,
and pairing it with the pole-based precession would tilt the of-date ecliptic
away from its own pole by up to ~6.5″ at −3000, surfacing as a spurious
ecliptic latitude on the Sun. The two realizations agree to <0.001″ across
1900–2100 (identically 0 at J2000); the direct series is not evaluated anywhere
in the library.

### Cached Nutation

For performance-critical applications, a cached variant avoids redundant calculations for the same Julian date via an LRU cache:

```python
from libephemeris.erfa_nutation import get_erfa_nutation_cached
dpsi, deps = get_erfa_nutation_cached(jd_tt)
```

## API Reference

LibEphemeris exposes the following PyERFA-based functions:

| Function | Description |
|----------|-------------|
| `has_erfa()` | Check if PyERFA is available |
| `get_erfa_nutation_nut00a(jd)` | IAU 2000A nutation (dpsi, deps in radians) |
| `get_erfa_nutation_nut06a(jd)` | IAU 2006/2000A combined nutation |
| `get_erfa_obliquity_iau2006(jd)` | Mean obliquity using IAU 2006 model |
| `get_erfa_pnm06a_matrix(jd)` | Precession-nutation matrix |
| `compare_nutation_models(jd)` | Compare all available models |
| `get_erfa_nutation_cached(jd)` | Cached nutation calculation |

## Precision and Validation

### Nutation Precision

- **IAU 2000A**: Full-precision nutation model (~0.2 mas accuracy)
- **IAU 2000B**: Truncated model (~1 mas accuracy), approximately 10x faster
- **nut06a**: IAU 2006 precession + IAU 2000A nutation (~0.05 mas accuracy)

The precision improvement is most significant for:
- Dates far from J2000.0 (before 1950 or after 2050)
- Applications requiring arcsecond or better accuracy
- Professional astronomical calculations

## Runtime Role

PyERFA is used for every standard installation and requires no user toggle.
Its `nut06a()` routine supplies nutation for the apparent-place and house
pipelines, while ERFA frame transformations keep the coordinate reductions
internally consistent. The Vondrák 2011 model supplies the separate long-term
precession behavior at remote epochs.

The alternative analytical helpers retained in the repository are useful for
tests and model comparisons only. They are not a faster selectable runtime
tier and are not used as a fallback when calculating positions.

## Installation

PyERFA is a **required runtime dependency** (declared in `pyproject.toml`) and is
installed automatically with LibEphemeris:

```bash
pip install libephemeris
```

The precession/nutation/obliquity reductions are built directly on PyERFA, so it
is always present at runtime. The pure-Python analytical helpers that remain in
`astrometry.py` (Lieske precession, the IAU 2006 obliquity polynomial, the
numpy nutation series) are **reference/test-only**: they are exercised by the
no-erfa unit tests (which monkeypatch the erfa flag off) but are not a runtime
fallback, since a real install always has PyERFA.

To confirm PyERFA is available:

```python
from libephemeris.erfa_nutation import has_erfa
assert has_erfa()  # always true in a normal install
```

## References

- IAU SOFA Library: [http://www.iausofa.org/](http://www.iausofa.org/)
- Capitaine, N. et al. "Expressions for IAU 2000 precession quantities" (2003), Astronomy & Astrophysics, 412, 567-586
- Mathews, P.M. et al. "Modeling of nutation and precession: New nutation series for nonrigid Earth and insights into the Earth's interior" (2002), Journal of Geophysical Research, 107(B4)
- IERS Conventions (2010), IERS Technical Note No. 36
