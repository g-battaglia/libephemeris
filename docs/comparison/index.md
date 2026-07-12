# Swiss Ephemeris Comparison

This section is the **single place** in the libephemeris documentation where the
library is compared head-to-head with Swiss Ephemeris (via its Python wrapper
`pyswisseph`). The methodology pages elsewhere in these docs explain *what*
libephemeris computes and *why*, on their own merits; this section is where the
direct, measured comparison lives.

libephemeris targets 1:1 API compatibility with Swiss Ephemeris, so existing code
ports with minimal changes. The two engines produce results that are extremely
close in the modern era but are not bit-identical, because they are built on
different foundations:

| | libephemeris | Swiss Ephemeris (pyswisseph) |
|---|---|---|
| Planetary ephemeris | NASA JPL DE440/DE441 (2021) | JPL DE431 (2013), repackaged + Moshier fallback |
| Lunar model | DE440 numerical integration | ELP/MPP02 analytical + DE431 |
| Precession / obliquity | Vondrák 2011 long-term (valid ±200 ka) via ERFA | IAU 2006 polynomial |
| Nutation | IAU 2006/2000A via pyerfa | IAU 2006/2000A (internal) |
| ΔT (TT − UT1) | IERS observed + SMH-2016, differs only in the deep-past / future extrapolation | IERS observed + SMH-2016 (≥ SE 2.06) |
| Outer-planet center | Body center (SPK + analytical), automatic | System barycenter by default |
| Velocities | Central finite difference (numerical) | Chebyshev derivative (analytical) |
| Implementation | Pure Python + Skyfield + pyerfa | C with Python wrapper |

Both engines produce scientifically valid results. Where they differ, the cause
is documented below, with the underlying numbers, so users know exactly what to
expect.

## Pages in this section

- **[Measured precision](precision.md)** — arcsecond-level precision tables for
  every category (planets, Moon, nodes, houses, fixed stars, eclipses,
  coordinate modes) plus independent triangulation against JPL Horizons and
  astropy/ERFA.
- **[Known differences](known-differences.md)** — the *why* behind each
  divergence (ephemeris generation, lunar model, ΔT, long-term sidereal time,
  asteroids), followed by a granular per-API catalog.
- **[Intentional divergences](intentional-divergences.md)** — the few cases where
  libephemeris deliberately departs from Swiss Ephemeris behaviour for physical
  correctness (`SIDEREAL | J2000` for lunar nodes/apsides; total-eclipse
  obscuration).
- **[API compatibility](api-compatibility.md)** — function-signature
  differences, the validation methodology, and validation results.

> All measurements here were taken with the reference API's Python binding
> used purely as a black-box oracle.
