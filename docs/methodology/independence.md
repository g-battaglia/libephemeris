# Independence and Methodology

LibEphemeris is an independent implementation: it shares an API surface with
the Swiss Ephemeris ecosystem so existing code can migrate without changes,
while its implementation is built from the independent sources and components
documented below. This page describes that stack and how behavioral parity is
achieved without sharing any code or implementation material.

## Different data sources

| Layer | LibEphemeris | Compatibility comparison |
|---|---|---|
| Planetary/lunar ephemeris | **NASA JPL DE440/DE441** through the local Skyfield path, precomputed LEB coefficients, or JPL Horizons | Numeric differences are measured from public API output only. |
| Asteroids / minor bodies | JPL **SBDB** osculating elements, JPL **Horizons** SPK kernels, optional ASSIST N-body propagation | Coverage and values are compared through public API calls. |
| Fixed stars | Catalog built from **Hipparcos** (van Leeuwen 2007), Tycho-2, VizieR and the IAU WGSN name list — public astrometry only | Names and returned astrometry are compared through public API calls. |
| Outer-planet centers | JPL satellite ephemerides (jup204, sat319, ura083, nep050, plu017) for barycenter→body-center correction | Center/barycenter distinctions are established from output comparisons. |
| Hypothetical bodies | Published primary sources (Witte & Lefeldt 1928, Neely 1988, Hoyt 1980, …) in an auditable CSV | Corresponding public body outputs are used only as interoperability targets. |

The local Skyfield path reads DE440/DE441 SPK kernels directly. The LEB path is
an intentional intermediate representation generated from the same JPL state
data, with its compression error bounded and verified separately; the Horizons
path receives JPL-produced state vectors over HTTPS.

## Independently sourced reduction chain

The apparent-place pipeline (geometric position → light-time → gravitational
deflection → aberration → frame rotation) is built on the **official IAU
routines** rather than on a private reimplementation:

- **Precession–nutation–obliquity:** IAU 2006/2000A via
  [pyerfa](https://github.com/liberfa/pyerfa) (the reference C library of the
  IAU Standards of Fundamental Astronomy), with the **Vondrák 2011**
  long-term precession model for epochs far from J2000 (valid ±200,000
  years). See [pyerfa-integration.md](pyerfa-integration.md) and
  [sidereal-time-longterm.md](sidereal-time-longterm.md).
- **ΔT (TT↔UT1):** the Stephenson–Morrison–Hohenkerk 2016 spline (with the
  Morrison et al. 2021 revision) blended with IERS observations — the most
  current published reconstruction of Earth rotation. See
  [delta-t.md](delta-t.md).
- **House systems:** each of the 25 systems is implemented from its
  published definition (Polich & Page 1964, Krusinski/Pisa/Goelzer, Savard,
  Makransky 1988, Knegt, Raman, …), driven by a shared spherical-geometry
  engine.
- **Lunar theory helpers:** ELP 2000-85 / Chapront series for the mean
  elements; osculating elements extracted from the DE440 lunar state for the
  true node and apsides.

## Own architecture

Several components are LibEphemeris-native features outside the compatibility
surface:

- **LEB / LEB2 binary format** — a precomputed Chebyshev ephemeris designed
  for this library (~14× faster than the general pipeline), with
  error-bounded lossy compression (4–10× smaller, <0.001″ deviation). See
  [the LEB guide](../leb/guide.md).
- **Horizons backend** — zero-install operation by querying the NASA JPL
  Horizons REST API when no local kernels are present.
- **Four calculation modes** — `auto`, `leb`, `horizons`, and `skyfield` expose
  a documented LEB → Horizons → Skyfield fallback contract.
- **Thread-safe contexts** — `EphemerisContext` isolates per-thread state
  alongside the classic global-state API.
- Pure Python throughout: inspectable algorithms, standard debugging and
  profiling, wheels on any platform.

## How parity is achieved

API compatibility extends beyond signatures to numeric behavior: flag
semantics, frame conventions, retflag echoes, edge-case handling. That
behavior is established **by measurement, not by reading source**: the
reference implementation is used strictly as a *black-box oracle* — inputs
in, outputs compared — from a separate validation repository that is not a
dependency of this package. Where the two engines legitimately disagree
(different underlying ephemeris, different ΔT, physical-model choices), the
difference is measured and documented in
[known-differences.md](../comparison/known-differences.md) and
[intentional-divergences.md](../comparison/intentional-divergences.md)
rather than papered over.

Two generated data sets are fitted against that oracle's *output* and are
disclosed as such (the interpolated lunar-apse residual tables, and four
fictitious-body element rows with no known publication) — see the
Calibration Data Disclosure in
[NOTICE.md](https://github.com/g-battaglia/libephemeris/blob/main/NOTICE.md).
They contain
computed positions, not source expression.

## What this means in practice

- Same function calls, same flags, same return shapes as the reference API.
- Positions come from JPL DE440/DE441 and IAU standard reductions; typical
  agreement with the reference engine is at the sub-arcsecond level, with
  every systematic difference measured and documented.
- No runtime dependency on any Swiss Ephemeris component; the package and
  its required dependencies are permissively licensed (Apache-2.0 core).

"Swiss Ephemeris" is a product of Astrodienst AG; the name is used here
nominatively only.
