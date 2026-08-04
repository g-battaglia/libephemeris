# Independence and Methodology

LibEphemeris shares an API surface with the Swiss Ephemeris ecosystem so
existing code can migrate without changes. Its implementation is built from
the sources and components documented below. This page describes that stack
and how behavioral parity is achieved.

## Different data sources

| Layer | LibEphemeris | Compatibility comparison |
|---|---|---|
| Planetary/lunar ephemeris | **NASA JPL DE440/DE441** through the local Skyfield path, precomputed LEB coefficients, or JPL Horizons | Numeric differences are measured from public API output only. |
| Asteroids / minor bodies | JPL **SBDB** osculating elements, JPL **Horizons** SPK kernels, optional ASSIST N-body propagation | Coverage and values are compared through public API calls. |
| Fixed stars | Catalog built from **Hipparcos** (van Leeuwen 2007), Tycho-2, VizieR and the IAU WGSN name list — public astrometry only | Names and returned astrometry are compared through public API calls. |
| Outer-planet centers | JPL satellite ephemerides (jup204, sat319, ura083, nep050, plu017) for barycenter→body-center correction | Center/barycenter distinctions are established from output comparisons. |
| Hypothetical bodies | Neely IDs 40–47, Harrington 50, Le Verrier 51, Adams 52, and Lowell 53 from page-level primary transcriptions; Selena 56 from Velichko–Larin's published seven-year rule, 1800–2000 endpoint unwrap, and independent 1879–2007 checkpoints; IDs 48, 49, 54, 55, 57, and 58 are recognised but unsupported | API shape and error behavior are checked ephemerally; unsupported source gaps are listed in [the missing-models inventory](missing-hypothetical-models.md). |

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
  cited public definition or independently documented construction, driven by
  a shared spherical-geometry engine. The four v3-remediated vector
  constructions have focused primary-source/provenance tests.
- **Lunar points:** true/osculating geometry comes from DE440/DE441 state
  vectors. Mean points use ERFA/IERS Delaunay arguments. Interpolated
  points use independently fitted Delaunay-series models and residual tables
  generated from DE440 apsis passages and pinned by SHA-256; see
  [lunar-apsides.md](lunar-apsides.md).

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
behavior is established by measuring public outputs; the reference API is
compared behaviorally, and its output is not used as a data source for
fixtures, coefficients, tables, or generated artifacts.
Where the two engines legitimately disagree
(different underlying ephemeris, different ΔT, physical-model choices), the
difference is measured and documented in
[known-differences.md](../comparison/known-differences.md) and
[intentional-divergences.md](../comparison/intentional-divergences.md)
rather than papered over.

Each currently registered model must have a model-specific derivation from
published JPL/IAU standards, primary literature, public catalogues, or a
permissively licensed source. The data-source inventory is stated in
[NOTICE.md](https://github.com/g-battaglia/libephemeris/blob/main/NOTICE.md).
The exhaustive source, derivation, project-choice, generator, and documentation
coverage is maintained in the machine-checked
[Algorithm and Data Provenance](algorithm-provenance.md) registry. A new Python
module or shipped scientific asset cannot pass that gate until it is classified
and its in-code provenance contract is present.

## What this means in practice

- Same function calls, same flags, same return shapes as the reference API.
- Positions come from JPL DE440/DE441 and IAU standard reductions; typical
  agreement with the reference engine is at the sub-arcsecond level, with
  every systematic difference measured and documented.
- No runtime dependency on a reference-product component. The package is
  licensed AGPL-3.0-only.

## Behavioral-compatibility inventory

For auditability, the public-API conventions adopted purely for
interoperability — observable semantics of a published API, never numeric
data — are concentrated in these families:

- **Call surface**: function names, argument order, tuple shapes and
  lengths, flag bit values, return-flag echo rules (which request bits are
  echoed, implied, or consumed), typed-error classes per input family.
- **Semantic conventions of documented output channels**: which frame a
  flag combination selects (e.g. the mean frame of a mode's t0 under the
  ECL_T0 projection), precedence between projection bits, which sidereal
  modes suppress a projection, the local-noon eclipse slot being defined
  only within a resolved window, the disc-area obscuration convention, the
  Hindu-rising bit's latitude-zeroing, and the interpolated-apsis point
  *concept* (its realization here is fitted to JPL DE440 apsis passages
  only).
- **Body-class behavior**: which ids compute, raise, or are consumed by a
  flag (e.g. the lunar points under the topocentric flag), and the
  trans-jovian barycentric re-referencing class.

Everything in this inventory is behavior observable from the public API of
any conforming implementation. Where an external implementation's behavior
is *not* derivable from published definitions and would alter numeric
output, it is deliberately not reproduced and is recorded in
[intentional divergences](../comparison/intentional-divergences.md).

"Swiss Ephemeris" is a product of Astrodienst AG; the name is used here
nominatively only.
