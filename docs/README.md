# LibEphemeris Documentation

For the source of every algorithm, generated table, public dataset, and
project-authored numerical choice, start with
[Algorithm and Data Provenance](methodology/algorithm-provenance.md). It is the
canonical, machine-checked map of the codebase; the pages below provide the
topic-level derivations.

## Getting Started

- **[Getting Started](guides/getting-started.md)** -- Installation, ephemeris tiers, first calculations, thread safety
- **[Migration from PySwissEph](guides/migration-guide.md)** -- API mapping and behavioral differences
- **[Optional Modules](guides/optional-modules.md)** -- Optional backends and extras (star catalog, n-body, SPK kernels)
- **[Precision Tuning](guides/precision-tuning.md)** -- Configuring optional dependencies for maximum precision
- **[Computation Tracing](guides/tracing.md)** -- Discover which backend computed each body (ContextVar + DEBUG logging)

## Architecture

- **[Calculation Backends](architecture/horizons-backend.md#calculation-modes)** -- Auto/LEB/Horizons/Skyfield mode selection
- **[Horizons API Backend](architecture/horizons-backend.md)** -- Zero-install ephemeris via NASA JPL Horizons REST API
- **[Architecture Overview](development/architecture-overview.md)** -- Module structure, data flow, backends, and state ownership

## API Reference

- **[Complete API Reference](api_reference.rst)** -- Every public function, class, and constant with signatures, parameters, and examples

## Reference

- **[Precision Report](reference/precision.md)** -- Measured accuracy across all body categories
- **[Flag Reference](reference/flags.md)** -- Complete calculation flag documentation with examples
- **[House Systems](reference/house-systems.md)** -- Mathematical documentation for all 25 house systems (26 codes with A/E alias)
- **[Ayanamsha Modes](reference/ayanamsha.md)** -- All predefined modes, source-audit status, and user mode
- **[Known Bugs](reference/known-bugs.md)** -- Active issues and Horizons limitations

## Methodology

- **[Algorithm and Data Provenance](methodology/algorithm-provenance.md)** -- Complete source, derivation, generator, and project-choice map for the codebase
- **[Classical Astrology Helpers](methodology/classical-astrology-helpers.md)** -- Per-function distinction between published definitions, mathematical identities, and project conventions
- **[Overview](methodology/overview.md)** -- Principal computational approaches
- **[Independence and Methodology](methodology/independence.md)** -- What differs from the reference stack: data sources, reduction chain, architecture, and how parity is measured
- **[Delta T (ΔT)](methodology/delta-t.md)** -- The multi-era TT−UT1 model and selector
- **[Long-term sidereal time & precession](methodology/sidereal-time-longterm.md)** -- Vondrák 2011, house cusps over ±13,000 years
- **[Planet Centers](methodology/planet-centers-spk.md)** -- Barycenter vs body center corrections for outer planets
- **[Lunar Apsides](methodology/lunar-apsides.md)** -- Perigee and apogee computation
- **[Interpolated Apogee](methodology/interpolated-apogee.md)** -- INTP_APOG analytical compatibility curve
- **[Interpolated Perigee](methodology/interpolated-perigee.md)** -- INTP_PERG analytical compatibility curve
- **[True Lilith](methodology/true-lilith.md)** -- Osculating lunar apogee calculation
- **[Hypothetical Bodies](methodology/hypothetical-bodies.md)** -- per-ID primary provenance, transformations, and fail-closed status
- **[Missing Hypothetical Models](methodology/missing-hypothetical-models.md)** -- field-by-field inventory of the six fail-closed IDs and the evidence needed to restore them
- **[pyerfa Integration](methodology/pyerfa-integration.md)** -- IAU standard nutation, precession, obliquity
- **[REBOUND Integration](methodology/rebound-integration.md)** -- N-body minor body orbit propagation

## Swiss Ephemeris Comparison

The single, centralized place for the head-to-head with Swiss Ephemeris (pyswisseph).

- **[Overview](comparison/index.md)** -- Foundational differences and navigation
- **[Measured precision](comparison/precision.md)** -- Arcsecond-level precision tables + independent triangulation
- **[Known differences](comparison/known-differences.md)** -- Root causes + granular per-API catalog
- **[Intentional divergences](comparison/intentional-divergences.md)** -- Where libephemeris deliberately departs for correctness
- **[API compatibility](comparison/api-compatibility.md)** -- Signature differences, validation methodology and results

## LEB Binary Ephemeris

- **[Quickstart](leb/quickstart.md)** -- Fast path to generating and using the LEB binary ephemeris
- **[Technical Guide](leb/guide.md)** -- Format specification, reader, fast-path pipeline, LEB2 compressed format
- **[Algorithms & Theory](leb/algorithms.md)** -- Chebyshev polynomials, Clenshaw, gravitational deflection, error analysis
- **[Comparison Testing](leb/testing.md)** -- LEB vs Skyfield comparison methodology

## Development

- **[Testing](development/testing.md)** -- Test suites, commands, expected failures
- **[Roadmap](development/roadmap.md)** -- Historical March 2026 planning snapshot
- **[Precision History](development/precision-history.md)** -- Record of precision fixes and investigations
- **[Keplerian Improvements](development/keplerian-improvements.md)** -- Fallback orbit propagation improvements
- **[Full Range Coverage](development/full-range-coverage.md)** -- Historical plan for the now-implemented minor-body coverage work

## Manuals

Beginner-friendly introductions to astronomical and astrological calculations.
No prior knowledge of astronomy or programming required.

- **[Manuale (Italiano)](manual/it/README.md)** -- prefazione + 15 capitoli
- **[Manual (English)](manual/en/README.md)** -- preface + 15 chapters
