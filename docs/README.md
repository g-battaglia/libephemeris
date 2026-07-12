# LibEphemeris Documentation

## Getting Started

- **[Getting Started](guides/getting-started.md)** -- Installation, ephemeris tiers, first calculations, thread safety
- **[Migration from PySwissEph](guides/migration-guide.md)** -- API mapping and behavioral differences
- **[Optional Modules](guides/optional-modules.md)** -- Optional backends and extras (star catalog, n-body, SPK kernels)
- **[Precision Tuning](guides/precision-tuning.md)** -- Configuring optional dependencies for maximum precision
- **[Computation Tracing](guides/tracing.md)** -- Discover which backend computed each body (ContextVar + DEBUG logging)

## Architecture

- **[Calculation Backends](architecture/horizons-backend.md#calculation-modes)** -- Auto/LEB/Horizons/Skyfield mode selection
- **[Horizons API Backend](architecture/horizons-backend.md)** -- Zero-install ephemeris via NASA JPL Horizons REST API
- **[Architecture Overview](development/architecture-overview.md)** -- Codebase metrics, module structure, performance analysis

## API Reference

- **[Complete API Reference](api_reference.rst)** -- Every public function, class, and constant with signatures, parameters, and examples

## Reference

- **[Precision Report](reference/precision.md)** -- Measured accuracy across all body categories
- **[Flag Reference](reference/flags.md)** -- Complete calculation flag documentation with examples
- **[House Systems](reference/house-systems.md)** -- Mathematical documentation for all 25 house systems (26 codes with A/E alias)
- **[Ayanamsha Modes](reference/ayanamsha.md)** -- Complete reference for all 47 predefined sidereal modes (+ user-defined)
- **[Known Bugs](reference/known-bugs.md)** -- Active issues and Horizons limitations

## Methodology

- **[Overview](methodology/overview.md)** -- Principal computational approaches
- **[Independence and Methodology](methodology/independence.md)** -- What differs from the reference stack: data sources, reduction chain, architecture, and how parity is measured
- **[Delta T (ΔT)](methodology/delta-t.md)** -- The multi-era TT−UT1 model and selector
- **[Long-term sidereal time & precession](methodology/sidereal-time-longterm.md)** -- Vondrák 2011, house cusps over ±13,000 years
- **[Planet Centers](methodology/planet-centers-spk.md)** -- Barycenter vs body center corrections for outer planets
- **[Lunar Apsides](methodology/lunar-apsides.md)** -- Perigee and apogee computation
- **[Interpolated Apogee](methodology/interpolated-apogee.md)** -- INTP_APOG and INTP_PERG
- **[Interpolated Perigee](methodology/interpolated-perigee.md)** -- Passage-interpolated harmonic fitting
- **[True Lilith](methodology/true-lilith.md)** -- Osculating lunar apogee calculation
- **[Hypothetical Bodies](methodology/hypothetical-bodies.md)** -- Hamburg-school Uranian planets and trans-Plutonian bodies
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
- **[Roadmap](development/roadmap.md)** -- Project status and open tasks
- **[Precision History](development/precision-history.md)** -- Record of precision fixes and investigations
- **[Keplerian Improvements](development/keplerian-improvements.md)** -- Fallback orbit propagation improvements
- **[Full Range Coverage](development/full-range-coverage.md)** -- Extending minor body coverage

## Manuals

Beginner-friendly introductions to astronomical and astrological calculations.
No prior knowledge of astronomy or programming required.

- **[Manuale (Italiano)](manual/it/)** -- 15 capitoli
- **[Manual (English)](manual/en/)** -- 15 chapters
