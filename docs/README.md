# LibEphemeris Documentation

## Getting Started

- **[Getting Started](guides/getting-started.md)** -- Installation, ephemeris tiers, first calculations, thread safety
- **[Migration from PySwissEph](guides/migration-guide.md)** -- API mapping and behavioral differences
- **[Precision Tuning](guides/precision-tuning.md)** -- Configuring optional dependencies for maximum precision
- **[Computation Tracing](guides/tracing.md)** -- Discover which backend computed each body (ContextVar + DEBUG logging)

## Architecture

- **[Calculation Backends](architecture/horizons-backend.md#calculation-modes)** -- Auto/LEB/Horizons/Skyfield mode selection
- **[Horizons API Backend](architecture/horizons-backend.md)** -- Zero-install ephemeris via NASA JPL Horizons REST API
- **[Architecture Overview](development/architecture-overview.md)** -- Codebase metrics, module structure, performance analysis

## Reference

- **[Flag Reference](reference/flags.md)** -- Complete calculation flag documentation with examples
- **[Precision Report](reference/precision.md)** -- Measured accuracy across all body categories
- **[Known Divergences](reference/divergences.md)** -- Documented differences vs reference implementations
- **[House Systems](reference/house-systems.md)** -- Mathematical documentation for all 25 house systems (26 codes with A/E alias)
- **[Ayanamsha Modes](reference/ayanamsha.md)** -- Complete reference for all 43 sidereal zodiac modes
- **[Known Bugs](reference/known-bugs.md)** -- Active issues and Horizons limitations
- **[Comparison Report](reference/swisseph-comparison.md)** -- 1,619 cross-validation tests

## Methodology

- **[Overview](methodology/overview.md)** -- Principal computational approaches
- **[Planet Centers](methodology/planet-centers-spk.md)** -- Barycenter vs body center corrections for outer planets
- **[Lunar Apsides](methodology/lunar-apsides.md)** -- Perigee and apogee computation
- **[Interpolated Apogee](methodology/interpolated-apogee.md)** -- INTP_APOG and INTP_PERG
- **[Interpolated Perigee](methodology/interpolated-perigee.md)** -- Passage-interpolated harmonic fitting
- **[True Lilith](methodology/true-lilith.md)** -- Osculating lunar apogee calculation
- **[pyerfa Integration](methodology/pyerfa-integration.md)** -- IAU standard nutation, precession, obliquity
- **[REBOUND Integration](methodology/rebound-integration.md)** -- N-body minor body orbit propagation

## LEB Binary Ephemeris

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
