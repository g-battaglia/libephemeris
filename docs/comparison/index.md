# Compatibility comparison

LibEphemeris targets 1:1 compatibility with the public reference API so
existing programs can migrate with minimal changes. The implementation itself
uses NASA JPL, IAU/ERFA, primary literature, and published
catalogue sources; the reference API serves as an external behavioral
comparison target for the public call surface.

## Model basis

| Area | LibEphemeris basis |
|---|---|
| Planetary and lunar states | NASA JPL DE440/DE441 through Skyfield, Horizons, or reviewed LEB approximations |
| Mean lunar arguments | ERFA implementation of IERS 2003 Delaunay arguments |
| Precession and nutation | IAU 2006/2000A and Vondrák 2011 through ERFA |
| Delta T | IERS observations plus the published Stephenson–Morrison–Hohenkerk model |
| Atmosphere | Published Sæmundsson/Bennett relations and ICAO atmosphere ray tracing |
| Photometry | Published almanac and peer-reviewed planetary magnitude models |
| Historical hypothetical bodies | IDs 40–47, 50–53, and 56 from cited publications; six unsupported IDs fail closed |

Numerical differences between the two stacks are interpreted in light of
LibEphemeris's own documented model choices, listed above and detailed in the
pages below.

## Pages in this section

- [Precision and validation](precision.md) explains the accuracy checks,
  cross-backend consistency, and how compatibility observations are
  summarized.
- [Known differences](known-differences.md) catalogs documented model and API
  differences.
- [Intentional divergences](intentional-divergences.md) records deliberate public
  behavior choices and the reason for each.
- [API compatibility](api-compatibility.md) covers signatures, return shapes,
  flags, errors, and the validation process.

See [Independence and Methodology](../methodology/independence.md) for the data
sources, the reduction chain, and how parity is measured.
