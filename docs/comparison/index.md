# Compatibility comparison

LibEphemeris targets 1:1 compatibility with the public PySwissEphemeris API so
existing programs can migrate with minimal changes. The implementation itself
uses NASA JPL, IAU/ERFA, primary literature, and published
catalogue sources.

PySwissEphemeris is used only as an external reference API for behavioral
comparison. A validation run may compare public return values in memory, but
those values are ephemeral. Raw output, per-date deltas, golden files, fitted
coefficients, recovered datasets, and encoded comparison artifacts may not be
written to this repository or its validation worktree.

## Independent model basis

| Area | LibEphemeris basis |
|---|---|
| Planetary and lunar states | NASA JPL DE440/DE441 through Skyfield, Horizons, or reviewed LEB approximations |
| Mean lunar arguments | ERFA implementation of IERS 2003 Delaunay arguments |
| Precession and nutation | IAU 2006/2000A and Vondrák 2011 through ERFA |
| Delta T | IERS observations plus the published Stephenson–Morrison–Hohenkerk model |
| Atmosphere | Published Sæmundsson/Bennett relations and ICAO atmosphere ray tracing |
| Photometry | Published almanac and peer-reviewed planetary magnitude models |
| Historical hypothetical bodies | IDs 40–47, 50–53, and 56 from cited publications; six unsupported IDs fail closed |

Current project policy prohibits inspecting or inferring the reference API's
internal models and data. Numerical differences are interpreted only in light
of LibEphemeris's own documented model choices. See the independence page for
the scope of historical remediation statements.

## Pages in this section

- [Precision and validation](precision.md) explains independent accuracy checks,
  cross-backend consistency, and how compatibility observations are summarized
  without retaining source output.
- [Known differences](known-differences.md) catalogs documented model and API
  differences without per-date reference data.
- [Intentional divergences](intentional-divergences.md) records deliberate public
  behavior choices and the independent reason for each.
- [API compatibility](api-compatibility.md) covers signatures, return shapes,
  flags, errors, and the validation process.

See [Independent Implementation Methodology](../methodology/independence.md) for
the project-wide clean-room rules.
