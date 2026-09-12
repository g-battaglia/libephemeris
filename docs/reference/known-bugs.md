# Known Bugs and Limitations

Documented issues in LibEphemeris. These are tracked and will be addressed
in future releases where possible.

## Active Bugs

No confirmed active defects are listed here. See
[Known Differences](../comparison/known-differences.md) for measured backend,
model, range, and compatibility limitations.

## Fixed Bugs (v1.0.0)

### True Node distance (FIXED)

The LEB fast path used mean orbital elements for True Node distance, causing
~3e-4 AU errors. Now uses full osculating orbit calculation from JPL DE data.
Error reduced from 3.25e-4 AU to ~1.5e-7 AU (2000x improvement).

### Sun heliocentric (FIXED)

Previously, `calc(jd, SUN, FLG_HELCTR)` returned ~126° longitude
(garbage from `arctan2` of near-zero values). Now correctly returns (0,0,0,0,0,0).

### Historical hypothetical-body availability (v3 provenance boundary)

IDs 40–48, 50–58 have documented runtime models. Vulcan (55) and Proserpina
(57) are explicit project conventions; the remaining modeled IDs use the
sources and constructions listed in
[Hypothetical bodies](../methodology/hypothetical-bodies.md). Only Nibiru (49)
raises `UnknownBodyError` because no source-complete model has been established;
see the exact
[missing-models field inventory](../methodology/missing-hypothetical-models.md).

## Horizons Backend Limitations

These are inherent to the Horizons API approach, not bugs:

| Limitation | Details |
|------------|---------|
| **Moon velocity** | ~0.001 deg/day difference vs analytical derivative (numerical differentiation) |
| **Heliocentric offset** | ~0.01-0.03" systematic difference (Horizons Sun center vs Skyfield SSB offset) |
| **Requires internet** | No offline fallback (use LEB or Skyfield for offline) |
| **Latency** | ~300-600ms first chart, cached afterward |
| **Unsupported bodies** | True Node (11), Oscu Apogee (13), Interp Apogee/Perigee (21,22), fixed stars, planetary moons — fall through to Skyfield |

Full details: [Horizons Backend Guide](../architecture/horizons-backend.md)
