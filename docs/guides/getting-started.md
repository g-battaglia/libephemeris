# Getting Started

## Installation

```bash
pip install libephemeris
```

The PyPI wheel includes the bundled LEB2 base-tier core (~10.23 MB) and
Hamburg-body companion (~46 KB), covering 1850–2150. The immutable `data-v3`
release supplies all five SHA-256-pinned groups for every tier. In `auto` or
`leb` mode, `libephemeris download auto` installs them cumulatively through the
configured tier (5, 10, or 15 files). Historical LEB1 files remain readable and
their standard tier filenames are auto-discovered.

For a traditional JPL/Skyfield workflow, download its DE kernel, planet-center
data and minor-body SPKs directly:

```bash
libephemeris download medium       # 1550-2650, >300 MB plus SPKs — recommended
libephemeris download base         # 1850-2150, lightweight
libephemeris download extended     # DE441 core span; minor-body ranges vary
```

### Optional Extras

```bash
pip install libephemeris[nbody]   # REBOUND/ASSIST n-body integration for TNOs
pip install libephemeris[stars]   # Star-catalog tooling (astropy + astroquery)
pip install libephemeris[all]     # Permissive runtime extras (currently astropy)
pip install libephemeris[all,nbody]  # ... plus the GPL-licensed N-body extra (provenance-ok)
```

**Requirements:** Python 3.12+ | skyfield >= 1.54 | pyerfa >= 2.0

## Ephemeris Tiers

| Tier | JPL File | Date Range | Size | Use Case |
|------|----------|------------|------|----------|
| `base` | de440s.bsp | 1849-2150 | ~31 MB | Modern-era, lightweight |
| `medium` | de440.bsp | 1550-2650 | ~114 MB | General purpose **(default)** |
| `extended` | de441.bsp | -13200 to +17191 | ~3.1 GB | Historical/far-future |

DE440s is a reduced-size subset of DE440 with no loss inside its stored range.
DE441 shares the same modern observational lineage but is the long-span
integration; it must not be described as numerically identical to DE440 over
their overlap. The best-by-date LEB router therefore prefers DE440s/DE440 and
uses DE441 for dates the shorter integration does not cover.

```python
import libephemeris as swe

# Select by tier
swe.set_precision_tier("extended")

# Or by file directly
swe.set_ephemeris_file("de441.bsp")

# Or via environment variable
# export LIBEPHEMERIS_PRECISION=extended
# export LIBEPHEMERIS_EPHEMERIS=de441.bsp
```

## First Calculation

```python
import libephemeris as swe
from libephemeris.constants import SUN, MOON, MARS, FLG_SPEED

jd = swe.julday(2024, 3, 26, 12.0)

# Planet positions with velocity
sun, _ = swe.calc_ut(jd, SUN, FLG_SPEED)
moon, _ = swe.calc_ut(jd, MOON, FLG_SPEED)
mars, _ = swe.calc_ut(jd, MARS, FLG_SPEED)

# Returns: (longitude, latitude, distance, speed_lon, speed_lat, speed_dist)
print(f"Sun:  {sun[0]:.4f} deg, speed {sun[3]:.4f} deg/day")
print(f"Moon: {moon[0]:.4f} deg, speed {moon[3]:.4f} deg/day")
print(f"Mars: {mars[0]:.4f} deg, speed {mars[3]:.4f} deg/day")
```

## House Cusps

```python
jd = swe.julday(2024, 11, 5, 18.0)

# Placidus houses for Rome (41.9N, 12.5E)
cusps, ascmc = swe.houses(jd, 41.9028, 12.4964, b"P")

print(f"ASC: {ascmc[0]:.4f}")
print(f"MC:  {ascmc[1]:.4f}")
for i, cusp in enumerate(cusps, 1):
    print(f"House {i:2d}: {cusp:.4f}")
```

## Thread Safety

The global API uses mutable global state (compatible with PySwissEph).
For concurrent calculations, use `EphemerisContext`:

```python
from libephemeris import EphemerisContext
from libephemeris.constants import SUN

ctx = EphemerisContext()
ctx.set_topo(12.5, 41.9, 0)
pos, _ = ctx.calc_ut(2451545.0, SUN, 0)
```

Each context has independent state (sidereal mode, topographic position, LEB reader).

## Minor Bodies

For asteroids and TNOs, JPL SPK kernels provide high precision.
Auto-download is supported:

```python
swe.set_auto_spk_download(True)
pos, _ = swe.calc_ut(2460000.0, swe.CHIRON, 0)
```

Fallback chain: SPK kernel → auto-download → strict precision check → REBOUND/ASSIST → Keplerian propagation.

With strict precision mode (enabled by default), `SPKRequiredError` is raised for mapped bodies when no SPK is available rather than silently falling back to low-precision Keplerian. Use `set_strict_precision(False)` to allow the full fallback chain.

See [Optional Modules and Calculation Backends](optional-modules.md) for details on the fallback chain, ASSIST n-body integration, and known limitations (e.g., Bennu).
