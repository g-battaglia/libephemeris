# Planet Centers SPK Generation

LibEphemeris generates a compact SPK (SPICE kernel) file containing planet center positions for the outer planets, correcting for the barycenter-to-center offset inherent in JPL Development Ephemerides.

## Table of Contents

- [Background](#background)
  - [Barycenters vs Planet Centers](#barycenters-vs-planet-centers)
  - [Barycenter Offsets](#barycenter-offsets)
- [Method](#method)
  - [JPL center segments](#jpl-center-segments)
  - [Explicit system-barycenter fallback](#explicit-system-barycenter-fallback)
  - [SPK File Format](#spk-file-format)
  - [Chaining with DE Ephemerides](#chaining-with-de-ephemerides)
- [SPK Generation Process](#spk-generation-process)
  - [Prerequisites](#prerequisites)
  - [Source Files](#source-files)
  - [Running the Script](#running-the-script)
  - [Output](#output)
- [Usage in LibEphemeris](#usage-in-libephemeris)
- [Precision and Validation](#precision-and-validation)
- [Troubleshooting](#troubleshooting)
- [References](#references)
- [Changelog](#changelog)

## Background

### Barycenters vs Planet Centers

JPL's Development Ephemerides (DE440, DE441, etc.) provide positions for planetary system **barycenters** rather than planet **centers** for the outer planets. This is because these planets have significant moon systems whose combined mass affects the barycenter location.

| NAIF ID | Body | Description |
|---------|------|-------------|
| 5 | Jupiter Barycenter | Center of mass of Jupiter system |
| 6 | Saturn Barycenter | Center of mass of Saturn system |
| 7 | Uranus Barycenter | Center of mass of Uranus system |
| 8 | Neptune Barycenter | Center of mass of Neptune system |
| 9 | Pluto Barycenter | Center of mass of Pluto-Charon system |
| 599 | Jupiter | Center of Jupiter itself |
| 699 | Saturn | Center of Saturn itself |
| 799 | Uranus | Center of Uranus itself |
| 899 | Neptune | Center of Neptune itself |
| 999 | Pluto | Center of Pluto itself |

### Barycenter Offsets

The barycenter can be significantly offset from the planet center:

| Planet | Typical Offset | Angular Size (from Earth) |
|--------|----------------|--------------------------|
| Jupiter | ~64 km | ~0.02 arcsec |
| Saturn | ~281 km | ~0.03 arcsec |
| Uranus | ~43 km | ~0.003 arcsec |
| Neptune | ~74 km | ~0.01 arcsec |
| Pluto | ~2131 km | ~0.15 arcsec |

For high-precision work, these offsets matter.

## Method

LibEphemeris uses JPL state data when a physical planet-center segment covers
the requested epoch. It does not synthesize a missing center from analytical
satellite theories.

This method applies to `auto`, `skyfield`, and explicit JPL workflows. In
sealed `leb` mode the planet-center kernel is deliberately disabled: the
runtime returns the Jupiter-through-Pluto system barycentres already stored in
the LEB core and never opens an auxiliary BSP.

### JPL center segments

The generator extracts planet-center segments from JPL satellite SPK files and
creates a compact file containing only those centers. At runtime the center is
resolved at the retarded epoch, so heliocentric and barycentric light time uses
the observer-to-body-center vector.

### Explicit system-barycenter fallback

When no JPL center segment covers the requested epoch, the stored planetary
system barycenter is returned as the target. This is a deliberate, traceable
fallback; no analytical COB correction is applied.

### SPK File Format

The output file uses SPK Type 2 (Chebyshev polynomial) format, the same format used in DE ephemerides and satellite SPKs. This format stores polynomial coefficients that can be evaluated at any time within the segment's coverage.

Type 2 SPK structure per segment:
- Chebyshev polynomial coefficients for X, Y, Z position
- Polynomial degree (typically 7-15)
- Record interval (typically 1-8 days)
- Time range (typically 1900-2100)

### Chaining with DE Ephemerides

Skyfield automatically chains SPK segments. When both DE440 and the planet centers SPK are loaded, Skyfield computes:

```
Earth -> SSB -> Jupiter Barycenter (from DE440)
                     |
              Jupiter Center (from planet_centers.bsp)
```

This chaining is transparent to the user. All segments use the J2000 reference frame (ICRF), consistent with DE ephemerides.

## SPK Generation Process

### Prerequisites

1. **spiceypy** >= 6.0.0 (Python wrapper for NAIF SPICE toolkit)
2. **Internet connection** (to download source SPK files)
3. **~6.5 GB temporary disk space** for the medium/extended source cache

### Source Files

`scripts/generate_planet_centers_spk.py` declares every input as an HTTPS URL
under the official JPL NAIF generic-kernel archive and extracts only the
listed barycenter-to-center segment with SPICE `spksub`. It never reads a
compatibility distribution or reference-API output.

| Tier | Jupiter | Saturn | Uranus | Neptune | Pluto |
|------|---------|--------|--------|---------|-------|
| base | `jup204.bsp` | `sat319.bsp` | `ura083.bsp` | `nep050.bsp` | `plu017.bsp` |
| medium | `jup365.bsp` | `sat441xl_part-1.bsp` + `sat441xl_part-2.bsp` | `ura111xl-799.bsp` | `nep097xl-899.bsp` | `plu060.bsp` |
| extended | `jup365.bsp` | `sat441xl_part-1.bsp` + `sat441xl_part-2.bsp` | `ura111xl-799.bsp` | `nep097xl-899.bsp` | `plu060.bsp` |

The base files are in NAIF's `a_old_versions` satellite directory; the medium
and extended files are in the current satellite directory. All downloads use
certificate-verified TLS and require a structurally valid BSP; medium/extended
inputs additionally require the reviewed source hashes recorded below. Every
matching DAF descriptor is copied: JPL products can split one target into
consecutive segments, so copying only the first descriptor truncates the
published coverage. The output is checked by reopening it and enumerating its
five expected targets.

The reviewed medium build used these exact JPL source bytes:

| Source | SHA-256 |
|--------|---------|
| `jup365.bsp` | `dbf016c01ba4d022154838000cf3f06962cf958ddc503a366f7fe8f81495c5cb` |
| `sat441xl_part-1.bsp` | `45e6f6687cf9a4d7c94c206c13650d35583740019df939c19f2ee0f0f7d2750d` |
| `sat441xl_part-2.bsp` | `d021819387e6045e2cd2819a7907ca62cc8cb3091ca925dfe894d2e0f4f7af1e` |
| `ura111xl-799.bsp` | `fa55665eb129613574a5d029bf2cc7b521d3f0d2edc7de5790add3cf950714bf` |
| `nep097xl-899.bsp` | `95e815cb8584e530d3c83f35e39af4bf67d0874f01b65538cbe9dbf279818df8` |
| `plu060.bsp` | `dfbb102491a26ed41ae08ca3f8963f22f0219df1d8f265ab87b9ad825a826fc6` |
| `naif0012.tls` | `678e32bdb5a744117a467cd9601cd6b373f0e9bc9bbde1371d5eee39600a039b` |

### Running the Script

```bash
# Using leph CLI (recommended)
leph generate planet-centers-medium

# Other tiers
leph generate planet-centers-base
leph generate planet-centers-extended
leph generate planet-centers-all       # All three tiers
```

### Output

The script generates one tier-specific file in `~/.libephemeris/` (the
default output directory):

```
~/.libephemeris/planet_centers_base.bsp
~/.libephemeris/planet_centers_medium.bsp
~/.libephemeris/planet_centers_extended.bsp
```

At runtime the loader looks for `planet_centers_{tier}.bsp` matching the
active precision tier. In the base tier only, it can fall back to the legacy
single-file name `planet_centers.bsp` when that file matches the pinned base
artifact. Each file contains five center targets; a target can have several
time-partitioned SPK segments:

| Center | Target | Description |
|--------|--------|-------------|
| 5 | 599 | Jupiter Barycenter -> Jupiter Center |
| 6 | 699 | Saturn Barycenter -> Saturn Center |
| 7 | 799 | Uranus Barycenter -> Uranus Center |
| 8 | 899 | Neptune Barycenter -> Neptune Center |
| 9 | 999 | Pluto Barycenter -> Pluto Center |

### Distributed output pins

The downloader accepts only these SHA-256-pinned outputs:

| File | SHA-256 |
|------|---------|
| `planet_centers_base.bsp` | `a9ec744ff412b095129166587ea0814f81c850faebf92586a738cb5dc103c92a` |
| `planet_centers_medium.bsp` | `d3c34f5efe9223ef588ec59a8c59c1bd6619b0eab5d5e0b35c353d675efe7b4d` |
| `planet_centers_extended.bsp` | `a07b046b89a9992fc7fda445b00e656341a3bab66a035adb8108de7d4bd69edc` |

For backward API compatibility, `download_planet_centers()` writes the pinned
base-tier bytes under the legacy destination name `planet_centers.bsp`; it no
longer accesses the former unverified legacy release object.

For the medium output, 50 position/velocity evaluations spanning the endpoints
and midpoints of all ten copied descriptors were compared directly with the
pinned JPL inputs: the maximum difference was exactly `0 km` and `0 km/day`.

## Usage in LibEphemeris

After generating the SPK file, LibEphemeris automatically loads it alongside
the main ephemeris. Covered outer-planet requests use the JPL center segment;
outside its coverage they use the explicit system-barycenter fallback.

### Automatic Loading

```python
import libephemeris as eph

# The pinned tier-specific file is loaded automatically if present.
# Base also accepts the pinned legacy planet_centers.bsp destination.
# Jupiter uses its JPL center segment when the epoch is covered
pos, _ = eph.calc_ut(jd, eph.JUPITER, 0)
```

### Verifying Planet Centers

```python
from skyfield.api import load

# Load both ephemerides
planets = load('de440.bsp')
centers = load('planet_centers.bsp')

# Get Jupiter center via chaining
ts = load.timescale()
t = ts.utc(2025, 1, 1)

earth = planets['earth']
jupiter_bary = planets['jupiter barycenter']
jupiter_center = centers['jupiter']

# Position of Jupiter center from Earth
astrometric = earth.at(t).observe(jupiter_bary + jupiter_center)
ra, dec, distance = astrometric.radec()
```

## Precision and Validation

### Time Coverage

The coverage depends on the source SPK files. A regenerated medium asset has:

| Planet | Start | End |
|--------|-------|-----|
| Jupiter | 1600 | 2200 |
| Saturn | before 1550 | after 2650 |
| Uranus | before 1550 | after 2650 |
| Neptune | before 1550 | after 2650 |
| Pluto | 1800 | 2200 |

JPL's published `jup365` comments define its 1600--2200 interval; the XL
products supply the wider Saturn/Uranus/Neptune intervals. See the official
[JPL NAIF satellite-kernel archive](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/).
No published JPL planet-center SPK currently fills the remaining Jupiter and
Pluto gaps, so runtime uses the system barycenter there.

### Position Accuracy

The extracted segments maintain the full precision of the source files:
- Position accuracy: sub-kilometer
- From Earth: sub-arcsecond (typically <0.001 arcsec)

## Troubleshooting

### SSL Certificate Errors

If SSL errors occur when downloading from JPL:

```
ssl.SSLError: certificate verify failed
```

The generator requires certificate-verified HTTPS through the certifi trust
store and fails closed on TLS errors. Check the local CA/certifi installation
or download the declared JPL sources manually, then verify their recorded
SHA-256 values before generation.

### spiceypy Not Found

```
ModuleNotFoundError: No module named 'spiceypy'
```

Install spiceypy:
```bash
pip install "spiceypy>=6.0.0"
# or
uv pip install "spiceypy>=6.0.0"
```

### Disk Space

Base generation needs roughly 500 MB of temporary space. Medium and extended
generation need about 6.5 GB for the full JPL satellite-kernel source set.
Temporary inputs are cleaned up after extraction unless an explicit cache is
requested.

## References

- [JPL NAIF Generic Kernels](https://naif.jpl.nasa.gov/naif/data_generic.html)
- [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
- [Skyfield SPK Support](https://rhodesmill.org/skyfield/planets.html)
- [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit.html)

## Changelog

- **2024**: Initial implementation
  - Created `generate_planet_centers_spk.py` script
  - Added CLI command (`leph generate planet-centers-medium`)
