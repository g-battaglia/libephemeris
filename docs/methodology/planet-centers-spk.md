# Planet Centers SPK Generation

LibEphemeris generates a compact SPK (SPICE kernel) file containing planet center positions for the outer planets, correcting for the barycenter-to-center offset inherent in JPL Development Ephemerides.

## Table of Contents

- [Background](#background)
  - [Barycenters vs Planet Centers](#barycenters-vs-planet-centers)
  - [Barycenter Offsets](#barycenter-offsets)
- [Method](#method)
  - [Analytical COB Correction (Current Approach)](#analytical-cob-correction-current-approach)
  - [Compact SPK Extraction (New Approach)](#compact-spk-extraction-new-approach)
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

LibEphemeris employs a three-tier correction strategy to convert barycenter positions to planet center positions.

### Analytical COB Correction (Current Approach)

The current approach uses analytical moon theories to compute the Center of Body (COB) offset from the barycenter:

- **Jupiter**: E5 theory (Galilean moons)
- **Saturn**: TASS 1.7 (Titan-dominated)
- **Neptune**: Triton Keplerian elements
- **Pluto**: Charon two-body solution
- **Uranus**: Keplerian elements for the 5 major moons (Titania/Oberon dominate)

This approach requires no extra files but has limited precision (~0.02-0.15 arcsec).

### Compact SPK Extraction (New Approach)

The new approach extracts planet center segments from JPL's satellite SPK files and creates a compact file containing only the planet centers. This provides maximum precision (<0.001 arcsec) with minimal file size (~5-10 MB).

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
3. **~500 MB temporary disk space** (source files are deleted after extraction)

### Source Files

The script downloads these satellite SPK files from JPL NAIF:

| File | Planet | Size | URL |
|------|--------|------|-----|
| jup204.bsp | Jupiter | 89 MB | [NAIF](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/jup204.bsp) |
| sat319.bsp | Saturn | 52 MB | [NAIF](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/sat319.bsp) |
| ura083.bsp | Uranus | 80 MB | [NAIF](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/ura083.bsp) |
| nep050.bsp | Neptune | 201 MB | [NAIF](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/nep050.bsp) |
| plu017.bsp | Pluto | 19 MB | [NAIF](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/plu017.bsp) |

These are older, more compact versions of the satellite ephemerides. Newer versions (jup365.bsp, etc.) are much larger (500 MB-1 GB each) but offer higher precision for satellite positions, which is not needed for planet centers.

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
active precision tier, falling back to the legacy single-file name
`planet_centers.bsp`. Each file contains 5 segments:

| Center | Target | Description |
|--------|--------|-------------|
| 5 | 599 | Jupiter Barycenter -> Jupiter Center |
| 6 | 699 | Saturn Barycenter -> Saturn Center |
| 7 | 799 | Uranus Barycenter -> Uranus Center |
| 8 | 899 | Neptune Barycenter -> Neptune Center |
| 9 | 999 | Pluto Barycenter -> Pluto Center |

## Usage in LibEphemeris

After generating the SPK file, LibEphemeris automatically loads it alongside the main ephemeris (DE440). The planet center positions are used instead of barycenter + COB correction.

### Automatic Loading

```python
import libephemeris as eph

# planet_centers.bsp is loaded automatically if present
# Jupiter now returns planet center, not barycenter
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

The coverage depends on the source SPK files, typically:

| Planet | Start | End |
|--------|-------|-----|
| Jupiter | ~1900 | ~2100 |
| Saturn | ~1900 | ~2100 |
| Uranus | ~1900 | ~2100 |
| Neptune | ~1900 | ~2100 |
| Pluto | ~1900 | ~2100 |

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

The script automatically falls back to unverified SSL for the trusted JPL servers. If this fails, the source files can be downloaded manually.

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

The script requires ~500 MB temporary space for downloading source files. These are automatically cleaned up after extraction.

## References

- [JPL NAIF Generic Kernels](https://naif.jpl.nasa.gov/naif/data_generic.html)
- [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
- [Skyfield SPK Support](https://rhodesmill.org/skyfield/planets.html)
- [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit.html)

## Changelog

- **2024**: Initial implementation
  - Created `generate_planet_centers_spk.py` script
  - Added CLI command (`leph generate planet-centers-medium`)
