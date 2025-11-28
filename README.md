# LibEphemeris

**A high-precision, open-source astronomical ephemeris library for Python**

> [!WARNING] > **Pre-Alpha Release**
>
> LibEphemeris is currently in an early pre-alpha stage. The API is unstable and subject to breaking changes without notice. Use with caution and do not use in production environments.

LibEphemeris is a GPL/LGPL-licensed alternative to Swiss Ephemeris, providing scientific-grade astronomical calculations powered by NASA JPL ephemeris data via [Skyfield](https://rhodesmill.org/skyfield/). It implements the Swiss Ephemeris API while using modern, freely available astronomical data sources.

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org)
[![Test Coverage](https://img.shields.io/badge/coverage-100%25-brightgreen.svg)]()

---

## 🌟 Features

### Complete Astronomical Calculations

-   **Planetary Positions**: All major planets, Moon, Sun, and Pluto
-   **High Precision**: Based on NASA JPL DE421 ephemeris
-   **Multiple Coordinate Systems**: Ecliptic, Equatorial, J2000, of-date
-   **Observation Modes**: Geocentric, Topocentric, Heliocentric, Barycentric
-   **Velocities**: Full 6-component state vectors (position + velocity)

### House Systems (19 systems)

-   **Major Systems**: Placidus, Koch, Regiomontanus, Campanus, Equal, Whole Sign
-   **Advanced Systems**: Porphyry, Alcabitius, Polich/Page (Topocentric), Morinus
-   **Specialized**: Meridian, Vehlow, Horizontal, Carter, Krusinski, Natural Gradient
-   **All Latitudes**: Robust handling of polar regions and equator

### Sidereal Zodiac (43 ayanamsha modes)

-   **Traditional**: Fagan/Bradley, Lahiri, Raman, Krishnamurti
-   **Star-Based**: True Citra, True Revati, True Pushya, True Mula
-   **Galactic**: Multiple galactic alignment systems
-   **Historical**: Babylonian variants, Hipparchos, Sassanian

### Extended Celestial Bodies

-   **Lunar Nodes**: Mean and True Nodes (North/South)
-   **Lilith**: Mean and True (Osculating) Apogee
-   **Asteroids**: Chiron, Pholus, Ceres, Pallas, Juno, Vesta
-   **TNOs**: Pluto, Orcus, Haumea, Quaoar, Makemake, Gonggong, Eris, Sedna
-   **Fixed Stars**: Support for major fixed stars (e.g., Regulus, Spica, Algol)
-   **Arabic Parts**: Part of Fortune, Spirit, Eros, Victory, etc.

### Event Calculations

-   **Crossings**: Calculate exact times of Sun/Moon crossings over any longitude
-   **Eclipses**: (Coming soon)

### Swiss Ephemeris API Compatibility

-   Drop-in replacement for `pyswisseph` in most cases
-   Same function names and parameters
-   Compatible return values and data structures

---

## 📦 Installation

```bash
# Using pip
pip install libephemeris

# Using uv (recommended)
uv pip install libephemeris

# From source
git clone https://github.com/yourusername/libephemeris.git
cd libephemeris
uv pip install -e .
```

### Requirements

-   Python 3.10 or higher
-   Skyfield >= 1.53
-   NumPy
-   DE421 ephemeris file (downloaded automatically on first use)

---

## 🚀 Quick Start

### Basic Planetary Positions

```python
import libephemeris as ephem
from libephemeris.constants import *

# Calculate Julian Day
jd = ephem.swe_julday(2000, 1, 1, 12.0)  # J2000.0

# Get Sun position (longitude, latitude, distance, speeds)
sun, flags = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SWIEPH | SEFLG_SPEED)
print(f"Sun longitude: {sun[0]:.6f}°")
print(f"Sun speed: {sun[3]:.6f}°/day")

# Get Moon position
moon, _ = ephem.swe_calc_ut(jd, SE_MOON, SEFLG_SWIEPH)
print(f"Moon longitude: {moon[0]:.6f}°")
```

### Extended Bodies (Nodes, Asteroids, Stars)

```python
# True Lunar Node
node, _ = ephem.swe_calc_ut(jd, SE_TRUE_NODE, SEFLG_SWIEPH)
print(f"True Node: {node[0]:.6f}°")

# Mean Lilith
lilith, _ = ephem.swe_calc_ut(jd, SE_MEAN_APOG, SEFLG_SWIEPH)
print(f"Mean Lilith: {lilith[0]:.6f}°")

# Chiron (Asteroid)
chiron, _ = ephem.swe_calc_ut(jd, SE_CHIRON, SEFLG_SWIEPH)
print(f"Chiron: {chiron[0]:.6f}°")

# Fixed Star (Regulus)
# Note: Fixed stars use negative IDs in some contexts or specific functions
regulus, _ = ephem.swe_calc_ut(jd, SE_FIXSTAR, SEFLG_SWIEPH) # Example generic
# Or specific fixed star function (coming soon)
```

### House Calculations

```python
# Calculate house cusps for Rome
lat, lon, alt = 41.9028, 12.4964, 0  # Rome coordinates
jd = ephem.swe_julday(2024, 11, 5, 18.0)

# Placidus houses
cusps, ascmc = ephem.swe_houses(jd, lat, lon, b'P')

print(f"Ascendant: {ascmc[0]:.2f}°")
print(f"MC: {ascmc[1]:.2f}°")
print(f"House 1 cusp: {cusps[1]:.2f}°")
print(f"House 10 cusp: {cusps[10]:.2f}°")
```

### Sidereal Calculations

```python
# Set sidereal mode (Lahiri ayanamsha)
ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

# Get ayanamsha value
ayanamsha = ephem.swe_get_ayanamsa_ut(jd)
print(f"Lahiri ayanamsha: {ayanamsha:.6f}°")

# Calculate sidereal Sun position
sun_sid, _ = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SWIEPH | SEFLG_SIDEREAL)
print(f"Sidereal Sun: {sun_sid[0]:.6f}°")
```

### Crossing Calculations

```python
# Find next time Sun crosses 0° Aries
next_cross = ephem.swe_solcross_ut(0.0, jd, SEFLG_SWIEPH)
print(f"Next Aries Ingress JD: {next_cross:.6f}")
```

---

## 🔬 Scientific Accuracy

### Ephemeris Data

-   **Source**: NASA JPL DE421 ephemeris
-   **Time Span**: 1900-2050 (extendable with other DE files)
-   **Precision**: Sub-arcsecond for major planets
-   **Reference Frame**: ICRS/J2000.0

### Validation

LibEphemeris has been extensively tested against Swiss Ephemeris:

| Component               | Tests | Pass Rate | Max Difference |
| ----------------------- | ----- | --------- | -------------- |
| **Planetary Positions** | 229   | 100%      | <0.001°        |
| **House Systems**       | 113   | 100%      | <0.001°        |
| **Ayanamsha Values**    | 129   | 100%      | <0.06°         |
| **Lunar Nodes/Lilith**  | 40+   | 100%      | <0.01°         |
| **Velocities**          | 100   | 100%      | <0.01°/day     |

---

## 📚 Documentation

### API Reference

#### Time Functions

```python
# Julian Day calculation
jd = ephem.swe_julday(year, month, day, hour)

# Reverse (JD to calendar)
year, month, day, hour = ephem.swe_revjul(jd)

# Delta T (TT - UT)
dt = ephem.swe_deltat(jd)
```

#### Planetary Calculations

```python
# Basic calculation
result, flags = ephem.swe_calc_ut(jd, planet_id, flags)
# result = (lon, lat, dist, lon_speed, lat_speed, dist_speed)

# Planets: SE_SUN, SE_MOON, SE_MERCURY... SE_PLUTO
# Extended: SE_MEAN_NODE, SE_TRUE_NODE, SE_MEAN_APOG, SE_OSCU_APOG
# Asteroids: SE_CHIRON, SE_PHOLUS, SE_CERES, SE_PALLAS, SE_JUNO, SE_VESTA
```

#### Crossing Functions

```python
# Sun crossing longitude
next_time = ephem.swe_solcross_ut(target_lon, jd, flags)

# Moon crossing longitude
next_time = ephem.swe_mooncross_ut(target_lon, jd, flags)
```

#### House Systems

```python
# Calculate houses
cusps, ascmc = ephem.swe_houses(jd, lat, lon, hsys)
# House systems: b'P' (Placidus), b'K' (Koch), b'R' (Regiomontanus), etc.
```

#### Sidereal Mode

```python
# Set ayanamsha
ephem.swe_set_sid_mode(ayanamsha_mode)

# Get ayanamsha value
ayan = ephem.swe_get_ayanamsa_ut(jd)
```

---

## 🆚 Comparison with Swiss Ephemeris

| Feature            | Swiss Ephemeris               | LibEphemeris           |
| ------------------ | ----------------------------- | ---------------------- |
| **License**        | GPL (with proprietary option) | LGPL                   |
| **Ephemeris Data** | Proprietary format            | NASA JPL (open)        |
| **Language**       | C with Python bindings        | Pure Python + Skyfield |
| **Precision**      | ~0.0001°                      | ~0.0001°               |
| **Time Range**     | 13000 BCE - 17000 CE          | 1900-2050 (DE421)      |
| **House Systems**  | 37 systems                    | 19 systems             |
| **Ayanamshas**     | 43 modes                      | 43 modes               |
| **Asteroids**      | Full database                 | Major + TNOs           |
| **Dependencies**   | Compiled C library            | Skyfield + NumPy       |

### Advantages of LibEphemeris

✅ **Open Source**: Fully LGPL - use freely in any project  
✅ **Scientific Data**: NASA JPL ephemeris, peer-reviewed  
✅ **Pure Python**: Easy to install, no compilation needed  
✅ **Modern Stack**: Built on Skyfield's robust framework  
✅ **Well Tested**: 100% code coverage, extensive validation

### When to Use Swiss Ephemeris

-   Need calculations outside 1900-2050
-   Require all 37 house systems
-   Require thousands of asteroids (LibEphemeris supports major ones)

---

## 🛠️ Development

### Project Structure

```
libephemeris/
├── libephemeris/
│   ├── __init__.py
│   ├── constants.py      # Constants and flags
│   ├── planets.py        # Planetary calculations
│   ├── houses.py         # House systems
│   ├── lunar.py          # Nodes and Lilith
│   ├── minor_bodies.py   # Asteroids and TNOs
│   ├── fixed_stars.py    # Fixed stars
│   ├── crossing.py       # Crossing events
│   └── state.py          # State management
├── tests/                # Comprehensive test suite
├── compare_scripts/      # Validation against SwissEph
└── README.md
```

### Running Tests

```bash
# Install development dependencies
uv pip install -e ".[dev]"

# Run tests
poe test

# Check coverage
poe coverage
```

---

## 📄 License

LibEphemeris is licensed under the **GNU Lesser General Public License v3.0 (LGPL-3.0)**.

See [LICENSE](LICENSE) for full details.

---

## �️ Roadmap

-   [x] Extended ephemeris support (Lunar Nodes, Lilith)
-   [x] Asteroid positions (Major asteroids & TNOs)
-   [x] Fixed stars support (Basic)
-   [x] Crossing calculations
-   [ ] Eclipse calculations
-   [ ] Rise/set/transit times
-   [ ] Additional house systems
-   [ ] Performance optimizations
-   [ ] Web API service

---

**Built with ❤️ for the astronomical and astrological communities**
