# Chapter 12 — Minor bodies: asteroids, centaurs, TNOs

## What you will learn

In this chapter, you will discover how the library handles thousands of celestial bodies beyond the main planets: main-belt asteroids, centaurs, trans-Neptunian objects (TNOs), and near-Earth asteroids. You will learn the "calculation chain" that the library uses to find the most precise position possible, and how to download and register SPK kernels for the bodies you are interested in.

---

## 12.1 The hierarchy of celestial bodies

Not all Solar System bodies are equal in terms of available precision:

**Planets (Mercury–Neptune)**: always available with sub-millimeter precision, calculated directly from the JPL DE440 ephemeris. You don't need to do anything special.

**Pluto and Chiron**: included in the JPL ephemeris like the planets. Same precision, same methods.

**Major asteroids (Ceres, Pallas, Juno, Vesta)**: have dedicated IDs (`CERES=17`, `PALLAS=18`, `JUNO=19`, `VESTA=20`). The library can download high-precision SPK kernels from NASA/JPL.

**Centaurs (Pholus, Nessus, Chariklo)**: bodies with orbits between Jupiter and Neptune. Available with downloadable SPK kernels or with the Keplerian fallback.

**Trans-Neptunians (Eris, Sedna, Haumea, Makemake, Quaoar)**: beyond Neptune. Same methods as the centaurs.

**Near-Earth (Apophis, Bennu, Eros)**: asteroids with orbits close to Earth. Particularly interesting for space missions.

```python
import libephemeris as ephem

# Dedicated IDs for the most important bodies
bodies = [
    (ephem.CHIRON,  "Chiron"),
    (ephem.PHOLUS,  "Pholus"),
    (ephem.CERES,   "Ceres"),
    (ephem.PALLAS,  "Pallas"),
    (ephem.JUNO,    "Juno"),
    (ephem.VESTA,   "Vesta"),
]

jd = ephem.julday(2024, 4, 8, 12.0)

signs = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

for body_id, name in bodies:
    pos, _ = ephem.calc_ut(jd, body_id, ephem.FLG_SPEED)
    sign = signs[int(pos[0] / 30)]
    degrees = pos[0] % 30
    print(f"{name:10s}  {degrees:5.1f}° {sign}")
```

```
Chiron       19.4° Ari
Pholus       10.3° Cap
Ceres        17.8° Cap
Pallas        8.2° Sgr
Juno          6.9° Vir
Vesta         2.4° Cnc
```

For TNOs and other bodies without a dedicated ID, you use the **catalog number + `AST_OFFSET`** (10000):

```python
import libephemeris as ephem

# Eris (number 136199)
ERIS = ephem.AST_OFFSET + 136199  # = 146199

jd = ephem.julday(2024, 4, 8, 12.0)
pos, _ = ephem.calc_ut(jd, ERIS, ephem.FLG_SPEED)
print(f"Eris: {pos[0]:.2f}°")
```

```
Eris: 24.74°
```

---

## 12.2 The calculation chain

When you request the position of a minor body, the library tries different methods in order, from most precise to least precise:

**1. SPK Kernel** — If a JPL binary file (SPK/BSP format) is registered for that body, it uses it. Precision: sub-arcsecond. It is the gold standard method.

**2. Automatic SPK download** — If automatic download is enabled and the kernel is not available locally, the library downloads it from JPL Horizons. Works for all 37 bodies in the SPK map.

**3. Strict precision check** — With strict precision mode (enabled by default), `SPKRequiredError` is raised for mapped bodies when no SPK is available, preventing silent precision loss. Use `set_strict_precision(False)` to allow lower-precision fallbacks.

**4. ASSIST n-body** — If `libephemeris[nbody]` is installed with data files, REBOUND/ASSIST provides sub-arcsecond precision for any body with orbital elements.

**5. Keplerian fallback** — Uses Kepler's laws with corrections for giant planet perturbations. Less precise (arcminutes over year scales), but works without an Internet connection.

```python
import libephemeris as ephem

# Enable automatic SPK download
ephem.set_auto_spk_download(True)

# Now calc_ut will download the SPK if necessary
jd = ephem.julday(2024, 4, 8, 12.0)
pos, _ = ephem.calc_ut(jd, ephem.CERES, ephem.FLG_SPEED)
print(f"Ceres (with automatic SPK): {pos[0]:.4f}°")
```

```
Ceres (with automatic SPK): 287.7759°
```

---

## 12.3 SPK Kernels: the gold standard

SPK (Spacecraft and Planet Kernel) kernels are NASA binary files containing precise trajectories in the form of Chebyshev polynomials. They are the same format used for the planets in the DE440 ephemeris.

### Downloading and registering an SPK

```python
import libephemeris as ephem

# Download the SPK for Ceres and register it
path = ephem.download_and_register_spk(
    "1;",             # identifier for JPL Horizons
    ephem.CERES,   # body ID in the library
    "2000-01-01",     # start date
    "2050-01-01",     # end date
)

print(f"Downloaded SPK: {path}")
```

```
Downloaded SPK: ~/.libephemeris/spk/1_200001_205001.bsp
```

### Managing registered SPKs

```python
import libephemeris as ephem

# Which bodies have an SPK loaded?
bodies = ephem.list_spk_bodies()
for body_id, (path, naif_id) in bodies.items():
    print(f"Body {body_id}: NAIF={naif_id}, file={path}")

# Does a specific body have the SPK?
if ephem.is_spk_available_for_body(ephem.CERES):
    print("Ceres: SPK available")
```

```
Body 17: NAIF=20000001, file=ceres_201603_203603.bsp
Ceres: SPK available
```

### "Major" asteroids with preconfigured SPK

The library knows the SPK download parameters for 37+ bodies. You can ensure an SPK is available with:

```python
import libephemeris as ephem

# Ensure Vesta's SPK is available
success = ephem.ensure_major_asteroid_spk(ephem.VESTA)
if success:
    print("Vesta's SPK ready")

# List all bodies with downloadable SPK
for body_id, name in ephem.list_major_asteroids():
    print(f"  {name} (ID: {body_id})")
```

```
Vesta's SPK ready
  Ceres (ID: 17)
  Pallas (ID: 18)
  Juno (ID: 19)
  Vesta (ID: 20)
  Chiron (ID: 15)
```

---

## 12.4 Searching for an asteroid by name or number

If you know the name of an asteroid but not its number, the library can search for it — first in the local database, then by querying the JPL SBDB (Small Body Database) service:

```python
import libephemeris as ephem

# Search by name
number = ephem.get_asteroid_number("Vesta")
print(f"Vesta = asteroid no. {number}")  # 4

number = ephem.get_asteroid_number("Apophis")
print(f"Apophis = asteroid no. {number}")  # 99942
```

```
Vesta = asteroid no. 4
Apophis = asteroid no. 99942
```

### Calculating by number

To calculate the position of any asteroid given its catalog number:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Calculate the position of Eros (no. 433)
lon, lat, dist = ephem.calc_asteroid_by_number(433, jd)
print(f"Eros: lon={lon:.2f}°, lat={lat:.2f}°, dist={dist:.4f} AU")
```

```
Eros: lon=91.13°, lat=5.97°, dist=1.1666 AU
```

The function automatically downloads the orbital elements from JPL SBDB if they are not already cached. You can also download them explicitly:

```python
import libephemeris as ephem

# Download the orbital elements of Apophis
elements = ephem.fetch_orbital_elements_from_sbdb(99942)
if elements:
    print(f"Name: {elements.name}")
    print(f"Semi-major axis: {elements.a:.4f} AU")
    print(f"Eccentricity: {elements.e:.6f}")
    print(f"Inclination: {elements.i:.4f}°")
```

> **Note**: this function requires Internet access to query the JPL SBDB (Small Body Database) service. The results are cached locally for subsequent requests.

---

## 12.5 The Keplerian fallback

When an SPK kernel is not available, the library calculates positions using **Kepler's laws** with corrections for giant planet perturbations.

The method works like this:

1. **Osculating orbital elements**: six numbers that describe the orbit's ellipse at a reference epoch (semi-major axis, eccentricity, inclination, ascending node, argument of perihelion, mean anomaly).

2. **Kepler's equation**: given the elapsed time, calculates where the body is located on the ellipse.

3. **Secular perturbations**: Jupiter, Saturn, Uranus, and Neptune slowly "push" the orbit, causing it to rotate and change shape. The library applies these corrections.

4. **Libration model**: for "plutinos" (bodies in a 2:3 resonance with Neptune, like Ixion and Orcus), the library corrects for the oscillation of the resonance argument.

The typical precision of the Keplerian fallback depends on the time elapsed from the epoch of the elements:

- **1 month**: ~7 arcseconds — excellent
- **1 year**: ~2 arcminutes — good for most uses
- **10 years**: ~30 arcminutes — acceptable for general purposes
- **50 years**: ~3.6° — orientation only

For modern astrology (dates from 1900 to 2100), the precision is generally sufficient. For research work or for dates far from the epoch of the elements, always use SPKs.

---

## Summary

In this chapter, we learned how to work with the minor bodies of the Solar System.

**Key concepts:**

- **Minor bodies** include asteroids, centaurs, TNOs, and near-Earth asteroids — thousands of bodies beyond the planets.
- The library uses a **calculation chain**: SPK kernel → auto-download → strict precision check → ASSIST n-body → Keplerian fallback.
- **SPK Kernels** are NASA binary files with precise trajectories — the gold standard for sub-arcsecond positions.
- **Strict precision** (default) raises `SPKRequiredError` for mapped bodies without SPK, preventing silent precision loss. Disable with `set_strict_precision(False)`.
- The **Keplerian fallback** works without the Internet but is only reached if strict precision is disabled or the body is not in the SPK map.
- For bodies without a dedicated ID, use `AST_OFFSET + catalog_number`.

**Introduced functions:**

- `calc_ut(jd, CERES, flag)` — calculates the position of a minor body with a dedicated ID, just like for the planets.
- `download_and_register_spk(body_id, jd_start, jd_end)` — downloads and registers an SPK kernel from JPL.
- `ensure_major_asteroid_spk(body_id)` — ensures the SPK is available, downloading it if necessary.
- `list_major_asteroids()` — lists the bodies with supported automatic SPK download.
- `list_spk_bodies()` — shows which bodies have a registered SPK.
- `set_auto_spk_download(True)` — enables automatic download of SPKs.
- `get_asteroid_number(name)` — looks up the catalog number of an asteroid by name.
- `calc_asteroid_by_number(number, jd)` — calculates the position of any asteroid given the number.
- `fetch_orbital_elements_from_sbdb(number)` — downloads orbital elements from JPL SBDB.
