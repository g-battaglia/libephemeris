# Chapter 4 — Ephemerides: what they are and how they work

## What you will learn

In this chapter, you will discover what is "inside" an ephemeris, where the data used by the library comes from, the three available precision tiers, and the difference between geocentric, heliocentric, barycentric, and topocentric positions.

---

## 4.1 What is an ephemeris

An **ephemeris** is a table that says "on this date, this celestial body is located here". The concept is as old as civilization: 3000-year-old Babylonian clay tablets contained lunar ephemerides to predict eclipses. The Nautical Almanac, published since 1767, provides stellar positions to sailors for navigation.

A typical printed ephemeris has a row for each day and a column for each planet. A digital ephemeris is the same idea, but with much greater precision and the ability to interpolate between tabulated points to obtain the position at any given instant.

Modern ephemerides are not simple tables: they are binary files containing
**Chebyshev polynomials** that interpolate a numerically integrated trajectory.
Their accuracy depends on the source, body, date, and interpolation budget; an
evaluated value is not literally an error-free physical position.

---

## 4.2 JPL DE440: our data source

LibEphemeris uses the ephemerides produced by NASA's **Jet Propulsion Laboratory** (JPL) — the same laboratory that guides space probes. JPL ephemerides are the gold standard: they are produced by numerically integrating the equations of motion for all bodies in the Solar System, taking into account hundreds of gravitational and relativistic effects.

The library supports three ephemeris files:

| File | Name | Range | Size |
|------|------|-----------|-----------|
| `de440s.bsp` | DE440 short | 1849 – 2150 | ~31 MB |
| `de440.bsp` | DE440 | 1550 – 2650 | ~114 MB |
| `de441.bsp` | DE441 | -13200 – +17191 | ~3.1 GB |

**DE440** (2020) is the reference ephemeris: it covers 1100 years with sub-milliarcsecond precision for the inner planets. **DE441** is the extended version up to 30000 years, useful for historical research but with gradually lower precision for very distant dates.

Internally, LibEphemeris uses **Skyfield** to read the JPL files and calculate positions. Skyfield is a high-quality Python astronomy library developed by Brandon Rhodes.

---

## 4.3 The three precision tiers

The choice of precision tier determines which ephemeris file is used:

- **`base`** (DE440s): lightweight (~31 MB), covers 1849–2150. Ideal for modern astrology and mobile applications.
- **`medium`** (DE440): balanced (~114 MB), covers 1550–2650. The **default**. Good for almost all uses.
- **`extended`** (DE441): maximum (~3.1 GB), covers -13200 to +17191. For historical research and ancient/future dates.

In reviewed LEB mode the tier is cumulative rather than a forced downgrade:
base is tried first for each body/date, then medium, then extended. The
configured maximum therefore expands coverage while preserving the narrower
DE440s/DE440 source wherever it applies.

```python
import libephemeris as ephem

# View the current tier
# The default is "medium"

# Change tier
ephem.set_precision_tier("extended")  # for ancient dates

# Download the required files (one-time)
ephem.download_for_tier("medium")

# Revert to default
ephem.set_precision_tier("medium")
```

### Which one to choose?

For **modern astrology** (natal charts from 1900 onwards): `base` is sufficient. For **historical astrology** (before 1849): at least `medium` is needed. For **ancient historical research** (Babylonian eclipses, medieval comets): `extended` is required.

---

## 4.4 Geocentric, heliocentric, and barycentric

The calculated positions depend on **where** you imagine observing from:

- **Geocentric** (default): viewed from the center of the Earth. It is the viewpoint of astrology — the sky as it appears from Earth.
- **Heliocentric**: viewed from the center of the Sun. Flag `FLG_HELCTR`. Useful for celestial mechanics.
- **Barycentric**: viewed from the barycenter (center of mass) of the Solar System. Flag `FLG_BARYCTR`. The Sun is not exactly at the barycenter: it "wobbles" by about 2 solar radii due to the gravitational pull of the giant planets.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Geocentric Mars (default)
geo, _ = ephem.calc_ut(jd, ephem.MARS, 0)

# Heliocentric Mars
helio, _ = ephem.calc_ut(jd, ephem.MARS, ephem.FLG_HELCTR)

print(f"Geocentric Mars:  {geo[0]:.4f}° (dist {geo[2]:.4f} AU)")
print(f"Heliocentric Mars: {helio[0]:.4f}° (dist {helio[2]:.4f} AU)")
```

```
Geocentric Mars:  342.8458° (dist 2.0618 AU)
Heliocentric Mars: 317.5543° (dist 1.3879 AU)
```

### Planetocentric: observing from another planet

With `calc_pctr` you can calculate the position of a body as seen from another planet:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# The Sun as seen from Jupiter
pos, _ = ephem.calc_pctr(jd, ephem.SUN, ephem.JUPITER, 0)
print(f"Sun seen from Jupiter: {pos[0]:.4f}°, dist {pos[2]:.4f} AU")
```

```
Sun seen from Jupiter: 234.7096°, dist 5.0060 AU
```

---

## 4.5 Topocentric: location matters

The **geocentric** position is calculated from the center of the Earth. But you are not at the center of the Earth — you are on its surface. The **topocentric** position takes your exact location into account.

For most celestial bodies, the difference is negligible (less than 1"). But for the **Moon**, which is very close, the difference can reach ~1° — the so-called **lunar parallax**. This means that a solar eclipse is total in one place and partial in another just a few hundred km away.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 18.0)

# Set the observer's location (Rome)
ephem.set_topo(12.4964, 41.9028, 50.0)

# Geocentric Moon (default)
moon_geo, _ = ephem.calc_ut(jd, ephem.MOON, 0)

# Topocentric Moon
moon_topo, _ = ephem.calc_ut(jd, ephem.MOON, ephem.FLG_TOPOCTR)

diff = moon_topo[0] - moon_geo[0]
print(f"Geocentric Moon:  {moon_geo[0]:.4f}°")
print(f"Topocentric Moon: {moon_topo[0]:.4f}°")
print(f"Parallax:         {diff * 3600:.1f}\" arcseconds")
```

```
Geocentric Moon:  19.1832°
Topocentric Moon: 18.2383°
Parallax:         -3401.7" arcseconds
```

### 🌍 Real life

The solar eclipse of April 8, 2024, was total in Dallas (Texas) but only partial in New York — the difference is entirely due to parallax: the Moon covers the Sun differently depending on where you are on the Earth's surface.

---

## Summary

- An **ephemeris** is a table of celestial positions. LibEphemeris uses NASA JPL DE440/DE441 ephemerides, read via Skyfield.
- Three tiers: `base` (1849–2150), `medium` (1550–2650, default), `extended` (-13200 to +17191).
- **Geocentric** (default): from the center of the Earth. **Heliocentric** (`FLG_HELCTR`): from the Sun. **Barycentric** (`FLG_BARYCTR`): from the center of mass of the Solar System.
- **Topocentric** (`FLG_TOPOCTR`): from your location on the Earth's surface. Crucial for the Moon (~1° of parallax).

### Functions and constants introduced

| Function / Constant | Use |
|---------------------|-----|
| `set_precision_tier(tier)` | Chooses `"base"`, `"medium"`, or `"extended"` |
| `download_for_tier(tier)` | Downloads the ephemeris files |
| `set_topo(lon, lat, alt)` | Sets the observer's location |
| `calc_pctr(jd, body, center, flag)` | Position as seen from another planet |
| `FLG_HELCTR` | Heliocentric coordinates |
| `FLG_BARYCTR` | Barycentric coordinates |
| `FLG_TOPOCTR` | Topocentric coordinates |
