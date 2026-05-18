# Chapter 3 — Where a planet is located: celestial coordinates

## What you will learn

In this chapter you will learn the three coordinate systems used in astronomy — ecliptic, equatorial and horizontal — and how to convert between them. You will also understand what reference frames are (J2000, ICRS), and the effects of precession, nutation and aberration.

---

## 3.1 Ecliptic coordinates: longitude and latitude

The most used coordinate system in astrology. The reference plane is the ecliptic (the Sun's path), the zero point is the vernal equinox.

- **Ecliptic longitude**: 0°–360° along the ecliptic. It is the "zodiac sign".
- **Ecliptic latitude**: distance from the ecliptic, from -90° to +90°.

The Sun always has a latitude of ~0° (by definition, it traces the ecliptic). The Moon can reach ±5.1° of latitude. The planets stay within a few degrees of the ecliptic, except Pluto which can reach ~17°.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)
pos, flag = ephem.calc_ut(jd, ephem.MOON, 0)

lon = pos[0]   # ecliptic longitude (degrees)
lat = pos[1]   # ecliptic latitude (degrees)
dist = pos[2]  # distance (AU)

sign_num = int(lon / 30)
signs = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]
print(f"Moon: {lon % 30:.2f}° {signs[sign_num]}, lat {lat:.2f}°")
# "Moon at 15° Cancer" = ecliptic longitude 105.x°
```

```
Moon: 15.43° Ari, lat -0.02°
```

---

## 3.2 Equatorial coordinates: right ascension and declination

The system used by telescopes. The reference plane is the celestial equator.

- **Right Ascension** (RA): 0°–360° along the celestial equator (sometimes expressed in hours: 0h–24h).
- **Declination** (Dec): distance from the celestial equator, from -90° to +90°.

Telescopes prefer this system because the celestial equator is aligned with the Earth's rotation: it is enough to rotate around a single axis to track an object.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Flag FLG_EQUATORIAL: returns RA and Dec instead of lon and lat
pos, flag = ephem.calc_ut(jd, ephem.MOON, ephem.FLG_EQUATORIAL)

ra = pos[0]    # Right Ascension in degrees (0–360)
dec = pos[1]   # Declination in degrees (-90 to +90)

# RA conversion from degrees to hours:minutes:seconds
ra_hours = ra / 15.0
ra_h = int(ra_hours)
ra_m = int((ra_hours - ra_h) * 60)
ra_s = ((ra_hours - ra_h) * 60 - ra_m) * 60

print(f"Moon: RA = {ra_h}h {ra_m}m {ra_s:.1f}s, Dec = {dec:+.2f}°")
```

```
Moon: RA = 0h 56m 52.2s, Dec = +6.06°
```

---

## 3.3 Horizontal coordinates: altitude and azimuth

The "practical" system — where should I look to see an object in the sky? This system depends on the observer's position and the time.

- **Altitude**: angle above the horizon. 0° = horizon, 90° = zenith. Negative values = below the horizon.
- **Azimuth**: direction on the horizon. In LibEphemeris: 0° = South, 90° = West, 180° = North, 270° = East.

The observer's position is required and, optionally, pressure and temperature for atmospheric refraction.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 9, 15, 21.0)

# Jupiter in equatorial coordinates
pos, _ = ephem.calc_ut(jd, ephem.JUPITER, ephem.FLG_EQUATORIAL)

# Rome (lon East, lat North, altitude in meters)
geopos = (12.4964, 41.9028, 50.0)

# From equatorial to horizontal (EQU2HOR = 1)
hor = ephem.azalt(jd, ephem.EQU2HOR, geopos, 1013.25, 15.0,
                  (pos[0], pos[1], pos[2]))

print(f"Azimuth: {hor[0]:.1f}° (from South)")
print(f"True altitude: {hor[1]:.1f}°")
print(f"Apparent altitude: {hor[2]:.1f}° (with refraction)")

# Inverse operation: from horizontal to equatorial
ra_dec = ephem.azalt_rev(jd, ephem.HOR2EQU, geopos, hor[0], hor[1])
print(f"Recovered RA: {ra_dec[0]:.4f}°, Dec: {ra_dec[1]:.4f}°")
```

```
Azimuth: 235.7° (from South)
True altitude: -3.2°
Apparent altitude: -3.2° (with refraction)
Recovered RA: 79.6564°, Dec: 22.3890°
```

### 🌍 Real life

"Is Jupiter visible tonight?" — calculate the altitude: if it is > 0°, it is above the horizon. If it is > 10°, it is high enough not to be disturbed by the atmosphere near the horizon. The azimuth tells you in which direction to look.

---

## 3.4 Reference frames: ICRS, J2000, equinox of date

The coordinates of a celestial object depend on the chosen **reference frame**. The zero point and axes can be defined in different ways:

- **Equinox of date**: the default in LibEphemeris. The coordinates take into account the current precession and nutation. The zero point is the vernal equinox *of that moment*. It is the system used in astrology.

- **J2000.0**: the coordinates are fixed to January 1, 2000, 12:00 TT. The zero point is the vernal equinox as it was in 2000. It is the standard for astronomical catalogs.

- **ICRS**: the modern system, based on the positions of very distant quasars (practically fixed). It is almost identical to J2000 — the difference ("frame bias") is ~0.02".

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# "Of date" ecliptic coordinates (default)
pos_date, _ = ephem.calc_ut(jd, ephem.MARS, 0)

# J2000 ecliptic coordinates
pos_j2000, _ = ephem.calc_ut(jd, ephem.MARS, ephem.FLG_J2000)

# ICRS equatorial coordinates
pos_icrs, _ = ephem.calc_ut(jd, ephem.MARS,
                            ephem.FLG_EQUATORIAL | ephem.FLG_ICRS)

print(f"Mars of-date: {pos_date[0]:.4f}°")
print(f"Mars J2000:   {pos_j2000[0]:.4f}°")
print(f"Difference:   {pos_date[0] - pos_j2000[0]:.4f}° (≈ precession)")
```

```
Mars of-date: 342.8458°
Mars J2000:   342.5082°
Difference:   0.3376° (≈ precession)
```

The difference between "of date" and J2000 is mainly the precession accumulated from 2000 to today — about 0.34° in 2024.

---

## 3.5 Precession, nutation and aberration

Three physical effects influence celestial coordinates. It is important to know they exist, even if in practice the library handles them automatically.

### Precession

The Earth's axis of rotation does not always point in the same direction: it draws a circle in the sky in about 26000 years, like a wobbling top. The effect is that the vernal point (0° Aries) shifts by about 50 arcseconds per year relative to the stars.

In 2000 years, precession has shifted the vernal point by about 28° — almost an entire zodiac sign. This is the reason why the tropical zodiac and constellations no longer coincide.

### Nutation

Superimposed on precession, the Earth's axis makes small oscillations with a main period of ~18.6 years (linked to the cycle of the lunar nodes). The amplitude is ~9" in longitude and ~17" in obliquity.

The `FLG_NONUT` flag calculates "mean" coordinates (without nutation) instead of "true" ones:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# With nutation (default = "true" coordinates)
pos_true, _ = ephem.calc_ut(jd, ephem.SUN, 0)

# Without nutation ("mean" coordinates)
pos_mean, _ = ephem.calc_ut(jd, ephem.SUN, ephem.FLG_NONUT)

diff_arcsec = (pos_true[0] - pos_mean[0]) * 3600
print(f"Nutation effect: {diff_arcsec:.2f}\"")
```

```
Nutation effect: -5.30"
```

### Aberration

Light takes time to arrive, and meanwhile the Earth moves. The effect is an apparent shift up to ~20" in the direction of the Earth's motion. The `FLG_NOABERR` flag disables it; `FLG_TRUEPOS` returns the pure geometric position (without correction for light-time nor for aberration).

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Apparent position (default: with aberration)
pos_app, _ = ephem.calc_ut(jd, ephem.MARS, 0)

# True geometric position (without aberration nor light-time)
pos_true, _ = ephem.calc_ut(jd, ephem.MARS, ephem.FLG_TRUEPOS)

diff = (pos_app[0] - pos_true[0]) * 3600
print(f"Aberration + light-time effect: {diff:.1f}\"")
```

```
Aberration + light-time effect: -33.3"
```

---

## 3.6 Cartesian coordinates (XYZ)

Sometimes positions are needed in Cartesian coordinates (X, Y, Z) instead of angular ones. Useful for 3D distance calculations, phase angles, or geometric problems.

The `FLG_XYZ` flag changes the meaning of the returned tuple: instead of (lon, lat, dist, vel_lon, vel_lat, vel_dist) you get (x, y, z, vx, vy, vz) in astronomical units and AU/day.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Equatorial Cartesian position
pos, _ = ephem.calc_ut(jd, ephem.MARS,
                       ephem.FLG_XYZ | ephem.FLG_EQUATORIAL)

x, y, z = pos[0], pos[1], pos[2]
import math
dist = math.sqrt(x**2 + y**2 + z**2)
print(f"Mars: X={x:.4f}, Y={y:.4f}, Z={z:.4f} AU")
print(f"Distance: {dist:.4f} AU = {dist * 149597870.7:.0f} km")
```

```
Mars: X=1.9696, Y=-0.5400, Z=-0.2829 AU
Distance: 2.0618 AU = 308436754 km
```

---

## 3.7 Converting between systems

LibEphemeris offers dedicated functions for conversions:

| Function | Conversion |
|----------|------------|
| `cotrans(coord, -obliquity)` | Ecliptic → Equatorial |
| `cotrans(coord, +obliquity)` | Equatorial → Ecliptic |
| `cotrans_sp(coord, speed, obl)` | As above, with speed |
| `azalt(jd, ECL2HOR, ...)` | Ecliptic → Horizontal |
| `azalt(jd, EQU2HOR, ...)` | Equatorial → Horizontal |
| `azalt_rev(jd, HOR2ECL, ...)` | Horizontal → Ecliptic |
| `azalt_rev(jd, HOR2EQU, ...)` | Horizontal → Equatorial |

The sign of the obliquity in `cotrans` controls the direction: **negative** = ecliptic→equatorial, **positive** = equatorial→ecliptic.

### Complete example: from ecliptic position to altitude in the sky

```python
import libephemeris as ephem

jd = ephem.julday(2024, 9, 15, 21.0)

# 1. Saturn's ecliptic position
pos, _ = ephem.calc_ut(jd, ephem.SATURN, 0)
lon, lat, dist = pos[0], pos[1], pos[2]

# 2. Get the obliquity
nut, _ = ephem.calc_ut(jd, ephem.ECL_NUT, 0)
obliquity = nut[0]

# 3. Ecliptic -> Equatorial (negative obliquity)
ra, dec, d = ephem.cotrans((lon, lat, dist), -obliquity)

# 4. Equatorial -> Horizontal (Rome)
geopos = (12.4964, 41.9028, 50.0)
hor = ephem.azalt(jd, ephem.EQU2HOR, geopos, 1013.25, 15.0,
                  (ra, dec, d))

print(f"Saturn: lon={lon:.2f}°, lat={lat:.2f}°")
print(f"        RA={ra:.2f}°, Dec={dec:+.2f}°")
print(f"        Altitude={hor[2]:.1f}°, Azimuth={hor[0]:.1f}° (from S)")
```

```
Saturn: lon=345.44°, lat=-2.20°
        RA=347.46°, Dec=-7.76°
        Altitude=35.5°, Azimuth=329.5° (from S)
```

> **Shortcut**: you can also pass ecliptic coordinates directly to `azalt` with `ECL2HOR`, avoiding manual conversion. The explicit conversion via `cotrans` is useful when you need intermediate values.

---

## Summary

- **Ecliptic** (lon, lat): system of astrology, based on the ecliptic. Default of `calc_ut`.
- **Equatorial** (RA, Dec): system of telescopes, based on the celestial equator. Flag `FLG_EQUATORIAL`.
- **Horizontal** (azimuth, altitude): practical system, depends on the location. Function `azalt`.
- **Of date** (default), **J2000** (`FLG_J2000`), **ICRS** (`FLG_ICRS`): three reference frames for the coordinate zero point.
- **Precession**: slow shift of the zero point (~50"/year). **Nutation**: oscillation of the axis (~9"). **Aberration**: shift due to Earth's motion (~20").
- `cotrans` converts between ecliptic and equatorial. `azalt` and `azalt_rev` convert from/to horizontal coordinates.

### Introduced functions and constants

| Function / Constant | Use |
|---------------------|-----|
| `cotrans(coord, obliquity)` | Ecliptic ↔ Equatorial |
| `cotrans_sp(coord, speed, obliquity)` | As above, with speed |
| `azalt_rev(jd, flag, geopos, az, alt)` | Horizontal → Ecl/Eq |
| `FLG_EQUATORIAL` | Equatorial coordinates |
| `FLG_XYZ` | Cartesian coordinates |
| `FLG_J2000`, `FLG_ICRS` | Reference frames |
| `FLG_NONUT`, `FLG_NOABERR`, `FLG_TRUEPOS` | Disable corrections |
| `ECL2HOR`, `EQU2HOR` | Direction for `azalt` |
| `HOR2ECL`, `HOR2EQU` | Direction for `azalt_rev` |
