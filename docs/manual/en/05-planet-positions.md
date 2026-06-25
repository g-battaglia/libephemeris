# Chapter 5 — Calculating Planetary Positions

## What You Will Learn

In this chapter, you will master `calc_ut` — the most used function in the library. You will learn what it returns, all available flags, how velocity and retrograde motion work, and how to calculate planetary phenomena such as magnitude, phase, and elongation.

---

## 5.1 The Main Function: `calc_ut`

The signature is simple:

```python
pos, flag = ephem.calc_ut(jd_ut, body, iflag)
```

- `jd_ut`: Julian Day in UT ("civil" time)
- `body`: celestial body identifier (`SUN`, `MOON`, etc.)
- `iflag`: calculation flags (combined with `|`)

The result consists of two parts: a `pos` tuple with 6 numbers and an integer `flag` confirming the calculation options used. The 6 numbers describe **where** the celestial body is located and **how it is moving**.

### The First Three: Where It Is

To understand these values, imagine the zodiac as a long circular road of 360 km. Every planet has a position along this road, a lateral distance from the road, and a distance from you.

- **`pos[0]` — Ecliptic Longitude** is the planet's position along this "zodiacal road". It ranges from 0° to 360° and corresponds to the 12 signs: the first 30° are Aries, from 30° to 60° is Taurus, and so on. If `pos[0]` is 105°, it means the planet is located at 15° of Cancer (because Cancer starts at 90°, and 105 - 90 = 15). It is the most important data point in astrology: when someone says "my Sun is in Leo", they mean that the Sun's ecliptic longitude at the time of their birth was between 120° and 150°.

- **`pos[1]` — Ecliptic Latitude** is the planet's distance from the ecliptic plane, measured in degrees. Think of the ecliptic as a floor: the latitude tells you how far the planet is above (+) or below (-) that floor. The Sun is always at ~0° (it defines the floor itself). The Moon can rise up to ±5.1°. Most planets stay within a few degrees of the ecliptic, but Pluto can reach up to ~17°.

- **`pos[2]` — Distance** is how far the celestial body is from Earth, measured in **astronomical units** (1 AU = ~150 million km, the average Earth-Sun distance). To give you an idea: the Moon is very close (~0.0026 AU, about 384,000 km), Mars varies greatly between 0.37 AU when it's on our side of the orbit and 2.68 AU when it's on the other side of the Sun, and Neptune is at about 30 AU.

### The Last Three: How It Is Moving

- **`pos[3]` — Speed in Longitude** tells you how many degrees the planet advances (or retreats) along the zodiac each day. The Sun moves at ~1°/day — almost imperceptible to the naked eye. The Moon races at ~13°/day, crossing an entire sign in just over two days. Saturn drags along at ~0.03°/day. If this value is **negative**, the planet is in **retrograde motion**: it appears to move backwards along the zodiac (we will talk about this in section 5.3).

- **`pos[4]` — Speed in Latitude** tells you how much the planet is approaching or moving away from the ecliptic each day. It is usually a small value — latitude changes slowly.

- **`pos[5]` — Speed in Distance** tells you whether the planet is approaching (negative value) or moving away from (positive) the Earth, in AU/day.

### How Values Change with Flags

The values described above are what you get with the default flag (`0`). But if you pass certain flags, the **meaning** of the 6 numbers changes completely — the structure of the tuple remains the same, but the numbers inside represent different things:

- With **`FLG_EQUATORIAL`**, `pos[0]` is no longer the ecliptic longitude but the **Right Ascension** (the position along the celestial equator, from 0° to 360°), and `pos[1]` becomes the **Declination** (the distance from the celestial equator, from -90° to +90°). This is the system used by telescopes. Distance and speeds keep the same meaning, but referred to the equatorial system.

- With **`FLG_XYZ`**, all 6 values become **Cartesian coordinates**: `pos[0]`, `pos[1]`, `pos[2]` are the X, Y, Z positions in astronomical units, and `pos[3]`, `pos[4]`, `pos[5]` are the corresponding speeds in AU/day. There are no more degrees — only distances and speeds in three-dimensional space.

- With **`FLG_RADIANS`**, angular values (longitude, latitude, speed) are expressed in **radians** instead of degrees (2π radians = 360°). Useful if you need to perform trigonometric calculations without converting.

If you don't use any of these flags (or pass `0`), you always get ecliptic coordinates in degrees — the system of astrology.

### Main Celestial Bodies

| Constant | Value | Body |
|----------|-------|------|
| `SUN` | 0 | Sun |
| `MOON` | 1 | Moon |
| `MERCURY` | 2 | Mercury |
| `VENUS` | 3 | Venus |
| `MARS` | 4 | Mars |
| `JUPITER` | 5 | Jupiter |
| `SATURN` | 6 | Saturn |
| `URANUS` | 7 | Uranus |
| `NEPTUNE` | 8 | Neptune |
| `PLUTO` | 9 | Pluto |
| `MEAN_NODE` | 10 | Mean Lunar Node |
| `TRUE_NODE` | 11 | True Lunar Node |
| `MEAN_APOG` | 12 | Mean Lunar Apogee (Lilith) |
| `OSCU_APOG` | 13 | Osculating Lunar Apogee |
| `CHIRON` | 15 | Chiron |

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# All planets at once
bodies = [
    (ephem.SUN, "Sun"), (ephem.MOON, "Moon"),
    (ephem.MERCURY, "Mercury"), (ephem.VENUS, "Venus"),
    (ephem.MARS, "Mars"), (ephem.JUPITER, "Jupiter"),
    (ephem.SATURN, "Saturn"), (ephem.URANUS, "Uranus"),
    (ephem.NEPTUNE, "Neptune"), (ephem.PLUTO, "Pluto"),
]

signs = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

for body_id, name in bodies:
    pos, _ = ephem.calc_ut(jd, body_id, ephem.FLG_SPEED)
    lon = pos[0]
    sign = signs[int(lon / 30)]
    degrees = lon % 30
    print(f"{name:10s} {degrees:5.1f}° {sign}  (vel: {pos[3]:+.4f}°/d)")
```

```text
Sun         19.1° Ari  (vel: +0.9831°/d)
Moon        15.4° Ari  (vel: +15.0292°/d)
Mercury     25.0° Ari  (vel: -0.6229°/d)
Venus        4.1° Ari  (vel: +1.2353°/d)
Mars        12.8° Psc  (vel: +0.7775°/d)
Jupiter     19.0° Tau  (vel: +0.2211°/d)
Saturn      14.4° Psc  (vel: +0.1078°/d)
Uranus      21.2° Tau  (vel: +0.0512°/d)
Neptune     28.2° Psc  (vel: +0.0359°/d)
Pluto        2.0° Aqr  (vel: +0.0114°/d)
```

---

## 5.2 Calculation Flags

The third argument of `calc_ut` is an integer that controls *how* the calculation is performed. Passing `0` means "use all default settings": ecliptic, geocentric coordinates, with aberration and nutation.

To change behavior, use the `FLG_*` constants and combine them with the `|` (bitwise OR) operator. For example, `ephem.FLG_EQUATORIAL | ephem.FLG_SPEED` requests equatorial coordinates with speed. You can combine as many as you want.

Here are the most important ones, grouped by category.

**Coordinates** — they change the meaning of the returned values:

- No flag (default): **ecliptic** coordinates — longitude, latitude, distance
- `FLG_EQUATORIAL`: **equatorial** coordinates — Right Ascension, Declination, distance
- `FLG_XYZ`: **Cartesian** coordinates — X, Y, Z in astronomical units
- `FLG_RADIANS`: angles in **radians** instead of degrees

**Center of observation** — where you are "looking" from:

- No flag (default): **geocentric** — from the center of the Earth
- `FLG_HELCTR`: **heliocentric** — from the center of the Sun
- `FLG_BARYCTR`: **barycentric** — from the center of mass of the Solar System
- `FLG_TOPOCTR`: **topocentric** — from your position on the Earth's surface (requires `set_topo()` first)

**Reference system** — which "zero point" to use:

- No flag (default): **equinox of date** — the system of astrology
- `FLG_J2000`: referred to J2000.0 — the system of astronomical catalogs
- `FLG_ICRS`: ICRS system — almost identical to J2000
- `FLG_SIDEREAL`: **sidereal** coordinates — requires `set_sid_mode()` first
- `FLG_NONUT`: without nutation — returns "mean" coordinates instead of "true" ones

**Corrections** — enables or disables physical effects:

- `FLG_SPEED`: also calculates **speeds** (pos[3], pos[4], pos[5]). Almost always useful.
- `FLG_TRUEPOS`: **true geometric** position, without light-time correction
- `FLG_NOABERR`: without **aberration** (the displacement due to Earth's motion)
- `FLG_NOGDEFL`: without **gravitational deflection** of light

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Common combination: equatorial with speed
pos, _ = ephem.calc_ut(jd, ephem.MARS,
                       ephem.FLG_EQUATORIAL | ephem.FLG_SPEED)

print(f"RA: {pos[0]:.4f}°, Dec: {pos[1]:+.4f}°")
print(f"Vel RA: {pos[3]:+.4f}°/d, Vel Dec: {pos[4]:+.4f}°/d")
```

```text
RA: 344.6682°, Dec: -7.8866°
Vel RA: +0.7256°/d, Vel Dec: +0.2960°/d
```

---

## 5.3 Speed and Retrograde Motion

The speed in longitude (`pos[3]`) indicates how much the planet moves each day:

- **Positive** = direct motion (the planet advances along the zodiac)
- **Negative** = retrograde motion (the planet appears to go backwards)
- **Zero** = station (the planet "stops" before reversing direction)

Retrograde motion is an optical effect: when the Earth overtakes an outer planet (or is overtaken by an inner one), that planet appears to move backward relative to the stars. It is like overtaking a car on the highway — for a moment it seems like the other car is going backward.

```python
import libephemeris as ephem
from libephemeris.crossing import is_retrograde

jd = ephem.julday(2024, 4, 8, 12.0)

# Check if Mercury is retrograde
retro = is_retrograde(ephem.MERCURY, jd)
print(f"Mercury retrograde: {retro}")

# Find Mercury's next station
jd_station, station_type = ephem.find_station_ut(ephem.MERCURY, jd)

year, month, day, hours = ephem.revjul(jd_station)
# station_type = "SR" (retrograde station) or "SD" (direct station)
print(f"Next station: {day}/{month}/{year} ({station_type})")
```

```text
Mercury retrograde: True
Next station: 25/4/2024 (SD)
```

### 🌍 Real Life

"Mercury retrograde" is one of the most popular astrological concepts. It happens about 3 times a year, for roughly 3 weeks each time. In reality, all planets (except the Sun and Moon) have retrograde periods — outer planets are retrograde for months.

---

## 5.4 Planetary Phenomena: Magnitude, Phase, Elongation

The `pheno_ut` function returns information on the physical appearance of a planet:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 21.0)

attr = ephem.pheno_ut(jd, ephem.JUPITER, 0)

phase_angle = attr[0]    # Sun-Planet-Earth angle (degrees)
phase = attr[1]          # illuminated fraction of the disk (0.0–1.0)
elongation = attr[2]     # angular distance from the Sun (degrees)
diameter = attr[3]       # apparent diameter (arcseconds)
magnitude = attr[4]      # apparent visual magnitude

print(f"Jupiter:")
print(f"  Elongation: {elongation:.1f}°")
print(f"  Phase: {phase:.2f} ({phase*100:.0f}% illuminated)")
print(f"  Magnitude: {magnitude:.1f}")
print(f"  Apparent diameter: {diameter:.1f}\"")
```

```text
Jupiter:
  Elongation: 29.6°
  Phase: 1.00 (100% illuminated)
  Magnitude: -2.0
  Apparent diameter: 0.0"
```

**Elongation** is the angular distance from the Sun — if it is small (< 10°), the planet is lost in the solar glare and not visible. **Magnitude** indicates brightness: lower numbers = brighter (Venus reaches -4.6, Jupiter -2.9).

---

## 5.5 Orbital Elements

The `get_orbital_elements` function returns the orbital parameters of a celestial body — the 6 numbers that describe an ellipse in space:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)
jd_tt = jd + ephem.deltat(jd)  # JD in TT is required

elem = ephem.get_orbital_elements(jd_tt, ephem.MARS, 0)

print(f"Mars - orbital elements:")
print(f"  Semi-major axis:     {elem[0]:.6f} AU")
print(f"  Eccentricity:        {elem[1]:.6f}")
print(f"  Inclination:         {elem[2]:.4f}°")
print(f"  Sidereal period:     {elem[10]:.2f} years")
print(f"  Synodic period:      {elem[13]:.1f} days")
print(f"  Perihelion dist.:    {elem[15]:.4f} AU")
print(f"  Aphelion dist.:      {elem[16]:.4f} AU")
```

```text
Mars - orbital elements:
  Semi-major axis:     1.523627 AU
  Eccentricity:        0.093278
  Inclination:         1.8479°
  Sidereal period:     1.88 years
  Synodic period:      779.9 days
  Perihelion dist.:    1.3815 AU
  Aphelion dist.:      1.6657 AU
```

To get the maximum and minimum distances with a single call:

```python
d_max, d_min, d_now = ephem.orbit_max_min_true_distance(
    jd, ephem.MARS, 0
)
print(f"Mars: min {d_min:.4f} AU, max {d_max:.4f} AU, now {d_now:.4f} AU")
```

```text
Mars: min 0.3648 AU, max 2.6825 AU, now 2.0618 AU
```

---

## 5.6 Nodes and Apsides

The **nodes** are the points where a planet's orbit intersects the ecliptic. The **apsides** are the points of maximum and minimum distance from the Sun (aphelion and perihelion).

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# nod_aps_ut returns 4 tuples of 6 values each:
# ascending node, descending node, perihelion, aphelion
nasc, ndsc, peri, aphe = ephem.nod_aps_ut(
    jd, ephem.MARS, 0, 0  # method=0: mean
)

print(f"Mars - ascending node: {nasc[0]:.4f}°")
print(f"Mars - descending node: {ndsc[0]:.4f}°")
print(f"Mars - perihelion: {peri[0]:.4f}° (dist {peri[2]:.4f} AU)")
print(f"Mars - aphelion:   {aphe[0]:.4f}° (dist {aphe[2]:.4f} AU)")
```

```text
Mars - ascending node: 37.4197°
Mars - descending node: 266.1986°
Mars - perihelion: 354.2716° (dist 2.2240 AU)
Mars - aphelion:   120.3522° (dist 1.1508 AU)
```

---

## Summary

- `calc_ut(jd, body, flag)` is the core function of the library. Given an instant and a celestial body, it returns 6 numbers: the first three tell where it is (longitude, latitude, distance), the last three how it is moving.
- The **flags** control the type of coordinates, the observation center, the reference system, and physical corrections. They are combined with `|`. Passing `0` yields geocentric ecliptic coordinates — the default for astrology.
- The **speed in longitude** (`pos[3]`) is positive when the planet advances along the zodiac (direct motion) and negative when it seems to go backwards (retrograde motion).
- `pheno_ut` provides information on the physical appearance of the planet: magnitude (brightness), phase (illuminated percentage), elongation (angular distance from the Sun), and apparent diameter.
- `get_orbital_elements` returns the parameters of the Keplerian orbit: semi-major axis, eccentricity, inclination, periods, and distances.
- `nod_aps_ut` returns the nodes (where the orbit crosses the ecliptic) and apsides (points of minimum and maximum distance from the Sun).

### Introduced Functions

- `pheno_ut(jd, body, flag)` — phenomena: phase, elongation, magnitude
- `is_retrograde(body, jd)` — is the planet currently retrograde? (import from `libephemeris.crossing`)
- `find_station_ut(body, jd)` — finds the next station (retrograde or direct)
- `get_orbital_elements(jd_tt, body, flag)` — Keplerian orbital elements
- `orbit_max_min_true_distance(jd, body, flag)` — min, max, and current distances from Earth
- `nod_aps_ut(jd, body, flag, method)` — orbit nodes and apsides
