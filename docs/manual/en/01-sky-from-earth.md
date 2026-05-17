# Chapter 1 — The sky seen from Earth

## What you will learn

In this chapter you will build a mental model of the sky. By the end you will know what the celestial sphere is, why the Sun traces a path called the ecliptic, how the zodiac divides that path into 12 signs, how your local horizon determines what you can see, and what the Ascendant and the Midheaven are.

These concepts are the foundation for the rest of the manual.

---

## 1.1 The celestial sphere

Imagine standing in an open field, at night, with no artificial lights. The stars seem painted on the inner surface of a huge sphere surrounding you on all sides. This is the **celestial sphere** — a geometric model that astronomers have used for millennia.

The celestial sphere doesn't physically exist: stars are at vastly different distances from each other. But to describe *where* an object appears in the sky, we don't need to know how far away it is — we just need to know in which direction to look. The celestial sphere is exactly this: a map of directions.

The sphere has two special points:

- The **north celestial pole**: the point in the sky where the Earth's axis of rotation points, extended northwards. The North Star (Polaris) is located less than 1° from this point — which is why it seems to stand still while the other stars turn.
- The **south celestial pole**: the opposite point, visible only from the southern hemisphere. It has no bright star marking it.

Halfway between the poles, the sphere is divided by the **celestial equator** — the extension of the Earth's equator projected into the sky. If you stand exactly on the Earth's equator, the celestial equator passes through your zenith (the point directly above your head).

### Why the stars "turn"

If you observe the sky for a few hours, you will see the stars moving slowly from east to west, tracing arcs parallel to the celestial equator. It's not the stars that move — it's the Earth rotating on its axis. The celestial sphere appears to turn in the opposite direction to the Earth's rotation: one complete rotation every 23 hours, 56 minutes, and 4 seconds (the **sidereal day** — we'll meet it in Chapter 2).

### 🌍 Real life

If you point a camera at the North Star with a long exposure, you will see the stars tracing concentric circles around it. Those circles are visual proof of the Earth's rotation.

---

## 1.2 The ecliptic: the path of the Sun

Stars appear fixed on the celestial sphere (they do move, but so slowly that to the naked eye it takes centuries to notice — we'll talk about this in Chapter 8). The Sun, on the other hand, moves: every day it shifts by about 1° relative to the background stars.

If we could see the stars during the day and marked the Sun's position every day for a year, the path would form a great circle on the celestial sphere. This circle is called the **ecliptic**.

The ecliptic does not coincide with the celestial equator: it is inclined to it by about 23.4°. This angle is called the **obliquity of the ecliptic** and is the reason seasons exist. When the Sun is in the part of the ecliptic that lies above the celestial equator, it is summer in the northern hemisphere — the Sun climbs higher in the sky and stays above the horizon longer.

The ecliptic and the celestial equator intersect at two points:

- The **vernal equinox** (or spring equinox, or First Point of Aries ♈): the Sun passes from the southern celestial hemisphere to the northern one. This is the zero point of the zodiac.
- The **autumnal equinox**: the Sun passes from the northern celestial hemisphere to the southern one.

### 💻 Code: obtaining the obliquity of the ecliptic

With LibEphemeris you can get the exact value of the obliquity for any date. The "pseudo-planet" `SE_ECL_NUT` returns information about obliquity and nutation:

```python
import libephemeris as ephem

# Vernal equinox 2024 (March 20, 3:06 UT)
jd = ephem.julday(2024, 3, 20, 3.1)

# SE_ECL_NUT returns obliquity and nutation
nut, flag = ephem.calc_ut(jd, ephem.ECL_NUT, 0)

true_obliquity = nut[0]    # true obliquity (with nutation)
mean_obliquity = nut[1]    # mean obliquity (without nutation)
delta_psi = nut[2]         # nutation in longitude
delta_eps = nut[3]         # nutation in obliquity

print(f"True obliquity:  {true_obliquity:.6f}°")
print(f"Mean obliquity:  {mean_obliquity:.6f}°")
```

```
True obliquity:  23.438703°
Mean obliquity:  23.436129°
```

The difference between true and mean obliquity is the **nutation in obliquity** (`delta_eps`): a slight wobble of the Earth's axis which we will discuss in Chapter 3.

### 🌍 Real life

In summer, the Sun climbs higher in the sky because it is located in the part of the ecliptic above the celestial equator. At the summer solstice, in Rome (lat. 41.9° N), the Sun reaches a maximum altitude of about 71.5° — almost at the zenith. At the winter solstice, barely 24.7°.

---

## 1.3 The zodiac: 12 sectors of 30°

The ecliptic is a 360° circle. For convenience, since ancient times it has been divided into **12 equal sectors of 30° each**, called **zodiac signs**:

| Sign | Range | Sign | Range |
|------|-------|------|-------|
| ♈ Aries | 0° – 30° | ♎ Libra | 180° – 210° |
| ♉ Taurus | 30° – 60° | ♏ Scorpio | 210° – 240° |
| ♊ Gemini | 60° – 90° | ♐ Sagittarius | 240° – 270° |
| ♋ Cancer | 90° – 120° | ♑ Capricorn | 270° – 300° |
| ♌ Leo | 120° – 150° | ♒ Aquarius | 300° – 330° |
| ♍ Virgo | 150° – 180° | ♓ Pisces | 330° – 360° |

The zero point — 0° Aries — coincides with the vernal equinox (First Point of Aries). When an astronomer says "the Sun is at 105° of ecliptic longitude", an astrologer says "the Sun is at 15° Cancer" (because 105 = 90 + 15, and Cancer starts at 90°).

### Signs and constellations: not the same thing

This is a fundamental point that generates perpetual confusion. Zodiac **signs** are geometric sectors of 30° each, defined starting from the vernal equinox. Zodiac **constellations** are groups of stars with different sizes (Virgo covers over 40° of the ecliptic, Scorpio less than 20°) and conventional boundaries established by the International Astronomical Union in 1930.

About 2000 years ago, signs and constellations roughly coincided. Not anymore: the vernal equinox has shifted by about 24° towards the stars of Pisces. This shift — about 50" of arc per year — is called the **precession of the equinoxes**, and it's due to the fact that the Earth's axis slowly wobbles like a spinning top (one complete wobble takes about 26000 years). We will talk about this in detail in Chapter 11, where we will look at the sidereal zodiac used in Vedic astrology.

### 🌍 Real life

When someone says "I'm an Aries", they mean that on the date of their birth the Sun was between 0° and 30° of tropical ecliptic longitude. This has nothing to do with the constellation of Aries: due to precession, that stretch of sky is today aligned with the stars of Pisces.

---

## 1.4 The local horizon

So far we have described the sky "in general" — as it would appear from any point on Earth. But the sky **you** see depends on **where** you are.

Three points define your local sky:

- **Zenith**: the point directly above your head, 90° from the horizon.
- **Nadir**: the point directly below your feet, opposite the zenith.
- **Horizon**: the 360° circle around you where sky and earth meet.

To locate an object in the sky from your viewing position, you need two coordinates:

- **Altitude** (or elevation): the angle above the horizon. 0° = on the horizon. 90° = at the zenith. Negative values = below the horizon.
- **Azimuth**: the direction along the horizon. In astronomy, the most widely used convention starts from the South: 0° = South, 90° = West, 180° = North, 270° = East.

> **Warning**: the azimuth convention varies. In navigation and daily life, 0° = North. In astronomy (and in LibEphemeris), 0° = South. If you pass data between different systems, check the convention.

### 💻 Code: altitude and azimuth of a celestial body

The `azalt` function converts ecliptic or equatorial coordinates into horizontal coordinates. It requires the observer's position (longitude, latitude, altitude):

```python
import libephemeris as ephem

# September 15, 2024, 21:00 UT
jd = ephem.julday(2024, 9, 15, 21.0)

# Position of Jupiter
pos, flag = ephem.calc_ut(jd, ephem.JUPITER, 0)

# Observer's position: Rome
# (longitude East, latitude North, altitude in meters)
geopos = (12.4964, 41.9028, 50.0)

# Convert from ecliptic to horizontal coordinates
# SE_ECL2HOR = from ecliptic to horizontal
# atpress = 1013.25 mbar (standard pressure)
# attemp = 15.0 °C (standard temperature)
hor = ephem.azalt(jd, ephem.ECL2HOR, geopos, 1013.25, 15.0,
                  (pos[0], pos[1], pos[2]))

azimuth = hor[0]           # from South, towards West
true_altitude = hor[1]     # without refraction
apparent_altitude = hor[2] # with atmospheric refraction

# Convert the azimuth from the astronomical convention (S=0)
# to the navigational convention (N=0) for readability
az_nav = (azimuth + 180.0) % 360.0

# Approximate direction
directions = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
dir_idx = int((az_nav + 22.5) % 360 / 45)
dir_name = directions[dir_idx]

if apparent_altitude > 0:
    print(f"Jupiter is visible! Altitude: {apparent_altitude:.1f}°, "
          f"direction: {dir_name} (azimuth {az_nav:.1f}° from N)")
else:
    print(f"Jupiter is below the horizon ({apparent_altitude:.1f}°)")
```

```
Jupiter is below the horizon (-3.2°)
```

Pressure and temperature are needed to calculate **atmospheric refraction**: the atmosphere bends light, making objects near the horizon appear a bit higher than they really are (about 34' at the horizon). If you pass `atpress = 0`, refraction is ignored.

---

## 1.5 The four angles: Ascendant, MC, Descendant, IC

The local horizon and meridian intersect the ecliptic at four points that are fundamental in astrology. These four points change continuously because the celestial sphere rotates — and it rotates fast: the Ascendant moves by about 1° every 4 minutes.

### The Ascendant (ASC)

The degree of the ecliptic that is currently rising in the east. It is the intersection point between the ecliptic and the eastern horizon. In a birth chart, the Ascendant defines the cusp of the first house.

The Ascendant depends on the **exact time** and **location**. Two people born on the same day but 10 minutes apart can have different Ascendants. This is why the time of birth is so important in astrology.

### The Midheaven (MC)

The degree of the ecliptic that is currently culminating — that is, crossing the upper meridian (the semicircle running from the north celestial pole to the zenith to the south celestial pole). The MC is the highest point that degree of the ecliptic will reach on its path across the sky. In a birth chart, the MC defines the cusp of the tenth house.

### Descendant (DSC) and Imum Coeli (IC)

The **Descendant** is opposite the Ascendant: the degree of the ecliptic that is setting in the west. The **Imum Coeli** (IC, "bottom of the sky") is opposite the MC: the degree of the ecliptic at the lower culmination, beneath the observer's feet. In a birth chart, DSC and IC define the cusp of the seventh and fourth houses, respectively.

> **Note**: MC and ASC are **not** always 90° apart. The distance depends on the observer's latitude and the time of day. Only at the equator, and only at certain times, are MC and ASC exactly 90° apart.

### 💻 Code: calculating the four angles

```python
import libephemeris as ephem

# April 8, 2024, 12:00 UT — solar eclipse
jd = ephem.julday(2024, 4, 8, 12.0)

# Rome: lat 41.9, lon 12.5
lat, lon = 41.9028, 12.4964

# ord('P') = Placidus, the most popular house system
cusps, ascmc = ephem.houses(jd, lat, lon, ord('P'))

asc     = ascmc[0]   # Ascendant
mc      = ascmc[1]   # Midheaven
armc    = ascmc[2]   # Right Ascension of the Midheaven (in degrees)
vertex  = ascmc[3]   # Vertex

signs = [
    "Aries", "Taurus", "Gemini", "Cancer",
    "Leo", "Virgo", "Libra", "Scorpio",
    "Sagittarius", "Capricorn", "Aquarius", "Pisces",
]

def format_position(degrees):
    """Converts decimal degrees to zodiacal format."""
    deg, min, sec, secfr, sign = ephem.split_deg(
        degrees, ephem.SPLIT_DEG_ZODIACAL | ephem.SPLIT_DEG_ROUND_SEC
    )
    return f"{deg}° {min}' {sec}\" {signs[sign]}"

print(f"Ascendant:   {format_position(asc)}")
print(f"Midheaven:   {format_position(mc)}")
print(f"Descendant:  {format_position((asc + 180) % 360)}")
print(f"Imum Coeli:  {format_position((mc + 180) % 360)}")
print(f"Vertex:      {format_position(vertex)}")
```

```
Ascendant:   13° 4' 44" Leo
Midheaven:   1° 54' 15" Taurus
Descendant:  13° 4' 44" Aquarius
Imum Coeli:  1° 54' 15" Scorpio
Vertex:      0° 46' 36" Capricorn
```

The `houses` function returns two tuples:

- `cusps`: the 12 house cusps (we will look at them in Chapter 7)
- `ascmc`: the principal angles — ASC, MC, ARMC, Vertex, plus four other specialized angles (Equatorial Ascendant, Koch co-Ascendant, Munkasey co-Ascendant, Polar Ascendant)

The `ord('P')` parameter indicates the Placidus house system, the most widely used in the West. There are over 20 different systems — we will explore them in Chapter 7.

### 🌍 Real life

When an astrologer asks "what is your Ascendant?", they are asking which degree of the ecliptic was rising at the exact moment of your birth, seen from the place where you were born. The Ascendant changes sign roughly every 2 hours (but not regularly — some signs rise faster than others, depending on latitude). This is why a precise date, time, and location are needed.

At very high latitudes (above the polar circle), at certain times of the year the ecliptic does not intersect the horizon at all, and the Ascendant cannot be calculated with some house systems. LibEphemeris handles this situation with automatic fallbacks — we'll talk about this in Chapter 7.

---

## Summary

- The **celestial sphere** is a geometric model: an imaginary sphere onto which we project the positions of celestial objects.
- The **ecliptic** is the Sun's annual path on the celestial sphere, inclined by ~23.4° relative to the celestial equator.
- The **zodiac** divides the ecliptic into 12 signs of 30° each, starting from the vernal equinox. Signs are geometric sectors, not constellations.
- The **local horizon** depends on your location. Altitude and azimuth describe where an object appears in your sky.
- The **Ascendant** is the degree of the ecliptic rising in the east; the **Midheaven** is the degree that culminates. They change rapidly and depend on time and place.

### Introduced functions and constants

| Function / Constant | Use |
|---------------------|-----|
| `calc_ut(jd, ephem.ECL_NUT, 0)` | Obliquity of the ecliptic and nutation |
| `azalt(jd, calc_flag, geopos, atpress, attemp, xin)` | Ecliptic/equatorial → horizontal coordinates |
| `houses(jd, lat, lon, ord('P'))` | House cusps and angles (ASC, MC, ...) |
| `split_deg(degrees, flags)` | Zodiacal formatting |
| `SE_ECL2HOR`, `SE_EQU2HOR` | Flags for `azalt`: ecliptic or equatorial input |
