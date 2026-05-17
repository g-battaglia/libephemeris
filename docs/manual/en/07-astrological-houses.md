# Chapter 7 — Astrological Houses

## What You Will Learn

In this chapter, you will discover what astrological houses are, why there are over 20 different systems to calculate them, how the Ascendant and Midheaven work, how to determine which house a planet falls into, and what happens at extreme latitudes where some systems fail to work.

---

## 7.1 What Are Houses

Imagine taking the celestial sphere and dividing it into 12 slices, like a pie. Each slice is a **house**. While the zodiac signs divide the ecliptic based on the Sun's annual path (and therefore depend only on the date), houses divide the sky based on the Earth's daily rotation — and therefore depend on the **date**, **time**, and **location**.

Each house "rules" a specific area of life according to astrological tradition:

- **1st house** — the self, physical appearance, temperament
- **2nd house** — resources, money, personal values
- **3rd house** — communication, siblings, short trips
- **4th house** — family, roots, physical home
- **5th house** — creativity, children, pleasure
- **6th house** — daily work, health, service
- **7th house** — relationships, partners, contracts
- **8th house** — transformations, inheritance, sexuality
- **9th house** — long journeys, philosophy, spirituality
- **10th house** — career, reputation, ambition
- **11th house** — friendships, groups, future projects
- **12th house** — the unconscious, isolation, transcendence

The fundamental difference between signs and houses is: two people born on the same day have their planets in the **same signs**, but if they were born at different times or in different places, their planets fall into **different houses**. The time of birth determines the houses.

The starting point of each house is called a **cusp**. The cusp of the 1st house is the Ascendant.

---

## 7.2 The Ascendant and Midheaven

We have already encountered these concepts in Chapter 1, but here we will delve into their role in the house system.

The **Ascendant** (ASC) is the degree of the ecliptic that is rising on the eastern horizon at a given time and place. It is the cusp of the 1st house in all house systems. It changes by about 1° every 4 minutes — which is why the exact time of birth is so important.

The **Midheaven** (MC) is the degree of the ecliptic that culminates at the meridian — the highest point that degree reaches in the sky. It is the cusp of the 10th house in most systems (but not all: the whole sign system and the equal house system do not use the MC as the 10th house cusp).

The `houses` function returns two sets of data:

- The **12 house cusps** (an ecliptic longitude for each house)
- The **angles**: Ascendant, Midheaven, ARMC, Vertex, and other special points

```python
import libephemeris as ephem

# April 8 2024, 14:30 UT — Rome
jd = ephem.julday(2024, 4, 8, 14.5)
lat, lon = 41.9028, 12.4964

# Calculate houses with Placidus system
cusps, ascmc = ephem.houses(jd, lat, lon, ord('P'))

# The 12 cusps
signs = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

for i in range(12):
    lon_c = cusps[i]
    sign = signs[int(lon_c / 30)]
    deg = lon_c % 30
    print(f"House {i+1:2d}: {deg:5.1f}° {sign}")

# Main angles
print(f"\nAscendant:  {ascmc[0]:.4f}°")
print(f"Midheaven:  {ascmc[1]:.4f}°")
print(f"ARMC:       {ascmc[2]:.4f}°")
print(f"Vertex:     {ascmc[3]:.4f}°")
```

```
House  1:  12.2° Vir
House  2:   6.3° Lib
House  3:   5.5° Sco
House  4:   9.0° Sgr
House  5:  13.5° Cap
House  6:  15.0° Aqr
House  7:  12.2° Psc
House  8:   6.3° Ari
House  9:   5.5° Tau
House 10:   9.0° Gem
House 11:  13.5° Cnc
House 12:  15.0° Leo

Ascendant:  162.2479°
Midheaven:  69.0379°
ARMC:       67.3366°
Vertex:     316.3110°
```

The **ARMC** (Right Ascension of the Midheaven) is the same concept as local sidereal time, expressed in degrees instead of hours (ARMC in degrees = sidereal time in hours × 15). It is the starting point for calculating all house systems.

The **Vertex** is the point where the prime vertical (the circle passing through East, Zenith, and West) intersects the ecliptic on the western side. In some astrological schools, it has a meaning related to fateful encounters and karmic relationships.

---

## 7.3 House Systems: Which One to Choose?

Why are there over 20 house systems? Because there is no single, unequivocal way to divide a three-dimensional sphere into 12 sectors. Each system uses a different criterion — dividing time, dividing space, dividing the equator — and each produces slightly different results.

Here are the most commonly used systems, with the letter to pass to `houses()`:

**Placidus** (`P`) has been the most widespread system in the West since the 17th century. It divides the **time** that each degree of the ecliptic takes to travel from the horizon to the meridian. It is based on the concept of "semi-arc": the time between the rising and culmination of a point. Its supporters consider it the most "natural" because it reflects the real motion of the sky. **Problem**: it does not work beyond the polar circle (~66.5° latitude).

**Koch** (`K`) is similar to Placidus but uses the semi-arc of the Midheaven instead of the point itself. Popular in Germany. Same problem at extreme latitudes.

**Regiomontanus** (`R`) divides the **celestial equator** into 12 equal parts of 30° each, then projects these divisions onto the ecliptic. It was the dominant system in the Middle Ages and the Renaissance.

**Campanus** (`C`) divides the **prime vertical** (the East-Zenith-West circle) into 12 equal parts. It is one of the oldest systems with a clear geometric basis.

**Equal from the Ascendant** (`E`) is the simplest: each house has exactly 30°, starting from the Ascendant. The 1st house goes from the Ascendant to Ascendant + 30°, the 2nd from +30° to +60°, and so on. The MC in this system is **not** necessarily the cusp of the 10th house — it can fall in any house between the 8th and the 11th. Widely used in ancient astrology and the Hellenistic tradition.

**Whole Sign** (`W`) is even simpler: each zodiac sign is a house. If the Ascendant is in Leo, the entire 1st house coincides with Leo (0°–30° Leo), the 2nd with Virgo, and so on. It is the oldest system, used in Greek and Indian astrology. It is currently experiencing a strong revival in modern astrology.

**Porphyry** (`O`) proportionally divides the four quadrants (ASC-IC, IC-DSC, DSC-MC, MC-ASC). It is a good compromise and works at **all latitudes** — which is why it is often used as a fallback when Placidus or Koch fail.

**Polich/Page** (`T`, Topocentric) is a topographical projection. Very popular in South America.

**Morinus** (`M`) divides the celestial equator from the MC. It works at all latitudes.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 14.5)
lat, lon = 41.9028, 12.4964

# Comparison: Placidus vs Equal vs Whole Sign
for letter, name in [('P', 'Placidus'), ('E', 'Equal'), ('W', 'Whole Sign')]:
    cusps, ascmc = ephem.houses(jd, lat, lon, ord(letter))
    print(f"\n{name} (1st cusp = ASC = {ascmc[0]:.1f}°):")
    for i in range(12):
        print(f"  House {i+1:2d}: {cusps[i]:.1f}°")

# Get system name from letter
print(ephem.house_name(ord('P')))  # "Placidus"
print(ephem.house_name(ord('W')))  # "Whole Sign"
```

```
Placidus (1st cusp = ASC = 162.2°):
  House  1: 162.2°
  House  2: 186.3°
  House  3: 215.5°
  House  4: 249.0°
  House  5: 283.5°
  House  6: 315.0°
  House  7: 342.2°
  House  8: 6.3°
  House  9: 35.5°
  House 10: 69.0°
  House 11: 103.5°
  House 12: 135.0°

Equal (1st cusp = ASC = 162.2°):
  House  1: 162.2°
  House  2: 192.2°
  House  3: 222.2°
  House  4: 252.2°
  House  5: 282.2°
  House  6: 312.2°
  House  7: 342.2°
  House  8: 12.2°
  House  9: 42.2°
  House 10: 72.2°
  House 11: 102.2°
  House 12: 132.2°

Whole Sign (1st cusp = ASC = 162.2°):
  House  1: 150.0°
  House  2: 180.0°
  House  3: 210.0°
  House  4: 240.0°
  House  5: 270.0°
  House  6: 300.0°
  House  7: 330.0°
  House  8: 0.0°
  House  9: 30.0°
  House 10: 60.0°
  House 11: 90.0°
  House 12: 120.0°

Placidus
Whole Sign
```

---

## 7.4 Which House Does a Planet Fall Into?

Knowing the house cusps is only half the job. The most frequent question is: **which house is my Sun in? And my Moon?** To answer this, you need the `house_pos` function.

### The Concept

Determining the house of a planet is not as trivial as it seems. It is not enough to simply compare the planet's longitude with the cusps. Why? Because many house systems (Placidus, Koch, Regiomontanus) work in **three dimensions**: they do not merely divide the ecliptic, but three-dimensional space. A planet with a non-zero ecliptic latitude could fall into a different house compared to a point with the same longitude but zero latitude.

The `house_pos` function performs the correct calculation for the chosen system. It returns a decimal number where:

- The **integer part** is the house number (1–12)
- The **decimal part** indicates how far "along" the planet is in the house (0.0 = exactly on the cusp, 0.99 = almost at the next cusp)

For example, a value of `7.50` means "exactly halfway through the 7th house".

### How to Use `house_pos`

The function requires:

- **ARMC**: the Right Ascension of the Midheaven (obtained from `houses()` as `ascmc[2]`)
- **Geographic latitude** of the observer
- **Obliquity of the ecliptic** (can be obtained with `calc_ut(jd, ECL_NUT)`)
- **House system** (the letter, such as `ord('P')`)
- **Ecliptic longitude and latitude** of the planet

```python
import libephemeris as ephem

# April 8 2024, 14:30 UT — Rome
jd = ephem.julday(2024, 4, 8, 14.5)
lat, lon = 41.9028, 12.4964

# 1. Calculate houses to get ARMC and cusps
cusps, ascmc = ephem.houses(jd, lat, lon, ord('P'))
armc = ascmc[2]  # ARMC in degrees

# 2. Get ecliptic obliquity
nut, _ = ephem.calc_ut(jd, ephem.ECL_NUT, 0)
obliquity = nut[0]  # true obliquity

# 3. Calculate Sun position
sun, _ = ephem.calc_ut(jd, ephem.SUN, ephem.FLG_SPEED)
sun_lon = sun[0]   # ecliptic longitude
sun_lat = sun[1]   # ecliptic latitude (near zero for the Sun)

# 4. Determine which house it falls in
pos = ephem.house_pos(armc, lat, obliquity, ord('P'), sun_lon, sun_lat)
house = int(pos)
position = pos - house  # how far along in the house (0.0 - 0.99)

print(f"The Sun is in house {house}")
print(f"Position within house: {position:.2%}")
```

```
The Sun is in house 8
Position within house: 46.37%
```

### All Planets in Their Houses

Here is a complete example showing each planet with its sign, degrees, and house:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 14.5)
lat, lon = 41.9028, 12.4964

# Calculate houses and obliquity
cusps, ascmc = ephem.houses(jd, lat, lon, ord('P'))
armc = ascmc[2]
nut, _ = ephem.calc_ut(jd, ephem.ECL_NUT, 0)
obliquity = nut[0]

signs = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

planets = [
    (ephem.SUN,      "Sun     "),
    (ephem.MOON,     "Moon    "),
    (ephem.MERCURY,  "Mercury "),
    (ephem.VENUS,    "Venus   "),
    (ephem.MARS,     "Mars    "),
    (ephem.JUPITER,  "Jupiter "),
    (ephem.SATURN,   "Saturn  "),
]

for body_id, name in planets:
    pos, _ = ephem.calc_ut(jd, body_id, ephem.FLG_SPEED)
    p_lon, p_lat = pos[0], pos[1]

    # Sign and degrees
    sign = signs[int(p_lon / 30)]
    deg = p_lon % 30

    # House
    hp = ephem.house_pos(armc, lat, obliquity, ord('P'), p_lon, p_lat)
    house = int(hp)

    print(f"{name}  {deg:5.1f}° {sign}  →  house {house:2d}")
```

```
Sun        19.2° Ari  →  house  8
Moon       17.0° Ari  →  house  8
Mercury    24.9° Ari  →  house  8
Venus       4.2° Ari  →  house  7
Mars       12.9° Psc  →  house  7
Jupiter    19.0° Tau  →  house  9
Saturn     14.4° Psc  →  house  7
```

### Does Ecliptic Latitude Matter?

For most planets, the ecliptic latitude is small (the Sun has an ecliptic latitude of almost zero by definition, and the planets stay within a few degrees of the ecliptic). But the Moon can reach ±5.1°, and Pluto up to ±17°. In these cases, ignoring the latitude can shift the planet into a different house, especially near the cusps.

If you pass `lat_body=0.0`, you are calculating the house based purely on longitude — as if the planet were exactly on the ecliptic. This is the traditional method used by many software programs. But if you want the correct three-dimensional calculation, use the true ecliptic latitude of the planet.

### Gauquelin Sectors

Michel Gauquelin (1928–1991) was a French psychologist and statistician who studied whether planetary positions at birth correlated with professions. He discovered that certain planets tend to be found immediately after rising or immediately after culmination in the charts of successful people (the "Mars effect" for athletes, the "Jupiter effect" for actors).

For this analysis, the space around the observer is divided into **36 sectors** instead of 12. You can calculate them using:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 14.5)
lat, lon = 41.9028, 12.4964

# Gauquelin sector of Mars
sector = ephem.gauquelin_sector(
    jd, ephem.MARS, 0,
    geopos=(lon, lat, 0.0)
)

num = int(sector)
print(f"Mars is in Gauquelin sector {num}")
print(f"Sector 1 = rising, 10 = culmination, 19 = setting, 28 = IC")
```

```
Mars is in Gauquelin sector 18
Sector 1 = rising, 10 = culmination, 19 = setting, 28 = IC
```

Sectors are counted clockwise from the Ascendant: sector 1 is the rising point (Ascendant), 10 is the upper culmination (MC), 19 is the setting point (Descendant), and 28 is the lower culmination (IC). The most significant "Gauquelin zones" are the sectors immediately after rising (1–3) and immediately after culmination (10–12).

---

## 7.5 Extreme Latitudes and the Polar Circle

### The Problem

Up until now, we have worked with Rome (41.9° N), a latitude where all house systems work without issues. But what happens if you need to calculate the birth chart for someone born in Tromsø (69.6° N), Reykjavík (64.1° N), or Murmansk (68.9° N)?

The problem is geometric. Systems like Placidus and Koch divide the sky based on the **time** a degree of the ecliptic takes to travel from the horizon to the meridian. But beyond the polar circle (~66.56° latitude), some parts of the ecliptic **never rise or set** — they remain either always above or always below the horizon. If a point never rises, you cannot measure "how long it takes to rise", and the calculation becomes impossible.

The exact threshold depends on the obliquity of the ecliptic (about 23.44° in the current epoch):

- **Polar threshold** = 90° − obliquity ≈ 90° − 23.44° = **66.56°**
- Above this latitude, Placidus (`P`), Koch (`K`), and Gauquelin (`G`) **do not work**

### Which Systems Work Where

Not all systems have this problem. Here is the situation:

**Systems that fail beyond the polar circle** — These are based on the semi-arc (the time between rising and culmination), which does not exist for circumpolar points:

- **Placidus** (`P`) — the most widespread, but does not work beyond ~66.56°
- **Koch** (`K`) — same problem
- **Gauquelin** (`G`) — the 36 sectors share the same limit

**Systems that are unstable at extreme latitudes (beyond 80°)** — They work technically but can yield numerically imprecise or "squashed" results:

- Campanus (`C`), Regiomontanus (`R`), Polich/Page (`T`), Alcabitius (`B`), Horizon (`H`), Krusinski (`U`), APC (`Y`), Carter (`F`)

**Systems stable at all latitudes** — These use geometric methods that do not depend on rising and setting:

- Equal (`E`), Whole Sign (`W`), Porphyry (`O`), Morinus (`M`), Meridian (`X`), Vehlow (`V`), Axial Rotation (`N`)

You can verify the situation for any latitude using `get_extreme_latitude_info`:

```python
import libephemeris as ephem

# Tromsø, Norway — inside the polar circle
info = ephem.get_extreme_latitude_info(69.6)

print(f"Latitude: {info['latitude']}°")
print(f"Is extreme (>80°)?     {info['is_extreme']}")
print(f"Is polar (>{info['polar_threshold']:.1f}°)?   {info['is_polar_circle']}")
print(f"Systems that fail:       {info['affected_systems']}")
print(f"Unstable systems:        {info['unstable_systems']}")
print(f"Always stable systems:   {info['stable_systems']}")
```

```
Latitude: 69.6°
Is extreme (>80°)?     False
Is polar (>66.6°)?   True
Systems that fail:       ['P', 'K', 'G']
Unstable systems:        []
Always stable systems:   ['E', 'W', 'O', 'M', 'X', 'V', 'N']
```

### What Happens If You Try Placidus at the Pole?

If you call `houses()` with Placidus at a polar latitude, the library raises a `PolarCircleError`:

```python
import libephemeris as ephem
from libephemeris.exceptions import PolarCircleError

jd = ephem.julday(2024, 6, 21, 12.0)  # Summer solstice

try:
    cusps, ascmc = ephem.houses(jd, 69.6, 19.0, ord('P'))
except PolarCircleError as e:
    print(f"Error: {e}")
    print(f"Latitude: {e.latitude}°")
    print(f"Polar threshold: {e.threshold:.2f}°")
    print(f"System: {e.house_system}")
```

```
Error: swe_houses: Placidus house system cannot be calculated at latitude
  69.60°N (within Northern polar circle). Polar threshold for obliquity
  23.44° is ±66.56°.
Latitude: 69.6°
Polar threshold: 66.56°
System: P
```

### The Solution: `houses_with_fallback`

In practice, you do not want your program to crash with an error. You want a reasonable result. This is why `houses_with_fallback` exists: it tries the system you requested, and if it fails, it automatically uses an alternative system (by default Porphyry, which works everywhere).

```python
import libephemeris as ephem

jd = ephem.julday(2024, 6, 21, 12.0)

# Tromsø — Placidus will fail, fallback uses Porphyry
cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
    jd, 69.6, 19.0,
    hsys=ord('P'),
    fallback_hsys=ord('O')   # Porphyry as alternative
)

if used_fallback:
    print(f"Warning: {warning}")
    print("Using Porphyry instead of Placidus")
else:
    print("Placidus calculated normally")

# Cusps are available regardless
signs = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

for i in range(12):
    lon_c = cusps[i]
    sign = signs[int(lon_c / 30)]
    deg = lon_c % 30
    print(f"House {i+1:2d}: {deg:5.1f}° {sign}")
```

```
Warning: Placidus house system unavailable at latitude 69.60°
  (polar circle threshold: 66.56°). Using Porphyry as fallback.
Using Porphyry instead of Placidus
House  1:   9.6° Lib
House  2:  12.3° Sco
House  3:  15.0° Sgr
House  4:  17.7° Cap
House  5:  15.0° Aqr
House  6:  12.3° Psc
House  7:   9.6° Ari
House  8:  12.3° Tau
House  9:  15.0° Gem
House 10:  17.7° Cnc
House 11:  15.0° Leo
House 12:  12.3° Vir
```

The function returns four values:

- **cusps** — the 12 cusps (calculated with the primary system or the fallback)
- **ascmc** — the 8 angles (Ascendant, MC, ARMC, Vertex, etc.)
- **used_fallback** — `True` if it had to use the alternative system
- **warning** — a warning message explaining what happened, or `None` if everything went well

### Practical Advice

If your software needs to work for users all over the world, always use `houses_with_fallback` instead of `houses`. This way:

- For normal latitudes (the vast majority of cases), you get exactly the system you requested
- For polar latitudes, you get a reasonable result with Porphyry, instead of a crash
- The `used_fallback` flag allows you to inform the user that the result uses a different system

Alternatively, if you work with Vedic or Hellenistic astrology, directly use the Whole Sign (`W`) or Equal (`E`) system — they work everywhere without needing a fallback.

---

## Summary

In this chapter, we explored astrological houses, from theory to practice.

**Key Concepts:**

- **Houses** divide the local sky into 12 sectors based on date, time, and location — unlike signs, which depend only on the date.
- The **Ascendant** (cusp of the 1st house) is the degree of the ecliptic rising in the East; the **Midheaven** (cusp of the 10th) is the degree that culminates at the meridian.
- There are over 20 house systems because there is no unequivocally correct way to divide a 3D sphere into 12 sectors; each system uses a different geometric criterion.
- At **polar latitudes** (beyond ~66.56°), Placidus, Koch, and Gauquelin do not work; systems like Porphyry, Equal, and Whole Sign work everywhere.
- A **planet's position in the houses** depends not only on its longitude but also on its ecliptic latitude.

**Functions Introduced:**

- `houses(jd, lat, lon, hsys)` — calculates the 12 cusps and angles for a given house system
- `house_pos(armc, lat, obliquity, hsys, lon, lat_body)` — determines which house a planet falls into (returns a decimal value where the integer part is the house number)
- `house_name(hsys)` — returns the readable name of a house system (e.g., `house_name(ord('P'))` → `"Placidus"`)
- `houses_with_fallback(jd, lat, lon, hsys, fallback_hsys)` — like `houses`, but with automatic fallback for polar latitudes
- `get_extreme_latitude_info(lat)` — returns a dictionary with information on which systems work at a given latitude
- `gauquelin_sector(jd, planet, method, geopos)` — calculates the Gauquelin sector (1–36) of a planet
