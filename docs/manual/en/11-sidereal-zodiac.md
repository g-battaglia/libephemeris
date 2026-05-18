# Chapter 11 — The Sidereal Zodiac and Ayanamshas

## What you will learn

In this chapter, you will discover why two different zodiacs exist (tropical and sidereal), what caused them to diverge over millennia, what the ayanamsha is, why there are over 40 different ways to calculate it, and how to use the library to work with the sidereal zodiac used in Indian Vedic astrology.

---

## 11.1 Two zodiacs: tropical and sidereal

### The same sky, two different measures

In Chapter 1, we saw that the zodiac is a division of the ecliptic into 12 sectors of 30° each. But where does the counting start? Here lies the problem — and one of the deepest divisions in the history of astrology.

The **tropical zodiac** (used in the West) starts the count from the **vernal equinox** (spring equinox): the point where the Sun crosses the celestial equator moving north, around March 20. This point is defined as 0° Aries, regardless of which stars are in that direction.

The **sidereal zodiac** (used in India and Vedic astrology) starts the count from a point linked to the **fixed stars**. The zero point is anchored to a reference star (usually Spica, fixed at 0° Libra = 180°).

About 2000 years ago, the two zodiacs coincided almost perfectly. But since then, they have shifted relative to each other, and today the difference is about **24 degrees**.

### Why they diverge: the precession of the equinoxes

The cause is **precession**: the Earth's axis slowly wobbles like a spinning top, completing a cycle in about 26,000 years. This causes the vernal equinox point to shift by about 50 arcseconds per year *relative to the stars*.

In practice:
- The tropical zodiac "follows" the equinox, and therefore shifts relative to the stars
- The sidereal zodiac "follows" the stars, and therefore shifts relative to the equinox

The practical effect is significant: if in your tropical (Western) natal chart the Sun is at 15° Aries, in the sidereal (Vedic) chart it will be at about 21° Pisces — almost an entire sign of difference.

```python
import libephemeris as ephem

# Comparison: Sun's position in tropical and sidereal
jd = ephem.julday(2024, 4, 8, 12.0)

# Tropical position (default)
pos_trop, _ = ephem.calc_ut(jd, ephem.SUN, ephem.FLG_SPEED)

# Lahiri Ayanamsha
ephem.set_sid_mode(ephem.SIDM_LAHIRI)
ayan = ephem.get_ayanamsa_ut(jd)

# Sidereal position = tropical - ayanamsha
lon_sid = (pos_trop[0] - ayan) % 360

signs = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

s_trop = signs[int(pos_trop[0] / 30)]
g_trop = pos_trop[0] % 30

s_sid = signs[int(lon_sid / 30)]
g_sid = lon_sid % 30

print(f"Tropical Sun: {g_trop:.1f}° {s_trop}")
print(f"Sidereal Sun:  {g_sid:.1f}° {s_sid}")
print(f"Ayanamsha (Lahiri): {ayan:.4f}°")
```

```
Tropical Sun: 19.1° Ari
Sidereal Sun:  24.9° Psc
Ayanamsha (Lahiri): 24.1961°
```

---

## 11.2 What is the ayanamsha

The **ayanamsha** (from Sanskrit *ayana* = solstice, *amsha* = portion) is simply the difference in degrees between the tropical zero point and the sidereal zero point at a given moment:

**sidereal longitude = tropical longitude − ayanamsha**

Today the ayanamsha is about 24°, and grows by about 50.3" per year (the precession rate). In a few millennia, the two zodiacs will coincide again — and then diverge again in the other direction.

### Calculating the ayanamsha

```python
import libephemeris as ephem

ephem.set_sid_mode(ephem.SIDM_LAHIRI)

# Ayanamsha in different epochs
for year in [0, 500, 1000, 1500, 2000, 2024, 2100]:
    jd = ephem.julday(year, 1, 1, 0.0)
    ayan = ephem.get_ayanamsa_ut(jd)
    print(f"Year {year:5d}: ayanamsha = {ayan:+7.2f}°")
```

```
Year     0: ayanamsha = +356.04°
Year   500: ayanamsha =   +2.97°
Year  1000: ayanamsha =   +9.92°
Year  1500: ayanamsha =  +16.88°
Year  2000: ayanamsha =  +23.86°
Year  2024: ayanamsha =  +24.19°
Year  2100: ayanamsha =  +25.25°
```

You will notice that the ayanamsha is negative in antiquity (the sidereal zodiac was "ahead" relative to the tropical one) and positive today (the tropical is "ahead").

### Using the FLG_SIDEREAL flag

Instead of calculating the ayanamsha and subtracting it manually, you can ask the library directly to return positions in sidereal coordinates:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Set the sidereal system
ephem.set_sid_mode(ephem.SIDM_LAHIRI)

# Calculate directly in sidereal
pos_sid, _ = ephem.calc_ut(
    jd, ephem.SUN,
    ephem.FLG_SIDEREAL | ephem.FLG_SPEED
)

signs = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

sign = signs[int(pos_sid[0] / 30)]
degrees = pos_sid[0] % 30

print(f"Sidereal Sun (Lahiri): {degrees:.2f}° {sign}")
```

```
Sidereal Sun (Lahiri): 24.95° Psc
```

The flag `FLG_SIDEREAL` (65536) can be combined with any other flag: `FLG_SPEED`, `FLG_EQUATORIAL`, etc.

---

## 11.3 The schools of ayanamsha

### The problem: nobody agrees

If the ayanamsha were a simple matter of fact, there would be no discussion. But the problem is: **where exactly was the zero point 2000 years ago?** There was no celestial GPS. Different astronomers and astrologers have proposed slightly different answers, based on different stars, different ancient texts, or different calculations. The difference between the most popular schools is a few degrees — enough to shift a planet into a different sign.

The library supports **47 ayanamsha systems**. Here are the most important ones:

### Indian Ayanamshas (the most used)

**Lahiri** (`SIDM_LAHIRI`, value 1) is the most widespread: it is the official standard of the Indian government, adopted in 1955 by the Indian Calendar Reform Committee (N.C. Lahiri). It fixes Spica (Citra in Sanskrit) at 180° of sidereal longitude. It is used by the vast majority of Vedic astrologers.

**Krishnamurti** (`SIDM_KRISHNAMURTI`, value 5) is used in the KP (Krishnamurti Paddhati) system, a very popular predictive method in South India. It is very similar to Lahiri, with a difference of only a few arcminutes.

**Raman** (`SIDM_RAMAN`, value 3) was proposed by B.V. Raman, one of the most influential Indian astrologers of the 20th century. It differs from Lahiri by about 1.5°.

**Lahiri variants**: there are also `SIDM_LAHIRI_1940` (43), `SIDM_LAHIRI_VP285` (44), and `SIDM_LAHIRI_ICRC` (46), which differ by a few arcseconds and reflect different interpretations of the original data.

### Western Sidereal Ayanamsha

**Fagan-Bradley** (`SIDM_FAGAN_BRADLEY`, value 0) was developed by Cyril Fagan and Donald Bradley in the 1950s for Western sidereal astrology. It differs from Lahiri by about 1°. It is rarely used outside a narrow circle of Western sidereal astrologers.

### Ayanamshas based on true stars

**True Citra** (`SIDM_TRUE_CITRA`, value 27) fixes the *true* position (with proper motion) of Spica at exactly 180°. Unlike Lahiri, which uses a polynomial formula calculated once and for all, True Citra recalculates the real position of Spica at each date, tracking its proper motion.

**True Revati** (`SIDM_TRUE_REVATI`, value 28) fixes the star Revati (zeta Piscium) at 359°50'.

**True Pushya** (`SIDM_TRUE_PUSHYA`, value 29) fixes the star Pushya (delta Cancri) at 106°.

### Galactic Ayanamshas

For those looking for a "cosmic" zero point, there are several systems based on the position of the galactic center or the galactic equator:

**Galactic Center 0 Sag** (`SIDM_GALCENT_0SAG`, value 17) places the galactic center at 0° Sagittarius.

### Babylonian Ayanamshas

For historical research on Babylonian astronomy: `SIDM_BABYL_KUGLER1` (9), `SIDM_BABYL_KUGLER2` (10), `SIDM_BABYL_KUGLER3` (11), `SIDM_BABYL_HUBER` (12), `SIDM_BABYL_ETPSC` (13), `SIDM_BABYL_BRITTON` (38).

### Comparison between ayanamshas

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

systems = [
    (ephem.SIDM_LAHIRI,         "Lahiri"),
    (ephem.SIDM_FAGAN_BRADLEY,  "Fagan-Bradley"),
    (ephem.SIDM_RAMAN,          "Raman"),
    (ephem.SIDM_KRISHNAMURTI,   "Krishnamurti"),
    (ephem.SIDM_TRUE_CITRA,     "True Citra"),
    (ephem.SIDM_TRUE_REVATI,    "True Revati"),
    (ephem.SIDM_GALCENT_0SAG,   "Galactic 0 Sag"),
]

print(f"Ayanamsha on {8}/04/2024:\n")
for mode, name in systems:
    ephem.set_sid_mode(mode)
    ayan = ephem.get_ayanamsa_ut(jd)
    print(f"  {name:20s}  {ayan:.4f}°")

# The name of a system
name = ephem.get_ayanamsa_name(ephem.SIDM_LAHIRI)
print(f"\nName: {name}")
```

```
Ayanamsha on 8/04/2024:

  Lahiri                24.1961°
  Fagan-Bradley         25.0793°
  Raman                 22.7498°
  Krishnamurti          24.0993°
  True Citra            24.1843°
  True Revati           20.3779°
  Galactic 0 Sag        27.1850°

Name: Lahiri
```

---

## 11.4 Custom Ayanamsha

If none of the 47 predefined systems meets your needs, you can define your own custom ayanamsha with `SIDM_USER`:

```python
import libephemeris as ephem

# Define a custom ayanamsha:
# "On January 1st, 2000 (J2000.0), the ayanamsha is 23.5°"
J2000 = 2451545.0

ephem.set_sid_mode(
    ephem.SIDM_USER,
    t0=J2000,       # reference epoch
    ayan_t0=23.5     # ayanamsha value at that epoch
)

# Now calculate the ayanamsha for any date
jd = ephem.julday(2024, 4, 8, 12.0)
ayan = ephem.get_ayanamsa_ut(jd)
print(f"Custom ayanamsha in 2024: {ayan:.4f}°")

# And use FLG_SIDEREAL for positions
pos, _ = ephem.calc_ut(jd, ephem.SUN,
    ephem.FLG_SIDEREAL | ephem.FLG_SPEED)
print(f"Custom sidereal Sun: {pos[0]:.4f}°")
```

```
Custom ayanamsha in 2024: 23.8390°
Custom sidereal Sun: 355.3029°
```

The library uses the IAU 2006 precession polynomials to calculate how the ayanamsha changes over time starting from the value you specified.

---

## 11.5 Practical calculations in sidereal

### Vedic birth chart

Here is how to calculate a complete birth chart in Lahiri sidereal — the format used in Vedic astrology (Jyotish):

```python
import libephemeris as ephem

# Date: August 15, 1947, 00:00 IST — Independence of India
# IST = UT + 5:30, so UT = 18:30 on August 14
jd = ephem.julday(1947, 8, 14, 18.5)
lat, lon = 28.6139, 77.2090  # New Delhi

# Set Lahiri sidereal
ephem.set_sid_mode(ephem.SIDM_LAHIRI)

signs = ["Mesha", "Vrishabha", "Mithuna", "Karka",
         "Simha", "Kanya", "Tula", "Vrischika",
         "Dhanus", "Makara", "Kumbha", "Meena"]

# Planets in sidereal
planets = [
    (ephem.SUN,      "Surya  "),
    (ephem.MOON,     "Chandra"),
    (ephem.MARS,     "Mangal "),
    (ephem.MERCURY,  "Budha  "),
    (ephem.JUPITER,  "Guru   "),
    (ephem.VENUS,    "Shukra "),
    (ephem.SATURN,   "Shani  "),
    (ephem.MEAN_NODE,"Rahu   "),
]

print("Sidereal birth chart (Lahiri):\n")
for body_id, name in planets:
    pos, _ = ephem.calc_ut(jd, body_id,
        ephem.FLG_SIDEREAL | ephem.FLG_SPEED)
    sign = signs[int(pos[0] / 30)]
    degrees = pos[0] % 30
    print(f"  {name}  {degrees:5.1f}°  {sign}")

# Ketu = opposite to Rahu
rahu_pos, _ = ephem.calc_ut(jd, ephem.MEAN_NODE,
    ephem.FLG_SIDEREAL | ephem.FLG_SPEED)
ketu_lon = (rahu_pos[0] + 180) % 360
sign_k = signs[int(ketu_lon / 30)]
degrees_k = ketu_lon % 30
print(f"  Ketu    {degrees_k:5.1f}°  {sign_k}")
```

```
Sidereal birth chart (Lahiri):

  Surya    28.0°  Karka
  Chandra   4.0°  Karka
  Mangal    7.5°  Mithuna
  Budha    13.7°  Karka
  Guru     25.9°  Tula
  Shukra   22.6°  Karka
  Shani    20.5°  Karka
  Rahu      5.1°  Vrishabha
  Ketu      5.1°  Vrischika
```

### The Nakshatras

In Vedic astrology, the zodiac is also divided into **27 Nakshatras** (lunar mansions) of 13°20' each. Each Nakshatra has a name, a meaning, and a ruling planet. The Moon takes about one day to traverse each of them.

```python
import libephemeris as ephem

nakshatra_names = [
    "Ashwini", "Bharani", "Krittika", "Rohini", "Mrigashira",
    "Ardra", "Punarvasu", "Pushya", "Ashlesha", "Magha",
    "Purva Phalguni", "Uttara Phalguni", "Hasta", "Chitra",
    "Swati", "Vishakha", "Anuradha", "Jyeshtha", "Mula",
    "Purva Ashadha", "Uttara Ashadha", "Shravana", "Dhanishtha",
    "Shatabhisha", "Purva Bhadrapada", "Uttara Bhadrapada", "Revati"
]

# Rulers in the Vimshottari Dasha sequence
rulers = [
    "Ketu", "Venus", "Sun", "Moon", "Mars",
    "Rahu", "Jupiter", "Saturn", "Mercury",
    "Ketu", "Venus", "Sun", "Moon", "Mars",
    "Rahu", "Jupiter", "Saturn", "Mercury",
    "Ketu", "Venus", "Sun", "Moon", "Mars",
    "Rahu", "Jupiter", "Saturn", "Mercury"
]

jd = ephem.julday(2024, 4, 8, 12.0)
ephem.set_sid_mode(ephem.SIDM_LAHIRI)

# Moon's Nakshatra
moon, _ = ephem.calc_ut(jd, ephem.MOON,
    ephem.FLG_SIDEREAL | ephem.FLG_SPEED)

nak_num = int(moon[0] / (360 / 27))  # 13°20' per Nakshatra
nak_pos = moon[0] % (360 / 27)       # position in the Nakshatra

# Pada (quarter) — each Nakshatra has 4 padas of 3°20'
pada = int(nak_pos / (360 / 108)) + 1  # 1-4

print(f"Sidereal Moon: {moon[0]:.2f}°")
print(f"Nakshatra: {nakshatra_names[nak_num]} (n. {nak_num + 1})")
print(f"Pada: {pada}")
print(f"Ruler: {rulers[nak_num]}")
```

```
Sidereal Moon: 351.24°
Nakshatra: Revati (n. 27)
Pada: 2
Ruler: Mercury
```

The Moon's Nakshatra at birth is fundamental in Vedic astrology: it determines the **Dasha** (the system of planetary periods that governs the person's life).

---

## 11.6 The extended version: `get_ayanamsa_ex_ut`

For advanced calculations where you need to specify additional flags or want the return flag, use `get_ayanamsa_ex_ut`:

```python
import libephemeris as ephem

ephem.set_sid_mode(ephem.SIDM_LAHIRI)

jd = ephem.julday(2024, 4, 8, 12.0)

# Simple version
ayan_simple = ephem.get_ayanamsa_ut(jd)

# Extended version (also returns the flag)
retflag, ayan_ex = ephem.get_ayanamsa_ex_ut(jd, ephem.FLG_SWIEPH)

print(f"Ayanamsha (simple): {ayan_simple:.6f}°")
print(f"Ayanamsha (extended):   {ayan_ex:.6f}°")
```

```
Ayanamsha (simple): 24.196111°
Ayanamsha (extended):   24.196111°
```

---

## Summary

In this chapter, we explored the sidereal zodiac, which is fundamental to Vedic astrology.

**Key concepts:**

- There are two zodiacs: the **tropical** (zero point = equinox) and the **sidereal** (zero point linked to the stars). Today they differ by about 24°.
- The divergence is caused by the **precession of the equinoxes**: the Earth's axis wobbles with a period of ~26,000 years, shifting the vernal point relative to the stars.
- The **ayanamsha** is the difference in degrees between the two zodiacs: sidereal longitude = tropical longitude − ayanamsha.
- There are over **47 ayanamsha systems** because there is no agreement on where the zero point was 2000 years ago.
- **Lahiri** is the most widely used (Indian government standard), **Fagan-Bradley** is used for Western sidereal astrology, and **True Citra** uses the real position of Spica.
- In Vedic astrology, the **27 Nakshatras** divide the sidereal zodiac into sectors of 13°20', each with a meaning and a ruling planet.

**Functions introduced:**

- `set_sid_mode(sid_mode, t0=0.0, ayan_t0=0.0)` — sets the ayanamsha system. Use `SIDM_USER` with `t0` and `ayan_t0` for a custom ayanamsha.
- `get_ayanamsa_ut(jd)` — returns the ayanamsha in degrees for a date in UT.
- `get_ayanamsa_ex_ut(jd, flags)` — extended version that also returns the return flag.
- `get_ayanamsa_name(sid_mode)` — returns the readable name of an ayanamsha system (e.g. `"Lahiri"`).
- `FLG_SIDEREAL` (65536) — flag to add to `calc_ut` to get positions directly in sidereal coordinates.