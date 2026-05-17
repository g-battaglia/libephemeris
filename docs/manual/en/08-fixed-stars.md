# Chapter 8 — Fixed Stars

## What you will learn

In this chapter, you will discover why "fixed" stars are not really fixed, how the library calculates their positions taking into account precession and proper motion, how to search for a star by name (even with spelling mistakes), and how to use fixed stars in astrological calculations.

---

## 8.1 Stars are not so fixed

The ancients called them "fixed stars" to distinguish them from planets ("wandering stars"), which visibly move from night to night. But stars move too — just much, much more slowly.

### Precession: the whole sky rotates

The largest effect is not the movement of individual stars, but that of our **reference point**. The Earth's axis wobbles slowly like a spinning top, completing a cycle in about 26,000 years. This phenomenon — the **precession of the equinoxes** — causes the vernal point (0° Aries) to shift by about 50 arcseconds per year relative to the stars.

In practice, all stars shift by about **1° every 72 years** in ecliptic longitude. Regulus, which 2000 years ago was in the heart of Leo, is today at almost 0° Virgo.

### Proper motion: every star has its path

Besides precession (which moves *all* stars in the same way), each star has its own **proper motion** — its actual displacement in space. The speed depends on the star's distance and actual velocity:

- **Sirius**: 1.33"/year — in 2000 years it has moved by almost half a degree
- **Arcturus**: 2.28"/year — one of the fastest proper motions among bright stars
- **Regulus**: 0.25"/year — relatively "still"
- **Barnard's Star**: 10.36"/year — the absolute record, but it is invisible to the naked eye

For the bright stars used in astrology, proper motion is almost always negligible over scales of a few centuries. But for historical calculations (ancient birth charts, stars in ancient Egypt) it becomes important.

The library accounts for both precession and proper motion when calculating a star's position. You can also manually propagate proper motion with `propagate_proper_motion`:

```python
import libephemeris as ephem

# Propagate Sirius's proper motion from the Hipparcos catalog (J1991.25) to J2000
# RA and Dec in degrees, proper motion in arcsec/year
ra_sirio = 101.2872  # RA J1991.25 in degrees
dec_sirio = -16.7161  # Dec J1991.25 in degrees
pm_ra = -0.5461   # proper motion in RA (includes cos(dec)), arcsec/year
pm_dec = -1.2232   # proper motion in Dec, arcsec/year

J1991_25 = 2448349.0625  # Hipparcos catalog epoch
J2000 = 2451545.0         # J2000.0 epoch

ra_2000, dec_2000 = ephem.propagate_proper_motion(
    ra_sirio, dec_sirio,
    pm_ra, pm_dec,
    J1991_25, J2000
)

print(f"Sirius at J2000.0: RA = {ra_2000:.4f}°, Dec = {dec_2000:.4f}°")
```

```
Sirius at J2000.0: RA = 101.2858°, Dec = -16.7191°
```

In daily practice, you will almost never need this function: `fixstar2_ut` does everything automatically — precession, proper motion, nutation, and aberration included.

---

## 8.2 The star catalog

The library includes an internal catalog of about 100 stars, based on data from ESA's **Hipparcos** satellite (1989–1993). Hipparcos measured the position, proper motion, and magnitude of over 118,000 stars with unprecedented precision (about 1 thousandth of an arcsecond).

The catalog covers all the important stars for astrology and observational astronomy: the 15 Behenian fixed stars, the four "royal" stars (Aldebaran, Regulus, Antares, Fomalhaut), the brightest stars in the sky, and the most used zodiacal stars.

### Calculating a star's position

The main function is `fixstar2_ut`. Pass it the star's name, the Julian Day, and the calculation flags:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Position of Regulus
star_name, pos, flag, err = ephem.fixstar2_ut("Regulus", jd, ephem.FLG_SPEED)

if not err:
    lon = pos[0]  # ecliptic longitude
    lat = pos[1]  # ecliptic latitude
    vel = pos[3]  # speed in longitude (degrees/day)

    signs = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
             "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]
    sign = signs[int(lon / 30)]
    deg = lon % 30

    print(f"{star_name}")
    print(f"Position: {deg:.2f}° {sign}")
    print(f"Latitude: {lat:.4f}°")
    print(f"Speed: {vel:.6f}°/day")
else:
    print(f"Error: {err}")
```

```
Regulus,alLeo
Position: 0.17° Vir
Latitude: 0.4657°
Speed: -0.000071°/day
```

The function returns four values:

- **star_name** — the full name of the star in the `"Name,Nomenclature"` format (e.g. `"Regulus,alLeo"` where `alLeo` stands for "alpha Leonis")
- **pos** — a tuple of 6 values: `(longitude, latitude, distance, lon_speed, lat_speed, dist_speed)`. The distance is fixed at 100,000 AU (stars are effectively at an infinite distance for our calculations)
- **flag** — the returned flags (usually identical to those passed)
- **err** — an error string, empty if everything went well

### The main stars of astrology

Here are the most used stars in astrological tradition, with their approximate position in 2024:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

stars = [
    "Aldebaran",    # Bull's Eye — "Watcher of the East"
    "Regulus",      # Lion's Heart — "Watcher of the North"
    "Antares",      # Scorpion's Heart — "Watcher of the West"
    "Fomalhaut",    # Mouth of the Southern Fish — "Watcher of the South"
    "Sirius",       # The brightest star in the sky
    "Algol",        # The "Demon Star" — beta Persei
    "Spica",        # The ear of wheat of Virgo
    "Betelgeuse",   # Orion's shoulder
    "Rigel",        # Orion's foot
    "Vega",         # The star of Lyra
    "Polaris",      # The North Star
    "Canopus",      # The second brightest in the sky
]

signs = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

for star in stars:
    star_name, pos, flag, err = ephem.fixstar2_ut(star, jd, 0)
    if not err:
        lon = pos[0]
        sign = signs[int(lon / 30)]
        deg = lon % 30
        short_name = star_name.split(",")[0]
        print(f"{short_name:12s}  {deg:5.1f}° {sign}")
```

```
Aldebaran      10.1° Gem
Regulus         0.2° Vir
Antares        10.1° Sgr
Fomalhaut       4.2° Psc
Sirius         14.4° Cnc
Algol          26.5° Tau
Spica          24.2° Lib
Betelgeuse     29.1° Gem
Rigel          17.2° Gem
Vega           15.7° Cap
Polaris        28.9° Gem
Canopus        15.3° Cnc
```

### Magnitude: how bright is it?

**Magnitude** indicates the apparent brightness of a star. The scale is counterintuitive: **lower** numbers mean **brighter** stars. Sirius has a magnitude of −1.46 (very bright), while stars barely visible to the naked eye have a magnitude of about 6.

```python
import libephemeris as ephem

# Magnitude of some stars
for star in ["Sirius", "Vega", "Polaris", "Algol"]:
    star_name, mag, err = ephem.fixstar2_mag(star)
    if not err:
        short_name = star_name.split(",")[0]
        print(f"{short_name:12s}  magnitude {mag:+.2f}")
```

```
Sirius        magnitude -1.46
Vega          magnitude +0.03
Polaris       magnitude +1.98
Algol         magnitude +2.12
```

### Equatorial coordinates

If you need Right Ascension and Declination (to point a telescope, for example), add the `FLG_EQUATORIAL` flag:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

nome, pos, flag, err = ephem.fixstar2_ut(
    "Sirius", jd, ephem.FLG_EQUATORIAL
)

if not err:
    ra = pos[0]   # Right Ascension in degrees
    dec = pos[1]  # Declination in degrees

    # Convert RA to hours, minutes, seconds
    ra_ore = ra / 15.0
    ore = int(ra_ore)
    minuti = int((ra_ore - ore) * 60)
    secondi = ((ra_ore - ore) * 60 - minuti) * 60

    print(f"Sirius: RA = {ore}h {minuti}m {secondi:.1f}s, Dec = {dec:+.4f}°")
```

```
Sirius: RA = 6h 46m 12.6s, Dec = -16.7520°
```

---

## 8.3 Fixed stars in astrology

In astrological tradition, fixed stars have a special role. They are not used like planets (with daily aspects and transits), but are considered significant when a planet or an angle of the birth chart is **in conjunction** — meaning at the same ecliptic longitude — with an important star.

### The four royal stars

The most powerful are the four "royal stars" or "watchers of the sky", one for each quadrant of the zodiac:

- **Aldebaran** (~10° Gemini in 2024) — the Bull's Eye, watcher of the vernal equinox. Associated with courage, military honor, and success through integrity.

- **Regulus** (~0° Virgo) — the Lion's Heart, watcher of the summer solstice. The "most astrological" star: associated with power, fame, and leadership. It is almost exactly on the ecliptic (latitude 0.46°), which makes it particularly significant.

- **Antares** (~10° Sagittarius) — the Scorpion's Heart, watcher of the autumnal equinox. "The rival of Mars" (anti-Ares), associated with intensity, passion, and destruction-regeneration.

- **Fomalhaut** (~4° Pisces) — the Mouth of the Southern Fish, watcher of the winter solstice. Associated with idealism, spirituality, and dreams.

### Conjunction with a planet

In stellar astrology, the **orb** (the maximum distance to consider an influence active) is much tighter than for planets: generally **1°** for first-magnitude stars, **30'** for others. Here is how to check if a planet is conjunct a star today:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Stars to check (the four royal stars + some important ones)
stars_to_check = [
    "Aldebaran", "Regulus", "Antares", "Fomalhaut",
    "Sirius", "Spica", "Algol", "Vega"
]

planets = [
    (ephem.SUN,     "Sun"),
    (ephem.MOON,    "Moon"),
    (ephem.MERCURY, "Mercury"),
    (ephem.VENUS,   "Venus"),
    (ephem.MARS,    "Mars"),
    (ephem.JUPITER, "Jupiter"),
    (ephem.SATURN,  "Saturn"),
]

# Calculate planets' positions
planet_lons = {}
for body_id, planet_name in planets:
    pos, _ = ephem.calc_ut(jd, body_id, ephem.FLG_SPEED)
    planet_lons[planet_name] = pos[0]  # ecliptic longitude

# Check conjunctions
orb = 1.0  # 1 degree of tolerance

for star in stars_to_check:
    star_name, pos_s, _, err = ephem.fixstar2_ut(star, jd, 0)
    if err:
        continue

    lon_s = pos_s[0]
    short_name = star_name.split(",")[0]

    for planet_name, lon_p in planet_lons.items():
        # Angular distance (handling the 360°→0° "jump")
        diff = abs(lon_p - lon_s)
        if diff > 180:
            diff = 360 - diff

        if diff <= orb:
            print(f"★ {planet_name} conjunct {short_name}! "
                  f"(distance: {diff:.2f}°)")
```

```
No conjunctions within 1° found on this date.
```

---

## 8.4 Searching for a star: names, designations, and fuzzy search

The library is very flexible in how it accepts star names. You can use:

**Traditional name** — the simplest and most intuitive way:

```python
ephem.fixstar2_ut("Sirius", jd, 0)
ephem.fixstar2_ut("Regulus", jd, 0)
ephem.fixstar2_ut("Betelgeuse", jd, 0)
```

**Bayer designation** — the Greek letter + the constellation, used in astronomical catalogs. Alpha Leonis = the brightest star in Leo:

```python
ephem.fixstar2_ut("Alpha Leonis", jd, 0)     # = Regulus
ephem.fixstar2_ut("Beta Persei", jd, 0)       # = Algol
ephem.fixstar2_ut("Alpha Canis Majoris", jd, 0)  # = Sirius
```

**Nomenclature code** — the abbreviated version of the Bayer designation, in the format used internally (e.g. `alLeo` = alpha Leonis):

```python
ephem.fixstar2_ut("alLeo", jd, 0)   # = Regulus
ephem.fixstar2_ut("bePer", jd, 0)   # = Algol
```

**HIP number** — the Hipparcos catalog number, for those working with astronomical data:

```python
ephem.fixstar2_ut("HIP 49669", jd, 0)  # = Regulus
ephem.fixstar2_ut("32349", jd, 0)       # = Sirius (HIP 32349)
```

**Partial search** — if you only remember the beginning of the name:

```python
ephem.fixstar2_ut("Reg", jd, 0)    # finds Regulus (if not ambiguous)
ephem.fixstar2_ut("Alde", jd, 0)   # finds Aldebaran
```

### Phonetic search: forgives spelling mistakes

Star names come from Arabic, Latin, and Greek, and transliterations vary. The library includes a **phonetic search** system that finds the right star even if you misspell it:

```python
# All these find Betelgeuse:
ephem.fixstar2_ut("Betelgeuse", jd, 0)   # correct spelling
ephem.fixstar2_ut("Betelgeuze", jd, 0)   # with z (German variant)
ephem.fixstar2_ut("Betelgeux", jd, 0)    # truncated

# And these find Fomalhaut:
ephem.fixstar2_ut("Fomalhaut", jd, 0)    # correct spelling
ephem.fixstar2_ut("Formalhaut", jd, 0)   # with an extra r
```

The phonetic system normalizes the name by removing double consonants, unifying similar vowels, and comparing the "consonant skeleton" of the name. It works well for common spelling variants, but returns an error if the search is ambiguous (multiple stars match the same pattern).

### The returned name format

When the search is successful, the first returned value is the full name in the `"Name,Nomenclature"` format:

```python
nome, pos, flag, err = ephem.fixstar2_ut("Sirius", jd, 0)
print(nome)  # "Sirius,alCMa" → alpha Canis Majoris

nome, pos, flag, err = ephem.fixstar2_ut("Algol", jd, 0)
print(nome)  # "Algol,bePer" → beta Persei
```

```
Sirius,alCMa
Algol,bePer
```

The part after the comma is the abbreviated Bayer designation: `al` = alpha, `be` = beta, `ga` = gamma, and so on. The constellation code follows the standard three-letter IAU abbreviation (Leo, Per, CMa, Vir, Ori...).

---

## Summary

In this chapter, we learned how to work with fixed stars.

**Key concepts:**

- "Fixed" stars move: **precession** shifts all stars by ~50"/year in ecliptic longitude, and **proper motion** adds an individual displacement for each star
- The library includes a catalog of ~100 stars based on Hipparcos satellite data, with sub-milliarcsecond precision
- In astrology, fixed stars matter when a planet or angle is in **conjunction** (within ~1°) with an important star
- The four royal stars (Aldebaran, Regulus, Antares, Fomalhaut) are considered the most powerful
- **Magnitude** indicates brightness: lower numbers = brighter star (Sirius = −1.46, naked eye limit ≈ 6)
- The search is flexible: traditional name, Bayer designation, HIP number, partial search, and even phonetic search for misspelled names

**Functions introduced:**

- `fixstar2_ut(nome, jd, flag)` — calculates the ecliptic (or equatorial with `FLG_EQUATORIAL`) position of a star, accepting names in many formats
- `fixstar2_mag(nome)` — returns the visual magnitude of a star
- `propagate_proper_motion(ra, dec, pm_ra, pm_dec, from_jd, to_jd)` — manually propagates a star's proper motion between two epochs
