# Chapter 13 — Hypothetical planets and Arabic parts

## What you will learn

In this chapter you will discover what Uranian planets are (mathematical points used in the Hamburg School), other hypothetical bodies like Transpluto and Vulcan, how to define custom fictitious orbits, and how to calculate Arabic parts — ancient formulas that combine planetary positions and angles.

---

## 13.1 The Uranian planets (Hamburg School)

In the 1920s, the German astrologer Alfred Witte founded the **Hamburg School** and postulated the existence of eight trans-Neptunian "planets". These bodies have never been observed — they do not exist physically. They are **mathematical points** with orbits defined by Keplerian elements, used exclusively in Uranian astrology and the technique of "midpoints".

The eight Uranian planets are:

- **Cupido** (`CUPIDO`, ID 40) — associated with family, art, and groups
- **Hades** (`HADES`, ID 41) — associated with the past, poverty, disease
- **Zeus** (`ZEUS`, ID 42) — associated with fire, machines, driving force
- **Kronos** (`KRONOS`, ID 43) — associated with authority, government, excellence
- **Apollon** (`APOLLON`, ID 44) — associated with expansion, science, commerce
- **Admetos** (`ADMETOS`, ID 45) — associated with depth, concentration, blocks
- **Vulkanus** (`VULKANUS`, ID 46) — associated with power, intensity, force
- **Poseidon** (`POSEIDON`, ID 47) — associated with mind, enlightenment, truth

### Calculating the position

Each Uranian planet has its dedicated function:

```python
import libephemeris as ephem

# Convert the date to TT (these work in TT, not UT)
jd_ut = ephem.julday(2024, 4, 8, 12.0)
delta_t = ephem.deltat(jd_ut)  # in days
jd_tt = jd_ut + delta_t

# Kronos position
pos = ephem.calc_kronos(jd_tt)
lon, lat, dist = pos[0], pos[1], pos[2]
vel = pos[3]  # velocity in degrees/day

signs = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]
sign = signs[int(lon / 30)]
degrees = lon % 30

print(f"Kronos: {degrees:.2f}° {sign}")
print(f"Velocity: {vel:.4f}°/day")
```

```
Kronos: 14.95° Cnc
Velocity: 0.0019°/day
```

You can also use the generic `calc_uranian_planet` function with the body ID:

```python
import libephemeris as ephem

jd_ut = ephem.julday(2024, 4, 8, 12.0)
jd_tt = jd_ut + ephem.deltat(jd_ut)

signs = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

uranians = [
    (ephem.CUPIDO,   "Cupido"),
    (ephem.HADES,    "Hades"),
    (ephem.ZEUS,     "Zeus"),
    (ephem.KRONOS,   "Kronos"),
    (ephem.APOLLON,  "Apollon"),
    (ephem.ADMETOS,  "Admetos"),
    (ephem.VULKANUS, "Vulkanus"),
    (ephem.POSEIDON, "Poseidon"),
]

for body_id, name in uranians:
    pos = ephem.calc_uranian_planet(body_id, jd_tt)
    sign = signs[int(pos[0] / 30)]
    degrees = pos[0] % 30
    print(f"{name:10s}  {degrees:5.1f}° {sign}")
```

```
Cupido        7.5° Cap
Hades        12.9° Cnc
Zeus         25.1° Lib
Kronos       14.9° Cnc
Apollon       6.5° Sco
Admetos       3.3° Gem
Vulkanus      3.8° Leo
Poseidon     16.3° Sco
```

The Uranian planets can also be calculated with `calc_ut` using their IDs, but the dedicated functions directly accept JD in TT.

---

## 13.2 Other hypothetical bodies

Besides the Uranians, the library includes other hypothetical bodies that have been proposed throughout history but never observationally confirmed:

**Transpluto / Isis** (`ISIS`, ID 48) — A hypothetical planet beyond Pluto, postulated before the discovery of Eris. The orbit used is based on Theodore Landscheidt's proposal.

```python
import libephemeris as ephem

jd_tt = ephem.julday(2024, 4, 8, 12.0) + ephem.deltat(ephem.julday(2024, 4, 8, 12.0))
pos = ephem.calc_transpluto(jd_tt)
print(f"Transpluto: {pos[0]:.2f}°")
```

```
Transpluto: 153.47°
```

**Vulcan** (`VULCAN`, ID 55) — A hypothetical planet between Mercury and the Sun, searched for throughout the 19th century to explain anomalies in Mercury's orbit. Einstein's general relativity later explained those anomalies without the need for an additional planet — but the concept remains in esoteric astrology.

```python
pos = ephem.calc_vulcan(jd_tt)
print(f"Vulcan: {pos[0]:.2f}°")
```

```
Vulcan: 45.53°
```

**White Moon / Selena** (`WHITE_MOON`, ID 56) — The point diametrically opposite to the Black Moon Lilith (mean lunar apogee). It is not a physical body, but a symbolic point used in some schools of astrology as the "luminous complement" of Lilith.

```python
pos = ephem.calc_white_moon_position(jd_tt)
print(f"White Moon: {pos[0]:.2f}°")
```

```
White Moon: 350.92°
```

**Waldemath's Moon** (ID 58) — A hypothetical second satellite of the Earth, "observed" in 1898 by Georg Waldemath. Never confirmed.

**Proserpina** (ID 57) — Another hypothetical trans-Plutonian.

**Pickering's Planet X** — William Pickering's 1919 prediction for a planet beyond Neptune.

All these bodies are calculated with the generic `calc_hypothetical_position` function:

```python
import libephemeris as ephem

jd_tt = ephem.julday(2024, 4, 8, 12.0) + ephem.deltat(ephem.julday(2024, 4, 8, 12.0))

# Any hypothetical body given its ID
pos = ephem.calc_hypothetical_position(ephem.WALDEMATH, jd_tt)
print(f"Waldemath's Moon: {pos[0]:.2f}°")
```

```
Waldemath's Moon: 61.36°
```

---

## 13.3 Custom fictitious orbits

If you need a hypothetical body not included in the library, you can define your own orbit. The library uses a file format with orbital elements (compatible with the Swiss Ephemeris `seorbel.txt` format).

### Loading predefined orbits

The library includes a file of predefined fictitious orbits:

```python
import libephemeris as ephem

# Load predefined orbits
orbits = ephem.load_bundled_fictitious_orbits()
print(f"Loaded {len(orbits)} fictitious orbits")

# Search for a body by name
body = ephem.get_orbital_body_by_name(orbits, "Cupido")
if body:
    print(f"Found: {body.name}")
    print(f"Semi-major axis: {body.semi_axis:.2f} AU")
```

```
Loaded 24 fictitious orbits
Found: Cupido
Semi-major axis: 41.00 AU
```

### Loading a custom file

```python
import libephemeris as ephem

# Load orbits from a custom file
orbits = ephem.parse_orbital_elements("/path/to/my/file.csv")

# Calculate the position of a body
jd_tt = ephem.julday(2024, 4, 8, 12.0) + ephem.deltat(ephem.julday(2024, 4, 8, 12.0))

body = ephem.get_orbital_body_by_name(orbits, "MyBody")
if body:
    pos = ephem.calc_orbital_position(body, jd_tt)
    print(f"MyBody: lon={pos[0]:.2f}°, lat={pos[1]:.2f}°")
```

---

## 13.4 Arabic parts (Lots)

The **Arabic parts** (in Greek *kleros*, "lots") are among the oldest techniques in astrology. They date back to Hellenistic astrology (2nd–3rd century AD) and were extensively developed by Persian and Arabic astrologers in the Middle Ages.

### How they work

Each Arabic part is a calculated point combining three positions: generally the Ascendant and two planets. The base formula is:

**Part = Ascendant + Planet A − Planet B**

The result (normalized to 0°–360°) is a point on the ecliptic with a specific meaning.

### The Part of Fortune

The most famous is the **Part of Fortune** (*Pars Fortunae*), associated with prosperity, material fortune, and physical well-being:

- **By day** (Sun above the horizon): Fortune = ASC + Moon − Sun
- **By night** (Sun below the horizon): Fortune = ASC + Sun − Moon

The formula is inverted between day and night because the "sect luminary" (the Sun by day, the Moon by night) acts as the starting point.

### Calculating all parts

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 14.5)
lat, lon = 41.9028, 12.4964  # Rome

# Calculate the necessary positions
cusps, ascmc = ephem.houses(jd, lat, lon, ord('P'))
sun, _ = ephem.calc_ut(jd, ephem.SUN, ephem.FLG_SPEED)
moon, _ = ephem.calc_ut(jd, ephem.MOON, ephem.FLG_SPEED)
mercury, _ = ephem.calc_ut(jd, ephem.MERCURY, ephem.FLG_SPEED)
venus, _ = ephem.calc_ut(jd, ephem.VENUS, ephem.FLG_SPEED)

# Prepare the positions
positions = {
    "Asc": ascmc[0],
    "Sun": sun[0],
    "Moon": moon[0],
    "Mercury": mercury[0],
    "Venus": venus[0],
}

# Calculate all Arabic parts
parts = ephem.calc_all_arabic_parts(
    positions,
    jd=jd,
    geo_lat=lat,
    geo_lon=lon,
)

signs = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

english_names = {
    "Pars_Fortunae": "Part of Fortune",
    "Pars_Spiritus": "Part of Spirit",
    "Pars_Amoris":   "Part of Love",
    "Pars_Fidei":    "Part of Faith",
}

for key, lon_part in parts.items():
    name = english_names.get(key, key)
    sign = signs[int(lon_part / 30)]
    degrees = lon_part % 30
    print(f"{name:22s}  {degrees:5.1f}° {sign}")
```

```
Part of Fortune           14.5° Vir
Part of Spirit            10.0° Vir
Part of Love              27.3° Leo
Part of Faith             20.2° Vir
```

The four calculated parts are:

- **Part of Fortune** (*Pars Fortunae*) — ASC + Moon − Sun (day) or ASC + Sun − Moon (night). Prosperity and material well-being.

- **Part of Spirit** (*Pars Spiritus*) — the inverse of the Part of Fortune. It represents the will, the spirit, and inner vocation.

- **Part of Love** (*Pars Amoris*) — ASC + Venus − Sun. Affectionate relationships and attraction.

- **Part of Faith** (*Pars Fidei*) — ASC + Mercury − Moon. Faith, trust, and beliefs.

---

## Summary

In this chapter we explored non-physical celestial bodies used in various astrological traditions.

**Key concepts:**

- The **Uranian planets** are eight mathematical points (Cupido, Hades, Zeus, Kronos, Apollon, Admetos, Vulkanus, Poseidon) with hypothetical orbits, used in the Hamburg School and in the midpoints technique
- Other **hypothetical bodies** (Transpluto, Vulcan, Waldemath's Moon) have been proposed historically but never confirmed; the library calculates them using defined orbital elements
- **Custom fictitious orbits** allow you to define any hypothetical body with its own orbital elements
- The **Arabic parts** are points calculated by the formula ASC + Planet A − Planet B, with the Part of Fortune being the most important

**Functions introduced:**

- `calc_cupido(jd_tt)`, `calc_hades(jd_tt)`, ... `calc_poseidon(jd_tt)` — calculate each Uranian planet, returning a tuple of 6 values (lon, lat, dist, vel_lon, vel_lat, vel_dist)
- `calc_uranian_planet(body_id, jd_tt)` — generic version for any Uranian planet
- `calc_transpluto(jd_tt)`, `calc_vulcan(jd_tt)`, `calc_waldemath(jd_tt)` — other hypothetical bodies
- `calc_white_moon_position(jd_tt)` — the White Moon (opposite of the mean Black Moon Lilith)
- `calc_hypothetical_position(body_id, jd_tt)` — generic function for any hypothetical body
- `load_bundled_fictitious_orbits()` — loads predefined fictitious orbits
- `parse_orbital_elements(filepath)` — loads fictitious orbits from a custom file
- `get_orbital_body_by_name(elements, name)` — searches for a body by name in the list of orbits
- `calc_orbital_position(elem, jd_tt)` — calculates the position of a body given its orbital elements
- `calc_all_arabic_parts(positions, jd=..., geo_lat=..., geo_lon=...)` — calculates the four main Arabic parts
