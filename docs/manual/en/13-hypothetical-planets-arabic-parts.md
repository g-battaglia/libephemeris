# Chapter 13 — Hypothetical planets and Arabic parts

## What you will learn

In this chapter you will discover what Uranian planets are (mathematical
points used in the Hamburg School), how built-in provenance is classified,
how to define custom fictitious orbits, and how
to calculate Arabic parts — ancient formulas that combine planetary positions
and angles.

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

### Availability

IDs 40–47 are built-in Keplerian calculations. The orbital elements are
transcribed from the Witte/Sieggruen publications.

---

## 13.2 Built-in hypothetical bodies

All historical IDs 40–58 compute. **Harrington** (`HARRINGTON`, ID 50) is
directly derived from Harrington's 1988 *Astronomical Journal* paper. All
models use built-in orbital elements; see the methodology page.

```python
import libephemeris as ephem

jd_tt = ephem.julday(2024, 4, 8, 12.0) + ephem.deltat(ephem.julday(2024, 4, 8, 12.0))
pos = ephem.calc_hypothetical_position(ephem.HARRINGTON, jd_tt)
```

No reference-distribution data file is bundled. The fixed-element dataset
contains a source/status field for every row.

---

## 13.3 Custom fictitious orbits

If you need a hypothetical body not included in the library, you can define
your own orbit using LibEphemeris's documented nine-field orbital-elements text
format. Supply values from sources you are entitled to use; reference-distribution
data files are not required or supported as bundled inputs.

### Loading the reviewed built-in orbit

The library includes a file of predefined fictitious orbits:

```python
import libephemeris as ephem

# Load predefined orbits
orbits = ephem.load_bundled_fictitious_orbits()
print(f"Loaded {len(orbits)} fictitious orbits")

# Search for the reviewed body by name
body = ephem.get_orbital_body_by_name(orbits, "Harrington")
if body:
    print(f"Found: {body.name}")
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
- All historical hypothetical IDs **40–58** calculate again
- Harrington has an independent primary derivation; all models use built-in
  orbital elements
- **Custom fictitious orbits** allow you to define any hypothetical body with its own orbital elements
- The **Arabic parts** are points calculated by the formula ASC + Planet A − Planet B, with the Part of Fortune being the most important

**Functions introduced:**

- `calc_hypothetical_position(body_id, jd_tt)` — calculate any historical
  hypothetical body ID 40–58
- `load_bundled_fictitious_orbits()` — loads predefined fictitious orbits
- `parse_orbital_elements(filepath)` — loads fictitious orbits from a custom file
- `get_orbital_body_by_name(elements, name)` — searches for a body by name in the list of orbits
- `calc_orbital_position(elem, jd_tt)` — calculates the position of a body given its orbital elements
- `calc_all_arabic_parts(positions, jd=..., geo_lat=..., geo_lon=...)` — calculates the four main Arabic parts
