# Chapter 9 — Solar and Lunar Eclipses

## What you will learn

In this chapter, you will discover how eclipses occur, why there isn't one every month, what types exist, and how to use the library to find the next eclipse, calculate its details, and determine if it will be visible from your observation location.

---

## 9.1 How an eclipse occurs

An eclipse is an almost perfect alignment between the Sun, Moon, and Earth. There are two fundamental types:

- **Solar eclipse**: the Moon passes in front of the Sun, casting its shadow on the Earth. It always occurs during a **New Moon** (Sun and Moon on the same side of the sky).

- **Lunar eclipse**: the Earth blocks the sunlight illuminating the Moon, and the Moon is obscured. It always occurs during a **Full Moon** (Sun and Moon on opposite sides).

### Why doesn't it happen every month?

If the Moon's orbit were exactly on the same plane as the Earth's orbit (the ecliptic), we would have a solar eclipse at every New Moon and a lunar eclipse at every Full Moon — two per month. But the lunar orbit is **inclined by about 5.1°** with respect to the ecliptic. The Moon almost always passes "above" or "below" the Sun (in solar eclipses) or "above" or "below" the Earth's shadow (in lunar eclipses).

An eclipse can only occur when the New or Full Moon falls **near a node** — one of the two points where the lunar orbit crosses the ecliptic (see Chapter 6). This happens during "eclipse seasons", windows of about 35 days that repeat every ~173 days (about every 6 months).

In a typical year there are 2–5 solar eclipses and 0–3 lunar eclipses.

---

## 9.2 Types of eclipses

### Solar eclipses

The type of solar eclipse depends on the Moon's distance from the Earth at that moment (the Moon moves closer and farther away because its orbit is elliptical):

- **Total** (`ECL_TOTAL`): the Moon is close enough to *completely* cover the solar disk. For a few minutes, the sky turns dark and the solar corona becomes visible — one of nature's most extraordinary spectacles. The path of totality is narrow: typically 100–250 km wide.

- **Annular** (`ECL_ANNULAR`): the Moon is too far away and its disk appears smaller than the Sun's. A bright **ring of fire** ("annulus") remains around the lunar disk. Spectacular, but the sky does not turn dark like in a total eclipse.

- **Partial** (`ECL_PARTIAL`): the Moon covers only a part of the Sun. It is the most common type to observe — you just need to be in the lunar penumbra, which covers a much wider area than the umbra.

- **Hybrid** (`ECL_ANNULAR_TOTAL`): a rare eclipse that is annular in some areas of the Earth and total in others. It happens when the Moon's shadow is at its limit — the apex of the umbral cone grazes the Earth's surface.

### Lunar eclipses

- **Total** (`ECL_TOTAL`): the Moon completely enters the Earth's umbral cone. It does not become invisible but takes on a **copper red** hue — sunlight filtered and refracted by the Earth's atmosphere dimly illuminates it. Every "Blood Moon" is a total lunar eclipse.

- **Partial** (`ECL_PARTIAL`): only a part of the Moon enters the umbra. A dark "bite" is seen on the lunar disk.

- **Penumbral** (`ECL_PENUMBRAL`): the Moon enters the Earth's penumbra (the partial shadow zone), but not the true umbra. The obscuration is so slight that it is often invisible to the naked eye.

### Magnitude and obscuration

Two numbers describe "how much" of an eclipse it is:

- The **magnitude** is the fraction of the Sun's (or Moon's) *diameter* covered. A magnitude of 0.5 means half the diameter is covered. For total eclipses, the magnitude is ≥ 1.0.

- The **obscuration** is the fraction of the *area* covered. For the same magnitude, the obscuration is always greater (because the area grows with the square of the radius). A magnitude of 0.5 corresponds to an obscuration of about 0.39.

---

## 9.3 Finding the next eclipse

### Global solar eclipse

The `sol_eclipse_when_glob` function searches for the next solar eclipse **anywhere in the world**. It doesn't tell you if it will be visible from your location — only that it will happen somewhere on Earth.

```python
import libephemeris as ephem

# Search for the next solar eclipse starting from today
jd_start = ephem.julday(2024, 1, 1, 0.0)

ecl_type, tret = ephem.sol_eclipse_when_glob(jd_start)

# Decode the type
tipi = []
if ecl_type & ephem.ECL_TOTAL:
    tipi.append("total")
if ecl_type & ephem.ECL_ANNULAR:
    tipi.append("annular")
if ecl_type & ephem.ECL_PARTIAL:
    tipi.append("partial")
if ecl_type & ephem.ECL_ANNULAR_TOTAL:
    tipi.append("hybrid")

# tret[0] = time of maximum eclipse
jd_max = tret[0]
y, m, d, h = ephem.revjul(jd_max)
ore = int(h)
minuti = int((h - ore) * 60)

print(f"Next solar eclipse: {tipi}")
print(f"Date: {d}/{m}/{y} at {ore:02d}:{minuti:02d} UT")

# Contact times
if tret[1] > 0:
    y1, m1, d1, h1 = ephem.revjul(tret[1])
    print(f"Start (1st contact): {d1}/{m1}/{y1} {int(h1):02d}:{int((h1%1)*60):02d} UT")
if tret[4] > 0:
    y4, m4, d4, h4 = ephem.revjul(tret[4])
    print(f"End (4th contact):   {d4}/{m4}/{y4} {int(h4):02d}:{int((h4%1)*60):02d} UT")
```

```
Next solar eclipse: ['total']
Date: 8/4/2024 at 18:17 UT
End (4th contact):   8/4/2024 16:41 UT
```

You can filter by eclipse type:

```python
# Search only for total eclipses
ecl_type, tret = ephem.sol_eclipse_when_glob(
    jd_start,
    ifltype=ephem.ECL_TOTAL
)
```

### The next 5 solar eclipses

```python
import libephemeris as ephem

jd = ephem.julday(2024, 1, 1, 0.0)

print("Next 5 solar eclipses:\n")

for i in range(5):
    ecl_type, tret = ephem.sol_eclipse_when_glob(jd)
    jd_max = tret[0]

    # Type
    if ecl_type & ephem.ECL_TOTAL:
        tipo = "Total"
    elif ecl_type & ephem.ECL_ANNULAR:
        tipo = "Annular"
    elif ecl_type & ephem.ECL_ANNULAR_TOTAL:
        tipo = "Hybrid"
    else:
        tipo = "Partial"

    y, m, d, h = ephem.revjul(jd_max)
    print(f"  {i+1}. {d:2.0f}/{m:02.0f}/{y:.0f}  {tipo}")

    # Advance past this eclipse
    jd = jd_max + 30
```

```
Next 5 solar eclipses:

  1.  8/04/2024  Total
  2.  2/10/2024  Annular
  3. 29/03/2025  Partial
  4. 21/09/2025  Partial
  5. 17/02/2026  Annular
```

### Solar eclipse visible from a location

The most common question is: "When will the next eclipse be visible from my home?" To answer, use `sol_eclipse_when_loc`:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 1, 1, 0.0)
lat, lon = 41.9028, 12.4964  # Rome

# Search for the next solar eclipse visible from Rome
ecl_type, times, attr = ephem.sol_eclipse_when_loc(
    jd, (lon, lat, 0.0)
)

jd_max = times[0]
y, m, d, h = ephem.revjul(jd_max)
hour = int(h)
minute = int((h - hour) * 60)

magnitude = attr[0]
obscuration = attr[2]
sun_alt = attr[4]

print(f"Next eclipse visible from Rome:")
print(f"Date: {d:.0f}/{m:.0f}/{y:.0f} at {hour:02d}:{minute:02d} UT")
print(f"Magnitude: {magnitude:.3f}")
print(f"Obscuration: {obscuration:.1%}")
print(f"Sun altitude at maximum: {sun_alt:.1f}°")
```

```
Next eclipse visible from Rome:
Date: 29/3/2025 at 11:03 UT
Magnitude: 0.073
Obscuration: 2.4%
Sun altitude at maximum: 51.6°
```

The difference between the two functions is crucial:

- `sol_eclipse_when_glob` finds the next eclipse *in the world* — there are 2–5 per year
- `sol_eclipse_when_loc` finds the next one visible *from your observation point* — years can pass between one and the next

---

## 9.4 Eclipse details

### Where is it visible at maximum?

Given the eclipse, `sol_eclipse_where` tells you where in the world the point of maximum falls:

```python
import libephemeris as ephem

# Total eclipse of April 8, 2024
jd_start = ephem.julday(2024, 4, 1, 0.0)
ecl_type, tret = ephem.sol_eclipse_when_glob(jd_start)
jd_max = tret[0]

# Where does the maximum fall?
ret, geopos, attr = ephem.sol_eclipse_where(jd_max)

lon_centro = geopos[0]
lat_centro = geopos[1]
larghezza = attr[3]  # path width in km

print(f"Center of the central line: {lat_centro:.2f}°N, {lon_centro:.2f}°E")
print(f"Width of the path of totality: {larghezza:.0f} km")
print(f"Magnitude at the center: {attr[0]:.4f}")
```

```
Center of the central line: 25.29°N, -104.15°E
Width of the path of totality: 207 km
Magnitude at the center: 1.0571
```

### Type and magnitude at a given time and location

If you know when an eclipse occurs and want to know how it appears from a specific location, use `sol_eclipse_how`:

```python
import libephemeris as ephem

# Eclipse of April 8, 2024 viewed from Dallas, Texas
jd_max = ephem.julday(2024, 4, 8, 18.0 + 42/60)  # ~18:42 UT
lat, lon = 32.78, -96.80  # Dallas

ret, attr = ephem.sol_eclipse_how(jd_max, (lon, lat, 0.0))

if ret > 0:
    magnitudine = attr[0]
    oscuramento = attr[2]
    alt_sole = attr[5]  # true altitude of the Sun
    az_sole = attr[4]   # azimuth of the Sun

    print(f"Magnitude: {magnitudine:.4f}")
    print(f"Obscuration: {oscuramento:.1%}")
    print(f"Sun: altitude {alt_sole:.1f}°, azimuth {az_sole:.1f}°")
else:
    print("No eclipse visible from this location at this time")
```

```
Magnitude: 1.0126
Obscuration: 111.6%
Sun: altitude 64.6°, azimuth 7.6°
```

### Complete details with `sol_eclipse_how_details`

For a full report with all contacts, position angles, and durations, use the extended version:

```python
import libephemeris as ephem

jd_start = ephem.julday(2024, 4, 1, 0.0)
ecl_type, tret = ephem.sol_eclipse_when_glob(jd_start)

# Complete details from Rome
info = ephem.sol_eclipse_how_details(
    tret[0], 41.9028, 12.4964
)

if info["is_visible"]:
    print(f"Type: {'total' if info['is_total'] else 'partial'}")
    print(f"Maximum magnitude: {info['max_magnitude']:.4f}")
    print(f"Maximum obscuration: {info['max_obscuration_percent']:.1f}%")
    print(f"Partial phase duration: {info['duration_partial_minutes']:.1f} min")

    if info["is_total"]:
        print(f"Totality duration: {info['duration_total_minutes']:.1f} min")
```

```
Not visible from Rome
```

---

## 9.5 Lunar eclipses

Lunar eclipses are found with similar functions.

### Finding the next lunar eclipse

```python
import libephemeris as ephem

jd = ephem.julday(2024, 1, 1, 0.0)

ecl_type, tret = ephem.lun_eclipse_when(jd)

# Type
if ecl_type & ephem.ECL_TOTAL:
    tipo = "Total"
elif ecl_type & ephem.ECL_PARTIAL:
    tipo = "Partial"
elif ecl_type & ephem.ECL_PENUMBRAL:
    tipo = "Penumbral"
else:
    tipo = "Unknown"

jd_max = tret[0]
y, m, d, h = ephem.revjul(jd_max)
print(f"Next lunar eclipse: {tipo}")
print(f"Date: {d:.0f}/{m:.0f}/{y:.0f}")
```

```
Next lunar eclipse: Penumbral
Date: 25/3/2024
```

The `tret` tuple for lunar eclipses contains 8 times:

- `tret[0]` — time of maximum eclipse
- `tret[1]` — start of partial phase (Moon enters umbra)
- `tret[2]` — start of totality (if total, otherwise 0)
- `tret[3]` — end of totality
- `tret[4]` — end of partial phase (Moon leaves umbra)
- `tret[5]` — start of penumbral phase
- `tret[6]` — end of penumbral phase

### Lunar eclipse visible from a location

A lunar eclipse is visible wherever the Moon is above the horizon during the event. But it might be partially visible — perhaps the Moon rises when the eclipse is already halfway through, or sets before it ends.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 1, 1, 0.0)
lat, lon = 41.9028, 12.4964  # Rome

ecl_type, tempi, attr = ephem.lun_eclipse_when_loc(jd, (lon, lat, 0.0))

jd_max = tempi[0]
y, m, d, h = ephem.revjul(jd_max)

# Umbral magnitude tells "how deep" the eclipse is
mag_umbrale = attr[0]

# Check visibility
visibile = ecl_type & ephem.ECL_VISIBLE

print(f"Lunar eclipse: {d:.0f}/{m:.0f}/{y:.0f}")
print(f"Umbral magnitude: {mag_umbrale:.3f}")

if visibile:
    alt_luna = attr[5]
    az_luna = attr[4]
    print(f"Visible from Rome!")
    print(f"Moon at maximum: altitude {alt_luna:.1f}°, azimuth {az_luna:.1f}°")
else:
    print("Not visible from Rome (Moon below horizon)")
```

```
Lunar eclipse: 25/3/2024
Umbral magnitude: 0.000
Visible from Rome!
Moon at maximum: altitude -22.7°, azimuth 109.2°
```

### Type and magnitude at a specific time

```python
import libephemeris as ephem

# Lunar eclipse: let's check a specific moment
jd = ephem.julday(2025, 3, 14, 6.0)  # example

ret, attr = ephem.lun_eclipse_how(jd, (12.4964, 41.9028, 0.0))

if ret > 0:
    mag_umbrale = attr[0]
    mag_penombrale = attr[1]
    saros = int(attr[9])
    saros_membro = int(attr[10])

    print(f"Umbral magnitude: {mag_umbrale:.3f}")
    print(f"Penumbral magnitude: {mag_penombrale:.3f}")
    if saros > 0:
        print(f"Saros series: {saros}, member {saros_membro}")
```

```
Umbral magnitude: 1.090
Penumbral magnitude: 2.000
Saros series: 123, member 53
```

---

## 9.6 Eclipse cycles: Saros and Inex

Eclipses do not happen randomly — they follow precise cycles known since antiquity.

### The Saros cycle

The **Saros** is a period of 6585.32 days, or roughly **18 years, 11 days, and 8 hours**. After one Saros, the Sun, Moon, and nodes return almost exactly to the same configuration, and a very similar eclipse repeats.

But there's that "almost": the 8 hours of difference mean the Earth has rotated about a third of a turn extra. So the next eclipse falls about **120° further west** on the Earth's surface.

Each "Saros series" is a family of related eclipses that is born as a partial eclipse at one pole, gradually becomes total, then returns to partial at the other pole, lasting a total of about 1200–1400 years and 70–80 eclipses.

```python
import libephemeris as ephem

# Total eclipse of April 8, 2024
jd_start = ephem.julday(2024, 4, 1, 0.0)
ecl_type, tret = ephem.sol_eclipse_when_glob(
    jd_start, ecltype=ephem.ECL_TOTAL
)
jd_ecl = tret[0]

# Which Saros series does it belong to?
saros = ephem.get_saros_number(jd_ecl, "solar")
print(f"Saros series: {saros}")

# The "sister" 18 years prior
jd_sorella = jd_ecl - ephem.SAROS_CYCLE_DAYS
y, m, d, h = ephem.revjul(jd_sorella)
print(f"Sister eclipse: ~{d:.0f}/{m:.0f}/{y:.0f}")
```

```
Saros series: 139
Sister eclipse: ~29/3/2006
```

### The Inex cycle

The **Inex** is a longer period of 10571.95 days, about **29 years**. It connects eclipses of different Saros series — when a series ends, the Inex connects it to the next series. The two cycles together form a grid covering all past and future eclipses.

```python
import libephemeris as ephem

jd_ecl = ephem.julday(2024, 4, 8, 18.0)

inex = ephem.get_inex_number(jd_ecl, "solar")
print(f"Inex number: {inex}")
```

```
Inex number: 50
```

---

## 9.7 Occultations

An **occultation** occurs when a celestial body passes in front of another, hiding it. The most common case is the Moon occulting a star or a planet.

### Occultation of a star

```python
import libephemeris as ephem

# Search for the next occultation of Regulus by the Moon
jd = ephem.julday(2024, 1, 1, 0.0)

ret, tret = ephem.lun_occult_when_glob(
    jd,
    body="Regulus",
)

if ret > 0:
    jd_max = tret[0]
    y, m, d, h = ephem.revjul(jd_max)
    print(f"Next occultation of Regulus: {d:.0f}/{m:.0f}/{y:.0f}")
```

```
Next occultation of Regulus: 16/10/2025
```

### Occultation of a planet

```python
import libephemeris as ephem

# Search for the next occultation of Jupiter by the Moon
jd = ephem.julday(2024, 1, 1, 0.0)

ret, tret = ephem.lun_occult_when_glob(
    jd,
    body=ephem.JUPITER,
)

if ret > 0:
    jd_max = tret[0]
    y, m, d, h = ephem.revjul(jd_max)
    print(f"Next occultation of Jupiter: {d:.0f}/{m:.0f}/{y:.0f}")
```

```
Next occultation of Jupiter: 8/9/2026
```

### Occultation visible from a location

```python
import libephemeris as ephem

jd = ephem.julday(2024, 1, 1, 0.0)

# Next lunar occultation of Spica visible from Rome
ret, tempi, attr = ephem.lun_occult_when_loc(
    jd,
    "Spica",
    geopos=(12.4964, 41.9028, 0.0)
)

if ret > 0:
    jd_inizio = tempi[1]  # first contact
    jd_fine = tempi[4]    # last contact
    y, m, d, h = ephem.revjul(tempi[0])

    durata = (jd_fine - jd_inizio) * 24 * 60  # in minutes

    print(f"Occultation of Spica: {d:.0f}/{m:.0f}/{y:.0f}")
    print(f"Duration: {durata:.0f} minutes")
    print(f"Altitude at maximum: {attr[5]:.1f}°")
else:
    print("No occultation found in search period")
```

```
Occultation of Spica: 1/3/2032
Duration: 34 minutes
Altitude at maximum: 28.3°
```

### Occultations between planets

Extremely rare but calculable: a planet occulting another planet or a star.

```python
import libephemeris as ephem

# Search for an occultation of Regulus by Venus
jd = ephem.julday(2024, 1, 1, 0.0)

ret, tret = ephem.planet_occult_when_glob(
    jd,
    occulting_planet=ephem.VENUS,
    occulted_planet=0,
    starname="Regulus"
)

if ret > 0:
    y, m, d, h = ephem.revjul(tret[0])
    print(f"Venus occults Regulus: {d:.0f}/{m:.0f}/{y:.0f}")
else:
    print("No occultation found (can search up to 150 years)")
```

> **Note**: This search can take a long time (minutes) because planetary occultations are extremely rare. The function searches up to 150 years into the future.

---

## Summary

In this chapter, we explored eclipses, from geometry to practice.

**Key concepts:**

- Eclipses occur only when a New Moon (solar eclipse) or Full Moon (lunar eclipse) falls near a lunar node — about 4–7 times a year in total
- Solar eclipses can be **total**, **annular**, **partial**, or **hybrid**, depending on the Moon's distance
- Lunar eclipses can be **total** (Blood Moon), **partial**, or **penumbral** (almost invisible)
- **Magnitude** measures the fraction of the diameter covered; **obscuration** measures the fraction of the area covered
- The **Saros** cycle (~18 years) connects similar eclipses; the **Inex** cycle (~29 years) connects different series
- **Occultations** are eclipses of stars or planets by the Moon or other planets

**Functions introduced:**

- `sol_eclipse_when_glob(jd, ifltype=0)` — finds the next solar eclipse in the world, returns type and contact times
- `sol_eclipse_when_loc(jd, geopos)` — finds the next solar eclipse visible from a location, with magnitude and obscuration
- `sol_eclipse_where(jd)` — given the time of maximum, returns the coordinates of the central line
- `sol_eclipse_how(jd, geopos)` — calculates type, magnitude, and Sun position for an already known eclipse
- `sol_eclipse_how_details(jd, lat, lon)` — extended version with all contacts, position angles, and durations
- `lun_eclipse_when(jd, ifltype=0)` — finds the next lunar eclipse (global)
- `lun_eclipse_when_loc(jd, geopos)` — finds the next lunar eclipse visible from a location
- `lun_eclipse_how(jd, geopos)` — details of a lunar eclipse at a time and location
- `lun_occult_when_glob(jd, body)` — finds the next lunar occultation of a star or planet
- `lun_occult_when_loc(jd, body, geopos)` — occultation visible from a location
- `planet_occult_when_glob(jd, occulting_planet, starname=...)` — occultation between planets
- `get_saros_number(jd, "solar"|"lunar")` — returns the Saros series of an eclipse
- `get_inex_number(jd, "solar"|"lunar")` — returns the Inex number
