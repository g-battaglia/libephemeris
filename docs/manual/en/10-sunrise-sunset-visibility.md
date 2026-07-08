# Chapter 10 — Sunrise, sunset and visibility

## What you will learn

In this chapter, you will discover how to calculate sunrise, sunset, and meridian transits of any celestial body, what the different types of twilight are, how atmospheric refraction works (which makes objects near the horizon "lift"), and how to determine whether a star or planet is visible to the naked eye at a given time.

---

## 10.1 Sunrise and sunset

### What does "the Sun rises" mean?

The definition seems obvious, but it is actually surprisingly complicated. The moment of sunrise is not when the *center* of the Sun touches the horizon, but when the **upper edge** (the "limb") of the solar disk appears on the horizon. Since the Sun's apparent diameter is about 32 arcminutes (0.53°), the center is still about 16' below the geometric horizon when we see the first ray.

But there's more: the Earth's atmosphere **bends light** (refraction), making objects near the horizon appear higher than they really are. At the horizon, refraction lifts the Sun by about **34 arcminutes**. Therefore, at the moment we "see" the Sun rise, its center is actually about 50' (almost a degree!) below the geometric horizon.

The library takes all of this into account automatically.

### The `rise_trans` function

```python
import libephemeris as ephem

# Sunrise and sunset of the Sun today in Milan
jd = ephem.julday(2024, 6, 21, 0.0)  # Summer solstice
lat, lon = 45.4642, 9.1900  # Milan

# Sunrise
ret, tret = ephem.rise_trans(
    jd, ephem.SUN,
    ephem.CALC_RISE,
    (lon, lat, 0.0),
    atpress=1013.25,
    attemp=15.0
)
jd_sunrise = tret[0]  # the event Julian Day

# Sunset
ret, tret = ephem.rise_trans(
    jd, ephem.SUN,
    ephem.CALC_SET,
    (lon, lat, 0.0)
)
jd_sunset = tret[0]

# Convert to local time (CEST = UT + 2)
tz_offset = 2  # summer daylight saving time

_, _, _, h_sunrise = ephem.revjul(jd_sunrise)
_, _, _, h_sunset = ephem.revjul(jd_sunset)

hours_sunrise = h_sunrise + tz_offset
hours_sunset = h_sunset + tz_offset
duration = (jd_sunset - jd_sunrise) * 24

print(f"Sunrise:           {int(hours_sunrise):02d}:{int((hours_sunrise % 1) * 60):02d} local time")
print(f"Sunset:            {int(hours_sunset):02d}:{int((hours_sunset % 1) * 60):02d} local time")
print(f"Daylight duration: {duration:.1f} hours")
```

```
Sunrise:           05:34 local time
Sunset:            21:15 local time
Daylight duration: 15.7 hours
```

The `rsmi` parameter accepts four base events:

- `CALC_RISE` (1) — rise (first ray at the horizon)
- `CALC_SET` (2) — set (last ray disappears)
- `CALC_MTRANSIT` (4) — upper meridian transit (highest point in the sky)
- `CALC_ITRANSIT` (8) — lower transit (lowest point, usually below the horizon)

You can add modifiers to these (using the `|` operator):

- `BIT_DISC_CENTER` (256) — use the center of the disc instead of the upper edge
- `BIT_DISC_BOTTOM` (8192) — use the lower edge of the disc
- `BIT_NO_REFRACTION` (512) — do not apply atmospheric refraction

### Moonrise and moonset

It works exactly like the Sun:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 0.0)
lat, lon = 41.9028, 12.4964  # Rome

ret, tret = ephem.rise_trans(
    jd, ephem.MOON, ephem.CALC_RISE,
    (lon, lat, 0.0)
)
jd_rise = tret[0]

ret, tret = ephem.rise_trans(
    jd, ephem.MOON, ephem.CALC_SET,
    (lon, lat, 0.0)
)
jd_set = tret[0]

if jd_rise > 0 and jd_set > 0:
    _, _, _, h_rise = ephem.revjul(jd_rise)
    _, _, _, h_set = ephem.revjul(jd_set)
    print(f"Moonrise: {int(h_rise):02d}:{int((h_rise % 1) * 60):02d} UT")
    print(f"Moonset:  {int(h_set):02d}:{int((h_set % 1) * 60):02d} UT")
```

```
Moonrise: 04:29 UT
Moonset:  17:36 UT
```

The Moon has a peculiarity: its rising is delayed by about 50 minutes each day. Therefore, sometimes it doesn't rise at all on a given day, or it rises but doesn't set (or vice versa).

### Circumpolar bodies

At high latitudes, some bodies never rise or never set. For example, at the North Pole, the Sun remains above the horizon for six months. When this happens, `rise_trans` returns a `retflag` of `-2`:

```python
import libephemeris as ephem

# Sun in Tromsø on June 21st — it never sets
jd = ephem.julday(2024, 6, 21, 0.0)
ret, tret = ephem.rise_trans(
    jd, ephem.SUN, ephem.CALC_SET,
    (19.0, 69.6, 0.0)
)

if ret == -2:
    print("The Sun never sets! (circumpolar)")
```

```
The Sun never sets! (circumpolar)
```

---

## 10.2 Meridian transits

The **upper transit** (or culmination) is the moment when a celestial body reaches its highest point in the sky, crossing the local meridian. For the Sun, this is **solar noon** — which almost never coincides with 12:00 PM on the clock.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 0.0)
lat, lon = 41.9028, 12.4964  # Rome

# Sun's meridian transit = solar noon
ret, tret = ephem.rise_trans(
    jd, ephem.SUN, ephem.CALC_MTRANSIT,
    (lon, lat, 0.0)
)
jd_transit = tret[0]

_, _, _, h = ephem.revjul(jd_transit)
print(f"Solar noon in Rome:  {int(h):02d}:{int((h % 1) * 60):02d} UT")
print(f"In local time (CET): {int(h+1):02d}:{int((h % 1) * 60):02d}")
```

```
Solar noon in Rome:  11:11 UT
In local time (CET): 12:11
```

Solar noon varies throughout the year due to the **equation of time** (see Chapter 2): it can be up to 16 minutes before or 14 minutes after 12:00 at the local meridian. And if your time zone doesn't exactly match your longitude, the difference is even greater.

---

## 10.3 Twilights

The sky doesn't become dark immediately after sunset. The transition from daylight to night happens gradually, and astronomers distinguish three phases of twilight, defined by the Sun's position below the horizon:

**Civil twilight** — the Sun is between 0° and −6° below the horizon. There is still enough light for most outdoor activities. You can read a book and distinguish the outlines of buildings. The brightest planets (Venus, Jupiter) begin to appear.

**Nautical twilight** — the Sun is between −6° and −12°. The sea horizon is still visible — hence the name: sailors could still use a sextant for navigation. The brightest stars are visible.

**Astronomical twilight** — the Sun is between −12° and −18°. The sky appears dark to the naked eye, but there is still a slight residual glow on the horizon. Faint stars are not yet visible.

**Astronomical night** — the Sun is below −18°. The sky is completely dark (aside from light pollution). All stars visible to the naked eye can be seen.

To calculate twilights, use `rise_trans` with the twilight flags:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 6, 21, 0.0)  # Summer solstice
lat, lon = 45.4642, 9.1900  # Milan
tz_offset = 2  # CEST

def local_time(jd_event, tz_offset):
    _, _, _, h = ephem.revjul(jd_event)
    h += tz_offset
    return f"{int(h):02d}:{int((h % 1) * 60):02d}"

# Sunset
_, tret = ephem.rise_trans(jd, ephem.SUN,
    ephem.CALC_SET, (lon, lat, 0.0))
jd_s = tret[0]

# End of civil twilight (Sun at -6°)
_, tret = ephem.rise_trans(jd, ephem.SUN,
    ephem.CALC_SET | ephem.BIT_CIVIL_TWILIGHT, (lon, lat, 0.0))
jd_civ = tret[0]

# End of nautical twilight (Sun at -12°)
_, tret = ephem.rise_trans(jd, ephem.SUN,
    ephem.CALC_SET | ephem.BIT_NAUTIC_TWILIGHT, (lon, lat, 0.0))
jd_nau = tret[0]

# End of astronomical twilight (Sun at -18°)
ret, tret = ephem.rise_trans(jd, ephem.SUN,
    ephem.CALC_SET | ephem.BIT_ASTRO_TWILIGHT, (lon, lat, 0.0))
jd_ast = tret[0]

print(f"Sunset:                       {local_time(jd_s, tz_offset)}")
print(f"End of civil twilight:        {local_time(jd_civ, tz_offset)}")
print(f"End of nautical twilight:     {local_time(jd_nau, tz_offset)}")

if ret != -2:
    print(f"End of astronomical twilight: {local_time(jd_ast, tz_offset)}")
else:
    print("Astronomical night: NEVER (the sky never gets completely dark)")
```

```
Sunset:                       21:15
End of civil twilight:        21:55
End of nautical twilight:     22:46
End of astronomical twilight: 23:57
```

At the summer solstice in Milan (45° N), astronomical twilight never completely ends — the Sun never drops below −18°. The white nights of northern latitudes are a consequence of this phenomenon.

---

## 10.4 Atmospheric refraction

**Refraction** is the bending of light caused by the Earth's atmosphere. The atmosphere is denser at the bottom and more rarefied at the top, and light follows a curve instead of a straight line. The effect is greatest at the horizon and decreases as you look higher up:

- At the horizon (0°): **~34'** of refraction — almost the diameter of the Sun!
- At 10° of altitude: ~5'
- At 45°: ~1'
- At the zenith (90°): 0' (light arrives vertically, no bending)

The `refrac` function converts between "true" (geometric) altitude and "apparent" altitude (what you see):

```python
import libephemeris as ephem

# From true to apparent altitude
true_alt = 0.0  # object exactly at the geometric horizon
apparent_alt = ephem.refrac(true_alt, calc_flag=ephem.TRUE_TO_APP)
print(f"True altitude: {true_alt:.2f}° → apparent: {apparent_alt:.2f}°")
# → about 0.57° (34 arcminutes above the horizon)

# From apparent to true
app_alt = 5.0
true_alt2 = ephem.refrac(app_alt, calc_flag=ephem.APP_TO_TRUE)
print(f"Apparent altitude: {app_alt:.2f}° → true: {true_alt2:.2f}°")
```

```
True altitude: 0.00° → apparent: 0.48°
Apparent altitude: 5.00° → true: 4.84°
```

Refraction depends on atmospheric conditions — pressure and temperature:

```python
import libephemeris as ephem

# Refraction at different temperatures (at the horizon)
for temp in [-20, 0, 15, 35]:
    r = ephem.refrac(0.0, pressure=1013.25, temperature=temp)
    print(f"T = {temp:+3d}°C → refraction at the horizon: {r:.2f}°")
```

```
T = -20°C → refraction at the horizon: 0.54°
T =  +0°C → refraction at the horizon: 0.50°
T = +15°C → refraction at the horizon: 0.48°
T = +35°C → refraction at the horizon: 0.45°
```

### Extended refraction: observers at high altitudes

If you observe from a mountain or an airplane, the horizon is lower than normal (the "dip of the horizon" or simply "dip"). The `refrac_extended` function takes this into account:

```python
import libephemeris as ephem

# Observation from the top of Mont Blanc (4810 m)
alt_obj = 0.5  # object half a degree above the geometric horizon

alt_result, details = ephem.refrac_extended(
    alt_obj,
    altitude_geo=4810.0,  # observer's altitude in meters
    pressure=550.0,       # reduced pressure at high altitude
    temperature=-10.0     # cold!
)

true_alt, app_alt, refraction, dip = details

print(f"True altitude:       {true_alt:.3f}°")
print(f"Apparent altitude:   {app_alt:.3f}°")
print(f"Refraction:          {refraction:.3f}°")
print(f"Dip of the horizon:  {dip:.3f}°")
```

```
True altitude:       0.500°
Apparent altitude:   0.744°
Refraction:          0.244°
Dip of the horizon:  -1.926°
```

### Custom horizon

If your horizon is not flat (mountains, buildings), you can use `rise_trans_true_hor` to specify a custom horizon altitude:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 0.0)
lat, lon = 45.4642, 9.1900  # Milan

# The Alps to the north raise the horizon by about 2 degrees
_, tret = ephem.rise_trans_true_hor(
    jd, ephem.SUN, ephem.CALC_RISE,
    (lon, lat, 0.0),
    horhgt=2.0,  # horizon raised 2° above the plane
)
jd_sunrise = tret[0]

_, _, _, h = ephem.revjul(jd_sunrise)
print(f"Sunrise with 2° horizon: {int(h):02d}:{int((h % 1) * 60):02d} UT")
```

```
Sunrise with 2° horizon: 05:01 UT
```

---

## 10.5 Heliacal visibility

**Heliacal rising** is one of the oldest concepts in astronomy. It is the first day of the year when a star (or a planet) becomes visible at dawn, after being hidden by the Sun's glare for weeks or months.

For the ancient Egyptians, the heliacal rising of **Sirius** (the brightest star) marked the beginning of the year and announced the impending flooding of the Nile — an event of vital importance. For the Babylonians, the heliacal rising of Venus was a fundamental omen.

Heliacal visibility depends on many factors: the brightness of the object, its distance from the Sun, the brightness of the sky at twilight, atmospheric conditions, and even the visual acuity of the observer.

### Finding the heliacal rising

The `heliacal_ut` function searches for the date of a heliacal event:

```python
import libephemeris as ephem

# When does Sirius become visible at dawn in 2024?
jd = ephem.julday(2024, 1, 1, 0.0)
lat, lon = 30.0, 31.2  # Cairo (where the Egyptians observed it)

jd_event, *_ = ephem.heliacal_ut(
    jd,
    geopos=(31.2, 30.0, 0.0),      # Cairo (lon, lat, alt)
    datm=(1013.25, 25.0, 30.0, 0.0),  # pressure, temp, humidity%, lapse
    dobs=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    object_name="Sirius",
    event_type=ephem.HELIACAL_RISING
)

if jd_event > 0:
    y, m, d, h = ephem.revjul(jd_event)
    print(f"Heliacal rising of Sirius in Cairo: {d:.0f}/{m:.0f}/{y:.0f}")
```

```
Heliacal rising of Sirius in Cairo: 4/8/2024
```

The four types of heliacal events are:

- `HELIACAL_RISING` (1) — **heliacal rising**: the first morning appearance. The object becomes visible at dawn after the period of invisibility near the Sun. Valid for planets and fixed stars (not the Sun or Moon).

- `HELIACAL_SETTING` (2) — **heliacal setting**: the last evening appearance. The object disappears into the Sun's glare at sunset. Valid for planets and fixed stars (not the Sun or Moon).

- `EVENING_FIRST` (3) — **evening first appearance**: the first day the body becomes visible in the evening after superior conjunction. Valid for the inner planets (Mercury, Venus) and the Moon.

- `MORNING_LAST` (4) — **morning last appearance**: the last day of visibility in the morning before inferior conjunction. Valid for the inner planets (Mercury, Venus) and the Moon.

### The reference-ephemeris-compatible API

For complete control over atmospheric conditions and the capabilities of the observer, use `heliacal_ut`:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 1, 1, 0.0)

geopos = (31.2, 30.0, 0.0)  # (lon, lat, alt) — attention: longitude FIRST

datm = (1013.25, 25.0, 30.0, 0.0)
#  pressure, temperature, relative humidity%, meteorological visibility

dobs = (36.0, 1.0, 0, 0, 0, 0)
#  age, visual acuity (Snellen), binocular, magnification, aperture, transmission

result = ephem.heliacal_ut(
    jd, geopos, datm, dobs,
    "Sirius",
    ephem.HELIACAL_RISING
)

jd_start = result[0]
if jd_start > 0:
    y, m, d, h = ephem.revjul(jd_start)
    print(f"Heliacal rising: {d:.0f}/{m:.0f}/{y:.0f}")
```

```
Heliacal rising: 4/8/2024
```

The `datm` tuple (atmospheric conditions):
- `[0]` — atmospheric pressure in mbar
- `[1]` — temperature in °C
- `[2]` — relative humidity as a percentage (0–100)
- `[3]` — meteorological visibility in km (0 = auto-calculate)

The `dobs` tuple (observer):
- `[0]` — observer's age in years (influences visual acuity)
- `[1]` — Snellen ratio (1.0 = normal vision, 1.5 = excellent vision)
- `[2]` — 0 = monocular, 1 = binocular (only with telescope)
- `[3]` — telescope magnification (0 = naked eye)
- `[4]` — telescope aperture in mm
- `[5]` — optical transmission coefficient

---

## 10.6 Atmospheric extinction and sky brightness

### How much is the light attenuated?

The atmosphere not only bends light (refraction), but also **absorbs** a part of it — especially near the horizon, where light travels through a thicker layer of air. This effect is called **atmospheric extinction**.

**Airmass** measures how much atmosphere the light has to travel through: at the zenith it is 1, at the horizon it's about 38:

```python
import libephemeris as ephem

# Airmass at different altitudes
for alt in [90, 60, 30, 10, 5, 1]:
    am = ephem.calc_airmass(float(alt))
    print(f"Altitude {alt:2d}° → airmass = {am:.2f}")
```

```
Altitude 90° → airmass = 1.00
Altitude 60° → airmass = 1.15
Altitude 30° → airmass = 1.99
Altitude 10° → airmass = 5.59
Altitude  5° → airmass = 10.31
Altitude  1° → airmass = 26.31
```

Total extinction depends on airmass and atmospheric conditions:

```python
import libephemeris as ephem

# How many magnitudes of light are lost at 10° altitude?
ext = ephem.calc_extinction_magnitude(10.0)
print(f"Extinction at 10° altitude: {ext:.2f} magnitudes")

# A star of magnitude 4.0 becomes:
app_mag = ephem.apparent_magnitude_with_extinction(4.0, 10.0)
print(f"Magnitude 4.0 → {app_mag:.2f} after extinction")
```

```
Extinction at 10° altitude: 1.58 magnitudes
Magnitude 4.0 → 5.58 after extinction
```

### The sky during twilight

The brightness of the sky during twilight determines which objects are visible. A bright sky "hides" faint stars:

```python
import libephemeris as ephem

# Sky brightness at twilight
for sun_alt in [0, -3, -6, -9, -12, -15, -18]:
    phase = ephem.get_twilight_phase(float(sun_alt))
    lim_mag = ephem.calc_limiting_magnitude_twilight(float(sun_alt))
    print(f"Sun at {sun_alt:+3d}° → {phase:14s}, "
          f"limiting magnitude: {lim_mag:.1f}")
```

```
Sun at  +0° → civil         , limiting magnitude: -0.7
Sun at  -3° → civil         , limiting magnitude: 0.6
Sun at  -6° → nautical      , limiting magnitude: 1.8
Sun at  -9° → nautical      , limiting magnitude: 3.7
Sun at -12° → astronomical  , limiting magnitude: 5.3
Sun at -15° → astronomical  , limiting magnitude: 5.9
Sun at -18° → night         , limiting magnitude: 6.4
```

The **limiting magnitude** is the magnitude of the faintest visible star. In the middle of the night, under perfect conditions, it's about 6.0–6.5. During civil twilight, it drops to 1–2 (only the brightest stars).

### Visibility of a planet

To determine whether an object is visible to the naked eye under specific conditions, use `vis_limit_mag`:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 12, 20, 17.0)  # evening, 18:00 CET

geopos = (12.4964, 41.9028, 0.0)  # Rome (lon, lat, alt)
datm = (1013.25, 10.0, 50.0, 0.0)
dobs = (36.0, 1.0, 0, 0, 0, 0)

result, dret = ephem.vis_limit_mag(
    jd, geopos, datm, dobs, "Venus"
)

limit_mag = dret[0]  # limiting magnitude
obj_mag = dret[7]  # apparent magnitude of the object
obj_alt = dret[1]  # altitude of the object

if result >= 0:  # not below the horizon
    visible = limit_mag > obj_mag
    print(f"Limiting magnitude: {limit_mag:.1f}")
    print(f"Magnitude of Venus: {obj_mag:.1f}")
    print(f"Altitude:           {obj_alt:.1f}°")
    print(f"Visible:            {'YES' if visible else 'NO'}")
else:
    print("Venus is below the horizon")
```

```
Limiting magnitude: 5.9
Magnitude of Venus: -3.5
Altitude:           20.0°
Visible:            YES
```

---

## Summary

In this chapter, we learned how to calculate when celestial bodies are visible in the sky.

**Key concepts:**

- **Sunrise** occurs when the upper edge of the solar disk appears on the horizon, taking atmospheric refraction into account (~34' of "lifting")
- **Twilights** are three progressive phases: civil (Sun at −6°), nautical (−12°), astronomical (−18°). Only below −18° is the sky truly dark
- Atmospheric **refraction** makes objects near the horizon appear higher than they actually are; the effect is greatest at the horizon (~34')
- Atmospheric **extinction** attenuates the light of objects low on the horizon; the airmass at the zenith is 1, at the horizon ~38
- **Heliacal rising** is the first appearance of a body at dawn after its period of invisibility near the Sun — a concept as old as astronomy itself
- The **limiting magnitude** depends on sky brightness, atmospheric conditions, and the visual acuity of the observer

**Introduced functions:**

- `rise_trans(jd, body, rsmi, geopos)` — finds the next rising, setting, or meridian transit of a celestial body
- `rise_trans_true_hor(jd, body, rsmi, geopos, horhgt=0.0)` — like `rise_trans` but with a custom horizon altitude
- `refrac(altitude, pressure, temperature, calc_flag)` — converts between true and apparent altitude (or vice versa), taking refraction into account
- `refrac_extended(altitude, altitude_geo, ...)` — extended refraction with dip of the horizon for observers at high altitudes
- `heliacal_ut(jd, geopos, datm, dobs, object_name, event_type)` — finds the date of a heliacal event (first/last visibility)
- `heliacal_ut(jd, geopos, datm, dobs, object_name, event_type)` — same function (bare name is now an alias)
- `vis_limit_mag(jd, geopos, atmo, observer, objname)` — determines whether an object is visible by comparing limiting magnitude and apparent magnitude
- `calc_airmass(altitude)` — calculates the air mass the light passes through
- `calc_extinction_magnitude(altitude)` — calculates how many magnitudes of light are lost due to atmospheric extinction
- `get_twilight_phase(sun_altitude)` — returns the twilight phase based on the Sun's altitude
