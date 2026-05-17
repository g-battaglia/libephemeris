# Chapter 2 — Measuring time in astronomy

## What you will learn

In this chapter, you will discover why time is much more complicated than it seems. You will learn what the Julian Day is and why astronomers use it, the difference between Universal Time (UT) and Terrestrial Time (TT), what Delta-T is and why it matters, and how sidereal time and local time work.

---

## 2.1 The Julian Day (JD)

### The problem with dates

How many days have passed between March 15, 44 BC (the assassination of Caesar) and July 20, 1969 (the moon landing)? To answer, you would have to handle the Julian calendar, the Gregorian calendar, the fact that 10 days "disappeared" in 1582, leap years with different rules in the two calendars, and the fact that there is no year zero (from 1 BC you go straight to 1 AD).

Astronomers solved this chaos in the simplest possible way: a single counter.

### The solution: a number for every instant

The **Julian Day** (JD) is a continuous count of days starting from noon on January 1, 4713 BC (Julian calendar). Every instant in history has a unique JD number:

| Event | JD |
|--------|-----|
| Start of the count (Jan 1, 4713 BC, noon) | 0.0 |
| Birth of Christ (conventional) | ~1 721 424 |
| Gregorian reform (Oct 15, 1582) | 2 299 161 |
| J2000.0 (Jan 1, 2000, TT noon) | 2 451 545.0 |
| Solar eclipse Apr 8, 2024, noon | 2 460 408.5 |

An important detail: the Julian day starts at **noon**, not at midnight. Therefore JD 2460408**.0** is midnight on April 8, 2024, and JD 2460408**.5** is noon. This is a historical remnant: astronomers worked at night, and starting the day at noon avoided changing the date in the middle of an observing session.

### 💻 Code

```python
import libephemeris as ephem

# Date -> Julian Day
jd = ephem.julday(2024, 4, 8, 12.0)
print(f"JD = {jd}")

# Julian Day -> Date
year, month, day, hour = ephem.revjul(jd)
print(f"{day}/{month}/{year} hour {hour:.1f}")

# Julian vs Gregorian calendar
# October 4, 1582 (Julian) is the day before October 15, 1582 (Gregorian)
jd_julian = ephem.julday(1582, 10, 4, 12.0, ephem.JUL_CAL)
jd_gregorian = ephem.julday(1582, 10, 15, 12.0, ephem.GREG_CAL)
print(f"Difference: {jd_gregorian - jd_julian} days")

# How many days between Caesar and the moon landing?
jd_caesar = ephem.julday(-43, 3, 15, 12.0, ephem.JUL_CAL)
jd_moon = ephem.julday(1969, 7, 20, 20.3)  # 20:17 UT
print(f"Days elapsed: {jd_moon - jd_caesar:.0f}")
```

```
JD = 2460409.0
8/4/2024 hour 12.0
Difference: 1.0 days
Days elapsed: 734997
```

> **Note on negative years**: LibEphemeris uses the astronomical convention where year 0 exists. 1 BC = year 0, 2 BC = year -1, and so on. So 44 BC = year -43.

---

## 2.2 UT, TT and Delta-T

### Two different times for two different purposes

**Universal Time** (UT, more precisely UT1) is based on the Earth's rotation. A UT day is the time it takes the Earth to make a complete rotation relative to the mean Sun. It is the time of our clocks, time zones, and daily life.

The problem is that the Earth is not a perfect clock. Its rotation gradually slows down (due to lunar tides) and has unpredictable irregularities. A UT day is not always the same.

**Terrestrial Time** (TT) is a perfectly uniform time, based on atomic clocks. A TT second is always exactly equal to another. Planets move according to the laws of physics, which operate in uniform time — this is why ephemerides are calculated in TT.

### Delta-T: the bridge between the two times

**Delta-T** (ΔT) is the difference between TT and UT:

    ΔT = TT − UT

Today Delta-T is about 69 seconds. In 1900 it was about 3 seconds. In 1800 about 14 seconds. Going back through the centuries, the uncertainty on Delta-T grows rapidly.

Why does it matter? Imagine calculating the position of the Moon for a date in 1800. The Moon moves by about 0.5" per second. If Delta-T has an error of 1 second, the position of the Moon has an error of 0.5". For dates in the Middle Ages, the uncertainty on Delta-T can be minutes, and the error on the Moon can be tens of arcminutes.

### The `_ut` functions vs without suffix

LibEphemeris offers two versions of many functions:

- `calc_ut(jd_ut, ...)`: accepts a Julian Day in **UT** — converts internally to TT by adding Delta-T
- `calc(jd_tt, ...)`: accepts a Julian Day in **TT** — for those who want to handle Delta-T manually

For the vast majority of uses, `calc_ut` is the right choice: you pass the date as you know it (in UT, i.e., "civil" time) and the library does the rest.

### 💻 Code

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Delta-T in days (the library's convention)
dt_days = ephem.deltat(jd)
dt_seconds = dt_days * 86400

print(f"Delta-T = {dt_seconds:.2f} seconds")

# The two equivalent calls:
pos_ut, _ = ephem.calc_ut(jd, ephem.MOON, 0)        # pass JD in UT
pos_tt, _ = ephem.calc(jd + dt_days, ephem.MOON, 0)  # pass JD in TT

print(f"Moon (via UT):  {pos_ut[0]:.6f}°")
print(f"Moon (via TT):  {pos_tt[0]:.6f}°")
```

```
Delta-T = 69.20 seconds
Moon (via UT):  15.429860°
Moon (via TT):  15.429860°
```

### 🌍 Real life

For a birth chart from 1800, Delta-T was about 14 seconds. The Moon moves by about 0.5" per second, so the error introduced by a wrong Delta-T of 1 second is ~0.5" on the Moon — negligible for astrology, but not for calculating historical eclipses. For ancient dates (before 1600), the uncertainty on Delta-T can be several minutes, making it impossible to calculate precise lunar positions.

---

## 2.3 TAI, UTC and leap seconds

Besides UT and TT, there are two other important times:

**UTC** (Coordinated Universal Time) is the time of our clocks and computers. It is based on atomic clocks, but is periodically adjusted with the addition of a **leap second** to stay within 0.9 seconds of UT1. Since 1972, 27 leap seconds have been added.

**TAI** (International Atomic Time) is pure atomic time, without adjustments. TAI runs continuously without ever stopping or skipping. The relationship is:

    TAI = UTC + N seconds    (where N is the number of accumulated leap seconds)
    TT  = TAI + 32.184 s     (fixed constant, forever)

Today (2024): TAI = UTC + 37 seconds, so TT = UTC + 69.184 seconds.

### 💻 Code

```python
import libephemeris as ephem

# Convert a precise UTC date to TAI JD
# April 8, 2024, 18:23:45.123 UTC
jd_tai = ephem.utc_to_tai_jd(2024, 4, 8, 18, 23, 45.123)
print(f"JD (TAI) = {jd_tai:.8f}")

# Convert from UTC to JD in TT and UT1
jd_tt, jd_ut = ephem.utc_to_jd(2024, 4, 8, 18, 23, 45.123)
print(f"JD (TT)  = {jd_tt:.8f}")
print(f"JD (UT1) = {jd_ut:.8f}")
print(f"TT - UT1 = {(jd_tt - jd_ut) * 86400:.2f} seconds")

# Convert back from JD TT to UTC components
year, month, day, hour, minute, sec = ephem.jdet_to_utc(jd_tt)
print(f"UTC: {day}/{month}/{year} {hour}:{minute}:{sec:.3f}")
```

```
JD (TAI) = 2460409.26692272
JD (TT)  = 2460409.26729522
JD (UT1) = 2460409.26649429
TT - UT1 = 69.20 seconds
UTC: 8/4/2024 18:23:45.123
```

### 🌍 Real life

Your smartphone shows the time in UTC (plus the local time zone). Every now and then — the last time on December 31, 2016 — a leap second is added: the clock shows 23:59:60 before moving to 00:00:00. Many computer systems do not handle this case well, and over the years leap seconds have caused server crashes, cloud service slowdowns, and bugs in GPS navigation systems.

---

## 2.4 Sidereal time

### The sidereal day

The **sidereal day** is the time it takes the Earth to make a complete rotation relative to the stars (not relative to the Sun). It lasts **23 hours, 56 minutes and 4 seconds** — about 4 minutes less than the solar day.

Why the difference? While the Earth rotates on its axis, it also moves along its orbit around the Sun. After a complete rotation relative to the stars, the Earth has traveled a short distance along its orbit, and must rotate another ~1° to "catch up" and realign with the Sun. That extra 1° takes about 4 minutes.

**Local sidereal time** measures how much rotation has occurred relative to the stars. In practice, it indicates which part of the celestial sphere is on the meridian at that moment. The **ARMC** (Right Ascension of the Midheaven) is the local sidereal clock, expressed in degrees (0°–360°).

Sidereal time is fundamental for calculating astrological houses: the houses depend on how the celestial sphere is oriented with respect to the local horizon, and this is exactly what sidereal time measures.

### 💻 Code

```python
import libephemeris as ephem

# April 8, 2024, 21:00 UT
jd = ephem.julday(2024, 4, 8, 21.0)

# Greenwich sidereal time (in hours)
st_greenwich = ephem.sidtime(jd)
print(f"Greenwich sidereal time: {st_greenwich:.4f} hours")

# Local sidereal time (Milan, lon 9.19° E)
st_milan = ephem.sidtime(jd, longitude=9.19)
print(f"Milan sidereal time:     {st_milan:.4f} hours")

# The difference is proportional to longitude:
# 9.19° / 15 = 0.6127 hours ≈ 36 minutes and 46 seconds
```

```
Greenwich sidereal time: 10.1738 hours
Milan sidereal time:     10.7865 hours
```

### 🌍 Real life

This is why constellations "rise earlier" every evening: every night, at the same time, the celestial sphere has rotated ~1° more (the 4 minutes difference between solar and sidereal day). After a month, the difference is about 2 hours. After 6 months, the night sky shows completely different constellations.

---

## 2.5 Local time, time zones and LMT

### Local solar time

Before time zones (introduced in 1884), every city used its own **local solar time**. Noon was when the Sun crossed the local meridian — and this depends on longitude.

**Local Mean Time** (LMT) is the mean solar time for a given longitude. The difference with respect to Greenwich is simply the longitude divided by 15 (because 360° / 24h = 15°/h):

    LMT = UT + longitude / 15    (in hours)

In Milan (longitude 9.19° E), the mean solar noon arrives about 37 minutes after Greenwich.

**Local Apparent Time** (LAT) is sundial time — the *true* solar time, not the mean. The difference between LAT and LMT is the **equation of time**: an oscillation of ±16 minutes during the year, caused by the eccentricity of the Earth's orbit and the obliquity of the ecliptic.

### 💻 Code

```python
import libephemeris as ephem

# April 8, 2024, 12:00 UT
jd = ephem.julday(2024, 4, 8, 12.0)

# Equation of time (returned in days)
eot_days = ephem.time_equ(jd)
eot_minutes = eot_days * 1440  # convert to minutes

print(f"Equation of time: {eot_minutes:.2f} minutes")
# If positive: the sundial is ahead of the clock

# LMT -> LAT conversion for Milan (lon 9.19° E)
jd_lmt = ephem.julday(2024, 4, 8, 12.0)  # LMT noon
jd_lat = ephem.lmt_to_lat(jd_lmt, 9.19)

# The difference is the equation of time
diff_minutes = (jd_lat - jd_lmt) * 1440
print(f"LAT - LMT = {diff_minutes:.2f} minutes")

# UTC to local time conversion
# Rome: UTC+2 in summer (CEST)
# utc_time_zone converts UTC -> local time by adding the offset
year, month, day, hour, minute, sec = ephem.utc_time_zone(
    2024, 4, 8, 14, 30, 0.0, 2.0  # 14:30 UTC, offset +2 for CEST
)
print(f"Rome local time: {hour}:{minute}:{sec:.0f}")  # 16:30 CEST
```

```
Equation of time: -0.42 minutes
LAT - LMT = -0.42 minutes
Rome local time: 16:30:0
```

### 🌍 Real life

In Milan (longitude 9.19° E, CET time zone = UTC+1), the mean solar noon falls around 12:23 CET in winter (UTC+1) and around 13:23 CEST in summer (UTC+2). But the equation of time adds a further oscillation: in early November the true solar noon is about 16 minutes before the mean solar noon, in mid-February about 14 minutes after.

---

## 2.6 IERS and observed Delta-T

### The prediction problem

Delta-T changes over time in a way that is not perfectly predictable, because it depends on the irregularities of Earth's rotation. For the recent past (from 1962 onwards), we have **precise measurements** of Delta-T thanks to the IERS (International Earth Rotation and Reference Systems Service).

By default, LibEphemeris uses Skyfield's Delta-T model, which combines historical data with predictions. But for maximum precision on recent dates, you can use the observed IERS data:

### 💻 Code

```python
import libephemeris as ephem

# Download IERS data (one-time, ~1 MB)
ephem.download_delta_t_data()

# Enable the use of observed IERS data
ephem.set_iers_delta_t_enabled(True)

# Verify that data is available for a date
jd = ephem.julday(2024, 4, 8, 12.0)
available = ephem.is_iers_data_available(jd)
print(f"IERS data available: {available}")

# Comparison: Calculated vs observed Delta-T
ephem.set_iers_delta_t_enabled(False)
dt_calculated = ephem.deltat(jd) * 86400  # in seconds

ephem.set_iers_delta_t_enabled(True)
dt_observed = ephem.deltat(jd) * 86400  # in seconds

print(f"Calculated Delta-T: {dt_calculated:.3f} s")
print(f"Observed Delta-T: {dt_observed:.3f} s")
print(f"Difference:         {abs(dt_observed - dt_calculated):.3f} s")

# For recent dates the difference is small (< 0.1 s)
# For dates far in the future, there is no IERS data
# and the library falls back to the calculated model
```

```
IERS data available: False
Calculated Delta-T: 69.200 s
Observed Delta-T: 69.199 s
Difference:         0.001 s
```

IERS data is updated weekly. The library can download it automatically:

```python
# Enable automatic download of IERS data
ephem.set_iers_auto_download(True)

# Or manually download everything
ephem.download_iers_finals()   # Earth Orientation data
ephem.download_leap_seconds()  # leap seconds table
ephem.download_delta_t_data()  # Delta-T historical series
```

---

## Summary

- The **Julian Day** (JD) is a continuous count of days, used to avoid calendar complications. `julday()` converts a date to JD, `revjul()` does the reverse.
- **UT** (Universal Time) is the time based on Earth's rotation — "civil" time. **TT** (Terrestrial Time) is the uniform time of atomic clocks.
- **Delta-T** = TT − UT, today ~69 seconds. The `_ut` functions accept JD in UT and convert automatically.
- **UTC** is the time of our clocks, with leap seconds. **TAI** is pure atomic time. TT = TAI + 32.184 s.
- **Sidereal time** measures rotation relative to the stars (~4 min shorter than the solar day). It is the basis for calculating astrological houses.
- **Local time** (LMT) depends on longitude. The **equation of time** is the difference between sundial time and mean time.
- **IERS data** provides observed Delta-T for maximum precision on recent dates.

### Introduced functions and constants

| Function / Constant | Usage |
|---------------------|-------|
| `julday(year, month, day, hour, gregflag)` | Date → Julian Day |
| `revjul(jd, gregflag)` | Julian Day → date (year, month, day, hour) |
| `deltat(jd)` | Delta-T in days |
| `calc_ut(jd_ut, body, flag)` | Calculation with JD in UT (converts internally to TT) |
| `calc(jd_tt, body, flag)` | Calculation with JD in TT |
| `utc_to_jd(year, month, day, hour, min, sec)` | UTC → JD in TT and UT1 |
| `utc_to_tai_jd(year, month, day, hour, min, sec)` | UTC → JD in TAI |
| `sidtime(jd, longitude)` | Local sidereal time in hours |
| `time_equ(jd)` | Equation of time in days |
| `lmt_to_lat(jd_lmt, longitude)` | Local Mean Time → Local Apparent Time |
| `utc_time_zone(year, month, day, hour, min, sec, offset)` | Time zone → UTC |
| `set_iers_delta_t_enabled(True/False)` | Enable/disable IERS data |
| `download_delta_t_data()` | Download observed Delta-T data |
| `SE_GREG_CAL`, `SE_JUL_CAL` | Gregorian / Julian calendar |
