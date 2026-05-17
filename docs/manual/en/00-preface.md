# Preface and how to use this manual

## What you will learn

In these few introductory pages you will discover who this manual is for, how it is organized, and you will perform your first astronomical calculation in three lines of Python.

---

## Who is it for

This manual is for three types of readers:

- **Developers** who need to integrate astronomical or astrological calculations into an application and want to understand what's behind the numbers, not just copy a line of code.
- **Astrologers** who want to switch from commercial software to an open and verifiable tool, without giving up precision.
- **Astronomy enthusiasts** who wonder how to calculate where a planet is, when the Sun rises, or why an eclipse is visible only in certain places.

You don't need to know anything about astronomy. You only need a bit of Python — the level of someone who knows what a function, a dictionary, and a `for` loop are.

## What you will learn reading this manual

This is not an API reference manual. It is a book that explains the **concepts** behind the calculations. The difference is important: an API reference tells you "the `calc_ut` function accepts three parameters"; this manual explains **why** those three parameters exist, what they mean physically, and what mistakes you could make if you don't understand them.

By the end of the manual you will know:

- How the celestial sphere works and why planets move the way they do
- Why time in astronomy is complicated and how to avoid traps
- What ecliptic, equatorial, and horizontal coordinates are
- How to calculate the position of any celestial body with sub-arcsecond precision
- How astrological houses work and why 20 different systems exist
- How to find eclipses, sunrises, sunsets, and visibility phenomena
- How to optimize performance for large-scale calculations

## How it is organized

Each chapter follows the same structure:

- **📐 Theory**: the astronomical or astrological concept, explained with analogies and without unnecessary formulas
- **🌍 Real life**: a concrete example connecting the concept to something observable or familiar
- **💻 Code**: a working example with LibEphemeris, which you can copy and paste

The chapters are designed to be read in order. The first three build the foundations (the sky, time, coordinates) and everything else rests on them. If you are already familiar with positional astronomy, you can skip straight to Chapter 5 and go back when needed.

Chapter 15 — the Cookbook — is designed for those in a hurry: copy-paste recipes for the most common calculations, with brief explanations.

## Prerequisites

- **Python 3.9** or higher
- **No astronomical knowledge** — this manual starts from scratch
- **No other packages** — LibEphemeris installs with a single command and has no heavy dependencies

## Installation

```bash
pip install libephemeris
```

On first use, the library automatically downloads the necessary ephemeris files (the DE440 file, about 114 MB). You can also force the download in advance:

```python
import libephemeris as ephem
ephem.download_for_tier("medium")
```

The three precision tiers — `base`, `medium`, `extended` — differ in the time range covered and file size. The `medium` tier (default) covers the years 1550–2650 and is fine for the vast majority of uses. We will talk about it in detail in Chapter 4.

## Your first calculation

Let's calculate where the Sun was during the total solar eclipse of April 8, 2024:

```python
import libephemeris as ephem

# Let's convert the date to Julian Day (noon UT)
jd = ephem.julday(2024, 4, 8, 12.0)

# Let's calculate the position of the Sun
# The third argument (0) means: ecliptic coordinates, without extra options
pos, flag = ephem.calc_ut(jd, ephem.SUN, 0)

longitude = pos[0]  # degrees along the ecliptic (0–360)
print(f"Sun at {longitude:.4f}° of ecliptic longitude")
```

```
Sun at 19.1404° of ecliptic longitude
```

What happened?

1. `julday(2024, 4, 8, 12.0)` converted the date "April 8, 2024, 12:00 UT" into a **Julian Day** — a unique number that identifies that instant. We will learn all about Julian Days in Chapter 2.

2. `calc_ut(jd, ephem.SUN, 0)` calculated the Sun's position for that instant. The result is a tuple of 6 numbers: longitude, latitude, distance, and their respective daily velocities. The third parameter (`0`) tells the library to use standard settings. Calculation flags are covered in Chapter 5.

3. The longitude `19.15°` tells us that the Sun was at about 19 degrees of the sign of Aries (the first sign goes from 0° to 30°). We will understand why in Chapter 1.

Let's take it one step further — let's show the position in zodiacal format:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)
pos, flag = ephem.calc_ut(jd, ephem.SUN, 0)

signs = [
    "Aries", "Taurus", "Gemini", "Cancer",
    "Leo", "Virgo", "Libra", "Scorpio",
    "Sagittarius", "Capricorn", "Aquarius", "Pisces",
]

# split_deg breaks down decimal degrees into degrees, minutes, seconds and sign
deg, min, sec, secfr, sign = ephem.split_deg(
    pos[0], ephem.SPLIT_DEG_ZODIACAL | ephem.SPLIT_DEG_ROUND_SEC
)

print(f"Sun at {deg}° {min}' {sec}\" {signs[sign]}")
```

```
Sun at 19° 8' 26" Aries
```

The `split_deg` function takes a value in decimal degrees and breaks it down into degrees, arc minutes, arc seconds, fraction of a second, and the number of the zodiac sign. It is the same format you find in astrology magazines or astronomical almanacs.

## Conventions used in this manual

Throughout the manual:

- **The import is always the same**: `import libephemeris as ephem`. All functions and constants are accessible via `ephem.name`.
- **The examples are complete**: every code block can be copied and pasted into a `.py` file and it works. No partial snippets.
- **The units are always indicated**: degrees (°), arc minutes ('), arc seconds ("), astronomical units (AU), kilometers (km). If you see a number without a unit, it's a mistake.
- **The dates in the examples are real events**: eclipses, conjunctions, equinoxes that actually happened. You can verify the results with any other astronomical software.

## Manual structure

| #  | Chapter | Topic |
|----|---------|-------|
| 1  | The sky seen from Earth | The celestial sphere, the ecliptic, the zodiac, the horizon |
| 2  | Measuring time | Julian Day, UT, TT, Delta-T, sidereal time |
| 3  | Celestial coordinates | Ecliptic, equatorial, horizontal, conversions |
| 4  | Ephemeris | JPL DE440, precision tiers, geocentric vs heliocentric |
| 5  | Planet positions | `calc_ut`, flags, velocities, retrogradation, phenomena |
| 6  | The Moon | Nodes, apogee, perigee, Lilith, crossings |
| 7  | Astrological houses | House systems, Ascendant, MC, extreme latitudes |
| 8  | Fixed stars | Catalog, proper motion, stars in astrology |
| 9  | Eclipses | Solar, lunar, Saros cycles, occultations |
| 10 | Sunrise and sunset | Rising, transits, twilights, refraction, heliacal visibility |
| 11 | Sidereal zodiac | Tropical vs sidereal, ayanamsha, Vedic calculations |
| 12 | Minor bodies | Asteroids, centaurs, TNOs, SPK kernels, Keplerian fallback |
| 13 | Hypothetical planets | Uranians, Transpluto, Arabic parts |
| 14 | Precision and performance | Comparisons, LEB mode, configuration |
| 15 | Cookbook | Birth chart, transits, eclipses, solar returns and more |

---

## Summary

- LibEphemeris is a Python library for high-precision astronomical and astrological calculations
- It installs with `pip install libephemeris` and uses NASA JPL DE440 data
- The base function is `calc_ut(jd, body, flag)`: given an instant and a celestial body, it returns the position
- `julday(year, month, day, hours)` converts a date to a Julian Day
- `split_deg(degrees, flag)` formats decimal degrees into degrees, minutes, seconds and zodiac sign
- Every chapter in this manual combines theory, real-life examples, and working code
