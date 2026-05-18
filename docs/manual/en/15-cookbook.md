# Chapter 15 — Cookbook: Practical Calculations from A to Z

## What You Will Learn

This chapter is a collection of **complete recipes**, ready to copy and paste. Each recipe solves a common practical problem: from the birth chart to transits, from retrogrades to eclipses. There is no theory to read — just working code with its output.

---

## Recipe 1 — Complete Birth Chart

**Problem**: Given the date, time, and place of birth, calculate all planetary positions and houses.

**What it is used for**: The birth chart is the foundation of any astrological analysis. This code produces the "photograph of the sky" at the moment of birth: the positions of the planets in the signs, the twelve houses, and the main angles (Ascendant and Midheaven). From here all interpretations begin — from personality (Sun position) to emotions (Moon), to the way of communicating (Mercury), up to the structure of life areas (houses). A concrete application: an astrological software that automatically generates the birth chart starting from the user's personal data.

```python
import libephemeris as ephem

# Rome, April 8, 2024, 14:30 UT
jd = ephem.julday(2024, 4, 8, 14.5)

signs = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

def zodiac_format(lon):
    sign_idx = int(lon / 30)
    degrees = lon % 30
    d, m, s, sf, si = ephem.split_deg(degrees, 0)
    return f'{d:2d}° {m:02d}\' {s:02d}" {signs[sign_idx]}'

bodies = [
    (ephem.SUN,       "Sun"),
    (ephem.MOON,      "Moon"),
    (ephem.MERCURY,   "Mercury"),
    (ephem.VENUS,     "Venus"),
    (ephem.MARS,      "Mars"),
    (ephem.JUPITER,   "Jupiter"),
    (ephem.SATURN,    "Saturn"),
    (ephem.URANUS,    "Uranus"),
    (ephem.NEPTUNE,   "Neptune"),
    (ephem.PLUTO,     "Pluto"),
    (ephem.MEAN_NODE, "North Node"),
    (ephem.CHIRON,    "Chiron"),
]

print("--- Birth Chart: Rome, April 8, 2024, 14:30 ---")
print()
for body_id, name in bodies:
    pos, _ = ephem.calc_ut(jd, body_id, ephem.FLG_SPEED)
    retro = " R" if pos[3] < 0 else ""
    print(f"{name:11s}  {zodiac_format(pos[0])}{retro}")

# Placidus Houses — Rome (lat 41.9, lon 12.5)
cusps, ascmc = ephem.houses(jd, 41.9, 12.5, ord("P"))
print()
print(f"Ascendant:   {zodiac_format(ascmc[0])}")
print(f"Midheaven:   {zodiac_format(ascmc[1])}")
print()
for i in range(12):
    print(f"House {i+1:2d}:    {zodiac_format(cusps[i])}")
```

```
--- Birth Chart: Rome, April 8, 2024, 14:30 ---

Sun          19° 14' 34" Ari
Moon         16° 59' 41" Ari
Mercury      24° 53' 58" Ari R
Venus         4° 14' 49" Ari
Mars         12° 55' 36" Psc
Jupiter      19° 00' 36" Tau
Saturn       14° 26' 16" Psc
Uranus       21° 09' 46" Tau
Neptune      28° 11' 03" Psc
Pluto         1° 57' 56" Aqr
North Node   15° 39' 19" Ari R
Chiron       19° 23' 45" Ari

Ascendant:   12° 15' 00" Vir
Midheaven:    9° 02' 28" Gem

House  1:    12° 15' 00" Vir
House  2:     6° 21' 02" Lib
House  3:     5° 32' 08" Sco
House  4:     9° 02' 28" Sgr
House  5:    13° 32' 25" Cap
House  6:    15° 02' 03" Aqr
House  7:    12° 15' 00" Psc
House  8:     6° 21' 02" Ari
House  9:     5° 32' 08" Tau
House 10:     9° 02' 28" Gem
House 11:    13° 32' 25" Cnc
House 12:    15° 02' 03" Leo
```

---

## Recipe 2 — Daily Transits

**Problem**: Given a birth chart, find which transiting planets are in aspect with natal planets today.

**What it is used for**: Transits are the most widespread method of astrological forecasting. The principle is simple: "today's" planets form significant angles (conjunction, square, trine...) with the positions of the birth chart, activating the themes represented by those planets. For example, Saturn squaring the natal Sun indicates a period of challenges and restructuring of one's identity — typically lasting a few weeks. A concrete application: an app that every morning shows the user the active transits on their chart, along with the orb (how close it is to the exact aspect).

```python
import libephemeris as ephem

# Natal positions: born March 15, 1990, 10:00 UT
jd_natal = ephem.julday(1990, 3, 15, 10.0)
# Transit date: April 8, 2024
jd_transit = ephem.julday(2024, 4, 8, 12.0)

planets = [
    (ephem.SUN, "Sun"), (ephem.MOON, "Moon"),
    (ephem.MERCURY, "Mercury"), (ephem.VENUS, "Venus"),
    (ephem.MARS, "Mars"), (ephem.JUPITER, "Jupiter"),
    (ephem.SATURN, "Saturn"),
]

# Calculate natal and transit positions
natal_pos = {}
for body_id, name in planets:
    pos, _ = ephem.calc_ut(jd_natal, body_id, 0)
    natal_pos[name] = pos[0]

transit_pos = {}
for body_id, name in planets:
    pos, _ = ephem.calc_ut(jd_transit, body_id, 0)
    transit_pos[name] = pos[0]

# Define aspects and orbs
aspects = {0: "Conjunction", 60: "Sextile", 90: "Square",
           120: "Trine", 180: "Opposition"}
orbs = {0: 8, 60: 6, 90: 7, 120: 7, 180: 8}

print("--- Transits of April 8, 2024, on March 15, 1990 chart ---")
print()
for t_name, t_lon in transit_pos.items():
    for n_name, n_lon in natal_pos.items():
        diff = abs(ephem.difdeg2n(t_lon, n_lon))
        for asp_deg, asp_name in aspects.items():
            orb = abs(diff - asp_deg)
            if orb <= orbs[asp_deg]:
                print(f"{t_name:10s} {asp_name:14s} {n_name:10s}"
                      f"  (orb: {orb:.1f}°)")
                break
```

```
--- Transits of April 8, 2024, on March 15, 1990 chart ---

Sun        Square         Saturn      (orb: 4.2°)
Mercury    Square         Saturn      (orb: 1.6°)
Venus      Sextile        Venus       (orb: 4.9°)
Venus      Sextile        Mars        (orb: 1.3°)
Venus      Square         Jupiter     (orb: 2.7°)
Mars       Trine          Moon        (orb: 4.1°)
Jupiter    Sextile        Sun         (orb: 5.6°)
Jupiter    Sextile        Mercury     (orb: 2.1°)
Jupiter    Trine          Saturn      (orb: 4.3°)
Saturn     Trine          Moon        (orb: 5.7°)
Saturn     Conjunction    Mercury     (orb: 6.6°)
```

---

## Recipe 3 — Sun Ingresses into Signs

**Problem**: Find the exact date and time when the Sun enters each zodiac sign during the year.

**What it is used for**: The Sun's ingresses into the signs mark the changes of season (Aries = spring equinox, Cancer = summer solstice, etc.) and are fundamental dates in mundane astrology. The chart erected for the moment of the Aries ingress (the so-called "chart of the year") is used to make predictions for the entire year for a nation. In astronomy, knowing the exact moment of the equinox is useful for calibrating sundials and almanacs. A concrete application: an astrological calendar that automatically shows the sign change dates for any given year.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 1, 1, 0.0)

signs = ["Aries", "Taurus", "Gemini", "Cancer",
         "Leo", "Virgo", "Libra", "Scorpio",
         "Sagittarius", "Capricorn", "Aquarius", "Pisces"]

print("--- Sun Ingresses into Signs (2024) ---")
print()
for i in range(12):
    target = i * 30.0
    jd_ingress = ephem.solcross_ut(target, jd)
    year, month, day, hour = ephem.revjul(jd_ingress)
    hours = int(hour)
    minutes = int((hour - hours) * 60)
    print(f"{signs[i]:12s}  {day:2d}/{month:02d}/{year}"
          f"  {hours:02d}:{minutes:02d} UT")
    jd = jd_ingress + 1
```

```
--- Sun Ingresses into Signs (2024) ---

Aries         20/03/2024  03:06 UT
Taurus        19/04/2024  13:59 UT
Gemini        20/05/2024  12:59 UT
Cancer        20/06/2024  20:50 UT
Leo           22/07/2024  07:44 UT
Virgo         22/08/2024  14:55 UT
Libra         22/09/2024  12:43 UT
Scorpio       22/10/2024  22:14 UT
Sagittarius   21/11/2024  19:56 UT
Capricorn     21/12/2024  09:20 UT
Aquarius      19/01/2025  20:00 UT
Pisces        18/02/2025  10:06 UT
```

---

## Recipe 4 — Mercury Retrograde

**Problem**: Find all the retrograde and direct stations of Mercury in the year, and check if it is retrograde on a specific date.

**What it is used for**: Mercury retrograde is one of the most popular — and feared — astrological events. When Mercury is retrograde (meaning it appears to move backward as seen from Earth), astrological tradition associates the period with misunderstandings in communication, problems with contracts and technology, and travel delays. In practice, many people avoid signing contracts or buying electronics during these periods. This code finds the exact start (retrograde station) and end (direct station) dates of each retrograde, along with the precise zodiacal degree. A concrete application: a digital planner that automatically highlights Mercury retrograde periods.

```python
import libephemeris as ephem
from libephemeris.crossing import find_station_ut, is_retrograde

jd = ephem.julday(2024, 1, 1, 0.0)
jd_end = ephem.julday(2025, 1, 1, 0.0)

signs = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

print("--- Mercury Retrogrades in 2024 ---")
print()
while jd < jd_end:
    jd_station, type_ = find_station_ut(
        ephem.MERCURY, jd, "any"
    )
    if jd_station >= jd_end:
        break
    year, month, day, hour = ephem.revjul(jd_station)
    hours = int(hour)
    minutes = int((hour - hours) * 60)
    pos, _ = ephem.calc_ut(jd_station, ephem.MERCURY, 0)
    sign = signs[int(pos[0] / 30)]
    degrees = pos[0] % 30
    label = ("Retrograde Station" if type_ == "SR"
             else "Direct Station    ")
    print(f"{label}  {day:2d}/{month:02d}/{year}"
          f"  {hours:02d}:{minutes:02d} UT  at {degrees:5.1f}° {sign}")
    jd = jd_station + 1

# Point-in-time check
print()
retro = is_retrograde(ephem.MERCURY,
                      ephem.julday(2024, 4, 8, 12.0))
print(f"Is Mercury retrograde on 8/4/2024? {retro}")
```

```
--- Mercury Retrogrades in 2024 ---

Direct Station       2/01/2024  03:07 UT  at  22.2° Sgr
Retrograde Station   1/04/2024  22:14 UT  at  27.2° Ari
Direct Station      25/04/2024  12:54 UT  at  16.0° Ari
Retrograde Station   5/08/2024  04:56 UT  at   4.1° Vir
Direct Station      28/08/2024  21:13 UT  at  21.4° Leo
Retrograde Station  26/11/2024  02:42 UT  at  22.7° Sgr
Direct Station      15/12/2024  20:56 UT  at   6.4° Sgr

Is Mercury retrograde on 8/4/2024? True
```

> **Note**: `find_station_ut` and `is_retrograde` are imported from `libephemeris.crossing` because they are not exposed at the main module level.

---

## Recipe 5 — New Moon and Full Moon

**Problem**: Find all the New and Full Moons of the year with date, time, and zodiac sign.

**What it is used for**: Lunar phases have set the rhythm of daily life for millennia. The New Moon (Sun-Moon conjunction) is traditionally the time to sow new projects, while the Full Moon (opposition) is the time of harvest and maximum energy — but also of emotional tension. In biodynamic agriculture, sowing still follows the lunar calendar. In astrology, the sign in which the New Moon falls indicates the "theme" of the month: a New Moon in Aries invites courageous action, one in Pisces to let go. A concrete application: a lunar calendar that shows the phases with their signs, useful for gardeners, fishermen, or anyone who follows lunar rhythms.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 1, 1, 0.0)
jd_end = ephem.julday(2025, 1, 1, 0.0)

signs = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

print("--- New and Full Moons in 2024 ---")
print()
while jd < jd_end:
    # Sun's current longitude
    pos_sun, _ = ephem.calc_ut(jd, ephem.SUN, 0)
    lon_sun = pos_sun[0]

    # New Moon: Moon conjuncts the Sun
    jd_new = ephem.mooncross_ut(lon_sun, jd)
    if jd_new < jd_end:
        year, month, day, hour = ephem.revjul(jd_new)
        hours = int(hour)
        minutes = int((hour - hours) * 60)
        pos_moon, _ = ephem.calc_ut(jd_new, ephem.MOON, 0)
        sign = signs[int(pos_moon[0] / 30)]
        degrees = pos_moon[0] % 30
        print(f"New Moon    {day:2d}/{month:02d}/{year}"
              f"  {hours:02d}:{minutes:02d} UT"
              f"  in {degrees:5.1f}° {sign}")

    # Full Moon: Moon opposite the Sun (180°)
    target_full = (lon_sun + 180.0) % 360.0
    jd_full = ephem.mooncross_ut(target_full, jd)
    if jd_full < jd_end and jd_full > jd:
        year, month, day, hour = ephem.revjul(jd_full)
        hours = int(hour)
        minutes = int((hour - hours) * 60)
        pos_moon, _ = ephem.calc_ut(jd_full, ephem.MOON, 0)
        sign = signs[int(pos_moon[0] / 30)]
        degrees = pos_moon[0] % 30
        print(f"Full Moon   {day:2d}/{month:02d}/{year}"
              f"  {hours:02d}:{minutes:02d} UT"
              f"  in {degrees:5.1f}° {sign}")

    print()
    # Advance by ~20 days to find the next pair
    jd = min(jd_new, jd_full) + 20
    if jd >= jd_end:
        break
```

```
--- New and Full Moons in 2024 ---

New Moon    10/01/2024  18:20 UT  in  10.0° Cap
Full Moon   23/01/2024  16:58 UT  in  10.0° Cnc

New Moon     9/02/2024  06:34 UT  in  10.3° Aqr
Full Moon   22/02/2024  10:11 UT  in  10.3° Leo

New Moon     9/03/2024  17:05 UT  in  10.2° Psc
Full Moon   23/03/2024  04:18 UT  in  10.2° Vir

New Moon     8/04/2024  02:31 UT  in   9.5° Ari
Full Moon   21/04/2024  22:20 UT  in   9.5° Lib

New Moon     7/05/2024  11:18 UT  in   8.3° Tau
Full Moon   21/05/2024  14:50 UT  in   8.3° Sco

New Moon     5/06/2024  19:59 UT  in   6.7° Gem
Full Moon   20/06/2024  04:55 UT  in   6.7° Sgr

New Moon     5/07/2024  05:23 UT  in   4.7° Cnc
Full Moon   19/07/2024  16:32 UT  in   4.7° Cap

New Moon     3/08/2024  16:25 UT  in   2.8° Leo
Full Moon   18/08/2024  02:22 UT  in   2.8° Aqr

New Moon     2/09/2024  05:50 UT  in   1.0° Vir
Full Moon   16/09/2024  11:18 UT  in   1.0° Psc

New Moon     1/10/2024  21:45 UT  in  29.7° Vir
Full Moon   15/10/2024  20:07 UT  in  29.7° Psc

New Moon    31/10/2024  15:26 UT  in  29.0° Lib
Full Moon   14/11/2024  05:21 UT  in  29.0° Ari

New Moon    30/11/2024  09:32 UT  in  28.8° Sco
Full Moon   13/12/2024  15:23 UT  in  28.8° Tau

New Moon    30/12/2024  02:45 UT  in  29.0° Sgr
```

---

## Recipe 6 — Next Solar Eclipse from My City

**Problem**: Find the next solar eclipses visible from a specific location, with contact times and magnitude.

**What it is used for**: Solar eclipses are spectacular but rare events for any given location — decades can pass between one total eclipse and the next. Knowing in advance when and how an eclipse will be visible from your city is essential for organizing observations (solar filters, travel, vacations). The magnitude indicates how much of the solar disk will be covered: 1.0 = total, 0.5 = half, 0.01 = barely perceptible. A concrete application: a website that, given the user's location, shows the next visible eclipse with a map of the contacts.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 1, 1, 0.0)
lat, lon, alt = 41.9, 12.5, 0.0  # Rome

print("--- Next Solar Eclipses Visible from Rome ---")
print()
found = 0
for _ in range(20):  # search up to 20 eclipses ahead
    ecl_type, times, attr = ephem.sol_eclipse_when_loc(
        jd, lat, lon, alt
    )
    if ecl_type > 0 and attr[0] > 0.01:
        jd_max = times[0]
        year, month, day, hour = ephem.revjul(jd_max)
        hours = int(hour)
        minutes = int((hour - hours) * 60)

        type_ = "Partial"
        if ecl_type & 4: type_ = "Annular"
        if ecl_type & 1: type_ = "Total"

        print(f"{type_:10s}  {day:2d}/{month:02d}/{year}"
              f"  {hours:02d}:{minutes:02d} UT"
              f"  magnitude: {attr[0]:.3f}")

        # Start and end times
        if times[1] > 0:
            _, _, _, h = ephem.revjul(times[1])
            print(f"  Start:   {int(h):02d}:"
                  f"{int((h - int(h)) * 60):02d} UT")
        if times[4] > 0:
            _, _, _, h = ephem.revjul(times[4])
            print(f"  End:     {int(h):02d}:"
                  f"{int((h - int(h)) * 60):02d} UT")
        print()
        found += 1
        if found >= 3:
            break

    jd = times[0] + 30 if times[0] > 0 else jd + 180
```

```
--- Next Solar Eclipses Visible from Rome ---

Partial     29/03/2025  11:03 UT  magnitude: 0.073
  Start:   10:35 UT
  End:     11:31 UT

Partial      2/08/2027  09:12 UT  magnitude: 0.786
  Start:   08:02 UT
  End:     10:26 UT

Partial      1/06/2030  05:04 UT  magnitude: 0.816
  Start:   04:02 UT
  End:     06:12 UT
```

---

## Recipe 7 — Sidereal Chart (Vedic)

**Problem**: Calculate planetary positions in the sidereal zodiac with Lahiri ayanamsha, including the Nakshatras.

**What it is used for**: Vedic astrology (Jyotish) uses the sidereal zodiac instead of the tropical one — positions are shifted by about 24° compared to Western ones. Nakshatras are 27 "lunar mansions" of 13°20' each, further divided into 4 Padas; the position of the Moon in the Nakshatra is fundamental for calculating the Dashas (planetary periods that govern life). For example, someone with the Moon in the Ashwini Nakshatra will live their first period under Ketu (7 years), then under Venus (20 years), etc. A concrete application: a Jyotish software that generates the sidereal chart, the Nakshatras, and the calculation of the Dashas.

```python
import libephemeris as ephem
from libephemeris.constants import (SIDM_LAHIRI,
                                    FLG_SIDEREAL, FLG_SPEED)

jd = ephem.julday(2024, 4, 8, 14.5)
ephem.set_sid_mode(SIDM_LAHIRI)

nakshatras = [
    "Ashwini", "Bharani", "Krittika", "Rohini",
    "Mrigashira", "Ardra", "Punarvasu", "Pushya",
    "Ashlesha", "Magha", "Purva Phalguni",
    "Uttara Phalguni", "Hasta", "Chitra", "Swati",
    "Vishakha", "Anuradha", "Jyeshtha", "Mula",
    "Purva Ashadha", "Uttara Ashadha", "Shravana",
    "Dhanishta", "Shatabhisha", "Purva Bhadra",
    "Uttara Bhadra", "Revati",
]

signs = ["Mesh", "Vrish", "Mithun", "Kark", "Simha",
         "Kanya", "Tula", "Vrischik", "Dhanu", "Makar",
         "Kumbh", "Meen"]

bodies = [
    (ephem.SUN,       "Surya"),
    (ephem.MOON,      "Chandra"),
    (ephem.MARS,      "Mangal"),
    (ephem.MERCURY,   "Budha"),
    (ephem.JUPITER,   "Guru"),
    (ephem.VENUS,     "Shukra"),
    (ephem.SATURN,    "Shani"),
    (ephem.MEAN_NODE, "Rahu"),
]

ayan = ephem.get_ayanamsa_ut(jd)
print("--- Sidereal Chart (Lahiri) ---")
print(f"Ayanamsha: {ayan:.4f}°")
print()

for body_id, name in bodies:
    pos, _ = ephem.calc_ut(jd, body_id,
                           FLG_SPEED | FLG_SIDEREAL)
    lon = pos[0]
    sign = signs[int(lon / 30)]
    degrees = lon % 30
    nak_idx = int(lon / (360 / 27))
    nak = nakshatras[nak_idx]
    pada = int((lon % (360 / 27)) / (360 / 108)) + 1
    print(f"{name:8s}  {degrees:5.1f}° {sign:10s}"
          f"  Nakshatra: {nak} (Pada {pada})")

# Ketu = opposite of Rahu
pos_rahu, _ = ephem.calc_ut(jd, ephem.MEAN_NODE,
                            FLG_SPEED | FLG_SIDEREAL)
ketu = (pos_rahu[0] + 180) % 360
sign_k = signs[int(ketu / 30)]
degrees_k = ketu % 30
nak_k = nakshatras[int(ketu / (360 / 27))]
pada_k = int((ketu % (360 / 27)) / (360 / 108)) + 1
print(f"Ketu      {degrees_k:5.1f}° {sign_k:10s}"
      f"  Nakshatra: {nak_k} (Pada {pada_k})")

ephem.set_sid_mode(0)  # reset
```

```
--- Sidereal Chart (Lahiri) ---
Ayanamsha: 24.1961°

Surya      25.0° Meen        Nakshatra: Revati (Pada 3)
Chandra    22.8° Meen        Nakshatra: Revati (Pada 2)
Mangal     18.7° Kumbh       Nakshatra: Shatabhisha (Pada 4)
Budha       0.7° Mesh        Nakshatra: Ashwini (Pada 1)
Guru       24.8° Mesh        Nakshatra: Bharani (Pada 4)
Shukra     10.1° Meen        Nakshatra: Uttara Bhadra (Pada 3)
Shani      20.2° Kumbh       Nakshatra: Purva Bhadra (Pada 1)
Rahu       21.5° Meen        Nakshatra: Revati (Pada 2)
Ketu       21.5° Kanya       Nakshatra: Hasta (Pada 4)
```

---

## Recipe 8 — Visibility of Planets Tonight

**Problem**: For every planet visible to the naked eye, calculate its altitude, direction, magnitude, and elongation from the Sun one hour after sunset.

**What it is used for**: Knowing which planets are visible tonight is the most practical question for stargazers. The altitude above the horizon tells whether the planet is high enough to be seen (above buildings and trees); the elongation from the Sun tells whether it is immersed in the twilight glare; the magnitude tells how bright it is (under 6 = visible to the naked eye, under 0 = very bright). For example, Jupiter at magnitude -2.0 and an altitude of 15° to the east is unmistakable — it is the "bright point" that many mistake for an airplane. A concrete application: an app for amateur astronomers that shows "what's in the sky tonight" with direction and brightness.

```python
import libephemeris as ephem

# Rome, April 8, 2024
lat, lon = 41.9, 12.5
jd_noon = ephem.julday(2024, 4, 8, 12.0)

# Find the exact sunset time
jd_set, _ = ephem.rise_trans(
    jd_noon, ephem.SUN, lat, lon, rsmi=2
)
year, month, day, hour = ephem.revjul(jd_set)
hours_s = int(hour)
minutes_s = int((hour - hours_s) * 60)
print(f"Sunset: {hours_s:02d}:{minutes_s:02d} UT")
print()

# One hour after sunset
jd_evening = jd_set + 1.0 / 24.0
ephem.set_topo(lon, lat, 0)

pos_sun, _ = ephem.calc_ut(jd_evening, ephem.SUN, 0)

planets = [
    (ephem.MERCURY, "Mercury"),
    (ephem.VENUS,   "Venus"),
    (ephem.MARS,    "Mars"),
    (ephem.JUPITER, "Jupiter"),
    (ephem.SATURN,  "Saturn"),
]

dirs = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]

print("--- Planet Visibility (Rome, 8/4/2024) ---")
print()
for body_id, name in planets:
    pos, _ = ephem.calc_ut(jd_evening, body_id, 0)

    # Altitude and azimuth
    az, alt_t, alt_app = ephem.azalt(
        jd_evening, 0, (lon, lat, 0),
        1013.25, 15.0, (pos[0], pos[1], pos[2])
    )

    # Elongation from the Sun
    elong = abs(ephem.difdeg2n(pos[0], pos_sun[0]))

    # Apparent magnitude
    pheno, _ = ephem.pheno_ut(jd_evening, body_id, 0)
    mag = pheno[4]

    # Cardinal direction
    dir_idx = int((az + 22.5) / 45) % 8
    direction = dirs[dir_idx]

    visible = ("VISIBLE" if alt_app > 5 and elong > 15
               else "not visible")
    print(f"{name:10s}  alt: {alt_app:5.1f}°"
          f"  az: {az:5.1f}° ({direction:2s})"
          f"  mag: {mag:+.1f}"
          f"  elong: {elong:5.1f}°"
          f"  {visible}")
```

```
Sunset: 17:43 UT

--- Planet Visibility (Rome, 8/4/2024) ---

Mercury     alt:  -5.5°  az: 111.8° (E )  mag: +4.8  elong:   5.4°  not visible
Venus       alt: -25.8°  az: 116.4° (SE)  mag: -3.9  elong:  15.0°  not visible
Mars        alt: -44.6°  az: 129.0° (SE)  mag: +1.2  elong:  36.4°  not visible
Jupiter     alt:  15.5°  az:  98.8° (E )  mag: -2.0  elong:  29.6°  VISIBLE
Saturn      alt: -43.6°  az: 127.4° (SE)  mag: +1.1  elong:  35.0°  not visible
```

> On this evening in April 2024, only Jupiter is visible — low in the east, very bright (mag -2.0). The other planets are below the horizon or too close to the Sun.

---

## Recipe 9 — Aspects Between Planets

**Problem**: Calculate the angular distance between all planets and identify the active aspects (conjunction, sextile, square, trine, opposition).

**What it is used for**: Aspects are the angular relationships between planets — the fundamental language of astrology. A conjunction (0°) merges the energies of the two planets; a trine (120°) harmonizes them; a square (90°) puts them in tension. For example, Sun conjunct Moon (as on April 8, 2024) indicates a New Moon — a time of emotional reset and new beginnings. Mars conjunct Saturn indicates a phase of frustration but also of constructive discipline. Knowing whether an aspect is "applying" (forming) or "separating" (dissolving) changes the interpretation: an applying aspect is more intense. A concrete application: an automatic interpretation module that reads the aspects and generates a descriptive text.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

planets = [
    (ephem.SUN,     "Sun"),
    (ephem.MOON,    "Moon"),
    (ephem.MERCURY, "Mercury"),
    (ephem.VENUS,   "Venus"),
    (ephem.MARS,    "Mars"),
    (ephem.JUPITER, "Jupiter"),
    (ephem.SATURN,  "Saturn"),
]

# Calculate positions and velocities
positions = {}
for body_id, name in planets:
    pos, _ = ephem.calc_ut(jd, body_id, ephem.FLG_SPEED)
    positions[name] = (pos[0], pos[3])  # longitude, velocity

aspects = {0: "Conjunction", 60: "Sextile",
           90: "Square", 120: "Trine",
           180: "Opposition"}
max_orbs = {0: 8, 60: 5, 90: 6, 120: 6, 180: 8}

print("--- Aspects Between Planets (April 8, 2024) ---")
print()
names = list(positions.keys())
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        n1, n2 = names[i], names[j]
        lon1, vel1 = positions[n1]
        lon2, vel2 = positions[n2]
        diff = abs(ephem.difdeg2n(lon1, lon2))

        for asp_deg, asp_name in aspects.items():
            orb = abs(diff - asp_deg)
            if orb <= max_orbs[asp_deg]:
                print(f"{n1:10s} {asp_name:14s} {n2:10s}"
                      f"  orb: {orb:4.1f}°")
                break
```

```
--- Aspects Between Planets (April 8, 2024) ---

Sun        Conjunction    Moon        orb:  3.7°
Sun        Conjunction    Mercury     orb:  5.8°
Mars       Conjunction    Saturn      orb:  1.6°
Jupiter    Sextile        Saturn      orb:  4.6°
```

---

## Recipe 10 — Solar Return

**Problem**: Find the exact moment when the Sun returns to the same longitude as at birth, and calculate the chart for that instant.

**What it is used for**: The solar return is one of the most widely used predictive techniques in astrology. Every year, the Sun returns exactly to the degree it was at the moment of birth — that precise instant generates a "chart of the year" that is superimposed on the birth chart. The Ascendant of the solar return indicates the "tone" of the year: for example, an Ascendant in Virgo suggests a year dedicated to work, organization, and health. The position of the Moon in the return indicates the emotional climate of the next twelve months. A concrete application: a software that automatically generates the solar return for every birthday and compares it with the birth chart.

```python
import libephemeris as ephem

# Born March 15, 1990, 10:00 UT, Rome
jd_natal = ephem.julday(1990, 3, 15, 10.0)
pos_natal, _ = ephem.calc_ut(jd_natal, ephem.SUN, 0)
lon_natal = pos_natal[0]

# Find the solar return for 2024
jd_search = ephem.julday(2024, 3, 1, 0.0)
jd_return = ephem.solcross_ut(lon_natal, jd_search)

year, month, day, hour = ephem.revjul(jd_return)
hours = int(hour)
minutes = int((hour - hours) * 60)
seconds = int(((hour - hours) * 60 - minutes) * 60)

signs = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

print("--- Solar Return 2024 ---")
print()
print(f"Natal Sun: {lon_natal:.6f}°")
print(f"Return:    {day:2d}/{month:02d}/{year}"
      f"  {hours:02d}:{minutes:02d}:{seconds:02d} UT")
print()

bodies = [
    (ephem.SUN,     "Sun"),
    (ephem.MOON,    "Moon"),
    (ephem.MERCURY, "Mercury"),
    (ephem.VENUS,   "Venus"),
    (ephem.MARS,    "Mars"),
    (ephem.JUPITER, "Jupiter"),
    (ephem.SATURN,  "Saturn"),
]

for body_id, name in bodies:
    pos, _ = ephem.calc_ut(jd_return, body_id, 0)
    sign = signs[int(pos[0] / 30)]
    degrees = pos[0] % 30
    d, m, s, sf, si = ephem.split_deg(degrees, 0)
    print(f"{name:10s}  {d:2d}° {m:02d}'  {sign}")

# Placidus Houses for Rome
cusps, ascmc = ephem.houses(
    jd_return, 41.9, 12.5, ord("P")
)
print()
print(f"SR Ascendant:   {ascmc[0]:.2f}°"
      f" ({signs[int(ascmc[0] / 30)]})")
print(f"SR Midheaven:   {ascmc[1]:.2f}°"
      f" ({signs[int(ascmc[1] / 30)]})")
```

```
--- Solar Return 2024 ---

Natal Sun: 354.556761°
Return:    14/03/2024  15:50:28 UT

Sun         24° 33'  Psc
Moon        23° 25'  Tau
Mercury      8° 25'  Ari
Venus        3° 24'  Psc
Mars        23° 31'  Aqr
Jupiter     13° 48'  Tau
Saturn      11° 34'  Psc

SR Ascendant:   158.76° (Vir)
SR Midheaven:   64.82° (Gem)
```

---

## Recipe 11 — Planetary Hours

**Problem**: Calculate the planetary hours of the day, dividing the time between sunrise and sunset into 12 parts and assigning a planet to each according to the Chaldean sequence.

**What it is used for**: Planetary hours are one of the oldest systems of time organization, dating back to Babylon. The idea is that every hour of the day is "ruled" by a planet, following the Chaldean sequence (Saturn → Jupiter → Mars → Sun → Venus → Mercury → Moon, and then it repeats). The first diurnal hour is ruled by the planet of the day: Monday = Moon, Tuesday = Mars, etc. A planetary hour does not last 60 minutes — it lasts 1/12 of the daylight period (in summer, diurnal hours are longer). In horary astrology and ceremonial magic, one chooses the hour of the appropriate planet to start an activity: the hour of Venus for matters of love, the hour of Mercury to sign contracts, the hour of Jupiter for financial affairs. A concrete application: an app that shows in real-time which planetary hour is currently active, with a countdown to the next one.

```python
import libephemeris as ephem

lat, lon = 41.9, 12.5  # Rome
jd_noon = ephem.julday(2024, 4, 8, 12.0)

# Find sunrise, sunset, and next sunrise
jd_rise, _ = ephem.rise_trans(
    jd_noon - 0.5, ephem.SUN, lat, lon, rsmi=1
)
jd_set, _ = ephem.rise_trans(
    jd_noon, ephem.SUN, lat, lon, rsmi=2
)
jd_rise_tmrw, _ = ephem.rise_trans(
    jd_set, ephem.SUN, lat, lon, rsmi=1
)

_, _, _, h_rise = ephem.revjul(jd_rise)
_, _, _, h_set = ephem.revjul(jd_set)
print(f"Sunrise:  {int(h_rise):02d}:"
      f"{int((h_rise - int(h_rise)) * 60):02d} UT")
print(f"Sunset:   {int(h_set):02d}:"
      f"{int((h_set - int(h_set)) * 60):02d} UT")

# Duration of a planetary hour
day_hour = (jd_set - jd_rise) / 12
night_hour = (jd_rise_tmrw - jd_set) / 12
print(f"Day hour duration:   {day_hour * 24 * 60:.1f} min")
print(f"Night hour duration: {night_hour * 24 * 60:.1f} min")
print()

# Chaldean sequence
chaldean = ["Saturn", "Jupiter", "Mars", "Sun",
            "Venus", "Mercury", "Moon"]
# April 8, 2024 = Monday -> Moon (index 6)
ruler = 6

print("--- Planetary Hours (Monday 8/4/2024, Rome) ---")
print()
print("Day hours:")
idx = ruler
for i in range(12):
    start = jd_rise + i * day_hour
    end = jd_rise + (i + 1) * day_hour
    _, _, _, hi = ephem.revjul(start)
    _, _, _, hf = ephem.revjul(end)
    print(f"  Hour {i+1:2d}  "
          f"{int(hi):02d}:{int((hi-int(hi))*60):02d}-"
          f"{int(hf):02d}:{int((hf-int(hf))*60):02d} UT"
          f"  {chaldean[idx % 7]}")
    idx += 1

print()
print("Night hours:")
for i in range(12):
    start = jd_set + i * night_hour
    end = jd_set + (i + 1) * night_hour
    _, _, _, hi = ephem.revjul(start)
    _, _, _, hf = ephem.revjul(end)
    print(f"  Hour {i+1:2d}  "
          f"{int(hi):02d}:{int((hi-int(hi))*60):02d}-"
          f"{int(hf):02d}:{int((hf-int(hf))*60):02d} UT"
          f"  {chaldean[idx % 7]}")
    idx += 1
```

```
Sunrise:  04:40 UT
Sunset:   17:43 UT
Day hour duration:   65.2 min
Night hour duration: 54.6 min

--- Planetary Hours (Monday 8/4/2024, Rome) ---

Day hours:
  Hour  1  04:40-05:45 UT  Moon
  Hour  2  05:45-06:51 UT  Saturn
  Hour  3  06:51-07:56 UT  Jupiter
  Hour  4  07:56-09:01 UT  Mars
  Hour  5  09:01-10:06 UT  Sun
  Hour  6  10:06-11:12 UT  Venus
  Hour  7  11:12-12:17 UT  Mercury
  Hour  8  12:17-13:22 UT  Moon
  Hour  9  13:22-14:27 UT  Saturn
  Hour 10  14:27-15:33 UT  Jupiter
  Hour 11  15:33-16:38 UT  Mars
  Hour 12  16:38-17:43 UT  Sun

Night hours:
  Hour  1  17:43-18:38 UT  Venus
  Hour  2  18:38-19:32 UT  Mercury
  Hour  3  19:32-20:27 UT  Moon
  Hour  4  20:27-21:22 UT  Saturn
  Hour  5  21:22-22:16 UT  Jupiter
  Hour  6  22:16-23:11 UT  Mars
  Hour  7  23:11-00:05 UT  Sun
  Hour  8  00:05-01:00 UT  Venus
  Hour  9  01:00-01:55 UT  Mercury
  Hour 10  01:55-02:49 UT  Moon
  Hour 11  02:49-03:44 UT  Saturn
  Hour 12  03:44-04:38 UT  Jupiter
```

> Note how the diurnal hours (65 min) are longer than the nocturnal ones (55 min) — in April in the Northern Hemisphere, the day is longer than the night.

---

## Recipe 12 — Printable Monthly Ephemeris

**Problem**: Generate a table with the position of all planets for each day of the month, in zodiacal format.

**What it is used for**: The ephemeris is the astrologer's working tool — a table showing "where the planets are" every day at noon. Before computers, astrologers consulted voluminous ephemeris books (such as the famous *Raphael's Ephemeris*) to calculate charts manually. Even today, a compact table is useful to get an overview of the month: you can see at a glance when a planet changes sign, when the Moon quickly transits through the zodiac, and when Mercury slows down before its retrograde. A concrete application: a generator of customized PDF ephemerides for any given year.

```python
import libephemeris as ephem

signs = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

def fmt(lon):
    """Format longitude in degrees + abbreviated sign."""
    s = signs[int(lon / 30)]
    d = lon % 30
    return f"{d:5.1f}{s}"

bodies = [
    (ephem.SUN,     "Sun"),
    (ephem.MOON,    "Moon"),
    (ephem.MERCURY, "Merc"),
    (ephem.VENUS,   "Venu"),
    (ephem.MARS,    "Mars"),
    (ephem.JUPITER, "Jupi"),
    (ephem.SATURN,  "Satu"),
]

print("--- Ephemeris April 2024 (Noon UT) ---")
print()
header = "Date     "
for _, name in bodies:
    header += f"  {name:>8s}"
print(header)
print("-" * len(header))

for day in range(1, 31):
    jd = ephem.julday(2024, 4, day, 12.0)
    row = f"{day:2d}/04/2024"
    for body_id, _ in bodies:
        pos, _ = ephem.calc_ut(jd, body_id, 0)
        row += f"  {fmt(pos[0]):>8s}"
    print(row)
```

```
--- Ephemeris April 2024 (Noon UT) ---

Date           Sun      Moon      Merc      Venu      Mars      Jupi      Satu
------------------------------------------------------------------------------
 1/04/2024   12.2Ari    4.4Cap   27.2Ari   25.5Psc    7.4Psc   17.5Tau   13.7Psc
 2/04/2024   13.2Ari   17.8Cap   27.2Ari   26.7Psc    8.2Psc   17.7Tau   13.8Psc
 3/04/2024   14.2Ari    1.7Aqr   27.1Ari   27.9Psc    9.0Psc   17.9Tau   13.9Psc
 4/04/2024   15.2Ari   15.9Aqr   26.8Ari   29.2Psc    9.7Psc   18.1Tau   14.0Psc
 5/04/2024   16.2Ari    0.5Psc   26.5Ari    0.4Ari   10.5Psc   18.3Tau   14.1Psc
 6/04/2024   17.2Ari   15.3Psc   26.1Ari    1.6Ari   11.3Psc   18.5Tau   14.2Psc
 7/04/2024   18.2Ari    0.4Ari   25.6Ari    2.9Ari   12.1Psc   18.8Tau   14.3Psc
 8/04/2024   19.1Ari   15.4Ari   25.0Ari    4.1Ari   12.8Psc   19.0Tau   14.4Psc
 9/04/2024   20.1Ari    0.4Tau   24.3Ari    5.4Ari   13.6Psc   19.2Tau   14.5Psc
10/04/2024   21.1Ari   15.1Tau   23.6Ari    6.6Ari   14.4Psc   19.4Tau   14.6Psc
11/04/2024   22.1Ari   29.4Tau   22.9Ari    7.8Ari   15.2Psc   19.7Tau   14.7Psc
12/04/2024   23.1Ari   13.4Gem   22.1Ari    9.1Ari   16.0Psc   19.9Tau   14.9Psc
13/04/2024   24.0Ari   26.8Gem   21.4Ari   10.3Ari   16.7Psc   20.1Tau   15.0Psc
14/04/2024   25.0Ari    9.9Cnc   20.6Ari   11.5Ari   17.5Psc   20.3Tau   15.1Psc
15/04/2024   26.0Ari   22.6Cnc   19.9Ari   12.8Ari   18.3Psc   20.6Tau   15.2Psc
16/04/2024   27.0Ari    4.9Leo   19.2Ari   14.0Ari   19.1Psc   20.8Tau   15.3Psc
17/04/2024   28.0Ari   17.0Leo   18.6Ari   15.2Ari   19.8Psc   21.0Tau   15.4Psc
18/04/2024   28.9Ari   28.9Leo   18.0Ari   16.5Ari   20.6Psc   21.2Tau   15.5Psc
19/04/2024   29.9Ari   10.8Vir   17.5Ari   17.7Ari   21.4Psc   21.5Tau   15.6Psc
20/04/2024    0.9Tau   22.5Vir   17.0Ari   18.9Ari   22.2Psc   21.7Tau   15.7Psc
21/04/2024    1.9Tau    4.4Lib   16.7Ari   20.2Ari   22.9Psc   21.9Tau   15.8Psc
22/04/2024    2.8Tau   16.3Lib   16.4Ari   21.4Ari   23.7Psc   22.2Tau   15.9Psc
23/04/2024    3.8Tau   28.3Lib   16.2Ari   22.6Ari   24.5Psc   22.4Tau   16.0Psc
24/04/2024    4.8Tau   10.5Sco   16.0Ari   23.9Ari   25.3Psc   22.6Tau   16.0Psc
25/04/2024    5.8Tau   22.9Sco   16.0Ari   25.1Ari   26.0Psc   22.8Tau   16.1Psc
26/04/2024    6.7Tau    5.5Sgr   16.0Ari   26.3Ari   26.8Psc   23.1Tau   16.2Psc
27/04/2024    7.7Tau   18.3Sgr   16.1Ari   27.6Ari   27.6Psc   23.3Tau   16.3Psc
28/04/2024    8.7Tau    1.3Cap   16.3Ari   28.8Ari   28.3Psc   23.5Tau   16.4Psc
29/04/2024    9.7Tau   14.6Cap   16.6Ari    0.0Tau   29.1Psc   23.8Tau   16.5Psc
30/04/2024   10.6Tau   28.1Cap   17.0Ari    1.3Tau   29.9Psc   24.0Tau   16.6Psc
```

> From the table it can be clearly read: Mercury moves backward from 27.2° Ari to 16.0° Ari (direct station on 25/04), then starts advancing again. The Moon traverses the entire zodiac in a month. The Sun advances about 1° per day through Aries and Taurus.

---

## Summary

This chapter has presented 12 complete recipes for the most common astronomical and astrological calculations:

- **Recipe 1** — Complete birth chart with planets, houses, and angles
- **Recipe 2** — Daily transits on a birth chart
- **Recipe 3** — Sun ingresses into the 12 signs (equinoxes and solstices)
- **Recipe 4** — Mercury retrogrades (stations and verification)
- **Recipe 5** — New and Full Moons of the year
- **Recipe 6** — Next solar eclipses visible from a city
- **Recipe 7** — Vedic sidereal chart with Nakshatras and Padas
- **Recipe 8** — Visibility of planets tonight (altitude, magnitude, direction)
- **Recipe 9** — Aspects between planets with orb
- **Recipe 10** — Solar Return
- **Recipe 11** — Planetary hours (Chaldean sequence)
- **Recipe 12** — Printable monthly ephemeris

Each recipe is self-contained and ready to use: just copy the code, modify the date and coordinates, and run it.

## Functions Used in the Recipes

- `julday(y, m, d, h)` / `revjul(jd)` — date ↔ Julian Day conversion
- `calc_ut(jd, body, flag)` — position of a celestial body
- `houses(jd, lat, lon, hsys)` — house cusps and angles
- `split_deg(deg, flag)` — splits decimal degrees into °, ', "
- `difdeg2n(p1, p2)` — normalized angular difference [-180, 180]
- `solcross_ut(x, jd)` — Sun crossing at a longitude
- `mooncross_ut(x, jd)` — Moon crossing at a longitude
- `find_station_ut(body, jd, type)` — next retrograde/direct station (from `libephemeris.crossing`)
- `is_retrograde(body, jd)` — retrograde verification (from `libephemeris.crossing`)
- `sol_eclipse_when_loc(jd, lat, lon, alt)` — local solar eclipse
- `set_sid_mode(mode)` / `get_ayanamsa_ut(jd)` — sidereal zodiac
- `rise_trans(jd, body, lat, lon, rsmi=...)` — sunrise and sunset
- `azalt(jd, flag, geopos, press, temp, xin)` — horizontal coordinates
- `pheno_ut(jd, body, flag)` — phenomena (magnitude, phase, elongation)
- `set_topo(lon, lat, alt)` — observer's position