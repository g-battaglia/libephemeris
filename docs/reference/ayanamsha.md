# Ayanamsha Definitions and Reference Systems

This document provides comprehensive documentation of all 47 predefined ayanamsha (sidereal zodiac) modes implemented in LibEphemeris. Each section explains the astronomical basis, historical origin, reference point, and the epoch when the ayanamsha was zero.

## Table of Contents

1. [Introduction to Ayanamsha](#introduction-to-ayanamsha)
2. [Western Sidereal Traditions](#western-sidereal-traditions)
3. [Indian Astrological Systems](#indian-astrological-systems)
4. [Babylonian Systems](#babylonian-systems)
5. [Star-Based (True) Ayanamshas](#star-based-true-ayanamshas)
6. [Galactic Alignment Systems](#galactic-alignment-systems)
7. [Historical Epoch Systems](#historical-epoch-systems)
8. [User-Defined Ayanamsha](#user-defined-ayanamsha)
9. [Technical Implementation](#technical-implementation)
10. [References](#references)

---

## Introduction to Ayanamsha

### What is Ayanamsha?

**Ayanamsha** (Sanskrit: ayanamsha, "portion of movement") is the angular difference between the tropical zodiac and the sidereal zodiac, measured in degrees. It represents the cumulative effect of the precession of the equinoxes since the two zodiacs were aligned.

```
Ayanamsha = Tropical Longitude - Sidereal Longitude
Sidereal Longitude = Tropical Longitude - Ayanamsha
```

### The Two Zodiacs

| Zodiac | Reference | Basis | Changes Over Time |
|--------|-----------|-------|-------------------|
| **Tropical** | Vernal Equinox (0° Aries) | Earth's seasons | Fixed to equinoxes |
| **Sidereal** | Fixed stars | Stellar positions | Fixed to stars |

### Precession of the Equinoxes

The Earth's axis wobbles like a spinning top over a ~25,772 year cycle (Platonic Year). This causes the vernal equinox to precess westward through the zodiac at approximately:

- **50.29 arcseconds per year** (modern value)
- **1.397° per century** (~5027.8 arcsec/century)

Because different traditions chose different starting points and epochs for when the zodiacs were aligned, there are numerous ayanamsha systems in use today.

---

## Western Sidereal Traditions

### SIDM_FAGAN_BRADLEY (0)

**Fagan-Bradley / Synetic Vernal Point**

| Property | Value |
|----------|-------|
| Constant | `SIDM_FAGAN_BRADLEY = 0` |
| J2000 Value | 24.7403° |
| Zero Epoch | 221 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
Cyril Fagan (1896-1970), an Irish astrologer, and Donald Bradley (1925-1974), an American statistician, developed this ayanamsha based on statistical analysis of horoscopes. They determined that the tropical and sidereal zodiacs coincided around 221 CE.

**Reference Point:**
The "Synetic Vernal Point" - determined empirically through statistical analysis rather than from stellar positions. Fagan believed this represented the true fiducial (fixed point) of the zodiac as used in ancient Babylon.

**Usage:**
Standard for Western sidereal astrology. Used by the American Federation of Astrologers (sidereal division).

---

### SIDM_DELUCE (2)

**De Luce**

| Property | Value |
|----------|-------|
| Constant | `SIDM_DELUCE = 2` |
| J2000 Value | 27.8158° |
| Zero Epoch | ~4 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
Robert De Luce (1877-1964) was an American astrologer who calculated this ayanamsha based on the star Spica being at 29° Virgo in the sidereal zodiac.

**Reference Point:**
Spica (Alpha Virginis) positioned at 29° Virgo.

---

### SIDM_DJWHAL_KHUL (6)

**Djwhal Khul / Alice Bailey**

| Property | Value |
|----------|-------|
| Constant | `SIDM_DJWHAL_KHUL = 6` |
| J2000 Value | 28.3597° |
| Zero Epoch | ~62 BCE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
Based on the esoteric teachings transmitted through Alice Bailey (1880-1949), attributed to the Tibetan master Djwhal Khul. This system places emphasis on the Age of Aquarius beginning around 2117 CE.

**Reference Point:**
Esoteric tradition placing the vernal point at a specific position relative to constellational boundaries.

---

## Indian Astrological Systems

### SIDM_LAHIRI (1)

**Lahiri / Chitrapaksha**

| Property | Value |
|----------|-------|
| Constant | `SIDM_LAHIRI = 1` |
| J2000 Value | 23.8571° |
| Zero Epoch | 285 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
Named after Dr. N.C. Lahiri (1906-1980), who chaired the Calendar Reform Committee appointed by the Government of India in 1952. This committee standardized the ayanamsha for the Indian national calendar (Rashtriya Panchang) in 1956.

**Reference Point:**
The star Spica (Chitra in Sanskrit, Alpha Virginis) is fixed at exactly 180° (0° Libra). "Chitrapaksha" means "Spica-wing" or "the side of Spica."

**Official Status:**
The official ayanamsha of the Government of India. Values are published annually in the Indian Astronomical Ephemeris (IAE).

**Zero Date:**
The tropical and sidereal zodiacs coincided around 285 CE, when the vernal equinox was at 0° Aries sidereal.

---

### SIDM_RAMAN (3)

**B.V. Raman**

| Property | Value |
|----------|-------|
| Constant | `SIDM_RAMAN = 3` |
| J2000 Value | 22.4108° |
| Zero Epoch | 389 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
B.V. Raman (1912-1998) was one of the most influential Indian astrologers of the 20th century. He calculated his own ayanamsha based on classical Sanskrit astronomical texts, particularly the Suryasiddhanta.

**Reference Point:**
Based on traditional Indian astronomical principles with adjustments for observed stellar positions.

---

### SIDM_USHASHASHI (4)

**Ushashashi**

| Property | Value |
|----------|-------|
| Constant | `SIDM_USHASHASHI = 4` |
| J2000 Value | 20.0575° |
| Zero Epoch | 559 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
Based on the traditional Indian calculation where the ayanamsha was 0° in 559 CE.

---

### SIDM_KRISHNAMURTI (5)

**K.S. Krishnamurti (KP System)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_KRISHNAMURTI = 5` |
| J2000 Value | 23.7602° |
| Zero Epoch | 291 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
K.S. Krishnamurti (1908-1972) developed the Krishnamurti Paddhati (KP) system, a simplified and systematized approach to Vedic astrology. His ayanamsha is a slight modification of Lahiri's, based on his own research.

**Reference Point:**
Similar to Lahiri but with minor adjustments based on Krishnamurti's observations and calculations.

---

### SIDM_YUKTESHWAR (7)

**Sri Yukteshwar**

| Property | Value |
|----------|-------|
| Constant | `SIDM_YUKTESHWAR = 7` |
| J2000 Value | 22.4788° |
| Zero Epoch | 499 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
Sri Yukteshwar Giri (1855-1936), the guru of Paramahansa Yogananda, presented his cosmological and astrological views in "The Holy Science" (1894). He proposed a unique theory of yuga cycles tied to the precession of the equinoxes.

**Reference Point:**
Based on Yukteshwar's calculation that the autumnal equinox coincided with the first point of Libra in 499 CE.

---

### SIDM_JN_BHASIN (8)

**J.N. Bhasin**

| Property | Value |
|----------|-------|
| Constant | `SIDM_JN_BHASIN = 8` |
| J2000 Value | 22.7621° |
| Zero Epoch | 364 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
J.N. Bhasin was an Indian astrologer who proposed this ayanamsha value based on his research into classical Indian astronomical texts.

---

### SIDM_SURYASIDDHANTA (21)

**Suryasiddhanta**

| Property | Value |
|----------|-------|
| Constant | `SIDM_SURYASIDDHANTA = 21` |
| J2000 Value | 20.8951° |
| Zero Epoch | 499 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
The Suryasiddhanta ("Sun Treatise") is one of the most important ancient Indian astronomical texts, dated to approximately 400 CE. It contains detailed instructions for calculating planetary positions.

**Reference Point:**
Based on the traditional calculation in the Suryasiddhanta text, where the tropical and sidereal zodiacs were aligned in 499 CE (during the Kaliyuga).

---

### SIDM_SURYASIDDHANTA_MSUN (22)

**Suryasiddhanta (Mean Sun)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_SURYASIDDHANTA_MSUN = 22` |
| J2000 Value | 20.6804° |
| Zero Epoch | ~514 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
A variant of the Suryasiddhanta calculation using the mean Sun position rather than the true Sun.

---

### SIDM_ARYABHATA (23)

**Aryabhata**

| Property | Value |
|----------|-------|
| Constant | `SIDM_ARYABHATA = 23` |
| J2000 Value | 20.8951° |
| Zero Epoch | 499 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
Aryabhata (476-550 CE) was one of the greatest mathematicians and astronomers of ancient India. His work "Aryabhatiya" (499 CE) contains important astronomical calculations.

**Reference Point:**
The epoch of the Aryabhatiya is midnight at the start of the Kali Yuga (3102 BCE) but calculated from 499 CE observations.

---

### SIDM_ARYABHATA_MSUN (24)

**Aryabhata (Mean Sun)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_ARYABHATA_MSUN = 24` |
| J2000 Value | 20.6574° |
| Zero Epoch | ~516 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
A variant using mean Sun position.

---

### SIDM_ARYABHATA_522 (37)

**Aryabhata 522**

| Property | Value |
|----------|-------|
| Constant | `SIDM_ARYABHATA_522 = 37` |
| J2000 Value | 20.5758° |
| Zero Epoch | 522 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
An alternative Aryabhata calculation with the zero epoch at 522 CE.

---

### SIDM_SS_REVATI (25)

**Suryasiddhanta Revati**

| Property | Value |
|----------|-------|
| Constant | `SIDM_SS_REVATI = 25` |
| J2000 Value | 20.1034° |
| Zero Epoch | 556 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
Based on the Suryasiddhanta with the star Revati (Zeta Piscium) as the reference point for 0° Aries.

**Reference Point:**
The star Revati (Zeta Piscium) is placed at 359°50' (or 29°50' Pisces), making it the defining star for the sidereal first point of Aries.

---

### SIDM_SS_CITRA (26)

**Suryasiddhanta Citra**

| Property | Value |
|----------|-------|
| Constant | `SIDM_SS_CITRA = 26` |
| J2000 Value | 23.0058° |
| Zero Epoch | 349 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
Based on the Suryasiddhanta with Spica (Citra) as the reference star at 180°.

---

## Babylonian Systems

### SIDM_BABYL_KUGLER1 (9)

**Babylonian (Kugler 1)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_BABYL_KUGLER1 = 9` |
| J2000 Value | 23.5336° |
| Zero Epoch | ~304 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
Franz Xaver Kugler (1862-1929) was a German Jesuit priest and Assyriologist who studied ancient Babylonian astronomical tablets. This is his first reconstruction of the Babylonian zodiac.

**Reference Point:**
Based on Kugler's analysis of Babylonian "goal-year" texts and astronomical diaries.

---

### SIDM_BABYL_KUGLER2 (10)

**Babylonian (Kugler 2)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_BABYL_KUGLER2 = 10` |
| J2000 Value | 24.9336° |
| Zero Epoch | ~204 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
Kugler's second reconstruction, differing in the placement of certain reference stars.

---

### SIDM_BABYL_KUGLER3 (11)

**Babylonian (Kugler 3)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_BABYL_KUGLER3 = 11` |
| J2000 Value | 25.7836° |
| Zero Epoch | ~144 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
Kugler's third variant.

---

### SIDM_BABYL_HUBER (12)

**Babylonian (Huber)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_BABYL_HUBER = 12` |
| J2000 Value | 24.7336° |
| Zero Epoch | ~220 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
Peter Huber's reconstruction of the Babylonian zodiac, based on detailed analysis of cuneiform tablets.

---

### SIDM_BABYL_ETPSC (13)

**Babylonian (Eta Piscium)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_BABYL_ETPSC = 13` |
| J2000 Value | 24.5225° |
| Zero Epoch | ~236 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
Based on the star Eta Piscium as a zodiacal reference point in Babylonian astronomy.

---

### SIDM_BABYL_BRITTON (38)

**Babylonian (Britton)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_BABYL_BRITTON = 38` |
| J2000 Value | 24.6158° |
| Zero Epoch | ~228 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
John Britton's (1939-2010) reconstruction based on his extensive study of Babylonian mathematical astronomy. Britton was a leading scholar of Babylonian astronomy at Yale University.

---

## Star-Based (True) Ayanamshas

These ayanamshas are calculated from the actual current position of reference stars, accounting for proper motion, precession, and nutation. They are recalculated for each date rather than using a fixed formula.

### SIDM_ALDEBARAN_15TAU (14)

**Aldebaran at 15° Taurus**

| Property | Value |
|----------|-------|
| Constant | `SIDM_ALDEBARAN_15TAU = 14` |
| J2000 Value | ~24.76° |
| Zero Epoch | ~223 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
This system places the bright star Aldebaran (Alpha Tauri, the "Eye of the Bull") at exactly 15° Taurus (45° ecliptic longitude) in the sidereal zodiac.

**Reference Point:**
- **Star:** Aldebaran (Alpha Tauri)
- **Hipparcos ID:** HIP 21421
- **Visual Magnitude:** 0.85 (13th brightest star)
- **Sidereal Position:** Fixed at 15°00' Taurus

---

### SIDM_TRUE_CITRA (27)

**True Citra / True Spica**

| Property | Value |
|----------|-------|
| Constant | `SIDM_TRUE_CITRA = 27` |
| J2000 Value | ~23.86° |
| Zero Epoch | ~285 CE |
| Calculation | Dynamic (star position) |

**Historical Basis:**
This is a dynamically calculated version of the Lahiri ayanamsha. Instead of using a fixed J2000 value and precession rate, it calculates the actual current position of Spica and places it at exactly 180°.

**Reference Point:**
- **Star:** Spica (Alpha Virginis) / Citra in Sanskrit
- **Hipparcos ID:** HIP 65474
- **Position:** Tropical longitude calculated at observation date
- **Sidereal Position:** Fixed at exactly 180° (0° Libra)

**Astronomical Data (Hipparcos):**
```
Right Ascension (J2000): 13h 25m 11.579s
Declination (J2000): -11° 09' 40.75"
Proper Motion (RA): -42.50 mas/yr
Proper Motion (Dec): -31.73 mas/yr
Parallax: 13.06 mas (distance: 250 ly)
```

---

### SIDM_TRUE_REVATI (28)

**True Revati**

| Property | Value |
|----------|-------|
| Constant | `SIDM_TRUE_REVATI = 28` |
| J2000 Value | ~20.05° |
| Zero Epoch | ~556 CE |
| Calculation | Dynamic (star position) |

**Historical Basis:**
Calculated from the actual position of the star Revati (Zeta Piscium), which traditionally marks the end of Pisces and beginning of Aries in Indian astronomy.

**Reference Point:**
- **Star:** Zeta Piscium (Revati)
- **Hipparcos ID:** HIP 5737
- **Sidereal Position:** Fixed at 359°50' (29°50' Pisces)

**Astronomical Data (Gaia DR3):**
```
Right Ascension (J2000): 01h 13m 43.883s
Declination (J2000): +07° 34' 31.274"
Proper Motion (RA): 109.57 mas/yr
Proper Motion (Dec): -16.20 mas/yr
Parallax: 21.85 mas (distance: 148 ly)
```

---

### SIDM_TRUE_PUSHYA (29)

**True Pushya**

| Property | Value |
|----------|-------|
| Constant | `SIDM_TRUE_PUSHYA = 29` |
| J2000 Value | ~22.72° |
| Zero Epoch | ~373 CE |
| Calculation | Dynamic (star position) |

**Historical Basis:**
Based on the nakshatra Pushya (Delta Cancri), which is significant in Vedic astrology as one of the most auspicious lunar mansions.

**Reference Point:**
- **Star:** Delta Cancri (Asellus Australis / Pushya)
- **Hipparcos ID:** HIP 42911
- **Sidereal Position:** Fixed at 106° (16° Cancer)

**Astronomical Data (Gaia DR3):**
```
Right Ascension (J2000): 08h 44m 41.099s
Declination (J2000): +18° 09' 15.509"
```

---

### SIDM_TRUE_MULA (35)

**True Mula**

| Property | Value |
|----------|-------|
| Constant | `SIDM_TRUE_MULA = 35` |
| J2000 Value | ~24.59° |
| Zero Epoch | ~240 CE |
| Calculation | Dynamic (star position) |

**Historical Basis:**
Based on the nakshatra Mula (Lambda Scorpii, Shaula), located at the "root" of Sagittarius near the Galactic Center.

**Reference Point:**
- **Star:** Lambda Scorpii (Shaula / Mula)
- **Hipparcos ID:** HIP 85927
- **Sidereal Position:** Fixed at 240° (0° Sagittarius)

**Astronomical Data (Hipparcos 2007):**
```
Right Ascension (J2000): 17h 33m 36.520s
Declination (J2000): -37° 06' 13.764"
Proper Motion (RA): -8.90 mas/yr
Proper Motion (Dec): -29.95 mas/yr
```

---

### SIDM_TRUE_SHEORAN (39)

**True Sheoran**

| Property | Value |
|----------|-------|
| Constant | `SIDM_TRUE_SHEORAN = 39` |
| J2000 Value | ~25.23° |
| Zero Epoch | ~186 CE |
| Calculation | Dynamic (star position) |

**Historical Basis:**
Named after its proponent, this ayanamsha uses Spica with a specific offset different from True Citra.

**Reference Point:**
- **Star:** Spica (Alpha Virginis)
- **Offset:** 178.607° from Spica's tropical longitude

---

### SIDM_VALENS_MOON (42)

**Valens Moon**

| Property | Value |
|----------|-------|
| Constant | `SIDM_VALENS_MOON = 42` |
| J2000 Value | ~22.80° |
| Zero Epoch | ~361 CE |
| Calculation | Dynamic (star position) |

**Historical Basis:**
Based on the work of Vettius Valens (120-175 CE), a Hellenistic astrologer from Antioch whose "Anthology" is the longest surviving astrological text from antiquity.

**Reference Point:**
- **Star:** Spica (Alpha Virginis)
- **Offset:** 181.0458° from Spica's tropical longitude

---

## Galactic Alignment Systems

These ayanamshas are based on the position of the Galactic Center (Sagittarius A*) or the Galactic Equator.

### SIDM_GALCENT_0SAG (17)

**Galactic Center at 0° Sagittarius**

| Property | Value |
|----------|-------|
| Constant | `SIDM_GALCENT_0SAG = 17` |
| J2000 Value | ~26.85° |
| Zero Epoch | ~69 CE |
| Calculation | Dynamic (GC position) |

**Historical Basis:**
This system places the Galactic Center (Sagittarius A*, the supermassive black hole at the center of the Milky Way) at exactly 0° Sagittarius (240° ecliptic longitude).

**Reference Point:**
- **Object:** Sagittarius A* (Galactic Center)
- **Position (Reid & Brunthaler 2004):**
```
Right Ascension (J2000): 17h 45m 40.0409s
Declination (J2000): -29° 00' 28.118"
```

---

### SIDM_GALCENT_RGILBRAND (30)

**Galactic Center (Gil Brand)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_GALCENT_RGILBRAND = 30` |
| J2000 Value | ~22.47° |
| Zero Epoch | ~381 CE |
| Calculation | Dynamic (GC position) |

**Historical Basis:**
Rafael Gil Brand's variant of the Galactic Center ayanamsha, with the Galactic Center positioned at approximately 4°23' Sagittarius rather than 0°.

---

### SIDM_GALCENT_MULA_WILHELM (36)

**Galactic Center at Mula (Wilhelm)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_GALCENT_MULA_WILHELM = 36` |
| J2000 Value | ~20.04° |
| Zero Epoch | ~558 CE |
| Calculation | Dynamic (GC position) |

**Historical Basis:**
Ernst Wilhelm's variant placing the Galactic Center at the nakshatra Mula (approximately 6°49' Sagittarius).

---

### SIDM_GALCENT_COCHRANE (40)

**Galactic Center (Cochrane)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_GALCENT_COCHRANE = 40` |
| J2000 Value | ~-3.15° |
| Zero Epoch | ~2231 CE (future) |
| Calculation | Dynamic (GC position) |

**Historical Basis:**
David Cochrane's variant placing the Galactic Center at 0° Capricorn (270°), based on his tropical/sidereal synthesis research.

---

### SIDM_GALEQU_IAU1958 (31)

**Galactic Equator (IAU 1958)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_GALEQU_IAU1958 = 31` |
| J2000 Value | ~30.11° |
| Zero Epoch | ~165 BCE |
| Calculation | Dynamic (galactic pole) |

**Historical Basis:**
Based on the IAU 1958 definition of the galactic coordinate system. The ascending node of the galactic equator on the ecliptic is used as a reference.

**Reference Point:**
- **Object:** Galactic North Pole (IAU 1958)
- **Galactic coordinates defined relative to B1950 epoch**
- **Node on ecliptic placed at 0° Sagittarius**

---

### SIDM_GALEQU_TRUE (32)

**Galactic Equator (True)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_GALEQU_TRUE = 32` |
| J2000 Value | ~30.11° |
| Zero Epoch | ~165 BCE |
| Calculation | Dynamic (galactic pole) |

**Historical Basis:**
Uses the true (current) position of the Galactic North Pole to determine the galactic equator intersection with the ecliptic.

---

### SIDM_GALEQU_MULA (33)

**Galactic Equator at Mula**

| Property | Value |
|----------|-------|
| Constant | `SIDM_GALEQU_MULA = 33` |
| J2000 Value | ~23.51° |
| Zero Epoch | ~310 CE |
| Calculation | Dynamic (galactic pole) |

**Historical Basis:**
A variant placing the galactic equator node at the nakshatra Mula.

---

### SIDM_GALALIGN_MARDYKS (34)

**Galactic Alignment (Mardyks)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_GALALIGN_MARDYKS = 34` |
| J2000 Value | 30.0178° |
| Reference Epoch | JD 2451079.771 (September 1998), ayanamsha exactly 30° |
| Calculation | Fixed-epoch frame mode |

**Historical Basis:**
Raymond Mardyks' system based on the alignment of the winter solstice with the Galactic Center, which he associates with the end of the Mayan calendar cycle in 2012.

**Frame behavior:** like the J2000/J1900/B1950 modes, sidereal *positions*
under this mode are frame projections onto the mean ecliptic and equinox of
the reference epoch (with a constant 30° longitude offset on the ecliptic/XYZ
channels; the equatorial channel carries no offset), not a scalar ayanamsha
subtraction. `get_ayanamsa()` still returns the scalar value above. See
`docs/comparison/known-differences.md` §10.2b.

---

### SIDM_GALEQU_FIORENZA (41)

**Galactic Equator (Fiorenza)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_GALEQU_FIORENZA = 41` |
| J2000 Value | 25.0000° |
| Zero Epoch | ~214 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
Nick Anthony Fiorenza's variant of the galactic equator-based ayanamsha, using a fixed formula rather than dynamic calculation.

---

## Historical Epoch Systems

### SIDM_HIPPARCHOS (15)

**Hipparchos**

| Property | Value |
|----------|-------|
| Constant | `SIDM_HIPPARCHOS = 15` |
| J2000 Value | 20.2478° |
| Zero Epoch | ~550 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
Named after Hipparchus of Nicaea (c. 190-120 BCE), the greatest astronomer of antiquity who discovered the precession of the equinoxes. Note: The zero epoch (~550 CE) is a conventional astronomical definition and differs from Hipparchus's historical era. The system is named in his honor as the discoverer of precession.

**Reference Point:**
Based on a conventional reconstruction of an ancient zodiac reference frame (see Meeus, "Astronomical Algorithms", Ch. 27).

---

### SIDM_SASSANIAN (16)

**Sassanian**

| Property | Value |
|----------|-------|
| Constant | `SIDM_SASSANIAN = 16` |
| J2000 Value | 19.9930° |
| Zero Epoch | 564 CE |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
Based on the astronomical tables of the Sassanian Empire (224-651 CE) in Persia. The Sassanians developed sophisticated astronomical traditions that later influenced both Islamic and Indian astronomy.

---

### SIDM_J2000 (18)

**J2000.0 (No Ayanamsha)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_J2000 = 18` |
| J2000 Value | 0.0° |
| Zero Epoch | 2000-01-01 12:00 TT |
| Calculation | Dynamic (precession only) |

**Historical Basis:**
This is not a traditional ayanamsha but rather places the sidereal zodiac to coincide with the tropical zodiac at the J2000.0 epoch. It effectively measures only the precession from J2000.0.

**Usage:**
Useful for astronomical calculations where a sidereal frame fixed to J2000.0 is needed.

**Formula:**
```
Ayanamsha = (5028.796195 * T + 1.1054348 * T²) / 3600
where T = Julian centuries from J2000.0
```

---

### SIDM_J1900 (19)

**J1900.0**

| Property | Value |
|----------|-------|
| Constant | `SIDM_J1900 = 19` |
| J2000 Value | 1.3966° |
| Zero Epoch | 1900-01-01 12:00 TT |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
Places the tropical and sidereal zodiacs aligned at the J1900.0 epoch. The value at J2000 represents 100 years of precession.

---

### SIDM_B1950 (20)

**B1950.0 (Besselian)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_B1950 = 20` |
| J2000 Value | 0.6984° |
| Zero Epoch | 1950-01-01 (Besselian) |
| Precession Rate | 5027.8 arcsec/century |

**Historical Basis:**
Based on the Besselian epoch B1950.0, which was the standard reference epoch for astronomical catalogs before J2000.0 was adopted.

---

## Additional Lahiri & Krishnamurti Variants

These predefined modes refine the Lahiri and Krishnamurti families with alternative
zero-point determinations.

### SIDM_LAHIRI_1940 (43)

**Lahiri (1940 Lahiri Commission value)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_LAHIRI_1940 = 43` |
| Basis | Lahiri ayanamsha using the 1940 Calendar Reform Committee value |

### SIDM_LAHIRI_VP285 (44)

**Lahiri (Vernal Point 285 CE)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_LAHIRI_VP285 = 44` |
| Basis | Lahiri family anchored to a vernal-point coincidence at 285 CE |

### SIDM_KRISHNAMURTI_VP291 (45)

**Krishnamurti (Vernal Point 291 CE)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_KRISHNAMURTI_VP291 = 45` |
| Basis | Krishnamurti (KP) family anchored to a vernal-point coincidence at 291 CE |

### SIDM_LAHIRI_ICRC (46)

**Lahiri (ICRC, Indian Calendar Reform Committee)**

| Property | Value |
|----------|-------|
| Constant | `SIDM_LAHIRI_ICRC = 46` |
| Basis | Lahiri value as adopted by the Indian Calendar Reform Committee |

---

## User-Defined Ayanamsha

### SIDM_USER (255)

**User-Defined Custom Ayanamsha**

| Property | Value |
|----------|-------|
| Constant | `SIDM_USER = 255` |
| J2000 Value | User-specified |
| Zero Epoch | User-specified |
| Precession Rate | 5027.8 arcsec/century |

**Usage:**
Allows users to define their own ayanamsha using:
- `t0`: Reference epoch (Julian Day)
- `ayan_t0`: Ayanamsha value at t0 (degrees)

**Formula:**
```python
ayanamsha = ayan_t0 + (precession * (tjd_tt - t0) / 36525) / 3600
```

**Example (mimicking Lahiri):**
```python
import libephemeris as leph
from libephemeris.constants import SIDM_USER

# Lahiri parameters: ayan_t0 = 23.857092 at J2000
leph.set_sid_mode(SIDM_USER, t0=2451545.0, ayan_t0=23.857092)
```

---

## Technical Implementation

### Calculation Method

LibEphemeris implements ayanamshas using two methods:

#### 1. Formula-Based Ayanamshas

For most ayanamshas, a fixed formula is used:

```
Ayanamsha = Ayan_J2000 + (Precession_Rate × T) / 3600
```

Where:
- `Ayan_J2000`: Ayanamsha value at J2000.0 epoch (degrees)
- `Precession_Rate`: Rate in arcseconds per century (typically 5027.8)
- `T`: Julian centuries from J2000.0 in TT (Terrestrial Time)

#### 2. Star-Based (True) Ayanamshas

For "True" modes, the actual stellar position is calculated:

1. Get the star's J2000.0 position (RA, Dec) from Hipparcos/Gaia
2. Apply proper motion to the observation date
3. Apply precession and nutation
4. Convert to ecliptic coordinates
5. Calculate: `Ayanamsha = Star_Longitude - Target_Sidereal_Position`

### Mean vs True Ayanamsha

- **Mean Ayanamsha**: `get_ayanamsa_ut()` returns the mean ayanamsha (without nutation)
- **True Ayanamsha**: For sidereal planetary positions, LibEphemeris adds nutation in longitude to get the true ayanamsha (IAU 2006/2000A model)

### Time System

All calculations use:
- **TT (Terrestrial Time)** for astronomical calculations
- **UT1 (Universal Time)** as input, converted internally to TT

---

## References

### Primary Sources

1. **Koch, Dieter & Treindl, Alois.** "Ayanamsha: The Sidereal Zodiac" (technical appendix on ayanamsha epoch definitions).

2. **Indian Astronomical Ephemeris** - Government of India, Positional Astronomy Centre

3. **Suryasiddhanta** - Ancient Indian astronomical text (~400 CE)

4. **Aryabhatiya** - Aryabhata (499 CE)

### Scholarly Works

5. **Kugler, F.X.** (1900-1935). "Sternkunde und Sterndienst in Babel" - Reconstruction of Babylonian astronomy

6. **Britton, J.P.** (2010). "Studies in Babylonian Lunar Theory" - Yale University

7. **Fagan, C. & Bradley, D.** (1971). "Sidereal Astrology" - Statistical basis for Fagan-Bradley ayanamsha

8. **Reid, M.J. & Brunthaler, A.** (2004). "The Proper Motion of Sagittarius A*" - Galactic Center position

### Star Catalogs

9. **Hipparcos Catalogue** (1997) - ESA, high-precision astrometry

10. **Gaia DR3** (2022) - ESA, microarcsecond precision for stellar positions

### Historical Texts

11. **Bailey, A.** (1922-1960). Esoteric writings - Djwhal Khul ayanamsha

12. **Yukteshwar, S.** (1894). "The Holy Science" - Yukteshwar's cosmological system

13. **Valens, V.** (~175 CE). "Anthology" - Hellenistic astrological text
