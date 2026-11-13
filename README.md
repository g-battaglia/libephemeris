# LibEphemeris

<p>
  <a href="https://pypi.org/project/libephemeris/"><img src="https://img.shields.io/pypi/v/libephemeris.svg" alt="PyPI version"></a>
  <a href="https://pypi.org/project/libephemeris/"><img src="https://img.shields.io/pypi/pyversions/libephemeris.svg" alt="Python versions"></a>
  <a href="https://pepy.tech/project/libephemeris"><img src="https://static.pepy.tech/badge/libephemeris/month" alt="Monthly downloads"></a>
  <a href="https://pepy.tech/project/libephemeris"><img src="https://static.pepy.tech/personalized-badge/libephemeris?period=total&units=INTERNATIONAL_SYSTEM&left_color=GREY&right_color=BLUE&left_text=downloads/total" alt="Total downloads"></a>
  <a href="https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/COPYING"><img src="https://img.shields.io/badge/License-AGPL--3.0-blue.svg" alt="License: AGPL-3.0"></a>
</p>

**A high-precision ephemeris engine for astrology: native Rust speed, a Python API, NASA/JPL DE440/DE441 data and IAU standards.**

**Version:** Python `4.0.0` · Rust core `1.0.0` · Python 3.12+ · [AGPL-3.0](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/COPYING)

[Get started](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/README.md#first-calculation)
· [Scientific precision](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/README.md#scientific-precision)
· [Documentation](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/docs/README.md)
· [English manual](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/beginners-manual/en/README.md)
· [Italian manual](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/beginners-manual/it/README.md)
· [Release notes](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/RELEASE_NOTES.md)

LibEphemeris computes the sky an astrology application is built on: positions
and speeds of the Sun, Moon and planets, houses and angles, sidereal zodiacs,
lunar nodes and apsides, eclipses, stations, risings and settings. Every value
comes from NASA/JPL ephemerides and published astronomical standards, evaluated
by a native Rust engine in microseconds and, with the bundled data, fully
offline.

## Hosted API

If you are building a commercial application, a SaaS, or prefer to keep the codebase closed-source, consider the hosted **[Astrologer API](https://www.kerykeion.net/astrologer-api/subscribe)** on RapidAPI.

Your app consumes LibEphemeris calculations as an external service rather than importing the AGPL library directly — no server setup, no copyleft concerns. Subscribing directly supports the ongoing development of this open-source project.

**[Full API Documentation](https://www.kerykeion.net/astrologer-api)**

## Why LibEphemeris

- **NASA/JPL ephemerides.** Positions come from DE440/DE441, the reference
  planetary and lunar ephemerides of NASA's Jet Propulsion Laboratory. The
  bundled coefficients reproduce direct DE440 evaluation to within
  **0.00024″** for the Sun, Moon and planets.
- **IAU standards, end to end.** IAU 2000A/2006 nutation and Earth orientation
  following the IERS Conventions (2010), through SOFA/ERFA-derived components,
  and Vondrák–Capitaine–Wallace (2011) long-term precession over ±200,000 years.
- **Modern Delta-T.** IERS observations for the modern era and the
  Stephenson–Morrison–Hohenkerk (2016) reconstruction of Earth's rotation for
  historical dates, with UTC, UT1 and TT always kept explicit.
- **Native speed.** About **9 µs** per planetary position and **11 µs** per
  Placidus house solution. A complete chart takes about **0.11 ms**, roughly
  **8× faster than v3**.
- **Offline and reproducible.** The base data ship inside the wheel. A sealed
  network policy, explicit data routes and SHA-256-identified datasets keep
  every result reproducible.
- **Deep time.** Coverage from −13200 to +17191 with DE441, for historical and
  long-range work.
- **Ready for servers.** `EphemerisContext` isolates observer and sidereal
  state per request, and native calls release the GIL.
- **Familiar API, two languages.** The `calc_ut`/`houses`-style functions,
  flags and result layouts that astrology code already uses, in Python, plus the
  `libephemeris-core` crate for Rust applications.

## What it calculates

| Area | Features |
| --- | --- |
| **Charts and transits** | Planetary positions, speeds, retrograde motion and angular separations |
| **Houses and angles** | 25 house systems (26 codes), including Placidus, Koch, Whole Sign, Equal, Regiomontanus, Campanus and Gauquelin sectors; Ascendant, Midheaven and house positions |
| **Sidereal astrology** | 47 predefined ayanamshas, including Lahiri, Fagan–Bradley, Babylonian, true-star and galactic-centre modes, plus a user-defined mode |
| **Lunar points** | Mean/true nodes, mean/osculating apogee and interpolated apsides, each with its own definition |
| **Events and returns** | Longitude crossings, stations, solar-return building blocks, eclipses and occultations |
| **Stars and minor bodies** | 1,447-star catalog; asteroids, centaurs and trans-Neptunian bodies where suitable data exist |
| **Local sky** | Topocentric positions, rise/set/transit, refraction, twilight and visibility models |

The [cookbook](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/beginners-manual/en/15-cookbook.md)
combines these building blocks into twelve practical recipes. Interpretation and
aspect orbs are left to the application.

## Installation

```sh
python -m pip install libephemeris
```

Wheels include the Rust engine and the base coefficient data, so the first
calculation needs no download. NumPy is the only Python runtime dependency.
Prebuilt wheels cover macOS ARM64, Linux x86_64/aarch64 and Windows x86_64;
other platforms can
[build from source](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/docs/installation.md#build-from-source).
Wider date ranges are separate datasets (see
[Data and offline operation](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/README.md#data-and-offline-operation)).

Rust applications can use [`libephemeris-core`](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/docs/rust.md)
directly, without Python: the `bundled-data` feature embeds the base data image
and `nbody` enables native numerical integration.

## First calculation

Tropical positions and Placidus angles for Rome on **8 April 2024 at 12:00 UTC**,
plus a sidereal Moon. The example runs offline with the bundled data:

```python
import libephemeris as ephe

# Fix the data route and clock policy before calculating.
ephe.set_calc_mode("leb")
ephe.set_precision_tier("base")
ephe.set_network_policy("sealed")
ephe.set_iers_delta_t_enabled(False)
ephe.set_auto_spk_download(False)

# utc_to_jd returns TT first, UT1 second. calc_ut and houses need UT1.
jd_tt, jd_ut1 = ephe.utc_to_jd(2024, 4, 8, 12, 0, 0.0)
signs = (
    "Aries", "Taurus", "Gemini", "Cancer", "Leo", "Virgo",
    "Libra", "Scorpio", "Sagittarius", "Capricorn", "Aquarius", "Pisces",
)

positions = {}
for body, name in ((ephe.SUN, "Sun"), (ephe.MOON, "Moon"), (ephe.MERCURY, "Mercury")):
    values, returned_flags = ephe.calc_ut(jd_ut1, body, ephe.FLG_SPEED)
    positions[name] = values
    longitude, speed = values[0], values[3]
    sign = signs[int(longitude // 30)]
    motion = "retrograde" if speed < 0 else "direct"
    print(f"{name}: {longitude % 30:.4f}° {sign}, {speed:+.4f}°/day ({motion})")

# houses takes latitude, then east-positive longitude; cusps[0] is house one.
cusps, angles = ephe.houses(jd_ut1, 41.9028, 12.4964, b"P")
print(f"Ascendant: {angles[0]:.4f}°; Midheaven: {angles[1]:.4f}°")

# A separate context keeps the sidereal choice out of the tropical calculations.
sidereal = ephe.EphemerisContext()
sidereal.set_sid_mode(ephe.SIDM_LAHIRI)
moon_sidereal, _ = sidereal.calc_ut(jd_ut1, ephe.MOON, ephe.FLG_SIDEREAL)
print(f"Moon, Lahiri sidereal longitude: {moon_sidereal[0]:.4f}°")
```

A few conventions to keep in mind:

- **Time scales.** `calc_ut` takes UT1 and `calc` takes TT. Convert local birth
  times to UTC before calling `utc_to_jd`.
- **Defaults.** Positions are apparent and geocentric, in degrees and AU.
- **Polar latitudes.** Some house systems are undefined there.
- **Concurrency.** Keep observer and sidereal settings in an
  `EphemerisContext` and configure shared resources before starting workers.
  Native calls release the GIL but remain synchronous.

[Python conventions](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/docs/python.md)
covers units, flags and time-conversion domains in detail.

## Scientific precision

LibEphemeris builds every calculation on identified data and published models:

| Layer | Source |
| --- | --- |
| **Ephemerides** | NASA/JPL DE440/DE441 (Park et al. 2021), read natively or through local LEB coefficient products |
| **Frames** | IAU 2000A/2006 nutation and Earth orientation, IERS Conventions (2010), SOFA/ERFA-derived components |
| **Precession** | Vondrák, Capitaine & Wallace (2011) long-term model, valid over ±200,000 years |
| **Clocks** | Explicit UTC/UT1/TT handling; IERS observations and Stephenson–Morrison–Hohenkerk (2016) historical Delta-T |
| **Stars** | ESA Hipparcos and XHIP astrometry with IAU/WGSN names |
| **Photometry** | Mallama & Hilton (2018) planetary magnitudes and IAU 2015 physical constants |

Measured results:

- **Bundled data vs. JPL.** Over 3,000 random epochs between 1850 and 2150, the
  bundled LEB coefficients matched direct DE440s evaluation within **0.00024″**
  in longitude for the Sun, Moon and planets (outer planets compared as system
  barycentres).
- **Sidereal time.** The long-term model agrees with a high-precision
  construction to **0.0024 microarcseconds** over 29,265 test pairs spanning
  ±2,000 Julian centuries.

The [source register](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/docs/sources.md)
documents where each model comes from and separates published models from
project conventions. Hypothetical bodies whose published models are incomplete
fail instead of returning invented orbits.

**Numerical precision is not physical accuracy.** Source uncertainty, epoch,
frame, centre and observer all matter, and for historical dates Delta-T can
dominate timing uncertainty. The documented migration criteria include
**0.001″ for angles** and **1 second for events** under aligned inputs; they are
not universal accuracy guarantees. See
[precision and limitations](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/docs/limitations.md)
for scope and qualifications. Scientific accuracy concerns the calculated
quantities, not astrological interpretation.

## Performance

Calculation kernels, JPL/SPK evaluation, LEB reading and N-body integration all
run in Rust. Immutable data and decoded caches are shared across contexts
instead of being rebuilt for each request.

| Operation | v4.0.0 | v3.2.2 |
| --- | --- | --- |
| Position with speed, LEB data | ~9 µs | ~82 µs |
| Position with speed, direct JPL route | ~75 µs | — |
| Placidus houses | ~11 µs | ~30 µs |
| Complete chart: Sun–Pluto, true node and houses | **~0.11 ms** | ~0.92 ms |

A single core computes roughly **9,000 complete charts per second** from the
bundled data. Figures are medians of warmed calls with the published wheels on
an Apple M5 Pro with CPython 3.12, in sealed LEB mode on identical data; actual
times depend on workload, data route, cache state and hardware.

## Data and offline operation

Every calculation follows one data route, selected with `set_calc_mode`:

| Mode | Route |
| --- | --- |
| `leb` | Local coefficients and models, the fastest route; never substitutes JPL or Horizons data implicitly |
| `skyfield` | Native Rust JPL/SPK reader (historical name; there is no Python Skyfield backend) |
| `auto` | Configured resources and provider routing |
| `horizons` | Explicit JPL Horizons network route |

The `sealed` network policy blocks any download during calculation, so
additional resources are provisioned up front:

```sh
libephemeris init
libephemeris download auto
libephemeris status -vv
```

Nominal JPL tiers cover **1850–2150** (`base`), **1550–2650** (`medium`) and
**−13200–17191** (`extended`). Check the actual coverage per file and body,
including the margin an event search needs around its interval; a wider tier
does not mean higher accuracy. For reproducible results, record data versions
and calculation settings.
See [data setup](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/docs/data.md)
and the [CLI reference](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/docs/cli.md).

## Documentation

- **Beginner manuals:** [English](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/beginners-manual/en/README.md) · [Italiano](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/beginners-manual/it/README.md), fifteen chapters and twelve recipes each.
- **Reference:** [Python](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/docs/python.md) · [Rust](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/docs/rust.md) · [Generated API](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/docs/api_reference.rst) · [Documentation index](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/docs/README.md).
- **Science:** [Sources and provenance](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/docs/sources.md) · [Precision and limitations](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/docs/limitations.md) · [Data format](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/docs/data-format.md).
- **Migration:** [API compatibility and changed behavior](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/MIGRATION.md). Familiar names and flags do not imply identical models or exception handling.
- **Development:** [Contributing](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/CONTRIBUTING.md) · [Rust core](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/core/README.md) · [Python bindings](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/python/README.md) · [Changelog](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/CHANGELOG.md).

## License

LibEphemeris is distributed under the [AGPL-3.0](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/COPYING);
contributions are accepted under the [CLA](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/CLA.md).
External components keep their [own terms](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/THIRD_PARTY_NOTICES.md),
and the manuals keep their [original licensing notice](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/beginners-manual/README.md#provenance-and-license).

For private licensing, contact Giacomo Battaglia at
[kerykeion.astrology@gmail.com](mailto:kerykeion.astrology@gmail.com).
Requests are evaluated case by case.

## Ecosystem

LibEphemeris is the ephemeris engine of the Kerykeion ecosystem. Each project
below builds on the previous one:

- **[Kerykeion](https://github.com/g-battaglia/kerykeion)**: the Python astrology
  library built on LibEphemeris, with chart data, aspects, reports and SVG charts.
- **[Astrologer API](https://www.kerykeion.net/astrologer-api)**: the hosted API
  powered by Kerykeion, for commercial, SaaS, mobile and closed-source
  applications.
  [Subscribe on RapidAPI](https://rapidapi.com/gbattaglia/api/astrologer/pricing)
  · [Coding-agent skill](https://github.com/g-battaglia/Astrologer-API/tree/v5/skills/astrologer-api)
- **[Astrologer Studio](https://www.astrologerstudio.com/)**:  Professional astrology software, built on LibEphemeris and Kerykeion.

Project home: [github.com/g-battaglia/libephemeris](https://github.com/g-battaglia/libephemeris).
