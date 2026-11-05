# LibEphemeris

**High-precision astrological ephemerides in Rust and Python, built on NASA/JPL data and published astronomical standards.**

**Beta:** Python `4.0.0b4` · Rust `1.0.0-beta.4` · Python 3.12+ · [AGPL-3.0](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/COPYING).

[Get started](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/README.md#first-calculation) · [Scientific precision](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/README.md#scientific-precision)
· [Italian manual](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/beginners-manual/it/README.md)
· [English manual](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/beginners-manual/en/README.md)
· [Documentation](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/docs/README.md)
· [v4 beta release notes](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/RELEASE_NOTES.md)

## Kerykeion ecosystem

LibEphemeris is the ephemeris engine in the
[**Kerykeion ecosystem**](https://github.com/g-battaglia/kerykeion).
Use LibEphemeris for low-level calculations, Kerykeion for chart data, aspects
and SVG rendering, or the hosted services below for ready-to-use applications.

## Astrologer API

For commercial closed-source projects, SaaS products and mobile apps,
[**Astrologer API**](https://www.kerykeion.net/content/astrologer-api/)
provides hosted astrology calculations and charts without a local installation.
Subscriptions support Kerykeion development.

[**Subscribe on RapidAPI**](https://rapidapi.com/gbattaglia/api/astrologer/pricing)
· [API documentation](https://www.kerykeion.net/content/astrologer-api/)
· [Coding-agent skill](https://github.com/g-battaglia/Astrologer-API/tree/v5/skills/astrologer-api)

## Built for astrology

| Area | Features |
| --- | --- |
| **Charts and transits** | Planetary positions, speeds, retrograde motion and angular separations |
| **Houses and angles** | Placidus, Whole Sign, Equal, Koch and other systems; Ascendant, Midheaven and house positions |
| **Sidereal astrology** | Lahiri, Fagan–Bradley and other ayanamshas, plus a user-defined mode |
| **Lunar points** | Mean/true nodes, mean/osculating apogee and interpolated apsides, with distinct definitions |
| **Events and returns** | Longitude crossings, stations, solar-return building blocks, eclipses and occultations |
| **Stars and minor bodies** | 1,447-star catalog; supported asteroids, centaurs and trans-Neptunian bodies with suitable data |
| **Local sky** | Topocentric positions, rise/set/transit, refraction, twilight and visibility models |

The [cookbook](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/beginners-manual/en/15-cookbook.md) turns these primitives into
twelve practical recipes. Interpretation and aspect orbs remain application choices.

## Scientific precision

- **Ephemerides:** NASA/JPL DE440/DE441-family kernels and local LEB coefficient products.
- **Frames:** IAU/IERS conventions, SOFA/ERFA-derived components and Vondrák–Capitaine–Wallace long-term precession.
- **Clocks:** explicit UTC/UT1/TT handling, Stephenson–Morrison–Hohenkerk (2016) historical Delta-T reconstruction and supported IERS inputs.
- **Validation:** identified JPL data, published standards, high-precision mathematical constructions and physical invariants—not output matching as a scientific source.

The [source register](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/docs/sources.md) documents provenance and separates published
models from project conventions. Source-incomplete hypothetical models fail
rather than return invented orbits.

**Numerical precision is not physical accuracy.** Source uncertainty, epoch,
frame, centre and observer all matter; historical Delta-T can dominate timing
uncertainty. The documented migration criteria include **0.001″ for angles** and
**1 second for events** under aligned inputs—not universal accuracy guarantees.
See [precision and limitations](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/docs/limitations.md) for scope and qualifications.
Scientific accuracy concerns the calculated quantities, not astrological interpretation.

## Installation

Once published, install the native beta explicitly:

```sh
python -m pip install "libephemeris==4.0.0b4"
```

Before publication, [build a local wheel](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/docs/installation.md#build-from-source).
Wheels include the Rust engine and modern base coefficients; NumPy is the Python
runtime dependency. Additional datasets are separate.

Rust applications use [`libephemeris-core`](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/docs/rust.md) without Python;
`bundled-data` embeds the base image and `nbody` enables native integration.

## First calculation

Tropical positions and Placidus angles for Rome, **8 April 2024 at 12:00 UTC**,
plus a sidereal Moon position. Runs offline with the bundled data:

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

Default positions are apparent and geocentric, in degrees and AU. `calc_ut`
accepts UT1; `calc` accepts TT. Convert local birth times to UTC first.
[Python conventions](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/docs/python.md) covers units, flags and time-conversion domains.
Some house systems are undefined at polar latitudes.

For concurrency, use contexts for observer/sidereal state and configure shared
resources before starting workers. Native calls release the GIL but remain synchronous.

## Data and offline operation

| Mode | Route |
| --- | --- |
| `leb` | Local coefficients/models; no implicit JPL or Horizons substitution |
| `skyfield` | Native Rust JPL/SPK reader; historical name, no Python Skyfield backend |
| `auto` | Configured resources and provider routing |
| `horizons` | Explicit JPL Horizons network route |

`sealed` blocks calculation-time acquisition. Provision additional resources explicitly:

```sh
libephemeris init
libephemeris download auto
libephemeris status -vv
```

Nominal JPL tiers: **1850–2150** (`base`), **1550–2650** (`medium`),
**−13200–17191** (`extended`). Check actual per-file/body coverage, including the
surrounding interval needed for event searches. A wider tier is not an accuracy
guarantee. Record data versions and calculation settings for reproducibility.
[Data setup](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/docs/data.md) · [CLI](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/docs/cli.md).

## Documentation

- **Complete beginner manuals:** [Italiano](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/beginners-manual/it/README.md) · [English](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/beginners-manual/en/README.md)—fifteen chapters and twelve recipes each.
- **Reference:** [Python](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/docs/python.md) · [Rust](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/docs/rust.md) · [Generated API](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/docs/api_reference.rst) · [Documentation index](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/docs/README.md).
- **Migration:** [API compatibility and changed behavior](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/MIGRATION.md). Familiar names and flags do not imply identical models or exception handling.
- **Development:** [Contributing](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/CONTRIBUTING.md) · [Rust core](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/core/README.md) · [Python bindings](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/python/README.md) · [Changelog](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/CHANGELOG.md).

## License

Code: [AGPL-3.0](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/COPYING). Contributions: [CLA](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/CLA.md).
External components retain their [own terms](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/THIRD_PARTY_NOTICES.md);
manuals retain their [original licensing notice](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b4/beginners-manual/README.md#provenance-and-license).

For inquiries about private licensing, contact Giacomo Battaglia at
[kerykeion.astrology@gmail.com](mailto:kerykeion.astrology@gmail.com).
Requests will be evaluated on a case-by-case basis.

Project home: [github.com/g-battaglia/libephemeris](https://github.com/g-battaglia/libephemeris).

## Astrologer Studio

Prefer a browser application? [**Astrologer Studio**](https://www.astrologerstudio.com/)
uses Kerykeion and Astrologer API for astrological calculations and charts—no
Python installation or ephemeris setup required.

[**Open Astrologer Studio**](https://www.astrologerstudio.com/)
