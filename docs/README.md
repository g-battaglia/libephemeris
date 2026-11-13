# Documentation

LibEphemeris supplies astrological ephemerides and astronomical calculations.
For higher-level chart data and rendering, use [Kerykeion](https://github.com/g-battaglia/kerykeion);
for hosted calculations, see [Astrologer API](https://www.kerykeion.net/astrologer-api).

Read the [v4 beta release notes](../RELEASE_NOTES.md) for the Rust rewrite,
API compatibility and data-v4 migration overview.

## Start here

1. [Install the library](installation.md) and run a first calculation offline.
2. Follow the complete [Italian](../beginners-manual/it/README.md) or
   [English](../beginners-manual/en/README.md) beginner manual.
3. Check [Python conventions](python.md#time-and-coordinates) or the
   [Rust API](rust.md) before integrating clocks, coordinates and result types.

## Build an astrological workflow

| Task | Worked example or explanation |
| --- | --- |
| Calculate a natal chart | [Birth chart recipe](../beginners-manual/en/15-cookbook.md#recipe-1--complete-birth-chart) |
| Follow transits and aspects | [Daily transits](../beginners-manual/en/15-cookbook.md#recipe-2--daily-transits) · [Aspects](../beginners-manual/en/15-cookbook.md#recipe-9--aspects-between-planets) |
| Choose houses or a zodiac | [House systems](../beginners-manual/en/07-astrological-houses.md) · [Sidereal zodiac](../beginners-manual/en/11-sidereal-zodiac.md) |
| Find a return or lunar phase | [Solar return](../beginners-manual/en/15-cookbook.md#recipe-10--solar-return) · [Lunar phases](../beginners-manual/en/15-cookbook.md#recipe-5--new-moon-and-full-moon) |

Recipes state their setup and data prerequisites. Both manuals cover the same
material; each chapter links to its translation.

## Technical reference

| Topic | Guide |
| --- | --- |
| Signatures, flags, units and results | [Python](python.md) · [Generated API](api_reference.rst) · [Rust](rust.md) |
| Datasets, coverage and offline setup | [Data and configuration](data.md) · [CLI](cli.md) |
| Accuracy, model limits and provenance | [Limitations](limitations.md) · [Scientific sources](sources.md) |
| Read or generate LEB files | [Binary format](data-format.md) |
| Upgrade an existing application | [Migration](../MIGRATION.md) |
| Build, test or contribute | [Contributing](../CONTRIBUTING.md) |
| Release on PyPI | [Releasing to PyPI](releasing.md) |
