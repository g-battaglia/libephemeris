# LibEphemeris for Beginners

A book for developers, astrologers and curious observers who want to understand
the sky before calculating it. You need basic Python, not prior astronomy.

Start at the preface and read chapters 1–3 in order. If you already know
positional astronomy, start at chapter 5; the cookbook offers complete practical
applications. The short [technical documentation](../../docs/README.md) is a
reference, not a substitute for the explanations here.

## Chapters

0. [Preface and first calculation](00-preface.md)
1. [The sky seen from Earth](01-sky-from-earth.md)
2. [Measuring time](02-measuring-time.md)
3. [Celestial coordinates](03-celestial-coordinates.md)
4. [Ephemerides](04-ephemeris.md)
5. [Planet positions](05-planet-positions.md)
6. [The Moon](06-moon.md)
7. [Astrological houses](07-astrological-houses.md)
8. [Fixed stars](08-fixed-stars.md)
9. [Eclipses and occultations](09-eclipses.md)
10. [Sunrise, sunset and visibility](10-sunrise-sunset-visibility.md)
11. [The sidereal zodiac](11-sidereal-zodiac.md)
12. [Minor bodies](12-minor-bodies.md)
13. [Hypothetical planets and Arabic parts](13-hypothetical-planets-arabic-parts.md)
14. [Precision, performance and data](14-precision-performance.md)
15. [Cookbook: twelve practical recipes](15-cookbook.md)

## Running the examples

Install the [native beta](../../docs/installation.md), then run this setup before
each complete example. Short fragments without an import continue the preceding
example and reuse its variables. This keeps the exercises independent of your
usual application's settings:

```python
import libephemeris as ephem

ephem.set_calc_mode("leb")
ephem.set_precision_tier("base")
ephem.set_network_policy("sealed")
ephem.set_iers_delta_t_enabled(False)
ephem.set_auto_spk_download(False)
```

The first calculation and most planetary examples work with the bundled core
alone. Unless a civil UTC conversion is shown explicitly, **UT means UT1** in
these exercises. Outputs are rounded examples, not astronomical error bounds.
The discussion of traditional astrology explains conventions, not scientific
evidence for predictions about personal events.

### Additional data and online examples

Lunar apsides, minor bodies and the Chiron natal-chart example can require
additional base groups. JPL/LEB comparisons in chapter 14 require a local JPL
kernel. Prepare those inputs explicitly, then keep the network sealed:

```sh
libephemeris download base
```

Use a dedicated data root with `LIBEPHEMERIS_DATA_DIR` if desired. Confirm body
coverage with `libephemeris status -vv`; do not assume that a tier includes every
asteroid. The [data guide](../../docs/data.md) explains files and configuration.

Blocks that download IERS/SPK/LEB data, query SBDB or select Horizons are
**optional online operations**, not offline exercises. They can write caches
or fetch large files. Read their prerequisites before running them. For a
network operation, deliberately change the policy with
`ephem.set_network_policy("allow")`; afterwards reapply the offline setup.
Path-selection examples require `LIBEPHEMERIS_DATA_DIR/leb/base_core.leb2`.
Online response illustrations are not locally verified output transcripts.

For maintainer checks, `tools/release/check_manual.py` inventories all examples,
executes offline/local-data examples without acquisition, and reports optional
online blocks separately. Never infer a successful online check from an offline
skip. See [contributing](../../CONTRIBUTING.md#documentation).

[Italiano](../it/README.md) · [Provenance and license](../README.md#provenance-and-license)

```{toctree}
:hidden:
:maxdepth: 1

00-preface
01-sky-from-earth
02-measuring-time
03-celestial-coordinates
04-ephemeris
05-planet-positions
06-moon
07-astrological-houses
08-fixed-stars
09-eclipses
10-sunrise-sunset-visibility
11-sidereal-zodiac
12-minor-bodies
13-hypothetical-planets-arabic-parts
14-precision-performance
15-cookbook
```
