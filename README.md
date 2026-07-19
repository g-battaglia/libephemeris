# LibEphemeris

<div align="left">
    <img src="https://static.pepy.tech/badge/libephemeris/month" alt="PyPI Downloads">
    <img src="https://static.pepy.tech/badge/libephemeris/week" alt="PyPI Downloads">
    <img src="https://static.pepy.tech/personalized-badge/libephemeris?period=total&units=INTERNATIONAL_SYSTEM&left_color=GREY&right_color=BLUE&left_text=downloads/total" alt="PyPI Downloads">
    <img src="https://img.shields.io/pypi/v/libephemeris.svg" alt="PyPI Version">
    <img src="https://img.shields.io/pypi/pyversions/libephemeris.svg" alt="Python Versions">
    <img src="https://img.shields.io/github/license/g-battaglia/libephemeris.svg" alt="License">
</div>

A high-precision astronomical ephemeris library for Python, powered by NASA JPL DE440/DE441 ephemerides and IAU 2006/2000A standards.

**Drop-in replacement for PySwissEph** - readable Python algorithms, standard debugging, easy deployment on the scientific Python stack (NumPy, Skyfield, pyerfa).

Licensed under the **[GNU Affero General Public License v3.0](LICENSE)**
([AGPL-3.0-only](LICENSING.md)). Built on NASA JPL DE440/DE441
ephemerides, IAU/ERFA standards, and cited primary literature. See
[NOTICE.md](NOTICE.md) for data sources and attribution.

---

## Public-source guarantee

**No astronomical result in LibEphemeris is accepted because it happens to
match another program, and no unexplained number is treated as scientific
authority.** Every model, coefficient, catalogue field, and generated table
must be traceable to public NASA/JPL, IERS, IAU, ESA, or other identified
scientific data; an independently citable standard, catalogue, or publication;
an explicitly identified permissively licensed upstream component; or a fully
stated project convention that is never misrepresented as an observed fact or
historical author's rule. Primary sources are required whenever the project
claims to reproduce a particular author's numerical model.

When LibEphemeris must choose a numerical tolerance, interpolation grid,
compression parameter, fallback, or branch convention that a publication does
not prescribe, the code labels it as a **project-authored choice** and supports
it with a derivation, error budget, convergence test, or mathematical/physical
invariant. If a historical model cannot be reconstructed completely from its
sources, the corresponding body fails closed instead of returning a guessed or
output-fitted orbit.

This is an enforced engineering contract, not a general aspiration:

- the [complete algorithm and data provenance map](docs/methodology/algorithm-provenance.md)
  explains every scientific chain and every project-authored boundary;
- the [machine-readable registry](docs/methodology/provenance-registry.toml)
  classifies every runtime module, generator, development script, and shipped
  scientific asset; and
- `uv run python scripts/check_algorithm_provenance.py` fails if a module is
  undocumented, a project-owned top-level function/class lacks a docstring, a
  scientific component has no source, a generated artifact has no generator,
  a vendored component lacks its upstream/license, or a citation is replaced
  by a vague placeholder.

Public API compatibility is kept separate from scientific provenance. Public
behavior may be measured ephemerally for regression testing, but those outputs
are never stored as model data or fitted into coefficients. The detailed policy
is in [Independence and Methodology](docs/methodology/independence.md).

“Public source” means independently identifiable and citable. Some standards,
books, and papers retain their publishers' copyright; LibEphemeris uses the
facts, equations, and definitions they publish and does not copy their prose or
another implementation's expression.

---

## Features

- **NASA JPL DE440/DE441** - modern planetary ephemerides via Skyfield, with full-range DE441 support for deep-history and far-future work
- **IAU + Vondrák 2011 standards** - long-term precession and of-date mean obliquity (Vondrák 2011, valid ±200,000 years), nutation (IAU 2006/2000A) via the official ERFA library
- **Latest-reconstruction Delta T (TT−UT1)** - IERS-observed values for the atomic-clock era and the most recent published reconstruction of Earth's rotation from ancient eclipse records (Stephenson, Morrison & Hohenkerk 2016 with the Morrison et al. 2021 update) for historical dates; the default realization keeps positions and house angles consistent, while explicit ΔT overrides have a documented Skyfield-mode exception ([details](https://github.com/g-battaglia/libephemeris/blob/main/docs/methodology/delta-t.md))
- **Source-based validation** - numerical checks use NASA JPL states, ERFA/IAU standards, cited literature, and mathematical invariants ([methodology](https://github.com/g-battaglia/libephemeris/blob/main/docs/methodology/independence.md))
- **Four backends, one API** - Skyfield, LEB (~14x speedup), Horizons API, and adaptive auto mode through the same `calc_ut()` interface
- **25 house systems (26 codes); all 47 predefined sidereal modes operational,
  plus the user-defined mode, with per-mode source/audit status**
- **Physical planet centers when covered** - outer planets use JPL center segments when available and the explicit system barycenter otherwise
- **Thread-safe contexts when you need them** - SwissEph-compatible globals for drop-in migration, `EphemerisContext` for concurrent workloads
- **15,000+ years of coverage** - `base`, `medium`, and `extended` precision tiers from modern use to -13200 / +17191 CE
- **Readable Python 3.12+** - the ephemeris algorithms are plain, inspectable Python; clean installs across CI, containers, and serverless from prebuilt scientific wheels

---

## Why LibEphemeris

Swiss Ephemeris is the industry standard for planetary calculations. But its Python binding (pyswisseph) wraps a large opaque C library - hard to build from source, hard to inspect or debug, tied to a single computation model.

LibEphemeris provides the **same API** with a modern foundation:

- **NASA JPL ephemerides** instead of semi-analytical theory - DE440/DE441 are the latest planetary ephemerides from the Jet Propulsion Laboratory, the same data used for spacecraft navigation.
- **IAU + Vondrák 2011 standards** - long-term precession and of-date mean obliquity (Vondrák, Capitaine & Wallace 2011, valid ±200,000 years instead of the IAU 2006 polynomial's few centuries), nutation (IAU 2006/2000A), all computed via the official ERFA library (the open-source implementation of IAU SOFA), not custom routines.
- **Up-to-date Earth-rotation timeline (ΔT)** - the TT↔UT1 conversion uses the *latest* published reconstruction of Earth's rotation from historical eclipse records (Stephenson-Morrison-Hohenkerk 2016 with the Morrison et al. 2021 revision), blended with IERS observations. The default realization is shared by positions and house angles; explicit model/user/IERS overrides affect LEB and Horizons positions but not forced Skyfield-mode positions.
- **Physical planet centers** - Jupiter, Saturn, Uranus, Neptune, and Pluto use
  JPL body-center segments where the published satellite kernels cover the
  requested epoch, with an explicit system-barycenter fallback outside those
  ranges.
- **Readable Python algorithms** - plain, inspectable source and standard debugging instead of an opaque C library. Installs from prebuilt wheels (NumPy/Skyfield/pyerfa) across any platform, CI, or serverless environment.

**Switching from pyswisseph?** Your existing code works with minimal changes. [Migration guide](https://github.com/g-battaglia/libephemeris/blob/main/docs/guides/migration-guide.md).

### Accuracy over deep time

Because house cusps derive from the long-term Vondrák 2011 model (valid ±200,000 years) and houses and bodies share one obliquity and one ΔT, charts stay correct and internally self-consistent across the whole ±13,000-year ephemeris range, where a truncated precession polynomial drifts by degrees. Cusp *speeds* are computed as the genuine dλ/dt of the full house solution, matching the real cusp motion to **< 0.005 °/day** — including the iteratively-solved Placidus and Koch systems near the polar circle, where an analytic speed approximation can be off by tens to hundreds of °/day.

Methodology: [Long-term sidereal time, precession & cusp speeds](https://github.com/g-battaglia/libephemeris/blob/main/docs/methodology/sidereal-time-longterm.md). Full head-to-head with Swiss Ephemeris: [Swiss Ephemeris Comparison](https://github.com/g-battaglia/libephemeris/blob/main/docs/comparison/index.md).

### Intentional divergences from the reference API

LibEphemeris aims for 1:1 API compatibility at the level of function names,
flags, and slot semantics. A small number of behaviors are deliberately
different, and all of them follow one criterion: **when a published definition
or a checkable physical invariant exists, LibEphemeris follows it** — an
ecliptic frame in which the Sun actually stays at latitude zero, a speed slot
that is the true derivative of the reported position, an obscuration that
respects its bounded definition, orbital elements taken from the page where
their author printed them. Where matching the reference would instead require
data fitted to another implementation's output, the value is derived from the
primary source. Every divergence is documented with its rationale in
[Intentional Divergences](https://github.com/g-battaglia/libephemeris/blob/main/docs/comparison/intentional-divergences.md);
the headline items are:

| Area | What LibEphemeris does | Why |
|------|----------------------|-----|
| **`SIDEREAL \| J2000` on lunar points** | Honors the flag uniformly for every lunar point (mean, true, osculating, interpolated). | Without it \|TrueNode − MeanNode\| grows without bound away from J2000 — far beyond the ±1.5° physical oscillation. |
| **Of-date mean obliquity** | Angle between the Vondrák ecliptic and equator poles (self-consistent frame). | A direction on the ecliptic (the Sun) stays at ~0° ecliptic latitude at every epoch; the separately fitted series drifts ~6″ at −3000. |
| **`deltat_ex` purity** | `deltat_ex()` never mutates `get_tid_acc()`. | Delta-T queries should not have side effects on library state. |
| **Total-eclipse obscuration** | `attr[2]` reports 1.0 during totality (bounded fraction). | Obscuration is a covered-area fraction \[0, 1\] by its published definition; the disc-area ratio stays derivable as `attr[1] ** 2`. |
| **Source-audited ayanamshas** | Every mode carries a per-mode evidence status: primary publication, public catalogue/frame geometry, secondary historical attribution, or explicit project convention. | No value is fitted to another implementation's output. Modes without a recovered primary statement are not mislabeled as author transcriptions; all epochs, anchors, propagation rules, and known offsets are documented in the [sidereal-mode inventory](docs/reference/ayanamsha.md). |
| **House-cusp speeds** | The derivative of the reported cusp position (centered finite difference). Whole Sign cusps report 0.0. | The speed slot is the derivative of the position slot, not a convention-dependent analytic approximation. |
| **Phase angle** | Light-time-consistent geometry. | Bounded 15–40″ difference on inner planets; elongation identical. |
| **Historical hypothetical bodies (IDs 40–58)** | Thirteen IDs compute from their primary published sources (Neely 1980, Harrington 1988, Le Verrier 1846, Adams 1846, Lowell 1915, Velichko–Larin 2007 for White Moon/Selena). IDs 48, 49, 54, 55, 57 and 58 keep their names and constants but raise `UnknownBodyError`. | A recognised name is not a numerical model: where no source-complete published definition could be recovered, failing closed is preferable to returning positions that cannot be traced to any source. The field-by-field inventory of what is missing per ID is in [Missing Hypothetical Models](https://github.com/g-battaglia/libephemeris/blob/main/docs/methodology/missing-hypothetical-models.md). |

Everywhere else the goal is indistinguishability, and it is met: modern-era
planetary and lunar positions, house cusps, eclipses, rise/set, and the
majority of sidereal modes agree with the reference API at the sub-arcsecond
level. The divergences above are documented choices with a stated rationale
and a measured bound — not accuracy limits. The full comparison methodology
and per-channel bounds are in
[Known Differences](https://github.com/g-battaglia/libephemeris/blob/main/docs/comparison/known-differences.md).

---

## Quick Start

```python
import libephemeris as swe
from libephemeris.constants import SUN, MOON, FLG_SPEED

jd = swe.julday(2000, 1, 1, 12.0)  # J2000.0

sun, _ = swe.calc_ut(jd, SUN, FLG_SPEED)
moon, _ = swe.calc_ut(jd, MOON, FLG_SPEED)

print(f"Sun:  {sun[0]:.4f} deg, speed {sun[3]:.4f} deg/day")
print(f"Moon: {moon[0]:.4f} deg, speed {moon[3]:.4f} deg/day")
```

```python
# House cusps (Placidus, Rome)
cusps, ascmc = swe.houses(swe.julday(2024, 11, 5, 18.0), 41.9028, 12.4964, b"P")
print(f"ASC: {ascmc[0]:.4f}, MC: {ascmc[1]:.4f}")
```

For concurrent or multi-threaded workloads, use `EphemerisContext` instead of the module-level global state.

More examples: [Getting Started](https://github.com/g-battaglia/libephemeris/blob/main/docs/guides/getting-started.md)

---

## Precision

[Precision report](https://github.com/g-battaglia/libephemeris/blob/main/docs/reference/precision.md) · [full comparison](https://github.com/g-battaglia/libephemeris/blob/main/docs/comparison/index.md).

Numerical correctness is established from independent NASA JPL/ERFA sources,
published defining conditions, explicitly labelled conventions, and
reference-free invariants.

---

## Four Backends, One API

Choose your trade-off between speed, locality, and setup. The same `calc_ut()` interface works across all four modes, from zero-install Horizons lookups to precomputed LEB throughput.

| Mode | Backend | Speed | Use case |
|------|---------|-------|----------|
| `"auto"` | LEB -> Horizons -> Skyfield | adaptive | **Default.** Best onboarding; resolves local or remote data transparently |
| `"skyfield"` | JPL DE440/DE441 via Skyfield | ~120 us | High-precision local JPL workflow |
| `"leb"` | Sealed LEB plus declared local models | ~5 us | Offline, source-pure repeated calculations; never opens JPL/BSP |
| `"horizons"` | NASA JPL Horizons REST API | ~300 ms | No local ephemeris files required |

```python
from libephemeris import set_calc_mode
set_calc_mode("leb")  # or via env: LIBEPHEMERIS_MODE=leb
```

---

## Installation

```bash
pip install libephemeris
```

Out of the box, the wheel includes the LEB2 base-tier core for the 14 core
bodies plus the lightweight Hamburg-body companion (1850–2150). Mean lunar
points come from ERFA/IERS arguments;
interpolated apsides are anchored to the actual JPL DE440 apsis passages
(documented in the [lunar methodology](https://github.com/g-battaglia/libephemeris/blob/main/docs/methodology/interpolated-perigee.md)).
The remaining versioned, SHA-256-pinned data-v3 groups are available through the
normal tier download commands, while local LEB1 files remain supported.

Recommended first-time setup:

```bash
libephemeris init                 # Optional but recommended interactive config
libephemeris download auto        # Download exactly what your config needs
libephemeris status               # Verify installed data and active setup
```

Prefer to install a tier directly? Use one of these:

```bash
libephemeris download base         # 1850-2150, lightweight
libephemeris download medium       # 1550-2650, ~200 MB (recommended)
libephemeris download extended     # DE441 core span; minor-body ranges vary
```

**Optional extras:** `pip install libephemeris[stars]` for star-catalog tooling, `[nbody]` for REBOUND/ASSIST n-body integration (GPL-3.0-or-later components — explicit opt-in), `[all]` for every permissive-licensed runtime extra. [Details](https://github.com/g-battaglia/libephemeris/blob/main/docs/guides/getting-started.md#optional-extras).

---

## Documentation

- [Getting Started](https://github.com/g-battaglia/libephemeris/blob/main/docs/guides/getting-started.md) - installation, ephemeris tiers, first calculations
- [Migration from PySwissEph](https://github.com/g-battaglia/libephemeris/blob/main/docs/guides/migration-guide.md) - API mapping, flag compatibility, known divergences
- [Optional Modules](https://github.com/g-battaglia/libephemeris/blob/main/docs/guides/optional-modules.md) - optional backends and extras (star catalog, n-body, SPK kernels)
- [Precision Tuning](https://github.com/g-battaglia/libephemeris/blob/main/docs/guides/precision-tuning.md) - configuring optional dependencies for maximum precision
- [Computation Tracing](https://github.com/g-battaglia/libephemeris/blob/main/docs/guides/tracing.md) - discover which backend computed each body
- [Complete API Reference](https://github.com/g-battaglia/libephemeris/blob/main/docs/api_reference.rst) - every public function, class, and constant with signatures and examples
- [Precision Report](https://github.com/g-battaglia/libephemeris/blob/main/docs/reference/precision.md) - models chosen and measured accuracy for every calculation
- [Compatibility Comparison](https://github.com/g-battaglia/libephemeris/blob/main/docs/comparison/index.md) - API semantics, known differences, intentional divergences, and clean-room validation policy
- [Long-term sidereal time, precession & cusp speeds](https://github.com/g-battaglia/libephemeris/blob/main/docs/methodology/sidereal-time-longterm.md) - why houses and cusp speeds stay correct over ±13,000 years
- [Delta T (ΔT)](https://github.com/g-battaglia/libephemeris/blob/main/docs/methodology/delta-t.md) - the multi-era ΔT model (IERS + Stephenson-Morrison-Hohenkerk 2016 / Morrison 2021), why it is piecewise, and the model selector
- [Flag Reference](https://github.com/g-battaglia/libephemeris/blob/main/docs/reference/flags.md) - all supported flags with examples
- [House Systems](https://github.com/g-battaglia/libephemeris/blob/main/docs/reference/house-systems.md) - all 25 systems (26 codes) with full methodology
- [Ayanamsha Modes](https://github.com/g-battaglia/libephemeris/blob/main/docs/reference/ayanamsha.md) - all predefined modes, source-audit status, and user mode
- [Known Bugs & Limitations](https://github.com/g-battaglia/libephemeris/blob/main/docs/reference/known-bugs.md) - active issues and backend limitations
- [LEB Binary Ephemeris](https://github.com/g-battaglia/libephemeris/blob/main/docs/leb/guide.md) - format, generation, LEB2 compression
- [Horizons Backend](https://github.com/g-battaglia/libephemeris/blob/main/docs/architecture/horizons-backend.md) - HTTP client, pipeline, precision
- [Architecture](https://github.com/g-battaglia/libephemeris/blob/main/docs/development/architecture-overview.md) - internal design and data flow
- [Methodology](https://github.com/g-battaglia/libephemeris/blob/main/docs/methodology/overview.md) - planet centers, lunar apsides, pyerfa integration
- [Algorithm and Data Provenance](https://github.com/g-battaglia/libephemeris/blob/main/docs/methodology/algorithm-provenance.md) - complete module-by-module source map, project-authored choices, generators, and enforcement gate
- [CLI Reference](https://github.com/g-battaglia/libephemeris/blob/main/CLI.md) - full command reference
- [Changelog](https://github.com/g-battaglia/libephemeris/blob/main/CHANGELOG.md) - release history

---

## Contributing

Every new or materially changed numerical implementation must include an
in-code `Provenance` section and a record in
`docs/methodology/provenance-registry.toml`. Cite the defining public source,
separate source values from project choices, state units/frame/time scale and
validity limits, and identify the generator for checked-in data. Agreement with
a compatibility target is never an acceptable source.

```bash
git clone https://github.com/g-battaglia/libephemeris.git
cd libephemeris && uv pip install -e ".[dev]"
poe lint                           # Ruff lint + auto-fix
uv run python scripts/check_algorithm_provenance.py
poe test:leb:fast                  # Recommended fast unit suite
poe test:skyfield:fast             # Skyfield backend unit suite
```

---

## Part of the Kerykeion Ecosystem

LibEphemeris is the computation engine behind:

- **[Astrologer Studio](https://www.astrologerstudio.com)** — professional online astrology software (in production)
- **[Kerykeion](https://github.com/g-battaglia/kerykeion)** — Python astrology library (v6 alpha)
- **[Astrologer API](https://www.kerykeion.net/astrologer-api)** — hosted REST API for astrology data and SVG charts (upcoming)

Learn more at [kerykeion.net](https://kerykeion.net).

---

## License

Licensed under the **[GNU Affero General Public License v3.0](LICENSE)**
(`AGPL-3.0-only`). See [LICENSING.md](LICENSING.md) for details.

> **Note:** the optional `libephemeris[nbody]` extra pulls in `rebound` and
> `assist` (GPL-3.0-or-later), which are not part of the core install and are
> never bundled. See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

### License history

Pre-releases `v3.0.0rc2` through `v3.0.0rc8` incorrectly declared
`Apache-2.0` in their PyPI metadata. The relicensing was premature: unresolved
components prevent an Apache-2.0 grant. Those releases are to be treated as
**AGPL-3.0-only** and have been yanked from PyPI. All versions from
`v3.0.0rc10` onward correctly declare `AGPL-3.0-only`.

LibEphemeris provides an API-compatible interface — see
[NOTICE.md](NOTICE.md). "Swiss Ephemeris" is a product of Astrodienst AG.
