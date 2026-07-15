# Algorithm and Data Provenance

## The project rule

**LibEphemeris does not accept an astronomical algorithm, coefficient, lookup
table, physical constant, or generated model merely because it resembles the
output of another program.** Every scientific component must instead be
traceable to at least one independently citable source:

1. public scientific data produced by NASA/JPL, IERS, ESA, IAU, or another
   identified scientific institution;
2. a publicly identifiable technical or astronomical standard, including
   IAU/SOFA, IERS, NAIF, and published coordinate/frame conventions;
3. a named primary paper, monograph, catalogue, or historical source with a
   sufficiently precise locator to repeat the derivation; or
4. an explicitly identified permissively licensed upstream component whose
   copyright, license, version, and local changes are preserved.

Project-authored numerical choices are allowed only when the source does not
dictate them. They must be labelled as project choices and justified by a
derivation, a mathematical invariant, an error budget, convergence evidence,
or a deliberately fail-closed behavior. A bare unexplained number is not an
acceptable source.

The complete machine-readable inventory is
[`provenance-registry.toml`](provenance-registry.toml). The command

```bash
uv run python scripts/check_algorithm_provenance.py
```

proves that every Python module in `libephemeris/`, every project script, and
the shipped scientific assets have a registry entry, a source classification,
and live documentation links. It also requires a docstring on every
project-owned top-level function and class, including private numerical
helpers. The check fails on missing modules, duplicate ownership, vague source
placeholders, missing documents, unlabelled generated files, or scientific
components with no source.

## What “public source” means

“Public” describes auditability, not ownership. NASA/JPL kernels, IERS tables,
catalogue records, facts of nature, mathematical identities, and standards can
be inspected independently of LibEphemeris. Some cited papers, books, and
standards remain copyrighted and are not necessarily open-licensed. The
project cites their scientific statements and writes its own code and prose;
it does not claim ownership of the publications or copy their expressive text.

Likewise, API compatibility is an output contract, not a scientific source.
During development another implementation's public API may be called for an
ephemeral pass/fail comparison. Those observations are not persisted as model
tables, fitted as coefficients, or used to choose an undocumented formula.
When a published definition and compatibility behavior differ, the published
definition controls and the measured divergence is documented.

## Provenance classes

| Class | Meaning | Required evidence |
|---|---|---|
| `scientific-model` | Evaluates a physical, astronomical, atmospheric, historical, or astrological model. | Named source for every model family; conventions, units, validity range, and project choices in the module or linked method document. |
| `scientific-pipeline` | Acquires or transforms source states without defining the underlying observations. | Data producer, frame/time standards, correction sequence, fallbacks, and convergence/error policy. |
| `scientific-data` | Stores or obtains source facts, identifiers, elements, or catalogue records. | Dataset or publication, field mapping, epoch/frame/units, and update procedure. |
| `generated-scientific-data` | Ships values produced by a deterministic generator. | Generator path, complete upstream chain, generation method, error statistics, and integrity pin where practical. |
| `numerical-method` | Encodes or evaluates scientific values using a project-native representation. | Published numerical methods plus explicit project format, precision, and corruption/error checks. |
| `public-api` | Defines names, flags, slots, or delegation. | Explicit statement that compatibility affects the interface only and supplies no scientific formula. |
| `infrastructure` / `development-tool` | Configuration, state, logging, CLI, packaging, or developer workflow. | Explicit statement that no astronomical model or coefficient is defined there. |
| `verification-tool` | Measures models or enforces policy. | Identification of the independent truth source or invariant and a statement that measured values do not become runtime coefficients. |
| `vendored` | Contains third-party source code. | Exact upstream/version, retained copyright and license, local-change list, and third-party notice. |

## End-to-end source chains

### Planetary and lunar states

The primary physical state is NASA/JPL DE440 or DE441, documented by Park et
al. (2021). The local Skyfield backend reads the JPL SPK kernel; the Horizons
backend asks the NASA/JPL service for states; the LEB backend evaluates a
project-generated Chebyshev representation sampled from the same registered
pipeline. LEB is therefore a storage/evaluation layer, not another planetary
theory.

The apparent-place chain is explicit:

```text
JPL ICRF barycentric state
  -> observer/center selection
  -> iterative light time
  -> gravitational deflection when requested
  -> stellar aberration when requested
  -> IAU/ERFA bias, precession, and nutation
  -> mean or true ecliptic/equator of the requested epoch
  -> optional topocentric or sidereal transformation
  -> requested spherical or Cartesian representation
```

Every bypass flag removes a named stage. Backend fallbacks are reported by the
tracing API; they do not silently substitute an unidentified model.

### Frames, precession, nutation, and Earth rotation

Modern precession-nutation follows IAU 2006/2000A through ERFA/pyerfa and the
IERS Conventions. Long-term precession and the shared of-date mean ecliptic use
Vondrak, Capitaine, and Wallace (2011), including the published corrigendum.
The geometric long-term sidereal-time construction combines those poles with
the Simon et al. (1994) mean-Earth expression and UT1 Earth rotation. Its
modern/long-term branch join is a project choice documented with continuity
tests in [sidereal-time-longterm.md](sidereal-time-longterm.md).

Delta T is not a single arbitrary polynomial. Modern values come from IERS
observations; historical values follow the Stephenson, Morrison, and Hohenkerk
(2016) reconstruction with the Morrison et al. (2021) update. Era selection,
blending, uncertainty, and user overrides are documented in
[delta-t.md](delta-t.md).

### Lunar nodes and apsides

Mean points use the IERS 2003 Delaunay arguments evaluated by ERFA. True and
osculating points are vector constructions from the active JPL Moon state:
angular momentum defines the orbital plane and node, while the eccentricity
vector defines the instantaneous apsidal direction. Interpolated apogee and
perigee are generated from every refined DE440 lunar-distance extremum across
the reviewed interval. Their Delaunay basis, residual grids, edge policy,
statistics, generator, and hashes are documented in
[lunar-apsides.md](lunar-apsides.md),
[interpolated-apogee.md](interpolated-apogee.md), and
[interpolated-perigee.md](interpolated-perigee.md).

### Houses and chart angles

House systems are spherical constructions over the horizon, equator, meridian,
prime vertical, and ecliptic. The common vector primitives are project-authored
from those definitions. Each named system has a definition, projection rule,
orientation convention, degeneracy policy, and precision record in
[house-systems.md](../reference/house-systems.md). Savard-A, Krusinski-Pisa,
Sunshine/Makransky, and APC additionally document their exact branch-selection
invariants in `house_constructions.py`.

Ascendant, Midheaven, Vertex, Equatorial Ascendant, and their antipodes are
merely named intersections returned by that common geometry; `angles.py` does
not contain a hidden alternate angle theory.

### Classical and historical astrology helpers

Lots, zodiac/nakshatra partitions, dignity tables, and convenience groupings
do not all have the same evidentiary status. The function-by-function matrix in
[classical-astrology-helpers.md](classical-astrology-helpers.md) distinguishes
published definitions from elementary modular identities and project
conventions. In particular, Fortune/Spirit map to al-Biruni section 476, while
the current Love/Faith formulas are disclosed project conventions; BPHS chapter
3 verses 49-56 support the exaltation and friendship rules. Unsourced aliases,
linear convenience scores, and group names are not promoted to ancient facts.

### Historical and astrological definitions

Historical names are not treated as self-authenticating algorithms. For the
small Lots API, the reversible Fortune/Spirit construction used here is
located specifically to al-Biruni, *The Book of Instruction in the Elements
of the Art of Astrology*, R. Ramsay Wright (1934), section 476. Broader
Hellenistic context is not presented as a direct transcription unless an
exact formula locator has been established. The code states the direction of
every arc explicitly as `ASC + A - B` and distinguishes the source rule from
modern presentation choices.

In particular, the fixed-direction Love and Faith convenience functions do
not claim to implement every sect-dependent historical variant. They preserve
their documented public API formulas; the module labels their fixed direction,
names, default sect behavior, and high-latitude day/night switch as project/API
conventions. These historical calculations are modular zodiacal arithmetic,
not physical astronomy and not coefficients inferred from another program.

The larger `contrib` namespace follows the same boundary: zodiac signs,
nakshatra divisions, classical aspects, and dignity tables identify their
declared tradition in code; purely structural compatibility names contain no
hidden scientific model, and unavailable database-backed functions fail
explicitly instead of fabricating records.

### Sidereal modes

`scripts/generate_ayanamsha_definitions.py` is the authoritative source
compiler. Each generated mode carries its source/status record, epoch, zero
point, and method. Published definitions and live public catalogue geometry
are distinguished from secondary attributions and disclosed project
conventions. No ayanamsha value is obtained by fitting another implementation.
The complete mode table and unresolved historical ambiguities are in
[ayanamsha.md](../reference/ayanamsha.md).

### Fixed stars

The shipped catalogue is generated from Hipparcos, the van Leeuwen (2007) new
reduction, identified VizieR cross-catalogues, and the IAU Working Group on Star
Names list. The generator records every field source. Runtime space motion,
proper motion, parallax/radial-velocity handling, and the IAU/Vondrak frame
chain are implemented in `fixed_stars.py`; no proprietary ephemeris star file
is an input.

### Minor and hypothetical bodies

Minor-body states come from public JPL Horizons SPKs when available. The
fallback uses JPL SBDB elements, classical two-body conversion, and explicitly
cited secular/short-period methods. `multi_epoch_data.py` is generated from
public SPK states: 36 bodies, 81 ten-year nodes each, with the exact historical
IAU-1976 generation constants and IAU-2012 au conversion disclosed for
reproducibility. The generator rejects a non-solar center or non-J2000 frame;
Bennu's separately sourced curated node is the declared exception. The
residual-Fourier experiment retained in
`generate_short_period_corrections.py` is labelled non-production: its output
is neither shipped nor imported.

Historical hypothetical bodies have a stricter rule. Every supported CSV field
has a page-level transcription and a documented frame/unit conversion in
[hypothetical-bodies.md](hypothetical-bodies.md). A known name is not enough to
invent an orbit. Models whose defining fields could not be recovered are
listed in [missing-hypothetical-models.md](missing-hypothetical-models.md) and
raise `UnknownBodyError`.

### Planetary moons and physical centers

Public JPL satellite kernels are preferred for moon positions and
planet-center corrections. Analytical fallbacks cite Lieske E5, TASS,
Jacobson/Brozovic orbit solutions, or public JPL satellite elements and state
their simplifications and precision limits. The Galilean implementation has a
standalone derivation in [galilean-e5-spec.md](galilean-e5-spec.md). The TASS
files are the intentionally retained MIT-licensed Stellarium adaptation, not
project-owned code; their notice and local boundary are in the repository's
[third-party notice](https://github.com/g-battaglia/libephemeris/blob/main/THIRD_PARTY_NOTICES.md).

### Eclipses, photometry, atmosphere, and heliacal visibility

Eclipse searches combine JPL states with published shadow-cone/Besselian
geometry and public physical radii. The event brackets, root tolerances, and
return-slot mapping are project decisions. Obscuration is kept within its
published definition as a covered-area fraction.

Planetary magnitudes use the cited Mallama/Hilton, Bowell H-G, and related
published phase relations. Nominal radii and constants use IAU Resolution B3
or an identified JPL source.

Atmospheric routines distinguish three kinds of input:

- published empirical relations, such as Schaefer visibility/extinction and
  Kasten-Young airmass;
- published closed forms, such as Bennett/Saemundsson refraction; and
- the project ray tracer over the ISO/ICAO standard atmosphere.

Quadrature order, numerical inversion tolerances, and user weather defaults are
project choices, labelled as such. They must not be described as measured laws
or silently promoted to universal constants.

### LEB and LEB2 numerical representation

LEB/LEB2 is original project architecture. Its scientific inputs retain the
source chains above. The representation uses:

- Chebyshev fitting and Clenshaw's published recurrence;
- a documented mapping of Julian dates to `[-1, 1]`;
- analytical Chebyshev derivatives with the time-domain scale restored;
- IEEE-754 binary64 mantissa truncation bounded per stored component;
- reversible coefficient ordering and byte shuffle; and
- Zstandard compression.

Degree, segment width, body target, chunk duration, and cache policy are
project parameters. They are accepted only after explicit error sweeps and are
stored or documented with the artifact. See [LEB algorithms](../leb/algorithms.md)
and the [bundled-core attestation](../leb/base-core-provenance.md).

## Complete runtime-module map

The detailed, path-exact source of truth is the TOML registry. This table is
the human navigation layer; every row names one registry component.

| Registry component | Runtime responsibility | Evidence boundary |
|---|---|---|
| `runtime-public-facade` | Top-level API, state, and isolated context dispatch | Architecture and API-compatibility policy |
| `runtime-infrastructure` | Configuration, cache, exceptions, logging, profiling, and tracing | Explicit no-model declarations and architecture |
| `production-cli-and-downloads` | CLI, wizard, hash-pinned data installation | CLI guide and asset attestations |
| `development-cli` | Development command routing | CLI reference; explicit no-model declaration |
| `public-identifiers` | Public IDs, flags, return masks, and NAIF mappings | Flag reference and NAIF standard |
| `coordinate-and-angle-utilities` | Coordinate transforms, chart angles, formatting, and refraction dispatch | Explicit equations, IERS/Meeus, house/refraction records |
| `historical-astrology-helpers` | Lots, zodiac/nakshatra/aspect/dignity arithmetic | Historical locators plus explicit project conventions |
| `astrometry-time-and-earth-orientation` | IAU/ERFA reductions, Vondrak precession, calendars, sidereal time, and Delta T | IERS, ERFA, Vondrak, Simon, and Earth-rotation publications |
| `house-systems` | All spherical house constructions | Per-system source and derivation table |
| `generated-sidereal-definitions` | Generated ayanamsha definitions | Per-mode source table and authoritative generator |
| `planetary-reduction-pipelines` | Skyfield, LEB, and Horizons apparent places | JPL states plus IAU/IERS reductions |
| `crossing-events` | Safeguarded angular event search | Brent/Newton methods and explicit error controls |
| `eclipse-geometry` | Global/local eclipse events and contacts | NASA canon, Meeus, JPL states, physical geometry |
| `atmospheric-refraction` | Closed-form and ray-traced refraction | Bennett/Saemundsson, ISO/ICAO, Auer/Standish |
| `atmospheric-visibility` | Extinction, contrast/visibility, and heliacal events | Schaefer, Kasten-Young, and named atmosphere publications |
| `fixed-star-runtime` | Catalogue lookup and three-dimensional space motion | Hipparcos/van Leeuwen/WGSN and ERFA/Vondrak |
| `generated-star-catalog` | Shipped public-catalogue rows | Catalog field map and authoritative generator |
| `lunar-points` | Mean, true, osculating, and interpolated nodes/apsides | IERS/ERFA and JPL vector derivations |
| `generated-lunar-apse-model` | Fitted DE440 passage model and residual grids | Generator, basis, hashes, and fit statistics |
| `historical-hypothetical-bodies` | Supported primary-source models and fail-closed IDs | Page-level source record and missing-field inventory |
| `minor-body-runtime` | Public body registry, SBDB elements, and Keplerian fallback | JPL services and published orbital mechanics |
| `generated-multi-epoch-elements` | JPL type-21 states converted to orbital elements | Horizons/NAIF field, frame, epoch, and generator record |
| `spk-access` | Horizons SPK download, cache, validation, and type dispatch | JPL Horizons and NAIF SPK standard |
| `optional-nbody` | REBOUND/ASSIST adapter | Upstream papers and optional-license boundary |
| `planetary-moon-spk-runtime` | JPL satellite-kernel positions and center corrections | JPL satellite/NAIF records and IERS frames |
| `project-moon-theories` | Project implementations of published satellite theories | Lieske E5 and named JPL system solutions |
| `vendored-tass` | MIT TASS adaptation and periodic-term table | Identified upstream, retained license, and TASS publication |
| `leb-representation` | Project LEB formats, compression, readers, and evaluation | Clenshaw, IEEE-754, Zstandard, measured error budgets |
| `astropy-evaluation` | Optional independent comparison bridge | Astropy upstream; diagnostics never become model data |
| `vendored-spk-type21` | MIT type-21 reader | Exact upstream/version/license/local edits |

## Generator and script map

Scripts are covered because provenance can be lost before a generated value
ever reaches runtime. The registry classifies every script, including rejected
experiments and retired tools.

| Registry component | Boundary |
|---|---|
| `script-namespace` | Namespace marker only; no executable or scientific content. |
| `manual-builders` | Documentation assembly only; no scientific output. |
| `star-catalog-generators` | Public catalogue inputs; `build_star_catalog_v2.py` is authoritative for the shipped table. |
| `provenance-and-packaging-gates` | Inventory, source, SPDX, and wheel checks; no scientific model. |
| `sidereal-generator` | Registered primary/geometric definitions, secondary attributions, and explicit conventions; no output fitting. |
| `lunar-apse-generator` | Physical DE440 apsis passages and IERS basis; deterministic artifact. |
| `jpl-minor-body-and-center-tools` | Named JPL endpoints and documented state/element transforms. |
| `experimental-short-period-fit` | Rejected residual-fit experiment; output not shipped or imported. |
| `leb-generators` | Registered JPL/IAU runtime channels plus measured compression error. |
| `leb-diagnostics` | Measurement only; thresholds cannot silently become generated coefficients. |
| `cross-backend-verification` | Ephemeral project/JPL backend reports; prohibited as model inputs. |
| `focused-verification` | State, precession, and Galilean invariants; no emitted runtime artifact. |
| `performance-tools` | Timing/profiling only; measurements cannot change numerical models. |
| `regression-fixture-helper` | Snapshots the project's reviewed core for regression only, never as scientific provenance. |
| `retired-release-tool` | Explicit failure-only legacy entry point; no upload or generation path. |

Unit and integration tests are not runtime models or artifact generators, so
they are not enumerated individually in the registry. They are evidence for
the documented algorithms. Any test helper that starts generating a shipped
table must first move into the registered generator boundary.

## Source register

The registry carries canonical citations and stable locators. The most central
sources are:

- Park et al. (2021), “The JPL Planetary and Lunar Ephemerides DE440 and
  DE441,” *AJ* 161:105, [DOI 10.3847/1538-3881/abd414](https://doi.org/10.3847/1538-3881/abd414).
- Petit and Luzum (eds.), *IERS Conventions (2010)*,
  [Technical Note 36](https://iers-conventions.obspm.fr/content/tn36.pdf).
- [ERFA](https://github.com/liberfa/erfa), the permissive IAU SOFA-derived
  standards implementation used through pyerfa.
- Vondrak, Capitaine, and Wallace (2011), “New precession expressions, valid
  for long time intervals,” [DOI 10.1051/0004-6361/201117274](https://doi.org/10.1051/0004-6361/201117274).
- Simon et al. (1994), “Numerical expressions for precession formulae and mean
  elements for the Moon and the planets,” *A&A* 282, 663-683.
- NASA/JPL [Horizons API](https://ssd-api.jpl.nasa.gov/doc/horizons.html),
  [SBDB API](https://ssd-api.jpl.nasa.gov/doc/sbdb.html), and NAIF
  [SPK specification](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html).
- Clenshaw (1955), “A note on the summation of Chebyshev series,”
  [DOI 10.1090/S0025-5718-1955-0071856-0](https://doi.org/10.1090/S0025-5718-1955-0071856-0).
- ESA Hipparcos (1997), van Leeuwen (2007), identified VizieR catalogues, and
  the IAU WGSN catalogue.

Domain-specific papers, historical editions, page locators, data licenses,
and permissive upstreams are recorded in the registry and the linked focused
method documents rather than hidden in release notes.

## Rules for future code

A contribution that introduces or changes a scientific result must include all
of the following in the same change:

1. a `Provenance` section in the module docstring;
2. a registry entry or update naming the source IDs;
3. a function-level explanation beside every non-obvious series, constant
   group, branch convention, and validity limit;
4. a focused methodology page when the derivation does not fit comfortably in
   code comments;
5. a test against a primary dataset, a standards implementation, a published
   checkpoint, or a mathematical/physical invariant;
6. a declaration of every project-authored parameter and why it is safe; and
7. deterministic generation plus integrity metadata for any shipped table or
   binary artifact.

“Matches expected output” is regression evidence, not provenance. “Standard
formula,” “reference documentation,” and “various sources” are deliberately
rejected by the automated gate because they cannot be independently audited.
