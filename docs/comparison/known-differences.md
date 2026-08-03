# Known differences from external implementations

LibEphemeris targets the public reference API while deriving every calculation
independently.

The entries below describe differences that follow from LibEphemeris's own
documented models or from public API behavior. Quantitative accuracy claims are
established against JPL, IAU/ERFA, primary literature, and internal
cross-backend checks.

## Why numerical differences occur

### Ephemeris states

Planetary and lunar states come from NASA JPL DE440/DE441, either directly
through Skyfield, through JPL Horizons, or through LEB Chebyshev approximations
generated from the same independent pipeline. An external implementation can
return a different valid state without exposing which data or reduction it
uses. LibEphemeris does not infer that hidden choice.

The three local tiers have different JPL coverage:

- `base`: DE440s, 1849–2150
- `medium`: DE440, 1550–2650
- `extended`: DE441, approximately −13200 to +17191

Apparent coordinates also depend on Earth orientation, precession, nutation,
light deflection, aberration, and time-scale conversion. Their published model
validity can be narrower than the underlying kernel coverage.

For the outer planets the body-center convention itself is visible: DE
kernels supply the system barycenter, and the local `skyfield` mode refines
it to the physical planet center with the `planet_centers` SPK data. The
geocentric-longitude effect is largest where the moon system is heaviest
relative to the planet: ~0.08" for Pluto/Charon, ~0.05" for Saturn/Titan,
~0.05" for Jupiter/Galileans, and a few milliarcseconds for Uranus and
Neptune. Implementations that stop at the system barycenter — including the
sealed LEB tier, whose channels store the barycentric state — differ from
the refined position by exactly that offset.

### Lunar points

Mean Node and Mean Apogee use ERFA/IERS Delaunay arguments and conventional
lunar-plane geometry. True Node and osculating apsides use the
active JPL Earth–Moon state.

Because the mean points are defined by a published mean-element theory, two
implementations can legitimately use polynomials with slightly different
secular coefficients. Measured against one external implementation, the
mean node and mean apogee longitudes agree at J2000 and drift apart
linearly by about 0.002 arcseconds per year (≈0.1" over the modern century,
≤0.35" at the 1850/2150 tier edges), identically on every backend and
entry point — the signature of a different published rate coefficient, not
of a frame or Delta-T effect. LibEphemeris keeps the IERS 2003 Delaunay
expressions and does not tune the rate toward external output. Interpolated apsides use separate
Delaunay-argument series plus a versioned, hash-pinned residual table,
both fitted at the actual JPL DE440 apsis passages — never to external
output. Their documented tolerances are verified by continuity, parity,
and rate checks.

See [Lunar apsides](../methodology/lunar-apsides.md) and
[Planetary nodes and apsides](../methodology/planetary-nodes-apsides.md).

### True Node under `FLG_ICRS`

For the osculating (true) lunar node, at least one external implementation
applies a longitude shift of about −0.26″ under `FLG_ICRS` that is neither
the frame bias applied to physical bodies (~−0.0065″) nor stable across
epochs (it changes sign around 1963). LibEphemeris applies no ICRS shift to
the node (the mean points behave identically in both implementations). The
divergence is bounded at ~0.35″ and affects only this niche combination.

### Delta T

LibEphemeris uses observed IERS values where available and the independently
published Stephenson–Morrison–Hohenkerk model outside that interval. Delta T
affects every UT-to-TT conversion, so two legitimate models can shift event
times and apparent positions. No sampled external Delta-T series is retained or
used to tune the implementation. Measured against one external implementation,
the difference is ~0.1–0.15 s across 2024–2025 (more recent IERS data on the
LibEphemeris side) and of order 20 s near 1582 (different historical models);
both are bounded model differences, not defects.

See [Delta T](../methodology/delta-t.md).

### UT1-UTC beyond the leap-second table

`utc_to_jd` and its inverses return a real UT1-UTC offset while the
leap-second table can describe the civil label, and the bare calendar JD
afterwards. LibEphemeris places that switch at 2035, following CGPM
Resolution 4 (27th CGPM, 2022), which decides the UT1-UTC tolerance will be
increased by or before that year. Both sides of the switch match an external
implementation exactly where observed IERS data exists (agreement to
≤ 1 ms across 1990-2020) and again once both have switched (exact from 2035
on). Between roughly 2024 and 2035 the two differ by up to ~0.6 s, for two
compounding reasons: UT1-UTC is a *predicted* quantity there, and the
external implementation switches earlier — measured at 2033-10, where its
own predicted offset reaches about -0.95 s and is discontinuously reset to
zero. Neither the switch epoch nor the predicted offset is reconstructed
from external output; the published CGPM decision fixes ours.

### `SIDBIT_SSY_PLANE` zero point

The invariable-plane projection takes its plane from Souami & Souchay
(2012). The plane *orientation* agrees with an external implementation to
better than 0.02″ — projected latitudes match on the calc, planet-centric,
node/apse and fixed-star surfaces alike — but the in-plane zero point
differs by a per-mode constant of about 27–32″ for the modes whose
ayanamsha is ~24° at J2000, shrinking toward zero for modes anchored at or
near a zero J2000 ayanamsha. The offset is stable in time (it varies by
under 0.5″ across 1850–2100), so it is a difference in the published pole
realization rather than a drift. The zero point is constructed
geometrically — the J2000 zero-point direction projected onto the plane —
and no per-mode correction is fitted to external output.

### Precession, nutation, and remote epochs

Modern reductions use IAU 2006/2000A through ERFA. Long-range frame evolution
uses the Vondrák–Capitaine–Wallace 2011 model. At remote epochs, different
published frame and Earth-rotation conventions can produce visible differences;
LibEphemeris reports its model basis instead of fitting a correction.

### Numerical derivatives

Some public speeds are obtained by differentiating the independently computed
position. A numerical derivative can differ from another implementation's
analytical derivative, especially near angle wraps, branch changes, and
singular geometries. LibEphemeris uses wrap-aware central differences and
documents cases in which a speed is undefined or intentionally zero.

### Asteroids

Major asteroids use NASA JPL SPK data over the kernel's coverage. Outside that
coverage, a Keplerian fallback may be available and has lower scientific
fidelity because it does not include the full perturbation history. Arbitrary
numbered asteroids require independently licensed state data; a separately
installed external data file is never consumed.

## Per-API catalog

### `calc*` and coordinate flags

- Geocentric, heliocentric, barycentric, topocentric, ecliptic, equatorial, and
  Cartesian results share the same JPL/ERFA reduction pipeline.
- Requested compatibility flags are echoed according to the public API
  contract, but they do not select or bundle proprietary ephemeris data.
- `FLG_MOSEPH` is accepted as an API flag; calculations remain based on the
  configured JPL backend.
- Deep-time results may differ because JPL kernels, Delta T, precession, and
  Earth-rotation models are independent choices.

### Houses and angles

House positions and cusps are implemented from spherical astronomy and
published house-system definitions. Polar-circle cases, exact horizon
singularities, and equatorial limiting cases require explicit conventions. The
library preserves the public return shape and error behavior where that can be
observed without reconstructing hidden numerical constants.

Cusp and angle speeds are numerical derivatives of the independently computed
geometry. Discontinuities at wraps and branch boundaries are handled as
coordinate discontinuities, not fitted toward external output.

For Sunshine-Makransky (system ``i``), solar declination zero is a regular
limit of Makransky's published construction: the ascensional difference is
zero, the diurnal and nocturnal semiarcs are both 90 degrees, and the eight
intermediate cusps vary continuously through the equinox. At geographic
latitude zero the house-circle poles are also zero, so the complete result is
independent of solar declination and remains cyclically ordered. One external
reference API instead exposes an exact-zero/equatorial branch discontinuity in
``houses_armc`` that disappears for arbitrarily small nonzero declinations.
LibEphemeris intentionally returns the source-defined continuous limit. This is
covered by provenance tests rather than by preserving the reference anomaly.

### Fixed stars

Fixed-star astrometry uses the permissively sourced LibEphemeris catalogue,
proper motions, and the same ERFA/Skyfield frame pipeline as planetary output.
Catalogue coverage, aliases, and ambiguous prefix resolution can differ from an
external catalogue. Numeric strings are interpreted as Hipparcos identifiers;
names and aliases follow LibEphemeris's documented resolver.

Sidereal star speeds are derived from the actual selected sidereal model. No
mode-dependent correction is recovered from external observations.

### Refraction and horizontal coordinates

The standard refraction path uses published Sæmundsson/Bennett-style formulae.
The extended atmospheric path uses an ICAO-standard atmosphere and physical
ray tracing. Below-horizon extrapolation and atmosphere cutoffs are inherently
conventional, so results can differ from other APIs without either result
identifying the other's implementation.

### Eclipses and occultations

Candidate events are found from JPL Sun, Moon, Earth, and target geometry.
Classification uses physical apparent radii, shadow cones, and contact
tangencies. Geographic circumstances use the WGS84 observer model. Event times
therefore inherit the configured JPL and Delta-T models.

No contact-time offset, magnitude scale, geographic correction, or event table
is fitted from compatibility output. Differences that cannot be resolved from
independent astronomical sources remain documented rather than calibrated away.

### `nod_aps*`

Osculating planetary nodes and apsides come from JPL state vectors and standard
two-body vector geometry. Mean lunar points use ERFA/IERS Delaunay arguments;
interpolated lunar apsides use separate analytical series and a
hash-pinned residual table, both fitted at the actual JPL DE440 apsis
passages (no value is fitted to any external implementation; see
scripts/generate_lunar_apse_model.py). Coordinate flags are then applied by the
shared reduction pipeline.

Some combinations, such as barycentric elements for a body without an
independent barycentric state, may return a documented fallback or unsupported
result. The implementation does not manufacture missing values from external
output.

`NODBIT_MEAN` uses the long-term mean-element polynomials of Simon et
al. (1994) for Mercury through Neptune; Pluto has no entry in that
mean-element set, so its mean request falls through to the osculating-state
path and the two methods coincide only for Pluto (see
[the planetary nodes/apsides methodology](../methodology/planetary-nodes-apsides.md)).
`NODBIT_OSCU` reduces the body's own JPL osculating state. Measured longitude
residuals against an
external osculating implementation (both reducing their own JPL state;
geocentric ecliptic of date, 1900–2024 sample) stay well inside the node
tolerance but are non-zero because osculating node and apse *longitudes* are
geometrically ill-conditioned:

- **Apse longitudes** lose definition as eccentricity falls: the
  perihelion/aphelion direction is poorly constrained for a near-circular
  orbit. Perihelion-longitude differences scale inversely with `e` — Mars
  (`e≈0.09`) agrees to ≲0.4″, Jupiter and Saturn (`e≈0.05`) to ≈0.5–1.3″, and
  Neptune (`e≈0.009`) to ≈15–20″. Pluto's apse longitudes differ ≈0.6″.
- **Node longitudes** are amplified by `1/sin i` for low-inclination orbits;
  because each node is reported as the *geocentric* projection of a distinct
  orbit point, the ascending and descending nodes are not exactly 180° apart
  and carry slightly different residuals. Jupiter's nodes differ ≈0.7″
  (ascending) and ≈3.0″ (descending), while the well-inclined Pluto (`i≈17°`)
  agrees to ≈0.03″.

These residuals do not grow monotonically toward remote epochs — Jupiter's
descending-node difference is ≈0.001″ near 1900 and ≈3″ near 2000 — which rules
out a precession, nutation, or epoch error. They are the inherent sensitivity
of the osculating decomposition to sub-milliarcsecond differences between two
independent ephemeris realizations, not a correctable offset; in the ecliptic
of date the underlying orbital-plane orientation agrees closely (node and apse
*latitudes* match to ≈1e-5°).

Under `FLG_J2000` the node/apse are still built on the ecliptic of date; only
their ecliptic *longitude* is precessed to J2000, after which the point is
re-projected geocentrically against the Earth referred to the J2000 ecliptic.
Because a point referred to the of-date plane is viewed against an Earth that
carries the of-date-to-J2000 ecliptic tilt, the J2000 *latitude* is a
distance-dependent parallax reduction of the of-date node — vanishing at J2000
and growing to roughly 130″ of latitude shift a century away — rather than a
rigid frame rotation of the of-date latitude (a rotation would be
distance-independent, but the reduction scales with 1/geocentric-distance).
LibEphemeris reproduces this reduction to ≈1″ for the nodes and to within a few
arcseconds for the high-latitude perihelia at century-scale separations from
J2000.

### `pheno*`

Phase angle, illuminated fraction, elongation, apparent diameter, and magnitude
use vector geometry plus published photometric models. Planetary magnitude
formulae follow Mallama and collaborators where applicable; lunar brightness
uses published almanac relations; angular diameters use independently sourced
IAU radii.

Photometric conventions are empirical and can differ between almanacs. Such
differences are reported as model choices, not used to reverse-engineer a
reference curve.

### Orbital elements

Osculating elements are derived from the active JPL state and standard
state-vector-to-elements equations. Element singularities, reference plane,
central body, and epoch convention can change the representation while
describing the same orbit. Arbitrary asteroids require an independent state
source.

### Sidereal calculations

Every predefined base `SIDM_*` ID 0--46 computes without fallback. Formula
modes share one defining table across backends and live catalogue modes can
differ slightly through catalogue and apparent-place choices. `SIDM_USER`
remains available for callers with their own independently sourced epoch and
zero point. See [Sidereal modes](../reference/ayanamsha.md) for the per-mode
source-audit status and its realization notes.

For a subset of the esoteric and historical modes, the primary statement
leaves a realization freedom (which instant "the year X" means, arcminute
precision of the published value, epoch-fixed versus star-tracking
convention, internal inconsistencies of the source) that different
implementations resolve differently. LibEphemeris realizes the cited
statement literally and does not tune the residual freedom toward external
output. Measured against one external implementation, the resulting stable
zero-point offsets are of order 10--30 arcseconds (De Luce, Djwhal Khul,
Huber, Lahiri-1940, Britton), a few arcseconds for the Siddhantic mean-Sun
and Revati modes, and larger for two modes whose published basis is
qualitatively different: Yukteshwar (~5.5 arcminutes, an unrecoverable
zero-point realization of the 1894 statement) and Aldebaran-15-Tau (offset
of order 2 arcminutes plus ~0.05"/yr, the star-tracking versus epoch-fixed
convention difference documented in the audit notes). These are documented
model choices, not defects; the modes with an unambiguous published anchor
(Fagan/Bradley, Lahiri, Raman, Krishnamurti, ICRC, the frame epochs, the
galactic and star-anchored modes) agree to under an arcsecond. Within that
envelope, several carry a measured *constant* residual with an exactly
matching precession rate: ~0.41" (Fagan/Bradley) and ~0.83" (Krishnamurti)
equal to the external implementation's own removable per-mode
precession-reconciliation offset (see the `SIDBIT_NO_PREC_OFFSET` note
below), and ~0.15–0.27" for the equinox-anchored historical rows, where
"the equinox of year X" admits both a computed-instant and a
conventional-date realization (~2 days apart at remote epochs) and the
measured external choices are not internally consistent across modes.
LibEphemeris keeps the literal published anchors and conventions.

The galactic-anchored modes use the canonical published coordinates — the
VLBI Sgr A* position with its apparent (solar-reflex) proper motion (Reid &
Brunthaler 2004) for the galactic-center modes, and the standard J2000
galactic north pole (192.85948°, +27.12825°) for the equator modes.
Measured against one external implementation, its galactic realizations
differ by a stable 0.08–0.36 arcseconds (center modes ~0.08", IAU-1958
equator modes ~0.19", the mid-Mula equator mode ~0.36"), unaffected by the
projection flags; the modes anchored to a fixed published value at a fixed
date (Fiorenza, Mardyks) agree exactly. LibEphemeris keeps the canonical
coordinates.

Of the `SIDBIT_*` projection flags, the frame projections (`ECL_T0`,
`SSY_PLANE`) and the ecliptic-of-date reference (`ECL_DATE`, realized from
the Vondrák 2011 ecliptic-pole geometry) are implemented for defining-pair
modes and suppressed for live star/galactic modes, whose value is already an
of-date longitude. The projections apply on the planetary, house and
fixed-star surfaces alike. For `SIDM_USER`, the `ECL_T0` projection plane is
the mean ecliptic of the `t0` passed to `set_sid_mode` taken literally (with
`SIDBIT_USER_UT` the stored `t0` is read as UT and converted to TT); the
plane realization uses the Vondrák 2011 long-term model, so agreement with
implementations built on the IAU 2006 developments degrades from the ~0.02"
near-J2000 convergence floor to arcminute level only for epochs several
millennia from J2000 (for example a literal `t0 = 0`, i.e. JD 0).

Equatorial output interacts with the projections as follows (measured
behavior, mirrored exactly): under `SIDBIT_ECL_T0` the RA/Dec channel is
the request reduced to the **mean equator and equinox of the mode's t0**
(no ayanamsha subtraction; position and speeds in the t0 mean frame), on
the calc and fixed-star paths alike; under `SIDBIT_SSY_PLANE` alone the
calc equatorial channel is unchanged, while the fixed-star equatorial
channel is the plain `FLG_J2000|FLG_NONUT` reduction. For Lahiri the
projection epoch is the mode's classical defining anchor — mean ayanamsha
23°15′ at the vernal equinox of 1956 (Report of the Calendar Reform
Committee, CSIR 1955) — while the ayanamsha *value* realization remains
the Indian Astronomical Ephemeris tabulation anchored at J2000; the two
roles are deliberately distinct. Two flags are accepted silently but reduce to the base
value: `SIDBIT_PREC_ORIG` (measured external effect below the ~0.1"
realization floor in the standard configuration) and `SIDBIT_NO_PREC_OFFSET`
(a time-constant per-mode reconciliation offset, non-zero for six modes and
bounded at ~0.83", whose defining values are internal metadata of the
external implementation's per-mode precession-model assignment and are not
reconstructible from published models without fitting output).

### Orbital-element period slots

The tropical-period slot of `get_orbital_elements` is the standard
equinox-referred definition 360°/(n + p) — mean motion plus the IAU
general precession in longitude (Explanatory Supplement, 3rd ed.,
glossary) — expressed in tropical years, and the synodic slot follows the
textbook sidereal-period relation 1/P_syn = 1/P_E − 1/P (e.g. Murray &
Dermott 1999). Measured
externally, one implementation multiplies the tropical-period slot by a
constant sidereal/tropical-year factor (~3.9 × 10⁻⁵ relative, identical
for every planet and flag); that unit convention has no published basis
and is not reproduced.

### Crossing functions

Crossing searches scan and bracket the independently computed longitude or
latitude function, then refine a root. The search window accounts for slow
bodies and wraparound. Returned events can shift with the chosen JPL and Delta-T
models; no reference event list is embedded.

### Heliacal events

Heliacal visibility uses published Schaefer-style visibility principles,
physical extinction, atmospheric refraction, object photometry, and local
Sun/target geometry. Weather, observer acuity, terrain, and twilight thresholds
make visibility predictions model-dependent. LibEphemeris uses documented
physical parameters rather than output-fitted seasonal or geographic terms.

A few auxiliary constants in the visibility chain are declared **project
calibrations** in `libephemeris/schaefer.py`: the humidity-based aerosol
fallback quadratic, the water-vapour V-band coefficient, and the twilight
phase-transition coefficients. Each is constrained by published physics
(Koschmieder 1924 for the visibility relation, Schaefer 1990/1993 and
Rozenberg 1966 for the twilight structure, the canonical ~6.5 mag dark-sky
naked-eye limit) and by continuity/monotonicity requirements — none is
recovered from another implementation's output, and the measured
0.3–0.9 mag limiting-magnitude differences from one external implementation
are the visible consequence of that independence.

### Hypothetical bodies

IDs 40–48 and 50–55, 57–58 calculate from independently transcribed
primary-source models: Neely, Harrington, Le Verrier, Adams, Lowell,
Pickering, and the historical Vulcan, Isis, Proserpina, and Waldemath
literature. ID 56 (Selena/White
Moon) uses Velichko and Larin's published seven-year uniform zodiac cycle,
unwrapped over their published January 1800–January 2000 endpoints and checked
against three unused rows through 2007; its compatibility-only radius is
derived from published IAU nominal constants. ID 49 (Nibiru)
remains named and recognised but raises `UnknownBodyError`; its
complete historical definition was not recovered. Source-exact epochs,
frames, and sexagesimal conversions can intentionally differ from an
undocumented legacy convention where reproducing that convention would require
unverified data.

See [Hypothetical bodies](../methodology/hypothetical-bodies.md) and
[the missing-models inventory](../methodology/missing-hypothetical-models.md).

### Constants, names, and return encodings

Public constants, function names, tuple shapes, flag encodings, and error types
are compatibility targets. They may be tested as ordinary API behavior. Numeric
astronomical datasets and sampled result series are not API metadata and are not
copied or persisted.

## Validation policy

Direct compatibility runs compare public return values and report pass/fail
status or a non-reconstructive aggregate bound; per-date comparison output is
not kept as project data.

Independent validation should be preferred:

- JPL Horizons or a local DE kernel for solar-system states
- ERFA/SOFA for IAU frame and time conventions
- primary literature for photometry, atmosphere, eclipses, and historical
  hypothetical elements
- LEB-versus-direct-JPL checks for binary approximation accuracy

Run `uv run python scripts/check_provenance.py` after provenance-sensitive
changes.
