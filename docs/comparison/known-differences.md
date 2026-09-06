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

### `SIDBIT_ECL_T0` plane for the Aldebaran mode

`SIDM_ALDEBARAN_15TAU` is a live star-anchored mode, so its ECL_T0
*projection plane* and its *value* come from different places. The plane is
the year −100 mean ecliptic of the Babylonian normal-star zodiac the mode
belongs to (Huber, Centaurus 5, 1958; Fagan, Zodiacs Old and New, 1950) —
the same plane already cited for the Kugler family — and with that plane the
equatorial channel matches an external implementation exactly.

Applying it is nevertheless blocked by data provisioning rather than by
sourcing. Because the mode's zero point is a live anchor evaluated *at* the
plane epoch, year −100 must be inside the loaded kernel; it is outside the
base and medium DE440 tiers. A catalog-only fallback was tried and removed:
it made the same public call return two longitudes 0.79° apart depending on
which tier happened to be installed, which is a worse defect than the
divergence it closed. The mode therefore keeps the J2000 fallback plane on
every tier, and the resulting difference from an external implementation
(~2′ in longitude) is bounded and stable. No zero-point correction is
fitted.

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

#### Vondrák precession through ERFA {#vondrak-erfa}

The long-term precession (Vondrák, Capitaine & Wallace 2011) is evaluated by
ERFA (`erfa.ltpecl`, `erfa.ltpequ`, `erfa.ltp`) instead of coefficient tables
kept in `sidereal_longterm.py`. ERFA is the reference realization of that
published model, so the of-date mean obliquity and the long-term sidereal time
now match it exactly; the former tables deviated from it by at most 1.7e-10
arcsec on the obliquity and 4.1e-10 arcsec on the sidereal time over
-13000..+17000, a summation-order effect. Every function that reads the
obliquity moves by that much or less, except three Sunshine-system charts
(`'I'`/`'i'`, latitude 70°, at the two ends of the DE441 span) where the cusp
geometry amplifies it to 1.4e-9 arcsec, and the cusp speeds of `houses_ex2`,
whose one-second centered difference divides a 5e-11 arcsec movement into
2e-6 arcsec/day: neither is visible at the 0.001 arcsec the library claims.

#### Sidereal time built from the definition {#sidereal-time-definition}

Greenwich mean sidereal time is the IAU 2006 expression in the era that
expression was fitted for, and the hour angle of the mean equinox of date --
Vondrák precession, the Simon et al. (1994) mean longitude of the Earth
retarded by the light time for unit distance, and UT1 -- everywhere else, the
two joined at the ends of the interval over which the international expression
is adopted.

Away from that interval the answer therefore carries the last places of an
accumulated mean longitude: an angle of 3e4 degrees a century out and 5e6
degrees at the ends of the DE441 span, whose last representable place is
already 1.3e-8 arcsec and 3.3e-6 arcsec respectively. Two evaluations of the
same published expressions that differ only in the order of their sums can
land on adjacent representable numbers there, and the difference is that whole
last place, not a fraction of it: the difference is below one representable
step of the accumulated mean longitude, which is the arithmetic floor of the
construction rather than a loss of precision within it. Evaluated with 50
significant digits at the epochs where this happens, the exact value falls
between the two, so the two implementations sit on opposite sides of it at a
distance of that same one step, both within 5.2e-8 arcsec of it, and neither is
the more nearly correct one in any sense the models support: the mean
longitude's own coefficients are published to fewer digits than that, and the
ΔT the epoch is placed with is uncertain by seconds at those dates. Nothing at this level is visible at the 0.001 arcsec the library claims,
and nothing at this level survives into the era of real use, where the answer
is the IAU 2006 expression itself.

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

#### Krusinski house position at the equator {#krusinski-equator}

The Krusinski–Pisa–Goelzer system (selector ``U``) divides the great circle
through the ascendant and the zenith into twelve arcs of 30 degrees, starting
at the ascendant and running towards the nadir, and carries the divisions to
the ecliptic along hour circles. A body is placed by the arc of the point of
that circle which shares its right ascension.

The construction has no singularity at the equator. At geographic latitude
zero the zenith lies on the celestial equator, the ascendant is an ordinary
point of the eastern horizon, and the circle through the two is an ordinary
great circle; the position varies smoothly through the configuration.
LibEphemeris evaluates it with unit vectors and a two-argument arctangent, so
the smoothness survives in floating point as well.

A formulation that carries the same geometry through a tangent ratio instead
loses digits as the latitude is driven towards zero — about 5e-13 relative at
0.001°, 4e-9 at 1e-8° and 8e-6 at zero — and shows the same stiffness a
thousandth of a degree from a pole, where the house circle is nearly an hour
circle. Measured against a 50-digit evaluation of the construction,
LibEphemeris stays within 8.5e-16 relative at every latitude, including
exactly zero and ±89.999. Output at those latitudes therefore differs from
implementations carrying the tangent form by up to 1.9e-5 of a house
(about 5 arcseconds of house-circle arc at the equator, 1.2e-7 relative at
±89.999); everywhere else the two agree to better than 1e-12 relative.

#### Krusinski house position at the geographic poles {#krusinski-poles}

At ``|latitude| = 90°`` the horizon is the celestial equator and the zenith is
a celestial pole, so the house circle through the ascendant and the zenith is
itself an hour circle. Every other hour circle meets it only at the two
celestial poles, and the projection therefore sends every body to the
quarter-turn or the three-quarter-turn mark of the scale: the position
collapses onto the 4th or the 10th cusp exactly, according to which side of
the ascendant the body's own hour circle falls on. Nothing in the arithmetic
is undefined — what is lost is the ability to tell bodies apart, which is a
property of the geometry at a pole and not of any particular formulation.

That collapse is also the limit of the position as the latitude approaches the
pole from inside, and LibEphemeris returns it for every polar chart. An
implementation that reaches the pole through a formulation whose intermediate
quantities cancel there can instead return chart-dependent values for a
minority of configurations that are geometrically indistinguishable from the
rest; such values differ from the limit by up to 3.2 houses. No exception is
raised at a pole and no body is refused; the answer is simply the degenerate
one the construction defines.
#### Topocentric house position: the exact root {#topocentric-root}

The Topocentric (Polich–Page) house position is the house angle of the one
position circle of the family that passes through the body. That angle is the
root of a transcendental equation — the trisected tangent of the geographic
latitude is a straight line in the house angle, set equal to a sinusoid — so,
unlike Regiomontanus, the system collapses to no closed form and the answer
must be found numerically. LibEphemeris brackets the root in the quadrant the
body occupies and converges to the precision of binary64: the residual of the
defining equation at the returned angle is at the arithmetic floor, below
1e-14 away from the geographic poles. A search halted on a fixed angular
tolerance instead answers an iterate whose distance from the root depends on
the schedule of the search rather than on the geometry; at the level such
searches usually stop, that is of the order of 1e-8 of a house typically and
about 2e-6 at worst.

Two identities that are exact in the geometry are therefore exact in the
arithmetic as well. At geographic latitude zero every pole height is zero,
every position circle is an hour circle and the family degenerates to the
equal division of the equator, so the Topocentric position coincides with the
Meridian system (selector `X`) to the last bits. And a body placed on one of
the twelve Topocentric cusps is answered with that cusp's own number to about
1e-14 of a house, outside the polar circles, where the cusp longitude names
the rising intersection.

At the geographic poles the construction has no answer at all: the tangent of
the latitude is unbounded, every circle whose pole share is non-zero collapses
onto the celestial equator, and a body of zero declination lies on all of them
at once. The family is read a thousandth of a degree inside the pole instead,
which keeps the latitude axis continuous and finite; the answers there pile up
on the meridian, and for a body near the equator they stay sensitive to the
last bits of its declination, which the collapsed family amplifies by the
tangent of that latitude.

### Fixed stars

Fixed-star astrometry uses the permissively sourced LibEphemeris catalogue,
proper motions, and the same ERFA/Skyfield frame pipeline as planetary output.
Catalogue coverage, aliases, and ambiguous prefix resolution can differ from an
external catalogue.

The reference-named `fixstar*` family resolves names with the reference's
exact search semantics: a leading-digit string is a **1-based sequential
number in the sorted catalogue** (so `"12"` is the twelfth catalogue entry,
not HIP 12, and `"32 Leonis"` is the 32nd entry, not the Flamsteed
designation), a leading comma keys an exact Bayer/Flamsteed nomenclature
match, and plain names match traditional names and curated alias-table
keys exactly first (`"Alpha Leonis"`, `"Betelgeux"` and `"Formalhaut"` are
literal alias keys and still resolve), then as a prefix in catalogue order
(`"Reg"` is Regulus). No fuzzy matching, no Hipparcos-number lookup
(`"HIP 49669"` raises), and no word-designation *parsing* happens on this
family: a designation absent from the alias table (`"Gamma Virginis"`)
raises. The richer lookups remain available in LibEphemeris's own helper
`resolve_star_name()` (`libephemeris.fixed_stars`), the migration path for
callers that relied on them.

Note that the sequential-number rule makes numeric-string lookups
inherently catalogue-dependent: two implementations sharing the rule but
sorting different catalogues will generally return **different stars** for
the same numeric string. Numeric lookups are only portable within one
catalogue.

Sidereal star speeds are derived from the actual selected sidereal model. No
mode-dependent correction is recovered from external observations.

### Refraction and horizontal coordinates

The standard refraction path uses published Sæmundsson/Bennett-style formulae.
The extended atmospheric path uses an ICAO-standard atmosphere and physical
ray tracing. Below-horizon extrapolation and atmosphere cutoffs are inherently
conventional, so results can differ from other APIs without either result
identifying the other's implementation.

#### Sinclair refraction crossover {#sinclair-crossover}

The rise/set path (`rise_trans`, `rise_trans_true_hor`) turns a true altitude
into an apparent one with a two-branch refraction model: Sinclair's rational
fit (34.46 + 4.23 h + 0.004 h²) / (1 + 0.505 h + 0.0845 h²) arcminutes at low
altitude, the first-order law 0.97 cot h above it, and Bennett's (1982)
pressure/temperature factor on both. The altitude at which the branches switch
is the root of 0.97 cot h = fit(h), solved at import by bisection down to two
adjacent doubles (17.90410464513974°), rather than a stored 12-digit literal.
The stored value sat 6.7e-9° below the root, so the refraction had a step of
1.3e-12° (4.5e-9″ in the standard atmosphere, 6.4e-9″ at 1100 mbar / −50 °C)
at the switch; the derived root brings it down to one rounding ulp
(1.4e-17°). Only apparent altitudes inside that 6.7e-9° band answer
differently, by at most the old step — orders of magnitude below one ulp of a
Julian Day at a rise or set time, so no rise or set time moves. `refrac` and
`refrac_extended` do not use this model and are unchanged.

#### Horizon dip {#horizon-dip}

The dip of the sea horizon — how far below the horizontal plane an elevated
observer sees the horizon — is the exact tangent geometry
`arccos(a_E / (a_E + h))` on the IERS Conventions (2010) equatorial radius
`a_E = 6 378 136.6 m`, scaled by `sqrt(1 - k)` for the coefficient of
terrestrial refraction

    k = a_E · c · P / T² · (g/R_d + dT/dh)

Every quantity in it is published: `c = (n-1)·T/P = 7.885e-5 K/mbar` from the
refractivity of air `n - 1 = 2.7726e-4` at 15 °C and 1013.25 mbar (Edlén 1966,
*Metrologia* 2, 71; Ciddor 1996, *Appl. Opt.* 35, 1566) in Gladstone–Dale
form, and `g/R_d = 9.80665 / 287.0528 = 0.03416 K/m`, the autoconvective lapse
rate of the ISO 2533/ICAO standard atmosphere. The product form is the
effective-radius reduction of the geodetic literature (Bomford 1980,
*Geodesy*, 4th ed.; Schaefer & Liller 1990, *PASP* 102, 796), applied to the
horizon of an elevated observer as in Thom (1971), *Megalithic Lunar
Observatories*.

Deriving the coefficient from those published values instead of a stored
composite makes the dip 4.3e-4 to 1.3e-3 larger in magnitude than in
LibEphemeris 3.2 and earlier: 0.58″ at 100 m, 1.8″ at 1000 m, 3.7″ at 4000 m
and 5.5″ at 9000 m under the standard atmosphere, and up to 10.6″ at 9000 m at
the cold, dense, steeply stratified corner of the model's range. Sea level,
zero pressure and the temperature and height edges are unchanged to the bit.
Because rise and set times read the dip through the `horhgt = -100` sentinel
and through the floor of the rise/set refraction inversion, they move too:
0 s at sea level and up to 4.4 s near the polar circle, where a degree of
horizon is worth over an hour. No event is gained or lost and no meridian
transit moves.

The published derivation is the closer one. Taking the refractivity of air
from ERFA's `eraRefco` at the sodium-D reference wavelength as an independent
oracle (dry air, 1013.25 mbar, 15 °C: 2.771965e-4) and putting it through the
same geometry, the former dip was 1.94″ from the oracle at 1000 m and 5.82″ at
9000 m; the derived one is 0.10″ and 0.31″. On rise and set times the former
horizon was up to 4.6 s from the oracle across a spread of sites, elevations
and events; the derived one stays within 0.26 s.

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
calibrations** in `libephemeris/extinction.py`: the humidity/altitude
aerosol fallback, the per-band water-vapour coefficients, and the twilight
sky-brightness anchor values at the phase boundaries. Each is constrained by
published physics (Koschmieder 1924 for the visibility relation, Patat 2003
and Rozenberg 1966 for the twilight structure) and by continuity/monotonicity
requirements — none is recovered from another implementation's output, and
the measured 0.3–0.9 mag limiting-magnitude differences from one external
implementation are the visible consequence of that independence.

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

#### AU conversion factors {#au-conversion-factors}

`AUNIT_TO_LIGHTYEAR` and `AUNIT_TO_PARSEC` are computed in one step from their
definitions — au / (c · Julian year) with the exact IAU 2012 astronomical unit
and the SI speed of light, and π / 648000 from the exact IAU 2015 parsec —
rather than stored as decimal literals. The stored literals were truncated
copies: they differed from the correctly rounded double of the definition by
5.9e-14 (light-year) and 1.8e-14 (parsec) in relative terms. The derived values
are exact to double precision; no computation inside the library reads either
factor, so positions, distances and magnitudes are unaffected.

### Heliacal phenomena

#### Heliacal visibility window {#heliacal-window}

`heliacal_pheno_ut` answers, in slots 12, 13 and 14 with the duration in slot
24, the interval of the requested date in which the sky can actually show the
object: the part of the four hours of twilight adjoining the Sun's rise or set
where the limiting magnitude of `vis_limit_mag` — atmospheric extinction,
twilight brightness, scattered moonlight and the observer's eye already folded
in — is fainter than the object's own magnitude. That condition is not solved
in closed form. The margin between the two magnitudes is evaluated, so the
three instants are the output of a search: LibEphemeris walks the twilight in
half-minute steps, refines the largest margin by golden section, and converges
each end of the interval by bisection to about a tenth of a second. Two
consequences are worth knowing.

The instants are reproducible only to the convergence of that search. An
implementation that stops its own search earlier answers the same interval
displaced by a second or two. That is far inside the uncertainty of the
underlying visibility model: a tenth of a magnitude near the detection limit
is minutes of twilight, so the interval is physically meaningful to minutes
and numerically determined to a fraction of a second.

The sky-brightness model is discontinuous where the Moon crosses the horizon.
The scattered-moonlight term is carried while the Moon is up and dropped once
it is down, so the limiting magnitude steps by a few tenths of a magnitude at
moonrise and at moonset. A search fine enough to sample those few minutes
reports what the model says there — an optimum that falls next to the step
rather than in the darkest part of the twilight, or a window a few minutes wide
that opens when the object rises and closes again when the Moon does. Both are
answers of the visibility model rather than of the window search: the optimum
next to the step really is where the limiting magnitude of `vis_limit_mag` is
furthest above the object, by a few hundredths to a quarter of a magnitude, and
the narrow window really is bounded by two sign changes of the margin. A coarser
search steps over them. Callers who need a window that is robust against that
step can compare the optimum with `vis_limit_mag` at the Moon's own rise and
set.

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
