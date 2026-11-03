# Coverage, precision and limitations

## Choose data for the question

Coverage belongs to each file, body and model, not just a precision-tier name.
An event search needs a surrounding interval; a file containing the starting
epoch may still be too short. Inspect the selected assets with
`libephemeris status -vv` and `libephemeris leb info FILE --json`.

- Missing/corrupt explicit LEB data are errors in forced LEB mode. Provision
  the required group or deliberately choose another supported provider.
- For minor bodies, enable `set_strict_precision(True)` when an approximate
  fallback must be refused. An osculating Kepler orbit is not a long-term
  perturbed trajectory.
- Some hypothetical-body names have no source-complete model. Unsupported
  calculations fail rather than inventing an orbit.
- Placidus and some other house definitions can fail at polar latitudes.
  Choose another system explicitly, or use the documented fallback helper
  and inspect which system was served.
- Eclipse and visibility searches have bounded numerical domains. Missing or
  unbracketed contact slots may be zero; do not turn every zero into a claim
  that an event is physically absent.

## Interpreting precision

More decimal places do not mean greater physical accuracy. Keep these separate:

1. **Source uncertainty:** observational and dynamical-model uncertainty.
2. **Interpolation/numerics:** how accurately stored or computed models are evaluated.
3. **Conventions:** clocks, reference frame, origin, apparent corrections,
   sidereal definition and observer location.

Align these before comparing programs. For historical events, uncertainty in
Earth's rotation (Delta-T) can dominate a very precise orbital interpolation.
LEB component error budgets are not blanket apparent-angle or event-time bounds.
See the [manual's precision chapter](../beginners-manual/en/14-precision-performance.md)
and [scientific sources](sources.md).

The v3.2.1 migration criteria were 0.001 arcsecond for angles, relative 1e-8
for distances, relative 1e-6 for rates (absolute 1e-6 degree/day near zero),
1 microsecond for clocks, 1 second for event times, 1e-6 for event geometry,
1e-5 degree for geographic coordinates and 0.001 magnitude for photometry.
These are comparison criteria under aligned, authenticated inputs, **not a
claim that every configuration meets them or an absolute-accuracy certificate**.
Flags, indices and public result layouts have separate compatibility checks.

Canonical GMST separates UT1 Earth rotation from a continuous TT celestial
origin using the VCW pole model and an IAU 2006 datum. Its numerical domain is
±2000 Julian centuries from J2000; the recorded 29,265-pair comparison reached
0.002420222 microarcseconds against a high-precision construction. That is a
finite numerical comparison, not observational accuracy of ancient precession.
See [the implementation](../core/src/time/canonical_sidereal.rs) for the model.

## Models with additional scope

**Stars:** the compiled 1,447-entry catalog combines astrometry and static
photometry. Proper-motion uncertainty grows away from the catalog epoch;
catalog magnitude does not predict a variable star's current brightness.
Spica's Johnson V magnitude is 0.98 (Hipparcos H5), correcting v3's 1.04.

**Atmosphere:** refraction and naked-eye visibility depend on pressure,
temperature, extinction, twilight and the observer. Positive geometric altitude
alone does not establish visibility. A model cannot account for local clouds
or terrain unless those inputs are explicitly represented.

**N-body:** integration is native Rust, not REBOUND/ASSIST. Named solar maps
can expose stationary-Sun/massless-target reductions rather than general
many-particle encounter algorithms. A prescribed JPL trajectory and the chosen
force model need not match: JPL's N16 contains trajectories from the larger
mutually perturbed N373 solution. Integrating only the other fifteen N16 bodies
is not the same physics. Selected-file constants, perturbers and pole policy
must be stated when interpreting a propagation result.

## Runtime and distribution

Global setters are process-wide; use contexts for observer/sidereal state and
configure shared resources before workers start. Native calls release the GIL
but remain synchronous. Cancelling an awaiting request does not stop a solver.
Indexed coefficient limits do not cap total process memory.

Use a wheel built for your OS/architecture. ABI3 targets supported GIL-enabled
CPython versions, not arbitrary free-threaded interpreters. Source installs
need the native toolchain. See [installation](installation.md) and
[migration](../MIGRATION.md) for supported interfaces and changed behavior.
