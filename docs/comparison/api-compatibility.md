# API compatibility and validation

LibEphemeris targets the public PySwissEphemeris interface: canonical function
and constant names, parameter order, tuple structure, flag encoding, state
transitions, and documented exceptions. Numerical calculations are independently
implemented from JPL, IAU/ERFA, primary literature, and permissively licensed
catalogues.

## Public interface

Most migration is an import change:

```python
import libephemeris as swe

position, retflag = swe.calc_ut(jd_ut, swe.SUN, swe.FLG_SPEED)
```

Compatibility tests cover:

- positional and keyword parameters
- native Python return types and tuple lengths
- flag implication, exclusivity, and return-flag echo behavior
- global state setters/getters and `close()`/`reset_session()` behavior
- body, star, house-system, and event acceptance
- exception class and unsupported-input behavior

`deltat_ex()` returns Delta T as a plain float. Eclipse, occultation,
`get_orbital_elements()`, `houses*()`, and `fixstar*()` entry points use the
canonical public tuple layouts documented in the API reference.

The strings the name lookups return (`house_name()`, `get_ayanamsa_name()`,
`get_planet_name()` and `contrib.house_system_name()`) belong to the
compatible surface as well: consumers store them in saved charts and printed
output, so they are held as one set of tables in `libephemeris/compat_names.py`
and returned unchanged, spelling and punctuation included.

## LibEphemeris-only extensions

The following convenience/performance APIs are outside the compatible surface:

| Function | Purpose |
|---|---|
| `reset_session` | Reset lightweight mutable state while preserving expensive readers and caches |
| `houses_with_fallback` | Apply an explicit polar-latitude house fallback |
| `houses_armc_with_fallback` | ARMC variant with explicit fallback |
| `sol_eclipse_max_time` | Refine maximum eclipse timing |
| `sol_eclipse_how_details` | Return named eclipse circumstances |
| `sol_eclipse_obscuration_at_loc` | Return a bounded covered-area fraction |
| `planet_occult_when_glob` | Search global planet–planet occultations |
| `planet_occult_when_loc` | Search local planet–planet occultations |
| `calc_angles` | Precompute angles used by Arabic-part calculations |
| `calc_eclipse_path_width` | Compute a path width from physical shadow geometry |
| `calc_eclipse_central_line` | Compute central-line coordinates |
| `calc_eclipse_northern_limit` | Compute northern-limit coordinates |
| `calc_eclipse_southern_limit` | Compute southern-limit coordinates |

These functions have LibEphemeris-defined contracts and should not be treated as
reference-compatible entry points.

## Validation methodology

Compatibility is checked behaviorally: the same public call is issued against
LibEphemeris and the reference API and the outputs are compared, yielding a
pass/fail status or an aggregate residual per surface.

The reference API's source, documentation prose, algorithms, data files, and
internal model choices are out of scope and must never enter this repository or
a validation worktree.

### Independent validation layers

Numerical correctness is established primarily through:

1. local DE-kernel results versus NASA JPL Horizons or direct SPK evaluation;
2. frame/time transforms versus ERFA/SOFA conventions;
3. LEB results versus the direct JPL/Skyfield path from which they are generated;
4. geometric identities for vectors, orbital elements, houses, eclipses, and
   coordinate transformations;
5. primary literature for Delta T, atmosphere, photometry, and historical
   hypothetical-body elements.

Compatibility observations then check the API contract without becoming a
source of shipped numerical behavior.

### Coverage strategy

Targeted suites sample:

- modern and remote epochs within the selected JPL tier;
- equatorial, tropical, temperate, and polar observers;
- individual flags and meaningful combined-flag grids;
- all supported house systems and sidereal-mode classifications;
- planets, lunar points, major asteroids, fixed stars, eclipses, crossings,
  rise/set, and time utilities;
- Skyfield, Horizons, and LEB paths where each is applicable.

Exact tolerances belong in reference-free tests and are justified by the
independent source's published uncertainty, numerical conditioning, or the LEB
approximation budget. They are not fitted from external residuals.

## Defects found through validation

Past validation exposed ordinary implementation defects including ignored output
flags, mixed coordinate frames, unstable phase-angle geometry, incomplete star
catalogue resolution, incorrect return ordering, polar house singularities,
eclipse classification inconsistencies, and crossing searches that failed near
stations. Fixes were derived from the public contract plus independent geometry
or published standards.

Historical per-case reference measurements were removed during review; they
are neither required to understand the fix nor permitted as persisted
validation data.

## References

1. Park, R. S. et al. (2021), “The JPL Planetary and Lunar Ephemerides DE440
   and DE441,” *Astronomical Journal* 161, 105.
2. Stephenson, F. R., Morrison, L. V. & Hohenkerk, C. Y. (2016), “Measurement
   of the Earth's rotation: 720 BC to AD 2015,” *Proceedings of the Royal
   Society A* 472, 20160404.
3. Vondrák, J., Capitaine, N. & Wallace, P. T. (2011), “New precession
   expressions, valid for long time intervals,” *A&A* 534, A22.
4. ESA (1997), *The Hipparcos and Tycho Catalogues*, ESA SP-1200.
5. IAU SOFA and ERFA documentation for standard frame and time routines.
