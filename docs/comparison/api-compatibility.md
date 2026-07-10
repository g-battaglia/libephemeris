# API compatibility and validation

libephemeris aims for 1:1 API compatibility with pyswisseph: same function names,
same flag and body constants, same return-value structures. This page documents the
few signature differences, the validation methodology behind the precision figures,
and the bugs found in libephemeris during that validation.

---

## 1. API signature differences

Developers migrating from pyswisseph should be aware of these.

### Return-value differences

| Function | pyswisseph returns | libephemeris returns |
|----------|-------------------|---------------------|
| `get_ayanamsa_ex_ut` | `(flags, ayanamsa)` | `(retflag, ayanamsa)` |
| `deltat_ex` | `(float, str)` | `float` — ΔT in days |
| `orbit_max_min_true_distance` | `(max, min, true)` | `(max, min, true)` |

### Parameter differences

| Function | pyswisseph signature | libephemeris signature |
|----------|---------------------|----------------------|
| `heliacal_ut` | `(jd, geopos, datm, dobs, name, event, flags)` → `(jd1, jd2, jd3)` | same signature |
| `lun_occult_when_loc` | body can be `int` or `str` | same |

> `heliacal_ut` / `heliacal_pheno_ut` share the reference's signature, return
> shape, body/event acceptance rules and output encodings. The visibility model
> is now Schaefer's VISLIMIT: `vis_limit_mag` limiting magnitude matches to
> ≈0.14 mag median / ≈0.23 mag max in the 3–9° object-altitude regime where
> heliacal events are actually decided (the disagreement is larger at low
> altitude or in daylight, where both sides sit far below the visibility
> threshold anyway); `kact`/`minTAV` and window widths track the reference, and
> 9 of 17 reference matrix event dates match exactly (the rest ±1 day at the
> marginal visibility transition) — see
> [Known differences §13](known-differences.md#13-heliacal-events).

### Structural differences

| Function | pyswisseph | libephemeris |
|----------|-----------|-------------|
| `get_orbital_elements` | Returns flat tuple (50 floats) | Returns flat tuple (50 floats) — identical shape |
| `houses_armc` | `ascmc[3]` = Vertex | `ascmc[3]` = Vertex — identical layout |
| `houses_ex2` | Returns cusp speeds (analytical for some systems) | Returns cusp speeds (numerical true derivative, always computed) |

### libephemeris-only extensions

These functions have no pyswisseph equivalent:

| Function | Description |
|----------|-------------|
| `houses_with_fallback` | Houses with automatic polar-latitude fallback |
| `houses_armc_with_fallback` | ARMC houses with automatic polar-latitude fallback |
| `sol_eclipse_max_time` | Precise maximum eclipse timing |
| `sol_eclipse_how_details` | Comprehensive eclipse circumstances (dict) |
| `sol_eclipse_obscuration_at_loc` | Eclipse obscuration at a geographic location |
| `planet_occult_when_glob` | Planet–planet occultation search (global) |
| `planet_occult_when_loc` | Planet–planet occultation search (local) |
| `calc_angles` | Pre-calculated angles for Arabic parts |
| `calc_eclipse_path_width` | Eclipse path width at a point |
| `calc_eclipse_central_line` | Eclipse central-line coordinates |
| `calc_eclipse_northern_limit` | Eclipse northern-limit coordinates |
| `calc_eclipse_southern_limit` | Eclipse southern-limit coordinates |

---

## 2. Validation methodology

The precision figures in [Measured precision](precision.md) come from a comparison
suite of **1,619 automated tests** across five suites:

| Suite | File | Tests | Focus |
|-------|------|-------|-------|
| 1 | `test_deep_validation.py` | 514 | Planetary positions (10 planets × 12 flag combos × 210 dates), houses, fixed stars, crossings, eclipses, rise/set, coordinates, utilities |
| 2 | `test_deep_validation_2.py` | 444 | Sidereal modes, topocentric, TT variants, nodal/apsides, orbital elements, phenomena, combined flags, eclipse details |
| 3 | `test_deep_validation_3.py` | 95 | Eclipse geography, occultations, Gauquelin sectors, heliacal events, ARMC ex2, orbit distances, cross_ut, helio_cross_ut |
| 4 | `test_deep_validation_4.py` | 63 | Polar fallback houses, eclipse max time/details/obscuration, planet occultations, heliacal_pheno_ut, calc_angles, state functions |
| 5 | `test_compare_helio_bary.py` | 503 | Heliocentric, barycentric, equatorial, XYZ modes (10 bodies × 10 dates × 4 modes), combined flags, return flag verification |

### Parameters

- **Mode:** Skyfield (LEB mode is validated separately).
- **Date range:** 1550–2650 CE (full DE440 range), concentrated around 1900–2100.
- **Sample density:** 100–210 dates per test (equinoxes, solstices, eclipses, random dates).
- **Locations:** 11 global locations including equator, tropics, mid-latitudes, Arctic, Antarctic.
- **Flags:** all individual flags and common combinations.
- **House systems:** all 25 supported systems.

### Tolerances

| Quantity | Tolerance | Notes |
|----------|-----------|-------|
| Ecliptic longitude (planets) | < 0.001° (3.6") | Sub-arcsecond for outer planets |
| Ecliptic longitude (Moon) | < 0.04° (135") | Model difference (see [known differences §2.2](known-differences.md#22-moon-precision)) |
| Ecliptic latitude | < 0.001° | Sub-arcsecond for all bodies |
| Equatorial RA/Dec | < 0.001° | Sub-arcsecond |
| Distance (AU) | 0.0001 AU | |
| Angular velocity | 0.003°/day | Lon speed; lat speed < 0.004°/day |
| True Node/Lilith longitude | < 0.001° (< 0.5") | Verified vs JPL Horizons |
| House cusps | < 0.02" | All 25 systems |
| Fixed star positions | < 0.51" | van Leeuwen 2007, verified vs SIMBAD |
| Eclipse/crossing times | < 8 seconds | Lunar eclipses; solar < 6 sec |
| Sidereal time | < 2 seconds | |
| Heliocentric/Barycentric | < 0.001° | Sub-arcsecond (except bary Sun) |
| XYZ Cartesian | < 0.00005 AU | Sub-arcsecond angular |

These tolerances were determined empirically and independently verified against JPL
Horizons, astropy/ERFA, and SIMBAD (see [Measured precision §2](precision.md#2-independent-verification)).

### Results

| Result | Count |
|--------|-------|
| Passed | 1,614 |
| Failed | 0 |
| Skipped | 5 (convergence-dependent tests) |
| Total | 1,619 |

All 87 public swe-compatible functions are covered by at least one test suite.

---

## 3. Bugs found and fixed in libephemeris

The validation effort uncovered 12 bugs in libephemeris, all since fixed and
verified. They are documented here as evidence of the validation's thoroughness.

### Output flag handling

| # | Bug | Fix | Commit |
|---|-----|-----|--------|
| 1 | `FLG_RADIANS` (8192) silently ignored — positions always in degrees | `_apply_output_flags()` in `planets.py` converts degrees→radians | `ef29a08` |
| 2 | `FLG_XYZ` (4096) silently ignored — positions always spherical | Same helper converts spherical→Cartesian (AU), incl. velocity Jacobian | `ef29a08` |

### Unit and threshold errors

| # | Bug | Fix | Commit |
|---|-----|-----|--------|
| 3 | `pheno_ut` returned apparent diameter in arcseconds; pyswisseph returns degrees | Divide by 3600.0 in `_calc_apparent_diameter()` | `ef29a08` |
| 4 | Lunar eclipse penumbral detection missed eclipses 12–18° from node | `ECLIPSE_LIMIT_LUNAR` 12.0° → 18.5° | `ef29a08` |

### House system

| # | Bug | Fix | Commit |
|---|-----|-----|--------|
| 5 | Sunshine house system (`'I'`) crashed/wrong at latitudes ≥ 67°N where MC falls below horizon | MC-under-horizon detection: flip MC 180°, compute, shift intermediate cusps 180° | `ef29a08` |

### Fixed stars

| # | Bug | Fix | Commit |
|---|-----|-----|--------|
| 6 | `fixstar_mag` returned 0.0 for most stars (only 2 magnitude entries) | Build from full `STAR_CATALOG` | `ef29a08` |
| 7 | 11 stars missing from catalog | Added all 11 from ESA Hipparcos via SIMBAD/CDS; fixed Tejat HIP + duplicate alias | `ef29a08` |
| 8 | `resolve_star_name()` substring matching caused false positives ("Al" → "Aldebaran") | Prefix matching | `ef29a08` |

Separately, fixed-star **flag handling** was hardened: `fixed_stars.py` previously
honored only 4 of 19 SEFLG flags and silently returned tropical coordinates for
`FLG_SIDEREAL`. All meaningful flags are now handled (`FLG_SIDEREAL`, `FLG_J2000`,
`FLG_NONUT`, `FLG_XYZ`, `FLG_RADIANS`, `FLG_TRUEPOS`, `FLG_MOSEPH`, `FLG_SPEED3`,
`FLG_TOPOCTR`).

### Lunar eclipse classification

| # | Bug | Fix | Commit |
|---|-----|-----|--------|
| 9 | Eclipse type misclassification — missing Danjon atmospheric shadow enlargement | `_SHADOW_ENLARGEMENT = 1 + 1/85 ≈ 1.0118` on umbra and penumbra radii | `e082ee5` |
| 10 | P1/P4 contact times inconsistent with `lun_eclipse_when` | Same Danjon factor in `_calc_lunar_eclipse_penumbral_separation()` | `e082ee5` |

The Danjon enlargement is a well-established correction (Earth's atmosphere refracts
sunlight around the limb, making the geometric shadow ~1.2% larger); the factor 1/85
is used by the Astronomical Almanac.

### Algorithmic improvements

| # | Bug | Fix | Commit |
|---|-----|-----|--------|
| 11 | `swe_cross_ut` diverged for slow outer planets near stations | Widened `STATION_SPEED_THRESHOLD` 0.001→0.01°/day; larger Brent windows; `dt_guess * 1.5` scaling; antipodal-point filter | `d39aaf8` |
| 12 | `swe_orbit_max_min_true_distance` returned 2 values in wrong order | Return 3-tuple `(max, min, true)` with true distance from `calc_ut` | `d39aaf8` |

---

## 4. References

1. Park, R.S. et al. (2021). "The JPL Planetary and Lunar Ephemerides DE440 and DE441." *Astronomical Journal*, 161(3), 105.
2. Folkner, W.M. et al. (2014). "The Planetary and Lunar Ephemerides DE430 and DE431." *Interplanetary Network Progress Report*, 196, 1–81.
3. Stephenson, F.R., Morrison, L.V. & Hohenkerk, C.Y. (2016). "Measurement of the Earth's rotation: 720 BC to AD 2015." *Proceedings of the Royal Society A*, 472, 20160404.
4. Espenak, F. & Meeus, J. (2006). "Five Millennium Canon of Solar Eclipses: −1999 to +3000." NASA/TP-2006-214141.
5. Danjon, A. (1951). "Les éclipses de Lune par la pénombre en 1951." *L'Astronomie*, 65, 51–53.
6. Chapront, J. et al. (2002). "A new determination of lunar orbital parameters." *Astronomy & Astrophysics*, 387, 700–709.
7. ESA (1997). *The Hipparcos and Tycho Catalogues*. ESA SP-1200.
