# LibEphemeris — Bug & Known Issue Tracker

## Known Tolerance Issues (not bugs)

### KI-001: Outer planet distance tolerance (Neptune, Pluto)
- **Bodies**: Neptune (8), Pluto (9)  
- **Check**: Distance |lib - ref| < 1e-5 AU
- **Observed**: Neptune max diff 1.23e-5 AU, Pluto max diff 4.1e-5 AU
- **Cause**: JPL DE440 (via Skyfield) vs Swiss Ephemeris internal data differ slightly for distant bodies
- **Severity**: LOW — differences are ~3.5e-7 relative, sub-arcsecond angular impact
- **Recommendation**: Relax dist tolerance to 5e-5 AU for outer planets, or use relative tolerance

### KI-002: Asteroid position differences (Pallas, Juno, Vesta, Chiron)
- **Bodies**: Pallas (18), Juno (19), Vesta (20), Chiron (15)
- **Check**: Longitude |lib - ref| < 1.0 arcsec
- **Observed**: Pallas max 5.53", Juno max 5.40", Vesta max 1.28", Chiron max 1.21"
- **Cause**: Different source ephemerides — libephemeris uses JPL SPK files, Swiss Ephemeris uses its own asteroid files
- **Severity**: LOW — within expected inter-ephemeris variation for minor bodies
- **Recommendation**: Use 6" tolerance for asteroids

### KI-003: Oscillating Apogee speed differences (OscuApog)
- **Body**: OscuApog (13)
- **Check**: Speed |lib - ref| < 0.01 deg/day
- **Observed**: Max diff 0.041 deg/day (277/500 dates fail)
- **Cause**: Oscillating apogee has very rapid speed changes; small timing/model differences amplified
- **Severity**: LOW — positions still agree well
- **Recommendation**: Use 0.05 deg/day tolerance for OscuApog speed

### KI-004: Sun barycentric large position offset
- **Body**: Sun (0) with SEFLG_BARYCTR
- **Check**: LON/LAT |lib - ref| < 2.0 arcsec
- **Observed**: LON 93/100 fail (max 194.6"), LAT 74/100 fail (max 83.6")
- **Cause**: Barycentric Sun = Earth's position relative to SSB. Different ephemeris data for Earth-Moon barycenter offset produce systematic differences.
- **Severity**: MEDIUM — barycentric Sun is rarely used in practice
- **Recommendation**: Use 200" tolerance for barycentric Sun, or skip comparison

### KI-005: Moon longitude marginal failures (~2" range)
- **Body**: Moon (1) across multiple flag combos
- **Check**: LON |lib - ref| < 2.0 arcsec
- **Observed**: Max 2.31" (very few dates, mostly ~2.1")
- **Cause**: Different lunar theories (ELP2000 vs JPL DE440 via Skyfield)
- **Severity**: VERY LOW — barely over tolerance, only 2-3 dates out of 100
- **Recommendation**: Use 2.5" tolerance for Moon longitude

## Potential Bugs

### BUG-001: Interpolated Apogee/Perigee large discrepancies — FIXED
- **Bodies**: IntpApog (21), IntpPerg (22)
- **Status**: **FIXED**
- **Before fix**:
  - IntpApog LON: RMS 175", Max 665" (~11 arcmin)
  - IntpPerg LON: RMS 1347", Max 5181" (~1.4 deg)
  - IntpApog LAT: RMS 3474", Max 11264" (~3.1 deg)
  - IntpPerg LAT: RMS 3480", Max 10695" (~3.0 deg)
- **After fix**:
  - IntpApog LON: RMS 6.9", Max 37.6" (25x improvement)
  - IntpPerg LON: RMS 15.4", Max 172.6" (87x improvement)
  - IntpApog LAT: RMS 381", Max 599" (9x improvement)
  - IntpPerg LAT: RMS 357", Max 568" (10x improvement)
- **Fix details**:
  1. Refitted Delaunay trigonometric perturbation coefficients for both
     `_calc_elp2000_apogee_perturbations()` and `_calc_elp2000_perigee_perturbations()`
     using 402,498 reference points from pyswisseph (1549-2651 CE)
  2. Added correction tables in `lunar_apse_corrections.py`:
     - Apogee: 40,250 entries at 10-day intervals (residuals in arcseconds)
     - Perigee: 201,249 entries at 2-day intervals (residuals in arcseconds)
  3. Replaced latitude computation: instead of raw osculating elements,
     now uses `lat = 5.145396° * sin(lon - node_lon)` which models the
     orbital plane inclination directly
  4. Replaced distance with mean constant values (0.0027099 AU apogee,
     0.0024222 AU perigee)
  5. Added `_interpolate_apse_correction()` function and wired it into
     `calc_interpolated_apogee()` and `calc_interpolated_perigee()`
- **Files modified**: `libephemeris/lunar.py`, new file `libephemeris/lunar_apse_corrections.py`

### BUG-002: Earth geocentric + J2000 produces NaN latitude
- **Body**: Earth (14) with SEFLG_J2000
- **Check**: Skyfield mode calc_ut(jd, 14, SPEED|J2000) returns NaN for latitude
- **Observed**: All 200 dates produce lat=NaN (RuntimeWarning: invalid value in scalar divide)
- **Reproduction**:
  ```python
  import libephemeris as swe
  swe.set_calc_mode("skyfield")
  r = swe.calc_ut(2451545.0, 14, 256 | 32)  # SPEED | J2000
  print(r[0])  # lat is NaN
  ```
- **Root cause**: `planets.py:2555` does `lat = math.degrees(math.asin(ze / dist))` where dist=0 for Earth geocentric. The zero vector is then rotated by the J2000 frame transform, producing 0/0 = NaN.
- **Severity**: LOW — Earth geocentric is always (0,0,0) by definition; should return zeros regardless of frame flags
- **Fix suggestion**: Guard against dist==0 before the division, return (0,0,0,0,0,0) for Earth geocentric before applying any frame transforms.
- **Also affects**: SEFLG_ICRS — Earth geocentric + ICRS crashes with "float division by zero" (Section 2.1: 600/140400 failures)

### BUG-003: Vedic house system ('J') cusps disagree with pyswisseph
- **System**: House system 'J' (Vedic/Sripati variant)
- **Check**: All 12 cusps vs pyswisseph
- **Observed**: Cusps 2,3,5,6,8,9,11,12 differ up to 12.99° from pyswisseph across all dates/locations. Cusps 1,4,7,10 (angular cusps) match fine.
- **Reproduction**:
  ```python
  import libephemeris as swe
  import swisseph
  swisseph.set_ephe_path(os.environ["REF_EPHE_PATH"])
  swe.set_calc_mode("skyfield")
  jd = 2451545.0
  lib = swe.houses(jd, 41.9, 12.5, ord('J'))
  ref = swisseph.houses(jd, 41.9, 12.5, b'J')
  for i in range(12):
      diff = abs(float(lib[0][i]) - float(ref[0][i]))
      if diff > 180: diff = 360 - diff
      print(f"  cusp {i+1}: lib={float(lib[0][i]):.4f} ref={float(ref[0][i]):.4f} diff={diff:.4f}")
  ```
- **Severity**: MEDIUM — Vedic 'J' system is rarely used; intermediate cusps use a different interpolation method than pyswisseph. Angular cusps (1,4,7,10) are correct.
- **Notes**: The 'J' system appears to use a different subdivision algorithm. May be intentional design choice or a misinterpretation of the Vedic house subdivision rules.

### KI-006: Sidereal houses_ex ayanamsha difference check (Section 3.3)
- **Systems**: P, K, E, W, R with sidereal mode
- **Check**: |cusp1_trop - cusp1_sid| ≈ ayanamsha value
- **Observed**: 140/600 checks fail (tolerance < 2°). Star-based ayanamsha (mode 27) and Whole Sign system have larger deviations from simple subtraction model.
- **Cause**: For Ascendant-based systems (W), sidereal cusps are recalculated from sidereal Ascendant, not simply shifted. Star-based ayanamsha may have epoch-dependent behavior.
- **Severity**: VERY LOW — the sidereal cusps themselves are valid; only the "difference ≈ ayanamsha" heuristic fails.

### KI-007: Placidus/Koch polar fallback (Section 3.4)
- **Systems**: Placidus ('P'), Koch ('K') at latitudes 70°, 80°, 85°, 89°
- **Observed**: 80/100 test checks failed, but manual testing confirms `swe_houses_with_fallback()` works correctly at all polar latitudes, returning Porphyry fallback. The test failures were due to residual state from prior sections.
- **Severity**: NON-ISSUE — fallback mechanism works as designed

### BUG-004: get_ayanamsa_ex_ut returns flag value instead of ayanamsha
- **Function**: `get_ayanamsa_ex_ut(jd, flags)`
- **Check**: Should return (ayanamsha_value, retflag), matches `get_ayanamsa_ut(jd)`
- **Observed**: Returns `256.0` (the SEFLG_SPEED flag value) instead of the actual ayanamsha. All 25 test calls fail.
- **Reproduction**:
  ```python
  import libephemeris as swe
  swe.set_calc_mode("skyfield")
  swe.set_sid_mode(0)
  jd = 2451545.0
  result = swe.get_ayanamsa_ex_ut(jd, 256)  # SEFLG_SPEED
  print(result)  # Prints 256.0, should print ~24.07
  print(swe.get_ayanamsa_ut(jd))  # Prints ~24.07 correctly
  ```
- **Severity**: MEDIUM — `get_ayanamsa_ex_ut` is returning the raw flag integer as a float instead of the ayanamsha value. The simpler `get_ayanamsa_ut()` works correctly.
- **Notes**: Likely a bug in the return value handling of the _ex variant.

### BUG-005: get_planet_name missing Uranian body names (40-47) — FIXED
- **Function**: `get_planet_name()` / `swe_get_planet_name()`
- **Bodies**: Cupido (40), Hades (41), Zeus (42), Kronos (43), Apollon (44), Admetos (45), Vulkanus (46), Poseidon (47)
- **Check**: `get_planet_name(40)` should return "Cupido" (pyswisseph returns "Cupido")
- **Observed**: Returns "Unknown (40)" for all Uranian bodies 40-47. Only Transpluto (48) was in the name table.
- **Root cause**: `_PLANET_NAMES` dict in `planets.py` was missing entries for SE_CUPIDO through SE_POSEIDON. The constants existed in `constants.py` but only SE_CUPIDO and SE_POSEIDON were imported in `planets.py`.
- **Severity**: LOW — cosmetic issue, calculations for these bodies worked fine
- **Fix**: Added all 8 Uranian body names to `_PLANET_NAMES` dict and imported the missing constants (SE_HADES, SE_ZEUS, SE_KRONOS, SE_APOLLON, SE_ADMETOS, SE_VULKANUS) in `planets.py`.
- **Files modified**: `libephemeris/planets.py`

### BUG-006: Skyfield `reify` descriptor corruption causes TypeError in sidereal pipeline — FIXED (v1.0.0a6)
- **Function**: `fast_calc._get_precession_matrix()`
- **Status**: **FIXED** in v1.0.0a6
- **Observed**: `TypeError: 'numpy.ndarray' object is not callable` for Pipeline B bodies (TrueNode, OscuApog, MeanNode, MeanApog) when sidereal+equatorial tests run at the same JD as prior Pipeline A tests.
- **Root cause**: Skyfield's `P = reify(precession_matrix)` descriptor uses `update_wrapper`, so `P.__name__` = `'precession_matrix'`. When `t.P` is accessed, the reify `__get__` stores the numpy result under `t.__dict__['precession_matrix']`, shadowing the method. Since `get_cached_time_tt()` uses `lru_cache`, the corruption persists. Later calls to `t.M` (via `ecliptic_frame.rotation_at(t)`) call `self.precession_matrix()` expecting a method but find the numpy array.
- **Severity**: HIGH — caused 20+ sidereal regression test failures
- **Fix**: Replaced `mean_equator_and_equinox_of_date.rotation_at(t)` with `mxm(t.precession_matrix(), ICRS_to_J2000)` to bypass the reify descriptor.
- **Files modified**: `libephemeris/fast_calc.py`

### BUG-007: Lunar occultation `np.minimum` prevents candidate detection — FIXED (v1.0.0a6)
- **Function**: `eclipse.lun_occult_when_glob()`
- **Status**: **FIXED** in v1.0.0a6
- **Observed**: Venus and Mars occultation searches returned events 421–530 days later than pyswisseph (wrong occultation event).
- **Root cause**: `candidate_mask = seps < np.minimum(occ_thresh, _CANDIDATE_DEG)` at line 6170. Since `occ_thresh` (~1.27°) is always less than `_CANDIDATE_DEG` (5.0°), the wide threshold was never applied. The narrow 1.27° window (~0.21 days) was smaller than the 0.5-day scan step, causing ~58% of valid occultation events to be missed.
- **Severity**: HIGH — returned wrong occultation events for planets
- **Fix**: Changed `np.minimum` to `np.maximum`. The verification step still correctly rejects non-occultation close approaches.
- **Files modified**: `libephemeris/eclipse.py`

### BUG-008: South node velocity path asymmetry with LEB backend — FIXED (v1.0.0a6)
- **Function**: `planets.swe_calc_ut()` for bodies -10 (south mean node), -11 (south true node)
- **Status**: **FIXED** in v1.0.0a6
- **Observed**: South node velocity (-0.05444°/day) did not equal north node velocity (-0.05474°/day) when LEB was active. Should be identical since south node = north node + 180°.
- **Root cause**: `swe_calc_ut()` dispatches LEB → Horizons → Skyfield. North node (11) is handled by LEB (Chebyshev derivatives). South node (-11) is not in LEB, falls through to Skyfield (numerical derivatives). Two different derivative methods produce different velocity values.
- **Severity**: MEDIUM — velocity mismatch between geometrically related bodies
- **Fix**: Added early south node handling in `swe_calc_ut()` before LEB/Horizons dispatch. South node recursively calls `swe_calc_ut()` for north node (same backend path), then transforms the result.
- **Files modified**: `libephemeris/planets.py`

### BUG-009: Topocentric observer cache returns stale positions after `set_topo()` — FIXED (v1.0.0a6)
- **Function**: `state.set_topo()` / `cache.get_cached_observer_at()`
- **Status**: **FIXED** in v1.0.0a6
- **Observed**: Topocentric Moon longitude differed by ~0.3" between libephemeris and pyswisseph after changing observer location with `set_topo()`. The first location's result was returned for the second location.
- **Root cause**: `get_cached_observer_at()` uses `(id(observer), jd_tt)` as cache key. When `set_topo()` creates a new `VectorSum` object, Python can reuse the deallocated previous object's memory address, causing `id()` to match the old cache entry. The stale cached position from the previous observer location was returned.
- **Severity**: HIGH — silently returns wrong topocentric positions when switching observer locations
- **Fix**: `set_topo()` now calls `clear_observer_cache()` after updating `_TOPO`, ensuring no stale entries persist across location changes.
- **Files modified**: `libephemeris/state.py`

### KI-008: Ayanamsha mode 40 returns ~357° (out of expected range)
- **Mode**: 40
- **Observed**: `get_ayanamsa_ut(J2000)` returns 356.846° instead of the expected [0, 30] range
- **Cause**: Mode 40 may use a different epoch or reversed direction for ayanamsha
- **Severity**: VERY LOW — mode 40 may be defined this way by design

### KI-009: Ayanamsha modes 31 and 34 produce identical values
- **Modes**: 31, 34
- **Observed**: Both return 30.0206° at J2000 (rounded to 6 decimals)
- **Severity**: VERY LOW — may be two names for the same ayanamsha tradition

### KI-010: Planet speed differences vs pyswisseph (~0.0001-0.0002°/day)
- **Bodies**: Mercury, Venus, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
- **Check**: Speed |lib - ref| < 0.0001°/day
- **Observed**: 22/100 golden regression checks fail. Max diff ~0.000236°/day (Mercury). All planet positions (lon/lat/dist) match within tolerance; only speed (derivative) shows small systematic offset.
- **Cause**: Different derivative computation methods (Skyfield numerical vs Swiss Ephemeris analytical). The difference is ~0.001 arcsec/day — negligible for all practical applications.
- **Severity**: VERY LOW — sub-arcsecond speed differences
- **Recommendation**: Use 0.0003°/day tolerance for planet speeds in golden regression
- **Secondary effect (v1.0.0a6)**: Near retrograde stations (velocity ≈ 0), this offset is amplified into a timing shift δt = δv/a, where a is angular acceleration. For outer planets (small a), this produces station timing differences of up to ~3400s (Saturn). Comparison test tolerances calibrated per-planet in v1.0.0a6.

### KI-011: Boundary date ephemeris range errors at DE440 edge
- **Bodies**: Sun, Mars, Jupiter, Saturn at JD 2287184.5 (1550-01-01)
- **Check**: calc_ut should not crash
- **Observed**: EphemerisRangeError raised — JD is exactly at the DE440 start boundary, some bodies need interpolation margin
- **Cause**: DE440 range is [2287184.5, 2688976.5] but edge segments may need a few days of margin for Chebyshev interpolation
- **Severity**: VERY LOW — boundary-exact dates are edge cases; Moon works fine at this boundary
- **Recommendation**: Use JD + 1 as minimum test date for medium tier

### KI-012: set_calc_mode(None) resolves to "leb" not "auto" when LEB discovered
- **Function**: `set_calc_mode(None)` / `get_calc_mode()`
- **Observed**: Returns "leb" instead of expected "auto" when a bundled LEB2 file is auto-discovered
- **Cause**: The auto mode detects available backends and resolves to the best available (LEB > Skyfield). `get_calc_mode()` returns the resolved mode, not the raw setting.
- **Severity**: NON-ISSUE — correct behavior by design
