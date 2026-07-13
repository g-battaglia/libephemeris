# LibEphemeris — Piano di Verifica Esaustiva

## Contesto Generale del Progetto

LibEphemeris e' un'implementazione clean-room in puro Python di un'API
compatibile, basata su NASA JPL DE440/DE441 e standard IAU/IERS. Non usa file di
distribuzione Swiss Ephemeris. Le chiamate pubbliche PySwissEph possono essere
confrontate solo in modo effimero: l'output non puo' essere salvato come
fixture, usato per fit/regressioni, trasformato in coefficienti o dati, oppure
committato.

> Questo documento e' uno snapshot storico del 2026-03-27. I punti relativi ai
> vecchi modelli lunari sono stati superati dalla revisione clean-room 2026-07.

---

## Open Issues from poe test:full (130 failures / 17376 passed)

Last updated: 2026-03-27

### Summary by Root Cause

| # | Root Cause | Failures | Priority | Status |
|---|-----------|----------|----------|--------|
| 1 | `get_ayanamsa_ut` returns numpy.float64, not native float | 2 | HIGH | PENDING |
| 2 | Eclipse local API signatures (`lun_eclipse_when_loc`, `sol_eclipse_when_loc`) | 8 | HIGH | PENDING |
| 3 | LEB files stale — derived lunar channels need independent regeneration | ~60 | MED | SUPERSEDED |
| 4 | Extended tier asteroids — "Invalid Time to evaluate" (SPK date range) | ~40 | MED | PENDING |
| 5 | TrueNode crosstier at boundary JD 2290867.5 (1560) — EphemerisRangeError | 2 | LOW | PENDING |
| 6 | Legacy fitted lunar perturbation consistency | 2 | LOW | RETIRED |
| 7 | Interpolated apogee edge case — range returns negative JD | 1 | MED | PENDING |
| 8 | Precision tier — `de441.bsp` != `de440s.bsp` (tier doesn't switch file) | 1 | LOW | PENDING |
| 9 | Sunshine house 'i' at lat > 58° — missing Makransky algorithm | 16 (verify) | LOW | KNOWN LIMITATION |

### Issue 1: numpy.float64 return type (2 failures)
- **Files**: `test_ayanamsha_comprehensive.py:43`, `test_state_management.py:190`
- **Error**: `assert <class 'numpy.float64'> is float`
- **Fix**: Wrap `get_ayanamsa_ut` return in `float()` in `ayanamsha.py`

### Issue 2: Eclipse local API signatures (8 failures)
- **Files**: `test_compare_leb_eclipses_lunar.py`, `test_compare_leb_eclipses_solar.py`
- **Errors**:
  - `swe_lun_eclipse_when_loc()` takes 2-4 positional args but 5 given
  - `sol_eclipse_when_loc` returns int where sequence expected (`object of type 'int' has no len()`)
- **Fix**: Check and fix API signatures in eclipse module

### Issue 3: LEB derived lunar channels (superseded)

Current generation evaluates mean points from ERFA/IERS arguments and smoothed
apsides from active JPL states. Any bundled LEB must be regenerated only from
those independent models; no persisted compatibility output may be used.

### Issue 4: Extended asteroids date range (~40 failures)
- **Files**: `test_extended_asteroids.py`, `test_compare_leb_asteroids.py`, `test_compare_leb_distances.py`, `test_compare_leb_velocities.py`
- **Error**: `ValueError: Invalid Time to evaluate` for Chiron(15), Ceres(17), Pallas(18), Juno(19), Vesta(20)
- **Cause**: Extended tier (de441.bsp, -13200 to +17191) tests at dates outside asteroid SPK kernel range
- **Fix**: Add date range validation for asteroid bodies, or catch ValueError and raise EphemerisRangeError

### Issue 5: TrueNode crosstier boundary (2 failures)
- **File**: `test_crosstier_medium_extended.py`
- **Error**: `EphemerisRangeError` at JD 2290867.5 (1560-01-31)
- **Cause**: TrueNode needs Moon position slightly outside medium tier range for internal computation
- **Fix**: Widen safety margin for True Node at tier boundaries

### Issue 6: Legacy lunar perturbation consistency (retired)

The output-fitted perturbation series and its expected-value tests were
removed. Current tests use IERS identities, orbital geometry, smoothing
invariants, and backend consistency rather than frozen oracle values.

### Issue 7: Interpolated apogee range (1 failure)
- **File**: `test_interpolated_apogee_edge_cases.py:88`
- **Error**: `assert -3100015.5 > 600000` — `test_ephemeris_range_returns_valid_values`
- **Cause**: Range computation returns very negative JD instead of valid range
- **Fix**: Fix range computation in interpolated apogee code

### Issue 8: Precision tier file (1 failure)
- **File**: `test_precision_tier.py:200`
- **Error**: `assert 'de441.bsp' == 'de440s.bsp'` — `test_tier_controls_ephemeris_file`
- **Cause**: `set_precision_tier("base")` doesn't switch to `de440s.bsp`
- **Fix**: Check precision tier / ephemeris file selection logic in `state.py`

### Issue 9: Sunshine house 'i' (known limitation)
- Treindl and Makransky algorithms use fundamentally different geometric constructions
- Both 'I' and 'i' currently dispatch to Treindl implementation
- Makransky diverges at lat > 58° when 2/3 × NSA > 90°
- Clean-room implementation requires Makransky's "Primary Directions" (1988) reference
