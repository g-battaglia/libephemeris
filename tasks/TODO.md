# LibEphemeris — Piano di Verifica Esaustiva

## Contesto Generale del Progetto

LibEphemeris e' una reimplementazione clean-room in puro Python della Swiss Ephemeris,
la libreria standard de facto per calcoli astronomici ed astrologici. Usa esclusivamente
dati NASA JPL (DE440/DE441) tramite Skyfield, anziche' i file di effemeridi proprietari
della Swiss Ephemeris. L'obiettivo e' la compatibilita' 1:1 con l'API PySwissEph
mantenendo una precisione sub-arcsecond rispetto ai dati JPL.

---

## Open Issues from poe test:full (130 failures / 17376 passed)

Last updated: 2026-03-27

### Summary by Root Cause

| # | Root Cause | Failures | Priority | Status |
|---|-----------|----------|----------|--------|
| 1 | `get_ayanamsa_ut` returns numpy.float64, not native float | 2 | HIGH | PENDING |
| 2 | Eclipse local API signatures (`lun_eclipse_when_loc`, `sol_eclipse_when_loc`) | 8 | HIGH | PENDING |
| 3 | LEB files stale — InterpApogee/InterpPerigee encoded with pre-fix values | ~60 | MED | PENDING (regen) |
| 4 | Extended tier asteroids — "Invalid Time to evaluate" (SPK date range) | ~40 | MED | PENDING |
| 5 | TrueNode crosstier at boundary JD 2290867.5 (1560) — EphemerisRangeError | 2 | LOW | PENDING |
| 6 | Lunar ELP2000 perturbation consistency after BUG-001 fix | 2 | LOW | PENDING |
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

### Issue 3: LEB stale InterpApogee/InterpPerigee (~60 failures)
- **Files**: All `test_leb/compare/` dirs for bodies 21 (InterpApogee), 22 (InterpPerigee)
- **Error**: Chebyshev polynomials encode pre-BUG-001 values; errors 500-10000+ arcseconds
- **Fix**: Regenerate LEB files: `poe leb:generate:base`, `poe leb:generate:medium`, `poe leb:generate:extended`
- **Note**: Data regeneration, not a code bug

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

### Issue 6: Lunar perturbation consistency (2 failures)
- **Files**: `test_elp2000_apogee_perturbations.py:137`, `test_elp2000_perigee_perturbations.py:238`
- **Error**: Decomposition mismatch — BUG-001 changed coefficients but test expectations weren't updated
- **Fix**: Update test expected values to match new perturbation coefficients

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
