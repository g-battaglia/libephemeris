# NEXT — Roadmap beyond v3

This document records the work planned after the v3.0.0 release: the path to
the highest achievable precision and coverage across every surface, including
the niche ones. It is a living plan, not a commitment to dates.

Items that shipped code in v3 already point here from the source ("planned for
a future release; see NEXT.md"): notably `nod_aps` for minor bodies (§Wave 1.3)
and the composite `SIDBIT` sidereal projections (§Wave 3.1). Those guards
currently raise or warn rather than return a plausible-looking wrong value; the
waves below turn them into full features.

---

## 1. Where we start

State at the v3.0.0 release:

- **Parity vs the reference engine: 89 GREEN + 22 CERTIFIED = 111/111, 0 RED.**
  The 22 CERTIFIED components are all deliberate, documented scientific
  divergences (see `docs/comparison/intentional-divergences.md`), not defects.
- **Unit suite: green** on both backends (Skyfield and the precomputed LEB
  Chebyshev backend).
- **Comparison suite: clean.**

### Binding philosophy

**We do not copy the reference: we stay independent, ideally more precise.**

- 1:1 parity applies to **behaviour, API and semantics**.
- Where we are more scientific / better aligned with modern astronomical
  standards (ΔT from IERS + a peer-reviewed long-term model, a continuous GMST
  join, true cusp derivatives, modern star catalogues, true rates in the speed
  slots) we **do not replicate the reference's data or its defects**.
- Precision is demonstrated against **independent arbiters** (JPL Horizons,
  IERS, astropy/ERFA, peer-reviewed literature), never against the reference
  itself.
- Published models (e.g. Schaefer 1993 for heliacal visibility, historical
  Ptolemaic/Babylonian criteria) are implemented **from the literature**.

---

## 2. Known scientific divergences (the state we defend, not "bugs to fix")

Documented in `docs/comparison/intentional-divergences.md` and
`docs/comparison/known-differences.md`. Summary:

| Topic | Divergence | Why we are the correct one |
|---|---|---|
| Modern ΔT | ±100–144 ms in 2024 | We use the real IERS signal; black-box reference output follows a different interpolation/extrapolation curve |
| Historic/future ΔT | up to 20 s @1582, 875 s @3000 | We use a newer peer-reviewed long-term model |
| GMST / houses from 2050 | ~2" (up to 9" @2340 via ΔT) | Black-box reference output has a ~1.9" discontinuity at 2050-01-01; we join continuously |
| Placidus/Koch cusp speed | up to 2392"/day | We report the true derivative |
| MEAN apside lat speed | different slot | The black-box slot equals the position at t+1d; we store the true rate |
| MEAN planetary apside speed | ~0.15°/day (speed channels only) | The reference's convention is unidentified; positions and nodes match exactly |
| Planetary nodes (scattered days) | transient spikes up to 3.5" | Isolated spikes in black-box reference output; internal cause unknown |
| Star radial-velocity slots | ~250 "AU" on Spica | We use a modern catalogue; angular positions are unaffected |
| Minor-body historic orbit solutions | 80–130" (Nessus/Amor/Bennu, 1900–1970) | We follow the SPK/Horizons solution; Horizons confirms us to ≤0.25" |
| Near-polar rise/set search window | Moon events at lat ±80–85° | The reference searches ~1 day and reports none; we search 2 |
| Heliacal visibility model | 3–14 days, apparition selection | Margin calibration; the physical rebuild lands in Wave 4 |

---

## 3. The six waves

### Wave 1 — Minor bodies: total coverage
*Goal: any numbered asteroid is first-class, surpassing the reference (which
needs a per-body file).*

- **1.1 Cartesian-state helper** `calc_minor_body_state(ipl, jd_tt) -> (r, v)`
  heliocentric ecliptic-J2000, with the chain SPK type21 → SPK type2/3 →
  ASSIST n-body → Keplerian elements. Today the public surface returns only
  spherical coordinates.
- **1.2 Universal ASSIST**: generalise the curated-asteroid ASSIST path to any
  small body via the SBDB element fetch; add an integration cache and an on-disk
  SBDB cache with a TTL.
- **1.3 `nod_aps` for minor bodies** (resolves the guard in `planets.py`):
  drop the raise-guard for asteroids/centaurs, source the state from 1.1,
  μ = GM_Sun, keep the rest of the machinery (aberration, J2000 reduction);
  MEAN falls back to osculating as Pluto already does. Validate on a grid of
  five bodies × three epochs × methods.
- **1.4 Universal `pheno`**: H, G from SBDB for bodies outside the curated
  table; diameter from SBDB when available.
- **1.5 Edge contract**: with universal ASSIST the out-of-SPK-coverage region
  becomes accurate. Enforce with a `strict_precision` opt-in that raises when
  neither SPK nor ASSIST is available, with an explicit opt-out for the
  degraded Keplerian fallback. Publish a source→precision-by-epoch matrix.

### Wave 2 — Planetary moons out of the box
*The infrastructure already exists (NAIF id parity, the full light-time /
aberration / sidereal pipeline, SPK registration with auto-detect).*

- **2.1 NAIF kernel auto-download** keyed by parent planet, modelled on the
  existing data downloaders, hooked into the moon dispatch.
- **2.2 `pheno` for moons**: a radius table for the ~20 major moons plus a
  dedicated branch (exact phase/elongation/diameter; best-effort magnitude from
  published models, documented as such).
- **2.3 Missing names/constants** for unmapped ids and planet-center ids.
- **2.4 Validation** against **Horizons** (an independent arbiter — the
  reference needs satellite files nobody ships): six moons × three epochs,
  target < 1" geocentric.

### Wave 3 — SIDBIT, rise_trans, remaining conventions

- **3.1 `SIDBIT`** (resolves the composite-mode guard in `set_sid_mode`): today
  the projection bits are accepted with a warning and the base mode is used.
  Implement the `ECL_T0` / `ECL_DATE` / `SSY_PLANE` / `PREC_ORIG` projections
  from their geometric definition, after a multi-date black-box characterisation
  (bit × mode × epoch).
- **3.2 `rise_trans`: full BIT_HINDU_RISING** — only the geocentric-no-ecliptic-
  latitude component is missing; the rest of the rsmi bits are already handled.
- **3.3 `FLG_JPLHOR` / `JPLHOR_APPROX` policy**: bits accepted and documented as
  no-ops (our positions are already DE440/ICRS-native); a true frame-fix
  emulation is optional and out of scope.
- **3.4 Close the MEAN planetary apside speed investigation**: one final
  time-boxed black-box session; if the convention does not emerge, keep the
  documented divergence and mark the channel EXPECTED in the validator.
- **3.5 Per-flag certification of minor surfaces**: micro-rounds for per-bit
  `rise_trans`, `refrac`/`refrac_extended` on a P/T/altitude grid, `azalt`
  round-trip, `gauquelin_sector` 0–5.

### Wave 4 — Heliacal: unified physics (the big build)

Today `heliacal.py` carries an empirical `SchaeferModel` with hard-coded margins.
The real physics already exists but is disconnected (`extinction.py`,
`schaefer.py`), the AVKIND criteria are ignored, and the search skeleton is
duplicated between the two backends. Target architecture:

```
extinction.py (physics source) → visibility.py NEW (vectorised VisibilityEngine)
→ heliacal_search.py NEW (single search skeleton, LEB/Skyfield PositionProvider)
→ heliacal.py (API + dispatch) → schaefer.py (deprecated shim)
```

Phased so each phase lands with a green suite: black-box characterisation +
golden captures; the engine in bit-identical legacy mode; physical
`vis_limit_mag` behind a transition flag; the unified search (twilight session
as the unit, boundary bisection, coherent jd1/jd2/jd3, conjunction tracking);
test migration behind a model flag; default flip; final validation against a
10-body × 6-location × 4-season grid and the literature (Schaefer 1993,
Krisciunas & Schaefer 1991, Yallop 1997, the Venus tablets of Ammisaduqa).

### Wave 5 — Certification: dual verdict

- **5.1 Dual-verdict validator**: alongside parity-vs-reference, an
  accuracy-vs-arbiter column (Horizons, astropy/ERFA, IERS, literature).
- **5.2 Divergence charter**: promote `intentional-divergences.md` to a
  versioned user document, with numeric arbiter evidence per entry.
- **5.3** — done in v3: the CERTIFIED verdict in the validation harness.

### Wave 6 — Standard-grade quality

- **6.1 Zero-setup**: `pip install` and it works to ±13000 years; unified
  downloaders; a `leph doctor` command that verifies/fetches everything and
  prints the precision-by-body-class matrix.
- **6.2 Continuous certification**: a CI validator profile plus a published
  per-release matrix. "Verified: N GREEN + M CERTIFIED / 111" badge.
- **6.3 Batch API**: `calc_ut_many(jds, planet, flags) -> ndarray (N, 6)`
  vectorising the LEB branch only, scalar fallback for Horizons/Skyfield.
  Target > 10× on the equivalent loop.
- **6.4 Hardening**: guard case for the Bennu-1980 boundary; cache hygiene
  (add the long-term sidereal LRU caches and the star-catalogue rebuild to
  `clear_caches()`).

---

## 4. Suggested order

Wave 1 → 2 → 3 → 5 → 4 (the long build) → 6, with 6.1 startable in parallel.

Each wave verifies with new targeted tests plus green `test:skyfield:fast` and
`test:leb:fast`, dedicated rounds/verify in the validation repo, a `./run.sh all`
at the end of the wave, and clean lint/format/typecheck on the touched files.
