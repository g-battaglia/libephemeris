# Delta T (ΔT) in libephemeris

This note documents how libephemeris computes **ΔT = TT − UT1**, the slowly
varying, partly unpredictable offset between uniform atomic time (TT) and
clock-from-Earth-rotation time (UT1). ΔT is required by every UT-based
calculation: it maps a civil (UT) instant to the dynamical time at which the
ephemeris and precession/nutation are evaluated. A wrong ΔT shifts **all** body
positions (∝ each body's mean motion) and, through sidereal time, the **house
cusps / Ascendant**.

Implementation: `libephemeris/time_utils.py` (`deltat`, `_deltat_espenak_meeus`)
and the model selector in `libephemeris/state.py`
(`set_delta_t_model` / `get_delta_t_model`).

> **In short.** libephemeris uses a **multi-era (piecewise)** ΔT model — IERS
> observed values for modern dates, Stephenson-Morrison-Hohenkerk 2016 for the
> eclipse era, and a long-term tidal model at the extremes — so a single, correct
> ΔT drives positions *and* house angles across the whole supported range. A
> second, fully self-contained model (the classic NASA Espenak-Meeus 2006
> polynomials) is selectable.

---

## 1. Why ΔT must be piecewise

ΔT is not given by one closed-form expression over all of history. Earth's
rotation is *measured* by atomic clocks since 1955, *reconstructed* from
telescopic observations back to ~1620, *inferred* from ancient eclipse/occultation
records back to ~720 BCE, and *unknowable* in the future (it depends on
unmodelled core-mantle coupling, post-glacial rebound, ice-mass redistribution,
etc.). A single polynomial applied outside its validated window diverges and
produces gross errors — for an Ascendant this can mean being off by a degree or
more. The robust solution is a **block architecture**, where each era uses the
source best constrained by data for that period:

| Era | Method | libephemeris source of truth |
|---|---|---|
| 1955 → present (+ short prediction) | IERS / USNO observed (atomic clocks) | Skyfield ΔT table (IERS finals) |
| ~1620 → 1955 | Telescopic-era tabulated values | Skyfield ΔT table / SMH-2016 |
| −720 → ~1620 | Stephenson-Morrison-Hohenkerk 2016 spline (eclipse records) | SMH-2016 |
| before −720 / far future | Long-term tidal model (integrated length-of-day) | SMH-2016 long-term tail |

libephemeris obtains this piecewise model from **Skyfield**, whose `delta_t`
implements Stephenson, Morrison & Hohenkerk (2016) blended with IERS observations
and a long-term tail — the modern scientific standard for ΔT reconstruction.

## 2. The multi-era model in numbers

ΔT in seconds, libephemeris default model, alongside the observed/consensus value
and — for illustration — the single long-term parabola `ΔT = −20 + 32·u²`
(u = (year−1820)/100) that a non-piecewise model would produce:

| year | libephemeris ΔT | observed / consensus | single parabola | active block |
|---|---|---|---|---|
| 2000 | 63.8 | 63.83 (IERS) | 83.7 | IERS observed |
| 2025 | 69.1 | ~69 (IERS) | 114.5 | IERS observed |
| 1900 | −2.7 | −2.8 | 0.5 | telescopic table |
| 0 | 10441 | ~10580 (SMH-2016) | 10580 | SMH-2016 spline |
| −500 | 16939 | ~16800 (SMH-2016) | 17204 | SMH-2016 spline |
| −6000 / +10000 | 198671 / 216871 | — | 195668 / 214100 | long-term tidal tail |

In the well-observed era libephemeris tracks the IERS value (e.g. 63.8 s at 2000,
the observed value, not the 83.7 s a single parabola would give). The block
boundaries are continuous — there are no steps across the 1600 → 2000 range.

## 3. The ΔT model selector

ΔT resolution order in `deltat()` is:

1. **user-defined** — `set_delta_t_userdef(value_in_days)` (highest priority);
2. **IERS observed** — if enabled via `set_iers_delta_t_enabled(True)`;
3. **the selected model** (this section).

The model is chosen with `set_delta_t_model(...)` / the
`LIBEPHEMERIS_DELTAT_MODEL` environment variable / the `deltat_model` TOML key:

| value | model | notes |
|---|---|---|
| `smh2016` *(default)* | Stephenson-Morrison-Hohenkerk 2016 via Skyfield | the modern scientific standard; piecewise as in §1 |
| `espenak_meeus` | classic NASA **Espenak & Meeus 2006** polynomial set | self-contained implementation; the multi-segment polynomials (−1999…+3000) used by NASA's eclipse predictions |

```python
import libephemeris as eph
eph.set_delta_t_model("espenak_meeus")   # or "smh2016"
eph.get_delta_t_model()                   # -> "espenak_meeus"
# or, by environment:  LIBEPHEMERIS_DELTAT_MODEL=espenak_meeus
```

Both models are **clean-room**: every coefficient comes from published sources (the
SMH-2016 / Skyfield implementation, and the NASA Espenak-Meeus polynomials). They
agree closely in the well-observed range and differ mainly at the extremes (deep
past / far future), where ΔT is itself uncertain. The `espenak_meeus` model
reproduces the classic NASA polynomials and exists for users and tools that
specifically expect them; the default `smh2016` is the recommended choice.

> **Scope — which calculations the selector (and `userdef` / IERS) affect.**
> All three ΔT priorities feed `deltat()` / `deltat_ex()`, and hence the TT used by
> the **LEB / fast / Horizons** position backends (the default is **LEB**), as well
> as eclipses, heliacal events and long-term sidereal time. The **`"skyfield"`**
> position backend is the exception: it derives TT from Skyfield's *own* internal
> ΔT model (SMH-2016 + IERS), so its planetary positions ignore
> `set_delta_t_model()`, `set_delta_t_userdef()` and the IERS toggle. Concretely,
> selecting `espenak_meeus` (or forcing a `userdef` value) shifts positions in the
> default LEB backend — e.g. a `userdef` step of 0.1 day moves the Sun by ~364″ —
> but leaves `"skyfield"`-mode positions unchanged. This is a pre-existing property
> of all three overrides, not specific to the model selector.

## 4. Independence

libephemeris is an independent implementation and imports no third-party ephemeris
engine at runtime; both ΔT models above are clean-room (published coefficients
only). When exact parity with another engine's ΔT is required for a controlled
comparison, that engine's ΔT can be injected from the *outside* via
`set_delta_t_userdef`, isolating the pure ephemeris-model difference from the ΔT
choice.

## 5. Validation

The ΔT model and its consequences are validated in the `validation/` repo:

- a dedicated **matched-ΔT sweep** that forces an external reference ΔT into
  libephemeris via `set_delta_t_userdef` and confirms positions match the reference
  to < 0.1″ (Moon to ~0.001″) across the supported range, isolating the pure
  ephemeris-model difference from ΔT;
- the position/house long-term rounds (`run.sh positions` / `run.sh houses`) bound
  the (expected, physical) ΔT-driven divergence at remote epochs and confirm the
  precession/sidereal-time model is an exact reproduction of the published physics
  at matched ΔT.

## 6. References

1. Stephenson, F.R., Morrison, L.V. & Hohenkerk, C.Y. (2016). *Measurement of the
   Earth's rotation: 720 BC to AD 2015.* Proc. Royal Society A, 472, 20160404.
2. Morrison, L.V., Stephenson, F.R., Hohenkerk, C.Y. & Zawilski, M. (2021).
   *Addendum 2020 to "Measurement of the Earth's rotation: 720 BC to AD 2015".*
   Proc. Royal Society A, 477, 20200776.
3. Espenak, F. & Meeus, J. (2006). *Five Millennium Canon of Solar Eclipses:
   −1999 to +3000.* NASA/TP-2006-214141 (polynomial ΔT expressions).
4. Morrison, L.V. & Stephenson, F.R. (2004). *Historical values of the Earth's
   clock error ΔT and the calculation of eclipses.* J. Hist. Astron. 35, 327.
5. IERS Earth Orientation Centre — observed ΔT / UT1 bulletins.

---

*For a head-to-head comparison of ΔT (and its effect on positions) with Swiss
Ephemeris, see [Swiss Ephemeris Comparison](../comparison/index.md).*
