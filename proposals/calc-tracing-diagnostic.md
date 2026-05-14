# Proposal: Detailed Calculation Tracing & Diagnostic

**Data**: 2026-05-10
**Stato**: Analysis complete, not implemented
**Branch**: feat/skyfield-to-leb-porting

---

## 1. Context e Motivazione

### Il problema

Quando un utente chiama `swe_calc_ut(jd, SE_JUPITER, flags)`, non è visibile:
- Da dove arrivano i dati (LEB? Skyfield? Horizons?)
- Se la correzione COB (center-of-body) usa il file SPK `planet_centers.bsp` o le teorie analitiche delle lune
- Se la nutation/precessione viene dal file LEB o da Skyfield
- Se Skyfield viene importato/toccato durante il calcolo
- Quante volte `eval_body()` viene chiamato e quante sono cache hit
- Quante iterazioni di light-time convergence servono

Il tracing esistente (`tracing.py`) registra solo `{body_id: "LEB"}` — il backend di primo livello. Non dà visibilità sulla catena interna di decisioni.

### Perché serve

1. **Debug di precisione**: capire se un calcolo usa SPK (precisione metri) o analitico (precisione km) per la COB
2. **Analisi dipendenze**: identificare esattamente quando/quando Skyfield viene toccato nella path LEB
3. **Ottimizzazione performance**: misurare cache hit rate, eval_body calls, light-time iterazioni
4. **Pianificazione Zig**: i dati del trace dicono esattamente quali funzioni portare in Zig e quali sono già self-contained
5. **Verifica utente**: dopo installazione, confermare che tutto funziona dalla fonte corretta

---

## 2. Stato Attuale: Cosa Esiste Già

### 2.1 `libephemeris/tracing.py` (47 righe)

Infrastructure ContextVar-based che registra quale backend ha calcolato ogni corpo.

```python
# API pubblica
start_tracing() -> Token          # Attiva tracing
get_trace_results() -> Dict[int, str]  # {body_id: "LEB"} ecc.

# Interna
_record(body_id, source)          # Chiamata da planets.py
```

**Limitazione**: registra solo il source di livello 1 (`"LEB"`, `"Skyfield"`, `"Horizons"`, `"SPK"`, `"ASSIST"`, `"Keplerian"`). Nessun dettaglio interno.

**Overhead**: `ContextVar.get(None)` quando inattivo = ~50ns per call. Zero overhead misurabile.

### 2.2 Logging DEBUG in `planets.py`

```python
get_logger().debug("body=%d jd=%.1f source=LEB", planet, tjdut)
get_logger().debug("body=%d jd=%.1f source=LEB->fallback reason=%s", ...)
get_logger().debug("body=%d jd=%.1f source=Horizons", planet, tjdut)
get_logger().debug("body=%d jd=%.1f source=Skyfield", planet, tjdut)
```

Attivabile con `LIBEPHEMERIS_LOG_LEVEL=DEBUG`. Registra il backend e il motivo di fallback, ma non i dettagli interni (COB, frame data, cache).

### 2.3 `leph diag positions-*` (script esistenti)

`scripts/_tier_diagnostic.py` calcola posizioni per tutti i corpi e mostra la fonte, ma usa una funzione statica `_get_source()` che **approssima** la fonte controllando la copertura SPK — non è un trace runtime.

### 2.4 `leph status` (CLI esistente)

Mostra lo stato dei file (LEB, DE kernel, planet_centers, SPK cache) ma non traccia cosa succede durante un calcolo.

---

## 3. Analisi: Cosa Succede Davvero in Ogni Calcolo

### 3.1 Matrice completa per `swe_calc_ut()` in modalità "auto" con LEB

Per ogni corpo, la pipeline passa attraverso questi stadi decisionali:

#### Corpi ICRS_BARY (Sun, Moon, Mercury-Mars, Earth, asteroidi)

```
swe_calc_ut(jd, body, flags)
  → get_leb_reader() → bundled base_core.leb2 ✓
  → fast_calc_ut()
    → reader.eval_body(body, jd)        ← LEB Chebyshev (cache miss)
    → reader.eval_body(SE_EARTH, jd)    ← LEB Chebyshev (cache miss)
    → Light-time: 3 iterazioni
      → reader.eval_body(body, jd-lt) × 3  ← LEB Chebyshev (cache hit dopo 1a)
    → _apply_gravitational_deflection() ← matematica pura
      → reader.eval_body(SE_SUN, jd)   ← LEB Chebyshev
      → reader.eval_body(SE_JUPITER, jd) ← LEB Chebyshev (attiva COB!)
      → reader.eval_body(SE_SATURN, jd)  ← LEB Chebyshev (attiva COB!)
    → _apply_aberration()               ← matematica pura
    → _frame_data()                     ← LEB nutation (se nel file)
    → Coordinate transforms             ← matematica pura
```

**Skyfield toccato?** NO per il corpo stesso. MA se la deflection chiama Jupiter/Saturn, la COB correction passa per `get_cached_time_tt()` che crea un oggetto Skyfield `Time`.

#### Corpi ICRS_BARY_SYSTEM (Jupiter-Pluto)

```
swe_calc_ut(jd, SE_JUPITER, flags)
  → get_leb_reader() → bundled base_core.leb2 ✓
  → fast_calc_ut()
    → reader.eval_body(SE_JUPITER, jd)  ← LEB Chebyshev (baricentro)
    → reader.eval_body(SE_EARTH, jd)    ← LEB Chebyshev
    → Light-time su baricentro: 3 iterazioni
    → _apply_cob_correction()
      → get_cached_time_tt(jd)          ← Skyfield Time object (inutile)
      → get_planet_center_segment(599)
        → get_planet_centers()
          → cerca planet_centers_base.bsp
          → NON ESISTE → ritorna None
      → SPK path: salta
      → get_cob_offset("jupiter", t)    ← Analitico puro (E5/Meeus)
        → t.tt → float                  ← Qui estrae il JD, Skyfield non serve
    → _apply_gravitational_deflection()
    → _apply_aberration()
    → _frame_data()                     ← LEB nutation
```

**Skyfield toccato?** SÌ — `get_cached_time_tt()` importa e usa Skyfield per creare un `Time` object da cui si estrae solo `.tt` (lo stesso float di partenza). Round-trip inutile.

### 3.2 Tabella riassuntiva: dipendenze per corpo

| Corpo | eval_body | COB | Frame data | Skyfield toccato? |
|-------|-----------|-----|------------|-------------------|
| Sun (0) | LEB | none | LEB nutation | NO |
| Moon (1) | LEB | none | LEB nutation | NO |
| Mercury (2) | LEB | none | LEB nutation | NO* |
| Venus (3) | LEB | none | LEB nutation | NO* |
| Mars (4) | LEB | none | LEB nutation | NO* |
| Jupiter (5) | LEB | analytical (E5) | LEB nutation | YES (via deflection + COB) |
| Saturn (6) | LEB | analytical (TASS) | LEB nutation | YES (via deflection + COB) |
| Uranus (7) | LEB | analytical (Kepl) | LEB nutation | YES (via COB) |
| Neptune (8) | LEB | analytical (Triton) | LEB nutation | YES (via COB) |
| Pluto (9) | LEB | analytical (Charon) | LEB nutation | YES (via COB) |
| Earth (14) | LEB | none | LEB nutation | NO |
| Mean Node (10) | LEB | none | LEB nutation | NO |
| True Node (11) | LEB | none | LEB nutation | NO |
| Mean Apogee (12) | LEB | none | LEB nutation | NO |
| Oscu Apogee (13) | LEB | none | LEB nutation | NO |
| IntpApog (21) | LEB | none | LEB nutation | NO |
| IntpPerig (22) | LEB | none | LEB nutation | NO |
| Chiron (15) | LEB | none | LEB nutation | NO* |
| Ceres (17) | LEB | none | LEB nutation | NO* |
| Pallas (18) | LEB | none | LEB nutation | NO* |
| Juno (19) | LEB | none | LEB nutation | NO* |
| Vesta (20) | LEB | none | LEB nutation | NO* |
| Uranians (40+) | LEB | none | LEB nutation | NO |

*NO* = no di per sé, ma se la deflection gravitazionale calcola Jupiter/Saturn come deflettori, la COB per quei corpi tocca Skyfield indirettamente.

### 3.3 COB: SPK vs Analitico

La correzione COB ha due path con fallback:

```python
# fast_calc.py:303-362
def _apply_cob_correction(pos, ipl, jd_tt):
    t = get_cached_time_tt(jd_tt)          # ← Skyfield sempre

    # Path 1: SPK (alta precisione, metri)
    seg = get_planet_center_segment(naif_id)
    if seg is not None:
        offset = seg.at(t).position.au     # ← Skyfield SPK evaluation
        return pos + offset

    # Path 2: Analitico (fallback, km)
    offset = get_cob_offset(bary_name, t)  # ← Python puro, ma t è Skyfield Time
    return pos + offset                     # usa solo t.tt (un float)
```

| Pianeta | Ampiezza COB | Precisione SPK | Precisione Analitica | Analitico sufficiente? |
|---------|-------------|----------------|---------------------|----------------------|
| Giove | ~742 km | metri | ~100-200 km (~0.05") | Sì per astrologia |
| Saturno | ~760 km | metri | ~50-100 km (~0.05") | Sì per astrologia |
| Urano | ~170 km | metri | ~100-200 km (~0.005") | Sì |
| Nettuno | ~320 km | metri | ~20-50 km (~0.003") | Sì |
| Plutone | ~2130 km | metri | ~1-5 km (~0.01") | Sì |

**Nota**: `planet_centers.bsp` non è mai auto-scaricato. La maggioranza degli utenti usa l'analitico senza saperlo.

### 3.4 Frame Data: LEB vs Skyfield

```python
# fast_calc.py:736-745
def _frame_data(jd_tt):
    if _active_reader_has_nutation:
        return _get_leb_frame_data(reader, jd_tt)   # ← Python puro, no Skyfield
    return _get_skyfield_frame_data(jd_tt)           # ← Skyfield t.M, t._nutation_angles
```

Se il file LEB ha dati di nutation (tutti i file generati li hanno), il path è **Python puro**. Se il file LEB non ha nutation, fallback a Skyfield.

### 3.5 Teorie Analitiche COB — Dettagli

| Pianeta | Teoria | File | Righe | Epoca riferimento | Validità stimata |
|---------|--------|------|-------|-------------------|-----------------|
| Giove | E5/Meeus (Lieske 1998) | `moon_theories/galilean.py` | ~450 | JDE 2443000.5 (~1976) | ~±200 anni |
| Saturno | TASS 1.7 (Vienne & Duriez) | `moon_theories/tass17.py` | ~400 | JD 2444240.0 (1980) | ~±500 anni |
| Urano | Kepleriano 5 lune (GUST86) | `moon_theories/uranian.py` | ~190 | J2000 | Degrada fuori epoca |
| Nettuno | Triton Kepleriano (NEP097) | `moon_theories/triton.py` | ~100 | J2000 | Degrada fuori epoca |
| Plutone | 2-body Charon (PLU060) | `moon_theories/charon.py` | ~100 | J2000 | Stabile |

Tutti i moduli sono **100% Python puro**, zero import Skyfield. L'unico tocco Skyfield è nel caller `fast_calc.py` che crea un oggetto `Time` inutile.

### 3.6 Copertura planet_centers.bsp per Tier

| Tier | File | Dimensione | Copertura |
|------|------|-----------|-----------|
| base | `planet_centers_base.bsp` | 25 MB | 1850-2150 (pieni) |
| medium | `planet_centers_medium.bsp` | 73 MB | Parziale: Giove 1600-2200, Plutone 1800-2200 |
| extended | `planet_centers_extended.bsp` | 223 MB | Parziale: Giove 1600-2200, Saturno -502/+4500, Urano/Nettuno pieni |

**Fuori copertura SPK**: sempre fallback analitico.

---

## 4. Piano di Implementazione

### 4.1 Estendere `tracing.py`

Aggiungere un secondo ContextVar per il detail trace:

```python
# Nuove aggiunte a tracing.py
_trace_detail: ContextVar[Optional[Dict[int, dict]]] = ContextVar(
    "libephemeris_trace_detail", default=None
)

TRACE_DETAIL_ENABLED: bool = os.environ.get("LIBEPHEMERIS_TRACE") == "1"

def start_tracing(detail: bool = False) -> Token:
    token = _trace_data.set({})
    if detail:
        _trace_detail.set({})
    return token

def get_trace_detail() -> Dict[int, dict]:
    return dict(_trace_detail.get() or {})

def _record_detail(body_id: int, **kwargs) -> None:
    """Registra dettagli interni per un body. No-op quando tracing detail è inattivo."""
    d = _trace_detail.get(None)
    if d is None:
        return
    entry = d.setdefault(body_id, {})
    entry.update(kwargs)
```

### 4.2 Strumentare `fast_calc.py`

Punti di instrumentazione (~5 chiamate `_record_detail` per pianeta):

| Punto | Cosa registrare |
|-------|----------------|
| `_fast_calc_core()` inizio | `source="LEB"`, `coord_type`, `leb_file` |
| `_apply_cob_correction()` | `cob="spk"` o `"analytical"`, `cob_theory`, `cob_precision` |
| `_frame_data()` | `frame_data="leb_nutation"` o `"skyfield"` |
| `_pipeline_icrs()` fine | `light_time_iters`, `deflection_applied`, `aberration_applied` |
| `skyfield_used=True/False` | Se qualsiasi chiamata Skyfield è stata fatta |

### 4.3 Strumentare `leb_reader.py`

Un punto in `eval_body()`:

| Punto | Cosa registrare |
|-------|----------------|
| `eval_body()` | `eval_cache_hit=True/False`, incrementare contatori `eval_body_calls` e `eval_cache_hits` |

### 4.4 Comando CLI: `leph diag trace`

Nuovo comando in `dev_cli/cmd_diag.py`:

```bash
leph diag trace [--jd JD] [--date "2000-01-01"] [--now] [--body BODY] [--flags FLAGS]
```

Output di esempio:

```
$ leph diag trace

libephemeris diagnostic trace — JD 2451545.0 (2000-01-01 12:00 TT)

Configuration:
  Calc mode:      auto
  Precision tier: base
  LEB file:       libephemeris/data/leb2/base_core.leb2 (10 MB)
  LEB2 format:    v2 (chunked)
  DE kernel:      not loaded
  Planet centers: not loaded
  Horizons:       not configured

Body Results:
  ID   Body          Source  COB          Frame      Skyfield  Evals  Cache  LT  Defl  Aber
  0    Sun           LEB     none         leb_nut    NO        1      0      1   yes   yes
  1    Moon          LEB     none         leb_nut    NO        2      0      3   no    yes
  2    Mercury       LEB     none         leb_nut    NO        5      1      3   yes   yes
  3    Venus         LEB     none         leb_nut    NO        4      1      3   yes   yes
  4    Mars          LEB     none         leb_nut    NO        4      1      3   yes   yes
  5    Jupiter       LEB     analytical   leb_nut    YES       5      1      3   yes   yes
  6    Saturn        LEB     analytical   leb_nut    YES       5      1      3   yes   yes
  7    Uranus        LEB     analytical   leb_nut    YES       4      1      3   yes   yes
  8    Neptune       LEB     analytical   leb_nut    YES       4      1      3   yes   yes
  9    Pluto         LEB     analytical   leb_nut    YES       4      1      3   yes   yes

COB Details:
  Jupiter:  galilean_e5 (E5/Meeus)     precision ~0.05"
  Saturn:   tass17 (TASS 1.7)          precision ~0.05"
  Uranus:   keplerian (5 moons)        precision ~0.005-0.01"
  Neptune:  triton (NEP097)            precision ~0.003"
  Pluto:    charon_2body (PLU060)       precision ~0.01"

Skyfield Touch Points:
  Jupiter  → get_cached_time_tt()  [COB wrapper]
  Saturn   → get_cached_time_tt()  [COB wrapper]
  ...

Summary:
  Total bodies: 15
  LEB:           15 (100%)
  Skyfield used: 5/15 (33%) — COB Time wrapper only (removable)
```

### 4.5 Environment Variable: `LIBEPHEMERIS_TRACE`

```bash
LIBEPHEMERIS_TRACE=1 python my_script.py
```

Quando attivo:
- Ogni `swe_calc_ut()` registra automaticamente trace detail
- Risultati accessibili via `get_trace_detail()`
- Il booleano `TRACE_DETAIL_ENABLED` è letto una volta all'import

### 4.6 API Python

```python
import libephemeris as swe

token = swe.start_tracing(detail=True)
result = swe.swe_calc_ut(2451545.0, swe.SE_JUPITER, swe.SEFLG_SPEED)
detail = swe.get_trace_detail()
# {5: {"source": "LEB", "cob": "analytical", "cob_theory": "galilean_e5", ...}}
token.var.reset(token)
```

---

## 5. Overhead

| Scenario | Overhead | Note |
|----------|----------|------|
| Tracing OFF (default, produzione) | +0.5% (~250ns) | 5 × ContextVar.get(None) → None → return |
| Tracing ON (debug) | +4% (~2µs) | dict.setdefault() + assignment |
| `LIBEPHEMERIS_TRACE` check | 0ns a runtime | Booleano letto una volta all'import |

---

## 6. File da Modificare

| File | Modifica | Righe nuove (stima) |
|------|----------|---------------------|
| `libephemeris/tracing.py` | Aggiungere `_trace_detail`, `start_tracing(detail=)`, `_record_detail()`, `get_trace_detail()`, `TRACE_DETAIL_ENABLED` | ~40 |
| `libephemeris/fast_calc.py` | Aggiungere `_record_detail()` in 5 punti: `_fast_calc_core`, `_apply_cob_correction`, `_frame_data`, `_pipeline_icrs`, e definizione `_record_detail` locale | ~35 |
| `libephemeris/leb_reader.py` | Aggiungere `_record_detail()` per cache hit/miss in `eval_body` | ~10 |
| `libephemeris/dev_cli/cmd_diag.py` | Nuovo comando `leph diag trace` | ~120 |
| `libephemeris/__init__.py` | Esporre `get_trace_detail()` nella public API | ~2 |

**Totale stimato**: ~200 righe nuove.

---

## 7. Verifica

1. `leph diag trace` — verifica output tabellare come nell'esempio
2. `leph diag trace --now --body jupiter` — dettaglio per singolo corpo
3. `LIBEPHEMERIS_TRACE=1 python -c "import libephemeris as s; s.swe_calc_ut(2451545.0, 5, 0); print(s.get_trace_detail())"` — env var
4. Con `planet_centers.bsp` presente: COB cambia da "analytical" a "spk"
5. `LIBEPHEMERIS_MODE=skyfield leph diag trace` — tutti corpi mostrano Source=Skyfield
6. `pytest tests/test_tracing.py` — test esistenti non si rompono
7. Benchmark: `pytest tests/test_performance_benchmarks.py -v` — nessun regressione >1%

---

## 8. Implicazione per il Futuro: Rimozione Dipendenza Skyfield

Il trace detail rivela che Skyfield è toccato nella path LEB **solo** per:
1. `get_cached_time_tt(jd_tt)` in `_apply_cob_correction` — crea un `Time` object per estrarne `.tt` (lo stesso float)
2. `get_cached_time_tt(jd_tt)` in `_get_skyfield_frame_data` — fallback solo se LEB non ha nutation (mai)

La fix per eliminare Skyfield dal runtime LEB è banale:

```python
# In fast_calc.py, _apply_cob_correction:
# PRIMA:
t = get_cached_time_tt(jd_tt)
offset = get_cob_offset(bary_name, t)

# DOPO:
offset = get_cob_offset(bary_name, jd_tt)  # passare il float direttamente
```

E in `moon_theories/constants.py`, cambiare `get_cob_offset` per accettare `float`:

```python
# PRIMA:
def get_cob_offset(planet_name: str, t) -> Tuple[float, float, float]:
    jd = t.tt  # <-- unico uso di t

# DOPO:
def get_cob_offset(planet_name: str, jd: float) -> Tuple[float, float, float]:
    # rimuovere jd = t.tt
```

Questo renderebbe la path LEB 100% Skyfield-independent per tutti i 30 corpi. Skyfield resterebbe dipendenza per generazione e fallback.

---

## Appendice A: Mappa Completa delle Fonti Dati

### File distribuiti con il pacchetto

| File | Dimensione | Contenuto | Richiede Skyfield? |
|------|-----------|-----------|-------------------|
| `libephemeris/data/leb2/base_core.leb2` | 10 MB | Chebyshev per Sun-Pluto, Earth, nodi, apogei | No |
| Skyfield (pip dependency) | ~30 MB | DE440 loader, timescale, nutation models | Sì |

### File scaricabili opzionali

| File | Dimensione | Contenuto | Richiede Skyfield? |
|------|-----------|-----------|-------------------|
| `ephemeris_base.leb` / `*_core.leb2` | 10-53 MB | LEB completo per tier | No |
| `planet_centers_base.bsp` | 25 MB | COB offsets SPK | Sì (per leggere .bsp) |
| `de440s.bsp` | 32 MB | JPL ephemeris (tier base) | Sì |
| `de440.bsp` | 114 MB | JPL ephemeris (tier medium) | Sì |
| `de441.bsp` | 3.1 GB | JPL ephemeris (tier extended) | Sì |

### Dati generati da Skyfield (offline, task sviluppatore)

- File LEB: generati da `scripts/generate_leb.py` usando Skyfield + numpy
- File planet_centers.bsp: generati da `scripts/generate_planet_centers_spk.py` usando spiceypy + JPL SPK sources

## Appendice B: Funzioni Coinvolte nel Dispatch

```
planets.py
  swe_calc_ut()                      # Entry point UT
    → state.get_leb_reader()          # Singleton LEB reader
    → fast_calc.fast_calc_ut()        # LEB path
    → state.get_horizons_client()     # Horizons path
    → _calc_body()                    # Skyfield path

fast_calc.py
  fast_calc_ut() → fast_calc_tt() → _fast_calc_core()
    → _pipeline_icrs() / _pipeline_ecliptic() / _pipeline_helio()
      → reader.eval_body()            [leb_reader.py]
      → _apply_cob_correction()       [fast_calc.py:303]
        → get_planet_center_segment() [state.py:1022] — SPK path
        → get_cob_offset()            [moon_theories/constants.py] — analitico
      → _apply_gravitational_deflection() [fast_calc.py:365]
      → _apply_aberration()           [fast_calc.py:201]
      → _frame_data()                 [fast_calc.py:736]
        → _get_leb_frame_data()       [fast_calc.py:669] — Python puro
        → _get_skyfield_frame_data()  [fast_calc.py:467] — Skyfield fallback

state.py
  get_calc_mode()                     # Risolve: programmatico → env → TOML → "auto"
  get_precision_tier()                # Risolve: programmatico → env → TOML → "medium"
  get_leb_reader()                    # Lazy load + auto-discover + auto-download
  get_planet_centers()                # Lazy load planet_centers_*.bsp
  get_planet_center_segment(naif_id)  # Cerca segmento per NAIF ID
  get_planets()                       # Lazy load DE kernel

tracing.py
  start_tracing()                     # Attiva accumulatore
  get_trace_results()                 # {body_id: source}
  _record(body_id, source)            # Chiamata da planets.py
```
