# Porting eclipse/fixed\_stars/heliacal da Skyfield a LEB

> **Documento storico di implementazione.** Questa analisi descrive il problema
> e il piano osservati nel maggio 2026; il porting è stato poi implementato,
> come indica lo stato seguente. Le metriche, i call site e le fasi future nel
> testo sono fotografie dell'epoca, non l'architettura corrente. Per lo stato
> attuale vedere [Architecture Overview](../development/architecture-overview.md)
> e [LEB Technical Guide](../leb/guide.md).

**Stato:** Superato da 3.0.0rc14 — i consumer vettoriali usano
`LEBVectorEphemeris`; in modalità `leb` non esiste più alcun fallback
Skyfield/JPL. Il fallback descritto nel testo storico riguarda soltanto la
precedente architettura e il routing non sigillato.
**Versione libephemeris:** 1.3.0+
**Data:** 2026-05-08 (analisi), 2026-05-08 (implementazione)
**Autore:** Giacomo Battaglia + Claude

---

## 1. Executive Summary

### Il problema

Tre moduli di libephemeris — `eclipse.py`, `fixed_stars.py`, `heliacal.py` —
bypassano il path LEB e chiamano direttamente l'API Skyfield tramite
`get_planets()`.  Quando vengono invocati, Skyfield carica il DE kernel
(DE441 per il tier extended: 3.1 GB su disco) in memoria.  Questo annulla
parte del risparmio RAM ottenuto con la rimozione di `MADV_WILLNEED`
(v1.3.0) e crea un limite di date per le feature avanzate in modalità
`leb`.

### Dimensione del problema

| Modulo | Righe | `get_planets()` call sites | Metodi Skyfield | Complessità porting |
|--------|------:|---------------------------:|----------------:|---------------------|
| `eclipse.py` | 15.067 | 32 | ~171 | Alta |
| `fixed_stars.py` | 5.937 | 2 | ~10 | Bassa |
| `heliacal.py` | 2.845 | 3 | ~19 | Media |
| **Totale** | **23.849** | **37** | **~200** | |

### Perche il DE kernel si carica

```
API endpoint (eclissi/stelle fisse/eliacali)
  → modulo libephemeris (eclipse.py / fixed_stars.py / heliacal.py)
    → get_planets()                         ← carica DE kernel via Skyfield
      → Skyfield SpiceKernel (DE441: 3.1 GB)
        → eph["earth"], eph["sun"], eph["moon"]
          → .at(t).observe(body).apparent()
            → .altaz(), .radec(), .separation_from(), .position.au
```

In modalità `leb`, `get_planets()` cerca un DE kernel leggero
(`de440s.bsp`, 31 MB, 1849-2150) ma questo limita il range di date
disponibili per le feature avanzate.

### Cosa i file LEB contengono già

I file LEB/LEB2 contengono posizioni pre-calcolate per Sole, Terra, Luna e
tutti i corpi planetari — sono letteralmente il risultato pre-computato dei
calcoli Skyfield.  `calc_ut()` legge questi dati senza toccare Skyfield.

Il porting consiste nel sostituire le chiamate
`eph["earth"].at(t).observe(body).apparent()` con chiamate a
`calc_ut()` + funzioni di trasformazione coordinate che **esistono già**
in libephemeris.

### Fattibilità

| Modulo | Porting possibile? | Effort stimato | Note |
|--------|-------------------|----------------|------|
| `fixed_stars.py` | **Sì, completamente** | 1-2 giorni | Tutti i building block esistono |
| `heliacal.py` | **Sì, completamente** | 2-3 giorni | `azalt()` già esiste |
| `eclipse.py` | **Parziale** | 2-4 settimane | Richiede supporto SEFLG\_XYZ in LEB mode |

---

## 2. Inventario dei building block esistenti

libephemeris ha già implementato quasi tutte le funzioni necessarie per
sostituire Skyfield.  Questa tabella elenca le equivalenze:

### Trasformazioni di coordinate

| Operazione Skyfield | Funzione libephemeris | File:linea | Stato |
|---------------------|----------------------|------------|-------|
| `.altaz()` — eclittica→orizzontale | `azalt(tjdut, ECL2HOR, geopos, atpress, attemp, xin)` | `utils.py:175` | Pronta |
| `.altaz()` — equatoriale→orizzontale | `azalt(tjdut, EQU2HOR, geopos, atpress, attemp, xin)` | `utils.py:175` | Pronta |
| Reverse: orizzontale→equatoriale | `azalt_rev(tjdut, flag, geopos, xin)` | `utils.py:334` | Pronta |
| Rifrazione atmosferica | `refrac(alt, atpress, attemp, flag)` | `utils.py:462` | Pronta |
| Rifrazione estesa | `refrac_extended(alt, alt_geo, atpress, attemp, lapse_rate, flag)` | `utils.py:547` | Pronta |
| Equatoriale↔eclittica | `cotrans_sp(lon, lat, obliquity)` | `utils.py:35` | Pronta |
| Equatoriale↔eclittica (vettore) | `cotrans(xpo, obliquity)` | `utils.py:649` | Pronta |

### Correzioni astrofisiche

| Operazione Skyfield | Funzione libephemeris | File:linea | Stato |
|---------------------|----------------------|------------|-------|
| `.apparent()` — aberrazione | `_apply_aberration(geo, earth_vel, light_time)` | `fast_calc.py:201` | Pronta |
| Precessione J2000↔data | `_precess_ecliptic(lon, lat, jd_tt, direction)` | `fast_calc.py:806` | Pronta |
| Precessione (alternativa) | `_precess_ecliptic(lon, lat, jd_from, jd_to)` | `astrometry.py:365` | Pronta |
| Nutazione dpsi/deps | `get_cached_nutation(jd_tt)` | `cache.py:49` | Pronta |
| Obliquità media/vera | `get_mean_obliquity(jd_tt)` / `get_true_obliquity(jd_tt)` | `cache.py:126` / `cache.py:132` | Pronta |
| Moto proprio stelle | `propagate_proper_motion(ra, dec, pm_ra, pm_dec, parallax, rv, jd)` | `fixed_stars.py:257` | Pronta |
| `.apparent()` — deflessione gravitazionale | — | — | **Non implementata** |

### Posizioni planetarie

| Operazione Skyfield | Funzione libephemeris | File:linea | Stato |
|---------------------|----------------------|------------|-------|
| `eph["sun"].at(t)` | `calc_ut(jd, SUN, FLG_SPEED)` | `planets.py:780` | Pronta (via LEB) |
| `eph["moon"].at(t)` | `calc_ut(jd, MOON, FLG_SPEED)` | `planets.py:780` | Pronta (via LEB) |
| `eph["earth"].at(t)` | `calc_ut(jd, EARTH, FLG_SPEED)` | `planets.py:780` | Pronta (via LEB) |
| Topocentric position | `calc_ut(jd, body, FLG_TOPOCTR)` | — | **KeyError in LEB mode** |
| Cartesian XYZ | `calc_ut(jd, body, FLG_XYZ)` | — | **KeyError in LEB mode** |
| Equatoriale (RA/Dec) | `calc_ut(jd, body, FLG_EQUATORIAL)` | `fast_calc.py:1105` | Pronta |

### Flag attualmente non supportati nel LEB path

`fast_calc.py:1508-1518` definisce i flag che causano **fallback a
Skyfield** in modalità LEB:

```python
if iflag & FLG_TOPOCTR:   raise KeyError(...)   # riga 1509
if iflag & FLG_XYZ:       raise KeyError(...)   # riga 1511
if iflag & FLG_RADIANS:   raise KeyError(...)   # riga 1513
if iflag & FLG_NONUT:     raise KeyError(...)   # riga 1515
if iflag & FLG_ICRS:      raise KeyError(...)   # riga 1517
```

Per il porting completo, sarà necessario implementare almeno
`FLG_TOPOCTR` e `FLG_XYZ` nel path LEB.

---

## 3. Analisi dettagliata: fixed\_stars.py

### Pipeline attuale (Skyfield)

```
fixstar_ut(star_name, jd, iflag)
  → calc_fixed_star_position()                          [riga 3904]
    → _calc_star_position_skyfield()                    [riga 3861]
      → earth = get_planets()["earth"]                  [riga 3897]  ← CARICA DE KERNEL
      → earth_at_t = earth.at(t)
      → _calc_star_position_from_observer(earth_at_t, star_data, ...)  [riga 3808]
        → star = Star(ra_hours, dec_degrees, pm_ra, pm_dec, parallax, rv)
        → astrometric = earth_at_t.observe(star)
        → apparent = astrometric.apparent()              # aberrazione + deflessione
        → ecl = apparent.frame_latlon(ecliptic_frame)    # → lon, lat, dist
    → _apply_fixstar_flags(lon, lat, dist, ...)          [riga 4068]
```

### Operazioni Skyfield e loro sostituti LEB

| Passo | Skyfield | Sostituto LEB |
|-------|----------|---------------|
| 1. Posizione Terra geocentrica | `earth.at(t)` | `calc_ut(jd, EARTH, FLG_SPEED)` |
| 2. Moto proprio J2000→data | `Star(pm_ra, pm_dec)` | `propagate_proper_motion()` (`fixed_stars.py:257`) |
| 3. Parallasse stellare | Interno a `Star()` | Formula: `dist_AU = 206265.0 / parallax_arcsec` |
| 4. Aberrazione annuale | `.apparent()` | `_apply_aberration(geo, earth_vel, lt)` (`fast_calc.py:201`) |
| 5. Deflessione gravitazionale | `.apparent()` | **Non implementata** — effetto < 0.004" per stelle lontane dal Sole, trascurabile per astrologia |
| 6. Equatoriale→eclittica | `.frame_latlon(ecliptic_frame)` | `cotrans_sp(lon, lat, obliquity)` (`utils.py:35`) |
| 7. Precessione | Implicita in `.apparent()` | `_precess_ecliptic()` (`fast_calc.py:806`) |
| 8. Nutazione | Implicita in `.apparent()` | `get_cached_nutation()` (`cache.py:49`) |

### Pipeline LEB proposta

```
fixstar_ut(star_name, jd, iflag)
  → calc_fixed_star_position()
    → _calc_star_position_leb(jd_tt, star_data, iflag)   # NUOVA FUNZIONE
      1. Leggere StarEntry dal LEB: reader.get_star(star_id)
      2. Propagare moto proprio: propagate_proper_motion(ra, dec, pm_ra, pm_dec, parallax, rv, jd)
      3. Equatoriale J2000 → eclittica J2000: cotrans_sp(ra, dec, obliquity_j2000)
      4. Precessione J2000 → data: _precess_ecliptic(lon, lat, jd_tt)
      5. Nutazione: applicare dpsi da get_cached_nutation(jd_tt)
      6. Aberrazione: _apply_aberration(geo, earth_vel, lt)
         - earth_vel da calc_ut(jd, EARTH, FLG_SPEED) → speed components
      7. Return (lon, lat, dist)
    → _apply_fixstar_flags(lon, lat, dist, ...)           # invariata
```

### Gap da colmare

1. **Deflessione gravitazionale** — Effetto massimo ~1.75" per stelle
   vicine al bordo solare, scende a < 0.004" per stelle a > 10° dal Sole.
   Per l'astrologia (precisione ~1') è trascurabile.  Se necessaria, si
   può implementare con la formula di Einstein usando la posizione del Sole
   da `calc_ut()`.

2. **Velocità della Terra** — `_apply_aberration()` richiede il vettore
   velocità della Terra in coordinate cartesiane ICRS.  Il path LEB
   restituisce velocità in eclittica.  Servono:
   - `calc_ut(jd, EARTH, FLG_SPEED)` → (lon, lat, dist, dlon, dlat, ddist)
   - Conversione delle velocità eclittiche in vettore ICRS (o usare
     direttamente il vettore dal LEB reader: `reader.eval_body(EARTH, jd)` 
     restituisce `((x, y, z), (vx, vy, vz))` in ICRS per corpi baricentici)

### Effort stimato: 1-2 giorni

- Creare `_calc_star_position_leb()` (~50-80 righe)
- Aggiungere dispatch in `calc_fixed_star_position()`: se LEB reader
  disponibile, usare path LEB; altrimenti fallback a Skyfield
- Adattare `batch_fixstars_ut()` per usare la nuova funzione
- Test di confronto: output LEB vs output Skyfield, target < 1"

---

## 4. Analisi dettagliata: heliacal.py

### Pipeline attuale (Skyfield)

Le tre funzioni pubbliche (`heliacal_ut`, `heliacal_pheno_ut`,
`vis_limit_mag`) seguono lo stesso pattern:

```
heliacal_ut(jd_start, geopos, body, ...)
  → eph = get_planets()                                 [riga 978]  ← CARICA DE KERNEL
  → observer = wgs84.latlon(lat, lon, alt)
  → observer_at = earth + observer                      # posizione topocentric
  → loop di ricerca (bisection/Newton):
      t = ts.ut1_jd(jd)
      sun_app = observer_at.at(t).observe(sun).apparent()
      body_app = observer_at.at(t).observe(body).apparent()
      moon_app = observer_at.at(t).observe(moon).apparent()
      sun_alt, _, _ = sun_app.altaz()                    # 11 chiamate .altaz()
      body_alt, body_az, _ = body_app.altaz()
      elongation = body_app.separation_from(sun_app)     # 7 chiamate .separation_from()
      → SchaeferModel.is_visible(sun_alt, body_alt, body_mag, moon_alt, moon_phase, ...)
```

### Operazioni Skyfield e loro sostituti LEB

| Passo | Skyfield | Sostituto LEB |
|-------|----------|---------------|
| 1. Posizione Sole | `observer_at.observe(sun).apparent()` | `calc_ut(jd, SUN, FLG_SPEED)` → eclittica |
| 2. Posizione corpo | `observer_at.observe(body).apparent()` | `calc_ut(jd, body_id, FLG_SPEED)` → eclittica |
| 3. Posizione Luna | `observer_at.observe(moon).apparent()` | `calc_ut(jd, MOON, FLG_SPEED)` → eclittica |
| 4. Altitudine/Azimut | `.altaz()` | `azalt(jd, ECL2HOR, geopos, atpress, attemp, (lon, lat, dist))` (`utils.py:175`) |
| 5. Separazione angolare | `.separation_from()` | Formula haversine su coordinate eclittiche |
| 6. Posizione stelle fisse | `Star(ra, dec, pm_ra, pm_dec)` | `fixstar_ut()` (dopo il porting del punto 3) |
| 7. Rifrazione | Custom (Bennet) | Già custom, non dipende da Skyfield |
| 8. Modello atmosferico | `SchaeferModel` | Già custom, non dipende da Skyfield |

### Pipeline LEB proposta

```
heliacal_ut(jd_start, geopos, body, ...)
  → set_topo(lon, lat, alt)                          # imposta osservatore
  → loop di ricerca:
      sun_pos = calc_ut(jd, SUN, FLG_SPEED)     # → (lon, lat, dist, ...)
      body_pos = calc_ut(jd, body_id, FLG_SPEED)
      moon_pos = calc_ut(jd, MOON, FLG_SPEED)
      _, sun_alt, _ = azalt(jd, ECL2HOR, geopos, 0, 0, sun_pos[:3])
      _, body_alt, body_az_app = azalt(jd, ECL2HOR, geopos, 0, 0, body_pos[:3])
      elongation = _angular_separation(sun_pos[0], sun_pos[1], body_pos[0], body_pos[1])
      → SchaeferModel.is_visible(...)                     # invariato
```

### Funzione helper necessaria

```python
def _angular_separation(lon1, lat1, lon2, lat2):
    """Separazione angolare tra due punti in coordinate sferiche (gradi)."""
    import math
    d2r = math.pi / 180.0
    lat1_r, lat2_r = lat1 * d2r, lat2 * d2r
    dlon_r = (lon2 - lon1) * d2r
    cos_sep = (math.sin(lat1_r) * math.sin(lat2_r) +
               math.cos(lat1_r) * math.cos(lat2_r) * math.cos(dlon_r))
    cos_sep = max(-1.0, min(1.0, cos_sep))
    return math.acos(cos_sep) / d2r
```

Questa formula è equivalente a `.separation_from()` di Skyfield.
La differenza è che Skyfield opera in ICRS mentre questa opera in
eclittica.  Per le separazioni angolari la differenza è trascurabile
(stessa sfera, coordinati ruotate).

### Gap da colmare

1. **Stelle fisse per heliacal** — `heliacal_ut()` gestisce anche stelle
   fisse (Sirio, Regolo, etc.) creando un oggetto `Star()` di Skyfield.
   Dopo il porting di fixed\_stars.py, si può usare `fixstar_ut()`.

2. **Coordinate geocentriche** — `heliacal_pheno_ut()` usa anche
   `.radec()` per calcolare l'altitudine geocentrica tramite angolo orario
   (riga 2165-2183).  Questo può essere ottenuto con
   `calc_ut(jd, body, FLG_EQUATORIAL)` che restituisce RA/Dec.

3. **Nota su `azalt()`** — La funzione `azalt()` in `utils.py:175` accetta
   coordinate eclittiche (con `ECL2HOR`) e include la rifrazione
   atmosferica.  È l'equivalente diretto di Skyfield `.altaz()`.

### Effort stimato: 2-3 giorni

- Creare versioni LEB delle tre funzioni (~200 righe)
- Dispatch: se LEB reader disponibile, usare path LEB
- Dipendenza: porting fixed\_stars.py per le stelle fisse in heliacal
- Test di confronto: risultati LEB vs Skyfield, target < 0.01° per altitudini

---

## 5. Analisi dettagliata: eclipse.py

### Panoramica

`eclipse.py` è il modulo più complesso (15.067 righe, 32 call sites).
Le funzioni pubbliche sono:

**Eclissi solari:**
- `sol_eclipse_when_glob()` / `sol_eclipse_when_glob()` — prossima eclissi globale
- `sol_eclipse_when_loc()` / `sol_eclipse_when_loc()` — prossima eclissi locale
- `sol_eclipse_where()` / `sol_eclipse_where()` — dove è visibile
- `sol_eclipse_how()` / `sol_eclipse_how()` — dettagli per location
- `sol_eclipse_how_details()` / `sol_eclipse_how_details()` — dettagli estesi

**Eclissi lunari:**
- `lun_eclipse_when()` / `lun_eclipse_when()` — prossima eclissi lunare
- `lun_eclipse_when_loc()` — eclissi lunare locale

**Occultazioni:**
- Funzioni per occultazioni planetarie e stellari (righe 5700+)

### Pattern di utilizzo Skyfield

I 32 call sites usano 5 pattern distinti:

#### Pattern A — Geocentrico apparente (più comune, ~15 siti)

```python
eph = get_planets()
earth, sun, moon = eph["earth"], eph["sun"], eph["moon"]
earth_at = earth.at(t)
sun_app = earth_at.observe(sun).apparent()
moon_app = earth_at.observe(moon).apparent()
sun_pos = sun_app.position.au        # vettore 3D ICRS (x, y, z)
moon_pos = moon_app.position.au
sun_dist = sun_app.distance().au
```

**Sostituto LEB:** `calc_ut(jd, body, FLG_SPEED)` per distanze e
coordinate eclittiche.  Per i vettori 3D ICRS, serve `FLG_XYZ`
(attualmente non supportato in LEB mode).

#### Pattern B — Topocentric con altaz (~8 siti)

```python
observer = wgs84.latlon(lat, lon, alt)
observer_at = earth + observer
sun_app = observer_at.at(t).observe(sun).apparent()
sun_alt, sun_az, _ = sun_app.altaz()
```

**Sostituto LEB:** `calc_ut(jd, body, FLG_SPEED)` + 
`azalt(jd, ECL2HOR, geopos, 0, 0, (lon, lat, dist))`.  Non serve
`FLG_TOPOCTR` se si usa `azalt()` che accetta coordinate geocentriche.

#### Pattern C — Separazioni angolari (~10 siti)

```python
sep = sun_app.separation_from(moon_app).degrees
```

**Sostituto LEB:** `_angular_separation()` su coordinate eclittiche
(stessa formula haversine di heliacal.py).

#### Pattern D — RA/Dec epoch="date" (~6 siti)

```python
ra, dec, dist = moon_app.radec(epoch="date")
```

**Sostituto LEB:**
`calc_ut(jd, body, FLG_EQUATORIAL)` restituisce RA/Dec (ma in J2000,
non epoch="date").  Per epoch="date" serve applicare la precessione IAU 2006
manualmente: `_precess_ecliptic()` + `cotrans_sp()`.

#### Pattern E — Vettori ICRS 3D per geometria del cono d'ombra (~14 siti)

```python
sun_pos = sun_app.position.au     # (x, y, z) AU
moon_pos = moon_app.position.au
shadow_dir = moon_pos - sun_pos   # vettore direzione ombra
```

**Questo è il pattern critico.** Gli algoritmi besseliani per le eclissi
solari operano su vettori 3D cartesiani ICRS.  Le coordinate eclittiche
(lon/lat/dist) non sono sufficienti — servono i componenti x, y, z.

**Sostituti possibili:**

1. **Implementare SEFLG\_XYZ in LEB mode** — Il LEB reader
   `eval_body()` restituisce già `((x, y, z), (vx, vy, vz))` per corpi
   con `coord_type = COORD_ICRS_BARY`.  Il flag SEFLG\_XYZ in
   `fast_calc.py` potrebbe restituire questi valori direttamente dopo le
   trasformazioni necessarie (precessione, aberrazione).

2. **Conversione eclittica→cartesiana** — Dalle coordinate eclittiche
   (lon, lat, dist) si possono ricavare i componenti cartesiani:
   ```
   x = dist * cos(lat) * cos(lon)
   y = dist * cos(lat) * sin(lon)
   z = dist * sin(lat)
   ```
   E poi ruotare dall'eclittica all'ICRS via matrice di obliquità.
   Meno preciso della lettura diretta ICRS, ma sufficiente per le eclissi.

### Stima di complessità per eclipse.py

| Sotto-task | Call sites | Difficoltà | Dipendenze |
|------------|----------:|------------|------------|
| Pattern B: topocentric altaz | 8 | Bassa | `azalt()` esiste |
| Pattern C: separazioni | 10 | Bassa | Formula haversine |
| Pattern D: RA/Dec epoch="date" | 6 | Media | Precessione + cotrans |
| Pattern A: geocentrico base | 15 | Media | `calc_ut()` |
| Pattern E: vettori 3D ICRS | 14 | **Alta** | Richiede SEFLG\_XYZ in LEB |
| Batch vettorizzati | 2 | Media | Loop scalare |
| Stelle in occultazioni | 4 | Bassa | Dopo porting fixed\_stars |

### Gap da colmare per eclipse.py

1. **SEFLG\_XYZ in LEB mode** — Attualmente causa `KeyError`.  Il LEB
   reader ha già i dati ICRS baricentici; serve implementare la
   trasformazione baricentric→geocentric→XYZ nel path LEB di
   `fast_calc.py`.  Questo è il prerequisito più grande.

2. **SEFLG\_TOPOCTR in LEB mode** — Attualmente causa `KeyError`.
   Necessario per le eclissi locali.  Richiede applicare la parallasse
   diurna alla posizione geocentrica usando la posizione dell'osservatore.

3. **Batch vettorizzato** — `calc_ut()` è scalare.  I 2 siti batch
   (righe ~6118-6165) dovranno usare un loop Python.  Impatto performance
   da valutare, ma i calcoli LEB per punto sono ~100x più veloci di
   Skyfield, quindi il loop potrebbe essere comunque più rapido.

### Effort stimato: 2-4 settimane

- Implementare SEFLG\_XYZ in fast\_calc.py (~200 righe)
- Implementare SEFLG\_TOPOCTR in fast\_calc.py (~150 righe)
- Refactor dei 32 call sites (~500 righe)
- Test di confronto estensivi (eclissi sono sensibili alla precisione)

---

## 6. Lacune critiche da implementare

### 6.1. SEFLG\_XYZ nel path LEB

**Stato attuale:** `fast_calc.py:1511` lancia `KeyError`.

**Cosa serve:** Restituire coordinate cartesiane (x, y, z) geocentriche
invece di (lon, lat, dist).

**Approccio proposto:**

Il LEB reader `eval_body()` restituisce `((x, y, z), (vx, vy, vz))` per
corpi con `coord_type = COORD_ICRS_BARY`.  Questi sono baricentici.
Per ottenere geocentrici:

```
geo_x = body_bary_x - earth_bary_x
geo_y = body_bary_y - earth_bary_y
geo_z = body_bary_z - earth_bary_z
```

Poi applicare precessione e nutazione se richiesto.

**Complessità:** Media.  Il calcolo esiste già nel path Skyfield; serve
replicarlo usando i dati dal LEB.

### 6.2. SEFLG\_TOPOCTR nel path LEB

**Stato attuale:** `fast_calc.py:1509` lancia `KeyError`.

**Cosa serve:** Aggiungere la parallasse diurna (spostamento dell'osservatore
rispetto al centro della Terra) alla posizione geocentrica.

**Approccio proposto:**

1. Ottenere posizione geocentrica eclittica: `fast_calc_ut(reader, jd, body, flags & ~FLG_TOPOCTR)`
2. Convertire in vettore cartesiano
3. Sottrarre il vettore posizione dell'osservatore (da `state._TOPO`)
4. Riconvertire in coordinate eclittiche

**Complessità:** Media.  Le formule sono standard in astrometria.
`state.py` ha già `_TOPO` per la posizione dell'osservatore.

### 6.3. Deflessione gravitazionale

**Stato attuale:** Non implementata nel path LEB.

**Effetto:** ~1.75" al bordo solare, < 0.01" a > 5° dal Sole.

**Per l'astrologia:** Trascurabile.  La precisione richiesta è ~1' (60"),
la deflessione a > 5° dal Sole è < 0.01".

**Per le eclissi:** Potenzialmente significativa per i calcoli di
occultazione precisi (il bordo lunare/solare è a ~0.25° dal corpo).

**Approccio:** Implementare la formula di Einstein per la deflessione
usando la posizione del Sole dal LEB.  ~30 righe di codice.

### 6.4. Separazione angolare

**Stato attuale:** Non esiste una funzione pubblica.

**Cosa serve:** Formula haversine/coseno sferico, ~10 righe:

```python
def _angular_separation(lon1, lat1, lon2, lat2):
    """Separazione angolare in gradi tra due punti in coord. sferiche."""
    d2r = math.pi / 180.0
    lat1_r, lat2_r = lat1 * d2r, lat2 * d2r
    dlon_r = (lon2 - lon1) * d2r
    cos_sep = (math.sin(lat1_r) * math.sin(lat2_r) +
               math.cos(lat1_r) * math.cos(lat2_r) * math.cos(dlon_r))
    return math.acos(max(-1.0, min(1.0, cos_sep))) / d2r
```

---

## 7. Roadmap di implementazione

### Fase 1: fixed\_stars.py (1-2 giorni)

Porting più semplice.  Nessuna dipendenza da flag mancanti.

1. Creare `_calc_star_position_leb()` in `fixed_stars.py`
2. Dispatch in `calc_fixed_star_position()`:
   LEB reader disponibile → path LEB; altrimenti → Skyfield
3. Adattare `batch_fixstars_ut()` per usare il nuovo path
4. Test: confronto output LEB vs Skyfield per ~100 stelle, target < 1"

### Fase 2: heliacal.py (2-3 giorni)

Dipendenza: Fase 1 (per le stelle fisse in heliacal).

1. Creare varianti LEB delle tre funzioni pubbliche
2. Usare `calc_ut()` + `azalt()` + `_angular_separation()`
3. Dispatch: LEB reader disponibile → path LEB
4. Test: confronto heliacal events LEB vs Skyfield, target < 0.01° altitudine

### Fase 3: SEFLG\_XYZ + SEFLG\_TOPOCTR in fast\_calc.py (3-5 giorni)

Prerequisito per Fase 4.

1. Implementare SEFLG\_XYZ: lettura ICRS baricentrici dal LEB →
   sottrazione Terra → precessione/nutazione → output cartesiano
2. Implementare SEFLG\_TOPOCTR: posizione geocentrica → parallasse
   diurna con `_TOPO`
3. Test: confronto con path Skyfield per tutti i pianeti, target < 0.001"

### Fase 4: eclipse.py (1-3 settimane)

Dipendenza: Fase 3.

1. Refactor incrementale: sostituire un pattern alla volta
   - Prima Pattern B (altaz) e C (separazioni) — più semplici
   - Poi Pattern D (RA/Dec) e A (geocentrico base)
   - Infine Pattern E (vettori 3D) — il più critico
2. Mantenere Skyfield fallback durante la transizione
3. Test: confronto eclissi LEB vs Skyfield per date storiche e moderne

### Fase 5: cleanup e rimozione fallback (1-2 giorni)

1. Rimuovere i fallback Skyfield dove non più necessari
2. Verificare che `get_planets()` non venga più chiamato da nessun path
   attivo in modalità LEB
3. Documentazione aggiornata

### Timeline totale stimata: 3-6 settimane

## 8. Rischi e tradeoff

### Precisione

| Effetto | Magnitudine | Impatto astrologia | Impatto eclissi |
|---------|-------------|-------------------|-----------------|
| Deflessione gravitazionale | < 0.01" a > 5° dal Sole | Trascurabile | Significativa near-limb |
| Differenza ICRS vs eclittica per separazioni | < 0.001" | Trascurabile | Trascurabile |
| Aberrazione LEB vs Skyfield | < 0.001" (stessa formula) | Nessuna | Nessuna |
| Precessione IAU 2006 | Identica (stessa libreria erfa) | Nessuna | Nessuna |

### Performance

| Operazione | Skyfield | LEB | Speedup atteso |
|------------|----------|-----|----------------|
| Posizione singolo corpo | ~120 μs | ~1.5 μs | ~80x |
| Star position (con aberrazione) | ~200 μs | ~20 μs (stimato) | ~10x |
| Altaz conversion | Skyfield built-in | `azalt()` pure Python | ~1x (simile) |
| Batch 1000 date | Skyfield vectorized | Loop scalare LEB | Da misurare |

L'unica incognita è il batch: Skyfield usa numpy vectorization, il path
LEB è scalare.  Tuttavia i calcoli LEB per punto sono ~80x più veloci,
quindi un loop di 1000 iterazioni potrebbe essere comparabile o più rapido
(1000 × 1.5 μs = 1.5 ms vs Skyfield batch ~10-50 ms).

### Rischio di regressione

Il rischio principale è nei calcoli eclissi (Fase 4), dove la geometria
del cono d'ombra richiede vettori 3D precisi.  Mitigazione:

- Test di regressione automatici contro output Skyfield
- Mantenere il fallback Skyfield come opzione di debug
- Confronto con effemeridi di riferimento (NASA eclipse predictions)

### Manutenibilità

Il porting elimina una dipendenza runtime dal DE kernel in modalità LEB,
semplificando l'architettura.  Le funzioni LEB usano building block già
testati e mantenuti (stessa base di `fast_calc.py`).

---

## 9. Appendice: inventario completo `get_planets()` call sites

### eclipse.py (32 siti)

| Riga | Funzione | Pattern | Skyfield methods usati |
|-----:|----------|---------|------------------------|
| 494 | `_solar_eclipse_max_time_impl` | A | `.observe().apparent().position.au`, `.distance()` |
| 610 | `_solar_eclipse_max_time_impl` | A | `.observe().apparent().position.au`, `.distance()` |
| 679 | `_solar_eclipse_max_time_impl` | A+C | `.position.au`, `.separation_from()` |
| 1226 | `sol_eclipse_when_glob` | A+C | `.apparent()`, `.separation_from()` |
| 1554 | `_calculate_local_eclipse_phases` | B | `wgs84`, `.altaz()`, `.separation_from()` |
| 1917 | `sol_eclipse_when_loc` | B+C | `wgs84`, `.altaz()`, `.apparent()`, `.separation_from()` |
| 2116 | `sol_eclipse_when_loc` | B+C | `wgs84`, `.altaz()`, `.radec()`, `.separation_from()` |
| 2779 | `sol_eclipse_where` | B+D | `wgs84`, `.radec(epoch="date")`, `.altaz()` |
| 3326 | `sol_eclipse_how` | B+D | `wgs84`, `.altaz()`, `.radec(epoch="date")` |
| 3701 | `sol_eclipse_how_details` | B+D+E | `wgs84`, `.position.au`, `.altaz()`, `.radec()` |
| 4874 | `lun_eclipse_when` | A+E | `.position.au`, `.distance()`, `.apparent()` |
| 5263 | `lun_eclipse_when_loc` | B+E | `wgs84`, `.position.au`, `.altaz()` |
| 5462 | `lun_eclipse_when_loc` | B | `wgs84`, `.altaz()` |
| 5719 | occultazione stellare | A+Star | `Star()`, `.observe()`, `.apparent()` |
| 6415 | occultazione planetaria | B | `wgs84`, `.altaz()`, `.radec()` |
| 7096 | occultazione planetaria loc | B | `wgs84`, `.altaz()`, `.separation_from()` |
| 7586 | occultazione stellare loc | B+Star | `Star()`, `wgs84`, `.altaz()` |
| 8109 | occultazione stellare loc | B+Star | `Star()`, `wgs84`, `.altaz()` |
| 8367 | occultazione details | B | `wgs84`, `.altaz()`, `.separation_from()` |
| 8849 | occultazione where | B+D | `wgs84`, `.radec(epoch="date")`, `.altaz()` |
| 9267 | occultazione how | B+D | `wgs84`, `.radec(epoch="date")`, `.altaz()` |
| 9707 | occultazione how | E | `.position.au` |
| 9867 | occultazione how | E | `.position.au` |
| 10031 | occultazione how | E | `.position.au`, `.separation_from()` |
| 10159 | occultazione how details | B+D+E | `wgs84`, `.position.au`, `.radec()`, `.altaz()` |
| 10316 | occultazione how details | E | `.position.au` |
| 10459 | occultazione how details | E | `.position.au`, `.separation_from()` |
| 12657 | eclissi lunare details | B | `wgs84`, `.altaz()` |
| 13661 | eclissi lunare where | B+D | `wgs84`, `.radec(epoch="date")` |
| 13859 | eclissi lunare how | B+D | `wgs84`, `.radec(epoch="date")`, `.altaz()` |
| 14489 | occultazione generalizzata | B+E | `wgs84`, `.position.au`, `.altaz()` |
| 14852 | occultazione generalizzata | B+Star | `Star()`, `wgs84`, `.altaz()` |

### fixed\_stars.py (2 siti)

| Riga | Funzione | Pattern |
|-----:|----------|---------|
| 3897 | `_calc_star_position_skyfield` | Star + earth.at().observe().apparent().frame_latlon() |
| 4308 | `batch_fixstars_ut` | Star + earth.at().observe().apparent().frame_latlon() |

### heliacal.py (3 siti)

| Riga | Funzione | Pattern |
|-----:|----------|---------|
| 978 | `heliacal_ut` | wgs84 + observe().apparent().altaz() + separation_from() |
| 2087 | `heliacal_pheno_ut` | wgs84 + observe().apparent().altaz() + radec() |
| 2622 | `vis_limit_mag` | wgs84 + observe().apparent().altaz() + separation_from() |
