# Handover — Indagine sui 16 fallimenti dei validation test

**Branch:** `fix/v3.0.0-review-round2`
**Data indagine:** 2026-07-04
**Stato:** SOLO INDAGINE. Nessuna implementazione/fix applicato. Nessun file di codice modificato.
**Scopo di questo doc:** passare a un altro agente il recap completo della diagnosi, così può decidere/implementare i fix.

---

## 0. Contesto di partenza

Una run dei validation test (repo esterno `validation/`, ~20k test vs pyswisseph) ha prodotto **16 failed, 19939 passed, 179 skipped, 14 xfailed, 29 xpassed**.

I 16 fallimenti si raggruppano in **4 aree indipendenti**. Conclusione generale:
**13 su 16 sono falsi positivi dei test di validation** (da correggere i test, NON la libreria); **3 su 16 sono un problema di prestazioni reale della libreria** (`planet_occult_when_glob`), non di correttezza.

Vincolo di progetto assoluto (da `CLAUDE.md`): **compatibilità 1:1 con pyswisseph**. Questo è il criterio dirimente per stabilire "chi ha ragione" tra test e libreria.

I 16 test falliti:
```
test_compare_minor_bodies.py::TestNearEarthAsteroids::test_nea_position[...111955-Bennu]  (x3)
test_fixed_stars_and_points.py::TestAngles::test_angles_require_location
test_fixed_stars_and_points.py::TestArabicParts::test_arabic_parts_require_precalc
test_hypothetical.py::TestWaldemathCalculations::* (x8)
test_deep_validation_4.py::TestPlanetOccultWhenGlob::test_return_structure          (timeout)
test_deep_validation_4.py::TestPlanetOccultWhenGlob::test_does_not_crash_with_various_planets  (timeout)
test_deep_validation_4.py::TestPlanetOccultWhenLoc::test_return_structure           (timeout)
```

---

## 0-bis. Scoperta trasversale: i test girano in modalità SKYFIELD, non LEB

Domanda emersa dall'utente: "i test non dovrebbero girare tutti in LEB di default?"

**Risposta: NO, di default girano in Skyfield.**

- `validation/compare_scripts/tests/conftest.py:24`
  `_COMPARE_MODE = os.environ.get("LIBEPHEMERIS_COMPARE_MODE", "skyfield").lower()` → **default `"skyfield"`**.
- ramo `else` (`conftest.py:33-35`) forza `LIBEPHEMERIS_MODE = "skyfield"`.
- Per girare in LEB serve esportare esplicitamente `LIBEPHEMERIS_COMPARE_MODE=leb` (allora `conftest.py:42-46` fa `set_calc_mode("auto")` e lascia intercettare LEB).
- `validation/run.sh:211` lancia `pytest compare_scripts/tests/` **senza** esportare `LIBEPHEMERIS_COMPARE_MODE` → eredita il default Skyfield.
- I traceback dei timeout confermano il path Skyfield (`skyfield/positionlib.py`, `jplephem/spk.py`), non il reader LEB.
- Quasi tutti gli script standalone `validation/compare_scripts/rounds/*.py` hanno `LIBEPHEMERIS_MODE=skyfield` **hardcoded** (unico outlier: `swiss_deltat_isolation.py` in LEB).

**IMPORTANTE:** passare il default a LEB **NON chiude i 16 fail**. Aiuta solo (parzialmente) il gruppo timeout. Vedi §4.

I file `.leb` esistono e sono presenti:
- `validation/data/leb/ephemeris_{base,medium,extended}.leb` (+ companion `.leb2`)
- copia su disco esterno `/run/media/giacomo/ExFAT/libephemeris/leb/`

---

## 1. Bennu / NEA — 3 fail — **PROBLEMA DEI TEST (setup) + limite by-design della lib**

### Sintomo
`test_compare_minor_bodies.py::TestNearEarthAsteroids::test_nea_position` per Bennu (body id **111955**) ritorna `source=Keplerian (fallback)` con diff di longitudine **4.77° / 5.46° / 44.87°** (tolleranza `NEA_TOL = 0.05°`). Warning da `planets.py:2974`.

### Causa radice
Nell'ambiente di test (modalità Skyfield pura), **ogni** ramo ad alta precisione per Bennu è saltato, e resta solo il Kepleriano:

| Ramo | Esito per Bennu | Riferimento |
|---|---|---|
| SPK type21 registrato | No (nessun kernel) | `planets.py:2869` |
| Auto-download SPK | **Bloccato** — Bennu è in `SPK_AUTO_DOWNLOAD_BLOCKED` | `constants.py:408-414`; gate `planets.py:2874` |
| ASSIST n-body | **Saltato** — `rebound`/`assist` NON installati nel venv validation | `planets.py:2923-2925`; `rebound_integration.py:610-621` |
| Keplerian (last resort) | **Eseguito** → warning | `planets.py:2974`; elementi `minor_bodies.py:1971-1981` |

- `conftest.py:24-35` forza Skyfield (no Horizons, no LEB).
- `conftest.py:68-70` `spk_auto.enable_common_bodies()` copre solo Chiron/Pholus/Ceres/Pallas/Juno/Vesta/Eris/Sedna — **NON Bennu**.
- Gli altri 8 NEA (Eros, Amor, Apophis, Itokawa, Ryugu, Toro, Toutatis, Icarus) **passano** perché non sono bloccati → auto-download SPK li copre.
- Bennu è l'**unico** NEA in `SPK_AUTO_DOWNLOAD_BLOCKED`, perché **JPL Horizons rifiuta di generare SPK per Bennu** (mission target OSIRIS-REx con kernel dedicati). Vedi `constants.py:408-414` (commento esplicito).

### PERCHÉ manca la "radice SPK" (domanda dell'utente, rimasta a metà)
Nota: la web-search di conferma sul motivo JPL è stata **interrotta** dall'utente prima di concludere. La spiegazione qui sotto viene dal codice/commenti della lib e va confermata con fonte JPL/NAIF se serve certezza pubblica.

- Il servizio Horizons "genera SPK on-demand" (`EPHEM_TYPE=SPK`) per corpi minori numerati. Per Bennu questo servizio è disabilitato lato JPL perché Bennu è un **mission target OSIRIS-REx**: la sua ephemeris di alta precisione è distribuita tramite **kernel SPICE dedicati della missione** (es. `bennu_refdrmc_v1.bsp`), non tramite la generazione automatica small-body. Fonte da verificare: NAIF OSIRIS-REx SPICE archive (`https://naif.jpl.nasa.gov/pub/naif/pds/pds4/orex/`).
- Conseguenza a cascata sul `.leb` (vedi §1-bis).

### È un bug della lib?
**No.** È comportamento previsto e documentato. La lib dichiara onestamente la precisione degradata (`planets.py:2974`) e non solleva eccezione (lo strict-check `planets.py:2906-2913` esenta i corpi bloccati). Il **vero problema è nel test**: pretende `0.05°` su Bennu in una config che non ha alcuna sorgente ad alta precisione per Bennu.

Nota complementare: sul lato riferimento, pyswisseph richiede `s101955s.se1` sotto `data/reference/`. Se quella dir è vuota, pyswisseph solleverebbe e il test andrebbe in **skip** (`try/except → pytest.skip`, `test_compare_minor_bodies.py:324-326`). I diff 4-44° presuppongono `data/reference/` popolato.

### File chiave
- `validation/compare_scripts/tests/test_compare_minor_bodies.py:311-332` (test), `:78,:330` (tol)
- `validation/compare_scripts/tests/conftest.py:24-35, 54-70`
- `libephemeris/planets.py:2860-2982` (catena fallback; warning 2974)
- `libephemeris/constants.py:408-414` (`SPK_AUTO_DOWNLOAD_BLOCKED`, Bennu)
- `libephemeris/minor_bodies.py:1971-1981` (elementi Bennu, epoch 2461000.5≈2027)
- `libephemeris/rebound_integration.py:610-654` (check ASSIST)
- `libephemeris/spk_auto.py:1021-1078` (`enable_common_bodies`, Bennu escluso)

---

## 1-bis. Perché Bennu NON è nei file `.leb` (approfondimento richiesto dall'utente)

L'utente ha chiesto: "ma io ho i `.leb`, non dovrebbero evitare ASSIST/rebound?". Risposta: **sì per 8 NEA su 9, NO per Bennu**.

- I `.leb` **non sono una sorgente indipendente**: sono una **compressione Chebyshev di un SPK di riferimento**. Docstring `exotic_bodies.py:11-12`: *"Each body is served from a JPL SPK kernel over its coverage window ... and falls back to the existing Keplerian path outside it."*
- Registry canonico degli exotic in `.leb`: `libephemeris/exotic_bodies.py`, lista `_REGISTRY` (righe 46-82). La sezione NEA (righe 74-81) contiene **8** NEA: Eros, Amor, Apophis, Itokawa, Ryugu, Toro, Toutatis, Icarus. **Bennu NON c'è.**
- Motivo esplicito, `exotic_bodies.py:21-22`: *"Bennu is intentionally absent: JPL Horizons blocks SPK generation for it (see constants.SPK_AUTO_DOWNLOAD_BLOCKED)."*
- Catena: `BODY_GROUPS["exotics"] = list(EXOTIC_IDS)` (`generate_leb.py:318`); `EXOTIC_IDS` deriva da `_REGISTRY` (`exotic_bodies.py:88`) → Bennu non nel registry ⇒ non in exotics ⇒ mai fittato nel `.leb`.
- Non poteva esserci: senza SPK di riferimento non c'è oracolo da comprimere in `.leb`. È lo **stesso** buco che manda Bennu in Kepleriano.
- Esclusione aggiuntiva sul tier **extended**: gli 8 NEA vengono droppati comunque (`generate_leb.py:2428-2440`, "chaotic NEA... integration diverges over millennia"). Secondario: per Bennu il blocco è già a monte in ogni tier.

**Conclusione:** per Bennu il `.leb` non può aiutare per costruzione. Le uniche vie che non passano da SPK di riferimento sono ASSIST e Horizons REST (vedi §1-ter).

---

## 1-ter. Come "sbloccare" Bennu — 3 opzioni concrete

### Opzione A — Registrare kernel SPK OSIRIS-REx manuale
- API: `register_spk_body(ipl, spk_file, naif_id)` (`__init__.py:454`, def `spk.py:542-644`) oppure `discover_local_spks(path)` (`spk_auto.py:1667-1800`).
- Intercettato a `planets.py:2869` **prima** di ASSIST/Kepler.
- ```python
  from libephemeris import register_spk_body, BENNU
  register_spk_body(BENNU, "/percorso/bennu_refdrmc_v1.bsp", 2101955)
  ```
- **LIMITE DECISIVO:** i kernel mission coprono solo ~2018-2021 (max 2011-2023). **NON coprono 1980 né 2000** → fuori range `EphemerisRangeError` (`planets.py:2891`) ⇒ ricade in ASSIST/Kepler. Da solo NON risolve i test storici.

### Opzione B — Attivare ASSIST (rebound + assist) ⟵ risposta alla domanda "se installo assist/rebound funziona?"
- Gate: `ipl in MINOR_BODY_ELEMENTS and check_assist_data_available()` (`planets.py:2923`). Bennu È in `MINOR_BODY_ELEMENTS`.
- `check_assist_available()` (`rebound_integration.py:610`) fa solo `import assist` → oggi fallisce (pacchetti non installati).
- **Data file GIÀ presenti** in `~/.libephemeris/assist/`: `linux_p1550p2650.440` (DE440 pianeti 1550-2650), `sb441-n16.bsp` (16 perturbatori), `linux_m13000p17000.441` (DE441 extended). **Coprono 1980/2000/2024.**
- Passo: `pip install rebound assist` → il path si attiva da solo (nessun'altra config).
- Cosa fa: `propagate_orbit_assist` (`rebound_integration.py:879-1048`) integra numericamente dall'epoch (2027) al jd, con Sun/Moon/8 pianeti + 16 asteroidi massicci + J2/J3/J4 + relatività; poi `_assist_position_at` (`planets.py:2025-2114`) riduce a geocentrico apparente.
- **Risposta secca "funziona poi?":** **Sì, in linea di principio Bennu smette di cadere in Kepleriano e usa ASSIST**, con precisione attesa buona (target 0.05° plausibile per 2024/2000; al 1980, 47 anni indietro, il drift è arcsec→arcmin, borderline ma probabilmente entro 0.05°). **DA VERIFICARE EMPIRICAMENTE** (non testato: verifica interrotta).
- **LIMITI:**
  1. Yarkovsky **non applicata** su questo path: `_assist_position_at` (`planets.py:2053`) chiama `propagate_orbit_assist` senza `include_non_gravitational=True` (default False, `rebound_integration.py:884`). `YARKOVSKY_DA_DT[BENNU]` (`minor_bodies.py:1527-1528`) è usato solo nel path Kepleriano. Per Bennu su 47 anni introduce drift along-track. Tuning auspicabile: passare `include_non_gravitational=True` + derivare A1/A2/A3.
  2. Single-epoch: propaga sempre da 2027 (anche all'indietro). `multi_epoch_data.py:24` dichiara "Bennu excluded".
  3. Performance: ogni `calc_ut` re-integra da zero (no cache). Lento su molte date. `FLG_SPEED` raddoppia (finite-diff ±1s).

### Opzione C — Estendere Horizons REST (VECTORS)
- Il blocco riguarda SOLO `EPHEM_TYPE=SPK` (generazione kernel). L'endpoint REST **VECTORS** (`horizons_backend.py:182`) serve state-vector on-demand per tutti i numerati, Bennu incluso, **non è bloccato**.
- Già 5 asteroidi mappati (`horizons_backend.py:39-57`). Aggiungere: `BENNU: "101955;"` (una riga).
- Attivo quando `get_horizons_client() is not None` → mode `horizons`, o `auto` senza DE440 locale. Dispatch già cablato (`planets.py:1130-1160`).
- **PRO:** copertura totale 1980/2000/2024 + precisione massima costante (orbit solution JPL), zero file locali, 1 riga.
- **CONTRO:** richiede internet; latenza per-test (5 fetch/calc, +3 con FLG_SPEED; cache LRU 4096); no `FLG_TOPOCTR`/`FLG_NONUT` (→ fallback Skyfield).

### Tabella copertura
| Data test | A (mission SPK) | B (ASSIST) | C (Horizons REST) |
|---|---|---|---|
| 1980 | ❌ fuori range | ✅ (drift, no Yarkovsky) | ✅ max precision |
| 2000 | ❌ fuori range | ✅ (drift modesto) | ✅ max precision |
| 2024 | ⚠️ solo se in finestra | ✅ (~3 anni, ottimo) | ✅ max precision |

**Raccomandazione:** C risolve tutti e 3 con una riga ma richiede modalità horizons/auto (non skyfield). B è il miglior fallback offline (data file già presenti), tuning Yarkovsky auspicabile per il 1980. A non risolve i test storici.

---

## 2. Angles + ArabicParts — 2 fail — **BUG DEI TEST (tipo eccezione errato)**

### Sintomo
- `test_angles_require_location` (`test_fixed_stars_and_points.py:174`): `pytest.raises(ValueError, match="observer location")` — ma la lib solleva `libephemeris.exceptions.Error` (`planets.py:3035`).
- `test_arabic_parts_require_precalc` (`test_fixed_stars_and_points.py:221`): `pytest.raises(ValueError, match="pre-calculated")` — ma la lib solleva `Error` (`planets.py:3051`).

### Causa radice
`Error(Exception)` (`exceptions.py:43`) deriva da `Exception`, **NON da `ValueError`**. `issubclass(Error, ValueError) == False`. Il `match=` combacia col messaggio, ma `pytest.raises(ValueError)` non cattura `Error` ⇒ l'eccezione si propaga ⇒ FAILURE.

### La lib è coerente con pyswisseph
- pyswisseph: `swe.Error = PyErr_NewException("swisseph.Error", NULL, NULL)` (`pyswisseph.c:5896`) deriva da `Exception`, non `ValueError`. Stessa gerarchia di `libephemeris.Error`.
- `ASCENDANT` (offset 2000000) e `PARS_FORTUNAE` (offset 3000000) sono **estensioni API di libephemeris**, non esposte da pyswisseph via `calc_ut` (l'asc si ottiene con `swe.houses()`). Non c'è "comportamento pyswisseph" di riferimento → sono smoke test interni.

### Verdetto
**Bug del test.** Fix: usare `pytest.raises(ephem.Error, ...)` (o `Exception`) invece di `ValueError`.

### File chiave
- `validation/compare_scripts/tests/test_fixed_stars_and_points.py:165-175, 212-222`
- `libephemeris/planets.py:3029-3053`
- `libephemeris/exceptions.py:43-65`

---

## 3. Waldemath — 8 fail — **IL COMMIT 17c43d2 È CORRETTO; I TEST SONO VECCHI**

### Contesto
Commit recente `17c43d2 fix(hypothetical): Waldemath uses the canonical Koch geocentric orbit` ha cambiato l'implementazione. I test validation sono rimasti all'aspettativa PRE-fix.

### Confronto parametri
| Campo | Test validation (PRE-fix, ERRATI) | 17c43d2 + pyswisseph reale | Fonte |
|---|---|---|---|
| epoch | 2451545.0 (J2000) | **2414290.958 (1898)** | seorbel.txt:68 |
| a | 0.0029833 | **0.0068400705250028** | seorbel.txt:68 |
| e | 0.0 (circolare) | **0.1587** (ellittica) | seorbel.txt:68 |
| i | 0.0 | **2.5°** | seorbel.txt:68 |
| lon @J2000 | 248.88° | **33.30°** | pyswisseph misurato |
| lat | 0.0 | **−1.24°** | pyswisseph |
| dist @J2000 | 0.0029833 (costante) | **0.0059506** (variabile) | pyswisseph |
| moto giorn. | 3.025°/g | **3.94°/g** | pyswisseph |

### Verifica empirica su pyswisseph 2.10.3.2 reale a J2000
```
pyswisseph        : lon=33.298377  lat=-1.244808  dist=0.0059506  dlon=3.9379
lib (post-17c43d2): lon=33.298512  lat=-1.244803  dist=0.0059500   → delta ~0.4 arcsec
```
La lib post-fix coincide **cifra per cifra** con pyswisseph. I numeri sono **esattamente** la riga 68 di `seorbel.txt` (sorgente Swiss Ephemeris per orbite fittizie), corpo dichiarato `geo`, ellittico, inclinato, rate secolari su M/ω/Ω, equinozio 1898.

### Chiarimenti sui "sintomi"
- `0.0059506` NON è `2× 0.0029833` (rapporto 1.9946): è la distanza istantanea a J2000 vicino al periastro (a(1−e)=0.00575, a(1+e)=0.00793). Nessun raddoppio: due orbite diverse.
- Distanza variabile tra date (0.0059506 vs 0.0060935) = orbita ellittica corretta.
- `M_rate = 109023.26 °/cy` ⇒ n = 2.985°/g ⇒ periodo ~120.6 g. Il vecchio "119 g / 3.024873" derivava da 360/119, non da Swiss Eph.
- Il commento test "from seorbel.txt #18" è impreciso: Waltemath è **#19**.

### Verdetto
`17c43d2` **NON è una regressione**: è la correzione 1:1 con pyswisseph (vincolo CLAUDE.md). I test validation sono da riscrivere.

### Fix
Riscrivere `TestWaldemathCalculations` (`validation/compare_scripts/tests/test_hypothetical.py:555-663`) sull'orbita Koch. **Modello già pronto:** `tests/test_hypothetical.py::TestWaldemathGeocentricMoon` (nel repo lib) con oracle a 5 date identico a pyswisseph.

### File chiave
- `libephemeris/hypothetical.py` (costante `WALDEMATH_ELEMENTS` ~1049-1110, `calc_waldemath` ~3169-3257)
- `libephemeris/data/fictitious_orbits.csv:110`
- `validation/compare_scripts/tests/test_hypothetical.py:555-663` (da aggiornare)
- `tests/test_hypothetical.py::TestWaldemathGeocentricMoon` (modello corretto)
- `seorbel.txt:68` (in pyswisseph sdist: `~/.cache/uv/sdists-v9/pypi/pyswisseph/2.10.3.2/.../src/libswe/seorbel.txt`)

---

## 4. Timeout occultazioni planetarie — 3 fail — **PRESTAZIONI DELLA LIB (non infinite-loop)**

Unico gruppo che richiede intervento sul CODICE della libreria. NON è un bug di correttezza.

### I test (`test_deep_validation_4.py`)
Tutti e 3 avvolgono la chiamata in `try/except Exception → pytest.skip` (righe 381-388, 397-403, 418-426). Falliscono **solo perché non ritornano entro 120s** (pytest-timeout). Se `Error` fosse sollevato in fretta, farebbero skip/pass.
- `TestPlanetOccultWhenGlob::test_return_structure`: `planet_occult_when_glob(2024-01-01, VENUS, JUPITER)`. Evento reale: 2065-11-22, ~15.301 giorni dopo.
- `test_does_not_crash_with_various_planets`: 2 chiamate `(VENUS,MARS)` e `(MARS,JUPITER)` da 2020 — nessun evento in 150 anni → scan completo x2.
- `TestPlanetOccultWhenLoc::test_return_structure`: `planet_occult_when_loc(2024-01-01, VENUS, MARS, 45, 10)` → delega a glob (`eclipse.py:13275`), timeout sulla prima chiamata interna.

### Il loop (`libephemeris/eclipse.py`)
- `MAX_SEARCH_YEARS=150`, `MAX_ITERATIONS=54750` (`eclipse.py:12780`).
- Loop `12956-13038`: step iniziale 1 giorno.
- Ogni iter: `_check_occultation(jd)` (`12964`) = **2 × `calc_ut(FLG_EQUATORIAL | FLG_SPEED)`** (`12793, 12843-12846`).
- Step adattivo (`13026-13033`): 0.01g se sep<0.5, ..., 1.0g nel "deserto".
- Refinement golden-section: se sep<0.2° → `_find_minimum_separation` = 50 iter × 2 calc_ut (`12976-12977, 12865-12890`).
- Uscite: return su hit (`13023`), o `raise Error` finale (`13046`). **Loop finito, NON infinito.**

### Costo/iter misurato (2 × calc_ut FLG_EQUATORIAL|FLG_SPEED, jd=2024)
| Modalità | ms/iter |
|---|---|
| auto/LEB | **4.83** |
| skyfield | **46.4** (~10×) |

### Proiezioni
| Test | LEB/auto | Skyfield |
|---|---|---|
| Glob return_structure (~15.301 iter) | ~74s + refinement → ~95-100s (borderline) | ~710s |
| Glob does_not_crash (2× 54.750 iter) | ~529s | ~5.081s |
| Loc return_structure (54.750 iter) | ~264s | ~2.540s |

### Tre inefficienze convergenti
1. **`FLG_SPEED` richiesto ma MAI usato nello scan.** `_check_occultation` (`12837-12853`) legge solo ra/dec/dist; i contact time usano finite-diff manuale (`dt_test=1/24`, `12900`). Ogni calc_ut FLG_SPEED triplica il costo per nulla. Rimuoverlo → ~3×.
2. **Nessun pre-filtro coarse.** Anche a sep 180° paga 2 riduzioni apparenti/step giornaliero. Swiss Eph nativo fa scan coarse su ephemeris read (Chebyshev diretto, no apparent place), raffina solo vicino congiunzione.
3. **Golden-section troppo eager.** Ogni sep<0.2° (`12976`) attiva 50 iter × 2 calc_ut anche per congiunzioni non-occultanti (decine in 150 anni).

### Dettaglio modalità
`planet_occult_when_glob` **non usa** `_call_with_leb_skyfield_fallback`, chiama il `calc_ut` globale (`12793`) ed eredita la modalità globale. Con default skyfield → 46 ms/iter → timeout garantito. In LEB → 4.8 ms ma **NON basta da solo** per 2 test su 3 (restano >120s). Serve comunque l'intervento sul codice (punti 1 e 2).

### Nota pyswisseph
`pyswisseph 2.10.03` **NON implementa** `planet_occult_when_glob/loc` (solo `lun_occult_*`). Questi 3 test sono **smoke test della sola estensione libephemeris**, non confronti con pyswisseph. Nessun reference behavior da confrontare.

### File chiave
- `validation/compare_scripts/tests/test_deep_validation_4.py:370-425`
- `libephemeris/eclipse.py:12672` (glob), `12780` (limiti), `12793` (`_get_planet_position`), `12837-12853` (`_check_occultation`), `12956-13038` (loop), `13046` (raise), `13072/13275` (loc → glob delega), `99-126` (`_call_with_leb_skyfield_fallback`)

---

## 5. Riepilogo azioni consigliate (per il prossimo agente)

| Area | Fail | Natura | Azione |
|---|---|---|---|
| Bennu/NEA | 3 | Setup test + limite by-design | Decidere: (C) 1 riga Horizons + modalità, o (B) `pip install rebound assist` + eventuale tuning Yarkovsky, o rilassare/xfail il test Bennu. **VERIFICARE empiricamente** che B/C portino <0.05°. |
| Angles/ArabicParts | 2 | Bug test | Cambiare `pytest.raises(ValueError)` → `pytest.raises(ephem.Error)` in `test_fixed_stars_and_points.py:174,221`. |
| Waldemath | 8 | Test vecchi | Riscrivere `TestWaldemathCalculations` sull'orbita Koch (modello: `tests/test_hypothetical.py::TestWaldemathGeocentricMoon`). NON toccare la lib. |
| Occultazioni | 3 | Prestazioni lib | Ottimizzare `planet_occult_when_glob`: (1) togliere `FLG_SPEED` dallo scan, (2) pre-filtro coarse, (3) golden-section meno eager. Opz. secondaria: default suite → LEB. |

**Priorità suggerita:** Waldemath + Angles/ArabicParts (10 fail, fix solo-test, zero rischio) → poi decidere Bennu (scelta di prodotto) → poi ottimizzazione occultazioni (unico intervento su codice lib, più delicato).

### Questioni ancora APERTE / da verificare
1. **Bennu + ASSIST: funziona davvero <0.05°?** Non testato empiricamente (verifica interrotta). Da provare: `pip install rebound assist` nel venv validation, poi rilanciare i 3 test Bennu.
2. **Motivo JPL del blocco SPK Bennu**: confermato dai commenti del codice, ma la conferma da fonte pubblica JPL/NAIF è rimasta a metà (web-search interrotta). Se serve certezza, verificare l'archivio NAIF OSIRIS-REx.
3. **Default validation → LEB?** Ortogonale ai 16 fail; migliora solo (parzialmente) le occultazioni. Decisione separata: cambiare `conftest.py:24` default o far esportare `LIBEPHEMERIS_COMPARE_MODE=leb` da `run.sh`.

---

*Documento generato durante indagine (no fix applicati). Tutti i file:riga verificati alla data indagine sul branch `fix/v3.0.0-review-round2`.*
