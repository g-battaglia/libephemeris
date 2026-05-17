# Capitolo 14 — Precisione, performance e modalità LEB

## Cosa imparerai

In questo capitolo scoprirai quanto è precisa LibEphemeris rispetto alle fonti astronomiche di riferimento, come si comporta la precisione su scale temporali diverse (dall'antichità al futuro lontano), cos'è la modalità LEB (efemeridi binarie precompilate) e come ottimizzare la velocità di calcolo, come configurare le modalità di calcolo e i livelli di precisione, e come gestire il download e la cache dei dati.

---

## 14.1 Quanto è precisa questa libreria?

Una domanda legittima prima di affidare i propri calcoli a qualsiasi software: quanto ci possiamo fidare dei risultati?

La risposta breve è: **per l'uso moderno (1900–2100), la precisione è sub-arcsecondo per tutti i corpi**. Questo significa che l'errore è inferiore a 1/3600 di grado — invisibile a qualsiasi scopo pratico, astrologico o astronomico amatoriale.

Vediamo il dettaglio per categoria.

### Pianeti principali

LibEphemeris usa direttamente le efemeridi **JPL DE440** (o DE441 per date estese), le stesse prodotte dal Jet Propulsion Laboratory della NASA per le missioni spaziali. Non stiamo "approssimando" le posizioni dei pianeti: stiamo leggendo gli stessi dati che la NASA usa per navigare le sonde.

La precisione è **sub-milliarcsecondo** per Sole, Luna e pianeti — il massimo ottenibile con la tecnologia attuale:

```python
import libephemeris as ephem

# Posizione del Sole al J2000.0 (1 gennaio 2000, mezzogiorno)
jd = ephem.julday(2000, 1, 1, 12.0)
pos, _ = ephem.calc_ut(jd, ephem.SUN, ephem.FLG_SPEED)
print(f"Sole al J2000.0:")
print(f"  Longitudine: {pos[0]:.9f}°")
print(f"  Latitudine:  {pos[1]:.9f}°")
print(f"  Distanza:    {pos[2]:.9f} UA")
print(f"  Velocità:    {pos[3]:.6f}°/giorno")
```

```
Sole al J2000.0:
  Longitudine: 280.368917535°
  Latitudine:  0.000227101°
  Distanza:    0.983328100 UA
  Velocità:    1.019428°/giorno
```

### Luna

La Luna è il corpo celeste più complesso da calcolare. Le sue perturbazioni gravitazionali (dovute al Sole, a Giove, alla forma non sferica della Terra) richiedono migliaia di termini. La precisione di LibEphemeris per la Luna è **sub-arcsecondo**, identica a quella delle efemeridi JPL.

### Stelle fisse

Le posizioni delle stelle fisse provengono dal catalogo **Hipparcos** (ESA, 1997), con una precisione originale di circa **1 milliarcsecondo**. La propagazione del moto proprio introduce un errore che cresce nel tempo, ma per pochi secoli dal catalogo (epoca J2000.0) la precisione rimane sotto **0.5"**.

### Case astrologiche

I calcoli delle cuspidi sono puramente geometrici (dipendono da ARMC, latitudine e obliquità dell'eclittica). La precisione è circa **0.02"** per tutti i sistemi di case — praticamente identica a qualsiasi altra implementazione.

### Ayanamsha

L'ayanamsha (precessione per lo zodiaco siderale) ha una precisione di circa **0.07"** rispetto ad altre implementazioni. Le differenze tra software diversi sono dovute a scelte di modello (quale formula di precessione usare, quale epoca di riferimento), non a errori di calcolo.

### Dove siamo più precisi

LibEphemeris usa fonti più recenti di molti altri software per due calcoli specifici:

- **Magnitudini planetarie**: implementazione basata su Mallama & Hilton (2018), più precisa dei modelli anni '80 usati altrove
- **Diametri apparenti**: valori IAU 2015, aggiornati rispetto alle costanti storiche

### Delta T: dove la precisione cala

Il **Delta T** (differenza TT − UT) è il tallone d'Achille di qualsiasi efemeride storica. Per date lontane dal presente, il Delta T è stimato con formule empiriche la cui incertezza cresce rapidamente:

```python
import libephemeris as ephem

# Delta T a diverse epoche
for anno in [1600, 1800, 1900, 1950, 2000, 2024]:
    jd = ephem.julday(anno, 1, 1, 12.0)
    dt = ephem.deltat(jd)
    print(f"Anno {anno}: Delta T = {dt * 86400:7.2f} secondi")
```

```
Anno 1600: Delta T =  109.11 secondi
Anno 1800: Delta T =   18.37 secondi
Anno 1900: Delta T =   -1.97 secondi
Anno 1950: Delta T =   28.93 secondi
Anno 2000: Delta T =   63.83 secondi
Anno 2024: Delta T =   69.18 secondi
```

Per l'anno 2000, il Delta T è noto con precisione sub-millisecondo (grazie ai dati IERS). Per l'anno 1600, l'incertezza è di diversi secondi. Per l'anno -500 (astronomia antica), l'incertezza può essere di **minuti** — il che si traduce in errori significativi sulla posizione della Luna.

---

## 14.2 Le scale temporali di errore

Non tutti i calcoli hanno la stessa validità temporale. Ogni componente della libreria ha un "orizzonte di precisione" oltre il quale l'errore cresce.

### Efemeridi JPL

Le efemeridi JPL sono generate integrando le equazioni del moto del Sistema Solare. Ogni versione copre un intervallo specifico:

- **DE440** (tier `medium`): **1550–2650** — precisione piena per oltre 1000 anni
- **DE440s** (tier `base`): **1849–2150** — versione leggera, sufficiente per uso moderno
- **DE441** (tier `extended`): **-13200 a +17191** — per ricerca storica e preistorica

Fuori dall'intervallo dell'efemeride scelta, la libreria solleva un errore: non produce risultati falsi.

### Delta T

Prima del 1600 circa, l'incertezza su Delta T diventa il fattore dominante dell'errore. Per la Luna (che si muove di ~0.5"/secondo), un errore di 30 secondi su Delta T produce un errore di ~15" sulla posizione — ancora accettabile per l'astrologia, ma non per eclissi storiche precise.

### Stelle fisse

Il moto proprio delle stelle è misurato con precisione dal catalogo Hipparcos (epoca J2000.0). La propagazione è lineare e precisa per **pochi secoli** dal catalogo. Per date molto antiche o future, il moto proprio reale potrebbe essere non lineare (a causa di compagni binari, perturbazioni gravitazionali).

### Corpi minori (Keplerian)

Per asteroidi e comete calcolati con propagazione kepleriana (senza file SPK), la precisione si degrada **rapidamente** — arcminuti in pochi mesi. Per questo la libreria supporta il download automatico di file SPK ad alta precisione da JPL Horizons.

### Regola pratica

Per **astrologia moderna** (1900–2100), la precisione è sub-arcsecondo per tutti i corpi principali, case e ayanamsha. Non c'è ragione di preoccuparsi della precisione per questo intervallo.

---

## 14.3 I livelli di precisione (Precision Tiers)

LibEphemeris organizza le efemeridi in tre **livelli di precisione** (tier), ciascuno con un diverso compromesso tra copertura temporale e dimensione dei file:

- **`base`** — usa `de440s.bsp` (~31 MB), copre **1849–2150**. Ideale per applicazioni moderne che non richiedono date storiche.

- **`medium`** (predefinito) — usa `de440.bsp` (~114 MB), copre **1550–2650**. Il miglior compromesso per la maggior parte degli usi.

- **`extended`** — usa `de441.bsp` (~3.1 GB), copre **-13200 a +17191**. Per ricerca storica, archeologia astronomica, calcoli su millenni.

### Selezionare un tier

```python
import libephemeris as ephem

# Vedere il tier corrente
print(f"Tier attivo: {ephem.get_precision_tier()}")

# Cambiare tier
ephem.set_precision_tier("base")
print(f"Nuovo tier: {ephem.get_precision_tier()}")

# Tornare al default
ephem.set_precision_tier("medium")
print(f"Tier ripristinato: {ephem.get_precision_tier()}")
```

```
Tier attivo: medium
Nuovo tier: base
Tier ripristinato: medium
```

Il tier può anche essere impostato tramite la variabile d'ambiente `LIBEPHEMERIS_PRECISION`:

```bash
export LIBEPHEMERIS_PRECISION=extended
```

### Informazioni sui tier disponibili

```python
import libephemeris as ephem
from libephemeris.state import list_tiers

for tier in list_tiers():
    print(f"{tier.name:10s}  {tier.ephemeris_file:14s}  {tier.description}")
```

```
base        de440s.bsp      Modern usage (1850-2150), ~31 MB
medium      de440.bsp       General purpose (1550-2650), ~114 MB
extended    de441.bsp       Extended range (-13200 to +17191), ~3.1 GB
```

---

## 14.4 Modalità LEB: efemeridi binarie precompilate

### Il problema

Skyfield (il motore di calcolo sottostante) è estremamente preciso, ma ogni singola chiamata a `calc_ut()` comporta diversi passaggi: lettura del file SPK, interpolazione dei polinomi di Chebyshev, rotazione dei frame di riferimento, correzione per aberrazione, nutazione, ecc. Per un singolo calcolo è impercettibile, ma per migliaia (un'efemeride mensile, una ricerca di transiti su un anno intero) il costo si accumula.

### La soluzione: LEB

**LEB** (LibEphemeris Binary) è un formato di efemeridi precompilate che memorizza approssimazioni con **polinomi di Chebyshev** per ogni corpo celeste, ogni intervallo temporale. Il file `.leb` contiene coefficienti ottimizzati che permettono di calcolare le posizioni con una sola valutazione polinomiale, saltando tutto il pipeline di Skyfield.

La precisione LEB è **identica** a quella di Skyfield — le differenze sono sotto il milliarcsecondo:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Calcolo con Skyfield
ephem.set_calc_mode("skyfield")
pos_sky, _ = ephem.calc_ut(jd, ephem.MARS, ephem.FLG_SPEED)

# Calcolo con LEB
ephem.set_calc_mode("leb")
pos_leb, _ = ephem.calc_ut(jd, ephem.MARS, ephem.FLG_SPEED)

diff_arcsec = abs(pos_sky[0] - pos_leb[0]) * 3600
print(f"Marte (Skyfield): {pos_sky[0]:.9f}°")
print(f"Marte (LEB):      {pos_leb[0]:.9f}°")
print(f"Differenza:       {diff_arcsec:.9f}\"")

ephem.set_calc_mode("auto")
```

```
Marte (Skyfield): 342.845795946°
Marte (LEB):      342.845795946°
Differenza:       0.000000000"
```

### Attivare LEB

Ci sono tre modi per attivare la modalità LEB:

**1. Scaricare il file LEB precompilato (consigliato):**

```python
import libephemeris as ephem

# Scarica il file LEB per il tier medium (~175 MB)
ephem.download_leb_for_tier("medium")
# Il file viene salvato in ~/.libephemeris/leb/ephemeris_medium.leb
# e attivato automaticamente per questa sessione
```

**2. Impostare il percorso manualmente:**

```python
import libephemeris as ephem

ephem.set_leb_file("/percorso/al/file/ephemeris_medium.leb")
```

**3. Tramite variabile d'ambiente:**

```bash
export LIBEPHEMERIS_LEB=/percorso/al/file/ephemeris_medium.leb
```

### Auto-discovery

Se hai scaricato un file LEB con `download_leb_for_tier()`, la libreria lo trova automaticamente nella posizione standard (`~/.libephemeris/leb/ephemeris_{tier}.leb`) senza bisogno di configurazione.

### Fallback automatico

Non tutti i corpi e le combinazioni di flag sono supportati da LEB. In particolare, LEB **non** supporta:

- `FLG_TOPOCTR` (posizioni topocentriche)
- `FLG_XYZ` (coordinate cartesiane)
- `FLG_RADIANS` (angoli in radianti)
- `FLG_NONUT` (senza nutazione)

Quando incontri uno di questi casi, la libreria passa automaticamente a Skyfield senza errori — il passaggio è trasparente.

### Tre tier LEB

Come le efemeridi JPL, anche i file LEB esistono in tre versioni:

- **`base`** (~53 MB) — copre 1850–2150
- **`medium`** (~175 MB) — copre 1550–2650
- **`extended`** (~1.6 GB) — copre -5000 a +5000

---

## 14.5 Modalità di calcolo e configurazione

La libreria supporta quattro modalità di calcolo, controllabili con `set_calc_mode()`:

### `"auto"` (predefinita)

Prova prima LEB se un file `.leb` è configurato (o auto-scoperto), poi l'API Horizons se non è disponibile un file DE440 locale, infine Skyfield. È la modalità raccomandata per l'uso normale:

```python
import libephemeris as ephem

ephem.set_calc_mode("auto")
print(f"Modalità: {ephem.get_calc_mode()}")
```

```
Modalità: auto
```

### `"skyfield"`

Forza sempre il percorso Skyfield, anche se un file LEB è disponibile. Utile per debug, confronti di precisione o quando vuoi essere sicuro di usare il pipeline completo:

```python
import libephemeris as ephem

ephem.set_calc_mode("skyfield")
print(f"Modalità: {ephem.get_calc_mode()}")

jd = ephem.julday(2024, 4, 8, 12.0)
pos, _ = ephem.calc_ut(jd, ephem.SUN, ephem.FLG_SPEED)
print(f"Sole (Skyfield forzato): {pos[0]:.6f}°")

ephem.set_calc_mode("auto")  # ripristina
```

```
Modalità: skyfield
Sole (Skyfield forzato): 19.140437°
```

### `"leb"`

Richiede un file LEB valido. La libreria prova il file configurato, poi l'auto-discovery, poi l'auto-download di LEB2 per il tier attivo. Solleva `RuntimeError` solo se nessun LEB può essere risolto. I corpi non presenti nel file LEB passano comunque a Skyfield:

```python
import libephemeris as ephem

ephem.set_calc_mode("leb")
print(f"Modalità: {ephem.get_calc_mode()}")

ephem.set_calc_mode("auto")  # ripristina
```

```
Modalità: leb
```

### `"horizons"`

Usa sempre l'API REST NASA JPL Horizons per i calcoli. Questa modalità richiede una connessione internet e non necessita di file di efemeridi locali. Supporta pianeti, asteroidi, Nodo Medio, Apogeo Medio e Uraniani. Corpi o flag non supportati da Horizons (es. `FLG_TOPOCTR`, stelle fisse) passano automaticamente a Skyfield:

```python
import libephemeris as ephem

ephem.set_calc_mode("horizons")
print(f"Modalità: {ephem.get_calc_mode()}")

jd = ephem.julday(2024, 4, 8, 12.0)
pos, _ = ephem.calc_ut(jd, ephem.SUN, ephem.FLG_SPEED)
print(f"Sole (Horizons): {pos[0]:.6f}°")

ephem.set_calc_mode("auto")  # ripristina
```

```
Modalità: horizons
Sole (Horizons): 19.140437°
```

### Variabile d'ambiente

La modalità può essere impostata anche tramite `LIBEPHEMERIS_MODE`:

```bash
export LIBEPHEMERIS_MODE=horizons
```

---

## 14.6 EphemerisContext: calcoli thread-safe

Per applicazioni multi-thread (server web, API, calcoli paralleli), la libreria offre `EphemerisContext` — un contesto che mantiene il proprio stato isolato (posizione dell'osservatore, modo siderale, cache degli angoli) condividendo le risorse costose (file delle efemeridi, timescale) in modo thread-safe.

### Uso base

```python
import libephemeris as ephem
from libephemeris import EphemerisContext
from libephemeris.constants import SUN, FLG_SPEED

ctx = EphemerisContext()
ctx.set_topo(12.5, 41.9, 0)  # Roma

jd = ephem.julday(2024, 4, 8, 12.0)
pos, flag = ctx.calc_ut(jd, SUN, FLG_SPEED)
print(f"Sole (da Roma, via contesto): {pos[0]:.4f}°")
```

```
Sole (da Roma, via contesto): 19.1404°
```

### Calcolo siderale nel contesto

Ogni contesto può avere il proprio modo siderale:

```python
import libephemeris as ephem
from libephemeris import EphemerisContext
from libephemeris.constants import SUN, FLG_SPEED, FLG_SIDEREAL
from libephemeris.constants import SE_SIDM_LAHIRI

ctx = EphemerisContext()
ctx.set_sid_mode(SE_SIDM_LAHIRI)

jd = ephem.julday(2024, 4, 8, 12.0)
pos, _ = ctx.calc_ut(jd, SUN, FLG_SPEED | FLG_SIDEREAL)
print(f"Sole siderale (Lahiri): {pos[0]:.4f}°")
```

```
Sole siderale (Lahiri): 354.9458°
```

### Case nel contesto

```python
import libephemeris as ephem
from libephemeris import EphemerisContext

ctx = EphemerisContext()
jd = ephem.julday(2024, 4, 8, 12.0)

cusps, ascmc = ctx.houses(jd, 41.9, 12.5, ord('P'))
print(f"Ascendente: {ascmc[0]:.4f}°")
print(f"Medio Cielo: {ascmc[1]:.4f}°")
```

```
Ascendente: 133.0806°
Medio Cielo: 31.9078°
```

### Calcoli paralleli con thread

```python
import libephemeris as ephem
from libephemeris import EphemerisContext
from libephemeris.constants import SUN, FLG_SPEED
import threading

risultati = {}

def calcola_per_citta(nome, lon, lat):
    ctx = EphemerisContext()
    ctx.set_topo(lon, lat, 0)
    jd = ephem.julday(2024, 4, 8, 12.0)
    pos, _ = ctx.calc_ut(jd, SUN, FLG_SPEED)
    risultati[nome] = pos[0]

citta = [
    ("Roma", 12.5, 41.9),
    ("New York", -74.0, 40.7),
    ("Tokyo", 139.7, 35.7),
]

threads = []
for nome, lon, lat in citta:
    t = threading.Thread(target=calcola_per_citta, args=(nome, lon, lat))
    threads.append(t)
    t.start()

for t in threads:
    t.join()

for nome, lon in risultati.items():
    print(f"{nome}: Sole a {lon:.4f}°")
```

```
Roma: Sole a 19.1404°
New York: Sole a 19.1404°
Tokyo: Sole a 19.1404°
```

> **Nota**: le posizioni geocentriche sono uguali per tutte le città perché il flag `FLG_TOPOCTR` non è stato usato. La differenza si vede nelle case e nelle posizioni topocentriche.

### LEB nel contesto

Ogni contesto può avere il proprio file LEB:

```python
from libephemeris import EphemerisContext

ctx = EphemerisContext()
ctx.set_leb_file("/percorso/al/file/ephemeris_medium.leb")
# I calcoli su questo contesto useranno LEB
```

### Chiusura delle risorse condivise

```python
from libephemeris import EphemerisContext

# Chiude file e risorse condivise da tutti i contesti
EphemerisContext.close()
```

---

## 14.7 Download e gestione dati

LibEphemeris scarica automaticamente i file necessari alla prima esecuzione, ma offre anche funzioni per gestire i download in modo esplicito.

### Scaricare dati per un tier

```python
import libephemeris as ephem

# Scarica tutto il necessario per il tier "medium":
# - de440.bsp (efemeridi JPL)
# - planet_centers_medium.bsp (centri dei pianeti)
# - SPK per tutti i corpi minori
ephem.download_for_tier("medium")
```

### Scaricare il file LEB

```python
import libephemeris as ephem

# Scarica le efemeridi precompilate LEB per il tier "medium"
# (~175 MB, attivato automaticamente dopo il download)
ephem.download_leb_for_tier("medium")
```

### Directory dei dati

Per impostazione predefinita, tutti i file vengono salvati in `~/.libephemeris/`. Puoi cambiare questa directory con la variabile d'ambiente:

```bash
export LIBEPHEMERIS_DATA_DIR=/dati/efemeridi
```

Oppure verificare la directory corrente:

```python
import libephemeris as ephem

print(f"Directory dati: {ephem.get_library_path()}")
```

```
Directory dati: ~/.libephemeris
```

### Reset completo

La funzione `close()` chiude tutti i file e resetta lo stato globale. Utile per liberare risorse o ricominciare con una configurazione pulita:

```python
import libephemeris as ephem

# Fai alcuni calcoli...
jd = ephem.julday(2024, 4, 8, 12.0)
pos, _ = ephem.calc_ut(jd, ephem.SUN, 0)
print(f"Prima di close: {pos[0]:.6f}°")

# Chiudi tutto
ephem.close()

# Il prossimo calcolo ricarica automaticamente le efemeridi
pos2, _ = ephem.calc_ut(jd, ephem.SUN, 0)
print(f"Dopo close: {pos2[0]:.6f}°")
print(f"Risultato identico: {abs(pos[0] - pos2[0]) < 0.000001}")
```

```
Prima di close: 19.140437°
Dopo close: 19.140437°
Risultato identico: True
```

### Riepilogo variabili d'ambiente

- `LIBEPHEMERIS_DATA_DIR` — directory base per tutti i dati (default: `~/.libephemeris`)
- `LIBEPHEMERIS_PRECISION` — tier di precisione: `base`, `medium`, `extended`
- `LIBEPHEMERIS_EPHEMERIS` — nome del file efemeride (es. `de441.bsp`)
- `LIBEPHEMERIS_LEB` — percorso al file `.leb` per la modalità binaria
- `LIBEPHEMERIS_MODE` — modalità di calcolo: `auto`, `skyfield`, `leb`, `horizons`
- `LIBEPHEMERIS_AUTO_SPK` — `1`/`0` per abilitare/disabilitare il download automatico SPK
- `LIBEPHEMERIS_SPK_DIR` — directory cache per file SPK dei corpi minori

---

## Riepilogo

- **Precisione sub-arcsecondo** per tutti i corpi nell'intervallo 1900–2100 — nessun compromesso per l'uso moderno
- **Tre tier di precisione**: `base` (1849–2150, 31 MB), `medium` (1550–2650, 114 MB), `extended` (-13200 a +17191, 3.1 GB)
- **LEB** (efemeridi binarie precompilate): velocizzano i calcoli mantenendo precisione identica a Skyfield; attivazione con `set_leb_file()` o `download_leb_for_tier()`
- **Quattro modalità di calcolo**: `auto` (predefinita — LEB, poi Horizons, poi Skyfield), `skyfield` (forza Skyfield), `leb` (richiede LEB), `horizons` (API NASA JPL Horizons, richiede connessione internet)
- **EphemerisContext**: contesto thread-safe per calcoli paralleli, con stato isolato per osservatore, modo siderale e cache
- **Download automatico** dei dati al primo utilizzo; gestione esplicita con `download_for_tier()` e `download_leb_for_tier()`
- **`close()`** per reset completo delle risorse e dello stato

## Funzioni introdotte

- `set_precision_tier(tier)` / `get_precision_tier()` — seleziona il livello di precisione (`"base"`, `"medium"`, `"extended"`)
- `set_calc_mode(mode)` / `get_calc_mode()` — imposta la modalità di calcolo (`"auto"`, `"skyfield"`, `"leb"`, `"horizons"`)
- `set_leb_file(filepath)` — attiva le efemeridi binarie precompilate
- `download_for_tier(tier)` — scarica tutti i dati per un tier
- `download_leb_for_tier(tier)` — scarica il file LEB precompilato
- `EphemerisContext()` — contesto isolato per calcoli thread-safe
- `EphemerisContext.calc_ut(jd, body, flag)` — calcolo posizione nel contesto
- `EphemerisContext.houses(jd, lat, lon, hsys)` — case nel contesto
- `EphemerisContext.set_topo(lon, lat, alt)` — osservatore nel contesto
- `EphemerisContext.set_sid_mode(mode)` — modo siderale nel contesto
- `EphemerisContext.set_leb_file(filepath)` — LEB per contesto
- `EphemerisContext.close()` — chiude risorse condivise
- `close()` — chiude tutte le risorse e resetta lo stato globale
- `get_library_path()` — percorso della directory dati
