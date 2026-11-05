# Capitolo 14 — Precisione, performance e modalità LEB

[Indice e configurazione](README.md) · [Precedente](13-pianeti-ipotetici-parti-arabe.md) · [Successivo](15-ricettario.md) · [English](../en/14-precision-performance.md)

## Cosa imparerai

In questo capitolo scoprirai quanto è precisa LibEphemeris rispetto alle fonti astronomiche di riferimento, come si comporta la precisione su scale temporali diverse (dall'antichità al futuro lontano), cos'è la modalità LEB (efemeridi binarie precompilate) e come ottimizzare la velocità di calcolo, come configurare le modalità di calcolo e i livelli di precisione, e come gestire il download e la cache dei dati.

---

## 14.1 Quanto è precisa questa libreria?

Una domanda legittima prima di affidare i propri calcoli a qualsiasi software: quanto ci possiamo fidare dei risultati?

La risposta utile parte da **quale corpo, data, dato e grandezza osservabile?** Un secondo d'arco è 1/3600 di grado, ma un errore numerico inferiore non garantisce un'accuratezza fisica universale. Incertezza delle sorgenti, rotazione terrestre, posizione dell'osservatore e atmosfera possono contare di più. Una longitudine zodiacale, un contatto di occultazione e la visibilità a occhio nudo hanno bilanci d'errore diversi.

Vediamo il dettaglio per categoria.

### Pianeti principali

LibEphemeris legge **JPL DE440/DE441** oppure coefficienti LEB locali. Le efemeridi JPL sono soluzioni orbitali numeriche adattate alle osservazioni, non traiettorie esatte. LEB aggiunge scelte di interpolazione e memorizzazione; le posizioni apparenti richiedono anche tempo-luce e riduzioni di coordinate. Consulta la pubblicazione e la copertura effettiva del file prima di attribuire un limite d'errore.

L'esempio stampa diverse cifre per illustrare il risultato, non per certificarle tutte. J2000.0 è mezzogiorno **TT**: usiamo `calc`, senza passare la stessa etichetta a `calc_ut`.

```python
import libephemeris as ephem

# J2000.0 è il 1 gennaio 2000 a mezzogiorno TT, non UT1.
jd_tt = ephem.julday(2000, 1, 1, 12.0)
pos, _ = ephem.calc(jd_tt, ephem.SUN, ephem.FLG_SPEED)
print(f"Sole al J2000.0:")
print(f"  Longitudine: {pos[0]:.9f}°")
print(f"  Latitudine:  {pos[1]:.9f}°")
print(f"  Distanza:    {pos[2]:.9f} UA")
print(f"  Velocità:    {pos[3]:.6f}°/giorno")
```

```
Sole al J2000.0:
  Longitudine: 280.368165285°
  Latitudine:  0.000227407°
  Distanza:    0.983327632 UA
  Velocità:    1.019434°/giorno
```

### Luna

La Luna è il corpo celeste più complesso da calcolare. Le sue perturbazioni gravitazionali (dovute al Sole, a Giove, alla forma non sferica della Terra) richiedono migliaia di termini. La traiettoria lunare scelta e le riduzioni apparenti contribuiscono entrambe al risultato. La vicinanza della Luna rende importanti posizione dell'osservatore e orientamento terrestre: un confronto geocentrico da solo non valida un contatto di eclissi locale.

### Stelle fisse

Il catalogo compilato combina **Hipparcos**, XHIP e fonti di identificazione. Le osservazioni Hipparcos hanno epoca J1991.25; le coordinate memorizzate a J2000 sono una rappresentazione propagata, non nuove osservazioni. Ogni stella ha la propria incertezza. L'incertezza del moto proprio cresce lontano dall'epoca di riferimento e un moto binario può rendere insufficiente la propagazione semplice. Non assegnare un unico limite a tutte le stelle e a tutti i secoli.

### Case astrologiche

I calcoli delle cuspidi sono puramente geometrici (dipendono da ARMC, latitudine e obliquità dell'eclittica). La riproducibilità richiede la stessa definizione di casa, UT1, latitudine e obliquità. Vicino alle singolarità polari, piccoli cambiamenti possono produrre grandi variazioni o rendere la costruzione impossibile. Una differenza tra sistemi non è di per sé un errore numerico.

### Ayanamsha

L'ayanamsha è una convenzione siderale scelta, non un angolo universale misurato direttamente. Confronta epoca di definizione, modello di precessione e ancoraggio stellare prima di interpretare una discrepanza. L'accordo con un altro programma misura la coerenza sotto quelle scelte, non dimostra che esista un'unica origine zodiacale fisicamente corretta.

### Identificare il modello fisico

Due esempi utili per capire perché conta la scelta delle fonti:

- **Magnitudini planetarie**: modelli basati su Mallama & Hilton (2018), con validità dipendente dal corpo e dalla fase
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
Anno 1800: Delta T =   12.60 secondi
Anno 1900: Delta T =   -2.70 secondi
Anno 1950: Delta T =   29.15 secondi
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

Per una carta moderna, parti dai dati base, dalla conversione corretta UT1/TT e da una convenzione esplicita per case e zodiaco siderale. Per lavori storici o scientifici, dichiara grandezza richiesta e incertezza tollerata. Confronta input allineati e osservazioni indipendenti quando disponibili: un intervallo familiare non elimina la necessità di verificare il calcolo concreto.

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
Tier attivo: base
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

Il percorso nativo Rust JPL/SPK compie diversi passaggi per ogni richiesta `calc_ut()`: lettura dei segmenti, interpolazione dei polinomi di Chebyshev, rotazione dei frame e correzioni della posizione apparente. Per un singolo calcolo è impercettibile, ma per migliaia (un'efemeride mensile, una ricerca di transiti su un anno intero) il costo si accumula.

### La soluzione: LEB

**LEB** (LibEphemeris Binary) è un formato di efemeridi precompilate che memorizza approssimazioni con **polinomi di Chebyshev** per ogni corpo celeste, ogni intervallo temporale. Un file `.leb` o `.leb2` memorizza canali polinomiali e dati ausiliari per un accesso locale efficiente. I canali grezzi richiedono ancora le riduzioni apparenti appropriate: LEB non elimina tutta la fisica dopo una sola valutazione. I blocchi LEB2 possono ridurre l'I/O; la lettura indicizzata limita lo spazio per i coefficienti decodificati.

Confronta LEB e JPL nativo alla stessa epoca, con gli stessi frame e impostazioni. Questo esempio richiede un kernel JPL locale oltre ai dati LEB. Misura una differenza di longitudine, non l'errore massimo su tutti i corpi e tutte le date:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# JPL/SPK nativo Rust (nome storico della modalità)
ephem.set_calc_mode("skyfield")
pos_sky, _ = ephem.calc_ut(jd, ephem.MARS, ephem.FLG_SPEED)

# Calcolo con LEB
ephem.set_calc_mode("leb")
pos_leb, _ = ephem.calc_ut(jd, ephem.MARS, ephem.FLG_SPEED)

diff_arcsec = abs(pos_sky[0] - pos_leb[0]) * 3600
print(f"Marte (JPL nativo): {pos_sky[0]:.9f}°")
print(f"Marte (LEB):      {pos_leb[0]:.9f}°")
print(f"Differenza:       {diff_arcsec:.9f}\"")

ephem.set_calc_mode("auto")
```

```
Marte (JPL nativo): 342.845795941°
Marte (LEB):      342.845795941°
Differenza:       0.000001281"
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

import os
from pathlib import Path

# Prima di avviare Python, imposta LIBEPHEMERIS_DATA_DIR sulla radice preparata.
path = Path(os.environ["LIBEPHEMERIS_DATA_DIR"]) / "leb" / "base_core.leb2"
assert path.is_file(), "Prepara prima il gruppo LEB2 base"
ephem.set_leb_file(path)
position, _ = ephem.calc_ut(2460409.0, ephem.SUN, ephem.FLG_SPEED)
assert len(position) == 6
```

**3. Tramite variabile d'ambiente:**

```bash
export LIBEPHEMERIS_LEB=/percorso/al/file/ephemeris_medium.leb
```

### Auto-discovery

Se hai scaricato un file LEB con `download_leb_for_tier()`, la libreria lo trova automaticamente nella posizione standard (`~/.libephemeris/leb/ephemeris_{tier}.leb`) senza bisogno di configurazione.

### Fallback automatico

LEB nativo supporta risultati topocentrici, cartesiani, in radianti e senza nutazione quando il corpo/modello fornisce gli input necessari. Questi flag non implicano più un passaggio a Skyfield Python.

La modalità forzata `leb` resta su coefficienti e modelli locali supportati. Corpi mancanti, file espliciti corrotti ed epoche non coperte sono errori, non un'autorizzazione a usare silenziosamente JPL o Horizons. In `auto`, l'instradamento è una politica separata. Un servizio può preferire un errore esplicito a un cambio di dati o a una richiesta di rete inattesi. Controlla la copertura e prepara i gruppi mancanti prima del calcolo.

### Tre tier LEB

Come le efemeridi JPL, anche i file LEB esistono in tre versioni:

- **`base`** — copertura moderna, nominalmente 1850–2150 per il nucleo.
- **`medium`** — intervallo storico più ampio, nominalmente 1550–2650.
- **`extended`** — prodotti a lungo termine; controlla i limiti dei corpi e dei dati ausiliari senza presumere l'intera estensione del kernel JPL.

LEB1 e LEB2 a blocchi hanno dimensioni diverse; i gruppi opzionali asteroidi,
apogeo ed esotici sono prodotti separati. `libephemeris leb info FILE --json`
mostra canali ed epoche del file effettivo. Un tier più grande non garantisce
che ogni corpo aggiuntivo sia presente o intrinsecamente più accurato.

---

## 14.5 Modalità di calcolo e configurazione

La libreria supporta quattro modalità di calcolo, controllabili con `set_calc_mode()`:

### `"auto"` (predefinita)

Usa risorse configurate e instradamento nativo. Disponibilità, supporto del corpo e politica di rete governano il risultato: non presumere una sequenza universale LEB → Horizons → Skyfield. Per un esercizio offline riproducibile, `leb` esplicito con `sealed` è più semplice da controllare:

```python
import libephemeris as ephem

ephem.set_calc_mode("auto")
print(f"Modalità: {ephem.get_calc_mode()}")
```

```
Modalità: auto
```

### `"skyfield"`

Seleziona JPL/SPK nativo Rust anche quando LEB è disponibile. La stringa è un nome di compatibilità, non una dipendenza Python. Prepara prima un kernel JPL locale adeguato; mantieni la rete chiusa per i confronti locali:

```python
import libephemeris as ephem

ephem.set_calc_mode("skyfield")
print(f"Modalità: {ephem.get_calc_mode()}")

jd = ephem.julday(2024, 4, 8, 12.0)
pos, _ = ephem.calc_ut(jd, ephem.SUN, ephem.FLG_SPEED)
print(f"Sole (JPL nativo): {pos[0]:.6f}°")

ephem.set_calc_mode("auto")  # ripristina
```

```
Modalità: skyfield
Sole (JPL nativo): 19.140437°
```

### `"leb"`

Seleziona coefficienti LEB e modelli locali supportati. Prepara esplicitamente il file necessario; il base incluso basta per i primi esercizi. File espliciti mancanti/corrotti e richieste non coperte producono errori. Non c'è sostituzione implicita con JPL/Horizons e il calcolo sealed non può scaricare gli input mancanti:

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

Seleziona le operazioni di rete JPL Horizons supportate. È un esempio online facoltativo, separato dagli esercizi offline. Verifica le combinazioni nel contratto attuale della funzione, senza presumere un equivalente remoto per ogni corpo e flag. Imposta deliberatamente la politica di rete e ripristinala dopo:

```python
import libephemeris as ephem

ephem.set_network_policy("allow")  # esempio online facoltativo
ephem.set_calc_mode("horizons")
print(f"Modalità: {ephem.get_calc_mode()}")

jd = ephem.julday(2024, 4, 8, 12.0)
pos, _ = ephem.calc_ut(jd, ephem.SUN, ephem.FLG_SPEED)
print(f"Sole (Horizons): {pos[0]:.6f}°")

ephem.set_calc_mode("leb")
ephem.set_network_policy("sealed")  # ripristina gli esercizi offline
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
from libephemeris.constants import SIDM_LAHIRI

ctx = EphemerisContext()
ctx.set_sid_mode(SIDM_LAHIRI)

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

# Ogni worker deve produrre un risultato; l'ordine di completamento può variare.
assert len(risultati) == len(citta)
for nome, lon in sorted(risultati.items()):
    print(f"{nome}: Sole a {lon:.4f}°")
```

```
New York: Sole a 19.1404°
Roma: Sole a 19.1404°
Tokyo: Sole a 19.1404°
```

> **Nota**: le posizioni geocentriche sono uguali per tutte le città perché il flag `FLG_TOPOCTR` non è stato usato. La differenza si vede nelle case e nelle posizioni topocentriche.

### LEB nel contesto

Ogni contesto può avere il proprio file LEB:

```python
from libephemeris import EphemerisContext

ctx = EphemerisContext()
import os
from pathlib import Path

path = Path(os.environ["LIBEPHEMERIS_DATA_DIR"]) / "leb" / "base_core.leb2"
assert path.is_file(), "Prepara prima il gruppo LEB2 base"
ctx.set_leb_file(path)
# La modalità è scelta nel setup: il percorso da solo non la cambia.
position, _ = ctx.calc_ut(2460409.0, ephem.SUN, ephem.FLG_SPEED)
assert len(position) == 6
```

### Chiusura delle risorse condivise

```python
from libephemeris import EphemerisContext

# Chiude file e risorse condivise da tutti i contesti
EphemerisContext.close()
```

---

## 14.7 Download e gestione dati

Prepara esplicitamente i dati necessari prima del calcolo. L'acquisizione dipende da modalità e politica di rete: gli esercizi sealed non devono scaricare file mancanti. I seguenti esempi di preparazione sono operazioni di rete facoltative, separate dal calcolo ordinario.

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

from pathlib import Path

data_directory = ephem.get_library_path()
print(f"Directory dati assoluta: {Path(data_directory).is_absolute()}")
```

```
Directory dati assoluta: True
```

### Reset completo

`close()` rilascia risorse predefinite condivise e ripristina lo stato transitorio dell'osservatore e alcune selezioni. Non cancella ogni impostazione né invalida gli handle conservati. È un'operazione del ciclo di vita, non la pulizia dopo ogni richiesta. Per confrontare i risultati prima/dopo, riapplica la stessa configurazione prima di calcolare:

```python
import libephemeris as ephem

# Fai alcuni calcoli...
jd = ephem.julday(2024, 4, 8, 12.0)
pos, _ = ephem.calc_ut(jd, ephem.SUN, 0)
print(f"Prima di close: {pos[0]:.6f}°")

# Chiudi tutto
ephem.close()

# Il prossimo calcolo ricarica automaticamente le efemeridi
ephem.set_precision_tier("base")
ephem.set_iers_delta_t_enabled(False)
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

- **La precisione dipende dalla grandezza osservabile**, dai dati, dall'epoca, dal modello fisico e dall'osservatore; le cifre stampate non la garantiscono.
- **Tre tier di precisione**: `base` (1849–2150, 31 MB), `medium` (1550–2650, 114 MB), `extended` (-13200 a +17191, 3.1 GB)
- **LEB**: coefficienti locali con copertura e bilancio di interpolazione propri; controlla il file scelto senza presumere risultati identici tra percorsi.
- **Quattro modalità di calcolo**: `auto` (instradamento configurato), `skyfield` (JPL/SPK nativo Rust), `leb` (coefficienti/modelli locali), `horizons` (rete esplicita).
- **EphemerisContext**: contesto thread-safe per calcoli paralleli, con stato isolato per osservatore, modo siderale e cache
- **Preparazione esplicita** con `download_for_tier()` e `download_leb_for_tier()`; i calcoli sealed non acquisiscono dati mancanti.
- **`close()`** rilascia risorse predefinite condivise e ripristina parte dello stato, non ogni impostazione esplicita.

## Funzioni introdotte

- `set_precision_tier(tier)` / `get_precision_tier()` — seleziona il livello di precisione (`"base"`, `"medium"`, `"extended"`)
- `set_calc_mode(mode)` / `get_calc_mode()` — imposta la modalità di calcolo (`"auto"`, `"skyfield"`, `"leb"`, `"horizons"`)
- `set_leb_file(filepath)` — seleziona un file di coefficienti; scegli separatamente la modalità di calcolo
- `download_for_tier(tier)` — scarica tutti i dati per un tier
- `download_leb_for_tier(tier)` — scarica il file LEB precompilato
- `EphemerisContext()` — contesto isolato per calcoli thread-safe
- `EphemerisContext.calc_ut(jd, body, flag)` — calcolo posizione nel contesto
- `EphemerisContext.houses(jd, lat, lon, hsys)` — case nel contesto
- `EphemerisContext.set_topo(lon, lat, alt)` — osservatore nel contesto
- `EphemerisContext.set_sid_mode(mode)` — modo siderale nel contesto
- `EphemerisContext.set_leb_file(filepath)` — LEB per contesto
- `EphemerisContext.close()` — chiude risorse condivise
- `close()` — rilascia risorse condivise e ripristina parte dello stato; handle conservati e impostazioni esplicite possono sopravvivere
- `get_library_path()` — percorso della directory dati
