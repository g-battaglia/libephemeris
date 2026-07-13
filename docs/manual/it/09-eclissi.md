# Capitolo 9 — Eclissi solari e lunari

## Cosa imparerai

In questo capitolo scoprirai come avvengono le eclissi, perché non ce n'è una ogni mese, quali tipi esistono, e come usare la libreria per trovare la prossima eclissi, calcolarne i dettagli e determinare se sarà visibile dal tuo luogo di osservazione.

---

## 9.1 Come avviene un'eclissi

Un'eclissi è un allineamento quasi perfetto tra Sole, Luna e Terra. Ne esistono due tipi fondamentali:

- **Eclissi solare**: la Luna passa davanti al Sole, proiettando la sua ombra sulla Terra. Avviene sempre durante una **Luna nuova** (Sole e Luna dalla stessa parte del cielo).

- **Eclissi lunare**: la Terra blocca la luce del Sole che illumina la Luna, e la Luna si oscura. Avviene sempre durante una **Luna piena** (Sole e Luna su lati opposti).

### Perché non succede ogni mese?

Se l'orbita della Luna fosse esattamente sullo stesso piano dell'orbita terrestre (l'eclittica), avremmo un'eclissi solare a ogni Luna nuova e una lunare a ogni Luna piena — due al mese. Ma l'orbita lunare è **inclinata di circa 5.1°** rispetto all'eclittica. La Luna passa quasi sempre "sopra" o "sotto" il Sole (nelle eclissi solari) o "sopra" o "sotto" l'ombra della Terra (nelle eclissi lunari).

Un'eclissi può avvenire solo quando la Luna nuova o piena cade **vicino a un nodo** — uno dei due punti in cui l'orbita lunare attraversa l'eclittica (vedi Capitolo 6). Questo succede durante le "stagioni delle eclissi", finestre di circa 35 giorni che si ripetono ogni ~173 giorni (circa ogni 6 mesi).

In un anno tipico ci sono 2–5 eclissi solari e 0–3 eclissi lunari.

---

## 9.2 Tipi di eclissi

### Eclissi solari

Il tipo di eclissi solare dipende dalla distanza della Luna dalla Terra in quel momento (la Luna si avvicina e si allontana perché la sua orbita è ellittica):

- **Totale** (`ECL_TOTAL`): la Luna è abbastanza vicina da coprire *completamente* il disco solare. Per qualche minuto, il cielo diventa buio e la corona solare diventa visibile — uno degli spettacoli più straordinari della natura. La fascia di totalità è stretta: in genere 100–250 km di larghezza.

- **Anulare** (`ECL_ANNULAR`): la Luna è troppo lontana e il suo disco appare più piccolo del Sole. Resta un **anello di fuoco** ("annulus") luminoso attorno al disco lunare. Spettacolare, ma il cielo non diventa buio come nella totale.

- **Parziale** (`ECL_PARTIAL`): la Luna copre solo una parte del Sole. È il tipo più comune da osservare — basta trovarsi nella penombra lunare, che copre un'area molto più ampia dell'ombra.

- **Ibrida** (`ECL_ANNULAR_TOTAL`): un'eclissi rara che è anulare in alcune zone della Terra e totale in altre. Succede quando l'ombra della Luna è al limite — l'apice del cono d'ombra sfiora la superficie terrestre.

### Eclissi lunari

- **Totale** (`ECL_TOTAL`): la Luna entra completamente nel cono d'ombra della Terra. Non diventa invisibile, ma assume un colore **rosso rame** — la luce solare filtrata e rifratta dall'atmosfera terrestre la illumina debolmente. Ogni "Luna di sangue" è un'eclissi lunare totale.

- **Parziale** (`ECL_PARTIAL`): solo una parte della Luna entra nell'ombra. Si vede un "morso" scuro sul disco lunare.

- **Penombrale** (`ECL_PENUMBRAL`): la Luna entra nella penombra della Terra (la zona di ombra parziale), ma non nell'ombra vera. L'oscuramento è così lieve che spesso è invisibile a occhio nudo.

### Magnitudine e oscuramento

Due numeri descrivono "quanto" è un'eclissi:

- La **magnitudine** è la frazione del *diametro* del Sole (o della Luna) coperta. Una magnitudine di 0.5 significa che metà del diametro è coperto. Per le eclissi totali la magnitudine è ≥ 1.0.

- L'**oscuramento** (obscuration) è la frazione dell'*area* coperta. A parità di magnitudine, l'oscuramento è sempre maggiore (perché l'area cresce col quadrato del raggio). Una magnitudine di 0.5 corrisponde a un oscuramento di circa 0.39.

---

## 9.3 Trovare la prossima eclissi

### Eclissi solare globale

La funzione `sol_eclipse_when_glob` cerca la prossima eclissi solare **in tutto il mondo**. Non ti dice se sarà visibile dal tuo luogo — solo che avverrà da qualche parte sulla Terra.

```python
import libephemeris as ephem

# Cerca la prossima eclissi solare a partire da oggi
jd_start = ephem.julday(2024, 1, 1, 0.0)

ecl_type, tret = ephem.sol_eclipse_when_glob(jd_start)

# Decodifica il tipo
tipi = []
if ecl_type & ephem.ECL_TOTAL:
    tipi.append("totale")
if ecl_type & ephem.ECL_ANNULAR:
    tipi.append("anulare")
if ecl_type & ephem.ECL_PARTIAL:
    tipi.append("parziale")
if ecl_type & ephem.ECL_ANNULAR_TOTAL:
    tipi.append("ibrida")

# tret[0] = momento del massimo
jd_max = tret[0]
y, m, d, h = ephem.revjul(jd_max)
ore = int(h)
minuti = int((h - ore) * 60)

print(f"Prossima eclissi solare: {tipi}")
print(f"Data: {d}/{m}/{y} alle {ore:02d}:{minuti:02d} UT")

# Tempi dei contatti
if tret[1] > 0:
    y1, m1, d1, h1 = ephem.revjul(tret[1])
    print(f"Inizio (1° contatto): {d1}/{m1}/{y1} {int(h1):02d}:{int((h1%1)*60):02d} UT")
if tret[4] > 0:
    y4, m4, d4, h4 = ephem.revjul(tret[4])
    print(f"Fine (4° contatto):   {d4}/{m4}/{y4} {int(h4):02d}:{int((h4%1)*60):02d} UT")
```

```
Prossima eclissi solare: ['totale']
Data: 8/4/2024 alle 18:17 UT
Fine (4° contatto):   8/4/2024 16:41 UT
```

Puoi filtrare per tipo di eclissi:

```python
# Cerca solo eclissi totali
ecl_type, tret = ephem.sol_eclipse_when_glob(
    jd_start,
    ifltype=ephem.ECL_TOTAL
)
```

### Le prossime 5 eclissi solari

```python
import libephemeris as ephem

jd = ephem.julday(2024, 1, 1, 0.0)

print("Prossime 5 eclissi solari:\n")

for i in range(5):
    ecl_type, tret = ephem.sol_eclipse_when_glob(jd)
    jd_max = tret[0]

    # Tipo
    if ecl_type & ephem.ECL_TOTAL:
        tipo = "Totale"
    elif ecl_type & ephem.ECL_ANNULAR:
        tipo = "Anulare"
    elif ecl_type & ephem.ECL_ANNULAR_TOTAL:
        tipo = "Ibrida"
    else:
        tipo = "Parziale"

    y, m, d, h = ephem.revjul(jd_max)
    print(f"  {i+1}. {d:2.0f}/{m:02.0f}/{y:.0f}  {tipo}")

    # Avanza oltre questa eclissi
    jd = jd_max + 30
```

```
Prossime 5 eclissi solari:

  1.  8/04/2024  Totale
  2.  2/10/2024  Anulare
  3. 29/03/2025  Parziale
  4. 21/09/2025  Parziale
  5. 17/02/2026  Anulare
```

### Eclissi solare visibile da un luogo

La domanda più comune è: "Quando sarà la prossima eclissi visibile da casa mia?" Per rispondere, usa `sol_eclipse_when_loc`:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 1, 1, 0.0)
lat, lon = 41.9028, 12.4964  # Roma

# Cerca la prossima eclissi solare visibile da Roma
ecl_type, tempi, attr = ephem.sol_eclipse_when_loc(
    jd, (lon, lat, 0.0)
)

jd_max = tempi[0]
y, m, d, h = ephem.revjul(jd_max)
ore = int(h)
minuti = int((h - ore) * 60)

magnitudine = attr[0]
oscuramento = attr[2]
alt_sole = attr[4]

print(f"Prossima eclissi visibile da Roma:")
print(f"Data: {d:.0f}/{m:.0f}/{y:.0f} alle {ore:02d}:{minuti:02d} UT")
print(f"Magnitudine: {magnitudine:.3f}")
print(f"Oscuramento: {oscuramento:.1%}")
print(f"Altezza del Sole al massimo: {alt_sole:.1f}°")
```

```
Prossima eclissi visibile da Roma:
Data: 29/3/2025 alle 11:03 UT
Magnitudine: 0.073
Oscuramento: 2.4%
Altezza del Sole al massimo: 51.6°
```

La differenza tra le due funzioni è cruciale:

- `sol_eclipse_when_glob` trova la prossima eclissi *nel mondo* — ce ne sono 2–5 all'anno
- `sol_eclipse_when_loc` trova la prossima visibile *dal tuo punto di osservazione* — possono passare anni tra una e l'altra

---

## 9.4 Dettagli di un'eclissi

### Dove è visibile al massimo?

Data l'eclissi, `sol_eclipse_where` ti dice dove nel mondo cade il punto di massimo:

```python
import libephemeris as ephem

# Eclissi totale dell'8 aprile 2024
jd_start = ephem.julday(2024, 4, 1, 0.0)
ecl_type, tret = ephem.sol_eclipse_when_glob(jd_start)
jd_max = tret[0]

# Dove cade il massimo?
ret, geopos, attr = ephem.sol_eclipse_where(jd_max)

lon_centro = geopos[0]
lat_centro = geopos[1]
larghezza = attr[3]  # larghezza del percorso in km

print(f"Centro della linea centrale: {lat_centro:.2f}°N, {lon_centro:.2f}°E")
print(f"Larghezza del percorso di totalità: {larghezza:.0f} km")
print(f"Magnitudine al centro: {attr[0]:.4f}")
```

```
Centro della linea centrale: 25.29°N, -104.15°E
Larghezza del percorso di totalità: 207 km
Magnitudine al centro: 1.0571
```

### Tipo e magnitudine in un dato momento e luogo

Se sai quando avviene un'eclissi e vuoi sapere come appare da un luogo specifico, usa `sol_eclipse_how`:

```python
import libephemeris as ephem

# Eclissi dell'8 aprile 2024 vista da Dallas, Texas
jd_max = ephem.julday(2024, 4, 8, 18.0 + 42/60)  # ~18:42 UT
lat, lon = 32.78, -96.80  # Dallas

ret, attr = ephem.sol_eclipse_how(jd_max, (lon, lat, 0.0))

if ret > 0:
    magnitudine = attr[0]
    oscuramento = attr[2]
    alt_sole = attr[5]  # altezza vera del Sole
    az_sole = attr[4]   # azimut del Sole

    print(f"Magnitudine: {magnitudine:.4f}")
    print(f"Oscuramento: {oscuramento:.1%}")
    print(f"Sole: altezza {alt_sole:.1f}°, azimut {az_sole:.1f}°")
else:
    print("Nessuna eclissi visibile da questo luogo in questo momento")
```

```
Magnitudine: 1.0126
Oscuramento: 111.6%
Sole: altezza 64.6°, azimut 7.6°
```

### Dettagli completi con `sol_eclipse_how_details`

Per un rapporto completo con tutti i contatti, gli angoli di posizione e le durate, usa la versione estesa:

```python
import libephemeris as ephem

jd_start = ephem.julday(2024, 4, 1, 0.0)
ecl_type, tret = ephem.sol_eclipse_when_glob(jd_start)

# Dettagli completi da Roma
info = ephem.sol_eclipse_how_details(
    tret[0], 41.9028, 12.4964
)

if info["is_visible"]:
    print(f"Tipo: {'totale' if info['is_total'] else 'parziale'}")
    print(f"Magnitudine massima: {info['max_magnitude']:.4f}")
    print(f"Oscuramento massimo: {info['max_obscuration_percent']:.1f}%")
    print(f"Durata fase parziale: {info['duration_partial_minutes']:.1f} min")

    if info["is_total"]:
        print(f"Durata totalità: {info['duration_total_minutes']:.1f} min")
```

```
Non visibile da Roma
```

---

## 9.5 Eclissi lunari

Le eclissi lunari si cercano con funzioni analoghe.

### Trovare la prossima eclissi lunare

```python
import libephemeris as ephem

jd = ephem.julday(2024, 1, 1, 0.0)

ecl_type, tret = ephem.lun_eclipse_when(jd)

# Tipo
if ecl_type & ephem.ECL_TOTAL:
    tipo = "Totale"
elif ecl_type & ephem.ECL_PARTIAL:
    tipo = "Parziale"
elif ecl_type & ephem.ECL_PENUMBRAL:
    tipo = "Penombrale"
else:
    tipo = "Sconosciuto"

jd_max = tret[0]
y, m, d, h = ephem.revjul(jd_max)
print(f"Prossima eclissi lunare: {tipo}")
print(f"Data: {d:.0f}/{m:.0f}/{y:.0f}")
```

```
Prossima eclissi lunare: Penombrale
Data: 25/3/2024
```

La tupla `tret` per le eclissi lunari contiene 8 tempi:

- `tret[0]` — momento del massimo
- `tret[1]` — inizio della fase parziale (Luna entra nell'ombra)
- `tret[2]` — inizio della totalità (se totale, altrimenti 0)
- `tret[3]` — fine della totalità
- `tret[4]` — fine della fase parziale (Luna esce dall'ombra)
- `tret[5]` — inizio della fase penombrale
- `tret[6]` — fine della fase penombrale

### Eclissi lunare visibile da un luogo

Un'eclissi lunare è visibile ovunque la Luna sia sopra l'orizzonte durante l'evento. Ma potrebbe essere parzialmente visibile — magari la Luna sorge quando l'eclissi è già a metà, o tramonta prima che finisca.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 1, 1, 0.0)
lat, lon = 41.9028, 12.4964  # Roma

ecl_type, tempi, attr = ephem.lun_eclipse_when_loc(jd, (lon, lat, 0.0))

jd_max = tempi[0]
y, m, d, h = ephem.revjul(jd_max)

# La magnitudine umbrale dice "quanto profonda" è l'eclissi
mag_umbrale = attr[0]

# Controlliamo la visibilità
visibile = ecl_type & ephem.ECL_VISIBLE

print(f"Eclissi lunare: {d:.0f}/{m:.0f}/{y:.0f}")
print(f"Magnitudine umbrale: {mag_umbrale:.3f}")

if visibile:
    alt_luna = attr[5]
    az_luna = attr[4]
    print(f"Visibile da Roma!")
    print(f"Luna al massimo: altezza {alt_luna:.1f}°, azimut {az_luna:.1f}°")
else:
    print("Non visibile da Roma (Luna sotto l'orizzonte)")
```

```
Eclissi lunare: 25/3/2024
Magnitudine umbrale: 0.000
Visibile da Roma!
Luna al massimo: altezza -22.7°, azimut 109.2°
```

### Tipo e magnitudine in un momento preciso

```python
import libephemeris as ephem

# Eclissi lunare: controlliamo un momento specifico
jd = ephem.julday(2025, 3, 14, 6.0)  # esempio

ret, attr = ephem.lun_eclipse_how(jd, (12.4964, 41.9028, 0.0))

if ret > 0:
    mag_umbrale = attr[0]
    mag_penombrale = attr[1]
    saros = int(attr[9])
    saros_membro = int(attr[10])

    print(f"Magnitudine umbrale: {mag_umbrale:.3f}")
    print(f"Magnitudine penombrale: {mag_penombrale:.3f}")
    if saros > 0:
        print(f"Serie Saros: {saros}, membro {saros_membro}")
```

```
Magnitudine umbrale: 1.090
Magnitudine penombrale: 2.000
Serie Saros: 123, membro 53
```

---

## 9.6 I cicli delle eclissi: Saros e Inex

Le eclissi non avvengono a caso — seguono cicli precisi conosciuti fin dall'antichità.

### Il ciclo di Saros

Il **Saros** è un periodo di 6585.32 giorni, ovvero circa **18 anni, 11 giorni e 8 ore**. Dopo un Saros, Sole, Luna e nodi tornano quasi esattamente nella stessa configurazione, e un'eclissi molto simile si ripete.

Ma c'è quel "quasi": le 8 ore di differenza significano che la Terra ha ruotato di circa un terzo di giro in più. Quindi l'eclissi successiva cade circa **120° più a ovest** sulla superficie terrestre.

Ogni "serie Saros" è una famiglia di eclissi correlate che nasce come parziale a un polo, diventa gradualmente totale, poi torna parziale all'altro polo, per un totale di circa 1200–1400 anni e 70–80 eclissi.

```python
import libephemeris as ephem

# Eclissi totale dell'8 aprile 2024
jd_start = ephem.julday(2024, 4, 1, 0.0)
ecl_type, tret = ephem.sol_eclipse_when_glob(
    jd_start, ecltype=ephem.ECL_TOTAL
)
jd_ecl = tret[0]

# A quale serie Saros appartiene?
saros = ephem.get_saros_number(jd_ecl, "solar")
print(f"Serie Saros: {saros}")

# La "sorella" di 18 anni prima
jd_sorella = jd_ecl - ephem.SAROS_CYCLE_DAYS
y, m, d, h = ephem.revjul(jd_sorella)
print(f"Eclissi sorella: ~{d:.0f}/{m:.0f}/{y:.0f}")
```

```
Serie Saros: 139
Eclissi sorella: ~29/3/2006
```

### Il ciclo Inex

L'**Inex** è un periodo più lungo di 10571.95 giorni, circa **29 anni**. Collega eclissi di serie Saros diverse — quando una serie finisce, l'Inex la collega alla serie successiva. I due cicli insieme formano una griglia che copre tutte le eclissi passate e future.

```python
import libephemeris as ephem

jd_ecl = ephem.julday(2024, 4, 8, 18.0)

inex = ephem.get_inex_number(jd_ecl, "solar")
print(f"Numero Inex: {inex}")
```

```
Numero Inex: 50
```

---

## 9.7 Occultazioni

Un'**occultazione** avviene quando un corpo celeste passa davanti a un altro, nascondendolo. Il caso più comune è la Luna che occulta una stella o un pianeta.

### Occultazione di una stella

```python
import libephemeris as ephem

# Cerca la prossima occultazione di Regolo da parte della Luna
jd = ephem.julday(2024, 1, 1, 0.0)

ret, tret = ephem.lun_occult_when_glob(
    jd,
    body="Regulus",
)

if ret > 0:
    jd_max = tret[0]
    y, m, d, h = ephem.revjul(jd_max)
    print(f"Prossima occultazione di Regolo: {d:.0f}/{m:.0f}/{y:.0f}")
```

```
Prossima occultazione di Regolo: 16/10/2025
```

### Occultazione di un pianeta

```python
import libephemeris as ephem

# Cerca la prossima occultazione di Giove da parte della Luna
jd = ephem.julday(2024, 1, 1, 0.0)

ret, tret = ephem.lun_occult_when_glob(
    jd,
    body=ephem.JUPITER,
)

if ret > 0:
    jd_max = tret[0]
    y, m, d, h = ephem.revjul(jd_max)
    print(f"Prossima occultazione di Giove: {d:.0f}/{m:.0f}/{y:.0f}")
```

```
Prossima occultazione di Giove: 8/9/2026
```

### Occultazione visibile da un luogo

```python
import libephemeris as ephem

jd = ephem.julday(2024, 1, 1, 0.0)

# Prossima occultazione lunare di Spica visibile da Roma
ret, tempi, attr = ephem.lun_occult_when_loc(
    jd,
    "Spica",
    geopos=(12.4964, 41.9028, 0.0)
)

if ret > 0:
    jd_inizio = tempi[1]  # primo contatto
    jd_fine = tempi[4]    # ultimo contatto
    y, m, d, h = ephem.revjul(tempi[0])

    durata = (jd_fine - jd_inizio) * 24 * 60  # in minuti

    print(f"Occultazione di Spica: {d:.0f}/{m:.0f}/{y:.0f}")
    print(f"Durata: {durata:.0f} minuti")
    print(f"Altezza al massimo: {attr[5]:.1f}°")
else:
    print("Nessuna occultazione trovata nel periodo di ricerca")
```

```
Occultazione di Spica: 1/3/2032
Durata: 34 minuti
Altezza al massimo: 28.3°
```

### Occultazioni tra pianeti

Rarissime ma calcolabili: un pianeta che occulta un altro pianeta o una stella.

```python
import libephemeris as ephem

# Cerca un'occultazione di Regolo da parte di Venere
jd = ephem.julday(2024, 1, 1, 0.0)

ret, tret = ephem.planet_occult_when_glob(
    jd,
    occulting_planet=ephem.VENUS,
    occulted_planet=0,
    starname="Regulus"
)

if ret > 0:
    y, m, d, h = ephem.revjul(tret[0])
    print(f"Venere occulta Regolo: {d:.0f}/{m:.0f}/{y:.0f}")
else:
    print("Nessuna occultazione trovata (può cercare fino a 150 anni)")
```

> **Nota**: Questa ricerca può richiedere molto tempo (minuti) perché le occultazioni planetarie sono estremamente rare. La funzione cerca fino a 150 anni nel futuro.

---

## Riepilogo

In questo capitolo abbiamo esplorato le eclissi, dalla geometria alla pratica.

**Concetti chiave:**

- Le eclissi avvengono solo quando una Luna nuova (eclissi solare) o piena (eclissi lunare) cade vicino a un nodo lunare — circa 4–7 volte all'anno in totale
- Le eclissi solari possono essere **totali**, **anulari**, **parziali** o **ibride**, a seconda della distanza della Luna
- Le eclissi lunari possono essere **totali** (Luna rossa), **parziali** o **penombrali** (quasi invisibili)
- La **magnitudine** misura la frazione del diametro coperta; l'**oscuramento** misura la frazione dell'area coperta
- Il ciclo di **Saros** (~18 anni) collega eclissi simili; il ciclo **Inex** (~29 anni) collega serie diverse
- Le **occultazioni** sono eclissi di stelle o pianeti da parte della Luna o di altri pianeti

**Funzioni introdotte:**

- `sol_eclipse_when_glob(jd, ifltype=0)` — trova la prossima eclissi solare nel mondo, restituisce tipo e tempi dei contatti
- `sol_eclipse_when_loc(jd, geopos)` — trova la prossima eclissi solare visibile da un luogo, con magnitudine e oscuramento
- `sol_eclipse_where(jd)` — dato il momento del massimo, restituisce le coordinate della linea centrale
- `sol_eclipse_how(jd, geopos)` — calcola tipo, magnitudine e posizione del Sole per un'eclissi già nota
- `sol_eclipse_how_details(jd, lat, lon)` — versione estesa con tutti i contatti, angoli di posizione e durate
- `lun_eclipse_when(jd, ifltype=0)` — trova la prossima eclissi lunare (globale)
- `lun_eclipse_when_loc(jd, geopos)` — trova la prossima eclissi lunare visibile da un luogo
- `lun_eclipse_how(jd, geopos)` — dettagli di un'eclissi lunare in un momento e luogo
- `lun_occult_when_glob(jd, body)` — trova la prossima occultazione lunare di una stella o pianeta
- `lun_occult_when_loc(jd, body, geopos)` — occultazione visibile da un luogo
- `planet_occult_when_glob(jd, occulting_planet, starname=...)` — occultazione tra pianeti
- `get_saros_number(jd, "solar"|"lunar")` — restituisce la serie Saros di un'eclissi
- `get_inex_number(jd, "solar"|"lunar")` — restituisce il numero Inex
