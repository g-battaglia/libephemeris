# Capitolo 7 — Le case astrologiche

## Cosa imparerai

In questo capitolo scoprirai cosa sono le case astrologiche, perché esistono oltre 20 sistemi diversi per calcolarle, come funzionano l'Ascendente e il Medio Cielo, come determinare in quale casa cade un pianeta, e cosa succede alle latitudini estreme dove alcuni sistemi non funzionano.

---

## 7.1 Cosa sono le case

Immagina di prendere la sfera celeste e di dividerla in 12 fette, come una torta. Ogni fetta è una **casa**. Mentre i segni zodiacali dividono l'eclittica in base al percorso annuale del Sole (e quindi dipendono solo dalla data), le case dividono il cielo in base alla rotazione quotidiana della Terra — e quindi dipendono dalla **data**, dall'**ora** e dal **luogo**.

Ogni casa "governa" un ambito della vita secondo la tradizione astrologica:

- **1a casa** — il sé, l'aspetto fisico, il temperamento
- **2a casa** — le risorse, il denaro, i valori personali
- **3a casa** — la comunicazione, i fratelli, i viaggi brevi
- **4a casa** — la famiglia, le radici, la casa fisica
- **5a casa** — la creatività, i figli, il piacere
- **6a casa** — il lavoro quotidiano, la salute, il servizio
- **7a casa** — le relazioni, i partner, i contratti
- **8a casa** — le trasformazioni, l'eredità, la sessualità
- **9a casa** — i viaggi lunghi, la filosofia, la spiritualità
- **10a casa** — la carriera, la reputazione, l'ambizione
- **11a casa** — le amicizie, i gruppi, i progetti futuri
- **12a casa** — l'inconscio, l'isolamento, la trascendenza

La differenza fondamentale tra segni e case: due persone nate lo stesso giorno hanno i pianeti negli **stessi segni**, ma se sono nate a ore o in luoghi diversi, i pianeti cadono in **case diverse**. L'ora di nascita determina le case.

Il punto di inizio di ogni casa si chiama **cuspide**. La cuspide della 1a casa è l'Ascendente.

---

## 7.2 L'Ascendente e il Medio Cielo

Abbiamo già incontrato questi concetti nel Capitolo 1, ma qui approfondiamo il loro ruolo nel sistema delle case.

L'**Ascendente** (ASC) è il grado dell'eclittica che sta sorgendo all'orizzonte orientale in un dato momento e luogo. È la cuspide della 1a casa in tutti i sistemi di case. Cambia di circa 1° ogni 4 minuti — ecco perché l'ora esatta di nascita è così importante.

Il **Medio Cielo** (MC) è il grado dell'eclittica che culmina al meridiano — il punto più alto che quel grado raggiunge nel cielo. È la cuspide della 10a casa nella maggior parte dei sistemi (ma non in tutti: il sistema a segni interi e il sistema uguale non usano il MC come cuspide della 10a).

La funzione `houses` restituisce due gruppi di dati:

- Le **12 cuspidi** delle case (una longitudine eclittica per ogni casa)
- Gli **angoli**: Ascendente, Medio Cielo, ARMC, Vertice e altri punti speciali

```python
import libephemeris as ephem

# 8 aprile 2024, ore 14:30 UT — Roma
jd = ephem.julday(2024, 4, 8, 14.5)
lat, lon = 41.9028, 12.4964

# Calcola le case con il sistema Placidus
cusps, ascmc = ephem.houses(jd, lat, lon, ord('P'))

# Le 12 cuspidi
segni = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

for i in range(12):
    lon_c = cusps[i]
    s = segni[int(lon_c / 30)]
    g = lon_c % 30
    print(f"Casa {i+1:2d}: {g:5.1f}° {s}")

# Gli angoli principali
print(f"\nAscendente:  {ascmc[0]:.4f}°")
print(f"Medio Cielo: {ascmc[1]:.4f}°")
print(f"ARMC:        {ascmc[2]:.4f}°")  # Asc. Retta del MC
print(f"Vertice:     {ascmc[3]:.4f}°")
```

```
Casa  1:  12.2° Vir
Casa  2:   6.3° Lib
Casa  3:   5.5° Sco
Casa  4:   9.0° Sgr
Casa  5:  13.5° Cap
Casa  6:  15.0° Aqr
Casa  7:  12.2° Psc
Casa  8:   6.3° Ari
Casa  9:   5.5° Tau
Casa 10:   9.0° Gem
Casa 11:  13.5° Cnc
Casa 12:  15.0° Leo

Ascendente:  162.2479°
Medio Cielo: 69.0379°
ARMC:        67.3366°
Vertice:     316.3110°
```

L'**ARMC** (Ascensione Retta del Medio Cielo) è lo stesso concetto del tempo siderale locale, espresso in gradi invece che in ore (ARMC in gradi = tempo siderale in ore × 15). È il dato di partenza per il calcolo di tutti i sistemi di case.

Il **Vertice** è il punto dove il primo verticale (il cerchio che passa per est, zenit e ovest) interseca l'eclittica dal lato ovest. In alcune scuole astrologiche ha un significato legato agli incontri fatali e alle relazioni karmiche.

---

## 7.3 I sistemi di case: quale scegliere?

Perché esistono oltre 20 sistemi di case? Perché non c'è un modo univoco di dividere una sfera tridimensionale in 12 settori. Ogni sistema usa un criterio diverso — dividere il tempo, dividere lo spazio, dividere l'equatore — e ognuno produce risultati leggermente diversi.

Ecco i sistemi più usati, con la lettera da passare a `houses()`:

**Placidus** (`P`) è il sistema più diffuso in Occidente dal XVII secolo. Divide il **tempo** che ogni grado dell'eclittica impiega per andare dall'orizzonte al meridiano. È basato sul concetto di "semi-arco": il tempo tra la levata e la culminazione di un punto. I suoi sostenitori lo considerano il più "naturale" perché riflette il moto reale del cielo. **Problema**: non funziona oltre il circolo polare (~66.5° di latitudine).

**Koch** (`K`) è simile a Placidus ma usa il semi-arco del Medio Cielo invece che del punto stesso. Popolare in Germania. Stesso problema alle latitudini estreme.

**Regiomontanus** (`R`) divide l'**equatore celeste** in 12 parti uguali di 30° ciascuna, poi proietta queste divisioni sull'eclittica. Fu il sistema dominante nel Medioevo e nel Rinascimento.

**Campanus** (`C`) divide il **primo verticale** (il cerchio est-zenit-ovest) in 12 parti uguali. È uno dei sistemi più antichi con una base geometrica chiara.

**Uguale dall'Ascendente** (`E`) è il più semplice: ogni casa ha esattamente 30°, partendo dall'Ascendente. La 1a casa va dall'Ascendente all'Ascendente + 30°, la 2a da +30° a +60°, e così via. Il MC in questo sistema **non** è necessariamente la cuspide della 10a casa — può cadere in qualsiasi casa tra l'8a e l'11a. Molto usato nell'astrologia antica e nella tradizione ellenistica.

**Segno intero** (`W`, Whole Sign) è ancora più semplice: ogni segno zodiacale è una casa. Se l'Ascendente è in Leone, tutta la 1a casa coincide con il Leone (0°–30° Leone), la 2a con la Vergine, e così via. È il sistema più antico, usato nell'astrologia greca e indiana. Sta vivendo una forte rinascita nell'astrologia moderna.

**Porfirio** (`O`) divide proporzionalmente i quattro quadranti (ASC-IC, IC-DSC, DSC-MC, MC-ASC). È un buon compromesso e funziona a **tutte le latitudini** — per questo è spesso usato come fallback quando Placidus o Koch falliscono.

**Polich/Page** (`T`, Topocentriche) è una proiezione topografica. Molto popolare in Sud America.

**Morino** (`M`) divide l'equatore celeste dal MC. Funziona a tutte le latitudini.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 14.5)
lat, lon = 41.9028, 12.4964

# Confronto: Placidus vs Uguale vs Segno intero
for lettera, nome in [('P', 'Placidus'), ('E', 'Uguale'), ('W', 'Segno intero')]:
    cusps, ascmc = ephem.houses(jd, lat, lon, ord(lettera))
    print(f"\n{nome} (cuspide 1a = ASC = {ascmc[0]:.1f}°):")
    for i in range(12):
        print(f"  Casa {i+1:2d}: {cusps[i]:.1f}°")

# Ottenere il nome del sistema da una lettera
print(ephem.house_name(ord('P')))  # "Placidus"
print(ephem.house_name(ord('W')))  # "Whole Sign"
```

```
Placidus (cuspide 1a = ASC = 162.2°):
  Casa  1: 162.2°
  Casa  2: 186.3°
  Casa  3: 215.5°
  Casa  4: 249.0°
  Casa  5: 283.5°
  Casa  6: 315.0°
  Casa  7: 342.2°
  Casa  8: 6.3°
  Casa  9: 35.5°
  Casa 10: 69.0°
  Casa 11: 103.5°
  Casa 12: 135.0°

Uguale (cuspide 1a = ASC = 162.2°):
  Casa  1: 162.2°
  Casa  2: 192.2°
  Casa  3: 222.2°
  Casa  4: 252.2°
  Casa  5: 282.2°
  Casa  6: 312.2°
  Casa  7: 342.2°
  Casa  8: 12.2°
  Casa  9: 42.2°
  Casa 10: 72.2°
  Casa 11: 102.2°
  Casa 12: 132.2°

Segno intero (cuspide 1a = ASC = 162.2°):
  Casa  1: 150.0°
  Casa  2: 180.0°
  Casa  3: 210.0°
  Casa  4: 240.0°
  Casa  5: 270.0°
  Casa  6: 300.0°
  Casa  7: 330.0°
  Casa  8: 0.0°
  Casa  9: 30.0°
  Casa 10: 60.0°
  Casa 11: 90.0°
  Casa 12: 120.0°

Placidus
Whole Sign
```

---

## 7.4 In quale casa cade un pianeta?

Sapere le cuspidi delle case è solo metà del lavoro. La domanda più frequente è: **in quale casa si trova il mio Sole? E la mia Luna?** Per rispondere serve la funzione `house_pos`.

### Il concetto

Determinare la casa di un pianeta non è banale come sembra. Non basta confrontare la longitudine del pianeta con le cuspidi. Perché? Perché molti sistemi di case (Placidus, Koch, Regiomontanus) lavorano in **tre dimensioni**: non dividono solo l'eclittica, ma lo spazio tridimensionale. Un pianeta con latitudine eclittica diversa da zero potrebbe cadere in una casa diversa rispetto a un punto con la stessa longitudine ma latitudine zero.

La funzione `house_pos` fa il calcolo corretto per il sistema scelto. Restituisce un numero decimale dove:

- La **parte intera** è il numero della casa (1–12)
- La **parte decimale** indica quanto il pianeta è "avanti" nella casa (0.0 = sulla cuspide, 0.99 = quasi alla cuspide successiva)

Per esempio, un valore di `7.50` significa "esattamente a metà della 7a casa".

### Come usare `house_pos`

La funzione ha bisogno di:

- **ARMC**: l'Ascensione Retta del Medio Cielo (la ottieni da `houses()` come `ascmc[2]`)
- **Latitudine geografica** dell'osservatore
- **Obliquità dell'eclittica** (la puoi ottenere con `calc_ut(jd, SE_ECL_NUT)`)
- **Sistema di case** (la lettera, come `ord('P')`)
- **Longitudine e latitudine eclittica** del pianeta

```python
import libephemeris as ephem

# 8 aprile 2024, ore 14:30 UT — Roma
jd = ephem.julday(2024, 4, 8, 14.5)
lat, lon = 41.9028, 12.4964

# 1. Calcola le case per ottenere ARMC e cuspidi
cusps, ascmc = ephem.houses(jd, lat, lon, ord('P'))
armc = ascmc[2]  # ARMC in gradi

# 2. Ottieni l'obliquità dell'eclittica
nut, _ = ephem.calc_ut(jd, ephem.ECL_NUT, 0)
obliquity = nut[0]  # obliquità vera

# 3. Calcola la posizione del Sole
sun, _ = ephem.calc_ut(jd, ephem.SUN, ephem.FLG_SPEED)
sun_lon = sun[0]   # longitudine eclittica
sun_lat = sun[1]   # latitudine eclittica (quasi zero per il Sole)

# 4. Determina in quale casa cade
pos = ephem.house_pos(armc, lat, obliquity, ord('P'), sun_lon, sun_lat)
casa = int(pos)
posizione = pos - casa  # quanto "avanti" nella casa (0.0 - 0.99)

print(f"Il Sole è nella {casa}a casa")
print(f"Posizione nella casa: {posizione:.2%}")
```

```
Il Sole è nella 8a casa
Posizione nella casa: 46.37%
```

### Tutti i pianeti nelle loro case

Ecco un esempio completo che mostra ogni pianeta con il segno, i gradi e la casa:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 14.5)
lat, lon = 41.9028, 12.4964

# Calcola case e obliquità
cusps, ascmc = ephem.houses(jd, lat, lon, ord('P'))
armc = ascmc[2]
nut, _ = ephem.calc_ut(jd, ephem.ECL_NUT, 0)
obliquity = nut[0]

segni = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

pianeti = [
    (ephem.SUN,      "Sole    "),
    (ephem.MOON,     "Luna    "),
    (ephem.MERCURY,  "Mercurio"),
    (ephem.VENUS,    "Venere  "),
    (ephem.MARS,     "Marte   "),
    (ephem.JUPITER,  "Giove   "),
    (ephem.SATURN,   "Saturno "),
]

for body_id, nome in pianeti:
    pos, _ = ephem.calc_ut(jd, body_id, ephem.FLG_SPEED)
    p_lon, p_lat = pos[0], pos[1]

    # Segno e gradi
    segno = segni[int(p_lon / 30)]
    gradi = p_lon % 30

    # Casa
    hp = ephem.house_pos(armc, lat, obliquity, ord('P'), p_lon, p_lat)
    casa = int(hp)

    print(f"{nome}  {gradi:5.1f}° {segno}  →  casa {casa:2d}")
```

```
Sole       19.2° Ari  →  casa  8
Luna       17.0° Ari  →  casa  8
Mercurio   24.9° Ari  →  casa  8
Venere      4.2° Ari  →  casa  7
Marte      12.9° Psc  →  casa  7
Giove      19.0° Tau  →  casa  9
Saturno    14.4° Psc  →  casa  7
```

### La latitudine eclittica conta?

Per la maggior parte dei pianeti la latitudine eclittica è piccola (il Sole ha latitudine quasi zero per definizione, i pianeti restano entro pochi gradi dall'eclittica). Ma la Luna può raggiungere ±5.1°, e Plutone fino a ±17°. In questi casi, ignorare la latitudine può spostare il pianeta in una casa diversa, specialmente vicino alle cuspidi.

Se passi `lat_body=0.0`, stai calcolando la casa basandoti solo sulla longitudine — come se il pianeta fosse esattamente sull'eclittica. Questo è il metodo tradizionale usato da molti software. Ma se vuoi il calcolo tridimensionale corretto, usa la latitudine eclittica vera del pianeta.

### Settori Gauquelin

Michel Gauquelin (1928–1991) fu uno psicologo e statistico francese che studiò se le posizioni planetarie alla nascita fossero correlate con la professione. Scoprì che certi pianeti tendono a trovarsi subito dopo la levata o subito dopo la culminazione nei temi di persone di successo (l'"effetto Marte" per gli sportivi, l'"effetto Giove" per gli attori).

Per questa analisi, lo spazio attorno all'osservatore è diviso in **36 settori** invece di 12. Puoi calcolarli con:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 14.5)
lat, lon = 41.9028, 12.4964

# Settore Gauquelin di Marte
settore = ephem.gauquelin_sector(
    jd, ephem.MARS, 0,
    geopos=(lon, lat, 0.0)
)

num = int(settore)
print(f"Marte è nel settore Gauquelin {num}")
print(f"Settore 1 = levata, 10 = culminazione, 19 = tramonto, 28 = IC")
```

```
Marte è nel settore Gauquelin 18
Settore 1 = levata, 10 = culminazione, 19 = tramonto, 28 = IC
```

I settori si contano in senso orario dall'Ascendente: il settore 1 è la levata (Ascendente), il 10 è la culminazione superiore (MC), il 19 è il tramonto (Discendente), il 28 è la culminazione inferiore (IC). Le "zone Gauquelin" più significative sono i settori subito dopo la levata (1–3) e subito dopo la culminazione (10–12).

---

## 7.5 Latitudini estreme e circolo polare

### Il problema

Fin qui abbiamo lavorato con Roma (41.9° N), una latitudine dove tutti i sistemi di case funzionano senza problemi. Ma cosa succede se devi calcolare il tema natale di qualcuno nato a Tromsø (69.6° N), a Reykjavík (64.1° N), o a Murmansk (68.9° N)?

Il problema è geometrico. Sistemi come Placidus e Koch dividono il cielo in base al **tempo** che un grado dell'eclittica impiega per andare dall'orizzonte al meridiano. Ma oltre il circolo polare (~66.56° di latitudine), alcune parti dell'eclittica **non sorgono e non tramontano mai** — restano sempre sopra o sempre sotto l'orizzonte. Se un punto non sorge mai, non puoi misurare "quanto tempo impiega a sorgere", e il calcolo diventa impossibile.

La soglia esatta dipende dall'obliquità dell'eclittica (circa 23.44° nell'epoca attuale):

- **Soglia polare** = 90° − obliquità ≈ 90° − 23.44° = **66.56°**
- Sopra questa latitudine, Placidus (`P`), Koch (`K`) e Gauquelin (`G`) **non funzionano**

### Quali sistemi funzionano dove

Non tutti i sistemi hanno questo problema. Ecco la situazione:

**Sistemi che falliscono oltre il circolo polare** — Si basano sul semi-arco (il tempo tra levata e culminazione), che non esiste per i punti circumpolari:

- **Placidus** (`P`) — il più diffuso, ma non funziona oltre ~66.56°
- **Koch** (`K`) — stesso problema
- **Gauquelin** (`G`) — i 36 settori hanno lo stesso limite

**Sistemi instabili alle latitudini estreme (oltre 80°)** — Funzionano tecnicamente, ma possono dare risultati numericamente imprecisi o "schiacciati":

- Campanus (`C`), Regiomontanus (`R`), Polich/Page (`T`), Alcabizio (`B`), Horizon (`H`), Krusinski (`U`), APC (`Y`), Carter (`F`)

**Sistemi stabili a tutte le latitudini** — Usano metodi geometrici che non dipendono dalla levata e dal tramonto:

- Uguale (`E`), Segno intero (`W`), Porfirio (`O`), Morino (`M`), Meridiano (`X`), Vehlow (`V`), Axial rotation (`N`)

Puoi verificare la situazione per qualsiasi latitudine con `get_extreme_latitude_info`:

```python
import libephemeris as ephem

# Tromsø, Norvegia — dentro il circolo polare
info = ephem.get_extreme_latitude_info(69.6)

print(f"Latitudine: {info['latitude']}°")
print(f"È estrema (>80°)?     {info['is_extreme']}")
print(f"È polare (>{info['polar_threshold']:.1f}°)? {info['is_polar_circle']}")
print(f"Sistemi che falliscono:  {info['affected_systems']}")
print(f"Sistemi instabili:       {info['unstable_systems']}")
print(f"Sistemi sempre stabili:  {info['stable_systems']}")
```

```
Latitudine: 69.6°
È estrema (>80°)?     False
È polare (>66.6°)? True
Sistemi che falliscono:  ['P', 'K', 'G']
Sistemi instabili:       []
Sistemi sempre stabili:  ['E', 'W', 'O', 'M', 'X', 'V', 'N']
```

### Cosa succede se provi Placidus al polo?

Se chiami `houses()` con Placidus a una latitudine polare, la libreria solleva un errore `PolarCircleError`:

```python
import libephemeris as ephem
from libephemeris.exceptions import PolarCircleError

jd = ephem.julday(2024, 6, 21, 12.0)  # Solstizio d'estate

try:
    cusps, ascmc = ephem.houses(jd, 69.6, 19.0, ord('P'))
except PolarCircleError as e:
    print(f"Errore: {e}")
    print(f"Latitudine: {e.latitude}°")
    print(f"Soglia polare: {e.threshold:.2f}°")
    print(f"Sistema: {e.house_system}")
```

```
Errore: swe_houses: Placidus house system cannot be calculated at latitude
  69.60°N (within Northern polar circle). Polar threshold for obliquity
  23.44° is ±66.56°.
Latitudine: 69.6°
Soglia polare: 66.56°
Sistema: P
```

### La soluzione: `houses_with_fallback`

Nella pratica, non vuoi che il tuo programma si blocchi con un errore. Vuoi un risultato ragionevole. Per questo esiste `houses_with_fallback`: prova il sistema che hai chiesto, e se fallisce, usa automaticamente un sistema alternativo (per default Porfirio, che funziona ovunque).

```python
import libephemeris as ephem

jd = ephem.julday(2024, 6, 21, 12.0)

# Tromsø — Placidus fallirà, il fallback userà Porfirio
cusps, ascmc, usato_fallback, avviso = ephem.swe_houses_with_fallback(
    jd, 69.6, 19.0,
    hsys=ord('P'),
    fallback_hsys=ord('O')   # Porfirio come alternativa
)

if usato_fallback:
    print(f"Attenzione: {avviso}")
    print("Usando Porfirio al posto di Placidus")
else:
    print("Placidus calcolato normalmente")

# Le cuspidi sono comunque disponibili
segni = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

for i in range(12):
    lon_c = cusps[i]
    s = segni[int(lon_c / 30)]
    g = lon_c % 30
    print(f"Casa {i+1:2d}: {g:5.1f}° {s}")
```

```
Attenzione: Placidus house system unavailable at latitude 69.60°
  (polar circle threshold: 66.56°). Using Porphyry as fallback.
Usando Porfirio al posto di Placidus
Casa  1:   9.6° Lib
Casa  2:  12.3° Sco
Casa  3:  15.0° Sgr
Casa  4:  17.7° Cap
Casa  5:  15.0° Aqr
Casa  6:  12.3° Psc
Casa  7:   9.6° Ari
Casa  8:  12.3° Tau
Casa  9:  15.0° Gem
Casa 10:  17.7° Cnc
Casa 11:  15.0° Leo
Casa 12:  12.3° Vir
```

La funzione restituisce quattro valori:

- **cusps** — le 12 cuspidi (calcolate con il sistema primario o con il fallback)
- **ascmc** — gli 8 angoli (Ascendente, MC, ARMC, Vertice, ecc.)
- **usato_fallback** — `True` se ha dovuto usare il sistema alternativo
- **avviso** — un messaggio che spiega cosa è successo, oppure `None` se tutto è andato bene

### Consiglio pratico

Se il tuo software deve funzionare per utenti di tutto il mondo, usa sempre `houses_with_fallback` al posto di `houses`. In questo modo:

- Per le latitudini normali (la stragrande maggioranza dei casi), ottieni esattamente il sistema che hai chiesto
- Per le latitudini polari, ottieni un risultato ragionevole con Porfirio, invece di un crash
- Il flag `usato_fallback` ti permette di informare l'utente che il risultato usa un sistema diverso

In alternativa, se lavori con l'astrologia vedica o ellenistica, usa direttamente il sistema a Segno intero (`W`) o Uguale (`E`) — funzionano ovunque senza bisogno di fallback.

---

## Riepilogo

In questo capitolo abbiamo esplorato le case astrologiche, dalla teoria alla pratica.

**Concetti chiave:**

- Le **case** dividono il cielo locale in 12 settori basati su data, ora e luogo — a differenza dei segni, che dipendono solo dalla data
- L'**Ascendente** (cuspide della 1a casa) è il grado dell'eclittica che sorge a est; il **Medio Cielo** (cuspide della 10a) è il grado che culmina al meridiano
- Esistono oltre 20 sistemi di case perché non c'è un modo univoco di dividere una sfera 3D in 12 settori; ogni sistema usa un criterio geometrico diverso
- Alle **latitudini polari** (oltre ~66.56°), Placidus, Koch e Gauquelin non funzionano; sistemi come Porfirio, Uguale e Segno intero funzionano ovunque
- La **posizione di un pianeta nelle case** dipende non solo dalla sua longitudine ma anche dalla sua latitudine eclittica

**Funzioni introdotte:**

- `houses(jd, lat, lon, hsys)` — calcola le 12 cuspidi e gli angoli per un dato sistema di case
- `house_pos(armc, lat, obliquity, hsys, lon, lat_body)` — determina in quale casa cade un pianeta (restituisce un valore decimale dove la parte intera è il numero della casa)
- `house_name(hsys)` — restituisce il nome leggibile di un sistema di case (es. `house_name(ord('P'))` → `"Placidus"`)
- `houses_with_fallback(jd, lat, lon, hsys, fallback_hsys)` — come `houses`, ma con fallback automatico per le latitudini polari
- `get_extreme_latitude_info(lat)` — restituisce un dizionario con informazioni su quali sistemi funzionano a quella latitudine
- `gauquelin_sector(jd, planet, method, geopos)` — calcola il settore Gauquelin (1–36) di un pianeta
