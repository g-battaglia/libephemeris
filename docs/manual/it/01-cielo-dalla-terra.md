# Capitolo 1 — Il cielo visto dalla Terra

## Cosa imparerai

In questo capitolo costruirai un modello mentale del cielo. Alla fine saprai cos'è la sfera celeste, perché il Sole traccia un percorso chiamato eclittica, come lo zodiaco divide quel percorso in 12 segni, come il tuo orizzonte locale determina cosa puoi vedere, e cosa sono l'Ascendente e il Medio Cielo.

Questi concetti sono il fondamento di tutto il resto del manuale.

---

## 1.1 La sfera celeste

Immagina di trovarti in un campo aperto, di notte, senza luci artificiali. Le stelle sembrano dipinte sulla superficie interna di un'enorme sfera che ti circonda da ogni lato. Questa è la **sfera celeste** — un modello geometrico che gli astronomi usano da millenni.

La sfera celeste non esiste fisicamente: le stelle sono a distanze diversissime tra loro. Ma per descrivere *dove* appare un oggetto nel cielo, non ci serve sapere quanto è lontano — ci basta sapere in quale direzione guardare. La sfera celeste è esattamente questo: una mappa delle direzioni.

La sfera ha due punti speciali:

- Il **polo nord celeste**: il punto nel cielo verso cui punta l'asse di rotazione della Terra, prolungato verso nord. La Stella Polare (Polaris) si trova a meno di 1° da questo punto — per questo sembra restare ferma mentre le altre stelle girano.
- Il **polo sud celeste**: il punto opposto, visibile solo dall'emisfero sud. Non ha una stella luminosa che lo segna.

A metà strada tra i poli, la sfera è divisa dall'**equatore celeste** — il prolungamento dell'equatore terrestre proiettato nel cielo. Se ti trovi esattamente sull'equatore terrestre, l'equatore celeste passa per il tuo zenit (il punto direttamente sopra la tua testa).

### Perché le stelle "girano"

Se osservi il cielo per qualche ora, vedrai le stelle muoversi lentamente da est a ovest, tracciando archi paralleli all'equatore celeste. Non sono le stelle a muoversi — è la Terra che ruota su se stessa. La sfera celeste sembra girare in senso opposto alla rotazione terrestre: una rotazione completa ogni 23 ore, 56 minuti e 4 secondi (il **giorno siderale** — lo incontreremo nel Capitolo 2).

### 🌍 Vita reale

Se punti una macchina fotografica verso la Stella Polare con un'esposizione lunga, vedrai le stelle tracciare cerchi concentrici attorno ad essa. Quei cerchi sono la prova visiva della rotazione terrestre.

---

## 1.2 L'eclittica: il cammino del Sole

Le stelle sembrano fisse sulla sfera celeste (si muovono, ma talmente lentamente che a occhio nudo servono secoli per notarlo — ne parleremo al Capitolo 8). Il Sole, invece, si muove: ogni giorno si sposta di circa 1° rispetto alle stelle circostanti.

Se potessimo vedere le stelle anche di giorno e segnassimo la posizione del Sole ogni giorno per un anno, il percorso formerebbe un grande cerchio sulla sfera celeste. Questo cerchio si chiama **eclittica**.

L'eclittica non coincide con l'equatore celeste: è inclinata rispetto ad esso di circa 23.4°. Questo angolo si chiama **obliquità dell'eclittica** ed è la ragione per cui esistono le stagioni. Quando il Sole si trova nella parte dell'eclittica che sta sopra l'equatore celeste, nell'emisfero nord è estate — il Sole sale più alto nel cielo e resta sopra l'orizzonte più a lungo.

L'eclittica e l'equatore celeste si intersecano in due punti:

- L'**equinozio di primavera** (o punto vernale, o punto gamma ♈): il Sole passa dall'emisfero celeste sud a quello nord. Questo è il punto zero dello zodiaco.
- L'**equinozio d'autunno**: il Sole passa dall'emisfero celeste nord a quello sud.

### 💻 Codice: ottenere l'obliquità dell'eclittica

Con LibEphemeris puoi ottenere il valore esatto dell'obliquità per qualsiasi data. La "pseudo-pianeta" `SE_ECL_NUT` restituisce le informazioni su obliquità e nutazione:

```python
import libephemeris as ephem

# Equinozio di primavera 2024 (20 marzo, ore 3:06 UT)
jd = ephem.julday(2024, 3, 20, 3.1)

# SE_ECL_NUT restituisce obliquità e nutazione
nut, flag = ephem.calc_ut(jd, ephem.ECL_NUT, 0)

obliquita_vera = nut[0]    # obliquità vera (con nutazione)
obliquita_media = nut[1]   # obliquità media (senza nutazione)
delta_psi = nut[2]         # nutazione in longitudine
delta_eps = nut[3]         # nutazione in obliquità

print(f"Obliquità vera:  {obliquita_vera:.6f}°")
print(f"Obliquità media: {obliquita_media:.6f}°")
```

```
Obliquità vera:  23.438703°
Obliquità media: 23.436129°
```

La differenza tra obliquità vera e media è la **nutazione in obliquità** (`delta_eps`): una piccola oscillazione dell'asse terrestre di cui parleremo al Capitolo 3.

### 🌍 Vita reale

D'estate il Sole sale più alto nel cielo perché si trova nella parte dell'eclittica sopra l'equatore celeste. Al solstizio d'estate, a Roma (lat. 41.9° N), il Sole raggiunge un'altezza massima di circa 71.5° — quasi allo zenit. Al solstizio d'inverno, appena 24.7°.

---

## 1.3 Lo zodiaco: 12 settori da 30°

L'eclittica è un cerchio di 360°. Per comodità, fin dall'antichità è stata divisa in **12 settori uguali da 30° ciascuno**, chiamati **segni zodiacali**:

| Segno | Intervallo | Segno | Intervallo |
|-------|-----------|-------|-----------|
| ♈ Ariete | 0° – 30° | ♎ Bilancia | 180° – 210° |
| ♉ Toro | 30° – 60° | ♏ Scorpione | 210° – 240° |
| ♊ Gemelli | 60° – 90° | ♐ Sagittario | 240° – 270° |
| ♋ Cancro | 90° – 120° | ♑ Capricorno | 270° – 300° |
| ♌ Leone | 120° – 150° | ♒ Acquario | 300° – 330° |
| ♍ Vergine | 150° – 180° | ♓ Pesci | 330° – 360° |

Il punto zero — 0° Ariete — coincide con l'equinozio di primavera (punto vernale). Quando un astronomo dice "il Sole è a 105° di longitudine eclittica", un astrologo dice "il Sole è a 15° Cancro" (perché 105 = 90 + 15, e il Cancro inizia a 90°).

### Segni e costellazioni: non sono la stessa cosa

Questo è un punto fondamentale che genera confusione perenne. I **segni** zodiacali sono settori geometrici di 30° ciascuno, definiti a partire dall'equinozio di primavera. Le **costellazioni** zodiacali sono gruppi di stelle con dimensioni diverse (la Vergine copre oltre 40° di eclittica, lo Scorpione meno di 20°) e confini convenzionali stabiliti dall'Unione Astronomica Internazionale nel 1930.

Circa 2000 anni fa segni e costellazioni coincidevano grossomodo. Oggi non più: l'equinozio di primavera si è spostato di circa 24° verso le stelle dei Pesci. Questo spostamento — circa 50" d'arco all'anno — si chiama **precessione degli equinozi**, ed è dovuto al fatto che l'asse terrestre oscilla lentamente come una trottola (un giro completo in circa 26000 anni). Ne parleremo in dettaglio nel Capitolo 11, dove vedremo lo zodiaco siderale usato nell'astrologia vedica.

### 🌍 Vita reale

Quando qualcuno dice "sono del segno dell'Ariete", intende che alla data della sua nascita il Sole si trovava tra 0° e 30° di longitudine eclittica tropicale. Questo non ha nulla a che fare con la costellazione dell'Ariete: per la precessione, quel tratto di cielo è oggi allineato con le stelle dei Pesci.

---

## 1.4 L'orizzonte locale

Fin qui abbiamo descritto il cielo "in generale" — come apparirebbe da qualsiasi punto della Terra. Ma il cielo che **tu** vedi dipende da **dove** ti trovi.

Tre punti definiscono il tuo cielo locale:

- **Zenit**: il punto direttamente sopra la tua testa, a 90° dall'orizzonte.
- **Nadir**: il punto direttamente sotto i tuoi piedi, opposto allo zenit.
- **Orizzonte**: il cerchio a 360° attorno a te dove cielo e terra si incontrano.

Per localizzare un oggetto nel cielo dal tuo punto di osservazione servono due coordinate:

- **Altezza** (o elevazione): l'angolo sopra l'orizzonte. 0° = sull'orizzonte. 90° = allo zenit. Valori negativi = sotto l'orizzonte.
- **Azimut**: la direzione sull'orizzonte. In astronomia la convenzione più diffusa parte da sud: 0° = Sud, 90° = Ovest, 180° = Nord, 270° = Est.

> **Attenzione**: la convenzione dell'azimut varia. In navigazione e nella vita quotidiana, 0° = Nord. In astronomia (e in LibEphemeris), 0° = Sud. Se passi dati tra sistemi diversi, verifica la convenzione.

### 💻 Codice: altezza e azimut di un corpo celeste

La funzione `azalt` converte coordinate eclittiche o equatoriali in coordinate orizzontali. Serve la posizione dell'osservatore (longitudine, latitudine, altitudine):

```python
import libephemeris as ephem

# 15 settembre 2024, ore 21:00 UT
jd = ephem.julday(2024, 9, 15, 21.0)

# Posizione di Giove
pos, flag = ephem.calc_ut(jd, ephem.JUPITER, 0)

# Posizione dell'osservatore: Roma
# (longitudine Est, latitudine Nord, altitudine in metri)
geopos = (12.4964, 41.9028, 50.0)

# Convertiamo da coordinate eclittiche a orizzontali
# SE_ECL2HOR = da eclittiche a orizzontali
# atpress = 1013.25 mbar (pressione standard)
# attemp = 15.0 °C (temperatura standard)
hor = ephem.azalt(jd, ephem.ECL2HOR, geopos, 1013.25, 15.0,
                  (pos[0], pos[1], pos[2]))

azimut = hor[0]            # da Sud, verso Ovest
altezza_vera = hor[1]      # senza rifrazione
altezza_apparente = hor[2] # con rifrazione atmosferica

# Convertiamo l'azimut dalla convenzione astronomica (S=0)
# alla convenzione navigazionale (N=0) per leggibilità
az_nav = (azimut + 180.0) % 360.0

# Direzione approssimativa
direzioni = ["N", "NE", "E", "SE", "S", "SO", "O", "NO"]
dir_idx = int((az_nav + 22.5) % 360 / 45)
dir_nome = direzioni[dir_idx]

if altezza_apparente > 0:
    print(f"Giove è visibile! Altezza: {altezza_apparente:.1f}°, "
          f"direzione: {dir_nome} (azimut {az_nav:.1f}° da N)")
else:
    print(f"Giove è sotto l'orizzonte ({altezza_apparente:.1f}°)")
```

```
Giove è sotto l'orizzonte (-3.2°)
```

La pressione e la temperatura servono per calcolare la **rifrazione atmosferica**: l'atmosfera curva la luce, facendo apparire gli oggetti vicini all'orizzonte un po' più alti di quanto siano realmente (circa 34' all'orizzonte). Se passi `atpress = 0`, la rifrazione viene ignorata.

---

## 1.5 I quattro angoli: Ascendente, MC, Discendente, IC

L'orizzonte e il meridiano del luogo intersecano l'eclittica in quattro punti che in astrologia sono fondamentali. Questi quattro punti cambiano continuamente perché la sfera celeste ruota — e ruota veloce: l'Ascendente si sposta di circa 1° ogni 4 minuti.

### L'Ascendente (ASC)

Il grado dell'eclittica che in quel momento sta sorgendo a est. È il punto di intersezione tra l'eclittica e l'orizzonte orientale. In un tema natale, l'Ascendente definisce la cuspide della prima casa.

L'Ascendente dipende dall'**ora esatta** e dal **luogo**. Due persone nate lo stesso giorno ma a 10 minuti di distanza possono avere Ascendenti diversi. Per questo, in astrologia, l'ora di nascita è così importante.

### Il Medio Cielo (MC)

Il grado dell'eclittica che in quel momento culmina — cioè attraversa il meridiano superiore (il semicerchio che va dal polo nord celeste allo zenit al polo sud celeste). Il MC è il punto più alto che quel grado dell'eclittica raggiungerà nel suo percorso attraverso il cielo. In un tema natale, il MC definisce la cuspide della decima casa.

### Discendente (DSC) e Imum Coeli (IC)

Il **Discendente** è opposto all'Ascendente: il grado dell'eclittica che sta tramontando a ovest. L'**Imum Coeli** (IC, "fondo del cielo") è opposto al MC: il grado dell'eclittica alla culminazione inferiore, sotto i piedi dell'osservatore. In un tema natale, DSC e IC definiscono rispettivamente la cuspide della settima e della quarta casa.

> **Nota**: MC e ASC **non** sono sempre a 90° di distanza tra loro. La distanza dipende dalla latitudine dell'osservatore e dal momento del giorno. Solo all'equatore, e solo in certi momenti, MC e ASC distano esattamente 90°.

### 💻 Codice: calcolare i quattro angoli

```python
import libephemeris as ephem

# 8 aprile 2024, ore 12:00 UT — eclissi solare
jd = ephem.julday(2024, 4, 8, 12.0)

# Roma: lat 41.9, lon 12.5
lat, lon = 41.9028, 12.4964

# ord('P') = Placidus, il sistema di case più diffuso
cusps, ascmc = ephem.houses(jd, lat, lon, ord('P'))

asc     = ascmc[0]   # Ascendente
mc      = ascmc[1]   # Medio Cielo
armc    = ascmc[2]   # Ascensione Retta del Medio Cielo (in gradi)
vertex  = ascmc[3]   # Vertice

segni = [
    "Ariete", "Toro", "Gemelli", "Cancro",
    "Leone", "Vergine", "Bilancia", "Scorpione",
    "Sagittario", "Capricorno", "Acquario", "Pesci",
]

def formatta_posizione(gradi):
    """Converte gradi decimali in formato zodiacale."""
    deg, min, sec, secfr, segno = ephem.split_deg(
        gradi, ephem.SPLIT_DEG_ZODIACAL | ephem.SPLIT_DEG_ROUND_SEC
    )
    return f"{deg}° {min}' {sec}\" {segni[segno]}"

print(f"Ascendente:  {formatta_posizione(asc)}")
print(f"Medio Cielo: {formatta_posizione(mc)}")
print(f"Discendente: {formatta_posizione((asc + 180) % 360)}")
print(f"Imum Coeli:  {formatta_posizione((mc + 180) % 360)}")
print(f"Vertex:      {formatta_posizione(vertex)}")
```

```
Ascendente:  13° 4' 44" Leone
Medio Cielo: 1° 54' 15" Toro
Discendente: 13° 4' 44" Acquario
Imum Coeli:  1° 54' 15" Scorpione
Vertex:      0° 46' 36" Capricorno
```

La funzione `houses` restituisce due tuple:

- `cusps`: le 12 cuspidi delle case (le vedremo nel Capitolo 7)
- `ascmc`: gli angoli principali — ASC, MC, ARMC, Vertex, più altri quattro angoli specialistici (Ascendente equatoriale, co-Ascendente Koch, co-Ascendente Munkasey, Ascendente polare)

Il parametro `ord('P')` indica il sistema di case Placidus, il più usato in Occidente. Esistono oltre 20 sistemi diversi — li esploreremo nel Capitolo 7.

### 🌍 Vita reale

Quando un astrologo chiede "qual è il tuo Ascendente?", sta chiedendo quale grado dell'eclittica stava sorgendo nel momento esatto della tua nascita, visto dal luogo dove sei nato. L'Ascendente cambia segno ogni 2 ore circa (ma non in modo regolare — alcuni segni sorgono più rapidamente di altri, a seconda della latitudine). Per questo servono data, ora e luogo precisi.

A latitudini molto alte (sopra il circolo polare), in certi momenti dell'anno l'eclittica non interseca affatto l'orizzonte, e l'Ascendente non è calcolabile con alcuni sistemi di case. LibEphemeris gestisce questa situazione con fallback automatici — ne parleremo nel Capitolo 7.

---

## Riepilogo

- La **sfera celeste** è un modello geometrico: una sfera immaginaria su cui proiettiamo le posizioni degli oggetti celesti.
- L'**eclittica** è il percorso annuale del Sole sulla sfera celeste, inclinata di ~23.4° rispetto all'equatore celeste.
- Lo **zodiaco** divide l'eclittica in 12 segni da 30° ciascuno, partendo dall'equinozio di primavera. I segni sono settori geometrici, non costellazioni.
- L'**orizzonte locale** dipende dalla tua posizione. Altezza e azimut descrivono dove appare un oggetto nel tuo cielo.
- L'**Ascendente** è il grado dell'eclittica che sorge a est; il **Medio Cielo** è il grado che culmina. Cambiano rapidamente e dipendono da ora e luogo.

### Funzioni e costanti introdotte

| Funzione / Costante | Uso |
|---------------------|-----|
| `calc_ut(jd, SE_ECL_NUT, 0)` | Obliquità dell'eclittica e nutazione |
| `azalt(jd, calc_flag, geopos, atpress, attemp, xin)` | Coordinate eclittiche/equatoriali → orizzontali |
| `houses(jd, lat, lon, ord('P'))` | Cuspidi delle case e angoli (ASC, MC, ...) |
| `split_deg(gradi, flags)` | Formattazione zodiacale |
| `SE_ECL2HOR`, `SE_EQU2HOR` | Flag per `azalt`: input eclittico o equatoriale |
