# Capitolo 4 — Le efemeridi: cosa sono e come funzionano

## Cosa imparerai

In questo capitolo scoprirai cosa c'è "dentro" un'efemeride, da dove vengono i dati che la libreria usa, i tre livelli di precisione disponibili, e la differenza tra posizioni geocentriche, eliocentriche, baricentriche e topocentriche.

---

## 4.1 Cos'è un'efemeride

Un'**efemeride** è una tabella che dice "a questa data, questo corpo celeste si trova qui". Il concetto è antico quanto la civiltà: le tavolette babilonesi di 3000 anni fa contenevano efemeridi lunari per prevedere le eclissi. L'Almanacco Nautico, pubblicato dal 1767, fornisce posizioni stellari ai marinai per la navigazione.

Un'efemeride cartacea tipica ha una riga per ogni giorno e una colonna per ogni pianeta. Un'efemeride digitale è la stessa idea, ma con precisione molto maggiore e la capacità di interpolare tra i punti tabulati per ottenere la posizione a qualsiasi istante.

Le efemeridi moderne non sono semplici tabelle: sono file binari che contengono **polinomi di Chebyshev** — funzioni matematiche che approssimano la traiettoria di ogni corpo celeste con errori sub-millimetrici. Dato un istante qualsiasi, la libreria valuta questi polinomi e ottiene la posizione esatta.

---

## 4.2 JPL DE440: la nostra sorgente dati

LibEphemeris usa le efemeridi prodotte dal **Jet Propulsion Laboratory** (JPL) della NASA — lo stesso laboratorio che guida le sonde spaziali. Le efemeridi JPL sono il gold standard: sono prodotte integrando numericamente le equazioni del moto di tutti i corpi del Sistema Solare, tenendo conto di centinaia di effetti gravitazionali e relativistici.

La libreria supporta tre file di efemeridi:

| File | Nome | Intervallo | Dimensione |
|------|------|-----------|-----------|
| `de440s.bsp` | DE440 short | 1849 – 2150 | ~31 MB |
| `de440.bsp` | DE440 | 1550 – 2650 | ~114 MB |
| `de441.bsp` | DE441 | -13200 – +17191 | ~3.1 GB |

**DE440** (2020) è l'efemeride di riferimento: copre 1100 anni con precisione sub-milliarcsecondo per i pianeti interni. **DE441** è la versione estesa a 30000 anni, utile per ricerca storica ma con precisione gradualmente inferiore per date molto lontane.

Internamente, LibEphemeris usa **Skyfield** per leggere i file JPL e calcolare le posizioni. Skyfield è una libreria astronomica Python di alta qualità, sviluppata da Brandon Rhodes.

---

## 4.3 I tre livelli di precisione

La scelta del livello di precisione determina quale file di efemeridi viene usato:

- **`base`** (DE440s): leggero (~31 MB), copre 1849–2150. Ideale per astrologia moderna e applicazioni mobile.
- **`medium`** (DE440): bilanciato (~114 MB), copre 1550–2650. Il **default**. Va bene per quasi tutti gli usi.
- **`extended`** (DE441): massimo (~3.1 GB), copre -13200 a +17191. Per ricerca storica e date antiche/future.

```python
import libephemeris as ephem

# Vedere il livello attuale
# Il default è "medium"

# Cambiare livello
ephem.set_precision_tier("extended")  # per date antiche

# Scaricare i file necessari (una tantum)
ephem.download_for_tier("medium")

# Tornare al default
ephem.set_precision_tier("medium")
```

### Quale scegliere?

Per **astrologia moderna** (temi natali dal 1900 in poi): `base` è sufficiente. Per **astrologia storica** (prima del 1849): serve almeno `medium`. Per **ricerca storica antica** (eclissi babilonesi, comete medievali): serve `extended`.

---

## 4.4 Geocentrico, eliocentrico e baricentrico

Le posizioni calcolate dipendono da **dove** ti immagini di osservare:

- **Geocentrico** (default): visto dal centro della Terra. È il punto di vista dell'astrologia — il cielo come appare dalla Terra.
- **Eliocentrico**: visto dal centro del Sole. Flag `FLG_HELCTR`. Utile per meccanica celeste.
- **Baricentrico**: visto dal baricentro (centro di massa) del Sistema Solare. Flag `FLG_BARYCTR`. Il Sole non sta esattamente al baricentro: "oscilla" di circa 2 raggi solari a causa dell'attrazione dei pianeti giganti.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Marte geocentrico (default)
geo, _ = ephem.calc_ut(jd, ephem.MARS, 0)

# Marte eliocentrico
helio, _ = ephem.calc_ut(jd, ephem.MARS, ephem.FLG_HELCTR)

print(f"Marte geocentrico:  {geo[0]:.4f}° (dist {geo[2]:.4f} UA)")
print(f"Marte eliocentrico: {helio[0]:.4f}° (dist {helio[2]:.4f} UA)")
```

```
Marte geocentrico:  342.8458° (dist 2.0618 UA)
Marte eliocentrico: 317.5543° (dist 1.3879 UA)
```

### Planetocentrico: osservare da un altro pianeta

Con `swe_calc_pctr` puoi calcolare la posizione di un corpo visto da un altro pianeta:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Il Sole visto da Giove
pos, _ = ephem.swe_calc_pctr(jd, ephem.SUN, ephem.JUPITER, 0)
print(f"Sole visto da Giove: {pos[0]:.4f}°, dist {pos[2]:.4f} UA")
```

```
Sole visto da Giove: 234.7096°, dist 5.0060 UA
```

---

## 4.5 Topocentrico: la posizione conta

La posizione **geocentrica** è calcolata dal centro della Terra. Ma tu non sei al centro della Terra — sei sulla sua superficie. La posizione **topocentrica** tiene conto della tua posizione esatta.

Per la maggior parte dei corpi celesti la differenza è trascurabile (meno di 1"). Ma per la **Luna**, che è vicinissima, la differenza può raggiungere ~1° — la cosiddetta **parallasse lunare**. Questo significa che un'eclissi solare è totale in un luogo e parziale in un altro a poche centinaia di km di distanza.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 18.0)

# Imposta la posizione dell'osservatore (Roma)
ephem.set_topo(12.4964, 41.9028, 50.0)

# Luna geocentrica (default)
luna_geo, _ = ephem.calc_ut(jd, ephem.MOON, 0)

# Luna topocentrica
luna_topo, _ = ephem.calc_ut(jd, ephem.MOON, ephem.FLG_TOPOCTR)

diff = luna_topo[0] - luna_geo[0]
print(f"Luna geocentrica:  {luna_geo[0]:.4f}°")
print(f"Luna topocentrica: {luna_topo[0]:.4f}°")
print(f"Parallasse:        {diff * 3600:.1f}\" d'arco")
```

```
Luna geocentrica:  19.1832°
Luna topocentrica: 18.2383°
Parallasse:        -3401.7" d'arco
```

### 🌍 Vita reale

L'eclissi solare dell'8 aprile 2024 era totale a Dallas (Texas) ma solo parziale a New York — la differenza è interamente dovuta alla parallasse: la Luna copre il Sole in modo diverso a seconda di dove ti trovi sulla superficie terrestre.

---

## Riepilogo

- Un'**efemeride** è una tabella di posizioni celesti. LibEphemeris usa le efemeridi NASA JPL DE440/DE441, lette tramite Skyfield.
- Tre livelli: `base` (1849–2150), `medium` (1550–2650, default), `extended` (-13200 a +17191).
- **Geocentrico** (default): dal centro della Terra. **Eliocentrico** (`FLG_HELCTR`): dal Sole. **Baricentrico** (`FLG_BARYCTR`): dal centro di massa del Sistema Solare.
- **Topocentrico** (`FLG_TOPOCTR`): dalla tua posizione sulla superficie terrestre. Cruciale per la Luna (~1° di parallasse).

### Funzioni e costanti introdotte

| Funzione / Costante | Uso |
|---------------------|-----|
| `set_precision_tier(tier)` | Sceglie `"base"`, `"medium"` o `"extended"` |
| `download_for_tier(tier)` | Scarica i file di efemeridi |
| `set_topo(lon, lat, alt)` | Imposta la posizione dell'osservatore |
| `swe_calc_pctr(jd, body, center, flag)` | Posizione vista da un altro pianeta |
| `FLG_HELCTR` | Coordinate eliocentriche |
| `FLG_BARYCTR` | Coordinate baricentriche |
| `FLG_TOPOCTR` | Coordinate topocentriche |
