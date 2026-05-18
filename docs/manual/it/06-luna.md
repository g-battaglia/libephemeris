# Capitolo 6 — La Luna: nodi, apogeo, perigeo e Lilith

## Cosa imparerai

In questo capitolo scoprirai perché la Luna è il corpo celeste più complesso da calcolare, cosa sono i nodi lunari e perché sono legati alle eclissi, come funzionano apogeo e perigeo (e cos'è la "superluna"), e cosa si intende per Lilith Nera in astrologia.

---

## 6.1 Perché la Luna è speciale

La Luna è il corpo celeste più vicino a noi e anche quello con il moto più complesso. Ecco perché merita un capitolo intero.

**È velocissima.** Mentre il Sole si sposta di circa 1° al giorno lungo lo zodiaco, la Luna percorre circa 12°–15° al giorno — attraversa un intero segno zodiacale in poco più di due giorni. Questo significa che l'ora di nascita conta moltissimo per la posizione della Luna: in sole 2 ore, la Luna si sposta di più di 1°.

**La sua orbita è instabile.** L'orbita lunare è un'ellisse, ma un'ellisse che cambia continuamente forma, inclinazione e orientamento. Il Sole "tira" la Luna alterandone l'orbita in modi complessi. I matematici hanno catalogato centinaia di perturbazioni lunari — è il problema a tre corpi più studiato della storia dell'astronomia.

**La parallasse conta.** Per tutti gli altri corpi celesti, la differenza tra posizione geocentrica (dal centro della Terra) e topocentrica (dalla tua posizione sulla superficie) è trascurabile — meno di 1 secondo d'arco. Per la Luna, questa differenza può raggiungere quasi 1°. È per questo che un'eclissi solare è totale in un luogo e parziale in un altro a poche centinaia di km di distanza.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 18.0)

# Posizione della Luna — notare la velocità
pos, _ = ephem.calc_ut(jd, ephem.MOON, 0)

print(f"Luna: {pos[0]:.4f}° di longitudine eclittica")
print(f"Latitudine: {pos[1]:+.4f}°")
print(f"Distanza: {pos[2]:.6f} UA ({pos[2] * 149597870.7:.0f} km)")
print(f"Velocità: {pos[3]:.4f}°/giorno")
```

```
Luna: 19.1832° di longitudine eclittica
Latitudine: +0.3292°
Distanza: 0.002405 UA (359781 km)
Velocità: 14.9963°/giorno
```

---

## 6.2 I nodi lunari: dove la Luna attraversa l'eclittica

L'orbita della Luna non giace sullo stesso piano dell'eclittica (il piano dell'orbita terrestre attorno al Sole). È inclinata di circa 5.1° rispetto ad essa. I due punti dove l'orbita lunare interseca l'eclittica si chiamano **nodi**.

Immagina l'eclittica come il pavimento di una stanza. L'orbita della Luna è un cerchio inclinato che attraversa il pavimento in due punti:

- Il **nodo ascendente** (o nodo nord) è il punto dove la Luna "sale" — passa dal sotto al sopra dell'eclittica, dalla latitudine negativa a quella positiva.
- Il **nodo discendente** (o nodo sud) è il punto opposto — la Luna "scende" sotto l'eclittica. È sempre a 180° dal nodo ascendente.

I nodi non sono fissi: si spostano lentamente all'indietro lungo lo zodiaco, completando un giro in circa **18.6 anni**. Questo moto retrogrado dei nodi è causato dall'attrazione gravitazionale del Sole sull'orbita lunare.

### Nodo medio e nodo vero

LibEphemeris offre due versioni dei nodi lunari:

- Il **nodo medio** (`MEAN_NODE`) è calcolato con un polinomio matematico regolare. Si muove in modo uniforme e prevedibile — nessuna oscillazione, nessun sussulto. È come l'orologio perfetto del moto nodale.

- Il **nodo vero** (`TRUE_NODE`) tiene conto di tutte le perturbazioni reali. Il nodo vero oscilla avanti e indietro attorno alla posizione media con un'ampiezza di circa ±1.5°. Nelle efemeridi dettagliate, il nodo vero può muoversi brevemente in moto diretto (in avanti), cosa che il nodo medio non fa mai.

In astrologia, la scelta tra nodo medio e vero è oggetto di dibattito. Molti astrologi usano il nodo medio per la sua regolarità; altri preferiscono il nodo vero perché riflette la realtà astronomica.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Nodo medio
pos_medio, _ = ephem.calc_ut(jd, ephem.MEAN_NODE, 0)

# Nodo vero
pos_vero, _ = ephem.calc_ut(jd, ephem.TRUE_NODE, 0)

diff = pos_vero[0] - pos_medio[0]
print(f"Nodo medio: {pos_medio[0]:.4f}°")
print(f"Nodo vero:  {pos_vero[0]:.4f}°")
print(f"Differenza: {diff:.4f}° ({diff * 60:.1f}' d'arco)")

# Funzioni dedicate (accettano JD in TT)
jd_tt = jd + ephem.deltat(jd)
lon_medio = ephem.calc_mean_lunar_node(jd_tt)
lon_vero, lat_vero, dist_vero = ephem.calc_true_lunar_node(jd_tt)

print(f"\nNodo medio (funzione dedicata): {lon_medio:.4f}°")
print(f"Nodo vero  (funzione dedicata): {lon_vero:.4f}°")
```

```
Nodo medio: 15.6610°
Nodo vero:  15.6269°
Differenza: -0.0340° (-2.0' d'arco)

Nodo medio (funzione dedicata): 15.6624°
Nodo vero  (funzione dedicata): 15.6269°
```

### Perché i nodi contano: le eclissi

Le eclissi avvengono **solo** quando il Sole si trova vicino a un nodo lunare. Perché? Un'eclissi solare richiede che la Luna passi davanti al Sole — ma se l'orbita lunare è inclinata di 5.1°, la Luna di solito passa sopra o sotto il Sole, mancandolo. Solo quando la Luna nuova avviene vicino a un nodo, la Luna è abbastanza allineata con l'eclittica da coprire il Sole (ne parleremo in dettaglio nel Capitolo 9).

---

## 6.3 Apogeo e perigeo: la Luna vicina e lontana

L'orbita della Luna non è un cerchio perfetto ma un'**ellisse** — una forma ovale. Questo significa che la distanza tra la Terra e la Luna cambia continuamente:

- Al **perigeo** (punto più vicino), la Luna dista circa 356000 km dalla Terra.
- All'**apogeo** (punto più lontano), la distanza sale a circa 407000 km.

La differenza è notevole: al perigeo la Luna appare circa il 14% più grande e il 30% più luminosa che all'apogeo. Quando una Luna piena coincide con il perigeo, i media la chiamano **"superluna"** — anche se la differenza a occhio nudo è difficile da percepire.

Ma c'è una complicazione: a causa delle perturbazioni del Sole, la posizione dell'apogeo e del perigeo **cambia continuamente**. L'ellisse dell'orbita lunare ruota lentamente nello spazio, completando un giro in circa 8.85 anni. Inoltre, l'ellisse stessa cambia forma di mese in mese.

Per questo LibEphemeris offre tre versioni diverse dell'apogeo (e analogamente del perigeo):

**Apogeo medio** — la posizione calcolata con un polinomio regolare, che avanza uniformemente lungo lo zodiaco. È il più usato in astrologia con il nome di "Lilith Nera" (ne parleremo nella sezione successiva).

**Apogeo osculante** — la posizione istantanea dell'apogeo, calcolata dall'orbita reale della Luna in quel preciso momento. Include tutte le perturbazioni e può oscillare di 20°–30° rispetto alla posizione media. "Osculante" viene dal latino *osculari* (baciare) — l'orbita osculante è l'ellisse che "bacia" la traiettoria reale in un dato istante.

**Apogeo interpolato** — una via di mezzo: le oscillazioni troppo rapide (artifici matematici dell'orbita osculante) vengono rimosse, mantenendo le variazioni fisicamente significative. È la versione più precisa dal punto di vista astronomico.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Apogeo medio (via calc_ut)
medio, _ = ephem.calc_ut(jd, ephem.MEAN_APOG, 0)

# Apogeo osculante (via calc_ut)
oscu, _ = ephem.calc_ut(jd, ephem.OSCU_APOG, 0)

# Apogeo interpolato (funzione dedicata, vuole JD in TT)
jd_tt = jd + ephem.deltat(jd)
lon_int, lat_int, dist_int = ephem.calc_interpolated_apogee(jd_tt)

print(f"Apogeo medio:       {medio[0]:.4f}°")
print(f"Apogeo osculante:   {oscu[0]:.4f}°")
print(f"Apogeo interpolato: {lon_int:.4f}°")
print(f"Differenza medio-osculante: {oscu[0] - medio[0]:.1f}°")
```

```
Apogeo medio:       170.9201°
Apogeo osculante:   182.7118°
Apogeo interpolato: 166.3793°
Differenza medio-osculante: 11.8°
```

---

## 6.4 Lilith: la Luna Nera

In astrologia, la **Lilith Nera** (o semplicemente "Lilith") è l'apogeo dell'orbita lunare — il punto dove la Luna è più lontana dalla Terra. Non è un corpo fisico: è un punto geometrico nello spazio, il fuoco vuoto dell'ellisse orbitale della Luna.

Il nome viene dalla mitologia: Lilith è la prima moglie di Adamo secondo alcune tradizioni ebraiche, associata al lato oscuro, selvaggio e indomabile della natura femminile. In astrologia, rappresenta gli istinti profondi, ciò che è nascosto o represso.

Esistono tre versioni di Lilith, che corrispondono alle tre versioni dell'apogeo descritte sopra:

**Lilith media** è la più usata in astrologia. Corrisponde all'apogeo medio e si muove in modo regolare, completando un giro dello zodiaco in circa 8.85 anni (poco meno di 9 anni). È la versione che trovi nella maggior parte dei software astrologici e nelle effemeridi stampate.

**Lilith vera (osculante)** è la posizione reale istantanea dell'apogeo. Può differire dalla Lilith media di 20°–30° e ha movimenti irregolari, incluse brevi retrogradazioni. Alcuni astrologi la preferiscono per la sua aderenza alla realtà astronomica.

**Lilith interpolata** è la versione lisciata — rimuove le oscillazioni artificiali dell'orbita osculante ma mantiene le variazioni fisicamente reali. È la più precisa dal punto di vista astronomico.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)
jd_tt = jd + ephem.deltat(jd)

# Lilith media (= apogeo medio)
lilith_media = ephem.calc_mean_lilith(jd_tt)

# Lilith vera/osculante
lon_vera, lat_vera, dist_vera = ephem.calc_true_lilith(jd_tt)

segni = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

def segno(lon):
    return f"{lon % 30:.1f}° {segni[int(lon / 30)]}"

print(f"Lilith media: {segno(lilith_media)}")
print(f"Lilith vera:  {segno(lon_vera)}")
print(f"Differenza:   {lon_vera - lilith_media:.1f}°")
```

```
Lilith media: 20.9° Vir
Lilith vera:  2.7° Lib
Differenza:   11.8°
```

### La Luna Bianca (Selena)

Il punto opposto a Lilith — cioè il **perigeo** medio dell'orbita lunare — è chiamato **Luna Bianca** o **Selena** in alcune tradizioni astrologiche. Se Lilith rappresenta il lato ombra, Selena rappresenta il lato luminoso.

```python
import libephemeris as ephem

jd_tt = ephem.julday(2024, 4, 8, 12.0) + ephem.deltat(
    ephem.julday(2024, 4, 8, 12.0)
)

# Selena: il punto opposto a Lilith
selena = ephem.calc_white_moon_position(jd_tt)
# Restituisce una tupla di 6 valori come calc_ut:
# (lon, lat, dist, vel_lon, vel_lat, vel_dist)

lilith = ephem.calc_mean_lilith(jd_tt)
print(f"Lilith:  {lilith:.4f}°")
print(f"Selena:  {selena[0]:.4f}°")
print(f"Differenza: {abs(selena[0] - lilith):.1f}°")  # ~180°
```

```
Lilith:  170.9216°
Selena:  350.9216°
Differenza: 180.0°
```

---

## 6.5 Perigeo interpolato e calibrazione

Così come esiste l'apogeo interpolato, esiste anche il **perigeo interpolato** — la versione lisciata del punto di massimo avvicinamento della Luna. È utile per calcoli di precisione, ad esempio per prevedere le maree o per determinare le "superlune" vere (Luna piena entro poche ore dal perigeo).

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)
jd_tt = jd + ephem.deltat(jd)

# Perigeo interpolato
lon, lat, dist = ephem.calc_interpolated_perigee(jd_tt)

# Distanza in km (1 UA = 149597870.7 km)
dist_km = dist * 149597870.7

print(f"Perigeo interpolato: {lon:.4f}°")
print(f"Distanza: {dist_km:.0f} km")
```

```
Perigeo interpolato: 4.1030°
Distanza: 358786 km
```

La precisione del perigeo interpolato in LibEphemeris è stata migliorata attraverso un processo di calibrazione contro le effemeridi JPL ad alta precisione. I dettagli tecnici di questo processo sono descritti nella documentazione di sviluppo del progetto.

---

## 6.6 Attraversamenti del nodo

In certi casi serve sapere **quando esattamente** la Luna attraversa il piano dell'eclittica — cioè quando la sua latitudine eclittica passa per zero. Questo è importante per il calcolo delle eclissi e per alcune tecniche astrologiche.

La funzione `mooncross_node_ut` cerca, a partire da una data iniziale, il prossimo momento in cui la Luna attraversa l'eclittica:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 1, 0.0)

# Trova il prossimo attraversamento del nodo
jd_cross, lon_cross, lat_cross = ephem.mooncross_node_ut(jd)

anno, mese, giorno, ore = ephem.revjul(jd_cross)
h = int(ore)
m = int((ore - h) * 60)

print(f"Prossimo attraversamento: {giorno}/{mese}/{anno} {h}:{m:02d} UT")
print(f"Longitudine al nodo: {lon_cross:.4f}°")

# lat_cross è ~0 per definizione (è il momento
# in cui la Luna è esattamente sull'eclittica)
```

```
Prossimo attraversamento: 8/4/2024 12:18 UT
Longitudine al nodo: 15.6269°
```

La Luna attraversa un nodo circa **due volte al mese** — una volta il nodo ascendente (latitudine da negativa a positiva) e una volta il nodo discendente (da positiva a negativa). Puoi distinguere i due casi dal segno della velocità in latitudine al momento dell'attraversamento: positiva = nodo ascendente, negativa = nodo discendente.

---

## Riepilogo

- La Luna è il corpo celeste più complesso da calcolare: è velocissima (~13°/giorno), la sua orbita è continuamente perturbata dal Sole, e la parallasse può raggiungere quasi 1°.
- I **nodi lunari** sono i punti dove l'orbita lunare interseca l'eclittica. Esistono in versione media (regolare) e vera (con oscillazioni). Le eclissi avvengono solo vicino ai nodi.
- **Apogeo** (punto più lontano) e **perigeo** (punto più vicino) esistono in tre versioni: medio, osculante e interpolato. L'apogeo medio è la "Lilith Nera" dell'astrologia.
- **Lilith** (Luna Nera) è l'apogeo dell'orbita lunare — un punto geometrico, non un corpo fisico. La Lilith media è la più usata; la vera può differire di 20°–30°.
- La **Luna Bianca** (Selena) è il punto opposto a Lilith, corrispondente al perigeo.
- `mooncross_node_ut` trova il momento esatto in cui la Luna attraversa l'eclittica.

### Funzioni introdotte

- `calc_mean_lunar_node(jd_tt)` — longitudine del nodo lunare medio
- `calc_true_lunar_node(jd_tt)` — longitudine, latitudine e distanza del nodo lunare vero
- `calc_mean_lilith(jd_tt)` — longitudine della Lilith media (apogeo medio)
- `calc_true_lilith(jd_tt)` — longitudine, latitudine e distanza della Lilith vera (apogeo osculante)
- `calc_interpolated_apogee(jd_tt)` — apogeo lunare interpolato (lisciato)
- `calc_interpolated_perigee(jd_tt)` — perigeo lunare interpolato
- `calc_white_moon_position(jd_tt)` — Luna Bianca (Selena), opposta a Lilith
- `mooncross_node_ut(jd_ut)` — prossimo attraversamento del nodo lunare
- `MEAN_NODE`, `TRUE_NODE` — nodi lunari via `calc_ut`
- `MEAN_APOG`, `OSCU_APOG` — apogeo lunare via `calc_ut`
