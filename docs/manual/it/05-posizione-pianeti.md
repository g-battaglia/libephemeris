# Capitolo 5 — Calcolare la posizione dei pianeti

## Cosa imparerai

In questo capitolo padroneggerai `calc_ut` — la funzione più usata della libreria. Imparerai cosa restituisce, tutti i flag disponibili, come funzionano velocità e moto retrogrado, e come calcolare fenomeni planetari come magnitudine, fase e elongazione.

---

## 5.1 La funzione principale: `calc_ut`

La firma è semplice:

```python
pos, flag = ephem.calc_ut(jd_ut, body, iflag)
```

- `jd_ut`: Giorno Giuliano in UT (il tempo "civile")
- `body`: identificatore del corpo celeste (`SE_SUN`, `SE_MOON`, ecc.)
- `iflag`: flag di calcolo (combinati con `|`)

Il risultato è composto da due parti: una tupla `pos` con 6 numeri e un intero `flag` che conferma le opzioni di calcolo usate. I 6 numeri descrivono **dove** si trova il corpo celeste e **come si sta muovendo**.

### I primi tre: dove si trova

Per capire questi valori, immagina lo zodiaco come una lunga strada circolare di 360 km. Ogni pianeta ha una posizione lungo questa strada, una distanza laterale dalla strada, e una distanza da te.

- **`pos[0]` — Longitudine eclittica** è la posizione del pianeta lungo questa "strada zodiacale". Va da 0° a 360° e corrisponde ai 12 segni: i primi 30° sono l'Ariete, da 30° a 60° il Toro, e così via. Se `pos[0]` vale 105°, significa che il pianeta si trova a 15° del Cancro (perché il Cancro inizia a 90°, e 105 - 90 = 15). È il dato più importante in astrologia: quando qualcuno dice "il mio Sole è in Leone", intende che la longitudine eclittica del Sole al momento della nascita era tra 120° e 150°.

- **`pos[1]` — Latitudine eclittica** è la distanza del pianeta dal piano dell'eclittica, misurata in gradi. Pensa all'eclittica come un pavimento: la latitudine dice quanto il pianeta si trova sopra (+) o sotto (-) quel pavimento. Il Sole è sempre a ~0° (è lui stesso che definisce il pavimento). La Luna può salire fino a ±5.1°. La maggior parte dei pianeti resta entro pochi gradi dall'eclittica, ma Plutone può arrivare a ~17°.

- **`pos[2]` — Distanza** è quanto è lontano il corpo celeste dalla Terra, misurato in **unità astronomiche** (1 UA = ~150 milioni di km, la distanza media Terra-Sole). Per dare un'idea: la Luna è vicinissima (~0.0026 UA, circa 384000 km), Marte varia moltissimo tra 0.37 UA quando è dalla nostra parte dell'orbita e 2.68 UA quando è dall'altra parte del Sole, Nettuno è a circa 30 UA.

### Gli ultimi tre: come si sta muovendo

- **`pos[3]` — Velocità in longitudine** dice di quanti gradi il pianeta avanza (o arretra) lungo lo zodiaco ogni giorno. Il Sole si muove di ~1°/giorno — quasi impercettibile a occhio nudo. La Luna corre a ~13°/giorno, attraversando un intero segno in poco più di due giorni. Saturno si trascina a ~0.03°/giorno. Se questo valore è **negativo**, il pianeta è in **moto retrogrado**: sembra muoversi all'indietro lungo lo zodiaco (ne parleremo nella sezione 5.3).

- **`pos[4]` — Velocità in latitudine** dice quanto il pianeta si avvicina o allontana dall'eclittica ogni giorno. Di solito è un valore piccolo — la latitudine cambia lentamente.

- **`pos[5]` — Velocità in distanza** dice se il pianeta si sta avvicinando (valore negativo) o allontanando (positivo) dalla Terra, in UA/giorno.

### Come cambiano i valori con i flag

I valori descritti sopra sono quelli che ottieni con il flag di default (`0`). Ma se passi certi flag, il **significato** dei 6 numeri cambia completamente — la struttura della tupla resta la stessa, ma i numeri dentro rappresentano cose diverse:

- Con **`FLG_EQUATORIAL`**, `pos[0]` non è più la longitudine eclittica ma l'**Ascensione Retta** (la posizione lungo l'equatore celeste, da 0° a 360°), e `pos[1]` diventa la **Declinazione** (la distanza dall'equatore celeste, da -90° a +90°). Questo è il sistema usato dai telescopi. La distanza e le velocità mantengono lo stesso significato, ma riferite al sistema equatoriale.

- Con **`FLG_XYZ`**, tutti e 6 i valori diventano **coordinate cartesiane**: `pos[0]`, `pos[1]`, `pos[2]` sono le posizioni X, Y, Z in unità astronomiche, e `pos[3]`, `pos[4]`, `pos[5]` sono le velocità corrispondenti in UA/giorno. Non ci sono più gradi — solo distanze e velocità nello spazio tridimensionale.

- Con **`FLG_RADIANS`**, i valori angolari (longitudine, latitudine, velocità) sono espressi in **radianti** invece che in gradi (2π radianti = 360°). Utile se devi fare calcoli trigonometrici senza convertire.

Se non usi nessuno di questi flag (o passi `0`), ottieni sempre coordinate eclittiche in gradi — il sistema dell'astrologia.

### I corpi celesti principali

| Costante | Valore | Corpo |
|----------|--------|-------|
| `SE_SUN` | 0 | Sole |
| `SE_MOON` | 1 | Luna |
| `SE_MERCURY` | 2 | Mercurio |
| `SE_VENUS` | 3 | Venere |
| `SE_MARS` | 4 | Marte |
| `SE_JUPITER` | 5 | Giove |
| `SE_SATURN` | 6 | Saturno |
| `SE_URANUS` | 7 | Urano |
| `SE_NEPTUNE` | 8 | Nettuno |
| `SE_PLUTO` | 9 | Plutone |
| `SE_MEAN_NODE` | 10 | Nodo lunare medio |
| `SE_TRUE_NODE` | 11 | Nodo lunare vero |
| `SE_MEAN_APOG` | 12 | Apogeo lunare medio (Lilith) |
| `SE_OSCU_APOG` | 13 | Apogeo lunare osculante |
| `SE_CHIRON` | 15 | Chirone |

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Tutti i pianeti in un colpo
corpi = [
    (ephem.SUN, "Sole"), (ephem.MOON, "Luna"),
    (ephem.MERCURY, "Mercurio"), (ephem.VENUS, "Venere"),
    (ephem.MARS, "Marte"), (ephem.JUPITER, "Giove"),
    (ephem.SATURN, "Saturno"), (ephem.URANUS, "Urano"),
    (ephem.NEPTUNE, "Nettuno"), (ephem.PLUTO, "Plutone"),
]

segni = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

for body_id, nome in corpi:
    pos, _ = ephem.calc_ut(jd, body_id, ephem.FLG_SPEED)
    lon = pos[0]
    segno = segni[int(lon / 30)]
    gradi = lon % 30
    print(f"{nome:10s} {gradi:5.1f}° {segno}  (vel: {pos[3]:+.4f}°/g)")
```

```
Sole        19.1° Ari  (vel: +0.9831°/g)
Luna        15.4° Ari  (vel: +15.0292°/g)
Mercurio    25.0° Ari  (vel: -0.6229°/g)
Venere       4.1° Ari  (vel: +1.2353°/g)
Marte       12.8° Psc  (vel: +0.7775°/g)
Giove       19.0° Tau  (vel: +0.2211°/g)
Saturno     14.4° Psc  (vel: +0.1078°/g)
Urano       21.2° Tau  (vel: +0.0512°/g)
Nettuno     28.2° Psc  (vel: +0.0359°/g)
Plutone      2.0° Aqr  (vel: +0.0114°/g)
```

---

## 5.2 I flag di calcolo

Il terzo argomento di `calc_ut` è un numero intero che controlla *come* viene fatto il calcolo. Passare `0` significa "usa tutte le impostazioni di default": coordinate eclittiche, geocentriche, con aberrazione e nutazione.

Per cambiare comportamento, usi le costanti `FLG_*` e le combini con l'operatore `|` (OR bit a bit). Ad esempio, `ephem.FLG_EQUATORIAL | ephem.FLG_SPEED` chiede coordinate equatoriali con velocità. Puoi combinarne quante ne vuoi.

Ecco i più importanti, raggruppati per categoria.

**Coordinate** — cambiano il significato dei valori restituiti:

- Senza flag (default): coordinate **eclittiche** — longitudine, latitudine, distanza
- `FLG_EQUATORIAL`: coordinate **equatoriali** — Ascensione Retta, Declinazione, distanza
- `FLG_XYZ`: coordinate **cartesiane** — X, Y, Z in unità astronomiche
- `FLG_RADIANS`: angoli in **radianti** invece di gradi

**Centro di osservazione** — da dove stai "guardando":

- Senza flag (default): **geocentrico** — dal centro della Terra
- `FLG_HELCTR`: **eliocentrico** — dal centro del Sole
- `FLG_BARYCTR`: **baricentrico** — dal centro di massa del Sistema Solare
- `FLG_TOPOCTR`: **topocentrico** — dalla tua posizione sulla superficie terrestre (richiede prima `set_topo()`)

**Sistema di riferimento** — quale "punto zero" usare:

- Senza flag (default): **equinozio del giorno** ("of date") — il sistema dell'astrologia
- `FLG_J2000`: riferito al J2000.0 — il sistema dei cataloghi astronomici
- `FLG_ICRS`: sistema ICRS — quasi identico a J2000
- `FLG_SIDEREAL`: coordinate **siderali** — richiede prima `set_sid_mode()`
- `FLG_NONUT`: senza nutazione — restituisce coordinate "medie" invece di "vere"

**Correzioni** — attiva o disattiva effetti fisici:

- `FLG_SPEED`: calcola anche le **velocità** (pos[3], pos[4], pos[5]). Quasi sempre utile.
- `FLG_TRUEPOS`: posizione **geometrica vera**, senza correzione per il tempo-luce
- `FLG_NOABERR`: senza **aberrazione** (lo spostamento dovuto al moto terrestre)
- `FLG_NOGDEFL`: senza **deflessione gravitazionale** della luce

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Combinazione comune: equatoriali con velocità
pos, _ = ephem.calc_ut(jd, ephem.MARS,
                       ephem.FLG_EQUATORIAL | ephem.FLG_SPEED)

print(f"RA: {pos[0]:.4f}°, Dec: {pos[1]:+.4f}°")
print(f"Vel RA: {pos[3]:+.4f}°/g, Vel Dec: {pos[4]:+.4f}°/g")
```

```
RA: 344.6682°, Dec: -7.8866°
Vel RA: +0.7256°/g, Vel Dec: +0.2960°/g
```

---

## 5.3 Velocità e moto retrogrado

La velocità in longitudine (`pos[3]`) indica quanto il pianeta si sposta ogni giorno:

- **Positiva** = moto diretto (il pianeta avanza lungo lo zodiaco)
- **Negativa** = moto retrogrado (il pianeta sembra andare all'indietro)
- **Zero** = stazione (il pianeta "si ferma" prima di invertire direzione)

Il moto retrogrado è un effetto ottico: quando la Terra sorpassa un pianeta esterno (o viene sorpassata da uno interno), quel pianeta sembra muoversi all'indietro rispetto alle stelle. È come sorpassare un'auto in autostrada — per un istante sembra che l'altra auto vada indietro.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Controllare se Mercurio è retrogrado
retro = ephem.is_retrograde(ephem.MERCURY, jd)
print(f"Mercurio retrogrado: {retro}")

# Trovare la prossima stazione di Mercurio
jd_stazione, tipo = ephem.swe_find_station_ut(ephem.MERCURY, jd)

anno, mese, giorno, ore = ephem.revjul(jd_stazione)
# tipo = "SR" (stazione retrograda) o "SD" (stazione diretta)
print(f"Prossima stazione: {giorno}/{mese}/{anno} ({tipo})")
```

```
Mercurio retrogrado: True
Prossima stazione: 25/4/2024 (SD)
```

### 🌍 Vita reale

"Mercurio retrogrado" è uno dei concetti astrologici più popolari. Succede circa 3 volte l'anno, per circa 3 settimane ogni volta. In realtà tutti i pianeti (tranne Sole e Luna) hanno periodi retrogradi — i pianeti esterni sono retrogradi per mesi.

---

## 5.4 Fenomeni planetari: magnitudine, fase, elongazione

La funzione `pheno_ut` restituisce informazioni sull'aspetto fisico di un pianeta:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 21.0)

attr, _ = ephem.pheno_ut(jd, ephem.JUPITER, 0)

angolo_fase = attr[0]    # angolo Sole-Pianeta-Terra (gradi)
fase = attr[1]           # frazione illuminata del disco (0.0–1.0)
elongazione = attr[2]    # distanza angolare dal Sole (gradi)
diametro = attr[3]       # diametro apparente (arcsecondi)
magnitudine = attr[4]    # magnitudine visuale apparente

print(f"Giove:")
print(f"  Elongazione: {elongazione:.1f}°")
print(f"  Fase: {fase:.2f} ({fase*100:.0f}% illuminato)")
print(f"  Magnitudine: {magnitudine:.1f}")
print(f"  Diametro apparente: {diametro:.1f}\"")
```

```
Giove:
  Elongazione: 29.6°
  Fase: 1.00 (100% illuminato)
  Magnitudine: -2.0
  Diametro apparente: 0.0"
```

L'**elongazione** è la distanza angolare dal Sole — se è piccola (< 10°), il pianeta è perso nel bagliore solare e non visibile. La **magnitudine** indica la luminosità: numeri più bassi = più luminoso (Venere arriva a -4.6, Giove a -2.9).

---

## 5.5 Elementi orbitali

La funzione `get_orbital_elements` restituisce i parametri dell'orbita di un corpo celeste — i 6 numeri che descrivono un'ellisse nello spazio:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)
jd_tt = jd + ephem.deltat(jd)  # serve JD in TT

elem = ephem.get_orbital_elements(jd_tt, ephem.MARS, 0)

print(f"Marte - elementi orbitali:")
print(f"  Semi-asse maggiore:  {elem[0]:.6f} UA")
print(f"  Eccentricità:        {elem[1]:.6f}")
print(f"  Inclinazione:        {elem[2]:.4f}°")
print(f"  Periodo siderale:    {elem[10]:.2f} anni")
print(f"  Periodo sinodico:    {elem[13]:.1f} giorni")
print(f"  Dist. perielio:      {elem[15]:.4f} UA")
print(f"  Dist. afelio:        {elem[16]:.4f} UA")
```

```
Marte - elementi orbitali:
  Semi-asse maggiore:  1.523627 UA
  Eccentricità:        0.093278
  Inclinazione:        1.8479°
  Periodo siderale:    1.88 anni
  Periodo sinodico:    779.9 giorni
  Dist. perielio:      1.3815 UA
  Dist. afelio:        1.6657 UA
```

Per ottenere le distanze massima e minima con un'unica chiamata:

```python
d_max, d_min, d_now = ephem.orbit_max_min_true_distance(
    jd, ephem.MARS, 0
)
print(f"Marte: min {d_min:.4f} UA, max {d_max:.4f} UA, ora {d_now:.4f} UA")
```

```
Marte: min 0.3648 UA, max 2.6825 UA, ora 2.0618 UA
```

---

## 5.6 Nodi e apsidi

I **nodi** sono i punti dove l'orbita di un pianeta interseca l'eclittica. Le **apsidi** sono i punti di massima e minima distanza dal Sole (afelio e perielio).

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# nod_aps_ut restituisce 4 tuple di 6 valori ciascuna:
# nodo ascendente, nodo discendente, perielio, afelio
nasc, ndsc, peri, aphe = ephem.nod_aps_ut(
    jd, ephem.MARS, 0, 0  # method=0: medio
)

print(f"Marte - nodo ascendente: {nasc[0]:.4f}°")
print(f"Marte - nodo discendente: {ndsc[0]:.4f}°")
print(f"Marte - perielio: {peri[0]:.4f}° (dist {peri[2]:.4f} UA)")
print(f"Marte - afelio:   {aphe[0]:.4f}° (dist {aphe[2]:.4f} UA)")
```

```
Marte - nodo ascendente: 37.4197°
Marte - nodo discendente: 266.1986°
Marte - perielio: 354.2716° (dist 2.2240 UA)
Marte - afelio:   120.3522° (dist 1.1508 UA)
```

---

## Riepilogo

- `calc_ut(jd, body, flag)` è la funzione centrale della libreria. Dato un istante e un corpo celeste, restituisce 6 numeri: i primi tre dicono dove si trova (longitudine, latitudine, distanza), gli ultimi tre come si sta muovendo.
- I **flag** controllano il tipo di coordinate, il centro di osservazione, il sistema di riferimento e le correzioni fisiche. Si combinano con `|`. Passando `0` si ottengono coordinate eclittiche geocentriche — il default dell'astrologia.
- La **velocità in longitudine** (`pos[3]`) è positiva quando il pianeta avanza lungo lo zodiaco (moto diretto) e negativa quando sembra andare all'indietro (moto retrogrado).
- `pheno_ut` fornisce informazioni sull'aspetto fisico del pianeta: magnitudine (luminosità), fase (percentuale illuminata), elongazione (distanza angolare dal Sole) e diametro apparente.
- `get_orbital_elements` restituisce i parametri dell'orbita kepleriana: semi-asse maggiore, eccentricità, inclinazione, periodi e distanze.
- `nod_aps_ut` restituisce i nodi (dove l'orbita attraversa l'eclittica) e le apsidi (punti di minima e massima distanza dal Sole).

### Funzioni introdotte

- `pheno_ut(jd, body, flag)` — fenomeni: fase, elongazione, magnitudine
- `is_retrograde(body, jd)` — il pianeta è retrogrado in questo momento?
- `swe_find_station_ut(body, jd)` — trova la prossima stazione (retrograda o diretta)
- `get_orbital_elements(jd_tt, body, flag)` — elementi orbitali kepleriani
- `orbit_max_min_true_distance(jd, body, flag)` — distanze min, max e attuale dalla Terra
- `nod_aps_ut(jd, body, flag, method)` — nodi e apsidi dell'orbita
