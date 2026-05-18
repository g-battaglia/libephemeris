# Capitolo 8 — Le stelle fisse

## Cosa imparerai

In questo capitolo scoprirai perché le stelle "fisse" non sono davvero fisse, come la libreria calcola le loro posizioni tenendo conto della precessione e del moto proprio, come cercare una stella per nome (anche sbagliando l'ortografia), e come usare le stelle fisse nei calcoli astrologici.

---

## 8.1 Le stelle non sono così fisse

Gli antichi le chiamavano "stelle fisse" per distinguerle dai pianeti ("stelle erranti"), che si muovono visibilmente di notte in notte. Ma anche le stelle si muovono — solo molto, molto più lentamente.

### La precessione: tutto il cielo ruota

L'effetto più grande non è il movimento delle singole stelle, ma quello del nostro **punto di riferimento**. L'asse della Terra oscilla lentamente come una trottola, completando un giro in circa 26.000 anni. Questo fenomeno — la **precessione degli equinozi** — fa sì che il punto vernale (0° Ariete) si sposti di circa 50 secondi d'arco all'anno rispetto alle stelle.

In pratica, tutte le stelle si spostano di circa **1° ogni 72 anni** in longitudine eclittica. Regolo, che 2000 anni fa era nel cuore del Leone, oggi è a quasi 0° Vergine.

### Il moto proprio: ogni stella ha il suo cammino

Oltre alla precessione (che muove *tutte* le stelle allo stesso modo), ogni stella ha un suo **moto proprio** — il suo reale spostamento nello spazio. La velocità dipende dalla distanza e dalla velocità reale della stella:

- **Sirio**: 1.33"/anno — in 2000 anni si è spostata di quasi mezzo grado
- **Arturo**: 2.28"/anno — una delle stelle con il moto proprio più veloce tra quelle luminose
- **Regolo**: 0.25"/anno — relativamente "ferma"
- **Stella di Barnard**: 10.36"/anno — il record assoluto, ma è invisibile a occhio nudo

Per le stelle luminose usate in astrologia, il moto proprio è quasi sempre trascurabile su scale di qualche secolo. Ma per calcoli storici (temi natali dell'antichità, stelle nell'antico Egitto) diventa importante.

La libreria tiene conto sia della precessione sia del moto proprio quando calcola la posizione di una stella. Puoi anche propagare manualmente il moto proprio con `propagate_proper_motion`:

```python
import libephemeris as ephem

# Propagare il moto proprio di Sirio dal catalogo Hipparcos (J1991.25) a J2000
# RA e Dec in gradi, moto proprio in arcsec/anno
ra_sirio = 101.2872  # RA J1991.25 in gradi
dec_sirio = -16.7161  # Dec J1991.25 in gradi
pm_ra = -0.5461   # moto proprio in RA (include cos(dec)), arcsec/anno
pm_dec = -1.2232   # moto proprio in Dec, arcsec/anno

J1991_25 = 2448349.0625  # epoca del catalogo Hipparcos
J2000 = 2451545.0         # epoca J2000.0

ra_2000, dec_2000 = ephem.propagate_proper_motion(
    ra_sirio, dec_sirio,
    pm_ra, pm_dec,
    J1991_25, J2000
)

print(f"Sirio a J2000.0: RA = {ra_2000:.4f}°, Dec = {dec_2000:.4f}°")
```

```
Sirio a J2000.0: RA = 101.2858°, Dec = -16.7191°
```

Nella pratica quotidiana non avrai quasi mai bisogno di questa funzione: `fixstar2_ut` fa tutto automaticamente — precessione, moto proprio, nutazione e aberrazione inclusi.

---

## 8.2 Il catalogo stellare

La libreria include un catalogo interno di circa 100 stelle, basato sui dati del satellite **Hipparcos** dell'ESA (1989–1993). Hipparcos ha misurato la posizione, il moto proprio e la magnitudine di oltre 118.000 stelle con una precisione senza precedenti (circa 1 millesimo di secondo d'arco).

Il catalogo copre tutte le stelle importanti per l'astrologia e l'astronomia osservativa: le 15 stelle di Agrippa, le quattro stelle "regali" (Aldebaran, Regolo, Antares, Fomalhaut), le stelle più luminose del cielo, e le stelle zodiacali più usate.

### Calcolare la posizione di una stella

La funzione principale è `fixstar2_ut`. Passale il nome della stella, il Giorno Giuliano e i flag di calcolo:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Posizione di Regolo
pos, nome, retflags = ephem.fixstar2_ut("Regulus", jd, ephem.FLG_SPEED)

lon = pos[0]  # longitudine eclittica
lat = pos[1]  # latitudine eclittica
vel = pos[3]  # velocità in longitudine (gradi/giorno)

segni = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]
segno = segni[int(lon / 30)]
gradi = lon % 30

print(f"{nome}")
print(f"Posizione: {gradi:.2f}° {segno}")
print(f"Latitudine: {lat:.4f}°")
print(f"Velocità: {vel:.6f}°/giorno")
```

```
Regulus,alLeo
Posizione: 0.17° Vir
Latitudine: 0.4657°
Velocità: -0.000071°/giorno
```

La funzione restituisce tre valori:

- **pos** — una tupla di 6 valori: `(longitudine, latitudine, distanza, vel_lon, vel_lat, vel_dist)`. La distanza è fissata a 100.000 UA (le stelle sono effettivamente a distanza infinita per i nostri calcoli)
- **nome** — il nome completo della stella nel formato `"Nome,Nomenclatura"` (es. `"Regulus,alLeo"` dove `alLeo` sta per "alpha Leonis")
- **retflags** — i flag restituiti (di solito identici a quelli passati)

Se la stella non viene trovata, la funzione lancia un'eccezione `Error`.

### Le stelle principali dell'astrologia

Ecco le stelle più usate nella tradizione astrologica, con la loro posizione approssimativa nel 2024:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

stelle = [
    "Aldebaran",    # Occhio del Toro — "Guardiano dell'Est"
    "Regulus",      # Cuore del Leone — "Guardiano del Nord"
    "Antares",      # Cuore dello Scorpione — "Guardiano dell'Ovest"
    "Fomalhaut",    # Bocca del Pesce Australe — "Guardiano del Sud"
    "Sirius",       # La stella più luminosa del cielo
    "Algol",        # La "stella del demonio" — beta Persei
    "Spica",        # La spiga della Vergine
    "Betelgeuse",   # Spalla di Orione
    "Rigel",        # Piede di Orione
    "Vega",         # La stella della Lira
    "Polaris",      # La Stella Polare
    "Canopus",      # La seconda più luminosa del cielo
]

segni = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

for star in stelle:
    pos, nome, _ = ephem.fixstar2_ut(star, jd, 0)
    lon = pos[0]
    segno = segni[int(lon / 30)]
    gradi = lon % 30
    nome_corto = nome.split(",")[0]
    print(f"{nome_corto:12s}  {gradi:5.1f}° {segno}")
```

```
Aldebaran      10.1° Gem
Regulus         0.2° Vir
Antares        10.1° Sgr
Fomalhaut       4.2° Psc
Sirius         14.4° Cnc
Algol          26.5° Tau
Spica          24.2° Lib
Betelgeuse     29.1° Gem
Rigel          17.2° Gem
Vega           15.7° Cap
Polaris        28.9° Gem
Canopus        15.3° Cnc
```

### Magnitudine: quanto è luminosa?

La **magnitudine** indica la luminosità apparente di una stella. La scala è controintuitiva: numeri più **bassi** significano stelle più **luminose**. Sirio ha magnitudine −1.46 (molto luminosa), mentre le stelle appena visibili a occhio nudo hanno magnitudine circa 6.

```python
import libephemeris as ephem

# Magnitudine di alcune stelle
for star in ["Sirius", "Vega", "Polaris", "Algol"]:
    mag, nome = ephem.fixstar2_mag(star)
    nome_corto = nome.split(",")[0]
    print(f"{nome_corto:12s}  magnitudine {mag:+.2f}")
```

```
Sirius        magnitudine -1.46
Vega          magnitudine +0.03
Polaris       magnitudine +1.98
Algol         magnitudine +2.12
```

### Coordinate equatoriali

Se hai bisogno dell'Ascensione Retta e della Declinazione (per puntare un telescopio, per esempio), aggiungi il flag `FLG_EQUATORIAL`:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

pos, nome, _ = ephem.fixstar2_ut(
    "Sirius", jd, ephem.FLG_EQUATORIAL
)

ra = pos[0]   # Ascensione Retta in gradi
dec = pos[1]  # Declinazione in gradi

# Converti RA in ore, minuti, secondi
ra_ore = ra / 15.0
ore = int(ra_ore)
minuti = int((ra_ore - ore) * 60)
secondi = ((ra_ore - ore) * 60 - minuti) * 60

print(f"Sirio: RA = {ore}h {minuti}m {secondi:.1f}s, Dec = {dec:+.4f}°")
```

```
Sirio: RA = 6h 46m 12.6s, Dec = -16.7520°
```

---

## 8.3 Stelle fisse in astrologia

Nella tradizione astrologica, le stelle fisse hanno un ruolo particolare. Non si usano come i pianeti (con aspetti e transiti quotidiani), ma si considerano significative quando un pianeta o un angolo del tema natale si trova **in congiunzione** — cioè alla stessa longitudine eclittica — con una stella importante.

### Le quattro stelle regali

Le più potenti sono le quattro "stelle regali" o "guardiani del cielo", una per ogni quadrante dello zodiaco:

- **Aldebaran** (~10° Gemelli nel 2024) — l'Occhio del Toro, guardiano dell'equinozio di primavera. Associata al coraggio, all'onore militare e al successo attraverso l'integrità.

- **Regolo** (~0° Vergine) — il Cuore del Leone, guardiano del solstizio d'estate. La stella "più astrologica": associata al potere, alla fama e alla leadership. È quasi esattamente sull'eclittica (latitudine 0.46°), il che la rende particolarmente significativa.

- **Antares** (~10° Sagittario) — il Cuore dello Scorpione, guardiano dell'equinozio d'autunno. "Il rivale di Marte" (anti-Ares), associata all'intensità, alla passione e alla distruzione-rigenerazione.

- **Fomalhaut** (~4° Pesci) — la Bocca del Pesce Australe, guardiano del solstizio d'inverno. Associata all'idealismo, alla spiritualità e ai sogni.

### Congiunzione con un pianeta

In astrologia stellare, l'**orbe** (la distanza massima per considerare attiva un'influenza) è molto più stretto che per i pianeti: in genere **1°** per le stelle di prima magnitudine, **30'** per le altre. Ecco come verificare se oggi un pianeta è congiunto a una stella:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Stelle da controllare (le quattro regali + alcune importanti)
stelle_da_cercare = [
    "Aldebaran", "Regulus", "Antares", "Fomalhaut",
    "Sirius", "Spica", "Algol", "Vega"
]

pianeti = [
    (ephem.SUN,     "Sole"),
    (ephem.MOON,    "Luna"),
    (ephem.MERCURY, "Mercurio"),
    (ephem.VENUS,   "Venere"),
    (ephem.MARS,    "Marte"),
    (ephem.JUPITER, "Giove"),
    (ephem.SATURN,  "Saturno"),
]

# Calcola le posizioni dei pianeti
pos_pianeti = {}
for body_id, nome_p in pianeti:
    pos, _ = ephem.calc_ut(jd, body_id, ephem.FLG_SPEED)
    pos_pianeti[nome_p] = pos[0]  # longitudine eclittica

# Controlla congiunzioni
orbe = 1.0  # 1 grado di tolleranza

for star in stelle_da_cercare:
    pos_s, nome_s, _ = ephem.fixstar2_ut(star, jd, 0)
    lon_s = pos_s[0]
    nome_corto = nome_s.split(",")[0]

    for nome_p, lon_p in pos_pianeti.items():
        # Distanza angolare (con gestione del "salto" 360°→0°)
        diff = abs(lon_p - lon_s)
        if diff > 180:
            diff = 360 - diff

        if diff <= orbe:
            print(f"★ {nome_p} congiunto a {nome_corto}! "
                  f"(distanza: {diff:.2f}°)")
```

```
Nessuna congiunzione entro 1° trovata in questa data.
```

---

## 8.4 Cercare una stella: nomi, designazioni e ricerca fuzzy

La libreria è molto flessibile nel modo in cui accetta i nomi delle stelle. Puoi usare:

**Nome tradizionale** — il modo più semplice e intuitivo:

```python
ephem.fixstar2_ut("Sirius", jd, 0)
ephem.fixstar2_ut("Regulus", jd, 0)
ephem.fixstar2_ut("Betelgeuse", jd, 0)
```

**Designazione Bayer** — la lettera greca + la costellazione, usata nei cataloghi astronomici. Alpha Leonis = la stella più luminosa del Leone:

```python
ephem.fixstar2_ut("Alpha Leonis", jd, 0)     # = Regulus
ephem.fixstar2_ut("Beta Persei", jd, 0)       # = Algol
ephem.fixstar2_ut("Alpha Canis Majoris", jd, 0)  # = Sirius
```

**Codice di nomenclatura** — la versione abbreviata della designazione Bayer, nel formato usato internamente (es. `alLeo` = alpha Leonis):

```python
ephem.fixstar2_ut("alLeo", jd, 0)   # = Regulus
ephem.fixstar2_ut("bePer", jd, 0)   # = Algol
```

**Numero HIP** — il numero del catalogo Hipparcos, per chi lavora con dati astronomici:

```python
ephem.fixstar2_ut("HIP 49669", jd, 0)  # = Regulus
ephem.fixstar2_ut("32349", jd, 0)       # = Sirius (HIP 32349)
```

**Ricerca parziale** — se ricordi solo l'inizio del nome:

```python
ephem.fixstar2_ut("Reg", jd, 0)    # trova Regulus (se non ambiguo)
ephem.fixstar2_ut("Alde", jd, 0)   # trova Aldebaran
```

### Ricerca fonetica: perdona gli errori di ortografia

I nomi delle stelle vengono dall'arabo, dal latino e dal greco, e le traslitterazioni variano. La libreria include un sistema di **ricerca fonetica** che trova la stella giusta anche se sbagli l'ortografia:

```python
# Tutti questi trovano Betelgeuse:
ephem.fixstar2_ut("Betelgeuse", jd, 0)   # ortografia corretta
ephem.fixstar2_ut("Betelgeuze", jd, 0)   # con la z (variante tedesca)
ephem.fixstar2_ut("Betelgeux", jd, 0)    # troncato

# E questi trovano Fomalhaut:
ephem.fixstar2_ut("Fomalhaut", jd, 0)    # ortografia corretta
ephem.fixstar2_ut("Formalhaut", jd, 0)   # con la r in più
```

Il sistema fonetico normalizza il nome rimuovendo le doppie consonanti, unificando le vocali simili, e confrontando lo "scheletro consonantico" del nome. Funziona bene per le varianti ortografiche comuni, ma restituisce un errore se la ricerca è ambigua (più stelle corrispondono allo stesso pattern).

### Il formato del nome restituito

Quando la ricerca ha successo, il primo valore restituito è il nome completo nel formato `"Nome,Nomenclatura"`:

```python
pos, nome, _ = ephem.fixstar2_ut("Sirius", jd, 0)
print(nome)  # "Sirius,alCMa" → alpha Canis Majoris

pos, nome, _ = ephem.fixstar2_ut("Algol", jd, 0)
print(nome)  # "Algol,bePer" → beta Persei
```

```
Sirius,alCMa
Algol,bePer
```

La parte dopo la virgola è la designazione Bayer abbreviata: `al` = alpha, `be` = beta, `ga` = gamma, e così via. Il codice della costellazione segue l'abbreviazione standard IAU a tre lettere (Leo, Per, CMa, Vir, Ori...).

---

## Riepilogo

In questo capitolo abbiamo imparato a lavorare con le stelle fisse.

**Concetti chiave:**

- Le stelle "fisse" si muovono: la **precessione** sposta tutte le stelle di ~50"/anno in longitudine eclittica, e il **moto proprio** aggiunge uno spostamento individuale per ogni stella
- La libreria include un catalogo di ~100 stelle basato sui dati del satellite Hipparcos, con precisione sub-millesimo di secondo d'arco
- In astrologia, le stelle fisse contano quando un pianeta o un angolo è in **congiunzione** (entro ~1°) con una stella importante
- Le quattro stelle regali (Aldebaran, Regolo, Antares, Fomalhaut) sono considerate le più potenti
- La **magnitudine** indica la luminosità: numeri più bassi = stella più luminosa (Sirio = −1.46, limite occhio nudo ≈ 6)
- La ricerca è flessibile: nome tradizionale, designazione Bayer, numero HIP, ricerca parziale, e persino ricerca fonetica per nomi scritti male

**Funzioni introdotte:**

- `fixstar2_ut(nome, jd, flag)` — calcola la posizione eclittica (o equatoriale con `FLG_EQUATORIAL`) di una stella, accettando nomi in molti formati
- `fixstar2_mag(nome)` — restituisce la magnitudine visuale di una stella
- `propagate_proper_motion(ra, dec, pm_ra, pm_dec, from_jd, to_jd)` — propaga manualmente il moto proprio di una stella tra due epoche
