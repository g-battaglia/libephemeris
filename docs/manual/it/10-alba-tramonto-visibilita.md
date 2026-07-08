# Capitolo 10 — Alba, tramonto e visibilità

## Cosa imparerai

In questo capitolo scoprirai come calcolare l'alba, il tramonto e i transiti al meridiano di qualsiasi corpo celeste, cosa sono i diversi tipi di crepuscolo, come funziona la rifrazione atmosferica (che fa "alzare" gli oggetti vicino all'orizzonte), e come determinare se una stella o un pianeta è visibile a occhio nudo in un dato momento.

---

## 10.1 Alba e tramonto

### Cosa significa "il Sole sorge"?

La definizione sembra ovvia, ma in realtà è sorprendentemente complicata. Il momento dell'alba non è quando il *centro* del Sole tocca l'orizzonte, ma quando il **bordo superiore** (il "lembo") del disco solare appare all'orizzonte. Dato che il diametro apparente del Sole è circa 32 minuti d'arco (0.53°), il centro è ancora circa 16' sotto l'orizzonte geometrico quando vediamo il primo raggio.

Ma c'è di più: l'atmosfera terrestre **curva la luce** (rifrazione), facendo apparire gli oggetti vicini all'orizzonte più alti di quanto siano realmente. All'orizzonte, la rifrazione solleva il Sole di circa **34 minuti d'arco**. Quindi, nel momento in cui "vediamo" il Sole sorgere, il suo centro è in realtà circa 50' (quasi un grado!) sotto l'orizzonte geometrico.

La libreria tiene conto di tutto questo automaticamente.

### La funzione `rise_trans`

```python
import libephemeris as ephem

# Alba e tramonto del Sole oggi a Milano
jd = ephem.julday(2024, 6, 21, 0.0)  # Solstizio d'estate
lat, lon = 45.4642, 9.1900  # Milano

# Alba (rise)
ret, tret = ephem.rise_trans(
    jd, ephem.SUN,
    ephem.CALC_RISE,
    (lon, lat, 0.0),
    atpress=1013.25,
    attemp=15.0
)
jd_alba = tret[0]  # il giorno giuliano dell'evento

# Tramonto (set)
ret, tret = ephem.rise_trans(
    jd, ephem.SUN,
    ephem.CALC_SET,
    (lon, lat, 0.0)
)
jd_tramonto = tret[0]

# Converti in ore locali (CEST = UT + 2)
fuso = 2  # ora legale estiva

_, _, _, h_alba = ephem.revjul(jd_alba)
_, _, _, h_tram = ephem.revjul(jd_tramonto)

ore_alba = h_alba + fuso
ore_tram = h_tram + fuso
durata = (jd_tramonto - jd_alba) * 24

print(f"Alba:     {int(ore_alba):02d}:{int((ore_alba % 1) * 60):02d} ora locale")
print(f"Tramonto: {int(ore_tram):02d}:{int((ore_tram % 1) * 60):02d} ora locale")
print(f"Durata del giorno: {durata:.1f} ore")
```

```
Alba:     05:34 ora locale
Tramonto: 21:15 ora locale
Durata del giorno: 15.7 ore
```

Il parametro `rsmi` accetta quattro eventi base:

- `CALC_RISE` (1) — levata (primo raggio all'orizzonte)
- `CALC_SET` (2) — tramonto (ultimo raggio scompare)
- `CALC_MTRANSIT` (4) — transito superiore al meridiano (punto più alto nel cielo)
- `CALC_ITRANSIT` (8) — transito inferiore (punto più basso, di solito sotto l'orizzonte)

A questi puoi aggiungere (con l'operatore `|`) dei modificatori:

- `BIT_DISC_CENTER` (256) — usa il centro del disco invece del bordo superiore
- `BIT_DISC_BOTTOM` (8192) — usa il bordo inferiore del disco
- `BIT_NO_REFRACTION` (512) — non applicare la rifrazione atmosferica

### Levata e tramonto della Luna

Funziona esattamente come per il Sole:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 0.0)
lat, lon = 41.9028, 12.4964  # Roma

ret, tret = ephem.rise_trans(
    jd, ephem.MOON, ephem.CALC_RISE,
    (lon, lat, 0.0)
)
jd_rise = tret[0]

ret, tret = ephem.rise_trans(
    jd, ephem.MOON, ephem.CALC_SET,
    (lon, lat, 0.0)
)
jd_set = tret[0]

if jd_rise > 0 and jd_set > 0:
    _, _, _, h_rise = ephem.revjul(jd_rise)
    _, _, _, h_set = ephem.revjul(jd_set)
    print(f"Sorgere della Luna: {int(h_rise):02d}:{int((h_rise % 1) * 60):02d} UT")
    print(f"Tramonto della Luna: {int(h_set):02d}:{int((h_set % 1) * 60):02d} UT")
```

```
Sorgere della Luna: 04:29 UT
Tramonto della Luna: 17:36 UT
```

La Luna ha una particolarità: il suo sorgere ritarda di circa 50 minuti al giorno. Quindi a volte in una giornata non sorge affatto, o sorge ma non tramonta (o viceversa).

### Corpi circumpolari

Alle latitudini alte, alcuni corpi non sorgono mai o non tramontano mai. Per esempio, al polo nord il Sole resta sopra l'orizzonte per sei mesi. Quando questo succede, `rise_trans` restituisce un `retflag` di `-2`:

```python
import libephemeris as ephem

# Sole a Tromsø il 21 giugno — non tramonta mai
jd = ephem.julday(2024, 6, 21, 0.0)
ret, tret = ephem.rise_trans(
    jd, ephem.SUN, ephem.CALC_SET,
    (19.0, 69.6, 0.0)
)

if ret == -2:
    print("Il Sole non tramonta! (circumpolare)")
```

```
Il Sole non tramonta! (circumpolare)
```

---

## 10.2 Transiti al meridiano

Il **transito superiore** (o culminazione) è il momento in cui un corpo celeste raggiunge il punto più alto nel cielo, attraversando il meridiano locale. Per il Sole, questo è il **mezzogiorno solare** — che non coincide quasi mai con le 12:00 dell'orologio.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 0.0)
lat, lon = 41.9028, 12.4964  # Roma

# Transito del Sole al meridiano = mezzogiorno solare
ret, tret = ephem.rise_trans(
    jd, ephem.SUN, ephem.CALC_MTRANSIT,
    (lon, lat, 0.0)
)
jd_transit = tret[0]

_, _, _, h = ephem.revjul(jd_transit)
print(f"Mezzogiorno solare a Roma: {int(h):02d}:{int((h % 1) * 60):02d} UT")
print(f"In ora locale (CET): {int(h+1):02d}:{int((h % 1) * 60):02d}")
```

```
Mezzogiorno solare a Roma: 11:11 UT
In ora locale (CET): 12:11
```

Il mezzogiorno solare varia durante l'anno a causa dell'**equazione del tempo** (vedi Capitolo 2): può essere fino a 16 minuti prima o 14 minuti dopo le 12:00 del meridiano locale. E se il tuo fuso orario non corrisponde esattamente alla tua longitudine, la differenza è ancora maggiore.

---

## 10.3 I crepuscoli

Il cielo non diventa buio immediatamente dopo il tramonto. La transizione dalla luce alla notte avviene gradualmente, e gli astronomi distinguono tre fasi di crepuscolo, definite dalla posizione del Sole sotto l'orizzonte:

**Crepuscolo civile** — il Sole è tra 0° e −6° sotto l'orizzonte. C'è ancora abbastanza luce per la maggior parte delle attività all'aperto. Si possono leggere un libro e distinguere i contorni degli edifici. I pianeti più luminosi (Venere, Giove) cominciano a comparire.

**Crepuscolo nautico** — il Sole è tra −6° e −12°. L'orizzonte marino è ancora visibile — da qui il nome: i marinai potevano ancora usare il sestante per la navigazione. Le stelle più luminose sono visibili.

**Crepuscolo astronomico** — il Sole è tra −12° e −18°. Il cielo sembra buio a occhio nudo, ma c'è ancora un leggero chiarore residuo all'orizzonte. Le stelle deboli non sono ancora visibili.

**Notte astronomica** — il Sole è sotto −18°. Il cielo è completamente buio (a parte l'inquinamento luminoso). Tutte le stelle visibili a occhio nudo sono visibili.

Per calcolare i crepuscoli, usa `rise_trans` con i flag di crepuscolo:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 6, 21, 0.0)  # Solstizio d'estate
lat, lon = 45.4642, 9.1900  # Milano
fuso = 2  # CEST

def ora_locale(jd_evento, fuso):
    _, _, _, h = ephem.revjul(jd_evento)
    h += fuso
    return f"{int(h):02d}:{int((h % 1) * 60):02d}"

# Tramonto
_, tret = ephem.rise_trans(jd, ephem.SUN,
    ephem.CALC_SET, (lon, lat, 0.0))
jd_t = tret[0]

# Fine crepuscolo civile (Sole a -6°)
_, tret = ephem.rise_trans(jd, ephem.SUN,
    ephem.CALC_SET | ephem.BIT_CIVIL_TWILIGHT, (lon, lat, 0.0))
jd_civ = tret[0]

# Fine crepuscolo nautico (Sole a -12°)
_, tret = ephem.rise_trans(jd, ephem.SUN,
    ephem.CALC_SET | ephem.BIT_NAUTIC_TWILIGHT, (lon, lat, 0.0))
jd_nau = tret[0]

# Fine crepuscolo astronomico (Sole a -18°)
ret, tret = ephem.rise_trans(jd, ephem.SUN,
    ephem.CALC_SET | ephem.BIT_ASTRO_TWILIGHT, (lon, lat, 0.0))
jd_ast = tret[0]

print(f"Tramonto:                {ora_locale(jd_t, fuso)}")
print(f"Fine crepuscolo civile:  {ora_locale(jd_civ, fuso)}")
print(f"Fine crepuscolo nautico: {ora_locale(jd_nau, fuso)}")

if ret != -2:
    print(f"Fine crepuscolo astron.: {ora_locale(jd_ast, fuso)}")
else:
    print("Notte astronomica: MAI (il cielo non diventa completamente buio)")
```

```
Tramonto:                21:15
Fine crepuscolo civile:  21:55
Fine crepuscolo nautico: 22:46
Fine crepuscolo astron.: 23:57
```

Al solstizio d'estate a Milano (45° N), il crepuscolo astronomico non finisce mai del tutto — il Sole non scende mai sotto −18°. Le notti bianche delle latitudini nordiche sono una conseguenza di questo fenomeno.

---

## 10.4 Rifrazione atmosferica

La **rifrazione** è la curvatura della luce causata dall'atmosfera terrestre. L'atmosfera è più densa in basso e più rarefatta in alto, e la luce segue una curva invece di una retta. L'effetto è massimo all'orizzonte e diminuisce man mano che si guarda più in alto:

- All'orizzonte (0°): **~34'** di rifrazione — quasi il diametro del Sole!
- A 10° di altezza: ~5'
- A 45°: ~1'
- Allo zenit (90°): 0' (la luce arriva verticalmente, nessuna curvatura)

La funzione `refrac` converte tra altezza "vera" (geometrica) e altezza "apparente" (quello che vedi):

```python
import libephemeris as ephem

# Da altezza vera ad apparente
alt_vera = 0.0  # oggetto esattamente all'orizzonte geometrico
alt_apparente = ephem.refrac(alt_vera, calc_flag=ephem.TRUE_TO_APP)
print(f"Altezza vera: {alt_vera:.2f}° → apparente: {alt_apparente:.2f}°")
# → circa 0.57° (34 arcminuti sopra l'orizzonte)

# Da apparente a vera
alt_app = 5.0
alt_true = ephem.refrac(alt_app, calc_flag=ephem.APP_TO_TRUE)
print(f"Altezza apparente: {alt_app:.2f}° → vera: {alt_true:.2f}°")
```

```
Altezza vera: 0.00° → apparente: 0.48°
Altezza apparente: 5.00° → vera: 4.84°
```

La rifrazione dipende dalle condizioni atmosferiche — pressione e temperatura:

```python
import libephemeris as ephem

# Rifrazione a diverse temperature (all'orizzonte)
for temp in [-20, 0, 15, 35]:
    r = ephem.refrac(0.0, pressure=1013.25, temperature=temp)
    print(f"T = {temp:+3d}°C → rifrazione all'orizzonte: {r:.2f}°")
```

```
T = -20°C → rifrazione all'orizzonte: 0.54°
T =  +0°C → rifrazione all'orizzonte: 0.50°
T = +15°C → rifrazione all'orizzonte: 0.48°
T = +35°C → rifrazione all'orizzonte: 0.45°
```

### Rifrazione estesa: osservatori in quota

Se osservi da una montagna o da un aereo, l'orizzonte è più basso del normale (l'"abbassamento dell'orizzonte" o "dip"). La funzione `refrac_extended` tiene conto di questo:

```python
import libephemeris as ephem

# Osservazione dalla cima del Monte Bianco (4810 m)
alt_obj = 0.5  # oggetto mezzo grado sopra l'orizzonte geometrico

alt_result, dettagli = ephem.refrac_extended(
    alt_obj,
    altitude_geo=4810.0,  # altitudine dell'osservatore in metri
    pressure=550.0,       # pressione ridotta in quota
    temperature=-10.0     # freddo!
)

alt_vera, alt_app, rifrazione, dip = dettagli

print(f"Altezza vera:     {alt_vera:.3f}°")
print(f"Altezza apparente: {alt_app:.3f}°")
print(f"Rifrazione:        {rifrazione:.3f}°")
print(f"Dip dell'orizzonte: {dip:.3f}°")
```

```
Altezza vera:     0.500°
Altezza apparente: 0.744°
Rifrazione:        0.244°
Dip dell'orizzonte: -1.926°
```

### Orizzonte personalizzato

Se il tuo orizzonte non è piatto (montagne, edifici), puoi usare `rise_trans_true_hor` per specificare un'altezza personalizzata dell'orizzonte:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 0.0)
lat, lon = 45.4642, 9.1900  # Milano

# Le Alpi a nord alzano l'orizzonte di circa 2 gradi
_, tret = ephem.rise_trans_true_hor(
    jd, ephem.SUN, ephem.CALC_RISE,
    (lon, lat, 0.0),
    horhgt=2.0,  # orizzonte alzato di 2° sopra il piano
)
jd_alba = tret[0]

_, _, _, h = ephem.revjul(jd_alba)
print(f"Alba con orizzonte a 2°: {int(h):02d}:{int((h % 1) * 60):02d} UT")
```

```
Alba con orizzonte a 2°: 05:01 UT
```

---

## 10.5 Visibilità eliacale

La **levata eliacale** è uno dei concetti più antichi dell'astronomia. È il primo giorno dell'anno in cui una stella (o un pianeta) diventa visibile all'alba, dopo essere stata nascosta dal bagliore del Sole per settimane o mesi.

Per gli antichi Egizi, la levata eliacale di **Sirio** (la stella più luminosa) segnava l'inizio dell'anno e annunciava l'imminente piena del Nilo — un evento di importanza vitale. Per i Babilonesi, la levata eliacale di Venere era un presagio fondamentale.

La visibilità eliacale dipende da molti fattori: la luminosità dell'oggetto, la sua distanza dal Sole, la brillanza del cielo al crepuscolo, le condizioni atmosferiche, e persino l'acutezza visiva dell'osservatore.

### Trovare la levata eliacale

La funzione `heliacal_ut` cerca la data di un evento eliacale:

```python
import libephemeris as ephem

# Quando Sirio diventa visibile all'alba nel 2024?
jd = ephem.julday(2024, 1, 1, 0.0)
lat, lon = 30.0, 31.2  # Il Cairo (dove gli Egizi la osservavano)

jd_evento, *_ = ephem.heliacal_ut(
    jd,
    geopos=(31.2, 30.0, 0.0),      # Il Cairo (lon, lat, alt)
    datm=(1013.25, 25.0, 30.0, 0.0),  # pressione, temp, umidità%, lapse
    dobs=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    object_name="Sirius",
    event_type=ephem.HELIACAL_RISING  # levata eliacale mattutina
)

if jd_evento > 0:
    y, m, d, h = ephem.revjul(jd_evento)
    print(f"Levata eliacale di Sirio al Cairo: {d:.0f}/{m:.0f}/{y:.0f}")
```

```
Levata eliacale di Sirio al Cairo: 4/8/2024
```

I quattro tipi di evento eliacale sono:

- `HELIACAL_RISING` (1) — **levata eliacale**: la prima apparizione mattutina. L'oggetto diventa visibile all'alba dopo il periodo di invisibilità vicino al Sole. Valido per pianeti e stelle fisse (non per Sole o Luna).

- `HELIACAL_SETTING` (2) — **tramonto eliaco**: l'ultima apparizione serale. L'oggetto scompare nel bagliore del Sole al tramonto. Valido per pianeti e stelle fisse (non per Sole o Luna).

- `EVENING_FIRST` (3) — **prima apparizione serale**: il primo giorno in cui il corpo diventa visibile la sera dopo la congiunzione superiore. Valido per i pianeti interni (Mercurio, Venere) e la Luna.

- `MORNING_LAST` (4) — **ultima apparizione mattutina**: l'ultimo giorno di visibilità al mattino prima della congiunzione inferiore. Valido per i pianeti interni (Mercurio, Venere) e la Luna.

### L'API compatibile con l'efemeride di riferimento

Per controllo completo sulle condizioni atmosferiche e le capacità dell'osservatore, usa `heliacal_ut`:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 1, 1, 0.0)

geopos = (31.2, 30.0, 0.0)  # (lon, lat, alt) — attenzione: longitudine PRIMA

datm = (1013.25, 25.0, 30.0, 0.0)
#  pressione, temperatura, umidità%, visibilità meteorologica

dobs = (36.0, 1.0, 0, 0, 0, 0)
#  età, acuità visiva (Snellen), binoculare, ingrandimento, apertura, trasmissione

result = ephem.heliacal_ut(
    jd, geopos, datm, dobs,
    "Sirius",
    ephem.HELIACAL_RISING
)

jd_inizio = result[0]
if jd_inizio > 0:
    y, m, d, h = ephem.revjul(jd_inizio)
    print(f"Levata eliacale: {d:.0f}/{m:.0f}/{y:.0f}")
```

```
Levata eliacale: 4/8/2024
```

La tupla `datm` (condizioni atmosferiche):
- `[0]` — pressione atmosferica in mbar
- `[1]` — temperatura in °C
- `[2]` — umidità relativa in percentuale (0–100)
- `[3]` — visibilità meteorologica in km (0 = calcola automaticamente)

La tupla `dobs` (osservatore):
- `[0]` — età dell'osservatore in anni (influenza l'acuità visiva)
- `[1]` — rapporto di Snellen (1.0 = vista normale, 1.5 = vista eccellente)
- `[2]` — 0 = monoculare, 1 = binoculare (solo con telescopio)
- `[3]` — ingrandimento del telescopio (0 = occhio nudo)
- `[4]` — apertura del telescopio in mm
- `[5]` — coefficiente di trasmissione ottica

---

## 10.6 Estinzione atmosferica e luminosità del cielo

### Quanto si attenua la luce?

L'atmosfera non solo curva la luce (rifrazione), ma ne **assorbe** una parte — soprattutto vicino all'orizzonte, dove la luce attraversa uno strato d'aria più spesso. Questo effetto si chiama **estinzione atmosferica**.

L'**airmass** misura quanto atmosfera la luce deve attraversare: allo zenit vale 1, all'orizzonte circa 38:

```python
import libephemeris as ephem

# Airmass a diverse altezze
for alt in [90, 60, 30, 10, 5, 1]:
    am = ephem.calc_airmass(float(alt))
    print(f"Altezza {alt:2d}° → airmass = {am:.2f}")
```

```
Altezza 90° → airmass = 1.00
Altezza 60° → airmass = 1.15
Altezza 30° → airmass = 1.99
Altezza 10° → airmass = 5.59
Altezza  5° → airmass = 10.31
Altezza  1° → airmass = 26.31
```

L'estinzione totale dipende da airmass e condizioni atmosferiche:

```python
import libephemeris as ephem

# Quante magnitudini di luce si perdono a 10° di altezza?
ext = ephem.calc_extinction_magnitude(10.0)
print(f"Estinzione a 10° di altezza: {ext:.2f} magnitudini")

# Una stella di magnitudine 4.0 diventa:
mag_app = ephem.apparent_magnitude_with_extinction(4.0, 10.0)
print(f"Magnitudine 4.0 → {mag_app:.2f} dopo estinzione")
```

```
Estinzione a 10° di altezza: 1.58 magnitudini
Magnitudine 4.0 → 5.58 dopo estinzione
```

### Il cielo durante il crepuscolo

La luminosità del cielo durante il crepuscolo determina quali oggetti sono visibili. Un cielo luminoso "nasconde" le stelle deboli:

```python
import libephemeris as ephem

# Brillanza del cielo al crepuscolo
for sun_alt in [0, -3, -6, -9, -12, -15, -18]:
    fase = ephem.get_twilight_phase(float(sun_alt))
    lim_mag = ephem.calc_limiting_magnitude_twilight(float(sun_alt))
    print(f"Sole a {sun_alt:+3d}° → {fase:14s}, "
          f"magnitudine limite: {lim_mag:.1f}")
```

```
Sole a  +0° → civil         , magnitudine limite: -0.7
Sole a  -3° → civil         , magnitudine limite: 0.6
Sole a  -6° → nautical      , magnitudine limite: 1.8
Sole a  -9° → nautical      , magnitudine limite: 3.7
Sole a -12° → astronomical  , magnitudine limite: 5.3
Sole a -15° → astronomical  , magnitudine limite: 5.9
Sole a -18° → night         , magnitudine limite: 6.4
```

La **magnitudine limite** è la magnitudine della stella più debole visibile. A notte piena, in condizioni perfette, è circa 6.0–6.5. Durante il crepuscolo civile, scende a 1–2 (solo le stelle più luminose).

### Visibilità di un pianeta

Per determinare se un oggetto è visibile a occhio nudo in condizioni specifiche, usa `vis_limit_mag`:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 12, 20, 17.0)  # sera, 18:00 CET

geopos = (12.4964, 41.9028, 0.0)  # Roma (lon, lat, alt)
datm = (1013.25, 10.0, 50.0, 0.0)
dobs = (36.0, 1.0, 0, 0, 0, 0)

result, dret = ephem.vis_limit_mag(
    jd, geopos, datm, dobs, "Venus"
)

mag_limite = dret[0]  # magnitudine limite
mag_oggetto = dret[7]  # magnitudine apparente dell'oggetto
alt_oggetto = dret[1]  # altezza dell'oggetto

if result >= 0:  # non sotto l'orizzonte
    visibile = mag_limite > mag_oggetto
    print(f"Magnitudine limite: {mag_limite:.1f}")
    print(f"Magnitudine di Venere: {mag_oggetto:.1f}")
    print(f"Altezza: {alt_oggetto:.1f}°")
    print(f"Visibile: {'SÌ' if visibile else 'NO'}")
else:
    print("Venere è sotto l'orizzonte")
```

```
Magnitudine limite: 5.9
Magnitudine di Venere: -3.5
Altezza: 20.0°
Visibile: SÌ
```

---

## Riepilogo

In questo capitolo abbiamo imparato a calcolare quando i corpi celesti sono visibili nel cielo.

**Concetti chiave:**

- L'**alba** avviene quando il bordo superiore del disco solare appare all'orizzonte, tenendo conto della rifrazione atmosferica (~34' di "alzamento")
- I **crepuscoli** sono tre fasi progressive: civile (Sole a −6°), nautico (−12°), astronomico (−18°). Solo sotto −18° il cielo è davvero buio
- La **rifrazione** atmosferica fa apparire gli oggetti vicini all'orizzonte più alti di quanto siano realmente; l'effetto è massimo all'orizzonte (~34')
- L'**estinzione** atmosferica attenua la luce degli oggetti bassi sull'orizzonte; l'airmass allo zenit è 1, all'orizzonte ~38
- La **levata eliacale** è la prima apparizione di un corpo all'alba dopo il periodo di invisibilità vicino al Sole — un concetto antico quanto l'astronomia stessa
- La **magnitudine limite** dipende dalla luminosità del cielo, dalle condizioni atmosferiche e dall'acutezza visiva dell'osservatore

**Funzioni introdotte:**

- `rise_trans(jd, body, rsmi, geopos)` — trova il prossimo sorgere, tramonto o transito al meridiano di un corpo celeste
- `rise_trans_true_hor(jd, body, rsmi, geopos, horhgt=0.0)` — come `rise_trans` ma con altezza personalizzata dell'orizzonte
- `refrac(altitude, pressure, temperature, calc_flag)` — converte tra altezza vera e apparente (o viceversa), tenendo conto della rifrazione
- `refrac_extended(altitude, altitude_geo, ...)` — rifrazione estesa con dip dell'orizzonte per osservatori in quota
- `heliacal_ut(jd, geopos, datm, dobs, object_name, event_type)` — trova la data di un evento eliacale (prima/ultima visibilità)
- `heliacal_ut(jd, geopos, datm, dobs, object_name, event_type)` — stessa funzione (il nome semplice è ora un alias)
- `vis_limit_mag(jd, geopos, atmo, observer, objname)` — determina se un oggetto è visibile, confrontando magnitudine limite e magnitudine apparente
- `calc_airmass(altitude)` — calcola la massa d'aria attraversata dalla luce
- `calc_extinction_magnitude(altitude)` — calcola quante magnitudini di luce si perdono per estinzione atmosferica
- `get_twilight_phase(sun_altitude)` — restituisce la fase del crepuscolo in base all'altezza del Sole
