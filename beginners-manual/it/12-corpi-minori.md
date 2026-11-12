# Capitolo 12 — Corpi minori: asteroidi, centauri, TNO

[Indice e configurazione](README.md) · [Precedente](11-zodiaco-siderale.md) · [Successivo](13-pianeti-ipotetici-parti-arabe.md) · [English](../en/12-minor-bodies.md)

## Cosa imparerai

In questo capitolo scoprirai come la libreria gestisce migliaia di corpi celesti oltre ai pianeti principali: asteroidi della fascia principale, centauri, oggetti trans-nettuniani (TNO) e asteroidi near-Earth. Imparerai la "catena di calcolo" che la libreria usa per trovare la posizione più precisa possibile, e come scaricare e registrare kernel SPK per i corpi che ti interessano.

---

## 12.1 La gerarchia dei corpi celesti

Non tutti i corpi del Sistema Solare sono uguali in termini di precisione disponibile:

**Pianeti (Mercurio–Nettuno)**: il nucleo incluso copre gli esercizi moderni; altre date e centri planetari richiedono dati adeguati. Incertezza della sorgente ed errore di interpolazione sono distinti: la disponibilità non implica un'accuratezza sub-millimetrica universale.

**Plutone e Chirone**: il sistema di Plutone appartiene alle efemeridi planetarie; Chirone è un corpo minore con dati propri, non un pianeta standard DE440. Gli esempi richiedono i gruppi LEB asteroidi/esotici appropriati o SPK registrati esplicitamente. Preparali seguendo [la configurazione degli esempi](README.md#dati-aggiuntivi-ed-esempi-online).

**Asteroidi maggiori (Cerere, Pallade, Giunone, Vesta)**: hanno ID dedicati (`CERES=17`, `PALLAS=18`, `JUNO=19`, `VESTA=20`). La libreria può scaricare kernel SPK ad alta precisione da NASA/JPL.

**Centauri (Pholus, Nessus, Chariklo)**: corpi con orbite tra Giove e Nettuno. Disponibili con kernel SPK scaricabili o con il fallback kepleriano.

**Trans-Nettuniani (Eris, Sedna, Haumea, Makemake, Quaoar)**: oltre Nettuno. Stessi metodi dei centauri.

**Near-Earth (Apophis, Bennu, Eros)**: asteroidi con orbite vicine alla Terra. Particolarmente interessanti per le missioni spaziali.

```python
import libephemeris as ephem

# Gli ID dedicati per i corpi più importanti
corpi = [
    (ephem.CHIRON,  "Chirone"),
    (ephem.PHOLUS,  "Pholus"),
    (ephem.CERES,   "Cerere"),
    (ephem.PALLAS,  "Pallade"),
    (ephem.JUNO,    "Giunone"),
    (ephem.VESTA,   "Vesta"),
]

jd = ephem.julday(2024, 4, 8, 12.0)

segni = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

for body_id, nome in corpi:
    pos, _ = ephem.calc_ut(jd, body_id, ephem.FLG_SPEED)
    segno = segni[int(pos[0] / 30)]
    gradi = pos[0] % 30
    print(f"{nome:10s}  {gradi:5.1f}° {segno}")
```

```
Chirone      19.4° Ari
Pholus       10.3° Cap
Cerere       17.8° Cap
Pallade       8.2° Sgr
Giunone       6.9° Vir
Vesta         2.4° Cnc
```

Per i TNO e gli altri corpi senza ID dedicato, usi il **numero del catalogo + `AST_OFFSET`** (10000):

```python
import libephemeris as ephem

# Eris (numero 136199)
ERIS = ephem.AST_OFFSET + 136199  # = 146199

jd = ephem.julday(2024, 4, 8, 12.0)
pos, _ = ephem.calc_ut(jd, ERIS, ephem.FLG_SPEED)
print(f"Eris: {pos[0]:.2f}°")
```

```
Eris: 24.74°
```

---

## 12.2 La catena di calcolo

I metodi disponibili dipendono da modalità e risorse selezionate. In `leb` forzato serve un gruppo di coefficienti coperto o un modello locale supportato: non si passa silenziosamente a un kernel JPL. Nei percorsi JPL nativo/automatici, queste sono le principali risorse e politiche da distinguere:

**1. Kernel SPK** — Un file JPL registrato fornisce una traiettoria sui segmenti effettivi. Contano centro, frame, copertura e incertezza della soluzione orbitale: il contenitore SPK da solo non garantisce un limite angolare fisso.

**2. Acquisizione SPK automatica** — I corpi mappati supportati possono essere acquisiti quando selezione automatica e politica di rete lo permettono entrambe. Un file assente in `sealed` non autorizza una connessione. In produzione prepara preferibilmente i dati prima del calcolo.

**3. Precisione strict** — `set_strict_precision(True)` rifiuta il fallback approssimato dove è richiesta una sorgente precisa supportata. Leggi l'errore dell'operazione invece di disabilitare il controllo per nascondere dati mancanti. Imposta `False` solo se l'approssimazione è adeguata alla tua applicazione.

**4. Propagazione N-body nativa** — Il motore Rust esegue le integrazioni supportate quando sono disponibili condizioni iniziali e file delle forze. Non serve un extra Python `[nbody]` per installare REBOUND/ASSIST. I nomi storici restano per compatibilità, ma l'integratore è indipendente. L'accuratezza dipende da condizioni iniziali, modello delle forze e intervallo, non dal solo possesso di elementi orbitali.

**5. Fallback kepleriano** — Un percorso approssimato supportato propaga gli elementi disponibili con le correzioni selezionate. Perturbazioni mancanti e incertezza degli elementi si accumulano nel tempo; gli incontri ravvicinati sono particolarmente delicati. Essere disponibile offline non dimostra che una traiettoria annuale sia sufficientemente accurata.

```python
import libephemeris as ephem

# Esempio online facoltativo: usciamo deliberatamente dagli esercizi LEB sealed.
ephem.set_calc_mode("skyfield")
ephem.set_network_policy("allow")
ephem.set_auto_spk_download(True)

# Ora calc_ut scaricherà l'SPK se necessario
jd = ephem.julday(2024, 4, 8, 12.0)
pos, _ = ephem.calc_ut(jd, ephem.CERES, ephem.FLG_SPEED)
print(f"Cerere (con SPK automatico): {pos[0]:.4f}°")
```

```
Cerere (con SPK automatico): 287.7759°
```

---

## 12.3 Kernel SPK: il gold standard

I kernel SPK (Spacecraft and Planet Kernel) sono file binari NASA che contengono traiettorie precise sotto forma di polinomi di Chebyshev. Sono lo stesso formato usato per i pianeti nelle efemeridi DE440.

### Scaricare e registrare un SPK

```python
import libephemeris as ephem

# Scarica l'SPK per Cerere e registralo
percorso = ephem.download_and_register_spk(
    "1;",             # identificativo per JPL Horizons
    ephem.CERES,   # ID del corpo nella libreria
    "2000-01-01",     # data inizio
    "2050-01-01",     # data fine
)

print(f"SPK scaricato: {percorso}")
```

```
SPK scaricato: ~/.libephemeris/spk/1_200001_205001.bsp
```

### Gestire gli SPK registrati

```python
import libephemeris as ephem

# Quali corpi hanno un SPK caricato?
bodies = ephem.list_spk_bodies()
for body_id, (path, naif_id) in bodies.items():
    print(f"Body {body_id}: NAIF={naif_id}, file={path}")

# Un corpo specifico ha l'SPK?
if ephem.is_spk_available_for_body(ephem.CERES):
    print("Cerere: SPK disponibile")
```

```

```

### Asteroidi "maggiori" con SPK preconfigurato

La libreria conosce i parametri di download SPK per 37+ corpi. Puoi assicurarti che un SPK sia disponibile con:

```python
import libephemeris as ephem

# Assicurati che l'SPK di Vesta sia disponibile
successo = ephem.ensure_major_asteroid_spk(ephem.VESTA)
if successo:
    print("SPK di Vesta pronto")

# Lista tutti i corpi con SPK scaricabile
for body_id, nome in ephem.list_major_asteroids():
    print(f"  {nome} (ID: {body_id})")
```

```
SPK di Vesta pronto
  Ceres (ID: 17)
  Pallas (ID: 18)
  Juno (ID: 19)
  Vesta (ID: 20)
  Chiron (ID: 15)
```

---

## 12.4 Cercare un asteroide per nome o numero

Se conosci il nome di un asteroide ma non il suo numero, la libreria può cercarlo — prima nel database locale, poi interrogando il servizio JPL SBDB (Small Body Database):

```python
import libephemeris as ephem

# Cerca per nome
numero = ephem.get_asteroid_number("Vesta")
print(f"Vesta = asteroide n. {numero}")  # 4

numero = ephem.get_asteroid_number("Apophis")
print(f"Apophis = asteroide n. {numero}")  # 99942
```

```
Vesta = asteroide n. 4
Apophis = asteroide n. 99942
```

### Calcolare per numero

Per calcolare la posizione di un asteroide qualsiasi dato il suo numero di catalogo:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Calcola la posizione di Eros (n. 433)
lon, lat, dist = ephem.calc_asteroid_by_number(433, jd)
print(f"Eros: lon={lon:.2f}°, lat={lat:.2f}°, dist={dist:.4f} UA")
```

```
Eros: lon=91.45°, lat=5.98°, dist=1.1666 UA
```

La funzione scarica automaticamente gli elementi orbitali da JPL SBDB se non sono già in cache. Puoi anche scaricarli esplicitamente:

```python
import libephemeris as ephem

# Scarica gli elementi orbitali di Apophis
elements = ephem.fetch_orbital_elements_from_sbdb(99942)
if elements:
    print(f"Nome: {elements.name}")
    print(f"Semi-asse maggiore: {elements.a:.4f} UA")
    print(f"Eccentricità: {elements.e:.6f}")
    print(f"Inclinazione: {elements.i:.4f}°")
```

> **Nota**: questa funzione richiede accesso a Internet per interrogare il servizio JPL SBDB (Small Body Database). I risultati vengono messi in cache locale per le richieste successive.

---

## 12.5 Il fallback kepleriano

Quando non è disponibile un kernel SPK, la libreria calcola le posizioni usando le **leggi di Keplero** con correzioni per le perturbazioni dei pianeti giganti.

Il metodo funziona così:

1. **Elementi orbitali osculanti**: sei numeri che descrivono l'ellisse dell'orbita a un'epoca di riferimento (semi-asse maggiore, eccentricità, inclinazione, nodo ascendente, argomento del perielio, anomalia media)

2. **Equazione di Keplero**: dato il tempo trascorso, calcola dove si trova il corpo sull'ellisse

3. **Perturbazioni secolari**: Giove, Saturno, Urano e Nettuno "spingono" lentamente l'orbita, facendola ruotare e cambiare forma. La libreria applica queste correzioni

4. **Modello di librazione**: per i "plutini" (corpi in risonanza 2:3 con Nettuno, come Ixion e Orcus), la libreria corregge per l'oscillazione dell'argomento di risonanza

La precisione tipica del fallback kepleriano dipende dal tempo trascorso dall'epoca degli elementi:

- **1 mese**: ~7 secondi d'arco — eccellente
- **1 anno**: ~2 minuti d'arco — buona per la maggior parte degli usi
- **10 anni**: ~30 minuti d'arco — accettabile per scopi generici
- **50 anni**: ~3.6° — solo orientativa

Per l'astrologia moderna (date dal 1900 al 2100), la precisione è generalmente sufficiente. Per lavori di ricerca o per date lontane dall'epoca degli elementi, usa sempre gli SPK.

---

## Riepilogo

In questo capitolo abbiamo imparato a lavorare con i corpi minori del Sistema Solare.

**Concetti chiave:**

- I **corpi minori** includono asteroidi, centauri, TNO e asteroidi near-Earth — migliaia di corpi oltre ai pianeti
- **L'instradamento dipende da modalità e politica**: gruppi LEB locali, SPK registrati, acquisizione permessa esplicitamente, strictness e propagazione nativa supportata sono scelte distinte.
- I **kernel SPK** sono file binari NASA con traiettorie precise — il gold standard per posizioni sub-arcsecondo
- La **strict precision** (default) solleva `SPKRequiredError` per i corpi mappati senza SPK, prevenendo perdite di precisione silenti. Disabilita con `set_strict_precision(False)`
- Il **fallback kepleriano** funziona senza Internet ma è raggiungibile solo se la strict precision è disabilitata o il corpo non è nella mappa SPK
- Per i corpi senza ID dedicato, usa `AST_OFFSET + numero_catalogo`

**Funzioni introdotte:**

- `calc_ut(jd, CERES, flag)` — calcola la posizione di un corpo minore con ID dedicato, proprio come per i pianeti
- `download_and_register_spk(body_id, jd_start, jd_end)` — scarica e registra un kernel SPK da JPL
- `ensure_major_asteroid_spk(body_id)` — si assicura che l'SPK sia disponibile, scaricandolo se necessario
- `list_major_asteroids()` — lista i corpi con download SPK automatico supportato
- `list_spk_bodies()` — mostra quali corpi hanno un SPK registrato
- `set_auto_spk_download(True)` — abilita il download automatico degli SPK
- `get_asteroid_number(nome)` — cerca il numero di catalogo di un asteroide per nome
- `calc_asteroid_by_number(numero, jd)` — calcola la posizione di qualsiasi asteroide dato il numero
- `fetch_orbital_elements_from_sbdb(numero)` — scarica gli elementi orbitali da JPL SBDB
