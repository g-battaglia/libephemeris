# Prefazione e come usare questo manuale

## Cosa imparerai

In queste poche pagine introduttive scoprirai a chi è rivolto questo manuale, come è organizzato, e farai il tuo primo calcolo astronomico in tre righe di Python.

---

## A chi è rivolto

Questo manuale è per tre tipi di lettori:

- **Sviluppatori** che devono integrare calcoli astronomici o astrologici in un'applicazione e vogliono capire cosa c'è dietro i numeri, non solo copiare una riga di codice.
- **Astrologi** che vogliono passare dal software commerciale a uno strumento aperto e verificabile, senza rinunciare alla precisione.
- **Curiosi di astronomia** che si chiedono come si calcola dove si trova un pianeta, quando sorge il Sole o perché un'eclissi è visibile solo in certi luoghi.

Non serve sapere nulla di astronomia. Serve solo un po' di Python — il livello di chi sa cos'è una funzione, un dizionario e un ciclo `for`.

## Cosa imparerai leggendo questo manuale

Questo non è un manuale di riferimento delle API. È un libro che spiega i **concetti** che stanno dietro ai calcoli. La differenza è importante: un riferimento API ti dice "la funzione `calc_ut` accetta tre parametri"; questo manuale ti spiega **perché** quei tre parametri esistono, cosa significano fisicamente, e quali errori potresti fare se non li capisci.

Alla fine del manuale saprai:

- Come funziona la sfera celeste e perché i pianeti si muovono come si muovono
- Perché il tempo in astronomia è complicato e come evitare le trappole
- Cosa sono le coordinate eclittiche, equatoriali e orizzontali
- Come calcolare la posizione di qualsiasi corpo celeste con precisione sub-arcsecondo
- Come funzionano le case astrologiche e perché esistono 20 sistemi diversi
- Come trovare eclissi, albe, tramonti e fenomeni di visibilità
- Come ottimizzare le performance per calcoli su larga scala

## Come è organizzato

Ogni capitolo segue la stessa struttura:

- **📐 Teoria**: il concetto astronomico o astrologico, spiegato con analogie e senza formule inutili
- **🌍 Vita reale**: un esempio concreto che collega il concetto a qualcosa di osservabile o familiare
- **💻 Codice**: un esempio funzionante con LibEphemeris, che puoi copiare e incollare

I capitoli sono pensati per essere letti in ordine. I primi tre costruiscono le basi (il cielo, il tempo, le coordinate) e tutto il resto si appoggia su di essi. Se hai già familiarità con l'astronomia posizionale, puoi saltare direttamente al Capitolo 5 e tornare indietro quando serve.

Il Capitolo 15 — il Ricettario — è pensato per chi ha fretta: ricette copia-incolla per i calcoli più comuni, con spiegazioni brevi.

## Prerequisiti

- **Python 3.9** o superiore
- **Nessuna conoscenza astronomica** — questo manuale parte da zero
- **Nessun altro pacchetto** — LibEphemeris si installa con un singolo comando e non ha dipendenze pesanti

## Installazione

```bash
pip install libephemeris
```

Al primo utilizzo, la libreria scarica automaticamente i file di efemeridi necessari (il file DE440, circa 114 MB). Puoi anche forzare il download in anticipo:

```python
import libephemeris as ephem
ephem.download_for_tier("medium")
```

I tre livelli di precisione — `base`, `medium`, `extended` — si differenziano per l'intervallo temporale coperto e la dimensione dei file. Il livello `medium` (default) copre gli anni 1550–2650 e va bene per la stragrande maggioranza degli usi. Ne parleremo nel dettaglio al Capitolo 4.

## Il tuo primo calcolo

Calcoliamo dove si trovava il Sole durante l'eclissi solare totale dell'8 aprile 2024:

```python
import libephemeris as ephem

# Convertiamo la data in Giorno Giuliano (mezzogiorno UT)
jd = ephem.julday(2024, 4, 8, 12.0)

# Calcoliamo la posizione del Sole
# Il terzo argomento (0) indica: coordinate eclittiche, senza opzioni extra
pos, flag = ephem.calc_ut(jd, ephem.SUN, 0)

longitudine = pos[0]  # gradi lungo l'eclittica (0–360)
print(f"Sole a {longitudine:.4f}° di longitudine eclittica")
```

```
Sole a 19.1404° di longitudine eclittica
```

Che cosa è successo?

1. `julday(2024, 4, 8, 12.0)` ha convertito la data "8 aprile 2024, ore 12:00 UT" in un **Giorno Giuliano** — un numero unico che identifica quell'istante. Impareremo tutto sui Giorni Giuliani nel Capitolo 2.

2. `calc_ut(jd, SE_SUN, 0)` ha calcolato la posizione del Sole per quell'istante. Il risultato è una tupla di 6 numeri: longitudine, latitudine, distanza, e le rispettive velocità giornaliere. Il terzo parametro (`0`) dice alla libreria di usare le impostazioni standard. I flag di calcolo sono trattati nel Capitolo 5.

3. La longitudine `19.15°` ci dice che il Sole si trovava a circa 19 gradi del segno dell'Ariete (il primo segno va da 0° a 30°). Capiremo il perché nel Capitolo 1.

Facciamo un passo in più — mostriamo la posizione in formato zodiacale:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)
pos, flag = ephem.calc_ut(jd, ephem.SUN, 0)

segni = [
    "Ariete", "Toro", "Gemelli", "Cancro",
    "Leone", "Vergine", "Bilancia", "Scorpione",
    "Sagittario", "Capricorno", "Acquario", "Pesci",
]

# split_deg scompone i gradi decimali in gradi, minuti, secondi e segno
deg, min, sec, secfr, segno = ephem.split_deg(
    pos[0], ephem.SPLIT_DEG_ZODIACAL | ephem.SPLIT_DEG_ROUND_SEC
)

print(f"Sole a {deg}° {min}' {sec}\" {segni[segno]}")
```

```
Sole a 19° 8' 26" Ariete
```

La funzione `split_deg` prende un valore in gradi decimali e lo scompone in gradi, minuti d'arco, secondi d'arco, frazione di secondo e numero del segno zodiacale. È lo stesso formato che trovi sulle riviste di astrologia o negli almanacchi astronomici.

## Convenzioni usate in questo manuale

In tutto il manuale:

- **L'import è sempre lo stesso**: `import libephemeris as ephem`. Tutte le funzioni e le costanti sono accessibili tramite `ephem.nome`.
- **Gli esempi sono completi**: ogni blocco di codice può essere copiato e incollato in un file `.py` e funziona. Niente frammenti parziali.
- **Le unità sono sempre indicate**: gradi (°), minuti d'arco ('), secondi d'arco ("), unità astronomiche (UA), chilometri (km). Se vedi un numero senza unità, è un errore.
- **Le date degli esempi sono eventi reali**: eclissi, congiunzioni, equinozi effettivamente avvenuti. Puoi verificare i risultati con qualsiasi altro software astronomico.

## Struttura del manuale

| #  | Capitolo | Argomento |
|----|----------|-----------|
| 1  | Il cielo visto dalla Terra | La sfera celeste, l'eclittica, lo zodiaco, l'orizzonte |
| 2  | Misurare il tempo | Giorno Giuliano, UT, TT, Delta-T, tempo siderale |
| 3  | Coordinate celesti | Eclittiche, equatoriali, orizzontali, conversioni |
| 4  | Le efemeridi | JPL DE440, livelli di precisione, geocentrico vs eliocentrico |
| 5  | Posizione dei pianeti | `calc_ut`, flag, velocità, retrogradazione, fenomeni |
| 6  | La Luna | Nodi, apogeo, perigeo, Lilith, attraversamenti |
| 7  | Case astrologiche | Sistemi di case, Ascendente, MC, latitudini estreme |
| 8  | Stelle fisse | Catalogo, moto proprio, stelle in astrologia |
| 9  | Eclissi | Solari, lunari, cicli di Saros, occultazioni |
| 10 | Alba e tramonto | Levata, transiti, crepuscoli, rifrazione, visibilità eliaca |
| 11 | Zodiaco siderale | Tropicale vs siderale, ayanamsha, calcoli vedici |
| 12 | Corpi minori | Asteroidi, centauri, TNO, kernel SPK, fallback kepleriano |
| 13 | Pianeti ipotetici | Uraniani, Transpluto, parti arabe |
| 14 | Precisione e performance | Confronti, modalità LEB, configurazione |
| 15 | Ricettario | Tema natale, transiti, eclissi, rivoluzioni solari e altro |

---

## Riepilogo

- LibEphemeris è una libreria Python per calcoli astronomici e astrologici di alta precisione
- Si installa con `pip install libephemeris` e usa i dati NASA JPL DE440
- La funzione base è `calc_ut(jd, corpo, flag)`: dato un istante e un corpo celeste, restituisce la posizione
- `julday(anno, mese, giorno, ore)` converte una data in Giorno Giuliano
- `split_deg(gradi, flag)` formatta i gradi decimali in gradi, minuti, secondi e segno zodiacale
- Ogni capitolo di questo manuale combina teoria, esempi di vita reale e codice funzionante
