# LibEphemeris per principianti

Un libro per sviluppatori, astrologi e osservatori curiosi che vogliono capire
il cielo prima di calcolarlo. Serve un minimo di Python, non di astronomia.

Parti dalla prefazione e leggi in ordine i capitoli 1–3. Se conosci già
l'astronomia posizionale, puoi iniziare dal capitolo 5; il ricettario propone
applicazioni pratiche complete. La breve [documentazione tecnica](../../docs/README.md)
serve per consultazione, non sostituisce le spiegazioni di questo libro.

## Capitoli

0. [Prefazione e primo calcolo](00-prefazione.md)
1. [Il cielo visto dalla Terra](01-cielo-dalla-terra.md)
2. [Misurare il tempo](02-misurare-il-tempo.md)
3. [Coordinate celesti](03-coordinate-celesti.md)
4. [Le efemeridi](04-efemeridi.md)
5. [Posizione dei pianeti](05-posizione-pianeti.md)
6. [La Luna](06-luna.md)
7. [Case astrologiche](07-case-astrologiche.md)
8. [Stelle fisse](08-stelle-fisse.md)
9. [Eclissi e occultazioni](09-eclissi.md)
10. [Alba, tramonto e visibilità](10-alba-tramonto-visibilita.md)
11. [Lo zodiaco siderale](11-zodiaco-siderale.md)
12. [Corpi minori](12-corpi-minori.md)
13. [Pianeti ipotetici e parti arabe](13-pianeti-ipotetici-parti-arabe.md)
14. [Precisione, prestazioni e dati](14-precisione-performance.md)
15. [Ricettario: dodici applicazioni pratiche](15-ricettario.md)

## Eseguire gli esempi

Installa la [beta nativa](../../docs/installation.md), poi esegui questa
configurazione prima di ogni esempio completo. I frammenti brevi senza import
proseguono l'esempio precedente e ne riutilizzano le variabili. Così gli esercizi
non dipendono dalle impostazioni abituali della tua applicazione:

```python
import libephemeris as ephem

ephem.set_calc_mode("leb")
ephem.set_precision_tier("base")
ephem.set_network_policy("sealed")
ephem.set_iers_delta_t_enabled(False)
ephem.set_auto_spk_download(False)
```

Il primo calcolo e quasi tutti gli esempi planetari funzionano con il solo
nucleo incluso. Salvo conversione esplicita da data civile UTC, negli esercizi
**UT significa UT1**. I risultati sono esempi arrotondati, non limiti d'errore
astronomico. La trattazione dell'astrologia tradizionale spiega convenzioni,
non prove scientifiche di previsioni sugli eventi personali.

### Dati aggiuntivi ed esempi online

Apsidi lunari, corpi minori e il tema natale con Chirone possono richiedere
altri gruppi base. Il confronto JPL/LEB del capitolo 14 richiede un kernel JPL
locale. Prepara esplicitamente questi input, poi mantieni la rete chiusa:

```sh
libephemeris download base
```

Puoi scegliere una directory dedicata con `LIBEPHEMERIS_DATA_DIR`. Controlla
la copertura dei corpi con `libephemeris status -vv`: un tier non contiene
necessariamente ogni asteroide. La [guida ai dati](../../docs/data.md) spiega
file e configurazione.

I blocchi che scaricano dati IERS/SPK/LEB, interrogano SBDB o selezionano
Horizons sono **operazioni online facoltative**, non esercizi offline.
Possono scrivere cache o acquisire file grandi: leggi i prerequisiti prima di
eseguirli. Per un'operazione di rete cambia deliberatamente la politica con
`ephem.set_network_policy("allow")`; poi riapplica la configurazione offline.
Gli esempi di selezione del file richiedono `LIBEPHEMERIS_DATA_DIR/leb/base_core.leb2`.
Le illustrazioni delle risposte online non sono trascrizioni verificate localmente.

Per i controlli di manutenzione, `tools/release/check_manual.py` censisce tutti
gli esempi, esegue quelli offline/con dati locali senza acquisizioni e riporta
a parte i blocchi online facoltativi. Un esempio saltato non è una verifica
online riuscita. Vedi [contribuire](../../CONTRIBUTING.md#documentation).

[English](../en/README.md) · [Provenienza e licenza](../README.md#provenance-and-license)

```{toctree}
:hidden:
:maxdepth: 1

00-prefazione
01-cielo-dalla-terra
02-misurare-il-tempo
03-coordinate-celesti
04-efemeridi
05-posizione-pianeti
06-luna
07-case-astrologiche
08-stelle-fisse
09-eclissi
10-alba-tramonto-visibilita
11-zodiaco-siderale
12-corpi-minori
13-pianeti-ipotetici-parti-arabe
14-precisione-performance
15-ricettario
```
