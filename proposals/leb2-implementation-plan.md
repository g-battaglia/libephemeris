# LEB2: piano di implementazione

**Stato:** implementato; piano storico sanificato dopo la review clean-room.

Le tabelle originali di dimensioni, rapporti di compressione, residui, nomi di
asset pubblicati e gruppi ipotetici sono state rimosse. Alcuni file storici
includevano modelli che non soddisfano più il requisito di provenienza e non
devono essere recuperati, convertiti o usati come riferimento.

## Confine attuale

- Il wheel contiene soltanto il `base_core.leb2` revisionato.
- Ogni altro LEB1/LEB2 viene generato localmente da kernel NASA JPL e modelli
  indipendenti con fonte citata.
- Gli asset LEB pubblicati in passato sono ritirati.
- Gli output dell'oracle di compatibilità non entrano in coefficienti, fixture,
  soglie, report o metadati distribuiti.
- Tra i corpi ipotetici integrati è supportato soltanto Harrington (ID 50),
  derivato dall'articolo di Harrington del 1988.

## Architettura implementata

1. Il formato mantiene header e directory di sezioni versionati.
2. L'indice per corpo dichiara intervallo temporale, coordinate e blocco dati.
3. I coefficienti sono compressi per corpo o per chunk temporale.
4. Il reader valida limiti, offset, lunghezze e checksum prima dell'uso.
5. La decompressione è lazy e la cache mutabile è protetta da lock.
6. `open_leb()` rileva LEB1/LEB2 e offre la stessa interfaccia al runtime.
7. Conversione e verifica accettano solo input generati localmente con
   provenienza completa.

## Verifica

- Round-trip e test strutturali usano fixture sintetiche.
- Gli error bound sono obiettivi numerici dichiarati, non soglie adattate
  all'output dell'oracle.
- I file locali sono confrontati con il percorso JPL/IAU indipendente.
- I gate di provenienza e packaging devono passare prima di considerare un
  artefatto distribuibile.

Per il workflow corrente vedere `docs/leb/guide.md`; per le decisioni sul
formato vedere `proposals/leb2-compressed-format.md`.
