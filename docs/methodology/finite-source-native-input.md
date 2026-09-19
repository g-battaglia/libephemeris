# Native finite-source input bridge

`libephemeris._finite_source_input_model` is a private, disconnected
arithmetic bridge. It implements the reviewed bounded input interpretation
in `validation/golden/specs/occultation_contact_residuals/finite-source-post-return-input-model.md`
at validation commit `aae8dbf0cd5c292e9b2f5830b92d5f0975380ccd`
(SHA-256 `e9d3937b290ca2bf7c355b762891cb7de3d3d072063b28638793938aaa837c36`).
The operation chains are pinned to product commit
`1147e6a907b35d6d52613e56c1986c094a5897c9`.

The bridge accepts synthetic native Python binary64 words for the Moon and
Sun apparent equatorial-of-date AU positions, their separately converted
kilometre positions, the UT word and the fixed radius, axis and flattening
chain. It checks each of the six AU-to-kilometre products by all 64 bits,
including signed zero. It replays every fixed native constant operation
and checks its intermediate and result words. Malformed types, nonfinite
words, overflowed products or mismatched bits are invalid input, rather
than physical shadow classifications.

It converts each accepted kilometre centre component individually with
`Fraction.from_float`, then derives the exact rational span `u=p-B`. The
existing native subtraction remains available only as a diagnostic because
it can round differently. From the accepted native Earth axis `a` and
flattening `f`, it derives `b=a(1-f)`,
`A=diag(1/a²,1/a²,1/b²)` and `R=a` exactly. The independent finite-source
checker owns the ordered theorem-domain outcomes and strict certificate
checks. The bridge also exposes the six exact differences between its
native-kilometre interpretation and an alternative interpretation using
the IAU 2012 exact astronomical unit. These differences are representation
facts, not bounds on physical or ephemeris error.

The module obtains no astronomical state or data asset. Numeric consistency
cannot authenticate a successful call, provider, time branch, frame, file
or coefficient. Coherent changes to an AU word and its kilometre product
can pass this bridge; a separately reviewed call-local source receipt must
reject them against its owned result. No runtime producer, witness
generator, public eclipse route, temporal contact or physical accuracy
claim follows from this module.
