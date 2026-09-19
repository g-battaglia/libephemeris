# Private Sun/Moon state diagnostic

`libephemeris._private_state_reducer` implements the reviewed bounded contract at
`validation/golden/specs/occultation_contact_residuals/private-sun-moon-state-reducer-contract.md`
(SHA-256 `ad2d7ae454f90ab7a2378e6c2d93aca4aecbf14566bae254f5fde4f66fb09299`).
It is disconnected from public eclipse and `calc_ut()` routing. Its results are
**diagnostic-only**. They are not source receipts, ordinary-call provenance,
contact geometry, or a whole-domain equivalence claim.

The private call owns fresh immutable base and medium LEB2-v2 snapshots and a
fresh three-asset default-time evaluator. It accepts only a native binary64
UT Julian date, the two exact equatorial Cartesian state request words, and
five asset locations. It evaluates private Delta-T once, forms TT with the
same ordered native addition as `fast_calc_ut()`, then calls the existing
`_fast_calc_core()` separately for Moon and Sun. Shared code also performs
Chebyshev decompression and evaluation; a private zstd decoder keeps this call
off the ordinary global decoder. The public decoder and mmap path retain
their existing behavior.

A private reader chooses base before medium independently for every target,
Earth, deflector, retarded split epoch, and nutation read. Each body record
retains exact time and result bits, its own ordered source history with both
candidate decisions, and control events for light-time, deflection,
aberration, and frame selection. The private frame computes selected LEB
nutation through the existing Vondrak/ERFA model without touching the
ordinary frame, precession, or mean-obliquity caches. A call-local provider
also prevents an ordinary `state.close()` generation change from sending
private frame work to Skyfield fallback. Private failures remain typed and
terminal, including corrupt deflector reads that the ordinary code could
otherwise skip as `ValueError`. The currently unresolved light-time early exit
and zero-geocentric deflection branches block with the private unsupported
outcome; they cannot yield a successful diagnostic pair.

Focused tests compare four fixed UT samples and both requests bit for bit
against a freshly controlled ordinary sealed-LEB path. Tests also exercise
mixed-tier source selection, two-part retarded epochs, asset changes before
and after admission, selected-child corruption, cache isolation, thread-local
restoration, and event/record alteration. These are finite regression samples;
they do not prove the complete dependency graph or a public call's source.
No astronomical vectors or reference API outputs are stored here or in tests.

The direct version guards cover the reviewed CPython, Skyfield, NumPy, pyerfa,
python-zstandard and libzstd versions, plus the earlier time component's
three direct Skyfield source-file hashes. A source-bearing claim still
requires independent review and exact binding of the *complete* transitive
numerical code, native binary build/ABI, and runtime constant inventory. In
particular, the deflection chain's `_DEFLECTORS` masses and `_GS`, `_C_MS`,
`_AU_M`, and `C_LIGHT_AU_DAY` values must be bound, along with relevant
branches and other loaded constants. The current control events show selected
body IDs and branch facts but cannot prove those constant values. A later
receipt schema must also specify stable encoding and independently verify
each result/dependency graph; a later same-call connector must authenticate
the ordinary calculation and chosen time branch.
