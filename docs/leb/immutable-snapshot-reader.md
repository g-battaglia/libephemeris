# Private immutable LEB2-v2 snapshot reader

`libephemeris._leb2_snapshot` is a project-authored, private byte-ownership
layer over the ordinary LEB2 parser and evaluator. It implements the bounded
contract at
`validation/golden/specs/occultation_contact_residuals/leb2-immutable-snapshot-reader-contract.md`
in validation commit `ec89c41ecd900fa92bbcbacf05784d448bc9e5df`
(SHA-256 `811f128845e87c47bb6ce6a44268873a496b988b2b7e8b94f316de891b985c88`).
The ordinary `LEB2Reader(path)` remains the same mmap-backed reader.

The production factory admits only the fixed `base` or `medium` LEB2-v2
core identity. Its expected byte count and SHA-256 come from the reviewed
allowlist, independently of the path and the bytes being read. It opens one
descriptor, copies no more than the expected count plus one byte into an
immutable `bytes` object, closes the descriptor, checks the exact count and
digest, and only then invokes the shared LEB2 parser. The base record is
10,232,283 bytes with SHA-256
`5d708bdbe3e799e0802ba575984e57a3c5e44720dbfa1b4a01cf826640e0cb82`;
the medium record is 37,276,175 bytes with SHA-256
`4d88ec9a79add7e3af9e75ac3ceabe5462a4af440447578a5ccba69a3e0a55b6`.
The digest check relies on the standard SHA-256 collision-resistance
assumption. A separate synthetic-fixture factory exercises the same copy,
hash and parse path but cannot issue a production-asset identity.

All subsequent metadata reads, chunk decompression and evaluations consume
the accepted immutable object. The snapshot has no mmap prefetch or cache
reclamation operation; `warm()` and `cool()` therefore perform no OS memory
advice. `close()` clears reader-owned caches and drops its reference to the
bytes; later evaluation and metadata queries fail with a typed closed-reader
error. A corruption detected during lazy evaluation closes the private
reader and propagates a typed corruption error. The first design is
call-local and does not permit concurrent close and evaluation.

For identical body and epoch arguments, this reader is intended to return
the same native Python float words as the ordinary LEB2 reader when the
ordinary reader's backing file stays byte-identical through its parse,
cache population and evaluation. Bitwise agreement under that condition
is an evaluator-equivalence result. The snapshot's verified asset identity
alone does not attest which tiered child, time branch, frame cache,
deflector, nutation or state words an ordinary `calc_ut()` call used. A
later private source computation must carry those dependencies into its
result before any B-16 astronomical source receipt can be claimed.
