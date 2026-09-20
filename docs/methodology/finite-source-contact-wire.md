# Private finite-source contact request bytes

`libephemeris/_finite_source_contact_wire.py` is an in-memory codec for one
exact, fixed-state request to the disconnected B-16 contact generator. Its
source is the project-authored audit candidate
`B16_CONTACT_WIRE_ENCODING_CANDIDATE.md` (SHA-256
`3b22f659e41868e8b09862bb605a77a344c0ea4bb3d456f1e9b8be968b9109b4`),
reviewed in `B16_CONTACT_WIRE_ENCODING_REVIEW_HIGH.md` (SHA-256
`e1cb80d805e425a4d3ca3f2563f298eea15f259060c001a7becaeaf40f41bb07`).
That review passed a **non-normative design candidate only**. This module is
an implementation candidate for its request half, not a production protocol
approval or an independent implementation review.

## Encoded request

The frame is the five ASCII bytes `LEFC1`, a one-byte kind (`1`), an
eight-byte unsigned big-endian body length, then exactly that many body bytes.
The body has three 32-byte SHA-256 identities in order (wire grammar,
generator contract, resource policy), a one-byte target (`1` penumbra, `2`
umbra, `3` antumbra), then 18 rational values in order:

```text
B[0:3], p[0:3], r, m, A[0:3][0:3] row-major, R
```

Each rational is a signed numerator followed by a positive denominator. A
signed integer is a one-byte sign (`0` positive or zero, `1` negative), a
four-byte unsigned big-endian magnitude length, and that many big-endian
magnitude bytes. Zero has sign `0` and zero magnitude length. A nonzero
magnitude has no leading zero byte. Fractions must be reduced; the only zero
rational is `0/1`. The decoder rejects extra bytes, noncanonical values,
unrecognized tags, wrong external identities, and inconsistent length fields.

`_ContactWireLimits` supplies a total-frame byte ceiling and an integer-bit
ceiling per call. The encoder checks the bit ceiling before making a magnitude
byte string. The decoder checks declared magnitude size, leading byte and
actual bit length before converting magnitude bytes into a Python integer.
It then reconstructs `_ValidatedGeometry` from the six base fields, so cached
axis, squared axis and inverse metric cannot be transmitted as authority.
Wire errors are distinct from invalid geometry and the three existing
out-of-domain outcomes. The SHA-256 of the whole frame is a request association
identifier, not authentication of astronomical states, code or a source
receipt.

## Boundary and unresolved work

The decoder receives a complete Python `bytes` object. It checks length
before parsing, but cannot prevent an upstream transport from allocating an
oversized body. Neither the limits' numeric values nor a supported-platform
whole-call deadline and termination policy have been approved. Python and
native exact arithmetic, domain reconstruction, the contact generator and
independent verifier are not confined by this codec. No response grammar,
worker, source receipt, astronomical state acquisition, temporal contact,
public route or physical result is implemented here. The arbitrary digest
arguments are equality checks only; no normative wire or resource-policy hash
is selected by this module.

Synthetic fixtures exercise the byte grammar, four exact geometries and
rejection paths without consulting or persisting a reference API output. A
production boundary still needs independently reviewed request and response
semantics, a bounded transport reader, concrete resource policy, semantic
replay, provenance, and complete preselected scientific validation. B-16
remains blocked.
