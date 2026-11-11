# Migrating from v3 to v4

See the [v4 release notes](RELEASE_NOTES.md) for the Rust rewrite,
independent implementation and published data-v4 inventory.

The import remains `import libephemeris`. Calculations now run in Rust through
PyO3; Python names, argument order, flags and ordinary result layouts remain
familiar. Python requires 3.12+. Use [installation](docs/installation.md) to
select 4.0.0 explicitly after publication, or build the prepared source locally.

| Previous assumption | v4 behavior / action |
| --- | --- |
| Skyfield/PyERFA Python calculation backend | Native Rust; `skyfield` remains a mode name for JPL/SPK |
| REBOUND/ASSIST dependencies | Independent Rust integration under compatibility names; external force files remain separate inputs |
| Python loader, timescale or LEB/SPK object internals | Retained native handles; use their public methods, not dependency-specific attributes |
| Arbitrary Python numeric-object arithmetic | Typed numeric/string/sequence conversion; values must fit declared widths and shapes |
| Dataclass/enum reflection on integration records | Native named values, not reconstructed Python class internals |
| Dependency-specific exceptions | Documented classes and useful fields, not dependency exception graphs or warning order |
| Mutating environment after configuration capture | Call `reload_configuration()` to reacquire settings; dotenv loading does not mutate `os.environ` |
| `close()` as per-request cleanup | Resource lifecycle operation; use independent observing contexts and bounded workers |
| Async cancellation interrupts calculation | Calls are synchronous even when the GIL is released |
| Filename establishes reviewed data identity | Inventory review requires an actual immutable-image SHA-256 match |
| Bit-for-bit equality with v3 | Compare physical quantities under aligned data/models; independently justified corrections remain |

`calc` takes TT; `calc_ut` takes UT1. `utc_to_jd` returns TT first, UT1 second.
Module-level `calc_pctr` takes TT, whereas the context method takes UT1.
House cusps remain zero-indexed. See [Python conventions](docs/python.md).

Explicit missing/corrupt LEB files are errors in forced LEB mode. In `auto`,
an unavailable optional image permits normal provider routing, not an arbitrary
replacement LEB. Resource-budget refusals and invalid reader settings remain
errors. Mutable catalog snapshots do not edit the compiled tables.

Spica's static Johnson V magnitude is corrected from 1.04 to Hipparcos H5 0.98.
Other model comparisons must align frame, clock, observer and actual input data;
see [precision and limitations](docs/limitations.md). The compatibility baseline
is v3.2.1 (`37ffb0555b7a9f5083b4f1b6dc1a0f10feeb85fa`), not a new audit of every
later v3 release.

## Rust callers

The package is `libephemeris-core`, imported as `libephemeris_core`.
`api` is the primary typed interface: pure functions and `Session` methods
return named records and `api::Error`. For example, `julday(...)?` returns a
record with `.value`, not a bare float.

The former broad SDK Session is now `native::Session`. Configuration and
acquisition remain in `context`/`config`; attach resources before building an
explicit primary Session. `nbody` enables native integration (`native-nbody`
is an alias); `python-internals` is not an application API. The C ABI/ctypes
package is retired. See [Rust](docs/rust.md) and the
[consumer contract](core/tests/common_api_contract.rs).
