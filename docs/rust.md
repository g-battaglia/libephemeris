# Rust API

Use `libephemeris_core::api` for the primary typed interface: pure functions
and configured `Session` methods return named records and `api::Error`.
Data setup belongs to `context` and `config`; `native` exposes the advanced SDK.

Build the reference locally:

```sh
cargo doc --locked -p libephemeris-core --no-deps --open
```

## First calculation

This example uses an **existing** local LEB file. Substitute its actual path;
a configured pathname does not open a file or download missing data.
See [data setup](data.md) to provision files.

```rust
use libephemeris_core::{
    api,
    context::{EphemerisContext, NetworkPolicy, ProviderMode},
    leb::LebFile,
};
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut context = EphemerisContext::new()
        .with_mode(ProviderMode::SealedLeb)
        .with_network_policy(Some(NetworkPolicy::Sealed))
        .with_iers_delta_t(Some(false));
    context.attach_leb_file(LebFile::open(Path::new("/path/to/base_core.leb2"))?);
    let session = context.primary_session_builder().build()?;
    let jd_ut1 = api::julday(2000, 1, 1, 12.0, api::Calendar::Gregorian)?;
    let sun = session.calc_ut(jd_ut1.value, api::SUN, api::FLG_SPEED)?;
    println!("Longitude: {} degrees", sun.values[0]);
    Ok(())
}
```

`julday` labels the supplied calendar time; it does not convert UTC to UT1.
Use the public UTC conversion when starting with civil fields. Position units
and frame flags follow the [shared conventions](python.md#time-and-coordinates).
Keep two-part epochs separate when retaining small offsets matters.

## Sessions and errors

A Session retains immutable settings and shared resources. Attach the native
resources before building an explicit session. For application configuration
and acquisition, use `config::Runtime`; `Runtime::primary_session()` exposes
its configured snapshot through the primary API.

Named results are not Python tuples: read `Position.values`, `returned_flags`
and diagnostics, or `JulianDay.value`. Handle failures with `Result` and `?`
rather than panicking on application input.

`api::Error::kind()` gives a broad class. `operation()`, `native_status()`,
`native_scalar()` and `diagnostic()` retain available details. Import
`std::error::Error` to inspect `source()`. Optional diagnostics remain `None`
when absent; do not infer a status from error text. Advanced SDK failures use
`errors::Error`.

## Features and concurrency

Default Cargo features are empty. `bundled-data` embeds the immutable base
image from the data-only crate; `nbody` enables native integration.
`native-nbody` is a compatibility alias; `python-internals` is not an
application API.

Calls are synchronous. Share retained resources and bound the number of native
jobs; cancelling an awaiting task does not interrupt a solver. The
[worker example](../core/examples/bounded_native_workers.rs) demonstrates this
pattern. Indexed LEB access can reduce coefficient memory; see [data](data.md).

For migration from the former large `api::Session` SDK, use `native::Session`
where appropriate. The [consumer contract](../core/tests/common_api_contract.rs)
checks the primary exported interface. Build and test instructions are in
[contributing](../CONTRIBUTING.md).
