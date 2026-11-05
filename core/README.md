# LibEphemeris Rust core

`libephemeris-core` provides native astronomical calculation. Public computation
lives in `libephemeris_core::api`; native resource configuration and advanced
integration live in `context`, `config` and `native`.

Version **1.0.0-beta.3** is a pre-release. Once available on the registry:

```sh
cargo add libephemeris-core@1.0.0-beta.3 --features bundled-data
```

From a source checkout's workspace root:

```sh
cargo check --locked -p libephemeris-core
cargo test --locked -p libephemeris-core
cargo build --release --locked -p libephemeris-core --features bundled-data,nbody
```

Default features are empty and require no Python. `bundled-data` statically
embeds the unchanged base image through the project-owned data-only package
`libephemeris-core-data`; `nbody` enables independent Rust propagation.
`native-nbody` is an alias. `python-internals` is an unstable binding bridge,
not an application interface. The declared MSRV is Rust 1.86.

The root workspace retains common build profiles and the lockfile. This
subdirectory owns the crate manifest, source, Rust integration tests/examples,
and actual upstream dependency notices. `data/` is the optional dataset crate,
not another calculation engine or shared API facade.

See the [primary API guide](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b3/docs/rust.md),
[getting started](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b3/docs/installation.md)
and [source register](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b3/docs/sources.md).
Project distribution terms are AGPL-3.0; complete texts accompany the archive.
For private licensing inquiries, see the [project contact](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b3/README.md#license);
requests are evaluated case by case.
