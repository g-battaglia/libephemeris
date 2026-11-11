# LibEphemeris bundled base data

This data-only crate supplies the exact published base-core LEB2 bytes to the
engine's optional `bundled-data` feature. It performs no calculations or runtime
I/O. Keeping it separate lets both registry archives remain below their size
limit without changing the data, removing a feature or adding a download.

The engine retains static, zero-copy ownership and validates the container.
Actual body/channel coverage comes from the image, not the package name.
The published data-v4 image is 10,232,271 bytes, SHA-256
`b891de02f82b68346046da67389fd3713aa81ab299287c632cb20cc396397964`.
Its stored core interval is JD 2396757.5–2506332.5; per-channel bounds remain
authoritative.
It is distinct from other independently identified local datasets.

Applications normally depend on `libephemeris-core` with `bundled-data`, not
this crate directly. Project distribution terms are AGPL-3.0; complete license
and upstream notices accompany the source package. See the canonical
[source documentation](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/docs/sources.md).
The stable data-crate version changes package identity, not the data-v4 image;
crates.io publication is separate from Python wheel publication.
For private licensing inquiries, see the [project contact](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0/README.md#license);
requests are evaluated case by case.
