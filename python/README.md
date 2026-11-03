# LibEphemeris Python stack

This directory owns the thin PyO3 binding crate, import-only runtime modules,
type stubs, installed-package tests and Python integration examples.
Astronomical calculation and configuration remain in the Rust engine.

The shared distribution metadata and Maturin entry point stay in the repository
root so ordinary source installs and source archives remain conventional. The
binding's `Cargo.toml` depends on `../core` and is not separately published to
crates.io. The root Cargo workspace owns common build profiles and the lockfile.

From the workspace root, with Python 3.12+:

```sh
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python maturin==1.9.4
.venv/bin/maturin build --release --locked --interpreter .venv/bin/python --out dist
uv pip install --python .venv/bin/python dist/*.whl
.venv/bin/python -m pytest python/tests
```

Python version **4.0.0b1** is a pre-release. See the shared
[beginner's manual](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b1/beginners-manual/en/README.md),
[migration guide](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b1/MIGRATION.md) and [release instructions](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b1/CONTRIBUTING.md).
For private licensing inquiries, see the [project contact](https://github.com/g-battaglia/libephemeris/blob/python-v4.0.0b1/README.md#license);
requests are evaluated case by case.
