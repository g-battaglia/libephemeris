# Installation

Python requires **3.12+**. This stable preparation uses Python `4.0.0` and Rust
core/data `1.0.0`. These instructions do not imply that the packages have already
been published to PyPI or crates.io; those are separate maintainer actions.

## Python wheel

After publication, select v4 explicitly:

```sh
python -m pip install "libephemeris==4.0.0"
```

Before publication, install a locally built wheel using the instructions below.
An unqualified stable install may select v3 instead. Check the imported package:

```python
import libephemeris
print(libephemeris.__version__)
print(libephemeris.__file__)
```

A wheel contains the native engine and modern base coefficients; NumPy is its
Python runtime dependency. It does not include all JPL kernels or data tiers.
The release wheel matrix covers macOS ARM64, manylinux2014 x86_64/aarch64 and
Windows x86_64. Other platforms, including Intel macOS, require a source build;
no native wheel is claimed for them. Source installs need the toolchain below;
an ABI3 wheel does not promise free-threaded Python support.

## First Python calculation

This works with the bundled core image, without additional data or downloads:

```python
import libephemeris as ephe

ephe.set_calc_mode("leb")
ephe.set_precision_tier("base")
ephe.set_network_policy("sealed")
ephe.set_iers_delta_t_enabled(False)

jd_tt, jd_ut1 = ephe.utc_to_jd(2024, 4, 8, 12, 0, 0.0)
values, flags = ephe.calc_ut(jd_ut1, ephe.MOON, ephe.FLG_SPEED)
print("Moon longitude:", values[0], "degrees")
cusps, angles = ephe.houses(jd_ut1, 41.9, 12.5, b"P")
print("Ascendant:", angles[0], "degrees")
```

For explanations, use the [beginner's manual](../beginners-manual/en/README.md).
For other datasets and explicit local files, see [data setup](data.md).

## Rust

After registry publication:

```toml
[dependencies]
libephemeris-core = { version = "1.0.0", features = ["bundled-data"] }
```

Before publication, use `path = "/path/to/libephemeris/core"` instead of the
registry version. The crate imports as `libephemeris_core`. Default features
are empty; `bundled-data` embeds the base image, while `nbody` enables native
integration. See the [Rust example](rust.md#first-calculation) for resource setup.

## Build from source

Install Rust (edition 2024, declared MSRV **1.86**), Python 3.12+, `uv`, and a C
compiler for Zstandard. From the repository root:

```sh
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python maturin==1.9.4
.venv/bin/maturin build --release --locked --interpreter .venv/bin/python --out dist
uv pip install --python .venv/bin/python dist/libephemeris-4.0.0-*.whl
cargo build --locked -p libephemeris-core --release
```

Install the wheel matching your platform if `dist/` contains multiple builds.
The Rust-only build needs no Python. See [contributing](../CONTRIBUTING.md) for
tests and release packaging; do not import a stale checkout extension to test
a newly built wheel.
