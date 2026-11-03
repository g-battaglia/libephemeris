# Contributing

Use small, well-commented changes with tests and identified scientific sources.
Read the [independent implementation policy](docs/sources.md#independent-implementation-policy)
before changing a numerical model. Contributions follow the [CLA](CLA.md);
project code is [AGPL-3.0](COPYING).

## Workspace

- `core/`: Rust engine, tests, examples and optional immutable data crate.
- `python/`: PyO3 argument/result conversion, imports, stubs and Python tests.
- `tools/generators/`: reproducible numerical tables and retained snapshots.
- `tools/leb_generation/`: optional offline ephemeris tools.
- `tools/release/`: packaging, license and documentation checks.

A call enters the Python binding or Rust `api::Session`, uses retained resources
and native models, and returns owned results/errors. Keep computation and
configuration in Rust; Python is not a second engine. Public contracts document
clocks, units, domains and errors. Use rustfmt and the
[Rust API Guidelines](https://rust-lang.github.io/api-guidelines/).

## Build and test locally

Follow [installation](docs/installation.md#build-from-source). Use an installed
wheel, not a stale extension in the checkout. Install the test/docs dependency
groups with `uv` as needed; select the Python interpreter explicitly.

```sh
export PYO3_PYTHON="$PWD/.venv/bin/python"
cargo fmt --all --check
cargo test --locked --workspace --all-features --no-fail-fast
cargo test --locked --workspace --all-features --release --no-fail-fast
cargo clippy --locked --workspace --all-targets --all-features -- -D warnings
RUSTDOCFLAGS='-D warnings' cargo doc --locked --workspace --all-features --no-deps
.venv/bin/python -m pytest -q -ra python/tests tools/tests python/examples/tests
.venv/bin/python tools/release/check_minimal.py
.venv/bin/python tools/release/package_licenses.py --check
```

Check Rust 1.86 MSRV, default/no-default features and `bundled-data,nbody`
consumer builds when changing the Rust surface. Repeat wheel checks in clean
Python 3.12 and 3.14 environments. Local checks remain the reference. The only
remote pipeline is [wheels.yml](https://github.com/g-battaglia/libephemeris/actions/workflows/wheels.yml): it builds and
tests Linux/Windows wheels, without publishing. Other remote pipelines remain
forbidden by `tools/release/check_minimal.py`.

Tests can discover existing assets in `~/.libephemeris`;
`LIBEPHEMERIS_TEST_DATA_DIR` and `LIBEPHEMERIS_TEST_LEB2` select explicit test
inputs. The historical topocentric-photometry regression retains its data-v3
reference words and requires the exact data-v3 base image via
`LIBEPHEMERIS_TEST_DATA_V3_LEB2` (SHA-256
`5d708bdbe3e799e0802ba575984e57a3c5e44720dbfa1b4a01cf826640e0cb82`).
Provision that published asset separately; do not substitute the data-v4 bundle.
The separate data-v4 regression checks observer and TT/UT contracts without
replacing historical expected values. Report optional-data skips with `-ra`.
Do not repin hashes, replace
production datasets, weaken tolerances, enable fast-math or fit outputs to make
a test pass. Preserve generator identities and upstream notices. Numerical
qualification remains separate from product build/tests.

## Documentation

The manuals are maintained in parallel under `beginners-manual/it` and `en`.
Preserve their teaching depth; update examples and translations together.
Run the beginner path with the bundled data and no network. Advanced examples
must identify their data prerequisites; never call missing-data skips a pass.

```sh
.venv/bin/python tools/release/check_docs.py
.venv/bin/python -m pytest -q -ra tools/tests/test_check_docs.py \
  tools/tests/test_documentation_build.py tools/tests/test_manual_examples.py
.venv/bin/python -m sphinx -E -W --keep-going -b html -c docs . docs/_build/html
```

`tools/release/check_manual.py --data-dir /path/to/data --report /tmp/manual.json`
runs the books with explicit local assets; `--refresh-output` refreshes only
successfully executed adjacent transcripts. For pytest, set
`LIBEPHEMERIS_MANUAL_DATA_DIR=/path/to/data` to include those cases. Online
acquisition examples are syntax/signature checked, not downloaded by these
checks.

Sphinx reads the installed wheel, and includes only docs and manuals from the
repository root. Build output is not source. For table regeneration, use the
relevant script under `tools/generators/` and its focused tests; generation is
not a prerequisite for installing the shipped assets.

## Release checklist

Use a clean committed revision, lockfile and explicit interpreter. Versions:
core/data crates `1.0.0-beta.1`, Python/offline tools `4.0.0b1`; the binding's
internal Cargo version is `4.0.0-beta.1`. Runtime and installed metadata must
agree. Package versions are independent of scientific model identifiers.

```sh
cargo package --locked --all-features -p libephemeris-core-data -p libephemeris-core
.venv/bin/maturin build --release --locked --interpreter .venv/bin/python --out dist
.venv/bin/maturin sdist --out dist
uv build --out-dir dist/offline-tools tools/leb_generation
.venv/bin/python -m twine check --strict dist/*.whl dist/*.tar.gz dist/offline-tools/*
.venv/bin/python tools/release/check_archives.py dist/*.whl dist/*.tar.gz target/package/*.crate
```

Each compressed crate must be below 10,485,760 bytes. Do not remove data or add
runtime downloads to fit this limit. Extract archives outside the checkout;
build the extracted core and rebuild a wheel from the sdist. Verify licenses,
notices, stubs, entry points, bundled data and beginner manuals. Caches, local
repositories and scientific comparison environments must not ship.

Publication needs explicit maintainer approval: data crate first, then its
exact dependent core, then wheels/sdist and optional tools. Verify registry
installs separately; dry runs do not establish registry acceptance. Use distinct
Rust/Python tags and new prerelease versions for corrections, never overwrite
published artifacts. Do not expose credentials in logs.

### Linux and Windows wheels

The [wheel workflow](https://github.com/g-battaglia/libephemeris/actions/workflows/wheels.yml) runs on pushes to `beta-v4`
and `main`. Manual dispatch becomes available when the workflow is present on
GitHub's default branch; until v4 is promoted, trigger it by pushing `beta-v4`.
It runs three native jobs: Ubuntu 24.04 x86_64, Ubuntu 24.04 aarch64 and Windows
x86_64. Linux builds use manylinux2014 (glibc 2.17); Windows builds use MSVC.
If the Linux toolchain cannot build in that image, review an explicit change
to manylinux_2_28 and document the higher glibc floor. Never just rename a
wheel to claim compatibility.

Each job builds one ABI3 wheel with Maturin 1.9.4, installs it in separate
Python 3.12 and 3.14 virtual environments, and runs `python/tests`. Tests use
SHA-verified data-v4 base groups through the product CLI and the pinned JPL
DE440 image. The latter is 119,799,808 bytes, SHA-256
`a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5`.
No private harness or local scientific data directory is required. On Windows,
Git Bash's `HOME` is removed from child processes so discovery problems seen
by ordinary PowerShell/cmd users are not hidden.

The workflow has read-only repository permissions and no registry credentials,
OIDC write permission, publication, repository mutations, execution of Rust
tests or shared caches. Windows additionally type-checks all core targets,
including ignored native tests, to catch Unix-only compilation assumptions.
macOS builds, source-install checks and native Rust numerical qualification
stay local. The workflow guard is a review aid, not a sandbox for arbitrary
scripts.

For a new release:

1. Run the local checks, update the candidate's versions and version-pinned
   documentation links, commit it and push `beta-v4`. Do not move existing tags
   or overwrite versions already published to a registry.
2. Find and watch the exact run, using the product repository explicitly:

   ```sh
   gh run list -R g-battaglia/libephemeris --workflow wheels.yml --branch beta-v4
   gh run watch RUN_ID -R g-battaglia/libephemeris --exit-status
   gh run view RUN_ID -R g-battaglia/libephemeris --json headSha,conclusion
   ```

3. Require all three jobs to pass. Review the `-ra` output and JUnit reports:
   skipped optional-data tests are not verified coverage. Compare every
   artifact's `reports/build.json` source SHA, version and wheel checksum with
   the selected run. Failed jobs retain diagnostic reports and any wheel they
   managed to build; artifact presence alone is not approval to release.
4. Create the new release tag on that exact successful SHA and push the tag.
   Download its artifacts locally:

   ```sh
   gh run download RUN_ID -R g-battaglia/libephemeris -D dist/ci
   git rev-parse 'NEW_TAG^{commit}'
   ```

   Replace `RUN_ID` and `NEW_TAG` with the reviewed values. The tag's commit,
   the run's `headSha` and all artifact source SHAs must agree.
5. Build and test the macOS arm64 wheel and sdist locally as above. For an
   Intel macOS wheel, add `x86_64-apple-darwin` with `rustup target add`, build
   with `maturin build --release --locked --target x86_64-apple-darwin`, and
   install/test it using a real x86_64 Python under Rosetta. A successful
   cross-build alone does not count as a runtime test.
6. Collect only approved wheel/sdist files in a dedicated upload directory.
   Run `twine check --strict` and `tools/release/check_archives.py` on all of
   them, including the downloaded Linux/Windows wheels. Keep the downloaded
   provenance/JUnit reports, but do not upload those reports to PyPI.
7. After explicit maintainer approval, publish a new data crate first, then its
   matching core crate; never repeat uploads of already published versions.
   Use one manual `uv publish` for the approved wheel and sdist files. In Bash,
   read a PyPI token without echoing it or saving it in shell history:

   ```sh
   read -rsp 'PyPI token: ' UV_PUBLISH_TOKEN; echo
   export UV_PUBLISH_TOKEN
   uv publish --trusted-publishing never dist/upload/*.whl dist/upload/*.tar.gz
   unset UV_PUBLISH_TOKEN
   ```

8. Install the explicit prerelease from PyPI in clean environments and run the
   installed-package smoke checks. Verify Rust registry consumers separately.

The first CI run is diagnostic. Keep `v4.0.0b1` on its existing commit and the
published `1.0.0-beta.1` crates unchanged. Do not upload a wheel from a later
commit under that frozen release merely because its version string matches.
Windows home-directory and test-portability fixes, if needed, belong in a
subsequent change and new `4.0.0b2` / `1.0.0-beta.2` candidates. Do not mask
failures with CI-only environment workarounds, test deselection or weakened
numerical tolerances.

## Issues and pull requests

Provide a minimal public-API reproduction, versions, OS/architecture, data
identity and configuration. Explain significant numerical differences with
independent sources or physical invariants. Keep API changes in
[MIGRATION.md](MIGRATION.md), not in a running refactoring diary. Preserve
unrelated work, consumer repositories and production data.
