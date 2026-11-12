# Releasing to PyPI

A PyPI release is five files built from **one commit**:

| Files | Built by |
|---|---|
| Linux x86_64, Linux aarch64 and Windows x86_64 wheels | CI, on every push to `main`, `beta-v4` and `release/*` |
| macOS arm64 wheel and sdist | a maintainer's Mac |

`tools/release/prepare_pypi_release.py` collects, tests and checks all five
files and writes a separate `publish.sh`. Preparation never tags or uploads;
the final publication command is an explicit maintainer action.

## Once per machine

```bash
gh auth login
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python --group release
```

You also need the Rust toolchain ([installation](installation.md#build-from-source))
and a PyPI API token scoped to the `libephemeris` project. Local tests read
data from `~/.libephemeris` or `LIBEPHEMERIS_TEST_DATA_DIR`, as usual.

## 1. Prepare the release commit

1. Choose unused registry versions. This preparation promotes Python `4.0.0`
   (binding Cargo `4.0.0`) and Rust core/data `1.0.0` from the qualified Python
   b4 source. Update the exact internal dependencies and workspace records in
   `Cargo.lock`; retain every external dependency and all scientific data bytes.
2. Pin Python and stable Rust README/project source links to `python-v4.0.0`.
   The unchanged optional tools keep their own existing version and source tag.
   The script and Sphinx read the source tag from the `Documentation` URL in
   `pyproject.toml`; release-link tests must agree.
3. Update the stable classifiers, `CHANGELOG.md`, release notes and migration
   instructions. Run the [local checks](../CONTRIBUTING.md#build-and-test-locally),
   commit and push the reviewed preparation. The push starts the existing wheel CI.
   Never substitute a local dirty build for the committed release candidate.

## 2. Prepare the files

```bash
.venv/bin/python tools/release/prepare_pypi_release.py
```

The script stops at the first problem. It:

1. requires a clean, pushed HEAD and a source tag that does not identify another commit;
2. builds macOS/sdist, provisions the same pinned public inputs as CI, runs
   `python/tests` on Python 3.12/3.14 and on a locked, offline sdist rebuild;
3. waits for the CI run of HEAD and accepts each CI wheel only if its report has
   the same commit, version and checksum, and passing tests;
4. runs `twine check --strict` and `tools/release/check_archives.py`;
5. writes five files to `target/pypi-release/<commit>/upload/`, plus checksums,
   provenance, JUnit reports and a separate `publish.sh`.

If CI or a test fails, fix it, commit, push and prepare the new candidate.
Preserve earlier artifacts and reports. Nothing has been published.

## 3. Publish (manual)

After every check passes, the maintainer freezes the **new** source tag on
that exact SHA. Never move an earlier tag. Publication then uses only the
verified upload directory, for example:

```bash
bash target/pypi-release/0123456789ab/publish.sh
```

The generated publisher verifies HEAD, the remote annotated tag and all five
file checksums, prompts for the PyPI token **locally without echoing it**, and
invokes `uv publish`. No token is stored. If the maintainer instead invokes
`uv publish` directly, first complete those same checks and name only the five
frozen files; never select a general `dist/` directory containing old versions.

The same files can resume a partial upload with `--check-url`; never rebuild
them after any file is published. Verify the installed release separately:

```bash
uv run --isolated --no-project --refresh --python 3.12 \
  --with libephemeris==4.0.0 python -I -c 'import libephemeris; print(libephemeris.__version__)'
```

Do not claim that a source preparation, tag or local wheel establishes registry
publication or availability of every platform artifact.

## From `main` or from a dedicated branch

The steps are the same; only the branch changes.

| | From `main` | From a dedicated branch |
|---|---|---|
| Use for | regular releases | preparation before merging or an older release line |
| Start | merge the qualified release PR, then `git switch main && git pull --ff-only` | create and push `release/4.0.0` from the commit to release |
| Afterwards | nothing | merge the qualified branch into `main` if its changes belong there |

Do not rebuild already published files after a merge. The immutable source tag
keeps the exact released commit reachable.

## Rules

- A published version is final. PyPI rejects a previously uploaded filename;
  a correction needs a new version. Never move source tags.
- After a partial upload, publish the same `upload/` files: the preparation
  script refuses to rebuild a version that is already on PyPI.
- Intel Macs have no wheel; pip builds from the sdist, which needs Rust.
- Rust crates have independent versions and registry publication. The Python
  wheel embeds its compiled engine and the sdist ships the Rust sources: neither
  requires stable crates to have been uploaded. See the separate
  [release checklist](../CONTRIBUTING.md#release-checklist) for crates.io.
