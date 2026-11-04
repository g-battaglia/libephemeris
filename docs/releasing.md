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

1. Choose unused registry versions. This candidate advances Python to `4.0.0b2`
   (Cargo `4.0.0-beta.2`) and retains Rust crates `1.0.0-beta.2`. Review the
   exact dependencies and `Cargo.lock`; scientific data-v4 bytes stay unchanged.
2. Pin Python README/project source links to `python-v4.0.0b2`, not the frozen
   `v4.0.0b1` or `python-v4.0.0b1`. Descriptions for unchanged Rust/tools packages
   retain their existing release-tag links. The script and Sphinx read the source tag from the `Documentation`
   URL in `pyproject.toml`; release-link tests must agree.
3. Update `CHANGELOG.md`, run the [local checks](../CONTRIBUTING.md#build-and-test-locally),
   commit and push. The push starts CI.

## 2. Prepare the files

```bash
.venv/bin/python tools/release/prepare_pypi_release.py
```

The script stops at the first problem. It:

1. requires a clean, pushed HEAD and a tag that does not exist on another commit;
2. builds macOS/sdist, provisions the same pinned public inputs as CI, runs
   `python/tests` on Python 3.12/3.14 and on a locked, offline sdist rebuild;
3. waits for the CI run of HEAD and accepts each CI wheel only if its report has
   the same commit, version and checksum, and passing tests;
4. runs `twine check --strict` and `tools/release/check_archives.py`;
5. writes five files to `target/pypi-release/<commit>/upload/`, plus checksums,
   provenance, JUnit reports and a separate `publish.sh`.

If CI or a test fails, fix it, commit, push and run the script again.
Nothing has been published.

## 3. Publish (manual)

After every check passes, the maintainer freezes the **new** source tag on
that exact SHA. Never move `v4.0.0b1`. The final operator then runs only the
command printed by preparation, for example:

```bash
bash target/pypi-release/0123456789ab/publish.sh
```

It verifies HEAD, the remote tag and all five file checksums, prompts for the
PyPI token **locally without echoing it**, and invokes `uv publish`. No token is
stored. The same files can resume a partial upload with `--check-url`; never
rebuild them after any file is published. Verify the installed release:

```bash
uv run --isolated --no-project --refresh --python 3.12 \
  --with libephemeris==4.0.0b2 python -I -c 'import libephemeris; print(libephemeris.__version__)'
```

## From `main` or from a dedicated branch

The steps are the same; only the branch changes.

| | From `main` | From a dedicated branch |
|---|---|---|
| Use for | regular releases | betas, fixes for an older line, releasing before merging |
| Start | merge the release PR, then `git switch main && git pull --ff-only` | `git switch -c release/4.0.0b2` from the commit to release, then `git push -u origin release/4.0.0b2` after step 1 |
| Afterwards | nothing | merge the branch into `main` with a PR if its changes belong there |

The tag keeps the released commit reachable even if the branch is deleted.

## Rules

- A published version is final. PyPI rejects a file name uploaded before, so a
  correction needs a new version. Never move a tag.
- After a partial upload, publish the same `upload/` files: the script refuses
  to rebuild a version that is already on PyPI.
- Intel Macs have no wheel; pip builds them from the sdist, which needs Rust.
- Rust crates have their own versions and are published separately; see the
  [release checklist](../CONTRIBUTING.md#release-checklist).
