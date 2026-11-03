# Releasing to PyPI

A PyPI release is five files built from **one commit**:

| Files | Built by |
|---|---|
| Linux x86_64, Linux aarch64 and Windows x86_64 wheels | CI, on every push to `main`, `beta-v4` and `release/*` |
| macOS arm64 wheel and sdist | a maintainer's Mac |

`tools/release/prepare_pypi_release.py` collects, tests and checks all five
files. It never tags or uploads: publication is manual.

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

1. Set the version in `python/Cargo.toml`, then run `cargo update --workspace`
   to update `Cargo.lock`. Cargo `4.0.0-beta.2` becomes Python `4.0.0b2`.
2. Point the version-pinned GitHub links to the new tag, `v<python-version>`:
   READMEs, `[project.urls]` in `pyproject.toml` and `TAG` in
   `tools/tests/test_release_links.py`. `git grep -nF <old-version>` finds them.
   The script reads the tag from the `Documentation` URL.
3. Update `CHANGELOG.md`, run the [local checks](../CONTRIBUTING.md#build-and-test-locally),
   commit and push. The push starts CI.

## 2. Prepare the files

```bash
.venv/bin/python tools/release/prepare_pypi_release.py
```

The script stops at the first problem. It:

1. requires a clean, pushed HEAD and a tag that does not exist on another commit;
2. builds the macOS wheel and the sdist, runs `python/tests` on the installed
   wheel with Python 3.12 and 3.14, and builds the extension again from the sdist;
3. waits for the CI run of HEAD and accepts each CI wheel only if its report has
   the same commit, version and checksum, and passing tests;
4. runs `twine check --strict` and `tools/release/check_archives.py`;
5. writes the files to `target/pypi-release/<commit>/upload/` with
   `SHA256SUMS`, then prints the commands for step 3.

If CI or a test fails, fix it, commit, push and run the script again.
Nothing has been published.

## 3. Publish (manual)

Copy the commands printed by the script; the values below are examples:

```bash
cd target/pypi-release/0123456789ab && shasum -a 256 -c SHA256SUMS
git tag -a v4.0.0b2 0123456789ab -m 'Python 4.0.0b2' && git push origin refs/tags/v4.0.0b2
read -rsp 'PyPI token: ' UV_PUBLISH_TOKEN && export UV_PUBLISH_TOKEN
uv publish --trusted-publishing never --check-url https://pypi.org/simple/ upload/*
unset UV_PUBLISH_TOKEN
```

`read -rsp` keeps the token out of the screen and the shell history.
`--check-url` skips files already on PyPI, so the same command resumes an
interrupted upload. Check the published release from a clean environment:

```bash
uv run --isolated --no-project --refresh --python 3.12 \
  --with libephemeris==4.0.0b2 python -c 'import libephemeris; print(libephemeris.__version__)'
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
