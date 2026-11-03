#!/usr/bin/env python3
"""Prepare the files of one PyPI release of HEAD; never tags or uploads.

Waits for the CI run of HEAD and accepts its Linux/Windows wheels only if their
reports match HEAD, builds and tests the macOS wheel and the sdist locally,
checks all files and writes them with SHA256SUMS to target/pypi-release/<commit>.
Publication stays manual: see docs/releasing.md.

Run on macOS from a clean checkout of a commit pushed to main, beta-v4 or release/*.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import tomllib
import urllib.error
import urllib.request

ROOT = Path(__file__).resolve().parents[2]
REPO = "g-battaglia/libephemeris"
WORKFLOW = "wheels.yml"
CI_PLATFORMS = ("linux-x86_64", "linux-aarch64", "windows-x86_64")
PYTHONS = ("3.12", "3.14")
VENV = ROOT / ".venv/bin"
# Import the installed package, never the checkout, and check its version.
SMOKE = (
    "import sys; from pathlib import Path; import libephemeris as e; "
    "assert e.__version__ == sys.argv[1], e.__version__; "
    "assert Path(e.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()); "
    "print('libephemeris', e.__version__, e.__file__)"
)


class ReleaseError(Exception):
    """A precondition failed; the message says what to do next."""


def run(command: list, **kwargs) -> None:
    command = [str(part) for part in command]
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True, **kwargs)


def output(command: list, check: bool = True) -> str:
    command = [str(part) for part in command]
    return subprocess.run(command, check=check, capture_output=True, text=True).stdout.strip()


def sha256(path: Path) -> str:
    with path.open("rb") as file:
        return hashlib.file_digest(file, "sha256").hexdigest()


def one(paths: list[Path], what: str) -> Path:
    if len(paths) != 1:
        raise ReleaseError(f"expected one {what}, found {len(paths)}")
    return paths[0]


def on_pypi(version: str) -> bool:
    try:
        urllib.request.urlopen(f"https://pypi.org/pypi/libephemeris/{version}/json", timeout=30)
        return True
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return False
        raise


def checked_head() -> str:
    """Return HEAD when it is committed, on a branch and pushed unchanged."""
    if output(["git", "-C", ROOT, "status", "--porcelain"]):
        raise ReleaseError("the checkout has local changes: commit or stash them first")
    branch = output(["git", "-C", ROOT, "branch", "--show-current"])
    if not branch:
        raise ReleaseError("detached HEAD: check out main, beta-v4 or a release/* branch")
    sha = output(["git", "-C", ROOT, "rev-parse", "HEAD"])
    remote = output(["git", "-C", ROOT, "ls-remote", "origin", f"refs/heads/{branch}"]).split()
    if remote[:1] != [sha]:
        raise ReleaseError(f"origin/{branch} is not HEAD: run `git push origin {branch}` first")
    return sha


def pinned_tag(pyproject: str) -> str:
    """Return the tag pinned by the Documentation URL (.../blob/<tag>/...)."""
    url = tomllib.loads(pyproject)["project"]["urls"]["Documentation"]
    if "/blob/" not in url:
        raise ReleaseError("pyproject Documentation URL must pin the release tag: .../blob/<tag>/...")
    return url.split("/blob/", 1)[1].split("/", 1)[0]


def release_tag(sha: str) -> str:
    """Return the pinned tag if it is new or already on HEAD; tags never move."""
    tag = pinned_tag((ROOT / "pyproject.toml").read_text())
    refs = output(["git", "-C", ROOT, "ls-remote", "origin", f"refs/tags/{tag}", f"refs/tags/{tag}^{{}}"])
    targets = dict(reversed(line.split("\t")) for line in refs.splitlines())
    remote = targets.get(f"refs/tags/{tag}^{{}}", targets.get(f"refs/tags/{tag}"))
    local = output(["git", "-C", ROOT, "rev-parse", "-q", "--verify", f"refs/tags/{tag}^{{commit}}"], check=False)
    for target in (remote, local):
        if target and target != sha:
            raise ReleaseError(
                f"tag {tag}, pinned by the package links, already exists on {target[:12]}: "
                "pin a new tag in the links, commit and push"
            )
    return tag


def find_run(sha: str) -> str:
    """Return the newest wheel workflow run of the commit; a fresh push may need a moment."""
    for _ in range(12):
        runs = json.loads(output([
            "gh", "run", "list", "-R", REPO, "--workflow", WORKFLOW, "--commit", sha,
            "--json", "databaseId", "-L", "20",
        ]))
        if runs:
            return str(max(run["databaseId"] for run in runs))
        time.sleep(10)
    raise ReleaseError(
        f"no {WORKFLOW} run for {sha[:12]}: CI builds pushes to main, beta-v4 and release/*"
    )


def report_problems(report: dict, wheel: Path, sha: str, version: str) -> list[str]:
    """Compare one CI build.json with the wheel next to it and with HEAD."""
    problems = []
    if report["source_sha"] != sha:
        problems.append(f"built from {report['source_sha'][:12]}, not HEAD")
    if report["version"] != version:
        problems.append(f"version {report['version']}, not {version}")
    if report["wheel"] != wheel.name or report["sha256"] != sha256(wheel):
        problems.append("wheel name or checksum differs from the report")
    for python in PYTHONS:
        code = report["tests"].get(python, {}).get("returncode")
        if code != 0:
            problems.append(f"Python {python} tests exited with {code}")
    return problems


def ci_wheels(run_id: str, sha: str, version: str, folder: Path) -> list[Path]:
    """Download the CI wheels and accept them only if their reports match HEAD."""
    run(["gh", "run", "download", run_id, "-R", REPO, "-D", folder])
    wheels = []
    for platform in CI_PLATFORMS:
        artifact = folder / f"wheels-{platform}"
        wheel = one(sorted((artifact / "dist").glob("*.whl")), f"{platform} wheel")
        report = json.loads((artifact / "reports/build.json").read_text())
        problems = report_problems(report, wheel, sha, version)
        if problems:
            raise ReleaseError(f"{platform} CI wheel rejected: {'; '.join(problems)}")
        wheels.append(wheel)
    return wheels


def new_venv(path: Path, python: str) -> Path:
    run(["uv", "venv", "--quiet", "--clear", "--python", python, path])
    return path / "bin/python"


def prepare() -> None:
    if sys.platform != "darwin":
        raise ReleaseError("run on macOS: CI already provides the Linux and Windows wheels")
    sha = checked_head()
    tag = release_tag(sha)
    work = ROOT / "target/pypi-release" / sha[:12]
    # Files of a (partly) published version must never be rebuilt: resume with them.
    for wheel in (work / "upload").glob("*.whl"):
        if on_pypi(wheel.name.split("-")[1]):
            raise ReleaseError(f"{work / 'upload'} is already (partly) on PyPI: upload those files, never rebuild")
    shutil.rmtree(work, ignore_errors=True)
    run_id = find_run(sha)
    print(f"CI run: https://github.com/{REPO}/actions/runs/{run_id}", flush=True)

    # Build and test locally while CI is still running.
    local = work / "local"
    run([VENV / "maturin", "build", "--release", "--locked", "--interpreter", VENV / "python", "--out", local])
    run([VENV / "maturin", "sdist", "--out", local])
    wheel = one(sorted(local.glob("*.whl")), "local wheel")
    sdist = one(sorted(local.glob("*.tar.gz")), "sdist")
    version = wheel.name.split("-")[1]
    if on_pypi(version):
        raise ReleaseError(f"libephemeris {version} is already on PyPI: bump the version")
    # Same bundled LEB2 image as CI; other data comes from the usual test locations.
    env = dict(os.environ, LIBEPHEMERIS_TEST_LEB2=str(ROOT / "core/data/leb2/base_core.leb2"))
    for python in PYTHONS:
        interpreter = new_venv(work / f"venv-{python}", python)
        run(["uv", "pip", "install", "--quiet", "--python", interpreter, "pytest", wheel])
        run([interpreter, "-I", "-c", SMOKE, version], cwd=work)
        run([interpreter, "-m", "pytest", "-q", "-ra", ROOT / "python/tests"], cwd=work, env=env)
    # Installing the sdist rebuilds the extension in a directory outside the checkout.
    interpreter = new_venv(work / "venv-sdist", PYTHONS[0])
    run(["uv", "pip", "install", "--quiet", "--no-cache", "--python", interpreter, sdist])
    run([interpreter, "-I", "-c", SMOKE, version], cwd=work)

    run(["gh", "run", "watch", run_id, "-R", REPO, "--exit-status", "--interval", "30"])
    upload = work / "upload"
    upload.mkdir()
    for path in [*ci_wheels(run_id, sha, version, work / "ci"), wheel, sdist]:
        shutil.copy2(path, upload)
    files = sorted(upload.iterdir())
    run([VENV / "python", "-m", "twine", "check", "--strict", *files])
    run([VENV / "python", ROOT / "tools/release/check_archives.py", *files])
    sums = "".join(f"{sha256(path)}  upload/{path.name}\n" for path in files)
    (work / "SHA256SUMS").write_text(sums)

    print(f"\nReady: libephemeris {version} from {sha[:12]} (CI run {run_id}), {len(files)} files:")
    print(sums)
    print(f"Publish manually (docs/releasing.md):\n"
          f"  cd {work.relative_to(ROOT)} && shasum -a 256 -c SHA256SUMS\n"
          f"  git tag -a {tag} {sha} -m 'Python {version}' && git push origin refs/tags/{tag}\n"
          f"  read -rsp 'PyPI token: ' UV_PUBLISH_TOKEN && export UV_PUBLISH_TOKEN\n"
          f"  uv publish --trusted-publishing never upload/*\n"
          f"  unset UV_PUBLISH_TOKEN")


def main() -> int:
    if sys.argv[1:]:
        print(__doc__)
        return 0 if sys.argv[1:] in (["-h"], ["--help"]) else 2
    try:
        prepare()
    except ReleaseError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    except subprocess.CalledProcessError as error:
        print(f"error: command failed with exit code {error.returncode}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
