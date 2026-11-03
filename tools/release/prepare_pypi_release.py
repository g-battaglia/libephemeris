#!/usr/bin/env python3
"""Prepare one verified PyPI release; never create tags or upload packages.

Run from a clean, pushed macOS checkout. The resulting publish.sh is a separate,
explicit maintainer action: it verifies the frozen files and asks for a token.
"""

from __future__ import annotations

from email.parser import BytesParser
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import tomllib
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[2]
REPO = "g-battaglia/libephemeris"
WORKFLOW = "wheels.yml"
CI_PLATFORMS = ("linux-x86_64", "linux-aarch64", "windows-x86_64")
TARGETS = {
    "linux-x86_64": ("x86_64-unknown-linux-gnu", "manylinux_2_17_x86_64"),
    "linux-aarch64": ("aarch64-unknown-linux-gnu", "manylinux_2_17_aarch64"),
    "windows-x86_64": ("x86_64-pc-windows-msvc", "win_amd64"),
}
PYTHONS = ("3.12", "3.14")
VENV = ROOT / ".venv/bin"
DE440_SHA256 = "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5"
SMOKE = (
    "import sys; from pathlib import Path; from importlib import metadata; import libephemeris as e; "
    "assert e.__version__ == metadata.version('libephemeris') == sys.argv[1]; "
    "assert Path(e.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()); "
    "print('libephemeris', e.__version__, e.__file__)"
)


class ReleaseError(Exception):
    """A release precondition failed; never continue to publication."""


def run(command: list, **kwargs) -> None:
    command = [str(part) for part in command]
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True, **kwargs)


def output(command: list, check: bool = True) -> str:
    return subprocess.run([str(part) for part in command], check=check,
                          capture_output=True, text=True).stdout.strip()


def sha256(path: Path) -> str:
    with path.open("rb") as file:
        return hashlib.file_digest(file, "sha256").hexdigest()


def one(paths: list[Path], what: str) -> Path:
    if len(paths) != 1:
        raise ReleaseError(f"expected one {what}, found {len(paths)}")
    return paths[0]


def on_pypi(version: str) -> bool:
    try:
        with urllib.request.urlopen(f"https://pypi.org/pypi/libephemeris/{version}/json", timeout=30):
            return True
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return False
        raise


def checked_head() -> str:
    """Require a clean checkout and the same commit on the remote branch."""
    if output(["git", "-C", ROOT, "status", "--porcelain"]):
        raise ReleaseError("commit or stash local changes before preparing a release")
    branch = output(["git", "-C", ROOT, "branch", "--show-current"])
    if branch not in ("main", "beta-v4") and not branch.startswith("release/"):
        raise ReleaseError("check out a branch built by CI: main, beta-v4 or release/*")
    sha = output(["git", "-C", ROOT, "rev-parse", "HEAD"])
    remote = output(["git", "-C", ROOT, "ls-remote", "origin", f"refs/heads/{branch}"]).split()
    if remote[:1] != [sha]:
        raise ReleaseError(f"origin/{branch} is not HEAD: push the reviewed commit first")
    return sha


def pinned_tag(pyproject: str) -> str:
    url = tomllib.loads(pyproject)["project"]["urls"]["Documentation"]
    if "/blob/" not in url:
        raise ReleaseError("Documentation must pin a source tag: .../blob/<tag>/...")
    return url.split("/blob/", 1)[1].split("/", 1)[0]


def release_tag(sha: str) -> str:
    """Accept a new tag or a tag already on HEAD; never move any existing tag."""
    tag = pinned_tag((ROOT / "pyproject.toml").read_text())
    run(["git", "check-ref-format", f"refs/tags/{tag}"])
    refs = output(["git", "-C", ROOT, "ls-remote", "origin", f"refs/tags/{tag}", f"refs/tags/{tag}^{{}}"])
    targets = dict(reversed(line.split("\t")) for line in refs.splitlines())
    remote = targets.get(f"refs/tags/{tag}^{{}}", targets.get(f"refs/tags/{tag}"))
    local = output(["git", "-C", ROOT, "rev-parse", "-q", "--verify", f"refs/tags/{tag}^{{commit}}"], check=False)
    if any(target and target != sha for target in (remote, local)):
        raise ReleaseError(f"{tag} already identifies another commit: choose and repin a new tag")
    return tag


def find_run(sha: str) -> str:
    for _ in range(12):
        runs = json.loads(output([
            "gh", "run", "list", "-R", REPO, "--workflow", WORKFLOW, "--commit", sha,
            "--json", "databaseId", "-L", "20",
        ]))
        if runs:
            return str(max(item["databaseId"] for item in runs))
        time.sleep(10)
    raise ReleaseError(f"no {WORKFLOW} run for {sha[:12]}: push the candidate to a CI branch")


def check_ci_run(record: dict, sha: str) -> None:
    if record["headSha"] != sha or record["status"] != "completed" or record["conclusion"] != "success":
        raise ReleaseError("CI is not a successful completed run of this exact commit")
    if {job["name"] for job in record["jobs"]} != set(CI_PLATFORMS):
        raise ReleaseError("CI does not contain the three required native jobs")
    if any(job["conclusion"] != "success" for job in record["jobs"]):
        raise ReleaseError("a required CI job did not succeed")


def wait_for_ci(run_id: str, sha: str) -> dict:
    # An API 500 can interrupt `watch` without stopping CI. Retry monitoring,
    # never rebuild/rerun the product merely to repair a monitoring failure.
    for _ in range(3):
        subprocess.run(["gh", "run", "watch", run_id, "-R", REPO, "--exit-status", "--interval", "30"])
        try:
            record = json.loads(output(["gh", "run", "view", run_id, "-R", REPO,
                                        "--json", "headSha,status,conclusion,jobs,url"]))
        except subprocess.CalledProcessError:
            time.sleep(10)
            continue
        if record["status"] == "completed":
            check_ci_run(record, sha)
            return record
        time.sleep(10)
    raise ReleaseError("CI monitoring is unfinished; inspect the run without treating it as a pass")


def report_problems(report: dict, wheel: Path, sha: str, version: str) -> list[str]:
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


def junit_summary(path: Path) -> dict:
    """Require executed tests with no failures; retain every optional-data skip."""
    cases = list(ET.parse(path).getroot().iter("testcase"))
    failed = [case for case in cases if case.find("failure") is not None or case.find("error") is not None]
    skipped = [{"test": case.get("classname", "") + "::" + case.get("name", ""),
                "reason": case.find("skipped").get("message", "")}
               for case in cases if case.find("skipped") is not None]
    passed = len(cases) - len(failed) - len(skipped)
    if failed or passed == 0:
        raise ReleaseError(f"{path}: failed, empty, or entirely skipped test suite")
    return {"passed": passed, "failed": 0, "skipped": skipped}


def ci_wheels(run_id: str, sha: str, version: str, folder: Path) -> list[Path]:
    run(["gh", "run", "download", run_id, "-R", REPO, "-D", folder])
    wheels = []
    for platform in CI_PLATFORMS:
        artifact = folder / f"wheels-{platform}"
        wheel = one(sorted((artifact / "dist").glob("*.whl")), f"{platform} wheel")
        reports = artifact / "reports"
        record = json.loads((reports / "build.json").read_text())
        problems = report_problems(record, wheel, sha, version)
        target, tag = TARGETS[platform]
        if record["target"] != target or "cp312-abi3-" not in wheel.name or tag not in wheel.name:
            problems.append("incorrect target or wheel ABI/platform tag")
        if (reports / "source.txt").read_text() != f"Commit: {sha}\nTarget: {target}\n":
            problems.append("source report differs from the candidate")
        with zipfile.ZipFile(wheel) as archive:
            name = one([Path(name) for name in archive.namelist() if name.endswith(".dist-info/METADATA")], "wheel metadata")
            metadata = BytesParser().parsebytes(archive.read(name.as_posix()))
            if metadata["Name"] != "libephemeris" or metadata["Version"] != version:
                problems.append("incorrect installed package metadata")
        if problems:
            raise ReleaseError(f"{platform} CI wheel rejected: {'; '.join(problems)}")
        for python in PYTHONS:
            junit_summary(reports / f"pytest-{python}.xml")
        wheels.append(wheel)
    return wheels


def new_venv(path: Path, python: str) -> Path:
    run(["uv", "venv", "--quiet", "--python", python, path])
    return path / "bin/python"


def provision_test_data(interpreter: Path, data: Path) -> None:
    """Use the same authenticated public inputs as CI, never a private harness."""
    run([interpreter.with_name("libephemeris"), "download", "leb2-base", "--data-dir", data, "--no-progress"])
    kernel = data / "de440.bsp"
    if not kernel.exists():
        partial = kernel.with_suffix(".part")
        with urllib.request.urlopen("https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440.bsp", timeout=120) as source:
            with partial.open("wb") as destination:
                shutil.copyfileobj(source, destination, length=1024 * 1024)
        if partial.stat().st_size != 119799808 or sha256(partial) != DE440_SHA256:
            raise ReleaseError("downloaded DE440 does not match the pinned CI input")
        partial.replace(kernel)
    if kernel.stat().st_size != 119799808 or sha256(kernel) != DE440_SHA256:
        raise ReleaseError("cached DE440 does not match the pinned CI input")


def installed_tests(interpreter: Path, work: Path, environment: dict, label: str, version: str) -> dict:
    run([interpreter, "-I", "-c", SMOKE, version], cwd=work, env=environment)
    report = work / "reports" / f"pytest-{label}.xml"
    run([interpreter, "-m", "pytest", "-q", "-ra", ROOT / "python/tests",
         f"--junitxml={report}"], cwd=work, env=environment)
    return junit_summary(report)


def write_publisher(work: Path, sha: str, tag: str) -> None:
    """Write, but never execute, the maintainer's final local publish command."""
    python = shlex.quote(str(VENV / "python"))
    root = shlex.quote(str(ROOT))
    script = f'''#!/usr/bin/env bash
set -euo pipefail
set +x
cd -- "$(dirname -- "$0")"
{python} - {root} {shlex.quote(sha)} {shlex.quote(tag)} <<'PY'
import hashlib, json, subprocess, sys
from pathlib import Path
root, sha, tag = sys.argv[1:]
record = json.loads(Path('release.json').read_text())
assert record['sha'] == sha and record['tag'] == tag
assert subprocess.check_output(['git', '-C', root, 'rev-parse', 'HEAD'], text=True).strip() == sha
remote = subprocess.check_output(['git', '-C', root, 'ls-remote', 'origin', 'refs/tags/' + tag + '^{{}}'], text=True).split()
assert remote[:1] == [sha], 'the reviewed release tag must already identify this commit'
actual = {{path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in Path('upload').iterdir()}}
assert actual == record['files'], 'the upload directory differs from the verified files'
print('Publishing libephemeris', record['version'], 'from', sha)
PY
if [ -z "${{UV_PUBLISH_TOKEN:-}}" ]; then
  read -rsp 'PyPI token (local, hidden): ' UV_PUBLISH_TOKEN
  printf '\\n'
fi
export UV_PUBLISH_TOKEN
trap 'unset UV_PUBLISH_TOKEN' EXIT
uv publish --trusted-publishing never --publish-url https://upload.pypi.org/legacy/ \\
  --check-url https://pypi.org/simple/ upload/*.whl upload/*.tar.gz
'''
    publisher = work / "publish.sh"
    publisher.write_text(script)
    publisher.chmod(0o700)


def prepare() -> None:
    if sys.platform != "darwin":
        raise ReleaseError("prepare on macOS; CI provides the native Linux/Windows wheels")
    sha = checked_head()
    tag = release_tag(sha)
    work = ROOT / "target/pypi-release" / sha[:12]
    if work.exists():
        raise ReleaseError(f"{work} already exists: preserve/review it; never rebuild release files in place")
    run_id = find_run(sha)
    print(f"CI run: https://github.com/{REPO}/actions/runs/{run_id}", flush=True)
    (work / "reports").mkdir(parents=True)
    local = work / "local"
    run([VENV / "maturin", "build", "--release", "--locked", "--interpreter", VENV / "python", "--out", local])
    run([VENV / "maturin", "sdist", "--out", local])
    wheel = one(sorted(local.glob("*.whl")), "local wheel")
    sdist = one(sorted(local.glob("*.tar.gz")), "sdist")
    version = wheel.name.split("-")[1]
    if on_pypi(version):
        raise ReleaseError(f"libephemeris {version} already exists on PyPI: choose a new version")
    data = ROOT / "target/pypi-release-inputs"
    environment = dict(os.environ, LIBEPHEMERIS_TEST_DATA_DIR=str(data),
                       LIBEPHEMERIS_TEST_LEB2=str(ROOT / "core/data/leb2/base_core.leb2"))
    tests = {}
    for python in PYTHONS:
        interpreter = new_venv(work / f"venv-{python}", python)
        run(["uv", "pip", "install", "--quiet", "--python", interpreter, "pytest", wheel])
        if python == PYTHONS[0]:
            provision_test_data(interpreter, data)
        tests[f"macos-{python}"] = installed_tests(interpreter, work, environment, f"macos-{python}", version)

    # Extract outside the checkout and build with locked, cached Rust inputs.
    rebuilt = work / "sdist-rebuilt"
    with tempfile.TemporaryDirectory(prefix="libephemeris-sdist-") as temporary:
        source = Path(temporary)
        with tarfile.open(sdist) as archive:
            archive.extractall(source, filter="data")
        extracted = one(list(source.iterdir()), "sdist root")
        build_env = dict(os.environ, CARGO_NET_OFFLINE="true", PYO3_PYTHON=str(VENV / "python"),
                         CARGO_TARGET_DIR=str(ROOT / "target/pypi-sdist-build"))
        run([VENV / "maturin", "build", "--release", "--locked", "--offline", "--out", rebuilt],
            cwd=extracted, env=build_env)
    interpreter = new_venv(work / "venv-sdist", PYTHONS[0])
    run(["uv", "pip", "install", "--quiet", "--python", interpreter, "pytest",
         one(list(rebuilt.glob("*.whl")), "rebuilt sdist wheel")])
    tests["sdist-3.12"] = installed_tests(interpreter, work, environment, "sdist-3.12", version)

    ci = wait_for_ci(run_id, sha)
    (work / "reports/run.json").write_text(json.dumps(ci, indent=2) + "\n")
    wheels = ci_wheels(run_id, sha, version, work / "ci")
    upload = work / "upload"
    upload.mkdir()
    for path in [*wheels, wheel, sdist]:
        shutil.copy2(path, upload)
    files = sorted(upload.iterdir())
    if len(files) != 5:
        raise ReleaseError("expected four native wheels and one source distribution")
    run([VENV / "python", "-m", "twine", "check", "--strict", *files])
    run([VENV / "python", ROOT / "tools/release/check_archives.py", *files])
    record = {"version": version, "sha": sha, "tag": tag, "ci_run": run_id,
              "files": {path.name: sha256(path) for path in files}, "local_tests": tests,
              "de440_sha256": DE440_SHA256}
    (work / "release.json").write_text(json.dumps(record, indent=2) + "\n")
    sums = "".join(f"{sha256(path)}  upload/{path.name}\n" for path in files)
    (work / "SHA256SUMS").write_text(sums)
    write_publisher(work, sha, tag)
    print(f"\nVerified {version}: {len(files)} files from {sha}; nothing has been uploaded.")
    print(f"Freeze the reviewed new tag {tag} on this SHA before publication.")
    print(f"Then run only: bash {work.relative_to(ROOT)}/publish.sh")


def main() -> int:
    if sys.argv[1:]:
        print(__doc__)
        return 0 if sys.argv[1:] in (["-h"], ["--help"]) else 2
    try:
        prepare()
    except (ReleaseError, subprocess.CalledProcessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
