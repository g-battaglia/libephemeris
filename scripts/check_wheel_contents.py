#!/usr/bin/env python3
"""Build the distribution artifacts and audit their contents.

Licensing/packaging gate: asserts neither the sdist nor the wheel contains
test code, the dev CLI, or Swiss Ephemeris reference data, and that the
wheel carries the expected payload (package, typing marker, bundled data,
vendored modules).

The wheel is built FROM THE SDIST (python -m build's default two-step), so
the audit is hermetic against stale in-tree ``build/`` directories — a
stale ``build/lib`` is exactly how ``dev_cli/`` once leaked into a wheel.

Usage:
    python scripts/check_wheel_contents.py            # build sdist+wheel
    python scripts/check_wheel_contents.py path.whl   # audit existing wheel
"""

from __future__ import annotations

import re
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

FORBIDDEN = (
    re.compile(r"(^|/)tests?/"),
    re.compile(r"(^|/)compare_scripts/"),
    re.compile(r"libephemeris/dev_cli/"),
    re.compile(r"\.se1$", re.IGNORECASE),
    re.compile(r"sefstars", re.IGNORECASE),
    re.compile(r"seorbel", re.IGNORECASE),
    re.compile(r"data/reference/"),
)
WHEEL_REQUIRED = (
    "libephemeris/__init__.py",
    "libephemeris/py.typed",
    "libephemeris/data/fictitious_orbits.csv",
    "libephemeris/data/leb2/base_core.leb2",
    "libephemeris/vendor/spktype21.py",
    "libephemeris/moon_theories/galilean.py",
    "libephemeris/moon_theories/tass17.py",
    "libephemeris/moon_theories/tass17_data.py",
)


def build_artifacts(outdir: Path) -> tuple[Path, Path]:
    """Build sdist + wheel (wheel from the sdist) into outdir."""
    result = subprocess.run(
        [sys.executable, "-m", "build", "--outdir", str(outdir)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        sys.stderr.write(result.stdout[-2000:] + "\n" + result.stderr[-2000:] + "\n")
        raise SystemExit("build failed")
    return next(outdir.glob("*.tar.gz")), next(outdir.glob("*.whl"))


def forbidden_in(names: list[str]) -> list[str]:
    return sorted(
        {name for name in names for pattern in FORBIDDEN if pattern.search(name)}
    )


def audit(sdist: Path | None, wheel: Path) -> int:
    problems = 0

    if sdist is not None:
        with tarfile.open(sdist) as tf:
            sdist_names = tf.getnames()
        bad = forbidden_in(sdist_names)
        for name in bad:
            print(f"FORBIDDEN in sdist: {name}")
        print(
            f"sdist audit: {sdist.name} — {len(sdist_names)} entries, "
            f"{len(bad)} forbidden"
        )
        problems += len(bad)

    names = zipfile.ZipFile(wheel).namelist()
    bad = forbidden_in(names)
    missing = [required for required in WHEEL_REQUIRED if required not in names]
    for name in bad:
        print(f"FORBIDDEN in wheel: {name}")
    for name in missing:
        print(f"MISSING from wheel: {name}")
    print(
        f"wheel audit: {wheel.name} — {len(names)} entries, "
        f"{len(bad)} forbidden, {len(missing)} missing"
    )
    problems += len(bad) + len(missing)
    return 1 if problems else 0


def main() -> int:
    if len(sys.argv) > 1:
        return audit(None, Path(sys.argv[1]))
    with tempfile.TemporaryDirectory() as tmp:
        sdist, wheel = build_artifacts(Path(tmp))
        return audit(sdist, wheel)


if __name__ == "__main__":
    sys.exit(main())
