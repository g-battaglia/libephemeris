"""Check the release guide and the offline parts of the preparation script.

Nothing here builds, tags, prompts for a token or uploads.
"""

import importlib.util
from pathlib import Path
import re
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[2]
GUIDE = (ROOT / "docs/releasing.md").read_text()
SCRIPT = ROOT / "tools/release/prepare_pypi_release.py"
SPEC = importlib.util.spec_from_file_location("prepare_pypi_release", SCRIPT)
PREPARE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PREPARE)


def test_release_guide_is_linked_from_the_book_and_contributing():
    assert "   releasing\n" in (ROOT / "docs/index.rst").read_text()
    assert "(releasing.md)" in (ROOT / "docs/README.md").read_text()
    assert "(docs/releasing.md)" in (ROOT / "CONTRIBUTING.md").read_text()


@pytest.mark.parametrize("block", re.findall(r"```bash\n(.*?)```", GUIDE, flags=re.DOTALL))
def test_guide_shell_blocks_parse(block):
    # -n only parses the commands; none of them runs.
    result = subprocess.run(["bash", "-n"], input=block, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr


def test_ci_builds_the_branches_named_by_the_guide():
    workflow = (ROOT / ".github/workflows/wheels.yml").read_text()
    assert "branches: [main, beta-v4, 'release/**']" in workflow
    assert "`main`, `beta-v4` and `release/*`" in GUIDE


def test_preparation_never_tags_pushes_or_uploads():
    commands = re.findall(r"(?:run|output)\(\[(.*?)\]", SCRIPT.read_text(), flags=re.DOTALL)
    assert commands
    for command in commands:
        assert not re.search(r'"(?:tag|push|publish|upload)"', command), command


def test_release_tag_comes_from_the_pinned_documentation_url():
    pinned = '[project.urls]\nDocumentation = "https://github.com/o/r/blob/v9.0.0b1/docs/README.md"\n'
    assert PREPARE.pinned_tag(pinned) == "v9.0.0b1"
    with pytest.raises(PREPARE.ReleaseError):
        PREPARE.pinned_tag('[project.urls]\nDocumentation = "https://example.org/docs"\n')


def test_ci_wheel_report_must_match_commit_version_checksum_and_tests(tmp_path):
    wheel = tmp_path / "libephemeris-9.0.0-cp312-abi3-win_amd64.whl"
    wheel.write_bytes(b"wheel")
    sha = "a" * 40
    report = {
        "source_sha": sha, "version": "9.0.0", "wheel": wheel.name, "sha256": PREPARE.sha256(wheel),
        "tests": {"3.12": {"returncode": 0}, "3.14": {"returncode": 0}},
    }
    assert PREPARE.report_problems(report, wheel, sha, "9.0.0") == []
    assert len(PREPARE.report_problems(report, wheel, "b" * 40, "9.0.1")) == 2
    wheel.write_bytes(b"replaced")
    report["tests"]["3.14"]["returncode"] = 1
    del report["tests"]["3.12"]
    assert len(PREPARE.report_problems(report, wheel, sha, "9.0.0")) == 3
