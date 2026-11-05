"""Registry descriptions must resolve before the beta is merged into main.

Check tagged links against the local candidate without requiring the tag to
have been pushed yet. The release process verifies the remote tag separately.
"""

import importlib.util
from pathlib import Path
import tomllib
from urllib.parse import unquote, urlsplit

import pytest

ROOT = Path(__file__).resolve().parents[2]
TAG = "python-v4.0.0b4"
UNCHANGED_PACKAGE_TAG = "python-v4.0.0b1"
PROJECT_PATH = "/g-battaglia/libephemeris"
SPEC = importlib.util.spec_from_file_location("release_link_checker", ROOT / "tools/release/check_docs.py")
CHECKER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECKER)


def assert_pinned_destination(target, expected_tag=TAG):
    """Resolve only this project's source links, never fetch third-party pages."""
    parsed = urlsplit(target)
    assert parsed.scheme in {"https", "mailto"}, f"relative or insecure link: {target}"
    if parsed.netloc != "github.com" or not parsed.path.startswith(PROJECT_PATH + "/"):
        # A fragment on the repository homepage would follow main's v3 README.
        assert not (parsed.netloc == "github.com" and parsed.path == PROJECT_PATH and parsed.fragment)
        return
    relative = parsed.path.removeprefix(PROJECT_PATH + "/").split("/", 2)
    if relative[0] not in {"blob", "tree"}:
        return  # Issues and historical releases are not candidate source files.
    kind, reference, path = relative
    assert reference == expected_tag, f"unpinned candidate source link: {target}"
    destination = (ROOT / unquote(path)).resolve()
    assert destination.is_relative_to(ROOT)
    assert destination.exists(), f"missing tagged destination: {target}"
    assert destination.is_dir() == (kind == "tree")
    if destination.suffix == ".md" and parsed.fragment:
        assert unquote(parsed.fragment) in CHECKER.anchors(destination.read_text()), target


@pytest.mark.parametrize("name,tag", [
    ("README.md", TAG), ("python/README.md", TAG),
    ("core/README.md", TAG),
    ("core/data/README.md", UNCHANGED_PACKAGE_TAG),
    ("tools/leb_generation/README.md", UNCHANGED_PACKAGE_TAG),
])
def test_registry_readme_links_are_absolute_and_pin_candidate_source(name, tag):
    content = CHECKER.prose((ROOT / name).read_text())
    links = [match.group(1) for match in CHECKER.LINK.finditer(content)]
    assert links, f"no registry description links found: {name}"
    for target in links:
        assert_pinned_destination(target, tag)


def test_registry_documentation_metadata_pins_the_release():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    tools = tomllib.loads((ROOT / "tools/leb_generation/pyproject.toml").read_text())["project"]
    core = tomllib.loads((ROOT / "core/Cargo.toml").read_text())["package"]
    assert_pinned_destination(project["urls"]["Documentation"])
    assert_pinned_destination(project["urls"]["Changelog"])
    assert_pinned_destination(tools["urls"]["Documentation"], UNCHANGED_PACKAGE_TAG)
    assert core["documentation"] == f"https://docs.rs/libephemeris-core/{core['version']}"
