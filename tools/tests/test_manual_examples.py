"""The full bilingual books remain intact and their offline examples execute.

Network examples are syntax/signature checked, never downloaded by the test
suite. Set LIBEPHEMERIS_MANUAL_DATA_DIR explicitly to run additional-data cases.
The initial reader path is always tested with an empty data directory.
"""
from __future__ import annotations

import ast
import importlib.util
import inspect
import math
import os
from pathlib import Path
import re
import subprocess
import sys

import pytest

ephe = pytest.importorskip("libephemeris")
ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("manual_checks", ROOT / "tools/release/check_manual.py")
MANUAL = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MANUAL
SPEC.loader.exec_module(MANUAL)
ITEMS = MANUAL.inventory()
NUMBER = re.compile(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


def test_both_books_keep_all_chapters_recipes_and_teaching_depth():
    for language in ("it", "en"):
        directory = ROOT / "beginners-manual" / language
        chapters = sorted(directory.glob("[0-9]*.md"))
        assert [p.name[:2] for p in chapters] == [f"{i:02d}" for i in range(16)]
        # The recovered books contained about 38,000 words apiece. This guard
        # rejects replacement by a short tutorial, without enforcing padding.
        assert sum(len(p.read_text().split()) for p in chapters) >= 37000
        assert len(re.findall(r"^## (?:Recipe|Ricetta) \d+", chapters[-1].read_text(), re.M)) == 12
        for chapter in chapters:
            assert "(README.md)" in chapter.read_text()
            assert ("../it/" if language == "en" else "../en/") in chapter.read_text()
    assert len(ITEMS) == 302
    it = [item for item in ITEMS if item.path.parent.name == "it"]
    en = [item for item in ITEMS if item.path.parent.name == "en"]
    assert [(item.key, item.category) for item in it] == [(item.key, item.category) for item in en]


@pytest.mark.parametrize("example", ITEMS, ids=lambda e: f"{e.path.parent.name}-{e.chapter:02d}-{e.number}")
def test_every_python_block_compiles_and_uses_current_call_signatures(example):
    tree = ast.parse(example.code)
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and isinstance(node.func.value, ast.Name) and node.func.value.id == "ephem"):
            continue
        function = getattr(ephe, node.func.attr)
        # Dynamic expansion needs runtime checks, not an invented static shape.
        if any(isinstance(arg, ast.Starred) for arg in node.args) or any(k.arg is None for k in node.keywords):
            continue
        try:
            signature = inspect.signature(function)
        except (ValueError, TypeError):
            continue
        signature.bind(*[None for _ in node.args], **{key.arg: None for key in node.keywords})


def assert_transcript(actual: str, expected: str) -> None:
    """Check displayed values without requiring bit identity of rounded floats."""
    actual, expected = actual.strip(), expected.strip()
    assert NUMBER.sub("#", actual) == NUMBER.sub("#", expected)
    actual_values = [float(word) for word in NUMBER.findall(actual)]
    expected_values = [float(word) for word in NUMBER.findall(expected)]
    assert len(actual_values) == len(expected_values)
    for actual_value, expected_value in zip(actual_values, expected_values, strict=True):
        assert math.isclose(actual_value, expected_value, rel_tol=1e-9, abs_tol=1e-6)


@pytest.mark.parametrize("example", ITEMS, ids=lambda e: f"{e.path.parent.name}-{e.chapter:02d}-{e.number}")
def test_manual_runtime_examples(example):
    configured = os.environ.get("LIBEPHEMERIS_MANUAL_DATA_DIR")
    data = Path(configured).resolve() if configured else None
    result = MANUAL.run_example(example, data)
    if result["status"] == "not-run":
        pytest.skip(result["reason"])
    assert result["status"] == "passed", result.get("stderr")
    # External data may differ from the recorded local snapshot.
    if example.expected is not None and example.category != "local-data":
        assert_transcript(result["stdout"], example.expected)


@pytest.mark.parametrize("language", ["it", "en"])
def test_cookbook_phase_search_solves_moving_sun_and_stays_inside_year(language, tmp_path):
    example = next(item for item in ITEMS if item.path.parent.name == language and item.key == (15, 5))
    # These checks validate the teaching algorithm, not just a successful print.
    program = MANUAL.setup(example.path) + "\n" + example.code + '''
assert len(events) == 25
assert sum(target == 0 for _, target in events) == 13
assert all(2460310.5 <= epoch < 2460676.5 for epoch, _ in events)
assert all(a[0] < b[0] for a, b in zip(events, events[1:]))
april_new = [epoch for epoch, target in events if target == 0 and 2460408.5 < epoch < 2460410.5]
assert len(april_new) == 1
assert abs(april_new[0] - 2460409.264) < 0.01
'''
    env = {key: value for key, value in os.environ.items() if not key.startswith("LIBEPHEMERIS_")}
    env.update(LIBEPHEMERIS_DATA_DIR=str(tmp_path), LIBEPHEMERIS_NETWORK_POLICY="sealed",
               LIBEPHEMERIS_CONFIG=str(tmp_path / "none.toml"),
               LIBEPHEMERIS_ENV_FILE=str(tmp_path / "none.env"))
    result = subprocess.run([sys.executable, "-c", program], cwd=tmp_path, env=env,
                            capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr


def test_transcript_check_rejects_wrong_units_and_significant_values():
    assert_transcript("Longitude: 19.1404375°", "Longitude: 19.140438°")
    with pytest.raises(AssertionError):
        assert_transcript("Distance: 1 AU", "Distance: 1 km")
    with pytest.raises(AssertionError):
        assert_transcript("Longitude: 20°", "Longitude: 19°")
