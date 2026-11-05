#!/usr/bin/env python3
"""Run the maintained manual examples with an installed wheel and sealed network.

Every Python fence is inventoried. Optional acquisition examples are compiled,
not executed: this command never authorizes downloads. Complete examples run in
fresh processes with the setup printed in their language's README. Explicit
continuations also execute the prerequisite example, with its output captured
separately. --refresh-output updates only adjacent output fences after success;
it never changes code, prose, or the output of an unexecuted example.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
PYTHON_BLOCK = re.compile(r"^```python\n(.*?)^```[ \t]*$", re.M | re.S)
OUTPUT_BLOCK = re.compile(r"\n\n```\n(.*?)^```[ \t]*$", re.M | re.S)
# These are instructions to acquire resources, not calculation examples.
ONLINE = {(0, 1), (2, 6), (2, 7), (4, 1), (12, 3), (12, 4), (12, 6),
          (12, 9), (14, 6), (14, 11), (14, 18), (14, 19)}
# Extra coefficient groups or an external JPL kernel, declared in the book index.
LOCAL_DATA = {(6, 3), (6, 4), (6, 6), (12, 1), (12, 2), (12, 5), (12, 8),
              (14, 5), (14, 7), (14, 9), (14, 16), (15, 1)}
CONTINUATIONS = {(5, 7): 6, **{(8, i): 5 for i in range(7, 14)},
                 (9, 2): 1, (13, 4): 3, (13, 5): 3}


@dataclass(frozen=True)
class Example:
    """One source fence, its resource class and an optional printed transcript."""
    path: Path
    chapter: int
    number: int
    code: str
    expected: str | None
    output_span: tuple[int, int] | None

    @property
    def key(self) -> tuple[int, int]:
        return self.chapter, self.number

    @property
    def category(self) -> str:
        if self.key in ONLINE:
            return "online"
        if self.key in LOCAL_DATA:
            return "local-data"
        return "continuation" if self.key in CONTINUATIONS else "bundled"


def examples(path: Path) -> list[Example]:
    """Extract fences and only their immediately adjacent, unlabelled output."""
    text = path.read_text(encoding="utf-8")
    result = []
    for number, match in enumerate(PYTHON_BLOCK.finditer(text), 1):
        output = OUTPUT_BLOCK.match(text, match.end())
        result.append(Example(path, int(path.name[:2]), number, match[1],
                              output[1] if output else None,
                              output.span(1) if output else None))
    return result


def inventory(root: Path = ROOT, language: str | None = None) -> list[Example]:
    """Keep chapter order and the two translations visible to callers/tests."""
    languages = [language] if language else ["it", "en"]
    return [example for lang in languages
            for path in sorted((root / "beginners-manual" / lang).glob("[0-9]*.md"))
            for example in examples(path)]


def setup(path: Path) -> str:
    """Use the reader-visible setup, not an undisclosed test-only configuration."""
    return PYTHON_BLOCK.search((path.parent / "README.md").read_text())[1]


def run_example(example: Example, data_dir: Path | None = None) -> dict:
    """Run one isolated example; compile and explicitly report non-runtime cases."""
    record = {"path": example.path.relative_to(ROOT).as_posix(),
              "block": example.number, "category": example.category}
    compile(example.code, f"{example.path}:{example.number}", "exec")
    if example.category == "online":
        return record | {"status": "not-run", "reason": example.category}
    if example.category == "local-data" and data_dir is None:
        return record | {"status": "not-run", "reason": "requires explicit --data-dir"}

    prerequisite = ""
    if example.key in CONTINUATIONS:
        previous = examples(example.path)[CONTINUATIONS[example.key] - 1]
        prerequisite = previous.code
    # Capture prerequisite output independently, so only this fence's transcript
    # is compared. The variables remain in the same namespace, as in the book.
    program = (setup(example.path) + "\nimport contextlib, io\n"
               "with contextlib.redirect_stdout(io.StringIO()):\n"
               f"    exec({prerequisite!r})\n" + example.code)
    with tempfile.TemporaryDirectory(prefix="libephemeris-manual-") as temporary:
        environment = {key: value for key, value in os.environ.items()
                       if not key.startswith("LIBEPHEMERIS_")}
        # Bundled exercises always use an empty data directory, even when the
        # caller supplied external assets for other examples in the same run.
        selected = data_dir.resolve() if example.category == "local-data" else Path(temporary)
        environment.update(LIBEPHEMERIS_DATA_DIR=str(selected),
                           LIBEPHEMERIS_CONFIG=f"{temporary}/absent.toml",
                           LIBEPHEMERIS_ENV_FILE=f"{temporary}/absent.env",
                           LIBEPHEMERIS_NETWORK_POLICY="sealed")
        try:
            process = subprocess.run([sys.executable, "-c", program], cwd=temporary,
                                     env=environment, text=True, capture_output=True,
                                     timeout=120)
        except subprocess.TimeoutExpired:
            return record | {"status": "failed", "stderr": "120-second timeout"}
    return record | {"status": "passed" if process.returncode == 0 else "failed",
                     "stdout": process.stdout, "stderr": process.stderr}


def refresh_outputs(items: list[Example], results: list[dict]) -> None:
    """Replace only successfully recomputed transcript spans, back to front."""
    changes: dict[Path, list[tuple[tuple[int, int], str]]] = {}
    for example, result in zip(items, results, strict=True):
        if result["status"] == "passed" and example.output_span is not None:
            changes.setdefault(example.path, []).append(
                (example.output_span, result["stdout"].rstrip() + "\n"))
    for path, replacements in changes.items():
        text = path.read_text(encoding="utf-8")
        for (start, end), output in sorted(replacements, reverse=True):
            text = text[:start] + output + text[end:]
        path.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--language", choices=("it", "en"))
    parser.add_argument("--data-dir", type=Path,
                        help="existing read-only assets for explicitly local-data examples")
    parser.add_argument("--refresh-output", action="store_true")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    if args.data_dir is not None and not args.data_dir.is_dir():
        parser.error("--data-dir must name an existing directory")
    items = inventory(language=args.language)
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(lambda item: run_example(item, args.data_dir), items))
    for result in results:
        print(f"{result['status']:7} {result['path']}:{result['block']} ({result['category']})")
        if result["status"] == "failed":
            print(result["stderr"])
    if args.refresh_output:
        refresh_outputs(items, results)
    if args.report:
        args.report.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n")
    counts = {status: sum(row["status"] == status for row in results)
              for status in ("passed", "failed", "not-run")}
    print(json.dumps(counts))
    return int(counts["failed"] > 0)


if __name__ == "__main__":
    raise SystemExit(main())
