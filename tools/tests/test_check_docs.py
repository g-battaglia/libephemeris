"""Local documentation checks must catch broken links without network access."""

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("check_docs", ROOT / "tools/release/check_docs.py")
CHECKER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECKER)


def test_anchors_preserve_code_labels_and_number_duplicates():
    text = "# `calc_ut` — Units\n# Repeated\n# Repeated\n<a id='explicit'></a>\n"
    assert CHECKER.anchors(text) == {"calc_ut--units", "repeated", "repeated-1", "explicit"}


def test_code_fences_do_not_create_links_or_headings(tmp_path):
    page = tmp_path / "README.md"
    page.write_text("# Guide\n```md\n[not a link](missing.md)\n# hidden\n```\n")
    assert CHECKER.document_errors(page, tmp_path) == []
    assert CHECKER.anchors(page.read_text()) == {"guide"}


def test_local_destinations_and_fragments_are_checked(tmp_path):
    page = tmp_path / "README.md"
    other = tmp_path / "other.md"
    other.write_text("# Actual heading\n")
    page.write_text("[good](other.md#actual-heading)\n[bad](other.md#missing)\n[absent](absent.md)\n")
    errors = CHECKER.document_errors(page, tmp_path)
    assert len(errors) == 2
    assert any("missing anchor" in error for error in errors)
    assert any("missing destination" in error for error in errors)


def test_external_urls_are_not_fetched(tmp_path):
    page = tmp_path / "README.md"
    page.write_text("[remote](https://example.invalid/page)\n[email](mailto:person@example.invalid)\n")
    assert CHECKER.document_errors(page, tmp_path) == []


def test_links_cannot_escape_or_follow_symlinks_outside_tree(tmp_path):
    root = tmp_path / "project"
    root.mkdir()
    secret = tmp_path / "outside.md"
    secret.write_text("# Outside")
    (root / "alias.md").symlink_to(secret)
    page = root / "README.md"
    page.write_text("[outside](../outside.md)\n[alias](alias.md)\n")
    assert len(CHECKER.document_errors(page, root)) == 2


def test_generated_output_and_process_archives_are_not_traversed(tmp_path):
    for relative in ("docs/_build", "docs/archive"):
        folder = tmp_path / relative
        folder.mkdir(parents=True)
        (folder / "old.md").write_text("[bad](missing.md)")
    (tmp_path / "README.md").write_text("# Current")
    assert CHECKER.check(tmp_path) == []


def test_source_move_rewrites_links_not_labels_fences_or_external_urls(tmp_path):
    (tmp_path / 'docs').mkdir()
    destination = tmp_path / 'core/src'
    destination.mkdir(parents=True)
    (destination / 'lib.rs').write_text('// current source')
    text = ('[../src/lib.rs](../src/lib.rs#L1)\n'
            '[external](https://example.invalid/src/lib.rs)\n'
            '```md\n[example](../src/lib.rs)\n```\n')
    page = tmp_path / 'docs/guide.md'
    rendered = CHECKER.relocate_links(text, page, tmp_path, [('src', 'core/src')])
    assert '[../src/lib.rs](../core/src/lib.rs#L1)' in rendered
    assert '[external](https://example.invalid/src/lib.rs)' in rendered
    assert '```md\n[example](../src/lib.rs)\n```' in rendered
    assert CHECKER.relocate_links(rendered, page, tmp_path, [('src', 'core/src')]) == rendered


def test_public_docs_cannot_depend_on_migration_history_or_process_files(tmp_path):
    (tmp_path / "PROGRESS.md").write_text("# Working notes")
    page = tmp_path / "README.md"
    page.write_text(
        "[old](https://github.com/g-battaglia/libephemeris-rust/pull/1)\n"
        "[checkpoint](PROGRESS.md)\n"
    )
    errors = CHECKER.document_errors(page, tmp_path)
    assert len(errors) == 2
    assert any("retired migration repository" in error for error in errors)
    assert any("process file" in error for error in errors)
