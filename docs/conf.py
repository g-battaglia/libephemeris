"""Build the current Python documentation against an installed native wheel.

Do not prepend the checkout to sys.path: that can load a stale local extension
and conceal a broken installation. Package metadata and runtime versions must
identify the same candidate before autodoc runs.
"""

from importlib import metadata
from pathlib import Path
import re
from urllib.parse import quote, urlsplit

import libephemeris

project = "LibEphemeris"
author = "Giacomo Battaglia"
copyright = "Giacomo Battaglia"
release = metadata.version("libephemeris")
if libephemeris.__version__ != release:
    raise RuntimeError("documentation runtime and installed metadata versions differ")
version = release

extensions = ["myst_parser", "sphinx.ext.autodoc"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
# Both books live beside docs/, so the source root is the repository root.
# An allowlist prevents worktrees, generated output and private files entering HTML.
root_doc = "docs/index"
include_patterns = ["docs/**", "beginners-manual/**"]
exclude_patterns = ["docs/_build/**", "**/.DS_Store"]
html_theme = "sphinx_rtd_theme"
html_theme_options = {"navigation_depth": 3}
myst_heading_anchors = 4
myst_enable_extensions = ["colon_fence", "deflist", "tasklist"]
# Native signatures/docstrings are authoritative; type stubs ship separately.
autodoc_typehints = "none"
add_module_names = False


def source_links(app, docname, source):
    """Resolve out-of-book links against the canonical source repository.

    GitHub-readable Markdown also links to root guides, source and numerical
    reports. Those are not Sphinx documents. Preserve internal book links for
    MyST, but pin other existing project files to the installed release's tag.
    No external page is fetched and no missing-link warning is suppressed.
    """
    root = Path(app.srcdir).resolve()
    book_roots = (root / "docs", root / "beginners-manual")
    origin = Path(app.env.doc2path(docname)).resolve().parent

    def replace(match):
        target = match.group(2)
        parsed = urlsplit(target)
        if parsed.scheme or parsed.netloc or not parsed.path:
            return match.group(0)
        destination = (origin / parsed.path).resolve()
        if any(destination.is_relative_to(book) for book in book_roots):
            return match.group(0)
        if not destination.is_relative_to(root) or not destination.exists():
            raise ValueError(f"invalid source link in {docname}: {target}")
        kind = "tree" if destination.is_dir() else "blob"
        relative = quote(destination.relative_to(root).as_posix())
        url = f"https://github.com/g-battaglia/libephemeris/{kind}/v{release}/{relative}"
        if parsed.fragment:
            url += "#" + parsed.fragment
        return match.group(1) + url + match.group(3)

    source[0] = re.sub(r"(!?\[[^\]]*\]\()([^\s)]+)(\))", replace, source[0])


def native_help_text(app, what, name, obj, options, lines):
    """Render native help verbatim rather than misparsing Markdown as RST.

    PyO3 help contains Markdown fences and mathematical notation intended for
    Python help() and Rust readers. Autodoc still indexes native signatures;
    a text block preserves every docstring word without altering runtime help.
    """
    if lines:
        lines[:] = [".. code-block:: text", ""] + ["    " + line for line in lines] + [""]


def setup(app):
    """Keep source links portable and native help faithful in the HTML book."""
    app.connect("source-read", source_links)
    app.connect("autodoc-process-docstring", native_help_text)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
