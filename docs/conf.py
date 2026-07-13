# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
import tomllib
from pathlib import Path

from docutils import nodes
from docutils.parsers.rst import Directive
from docutils.statemachine import StringList

# Add the project root to the path so autodoc can find the modules
PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

# -- Project information -----------------------------------------------------
project = "libephemeris"
copyright = "2024-2026, Giacomo Battaglia"
author = "Giacomo Battaglia"

# pyproject.toml is the packaging source of truth, including in an uninstalled
# source checkout. Avoid a hard-coded fallback that can drift between releases.
with (PROJECT_ROOT / "pyproject.toml").open("rb") as pyproject_file:
    release = tomllib.load(pyproject_file)["project"]["version"]
version = ".".join(release.split(".")[:2])

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",  # Support for Google/NumPy style docstrings
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "myst_parser",  # Support for Markdown files
]

# Napoleon settings for Google/NumPy style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_use_ivar = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# Autodoc settings. The API reference derives its public names and signatures
# from libephemeris.__all__ and the runtime objects.
autodoc_default_options = {
    "member-order": "alphabetical",
    "exclude-members": "__weakref__",
}
autodoc_typehints = "description"
autodoc_typehints_description_target = "documented"

# Autosummary settings
autosummary_generate = True

# Intersphinx settings (link to external documentation)
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Source file suffixes
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
master_doc = "index"

# MyST parser settings
myst_enable_extensions = [
    "colon_fence",  # ::: directive syntax
    "deflist",  # Definition lists
    "tasklist",  # Task lists with checkboxes
]
myst_heading_anchors = 3  # Add anchors to headings up to level 3

# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"

# Theme options
html_theme_options = {
    "navigation_depth": 4,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "includehidden": True,
    "titles_only": False,
}

# -- Options for autodoc -----------------------------------------------------
# Document classes and functions
add_module_names = False


class PublicDataDirective(Directive):
    """Render every non-callable export in ``libephemeris.__all__``."""

    has_content = False

    def run(self):
        import libephemeris

        source = self.state.document.current_source
        rst_lines: list[str] = []

        for name in libephemeris.__all__:
            value = getattr(libephemeris, name)
            if callable(value):
                continue

            rst_lines.append(f".. py:data:: {name}")
            rst_lines.append(f"   :type: {type(value).__name__}")
            if isinstance(value, (str, int, float, bool, type(None))):
                rst_lines.append(f"   :value: {value!r}")
            rst_lines.append("")

        container = nodes.container()
        self.state.nested_parse(
            StringList(rst_lines, source=source), self.content_offset, container
        )
        return [container]


def setup(app):
    app.add_directive("public-data", PublicDataDirective)
    return {"parallel_read_safe": True}


# -- Options for man pages ---------------------------------------------------
man_pages = [(master_doc, "libephemeris", "libephemeris Documentation", [author], 1)]
