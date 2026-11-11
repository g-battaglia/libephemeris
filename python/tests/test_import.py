"""Smoke checks for the installed native package."""

from importlib import metadata

import libephemeris
from libephemeris import _core


def test_package_exports_native_version():
    assert libephemeris.__version__ == "4.0.0"
    assert libephemeris.__version__ == _core.__version__ == libephemeris.version
    assert libephemeris.__version__ == metadata.version("libephemeris")


def test_public_exports_exist():
    for name in libephemeris.__all__:
        assert hasattr(libephemeris, name)
