"""Help text exposed by the installed native package, not a source copy.

Presence checks prevent accidental loss of Python documentation. They supplement,
not replace, human review of units, domains, clocks, outputs and error behavior.
"""
import inspect

import libephemeris as engine
import pytest


FUNCTIONS = [name for name in engine.__all__ if inspect.isbuiltin(getattr(engine, name))]
CLASSES = [name for name in engine.__all__ if isinstance(getattr(engine, name), type)]


def public_members():
    """Visit each native class once when several root names alias its type."""
    seen = set()
    for exported_name in CLASSES:
        cls = getattr(engine, exported_name)
        if id(cls) in seen:
            continue
        seen.add(id(cls))
        for name, value in vars(cls).items():
            if name.startswith('_'):
                continue
            if (inspect.ismethoddescriptor(value) or inspect.isdatadescriptor(value)
                    or inspect.isbuiltin(value) or isinstance(value, (classmethod, staticmethod))):
                yield cls.__name__ + '.' + name, value


MEMBERS = list(public_members())


def test_documentation_checks_cover_real_native_entry_points():
    assert {'calc', 'calc_ut', 'julday', 'houses_ex', 'fixstar_ut', 'rise_trans'} <= set(FUNCTIONS)
    assert {'EphemerisContext', 'OrbitalElements', 'StateVector'} <= set(CLASSES)
    assert MEMBERS


@pytest.mark.parametrize('name', FUNCTIONS)
def test_native_public_function_has_help_text(name):
    text = getattr(engine, name).__doc__
    assert isinstance(text, str) and text.strip(), name


@pytest.mark.parametrize('name', CLASSES)
def test_native_public_class_has_help_text(name):
    text = getattr(engine, name).__doc__
    assert isinstance(text, str) and text.strip(), name


@pytest.mark.parametrize('name,member', MEMBERS, ids=[name for name, _ in MEMBERS])
def test_native_public_member_has_help_text(name, member):
    text = member.__doc__
    assert isinstance(text, str) and text.strip(), name
