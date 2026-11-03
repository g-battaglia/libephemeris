"""Keep build and distribution paths aligned with the core/Python split."""

from pathlib import Path
import tomllib

ROOT = Path(__file__).resolve().parents[2]


def test_workspace_uses_core_directory():
    workspace = tomllib.loads((ROOT / 'Cargo.toml').read_text())['workspace']
    expected = {'core', 'core/data', 'python'}
    # Maturin retains only the binding as an explicit sdist member. Cargo
    # discovers both core crates through their same-workspace path dependencies.
    assert 'python' in workspace['members']
    assert set(workspace['members']).issubset(expected)
    assert workspace['default-members'] == ['core']
    for member in expected:
        assert (ROOT / member / 'Cargo.toml').is_file()
    assert not (ROOT / 'rust').exists()


def test_bindings_resolve_actual_engine_manifest():
    binding = tomllib.loads((ROOT / 'python/Cargo.toml').read_text())
    dependency = binding['dependencies']['libephemeris-core']
    engine = (ROOT / 'python' / dependency['path']).resolve()
    assert engine == ROOT / 'core'
    manifest = tomllib.loads((engine / 'Cargo.toml').read_text())
    assert manifest['package']['name'] == 'libephemeris-core'


def test_owner_retired_process_documents_are_absent():
    # The owner removed these development directives from the product tree.
    # Check the actual checkout too, not merely source-package export filters.
    for name in ('AGENTS.md', 'PLAN.md', 'PROGRESS.md'):
        assert not (ROOT / name).exists()


def test_license_packaging_uses_actual_core_notices():
    project = tomllib.loads((ROOT / 'pyproject.toml').read_text())['project']
    assert 'core/licenses/**/*' in project['license-files']
    assert (ROOT / 'core/licenses/index.json').is_file()
