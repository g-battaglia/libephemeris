"""Explicit optional test data root; never changes engine configuration.

LEB identities follow the owner-selected current comparison distribution.
The bundled published image and non-LEB pins retain their own identities.
"""
import hashlib
import os
from pathlib import Path

# Product test input pin; independent of the private qualification harness.
CURRENT_CORE_SHA256 = 'b891de02f82b68346046da67389fd3713aa81ab299287c632cb20cc396397964'

# All product tests authenticate the same optional local input images. These
# identities are fixtures, not a dependency on a comparison engine or harness.
LEB_HASHES = {
    'base': {
        'core': CURRENT_CORE_SHA256,
        'asteroids': '9364f821f09076db3cd3cfc428379e136a0eea80917af692681be405601652d8',
        'apogee': 'f76ca7354d18e21e9f4e96f859149f29b0a2f0c19c19e305db1b7dba81f9c616',
        'exotics': 'e0b5bf4b14b10495c704d553557922f91f3b9db71117e21f1eb4fc0a134f4d1a',
    },
    'medium': {
        'core': 'b74970cc37f8d0316c6a6fc0c72a5690d2b1f38f65352590a70ecb6500c9572d',
        'asteroids': 'f7a3784a1e29a1fcc42ad9736f0333bc690535befea008e86a7b799a4a7652d9',
        'apogee': 'a2f8970cfd08150394556828862f497401f4b47fc7f75e1db7e9c63057324a3a',
        'exotics': 'c549f12032dd2d74ad547ddad43c8335f3fe975121b71d2c1b5f1c4b37ed16b5',
    },
    'extended': {
        'core': '2122b11a6b5d19178e54de41cb1d956e4e0fb54331f2c5b7e0dd8991ffb572cf',
        'asteroids': '2c9c9f999e22857bc24c0bbfd55e95192739f4984da56a8d19156d8c03e6907c',
        'apogee': '9426373afabdfaa1d4fa08c8eb580ff532508f4ff07aa7c4a9ef08ca7ddb193a',
        'exotics': '8b5094beaacd2a14807f60ddb6867922bdfa757624b6c7404b8fe02baa341d94',
    },
}

DATA_ROOT = Path(os.environ.get('LIBEPHEMERIS_TEST_DATA_DIR', str(Path.home() / '.libephemeris'))).expanduser().resolve()

# Numerical test inputs always use the owner-selected production root.
CORE_LEB = Path(os.environ.get('LIBEPHEMERIS_TEST_LEB2', str(DATA_ROOT / 'leb/base_core.leb2'))).expanduser()
# Manifest regressions require the exact published data-v4 distribution. A
# same-named owner-selected file must independently match its actual digest.
PUBLISHED_CORE_LEB = Path(__file__).resolve().parents[2] / 'core/data/leb2/base_core.leb2'


def is_published_core(path):
    """Check real file bytes against the public distribution pin, not its name."""
    from libephemeris.download import DATA_FILES
    with Path(path).open('rb') as stream:
        actual = hashlib.file_digest(stream, 'sha256').hexdigest()
    return actual == DATA_FILES['base_core.leb2']['sha256']


ADMITTED_CORE_HASHES = frozenset((
    '5d708bdbe3e799e0802ba575984e57a3c5e44720dbfa1b4a01cf826640e0cb82',
    CURRENT_CORE_SHA256,
))
