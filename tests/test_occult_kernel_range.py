"""Regression guard for the occultation-search kernel range clamp.

DE441 (the ``extended`` precision tier) splits every body into two SPK
segments at JD 2440432.5 (1969-07-30). The forward search bound in
``lun_occult_when_glob`` must be the *maximum* segment ``end_jd`` (the
kernel's true end, ~JD 8000016.5 / year 17191), not the *minimum* (the 1969
split point). Using the minimum clamped every global lunar-occultation
search to 1969, so any search starting after that date raised "no
occultation within 20 years" even where events exist.

These tests pin the bound via ``_kernel_jd_bounds`` with synthetic segments,
so they run fast and deterministically without loading the 3.1 GB DE441 file.
"""

from __future__ import annotations

import types

from libephemeris.eclipse import _kernel_jd_bounds


def _fake_eph(segment_bounds):
    """Build a stand-in ephemeris exposing ``.segments[*].spk_segment.{start,end}_jd``."""
    segments = [
        types.SimpleNamespace(
            spk_segment=types.SimpleNamespace(start_jd=start, end_jd=end)
        )
        for start, end in segment_bounds
    ]
    return types.SimpleNamespace(segments=segments)


def test_split_kernel_uses_max_end_jd():
    """DE441-like kernel: bound must be the latest end_jd, not the 1969 split."""
    eph = _fake_eph(
        [
            (-3100015.5, 2440432.5),  # body A, pre-1969 segment
            (2440432.5, 8000016.5),  # body A, post-1969 segment
            (-3100015.5, 2440432.5),  # body B, pre-1969 segment
            (2440432.5, 8000016.5),  # body B, post-1969 segment
        ]
    )
    lo, hi = _kernel_jd_bounds(eph)
    assert lo == -3100015.5
    assert hi == 8000016.5  # NOT 2440432.5 (the 1969 split point)


def test_single_segment_kernel_unchanged():
    """DE440-like kernel (no split): bounds are just that segment's range."""
    eph = _fake_eph([(2287184.5, 2688976.5)])
    assert _kernel_jd_bounds(eph) == (2287184.5, 2688976.5)
