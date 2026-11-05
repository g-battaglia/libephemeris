"""Stored bounds remain available when a dated query cannot use them."""

from data_paths import CORE_LEB, DATA_ROOT, is_published_core

import json
import os
from pathlib import Path
import subprocess
import sys

import libephemeris as e
from libephemeris.inventory import get_reader_body_coverage
import pytest

FIXTURE = CORE_LEB


def test_explicit_native_reader_reports_out_of_range_bounds_without_containment():
    e.set_calc_mode("leb")
    e.set_precision_tier("base")
    e.set_leb_file(str(FIXTURE))
    reader = e.get_leb_reader()
    known = get_reader_body_coverage(reader, e.SUN)
    assert known is not None
    assert known.reviewed is is_published_core(FIXTURE)
    for jd in (known.jd_start - 1, known.jd_end + 1):
        record = get_reader_body_coverage(reader, e.SUN, jd)
        assert record.to_dict() == known.to_dict()
        assert not record.contains(jd)
        # The explicit reader retains its stored file label. The public runtime
        # answers an uncovered epoch with an unlabelled diagnostic envelope;
        # it must not pretend that a file served that date.
        expected_union = dict(record.to_dict(), data_file=None, group=None)
        assert e.get_body_coverage(e.SUN, jd).to_dict() == expected_union
    assert get_reader_body_coverage(reader, -999999) is None
    assert get_reader_body_coverage(reader, e.CUPIDO) is None
    reader.close()
    with pytest.raises(e.ConfigurationError):
        get_reader_body_coverage(reader, e.SUN)
