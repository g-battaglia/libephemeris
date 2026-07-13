# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Full-mock branch-coverage tests for ``libephemeris.iers_data``.

This module exercises the network/file-IO/optional-dependency paths of the
IERS data backend without ever touching the real network. All HTTP access is
mocked at the ``urllib.request.urlopen`` seam (or by monkeypatching the
``_download_file`` helper), file ages are forced via ``time``/``os.path``
patches, and ``time.sleep`` is neutralised. The goal is to drive every
reachable line and branch to 100% coverage.
"""

from __future__ import annotations

import os
from contextlib import contextmanager
from unittest import mock

import pytest

from libephemeris import iers_data


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _isolated_iers_state(tmp_path):
    """Isolate cache directory and in-memory state per test."""
    iers_data.clear_iers_cache()
    iers_data.set_iers_cache_dir(str(tmp_path / "cache"))
    iers_data.set_iers_auto_download(False)
    yield
    iers_data.clear_iers_cache()
    iers_data.set_iers_cache_dir(None)
    iers_data.set_iers_auto_download(None)


@contextmanager
def _fake_urlopen(content: bytes):
    """Patch urllib.request.urlopen to return a context manager yielding bytes.

    Mimics the http.client response contract used by the shared
    download.download_file() pipeline: .headers and chunked .read(size).
    """

    class _Resp:
        headers = {"Content-Length": str(len(content))}

        def __init__(self):
            self._data = content

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def read(self, size=-1):
            data, self._data = self._data, b""
            return data

    with mock.patch("urllib.request.urlopen", return_value=_Resp()) as m:
        yield m


# ---------------------------------------------------------------------------
# Configuration getters/setters
# ---------------------------------------------------------------------------


class TestConfig:
    """Cover the configuration getters and the auto-download resolution."""

    def test_get_cache_dir_override(self):
        """get_iers_cache_dir returns the explicit override (line 194)."""
        iers_data.set_iers_cache_dir("/tmp/some-dir")
        assert iers_data.get_iers_cache_dir() == "/tmp/some-dir"

    def test_auto_download_explicit_true(self):
        """Explicit True short-circuits before env/toml lookup."""
        iers_data.set_iers_auto_download(True)
        assert iers_data.get_iers_auto_download() is True

    def test_auto_download_env_enabled(self, monkeypatch):
        """Env var resolving to an enabled token returns True (line 231)."""
        iers_data.set_iers_auto_download(None)
        monkeypatch.setenv(iers_data._IERS_AUTO_DOWNLOAD_ENV_VAR, "yes")
        assert iers_data.get_iers_auto_download() is True

    def test_auto_download_env_disabled(self, monkeypatch):
        """Env var resolving to a disabled token returns False (line 233)."""
        iers_data.set_iers_auto_download(None)
        monkeypatch.setenv(iers_data._IERS_AUTO_DOWNLOAD_ENV_VAR, "off")
        assert iers_data.get_iers_auto_download() is False

    def test_auto_download_toml_true(self, monkeypatch):
        """No env var -> TOML fallback returns its value (line 240)."""
        iers_data.set_iers_auto_download(None)
        monkeypatch.delenv(iers_data._IERS_AUTO_DOWNLOAD_ENV_VAR, raising=False)
        with mock.patch(
            "libephemeris._config_toml.get_bool", return_value=True
        ) as toml_bool:
            assert iers_data.get_iers_auto_download() is True
        toml_bool.assert_called_once_with("iers_auto_download")

    def test_auto_download_default_false(self, monkeypatch):
        """No env var and no TOML value -> default False."""
        iers_data.set_iers_auto_download(None)
        monkeypatch.delenv(iers_data._IERS_AUTO_DOWNLOAD_ENV_VAR, raising=False)
        with mock.patch("libephemeris._config_toml.get_bool", return_value=None):
            assert iers_data.get_iers_auto_download() is False


# ---------------------------------------------------------------------------
# Cache directory creation
# ---------------------------------------------------------------------------


class TestCacheDirCreation:
    """Cover the makedirs branch in _get_cache_dir (line 260)."""

    def test_get_cache_dir_creates_missing_dir(self, tmp_path):
        """A non-existent cache dir is created on demand."""
        target = tmp_path / "fresh_cache"
        iers_data.set_iers_cache_dir(str(target))
        assert not target.exists()
        result = iers_data._get_cache_dir()
        assert os.path.isdir(result)
        assert os.path.abspath(str(target)) == result

    def test_get_cache_dir_default_location(self, tmp_path):
        """No override -> falls back to the state data dir (line 257)."""
        iers_data.set_iers_cache_dir(None)
        with mock.patch(
            "libephemeris.state._get_data_dir", return_value=str(tmp_path / "data")
        ):
            result = iers_data._get_cache_dir()
        assert result.endswith(iers_data.DEFAULT_CACHE_DIR)
        assert os.path.isdir(result)


# ---------------------------------------------------------------------------
# _download_file
# ---------------------------------------------------------------------------


class TestDownloadFile:
    """Cover _download_file success (incl. makedirs) and failure arms."""

    def test_download_success_creates_dirs(self, tmp_path):
        """Success path writes atomically and creates the parent dir (line 309)."""
        out = tmp_path / "nested" / "newdir" / "out.dat"
        with _fake_urlopen(b"payload-bytes"):
            ok = iers_data._download_file("https://example/x", str(out))
        assert ok is True
        assert out.read_bytes() == b"payload-bytes"
        assert not (tmp_path / "nested" / "newdir" / "out.dat.tmp").exists()

    def test_download_failure_url_error(self, tmp_path):
        """URLError -> returns False (except arm at line 317-318)."""
        import urllib.error

        out = tmp_path / "out.dat"
        with mock.patch(
            "urllib.request.urlopen", side_effect=urllib.error.URLError("boom")
        ):
            ok = iers_data._download_file("https://example/x", str(out))
        assert ok is False
        assert not out.exists()

    def test_download_failure_timeout(self, tmp_path):
        """TimeoutError is also caught -> returns False."""
        out = tmp_path / "out.dat"
        with mock.patch("urllib.request.urlopen", side_effect=TimeoutError("slow")):
            ok = iers_data._download_file("https://example/x", str(out))
        assert ok is False


# ---------------------------------------------------------------------------
# download_* dispatchers (cache-fresh / primary / backup / both-fail)
# ---------------------------------------------------------------------------


class TestDownloadDispatchers:
    """Cover the three download_* dispatchers' branch structure."""

    @pytest.fixture(autouse=True)
    def _no_sleep(self):
        """Neutralise any retry sleeps."""
        with mock.patch("time.sleep", lambda *a, **k: None):
            yield

    # ---- finals --------------------------------------------------------

    def test_finals_cache_fresh_skip(self):
        """A fresh cached file is returned without downloading (lines 343-345)."""
        path = iers_data._get_finals_cache_path()
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            f.write("cached")
        with (
            mock.patch("time.time", return_value=1_000_000.0),
            mock.patch("os.path.getmtime", return_value=1_000_000.0),
            mock.patch.object(iers_data, "_download_file") as dl,
        ):
            result = iers_data.download_iers_finals(force=False)
        assert result == path
        dl.assert_not_called()

    def test_finals_stale_triggers_download(self):
        """A stale finals cache falls through to download (344->347)."""
        path = iers_data._get_finals_cache_path()
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            f.write("stale")
        with (
            mock.patch("time.time", return_value=10_000_000.0),
            mock.patch("os.path.getmtime", return_value=0.0),
            mock.patch.object(iers_data, "_download_file", return_value=True) as dl,
        ):
            result = iers_data.download_iers_finals(force=False)
        assert result == path
        assert dl.call_count == 1

    def test_leap_stale_triggers_download(self):
        """A stale leap-seconds cache falls through to download (390->393)."""
        path = iers_data._get_leap_seconds_cache_path()
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            f.write("stale")
        with (
            mock.patch("time.time", return_value=10_000_000.0),
            mock.patch("os.path.getmtime", return_value=0.0),
            mock.patch.object(iers_data, "_download_file", return_value=True) as dl,
        ):
            result = iers_data.download_leap_seconds(force=False)
        assert result == path
        assert dl.call_count == 1

    def test_finals_primary_success(self):
        """Primary URL succeeds (lines 352-353)."""
        path = iers_data._get_finals_cache_path()
        with mock.patch.object(iers_data, "_download_file", return_value=True) as dl:
            result = iers_data.download_iers_finals(force=True)
        assert result == path
        assert dl.call_count == 1

    def test_finals_backup_success(self):
        """Primary fails, backup succeeds."""
        path = iers_data._get_finals_cache_path()
        with mock.patch.object(
            iers_data, "_download_file", side_effect=[False, True]
        ) as dl:
            result = iers_data.download_iers_finals(force=True)
        assert result == path
        assert dl.call_count == 2

    def test_finals_both_fail_raises(self):
        """Both sources fail -> ConnectionError (line 361)."""
        with mock.patch.object(iers_data, "_download_file", side_effect=[False, False]):
            with pytest.raises(ConnectionError):
                iers_data.download_iers_finals(force=True)

    # ---- leap seconds --------------------------------------------------

    def test_leap_cache_fresh_skip(self):
        """Fresh leap-seconds cache is returned without download (389-391)."""
        path = iers_data._get_leap_seconds_cache_path()
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            f.write("cached")
        with (
            mock.patch("time.time", return_value=1_000_000.0),
            mock.patch("os.path.getmtime", return_value=1_000_000.0),
            mock.patch.object(iers_data, "_download_file") as dl,
        ):
            result = iers_data.download_leap_seconds(force=False)
        assert result == path
        dl.assert_not_called()

    def test_leap_backup_success(self):
        """Primary fails, backup succeeds (lines 402-405)."""
        path = iers_data._get_leap_seconds_cache_path()
        with mock.patch.object(
            iers_data, "_download_file", side_effect=[False, True]
        ) as dl:
            result = iers_data.download_leap_seconds(force=True)
        assert result == path
        assert dl.call_count == 2

    def test_leap_both_fail_raises(self):
        """Both leap-second sources fail -> ConnectionError (line 407)."""
        with mock.patch.object(iers_data, "_download_file", side_effect=[False, False]):
            with pytest.raises(ConnectionError):
                iers_data.download_leap_seconds(force=True)

    # ---- delta t -------------------------------------------------------

    def test_delta_t_cache_fresh_skip(self):
        """Fresh delta-t cache returned without download (444->447 fall-through)."""
        path = iers_data._get_delta_t_cache_path()
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            f.write("cached")
        with (
            mock.patch("time.time", return_value=1_000_000.0),
            mock.patch("os.path.getmtime", return_value=1_000_000.0),
            mock.patch.object(iers_data, "_download_file") as dl,
        ):
            result = iers_data.download_delta_t_data(force=False)
        assert result == path
        dl.assert_not_called()

    def test_delta_t_stale_triggers_download(self):
        """A stale cached file falls through to the download path (444->447)."""
        path = iers_data._get_delta_t_cache_path()
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            f.write("stale")
        # mtime far in the past relative to time.time -> age > max -> download
        with (
            mock.patch("time.time", return_value=10_000_000.0),
            mock.patch("os.path.getmtime", return_value=0.0),
            mock.patch.object(iers_data, "_download_file", return_value=True) as dl,
        ):
            result = iers_data.download_delta_t_data(force=False)
        assert result == path
        assert dl.call_count == 1

    def test_delta_t_backup_success(self):
        """Primary fails, backup succeeds (lines 456-459)."""
        path = iers_data._get_delta_t_cache_path()
        with mock.patch.object(
            iers_data, "_download_file", side_effect=[False, True]
        ) as dl:
            result = iers_data.download_delta_t_data(force=True)
        assert result == path
        assert dl.call_count == 2

    def test_delta_t_both_fail_raises(self):
        """Both delta-t sources fail -> ConnectionError (line 461)."""
        with mock.patch.object(iers_data, "_download_file", side_effect=[False, False]):
            with pytest.raises(ConnectionError):
                iers_data.download_delta_t_data(force=True)


# ---------------------------------------------------------------------------
# _parse_finals_data
# ---------------------------------------------------------------------------


class TestParseFinalsData:
    """Cover the finals2000A fixed-width parser branches."""

    def test_parse_finals_full_and_edge_cases(self, tmp_path):
        """Exercise short lines, empty MJD, missing flag, error parsing, skips."""
        # Build lines exercising distinct branches:
        # 1) too-short line (<68) -> continue (line 497)
        short = "73 1 1"
        # 2) empty MJD field -> continue (line 515): cols 0-6 date, 7-15 mjd blank
        #    Provide a line long enough (>=68) but with blank MJD slice.
        empty_mjd = "730101 " + " " * 8 + " " * 60

        # 3) full observed row with UT1-UTC and a valid error -> stored,
        #    is_prediction False (flag 'I'), error parsed (lines 533-535, store).
        # Column layout (Python slices):
        #   [0:2] yy, [2:4] mm, [4:6] dd, [7:15] mjd, col57 flag, [58:68] ut1,
        #   [68:78] err
        def make_row(yy, mm, dd, mjd, flag, ut1, err):
            line = list(" " * 80)
            line[0:2] = list(yy)
            line[2:4] = list(mm)
            line[4:6] = list(dd)
            line[7:15] = list(mjd.rjust(8))
            line[57] = flag
            ut1s = ut1.rjust(10)
            line[58:68] = list(ut1s)
            errs = err.rjust(10)
            line[68:78] = list(errs)
            return "".join(line)

        observed = make_row("00", "01", "01", "51544", "I", "0.355000", "0.000020")
        predicted = make_row("99", "12", "31", "51543", "P", "0.357000", "0.000030")
        # 4) row with a non-numeric error field -> ValueError caught (538-539),
        #    err stays None but the point is still stored.
        bad_err = make_row("01", "01", "01", "51910", "I", "0.350000", "garbage..")
        # 5) row with empty UT1 field -> continue (line 521)
        empty_ut1 = make_row("01", "06", "01", "52061", "I", "          ", "0.0001")
        # 6) row whose year/month is garbage -> ValueError (line 551-553)
        malformed = "XX01 1  51544.0" + " " * 60

        content = "\n".join(
            [short, empty_mjd, observed, predicted, bad_err, empty_ut1, malformed]
        )
        fp = tmp_path / "finals.data"
        fp.write_text(content + "\n")

        data = iers_data._parse_finals_data(str(fp))
        # observed, predicted, bad_err -> 3 stored points
        assert len(data) == 3
        assert data[51544.0].is_prediction is False
        assert data[51544.0].ut1_utc_err == pytest.approx(0.000020)
        assert data[51543.0].is_prediction is True
        assert data[51910.0].ut1_utc_err is None

    def test_parse_finals_no_error_column(self, tmp_path):
        """A line >=68 but <78 chars: error column absent (533->541 skip)."""
        line = list(" " * 70)
        line[0:2] = list("00")
        line[2:4] = list("01")
        line[4:6] = list("01")
        line[7:15] = list("51544".rjust(8))
        line[57] = "I"
        line[58:68] = list("0.355000".rjust(10))
        content = "".join(line)
        fp = tmp_path / "finals_noerr.data"
        fp.write_text(content + "\n")
        data = iers_data._parse_finals_data(str(fp))
        assert len(data) == 1
        assert data[51544.0].ut1_utc_err is None

    def test_parse_finals_empty_error_field(self, tmp_path):
        """Error column present but blank -> err_str falsy (535->541 skip)."""
        line = list(" " * 80)
        line[0:2] = list("00")
        line[2:4] = list("01")
        line[4:6] = list("01")
        line[7:15] = list("51544".rjust(8))
        line[57] = "I"
        line[58:68] = list("0.355000".rjust(10))
        line[68:78] = list(" " * 10)
        content = "".join(line)
        fp = tmp_path / "finals_blankerr.data"
        fp.write_text(content + "\n")
        data = iers_data._parse_finals_data(str(fp))
        assert len(data) == 1
        assert data[51544.0].ut1_utc_err is None

    def test_parse_finals_no_flag_column(self, tmp_path):
        """Line long enough for UT1 but short of flag col -> 527->532 skip.

        UT1 lives at [58:68]; the prediction flag at index 57. A line that
        reaches col 68 necessarily includes index 57, so to take the
        ``len(line) > 57`` False arm we craft a line that has the UT1 slice
        but is built so index 57 is missing is impossible; instead this
        covers the inverse: ensure the flag branch is taken with default.
        """
        # This row keeps is_prediction default-False via a non-'P' flag.
        line = list(" " * 80)
        line[0:2] = list("85")
        line[2:4] = list("07")
        line[4:6] = list("01")
        line[7:15] = list("46247".rjust(8))
        line[57] = " "
        line[58:68] = list("0.100000".rjust(10))
        content = "".join(line)
        fp = tmp_path / "finals_flag.data"
        fp.write_text(content + "\n")
        data = iers_data._parse_finals_data(str(fp))
        assert len(data) == 1
        assert data[46247.0].is_prediction is False


# ---------------------------------------------------------------------------
# _parse_leap_seconds (IERS + IETF + fallback)
# ---------------------------------------------------------------------------


class TestParseLeapSeconds:
    """Cover the leap-seconds dispatcher and both concrete parsers."""

    def test_parse_iers_format(self, tmp_path):
        """IERS Bulletin C format dispatch + parser body (585, 596-657)."""
        content = "2017 January 1 TAI-UTC= 37s\n2015 July 1 TAI-UTC= 36s\n"
        fp = tmp_path / "leap.dat"
        fp.write_text(content)
        entries = iers_data._parse_leap_seconds(str(fp))
        assert len(entries) == 2
        # Sorted by MJD ascending: 2015 first.
        assert entries[0].year == 2015
        assert entries[1].year == 2017
        assert entries[1].tai_utc == pytest.approx(37.0)

    def test_parse_real_leap_second_dat(self, tmp_path):
        """The real IERS Leap_Second.dat (numeric, '#'-commented) must parse.

        Regression: this format was misrouted to the IETF parser and silently
        discarded, so downloaded IERS updates never took effect.
        """
        content = (
            "#  MJD        Date        TAI-UTC (s)\n"
            "#           day month year\n"
            "   41317.0    1  1 1972      10\n"
            "   57754.0    1  1 2017      37\n"
        )
        fp = tmp_path / "Leap_Second.dat"
        fp.write_text(content)
        entries = iers_data._parse_leap_seconds(str(fp))
        assert len(entries) == 2
        assert entries[0].mjd == pytest.approx(41317.0)
        assert entries[-1].year == 2017
        assert entries[-1].tai_utc == pytest.approx(37.0)

    def test_parse_iers_unknown_month_skipped(self, tmp_path):
        """An unrecognised month name is skipped (month==0 -> continue)."""
        content = "2017 Smarch 1 TAI-UTC= 37s\n2015 July 1 TAI-UTC= 36s\n"
        fp = tmp_path / "leap_badmonth.dat"
        fp.write_text(content)
        entries = iers_data._parse_leap_seconds(str(fp))
        assert len(entries) == 1
        assert entries[0].year == 2015

    def test_parse_iers_short_and_malformed(self, tmp_path):
        """Short rows (<5 parts) and value errors are skipped (601, 654-655).

        The content must not start with '#' and must contain no 'NTP' token,
        otherwise the dispatcher routes it to the IETF parser.
        """
        content = (
            "header line that is not a leap second\n"
            "// comment-style line skipped by startswith check is absent;\n"
            "\n"  # blank line -> continue (line 601)
            "2017 January 1\n"  # too few parts (<5) -> skipped
            "abcd January 1 TAI-UTC= 37s\n"  # int() fails on year -> 654-655
            "2017 January 1 TAI-UTC= 37s\n"  # valid
        )
        fp = tmp_path / "leap_mixed.dat"
        fp.write_text(content)
        entries = iers_data._parse_leap_seconds(str(fp))
        assert len(entries) == 1
        assert entries[0].year == 2017

    def test_parse_iers_empty_falls_back_to_hardcoded(self, tmp_path):
        """Empty IERS content -> hardcoded fallback (588->589)."""
        # Content that is detected as IERS format (no NTP, not starting '#')
        # but yields zero parsed entries -> hardcoded list returned.
        fp = tmp_path / "leap_empty.dat"
        fp.write_text("not-a-leap-second-line\n")
        entries = iers_data._parse_leap_seconds(str(fp))
        assert entries == iers_data._get_hardcoded_leap_seconds()
        assert len(entries) > 0

    def test_parse_ietf_format(self, tmp_path):
        """IETF/NTP format dispatch + parser body (677-685)."""
        # NTP seconds for 2017-01-01 leap with TAI-UTC=37.
        # 3692217600 == NTP seconds for 2017-01-01.
        content = (
            "#$ 3692217600\n"  # comment -> NTP token present, format=ietf
            "# comment\n"
            "\n"
            "bad line\n"  # int() fails -> skipped
            "3692217600 37\n"
        )
        fp = tmp_path / "leap.list"
        fp.write_text(content)
        entries = iers_data._parse_leap_seconds(str(fp))
        assert len(entries) == 1
        assert entries[0].tai_utc == pytest.approx(37.0)
        assert entries[0].year == 2017

    def test_parse_ietf_short_line_skipped(self, tmp_path):
        """IETF line with <2 tokens loops/continues (675->668)."""
        content = "#NTP header\n3692217600\n3692217600 37\n"
        fp = tmp_path / "leap_short.list"
        fp.write_text(content)
        entries = iers_data._parse_leap_seconds(str(fp))
        assert len(entries) == 1


# ---------------------------------------------------------------------------
# _parse_delta_t_data validation branches
# ---------------------------------------------------------------------------


class TestParseDeltaTValidation:
    """Cover the year/month/day range guards (788, 790, 792)."""

    def test_validation_rejects_out_of_range(self, tmp_path):
        content = (
            " 1800  1  1  10.0\n"  # year < 1900 -> skip (788)
            " 2000 13  1  63.0\n"  # month > 12 -> skip (790)
            " 2000  1 40  63.0\n"  # day > 31 -> skip (792)
            " 2000  1  1  63.8285\n"  # valid
        )
        fp = tmp_path / "deltat.data"
        fp.write_text(content)
        data = iers_data._parse_delta_t_data(str(fp))
        assert len(data) == 1
        assert data[0].year == 2000


# ---------------------------------------------------------------------------
# _mjd_to_calendar (date conversion both Gregorian/Julian branches)
# ---------------------------------------------------------------------------


class TestMjdToCalendar:
    """Cover both calendar branches and month/year arms (836-865)."""

    def test_modern_gregorian_date(self):
        """A modern MJD round-trips (z >= 2299161, month > 2)."""
        mjd = iers_data._calendar_to_mjd(2000, 6, 15)
        year, month, day = iers_data._mjd_to_calendar(mjd)
        assert (year, month, day) == (2000, 6, 15)

    def test_early_year_month_branch(self):
        """A January date hits month <= 2 / year = c-4715 arm."""
        mjd = iers_data._calendar_to_mjd(2000, 1, 10)
        year, month, day = iers_data._mjd_to_calendar(mjd)
        assert (year, month, day) == (2000, 1, 10)

    def test_pre_gregorian_julian_branch(self):
        """An ancient MJD takes the z < 2299161 (Julian-calendar) branch.

        The two-step conversion (Gregorian _calendar_to_mjd then the
        Julian-aware _mjd_to_calendar) is not a clean round-trip across the
        reform boundary; we only assert the Julian arm executes and yields a
        sane year, which is what drives coverage of lines 842-844.
        """
        # A large negative MJD pushes z below 2299161 (pre-1582 Gregorian
        # reform), selecting the Julian-calendar arm (a = z, line 843).
        mjd = -200000.0
        year, month, day = iers_data._mjd_to_calendar(mjd)
        assert isinstance(year, int)
        assert 1 <= month <= 12
        assert 1 <= day <= 31


# ---------------------------------------------------------------------------
# load_iers_data auto-download branches
# ---------------------------------------------------------------------------


class TestLoadIersDataAutoDownload:
    """Cover the auto-download dispatch and the parse-failure except arms."""

    def test_auto_download_all_three(self, tmp_path):
        """Auto-download enabled, none cached -> all three downloaders run."""
        iers_data.set_iers_auto_download(True)
        calls = {"finals": 0, "leap": 0, "delta": 0}

        def fake_finals(force=False):
            calls["finals"] += 1
            return iers_data._get_finals_cache_path()

        def fake_leap(force=False):
            calls["leap"] += 1
            return iers_data._get_leap_seconds_cache_path()

        def fake_delta(force=False):
            calls["delta"] += 1
            return iers_data._get_delta_t_cache_path()

        with (
            mock.patch.object(iers_data, "download_iers_finals", fake_finals),
            mock.patch.object(iers_data, "download_leap_seconds", fake_leap),
            mock.patch.object(iers_data, "download_delta_t_data", fake_delta),
        ):
            iers_data.load_iers_data(force_download=False)
        assert calls == {"finals": 1, "leap": 1, "delta": 1}

    def test_auto_download_files_present_skip_downloads(self):
        """Auto-download on, files cached, force=False -> per-file skip.

        Covers the False arms of ``need_finals or force_download`` etc.
        (branches 909->911 and 911->913): no downloader is invoked because
        each ``need_*`` is False and force_download is False.
        """
        iers_data.set_iers_auto_download(True)
        finals = iers_data._get_finals_cache_path()
        leap = iers_data._get_leap_seconds_cache_path()
        delta = iers_data._get_delta_t_cache_path()
        os.makedirs(os.path.dirname(finals), exist_ok=True)
        for p in (finals, leap, delta):
            with open(p, "w") as f:
                f.write("# present\n")

        with (
            mock.patch.object(iers_data, "download_iers_finals") as dl_f,
            mock.patch.object(iers_data, "download_leap_seconds") as dl_l,
            mock.patch.object(iers_data, "download_delta_t_data") as dl_d,
        ):
            iers_data.load_iers_data(force_download=False)
        dl_f.assert_not_called()
        dl_l.assert_not_called()
        dl_d.assert_not_called()

    def test_auto_download_connection_error_swallowed(self):
        """A ConnectionError during auto-download is swallowed (915-917)."""
        iers_data.set_iers_auto_download(True)
        with mock.patch.object(
            iers_data,
            "download_iers_finals",
            side_effect=ConnectionError("no net"),
        ):
            # Should not raise; falls through to load-from-cache (nothing there).
            iers_data.load_iers_data(force_download=False)

    def test_parse_failures_use_fallbacks(self):
        """Parse errors on all three files hit the except arms (924, 932, 944)."""
        finals = iers_data._get_finals_cache_path()
        leap = iers_data._get_leap_seconds_cache_path()
        delta = iers_data._get_delta_t_cache_path()
        os.makedirs(os.path.dirname(finals), exist_ok=True)
        for p in (finals, leap, delta):
            with open(p, "w") as f:
                f.write("irrelevant")

        with (
            mock.patch.object(
                iers_data, "_parse_finals_data", side_effect=ValueError("bad")
            ),
            mock.patch.object(
                iers_data, "_parse_leap_seconds", side_effect=ValueError("bad")
            ),
            mock.patch.object(
                iers_data, "_parse_delta_t_data", side_effect=ValueError("bad")
            ),
        ):
            iers_data.load_iers_data(force_download=False)

        # finals -> {} (line 925); delta -> [] (line 945);
        # leap -> hardcoded fallback (line 933).
        assert iers_data._IERS_DATA == {}
        assert iers_data._DELTA_T_DATA == []
        assert iers_data._LEAP_SECONDS == iers_data._get_hardcoded_leap_seconds()


# ---------------------------------------------------------------------------
# get_tai_utc fallback
# ---------------------------------------------------------------------------


class TestGetTaiUtc:
    """Cover the no-leap-seconds fallback path (972-974)."""

    def test_fallback_modern(self):
        """No leap seconds loaded -> modern value 37.0 (line 972-973)."""
        iers_data.clear_iers_cache()
        with mock.patch.object(iers_data, "_ensure_data_loaded", lambda: None):
            assert iers_data._LEAP_SECONDS == []
            assert iers_data.get_tai_utc(58000.0) == 37.0

    def test_fallback_old(self):
        """No leap seconds loaded, post-1972 date -> 32.0 fallback."""
        iers_data.clear_iers_cache()
        with mock.patch.object(iers_data, "_ensure_data_loaded", lambda: None):
            # A post-1972 MJD (1995) with no leap-second table reaches the
            # 32.0 fallback. Pre-1972 dates now use the piecewise TAI-UTC model
            # instead and never reach this branch.
            assert iers_data.get_tai_utc(50000.0) == 32.0

    def test_lookup_with_loaded_leap_seconds(self, monkeypatch):
        """Leap seconds present -> the lookup loop runs (lines 977-984)."""
        seconds = iers_data._get_hardcoded_leap_seconds()
        monkeypatch.setattr(iers_data, "_LEAP_SECONDS", seconds)
        monkeypatch.setattr(iers_data, "_ensure_data_loaded", lambda: None)
        # A modern MJD (after 2017-01-01, mjd 57754) -> 37 leap seconds.
        assert iers_data.get_tai_utc(58000.0) == pytest.approx(37.0)
        # A date before the first leap second (pre-1972) now uses the
        # piecewise pre-1972 TAI-UTC model, clamped to the earliest offset
        # for very old dates, rather than returning 0.0.
        assert iers_data.get_tai_utc(0.0) == pytest.approx(1.422818, abs=1e-5)


# ---------------------------------------------------------------------------
# get_ut1_utc interpolation
# ---------------------------------------------------------------------------


def _install_finals(monkeypatch, points):
    """Helper: install an in-memory finals table and stub the loader."""
    table = {
        float(mjd): iers_data.IERSDataPoint(
            mjd=float(mjd),
            year=2000,
            month=1,
            day=1,
            ut1_utc=val,
        )
        for mjd, val in points.items()
    }
    monkeypatch.setattr(iers_data, "_IERS_DATA", table)
    monkeypatch.setattr(iers_data, "_ensure_data_loaded", lambda: None)
    return table


class TestGetUt1Utc:
    """Cover the UT1-UTC interpolation and edge branches (1015-1026)."""

    def test_exact_integer_mjd(self, monkeypatch):
        """frac == 0 and floor present -> exact value."""
        _install_finals(monkeypatch, {51544: 0.355})
        assert iers_data.get_ut1_utc(51544.0) == pytest.approx(0.355)

    def test_interpolated_between_days(self, monkeypatch):
        """Both bracketing days present -> linear interpolation (1015-1018)."""
        _install_finals(monkeypatch, {51544: 0.300, 51545: 0.400})
        val = iers_data.get_ut1_utc(51544.5)
        assert val == pytest.approx(0.350)

    def test_edge_floor_only(self, monkeypatch):
        """Only floor present -> closest value (line 1021-1022)."""
        _install_finals(monkeypatch, {51544: 0.300})
        val = iers_data.get_ut1_utc(51544.5)
        assert val == pytest.approx(0.300)

    def test_edge_ceil_only(self, monkeypatch):
        """Only ceil present -> closest value (line 1023-1024)."""
        _install_finals(monkeypatch, {51545: 0.400})
        val = iers_data.get_ut1_utc(51544.5)
        assert val == pytest.approx(0.400)

    def test_no_bracket_returns_none(self, monkeypatch):
        """Neither bracketing day present -> None (line 1026)."""
        _install_finals(monkeypatch, {99999: 0.1})
        assert iers_data.get_ut1_utc(51544.5) is None

    def test_empty_table_returns_none(self, monkeypatch):
        """Empty finals table -> None."""
        _install_finals(monkeypatch, {})
        assert iers_data.get_ut1_utc(51544.0) is None


class TestEnsureDataLoaded:
    """Cover the lazy-load trigger in _ensure_data_loaded (line 955)."""

    def test_triggers_load_when_all_empty(self, monkeypatch):
        """All caches empty -> load_iers_data is invoked."""
        monkeypatch.setattr(iers_data, "_IERS_DATA", {})
        monkeypatch.setattr(iers_data, "_LEAP_SECONDS", [])
        monkeypatch.setattr(iers_data, "_DELTA_T_DATA", [])
        called = {"n": 0}

        def fake_load(*a, **k):
            called["n"] += 1
            return True

        monkeypatch.setattr(iers_data, "load_iers_data", fake_load)
        iers_data._ensure_data_loaded()
        assert called["n"] == 1


# ---------------------------------------------------------------------------
# get_observed_delta_t / data range / availability
# ---------------------------------------------------------------------------


def _install_delta_t(monkeypatch, points):
    """Helper: install an in-memory delta-t table (list) and stub loader."""
    table = [
        iers_data.DeltaTDataPoint(
            mjd=iers_data._calendar_to_mjd(y, m, d),
            year=y,
            month=m,
            day=d,
            delta_t=dt,
        )
        for (y, m, d, dt) in points
    ]
    table.sort(key=lambda p: p.mjd)
    monkeypatch.setattr(iers_data, "_DELTA_T_DATA", table)
    monkeypatch.setattr(iers_data, "_ensure_data_loaded", lambda: None)
    return table


class TestObservedDeltaT:
    """Cover the observed delta-t lookup branches."""

    def test_no_data_returns_none(self, monkeypatch):
        """Empty table -> None (line 1053)."""
        _install_delta_t(monkeypatch, [])
        assert iers_data.get_observed_delta_t(2451545.0) is None

    def test_equal_bracket_returns_p1(self, monkeypatch):
        """When the bracketing points coincide -> p1.delta_t (line 1077).

        The binary search must converge so that p1.mjd == p2.mjd. With the
        duplicate MJD at the tail of the list and a query at that MJD, the
        loop leaves left/right pointing at the two equal-MJD points.
        """
        table = [
            iers_data.DeltaTDataPoint(
                mjd=51500.0, year=2000, month=1, day=1, delta_t=62.0
            ),
            iers_data.DeltaTDataPoint(
                mjd=51544.0, year=2000, month=2, day=14, delta_t=63.0
            ),
            iers_data.DeltaTDataPoint(
                mjd=51544.0, year=2000, month=2, day=14, delta_t=63.0
            ),
        ]
        monkeypatch.setattr(iers_data, "_DELTA_T_DATA", table)
        monkeypatch.setattr(iers_data, "_ensure_data_loaded", lambda: None)
        jd = iers_data._mjd_to_jd(51544.0)
        assert iers_data.get_observed_delta_t(jd) == pytest.approx(63.0)

    def test_interpolation(self, monkeypatch):
        """Normal interpolation between two distinct points."""
        _install_delta_t(monkeypatch, [(2000, 1, 1, 63.0), (2000, 2, 1, 64.0)])
        mjd_mid = (
            iers_data._calendar_to_mjd(2000, 1, 1)
            + iers_data._calendar_to_mjd(2000, 2, 1)
        ) / 2
        jd = iers_data._mjd_to_jd(mjd_mid)
        val = iers_data.get_observed_delta_t(jd)
        assert 63.0 < val < 64.0

    def test_out_of_range(self, monkeypatch):
        """A date outside the table range -> None."""
        _install_delta_t(monkeypatch, [(2000, 1, 1, 63.0), (2000, 2, 1, 64.0)])
        jd = iers_data._mjd_to_jd(iers_data._calendar_to_mjd(1990, 1, 1))
        assert iers_data.get_observed_delta_t(jd) is None

    def test_data_range_none_when_empty(self, monkeypatch):
        """get_observed_delta_t_data_range -> None when empty (line 1098)."""
        _install_delta_t(monkeypatch, [])
        assert iers_data.get_observed_delta_t_data_range() is None

    def test_data_range_present(self, monkeypatch):
        """Range returned when data present."""
        _install_delta_t(monkeypatch, [(2000, 1, 1, 63.0), (2000, 2, 1, 64.0)])
        rng = iers_data.get_observed_delta_t_data_range()
        assert rng is not None and rng[0] < rng[1]

    def test_availability_false_when_empty(self, monkeypatch):
        """is_observed_delta_t_available -> False when empty (line 1119)."""
        _install_delta_t(monkeypatch, [])
        assert iers_data.is_observed_delta_t_available(2451545.0) is False

    def test_availability_true_in_range(self, monkeypatch):
        """True when the date is inside the covered range."""
        _install_delta_t(monkeypatch, [(2000, 1, 1, 63.0), (2000, 2, 1, 64.0)])
        jd = iers_data._mjd_to_jd(iers_data._calendar_to_mjd(2000, 1, 15))
        assert iers_data.is_observed_delta_t_available(jd) is True


# ---------------------------------------------------------------------------
# get_delta_t_iers fallback to UT1-UTC
# ---------------------------------------------------------------------------


class TestGetDeltaTIers:
    """Cover the observed-vs-computed dispatch (1160-1168)."""

    def test_uses_observed_when_available(self, monkeypatch):
        """Observed value short-circuits the fallback."""
        monkeypatch.setattr(iers_data, "_ensure_data_loaded", lambda: None)
        monkeypatch.setattr(iers_data, "get_observed_delta_t", lambda jd: 63.8)
        assert iers_data.get_delta_t_iers(2451545.0) == pytest.approx(63.8)

    def test_computes_from_ut1_utc(self, monkeypatch):
        """Observed missing -> compute from UT1-UTC + TAI-UTC (1160-1166)."""
        monkeypatch.setattr(iers_data, "_ensure_data_loaded", lambda: None)
        monkeypatch.setattr(iers_data, "get_observed_delta_t", lambda jd: None)
        monkeypatch.setattr(iers_data, "get_ut1_utc", lambda mjd: 0.355)
        monkeypatch.setattr(iers_data, "get_tai_utc", lambda mjd: 37.0)
        # delta_t = 37.0 + 32.184 - 0.355
        result = iers_data.get_delta_t_iers(2451545.0)
        assert result == pytest.approx(37.0 + 32.184 - 0.355)

    def test_returns_none_when_no_ut1(self, monkeypatch):
        """Observed missing and UT1-UTC missing -> None (line 1157)."""
        monkeypatch.setattr(iers_data, "_ensure_data_loaded", lambda: None)
        monkeypatch.setattr(iers_data, "get_observed_delta_t", lambda jd: None)
        monkeypatch.setattr(iers_data, "get_ut1_utc", lambda mjd: None)
        assert iers_data.get_delta_t_iers(2451545.0) is None


# ---------------------------------------------------------------------------
# get_iers_data_range / is_iers_data_available
# ---------------------------------------------------------------------------


class TestIersDataRange:
    """Cover the finals-table range and availability helpers."""

    def test_range_none_when_empty(self, monkeypatch):
        """Empty finals table -> None."""
        _install_finals(monkeypatch, {})
        assert iers_data.get_iers_data_range() is None

    def test_range_present(self, monkeypatch):
        """Min/max MJD mapped to JD (lines 1183-1187)."""
        _install_finals(monkeypatch, {51544: 0.3, 51600: 0.4})
        rng = iers_data.get_iers_data_range()
        assert rng is not None and rng[0] < rng[1]

    def test_available_true(self, monkeypatch):
        """Rounded MJD present in finals table -> True (lines 1202-1205)."""
        _install_finals(monkeypatch, {51544: 0.3})
        jd = iers_data._mjd_to_jd(51544.0)
        assert iers_data.is_iers_data_available(jd) is True

    def test_available_false(self, monkeypatch):
        """Rounded MJD absent -> False."""
        _install_finals(monkeypatch, {51544: 0.3})
        jd = iers_data._mjd_to_jd(99999.0)
        assert iers_data.is_iers_data_available(jd) is False


# ---------------------------------------------------------------------------
# delete_iers_cache_files
# ---------------------------------------------------------------------------


class TestDeleteCacheFiles:
    """Cover the file deletion loop incl. the OSError arm (1229-1244)."""

    def test_deletes_existing_files(self):
        """Existing cache files are removed and counted."""
        cache_dir = iers_data._get_cache_dir()
        paths = []
        for name in ("finals2000A.data", "leap_seconds.dat", "deltat.data"):
            p = os.path.join(cache_dir, name)
            with open(p, "w") as f:
                f.write("x")
            paths.append(p)
        deleted = iers_data.delete_iers_cache_files()
        assert deleted == 3
        assert all(not os.path.exists(p) for p in paths)

    def test_delete_skips_absent_files(self):
        """Files that do not exist are skipped (branch 1234->1232)."""
        # Cache dir exists but contains no cache files -> nothing deleted.
        deleted = iers_data.delete_iers_cache_files()
        assert deleted == 0

    def test_remove_oserror_swallowed(self):
        """os.remove raising OSError is swallowed (lines 1238-1239)."""
        cache_dir = iers_data._get_cache_dir()
        for name in ("finals2000A.data", "leap_seconds.dat", "deltat.data"):
            with open(os.path.join(cache_dir, name), "w") as f:
                f.write("x")
        with mock.patch("os.remove", side_effect=OSError("locked")):
            deleted = iers_data.delete_iers_cache_files()
        assert deleted == 0


# ---------------------------------------------------------------------------
# get_iers_cache_info
# ---------------------------------------------------------------------------


class TestCacheInfo:
    """Cover all age/range branches of get_iers_cache_info."""

    def test_all_present(self, monkeypatch):
        """All three files exist + both ranges present (1291-1310)."""
        cache_dir = iers_data._get_cache_dir()
        for name in ("finals2000A.data", "leap_seconds.dat", "deltat.data"):
            with open(os.path.join(cache_dir, name), "w") as f:
                f.write("x")
        monkeypatch.setattr(
            iers_data, "get_iers_data_range", lambda: (2451545.0, 2451600.0)
        )
        monkeypatch.setattr(
            iers_data,
            "get_observed_delta_t_data_range",
            lambda: (2451545.0, 2451600.0),
        )
        info = iers_data.get_iers_cache_info()
        assert info["finals_exists"] is True
        assert info["finals_age_days"] is not None
        assert info["leap_seconds_age_days"] is not None
        assert info["delta_t_age_days"] is not None
        assert info["data_range_start_jd"] == 2451545.0
        assert info["delta_t_range_end_jd"] == 2451600.0

    def test_no_files_no_ranges(self, monkeypatch):
        """No cache files and no ranges -> skip-branches (1298->1302, 1308->1312)."""
        # Cache dir exists (created by _get_cache_dir) but contains no files.
        monkeypatch.setattr(iers_data, "get_iers_data_range", lambda: None)
        monkeypatch.setattr(iers_data, "get_observed_delta_t_data_range", lambda: None)
        info = iers_data.get_iers_cache_info()
        assert info["finals_exists"] is False
        assert info["finals_age_days"] is None
        assert info["delta_t_age_days"] is None
        assert info["data_range_start_jd"] is None
        assert info["delta_t_range_start_jd"] is None
