"""Tests for the built-in .env file loader (libephemeris/_dotenv.py)."""

from __future__ import annotations

import os
import textwrap
from pathlib import Path
from unittest import mock


from libephemeris._dotenv import _parse_line, _find_env_file, load_dotenv


# =============================================================================
# _parse_line tests
# =============================================================================


class TestParseLine:
    """Unit tests for single-line parsing."""

    def test_simple_key_value(self) -> None:
        assert _parse_line("FOO=bar") == ("FOO", "bar")

    def test_key_value_with_spaces_around_eq(self) -> None:
        assert _parse_line("FOO = bar") == ("FOO", "bar")

    def test_double_quoted_value(self) -> None:
        assert _parse_line('FOO="bar baz"') == ("FOO", "bar baz")

    def test_single_quoted_value(self) -> None:
        assert _parse_line("FOO='bar baz'") == ("FOO", "bar baz")

    def test_empty_value(self) -> None:
        assert _parse_line("FOO=") == ("FOO", "")

    def test_empty_quoted_value(self) -> None:
        assert _parse_line('FOO=""') == ("FOO", "")

    def test_export_prefix(self) -> None:
        assert _parse_line("export FOO=bar") == ("FOO", "bar")

    def test_export_with_quotes(self) -> None:
        assert _parse_line('export FOO="hello world"') == ("FOO", "hello world")

    def test_comment_line(self) -> None:
        assert _parse_line("# this is a comment") is None

    def test_empty_line(self) -> None:
        assert _parse_line("") is None

    def test_whitespace_only_line(self) -> None:
        assert _parse_line("   ") is None

    def test_inline_comment_unquoted(self) -> None:
        assert _parse_line("FOO=bar # this is a comment") == ("FOO", "bar")

    def test_inline_comment_not_stripped_in_quotes(self) -> None:
        assert _parse_line('FOO="bar # not a comment"') == (
            "FOO",
            "bar # not a comment",
        )

    def test_no_equals_sign(self) -> None:
        assert _parse_line("JUSTAKEYWITHNOEQ") is None

    def test_empty_key(self) -> None:
        assert _parse_line("=value") is None

    def test_value_with_equals(self) -> None:
        assert _parse_line("FOO=bar=baz") == ("FOO", "bar=baz")

    def test_value_with_path(self) -> None:
        assert _parse_line("DATA_DIR=/Volumes/Data/libephemeris") == (
            "DATA_DIR",
            "/Volumes/Data/libephemeris",
        )

    def test_leading_trailing_whitespace(self) -> None:
        assert _parse_line("  FOO=bar  ") == ("FOO", "bar")

    def test_quoted_value_with_equals(self) -> None:
        assert _parse_line('FOO="a=b=c"') == ("FOO", "a=b=c")

    def test_hash_in_value_without_space(self) -> None:
        # "FOO=bar#baz" — the # is NOT preceded by a space, so no inline comment
        assert _parse_line("FOO=bar#baz") == ("FOO", "bar#baz")


# =============================================================================
# load_dotenv tests
# =============================================================================


class TestLoadDotenv:
    """Integration tests for load_dotenv."""

    def test_load_from_explicit_path(self, tmp_path: Path) -> None:
        env_file = tmp_path / ".env"
        env_file.write_text(
            "LIBEPHEMERIS_TEST_DOTENV_A=hello\nLIBEPHEMERIS_TEST_DOTENV_B=world\n"
        )

        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("LIBEPHEMERIS_TEST_DOTENV_A", None)
            os.environ.pop("LIBEPHEMERIS_TEST_DOTENV_B", None)

            result = load_dotenv(env_file)

            assert result is True
            assert os.environ["LIBEPHEMERIS_TEST_DOTENV_A"] == "hello"
            assert os.environ["LIBEPHEMERIS_TEST_DOTENV_B"] == "world"

        # Cleanup
        os.environ.pop("LIBEPHEMERIS_TEST_DOTENV_A", None)
        os.environ.pop("LIBEPHEMERIS_TEST_DOTENV_B", None)

    def test_does_not_override_existing(self, tmp_path: Path) -> None:
        env_file = tmp_path / ".env"
        env_file.write_text("LIBEPHEMERIS_TEST_DOTENV_X=from_file\n")

        with mock.patch.dict(
            os.environ, {"LIBEPHEMERIS_TEST_DOTENV_X": "from_env"}, clear=False
        ):
            load_dotenv(env_file)
            assert os.environ["LIBEPHEMERIS_TEST_DOTENV_X"] == "from_env"

    def test_override_existing_when_flag_set(self, tmp_path: Path) -> None:
        env_file = tmp_path / ".env"
        env_file.write_text("LIBEPHEMERIS_TEST_DOTENV_Y=from_file\n")

        with mock.patch.dict(
            os.environ, {"LIBEPHEMERIS_TEST_DOTENV_Y": "from_env"}, clear=False
        ):
            load_dotenv(env_file, override=True)
            assert os.environ["LIBEPHEMERIS_TEST_DOTENV_Y"] == "from_file"

    def test_returns_false_when_file_missing(self, tmp_path: Path) -> None:
        result = load_dotenv(tmp_path / "nonexistent.env")
        assert result is False

    def test_returns_false_when_no_file_found(self) -> None:
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("LIBEPHEMERIS_ENV_FILE", None)
            # Point search to non-existent dirs
            with mock.patch(
                "libephemeris._dotenv._SEARCH_PATHS",
                (lambda: Path("/nonexistent/path/.env"),),
            ):
                result = load_dotenv()
                assert result is False

    def test_skips_comments_and_blank_lines(self, tmp_path: Path) -> None:
        env_file = tmp_path / ".env"
        env_file.write_text(
            textwrap.dedent("""\
                # This is a comment
                LIBEPHEMERIS_TEST_DOTENV_C=value1

                # Another comment
                LIBEPHEMERIS_TEST_DOTENV_D=value2
            """)
        )

        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("LIBEPHEMERIS_TEST_DOTENV_C", None)
            os.environ.pop("LIBEPHEMERIS_TEST_DOTENV_D", None)

            load_dotenv(env_file)

            assert os.environ["LIBEPHEMERIS_TEST_DOTENV_C"] == "value1"
            assert os.environ["LIBEPHEMERIS_TEST_DOTENV_D"] == "value2"

        os.environ.pop("LIBEPHEMERIS_TEST_DOTENV_C", None)
        os.environ.pop("LIBEPHEMERIS_TEST_DOTENV_D", None)

    def test_handles_export_prefix(self, tmp_path: Path) -> None:
        env_file = tmp_path / ".env"
        env_file.write_text("export LIBEPHEMERIS_TEST_DOTENV_E=exported_val\n")

        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("LIBEPHEMERIS_TEST_DOTENV_E", None)
            load_dotenv(env_file)
            assert os.environ["LIBEPHEMERIS_TEST_DOTENV_E"] == "exported_val"

        os.environ.pop("LIBEPHEMERIS_TEST_DOTENV_E", None)

    def test_handles_quoted_values_with_spaces(self, tmp_path: Path) -> None:
        env_file = tmp_path / ".env"
        env_file.write_text(
            "LIBEPHEMERIS_TEST_DOTENV_F=\"hello world\"\nLIBEPHEMERIS_TEST_DOTENV_G='single quoted'\n"
        )

        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("LIBEPHEMERIS_TEST_DOTENV_F", None)
            os.environ.pop("LIBEPHEMERIS_TEST_DOTENV_G", None)

            load_dotenv(env_file)

            assert os.environ["LIBEPHEMERIS_TEST_DOTENV_F"] == "hello world"
            assert os.environ["LIBEPHEMERIS_TEST_DOTENV_G"] == "single quoted"

        os.environ.pop("LIBEPHEMERIS_TEST_DOTENV_F", None)
        os.environ.pop("LIBEPHEMERIS_TEST_DOTENV_G", None)

    def test_real_world_libephemeris_config(self, tmp_path: Path) -> None:
        """Test with realistic libephemeris environment variables."""
        env_file = tmp_path / ".env"
        env_file.write_text(
            textwrap.dedent("""\
                # Libephemeris configuration
                LIBEPHEMERIS_DATA_DIR=/Volumes/Data/libephemeris
                LIBEPHEMERIS_PRECISION=extended
                LIBEPHEMERIS_LOG_LEVEL=INFO
            """)
        )

        keys = [
            "LIBEPHEMERIS_DATA_DIR",
            "LIBEPHEMERIS_PRECISION",
            "LIBEPHEMERIS_LOG_LEVEL",
        ]

        with mock.patch.dict(os.environ, {}, clear=False):
            for k in keys:
                os.environ.pop(k, None)

            load_dotenv(env_file)

            assert os.environ["LIBEPHEMERIS_DATA_DIR"] == "/Volumes/Data/libephemeris"
            assert os.environ["LIBEPHEMERIS_PRECISION"] == "extended"
            assert os.environ["LIBEPHEMERIS_LOG_LEVEL"] == "INFO"

        for k in keys:
            os.environ.pop(k, None)

    def test_non_libephemeris_keys_are_not_exported(self, tmp_path):
        """Keys outside LIBEPHEMERIS_* never reach os.environ.

        A project .env regularly holds unrelated secrets; importing the
        library must not export them into the process.
        """
        import os

        env_file = tmp_path / ".env"
        env_file.write_text("SECRET_API_KEY=hunter2\nLIBEPHEMERIS_TEST_DOTENV_OK=1\n")
        os.environ.pop("SECRET_API_KEY", None)
        os.environ.pop("LIBEPHEMERIS_TEST_DOTENV_OK", None)
        try:
            assert load_dotenv(env_file) is True
            assert "SECRET_API_KEY" not in os.environ
            assert os.environ["LIBEPHEMERIS_TEST_DOTENV_OK"] == "1"
        finally:
            os.environ.pop("LIBEPHEMERIS_TEST_DOTENV_OK", None)

    def test_string_path_argument(self, tmp_path: Path) -> None:
        env_file = tmp_path / ".env"
        env_file.write_text("LIBEPHEMERIS_TEST_DOTENV_H=strpath\n")

        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("LIBEPHEMERIS_TEST_DOTENV_H", None)
            result = load_dotenv(str(env_file))
            assert result is True
            assert os.environ["LIBEPHEMERIS_TEST_DOTENV_H"] == "strpath"

        os.environ.pop("LIBEPHEMERIS_TEST_DOTENV_H", None)


# =============================================================================
# _find_env_file tests
# =============================================================================


class TestFindEnvFile:
    """Tests for .env file discovery."""

    def test_env_var_override(self, tmp_path: Path) -> None:
        env_file = tmp_path / "custom.env"
        env_file.write_text("X=1\n")

        with mock.patch.dict(
            os.environ, {"LIBEPHEMERIS_ENV_FILE": str(env_file)}, clear=False
        ):
            found = _find_env_file()
            assert found == env_file

    def test_env_var_override_missing_file(self) -> None:
        with mock.patch.dict(
            os.environ,
            {"LIBEPHEMERIS_ENV_FILE": "/nonexistent/path/.env"},
            clear=False,
        ):
            found = _find_env_file()
            assert found is None

    def test_cwd_env_file(self, tmp_path: Path) -> None:
        env_file = tmp_path / ".env"
        env_file.write_text("X=1\n")

        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("LIBEPHEMERIS_ENV_FILE", None)
            with mock.patch(
                "libephemeris._dotenv._SEARCH_PATHS",
                (lambda: env_file, lambda: Path("/nonexistent/.env")),
            ):
                found = _find_env_file()
                assert found == env_file

    def test_home_env_file(self, tmp_path: Path) -> None:
        home_env = tmp_path / ".libephemeris" / ".env"
        home_env.parent.mkdir(parents=True)
        home_env.write_text("X=1\n")

        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("LIBEPHEMERIS_ENV_FILE", None)
            with mock.patch(
                "libephemeris._dotenv._SEARCH_PATHS",
                (lambda: Path("/nonexistent/.env"), lambda: home_env),
            ):
                found = _find_env_file()
                assert found == home_env

    def test_cwd_takes_priority_over_home(self, tmp_path: Path) -> None:
        cwd_env = tmp_path / "project" / ".env"
        cwd_env.parent.mkdir(parents=True)
        cwd_env.write_text("SOURCE=cwd\n")

        home_env = tmp_path / "home" / ".libephemeris" / ".env"
        home_env.parent.mkdir(parents=True)
        home_env.write_text("SOURCE=home\n")

        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("LIBEPHEMERIS_ENV_FILE", None)
            with mock.patch(
                "libephemeris._dotenv._SEARCH_PATHS",
                (lambda: cwd_env, lambda: home_env),
            ):
                found = _find_env_file()
                assert found == cwd_env

    def test_no_env_file_found(self) -> None:
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("LIBEPHEMERIS_ENV_FILE", None)
            with mock.patch(
                "libephemeris._dotenv._SEARCH_PATHS",
                (lambda: Path("/nonexistent/.env"),),
            ):
                found = _find_env_file()
                assert found is None


# =============================================================================
# Public API tests
# =============================================================================


class TestPublicAPI:
    """Test that load_dotenv is accessible from the public package API."""

    def test_importable_from_package(self) -> None:
        from libephemeris import load_dotenv as fn

        assert callable(fn)

    def test_in_all(self) -> None:
        import libephemeris

        assert "load_dotenv" in libephemeris.__all__


class TestErrorBranches:
    """Cover the defensive except arms."""

    def test_find_env_file_skips_raising_search_path(self, monkeypatch) -> None:
        import libephemeris._dotenv as dotenv

        def _boom() -> Path:
            raise OSError("unreadable home dir")

        monkeypatch.delenv(dotenv._ENV_FILE_VAR, raising=False)
        monkeypatch.setattr(dotenv, "_SEARCH_PATHS", (_boom,))

        # The raising search path is swallowed; no candidate -> None.
        assert dotenv._find_env_file() is None

    def test_load_dotenv_returns_false_on_unreadable_file(
        self, tmp_path, monkeypatch
    ) -> None:
        # A file that exists but cannot be decoded as UTF-8 must yield False.
        bad = tmp_path / "bad.env"
        bad.write_bytes(b"\xff\xfe\x00invalid")
        assert load_dotenv(str(bad)) is False
