"""Regression tests for the LEB sealed-network invariant."""

from __future__ import annotations

import ast
import socket
from pathlib import Path

import pytest

import libephemeris as eph
from libephemeris import state
from libephemeris.exceptions import NetworkSealedError


@pytest.fixture(autouse=True)
def _restore_process_policy():
    eph.close()
    eph.set_calc_mode("leb")
    eph.set_network_policy("sealed")
    yield
    eph.close()
    eph.set_calc_mode(None)
    eph.set_network_policy(None)


def test_runtime_network_imports_and_urlopen_have_one_choke_point() -> None:
    """No package module except net.py may import networking or call urlopen."""
    package_dir = Path(eph.__file__).parent
    offenders: list[str] = []
    restricted_modules = ("urllib", "http.client", "socket", "requests", "httpx")
    for path in package_dir.rglob("*.py"):
        if path.name == "net.py":
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name in restricted_modules or alias.name.startswith(
                        tuple(f"{name}." for name in restricted_modules)
                    ):
                        offenders.append(
                            f"{path.relative_to(package_dir)}:{alias.name}"
                        )
            elif isinstance(node, ast.ImportFrom):
                module = node.module or ""
                if module in restricted_modules or module.startswith(
                    tuple(f"{name}." for name in restricted_modules)
                ):
                    offenders.append(f"{path.relative_to(package_dir)}:{module}")
                if module == "skyfield.api" and any(
                    alias.name == "Loader" for alias in node.names
                ):
                    offenders.append(
                        f"{path.relative_to(package_dir)}:raw Skyfield Loader"
                    )
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            if isinstance(func, ast.Attribute) and func.attr == "urlopen":
                offenders.append(str(path.relative_to(package_dir)))
            elif isinstance(func, ast.Name) and func.id == "urlopen":
                offenders.append(str(path.relative_to(package_dir)))
    assert offenders == []


def test_sealed_open_url_stops_before_socket(monkeypatch: pytest.MonkeyPatch) -> None:
    attempts: list[tuple[str, object]] = []

    def _blocked_connect(*args, **kwargs):  # type: ignore[no-untyped-def]
        attempts.append(("connect", args))
        raise AssertionError("socket connect was attempted")

    def _blocked_dns(*args, **kwargs):  # type: ignore[no-untyped-def]
        attempts.append(("dns", args))
        raise AssertionError("DNS resolution was attempted")

    monkeypatch.setattr(socket.socket, "connect", _blocked_connect)
    monkeypatch.setattr(socket, "getaddrinfo", _blocked_dns)

    from libephemeris.net import open_url

    with pytest.raises(NetworkSealedError, match="sealed-network smoke test"):
        open_url(
            "https://example.invalid/",
            purpose="sealed-network smoke test",
            timeout=0.1,
        )
    assert attempts == []


def test_leb_loader_cannot_implicitly_download_kernel(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Sealed LEB mode rejects JPL before consulting Skyfield's Loader."""
    attempts: list[str] = []

    def _blocked_connect(*args, **kwargs):  # type: ignore[no-untyped-def]
        attempts.append("connect")
        raise AssertionError("socket connect was attempted")

    def _blocked_dns(*args, **kwargs):  # type: ignore[no-untyped-def]
        attempts.append("dns")
        raise AssertionError("DNS resolution was attempted")

    monkeypatch.setattr(socket.socket, "connect", _blocked_connect)
    monkeypatch.setattr(socket, "getaddrinfo", _blocked_dns)
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))

    with pytest.raises(RuntimeError, match="JPL/SPICE ephemeris access is disabled"):
        state.get_planets()
    assert attempts == []


def test_context_loader_cannot_bypass_sealed_policy(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The context API cannot bypass the sealed LEB source boundary."""
    attempts: list[str] = []

    def _blocked_connect(*_args, **_kwargs):  # type: ignore[no-untyped-def]
        attempts.append("connect")
        raise AssertionError("socket connect was attempted")

    def _blocked_dns(*_args, **_kwargs):  # type: ignore[no-untyped-def]
        attempts.append("dns")
        raise AssertionError("DNS resolution was attempted")

    monkeypatch.setattr(socket.socket, "connect", _blocked_connect)
    monkeypatch.setattr(socket, "getaddrinfo", _blocked_dns)
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))

    eph.EphemerisContext.close()
    with pytest.raises(RuntimeError, match="JPL/SPICE ephemeris access is disabled"):
        eph.EphemerisContext(ephe_file="missing-context-kernel.bsp").get_planets()
    assert attempts == []


def test_explicit_download_cli_boundary_allows_network_policy(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Provisioning can explicitly override a sealed calculation config."""
    from click.testing import CliRunner

    from libephemeris.cli import cli
    from libephemeris import download

    monkeypatch.setattr(download, "download_leb2_for_tier", lambda **_kwargs: [])

    result = CliRunner().invoke(cli, ["download", "leb2-base", "--quiet"])
    assert result.exit_code == 0
    assert eph.get_network_policy() == "allow"


def test_explicit_iers_cli_uses_provisioning_boundary(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from click.testing import CliRunner

    from libephemeris import iers_data
    from libephemeris.cli import cli

    calls: list[bool] = []

    def _download(*, force: bool) -> str:
        calls.append(force)
        return "/tmp/iers"

    monkeypatch.setattr(iers_data, "download_iers_finals", _download)
    monkeypatch.setattr(iers_data, "download_leap_seconds", _download)
    monkeypatch.setattr(iers_data, "download_delta_t_data", _download)

    result = CliRunner().invoke(cli, ["download", "iers", "--force", "--quiet"])

    assert result.exit_code == 0
    assert calls == [True, True, True]
    assert eph.get_network_policy() == "allow"


def test_invalid_env_network_policy_fails_loudly(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A typo in the env policy must not silently widen the boundary."""
    from libephemeris import net

    eph.set_network_policy(None)
    monkeypatch.setenv(net.NETWORK_POLICY_ENV, "seales")

    with pytest.raises(ValueError, match="'seales'") as excinfo:
        eph.get_configured_network_policy()
    assert "['auto', 'allow', 'sealed']" in str(excinfo.value)
    assert net.NETWORK_POLICY_ENV in str(excinfo.value)


def test_invalid_env_policy_propagates_through_effective_resolution(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """No caller between config and enforcement may swallow the error."""
    from libephemeris import net

    eph.set_network_policy(None)
    monkeypatch.setenv(net.NETWORK_POLICY_ENV, "seales")

    with pytest.raises(ValueError, match="'seales'"):
        eph.get_network_policy()


def test_invalid_toml_network_policy_fails_loudly(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from libephemeris import _config_toml, net

    eph.set_network_policy(None)
    monkeypatch.delenv(net.NETWORK_POLICY_ENV, raising=False)
    monkeypatch.setattr(
        _config_toml,
        "get_str",
        lambda key: "sealde" if key == "network_policy" else None,
    )

    with pytest.raises(ValueError, match="'sealde'") as excinfo:
        eph.get_configured_network_policy()
    assert "network_policy" in str(excinfo.value)


def test_absent_and_blank_configured_policy_resolve_to_auto(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from libephemeris import _config_toml, net

    eph.set_network_policy(None)
    monkeypatch.delenv(net.NETWORK_POLICY_ENV, raising=False)
    monkeypatch.setattr(_config_toml, "get_str", lambda _key: None)
    assert eph.get_configured_network_policy() == "auto"

    monkeypatch.setenv(net.NETWORK_POLICY_ENV, "   ")
    monkeypatch.setattr(_config_toml, "get_str", lambda _key: "")
    assert eph.get_configured_network_policy() == "auto"


def test_configured_policy_normalizes_case_and_whitespace(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from libephemeris import _config_toml, net

    eph.set_network_policy(None)
    monkeypatch.setenv(net.NETWORK_POLICY_ENV, "  SEALED ")
    assert eph.get_configured_network_policy() == "sealed"

    monkeypatch.delenv(net.NETWORK_POLICY_ENV, raising=False)
    monkeypatch.setattr(
        _config_toml,
        "get_str",
        lambda key: " Allow" if key == "network_policy" else None,
    )
    assert eph.get_configured_network_policy() == "allow"


def test_extended_auto_spk_padding_is_clamped_to_horizons_limits(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from libephemeris import planets, spk

    captured: dict[str, object] = {}

    def _capture_download(**kwargs):  # type: ignore[no-untyped-def]
        captured.update(kwargs)

    monkeypatch.setattr(
        state,
        "get_spk_date_range_for_tier",
        lambda: ("1600-01-01", "2500-01-01"),
    )
    monkeypatch.setattr(spk, "download_and_register_spk", _capture_download)
    monkeypatch.setattr(spk, "calc_spk_body_position", lambda *_args: (1.0,) * 6)
    eph.set_network_policy("allow")

    result = planets._try_auto_spk_download(object(), eph.CHIRON, 0)

    assert result == (1.0,) * 6
    assert captured["start"] == "1600-01-01"
    assert captured["end"] == "2500-01-01"
