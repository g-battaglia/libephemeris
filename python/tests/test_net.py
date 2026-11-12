"""Explicit network settings and the provider-dependent auto gate."""
from data_paths import CORE_LEB

import libephemeris as eph
import pytest
from libephemeris.net import NetworkSealedError, require_network


def test_network_policy_resolution_and_failed_setter_is_atomic():
    try:
        for mode, expected in (("leb", "sealed"), ("auto", "allow"), ("skyfield", "allow"), ("horizons", "allow")):
            eph.set_calc_mode(mode)
            eph.set_network_policy("auto")
            assert eph.get_configured_network_policy() == "auto"
            assert eph.get_network_policy() == expected
        eph.set_network_policy(" SEALED ")
        assert eph.get_configured_network_policy() == eph.get_network_policy() == "sealed"
        with pytest.raises(ValueError):
            eph.set_network_policy("sometimes")
        assert eph.get_network_policy() == "sealed"
    finally:
        eph.set_calc_mode("leb")
        eph.set_network_policy(None)


def test_public_network_gate_refuses_without_opening_data_or_sockets(monkeypatch, tmp_path):
    import socket
    from pathlib import Path

    def unexpected(*args, **kwargs):
        raise AssertionError("a sealed gate must not access a socket or DNS")

    monkeypatch.setattr(socket, "socket", unexpected)
    monkeypatch.setattr(socket, "getaddrinfo", unexpected)
    fixture = str(CORE_LEB)
    try:
        eph.set_calc_mode("leb")
        eph.set_leb_file(str(tmp_path / "not-opened.leb2"))
        eph.set_network_policy("sealed")
        with pytest.raises(NetworkSealedError) as error:
            require_network("sealed regression probe")
        assert error.value.purpose == "sealed regression probe"
        assert issubclass(NetworkSealedError, RuntimeError)
        eph.set_network_policy("allow")
        assert require_network("explicitly permitted operation") is None
        eph.set_network_policy("auto")
        with pytest.raises(NetworkSealedError):
            require_network("coefficient-only auto policy")
        eph.set_calc_mode("skyfield")
        assert require_network("native JPL auto policy") is None
    finally:
        eph.set_calc_mode("leb")
        eph.set_leb_file(fixture)
        eph.set_network_policy(None)
