"""Owned local downloads verify before atomic publication; no remote data requests."""

from data_paths import CORE_LEB
from contextlib import contextmanager
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import inspect
import json
from pathlib import Path
import subprocess
from threading import Thread

import libephemeris as eph
import libephemeris.download as download
from libephemeris.net import NetworkSealedError
import pytest

FIXTURE = CORE_LEB


@pytest.fixture
def isolated_data(tmp_path):
    settings = tmp_path / "settings.env"
    data = tmp_path / "owned-data"
    settings.write_text(
        f"LIBEPHEMERIS_DATA_DIR={data}\nLIBEPHEMERIS_SPK_DIR={data / 'spk'}\n"
    )
    eph.load_dotenv(settings, override=True)
    eph.set_calc_mode("leb")
    eph.set_network_policy("sealed")
    try:
        yield data
    finally:
        eph.set_leb_file(str(FIXTURE))
        eph.set_network_policy(None)
        eph.reload_configuration()


@contextmanager
def local_payload(payload):
    requests = []

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            requests.append((self.path, self.headers.get("User-Agent")))
            self.send_response(200)
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)

        def log_message(self, *args):
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{server.server_port}/owned-payload", requests
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)
        assert not thread.is_alive()


def test_status_paths_are_local_metadata_without_installing(isolated_data, capfd):
    data = isolated_data
    assert download.get_data_dir() == data
    assert not data.exists()
    status = download.check_data_status()
    assert set(status) == set(download.DATA_FILES)
    assert all(
        not row["exists"] and row["path"] is None and row["size"] is None
        for row in status.values()
    )
    assert not data.exists()
    leb = download.get_leb_dir()
    assert leb == data / "leb" and leb.is_dir()
    assert (
        download.get_leb_path_for_tier("future-label")
        == leb / "ephemeris_future-label.leb"
    )
    sample = data / "planet_centers.bsp"
    sample.write_bytes(b"metadata does not claim structural validity")
    row = download.check_data_status()[sample.name]
    assert (
        row["exists"] and row["path"] == sample and row["size"] == sample.stat().st_size
    )
    assert download.print_data_status(as_json=True) is None
    snapshot = json.loads(capfd.readouterr().out)
    assert snapshot["data_directory"] == str(data)
    assert snapshot["calc_mode"] == "leb"
    assert snapshot["spk_cache"]["file_count"] == 0
    assert not (data / "spk").exists()


def test_publish_replaces_only_the_requested_destination(tmp_path):
    temporary = tmp_path / "complete.tmp"
    destination = tmp_path / "installed.bin"
    sibling = tmp_path / "unrelated.tmp"
    temporary.write_bytes(b"new")
    destination.write_bytes(b"old")
    sibling.write_bytes(b"owned elsewhere")
    assert download.publish_temp_file(temporary, destination) is None
    assert destination.read_bytes() == b"new" and not temporary.exists()
    assert sibling.read_bytes() == b"owned elsewhere"
    with pytest.raises(OSError):
        download.publish_temp_file(temporary, destination)
    assert destination.read_bytes() == b"new"


def test_stream_hash_and_format_refusals_preserve_destination(isolated_data, tmp_path):
    payload = b"actual owned local bytes" * 8000
    digest = hashlib.sha256(payload).hexdigest()
    destination = tmp_path / "atomic" / "result.bin"
    eph.set_network_policy("allow")
    with local_payload(payload) as (url, requests):
        assert download.download_file(url, destination, expected_sha256=digest)
        assert destination.read_bytes() == payload
        assert requests == [("/owned-payload", "libephemeris-download/1.0")]
        with pytest.raises(ValueError):
            download.download_file(url, destination, expected_sha256="0" * 64)
        assert destination.read_bytes() == payload
        with pytest.raises(ValueError):
            download.download_file(
                url, destination, expected_sha256=digest, validator="bsp"
            )
        assert destination.read_bytes() == payload
        with pytest.raises(ValueError):
            download.download_file(url, destination, validator="unknown")
        with pytest.raises(TypeError):
            download.download_file(url, destination, validator=lambda path: True)
    assert list(destination.parent.iterdir()) == [destination]


def test_sealed_and_bundled_group_installation_are_offline(isolated_data, tmp_path):
    destination = tmp_path / "sealed" / "result.bin"
    destination.parent.mkdir()
    destination.write_bytes(b"keep")
    with pytest.raises(NetworkSealedError):
        download.download_file("not a URL", destination)
    assert destination.read_bytes() == b"keep"
    assert list(destination.parent.iterdir()) == [destination]
    installed = download.download_leb2_for_tier(
        "base", ["unknown", "core"], quiet=True, activate=False
    )
    assert installed == [isolated_data / "leb/base_core.leb2"]
    assert (
        hashlib.sha256(installed[0].read_bytes()).hexdigest()
        == download.DATA_FILES["base_core.leb2"]["sha256"]
    )
    assert (
        download.download_leb2_for_tier("base", ["core"], quiet=True, activate=False)
        == installed
    )
    assert download.download_leb2_for_tier("base", [], quiet=True) == []
    with pytest.raises(ValueError):
        download.download_leb2_for_tier("unknown", quiet=True)
    summary = download.init_all(
        cache_dir=str(tmp_path / "empty-cache"),
        quiet=True,
        start_year=2000,
        end_year=2000,
    )
    assert (
        summary["spk_success"]
        == summary["spk_failed"]
        == summary["spk_skipped"]
        == summary["spk_out_of_range"]
        == 0
    )
    assert summary["failed_details"] == []
