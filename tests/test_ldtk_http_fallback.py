from pathlib import Path

import pytest

pytest.importorskip("ldtk")

from exotic.api import gael_ld  # noqa: E402


class DummyResponse:
    def __init__(self, chunks=None, text=""):
        self._chunks = chunks
        self.text = text
        self.closed = False

    def raise_for_status(self):
        return None

    def iter_content(self, chunk_size):
        return iter(self._chunks)

    def close(self):
        self.closed = True


class DummyLDTkFile:
    def __init__(self, cache_path):
        self.name = "lte02300+0.00+0.5.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits"
        self._zstr = "Z+0.5"
        self.local_path = str(Path(cache_path) / self._zstr / self.name)

    @property
    def local_exists(self):
        return Path(self.local_path).exists()


class DummyClient:
    def __init__(self, cache_path):
        self.edir = "SpecInt50FITS/PHOENIX-ACES-AGSS-COND-SPECINT-2011"
        self.files = [DummyLDTkFile(cache_path)]
        self.not_cached = len(self.files)
        self.checked_paths = None

    def check_file_corruption(self, paths):
        self.checked_paths = paths
        return False


def test_ldtk_http_fallback_downloads_missing_file_to_ldtk_cache(monkeypatch, tmp_path):
    client = DummyClient(tmp_path / "cache_vis-lowres")
    captured = {}

    def fake_get(url, stream, timeout):
        captured["url"] = url
        captured["stream"] = stream
        captured["timeout"] = timeout
        return DummyResponse([b"phoenix", b"", b"-fits"])

    monkeypatch.setenv(gael_ld._LDTK_HTTP_FALLBACK_ENV, "https://mirror.example/PHOENIX/")
    monkeypatch.setattr(gael_ld.requests, "get", fake_get)

    assert gael_ld._download_ldtk_uncached_files_from_http(client) is False

    assert captured == {
        "url": (
            "https://mirror.example/PHOENIX/"
            "SpecInt50FITS/PHOENIX-ACES-AGSS-COND-SPECINT-2011/"
            "Z%2B0.5/lte02300%2B0.00%2B0.5.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits"
        ),
        "stream": True,
        "timeout": gael_ld._LDTK_DOWNLOAD_TIMEOUT,
    }
    assert Path(client.files[0].local_path).read_bytes() == b"phoenix-fits"
    assert client.checked_paths == [client.files[0].local_path]
    assert client.not_cached == 0


def test_ldtk_http_fallback_tries_gwdg_before_nextastro(monkeypatch, tmp_path):
    client = DummyClient(tmp_path / "cache_vis-lowres")
    requested_urls = []

    def fake_get(url, stream, timeout):
        requested_urls.append(url)
        if url.startswith("https://ftp.gwdg.de/"):
            raise RuntimeError("gwdg unavailable")
        return DummyResponse([b"nextastro"])

    monkeypatch.delenv(gael_ld._LDTK_HTTP_FALLBACK_ENV, raising=False)
    monkeypatch.setattr(gael_ld.requests, "get", fake_get)

    assert gael_ld._download_ldtk_uncached_files_from_http(client) is False

    assert requested_urls[0].startswith("https://ftp.gwdg.de/pub/misc/phoenix/")
    assert requested_urls[1].startswith("https://downloads.nextastro.org/PHOENIX/")
    assert Path(client.files[0].local_path).read_bytes() == b"nextastro"


def test_ldtk_download_wrapper_tries_original_before_http_fallback(monkeypatch, tmp_path):
    from ldtk.client import Client

    client = DummyClient(tmp_path / "cache_vis-lowres")
    calls = []

    def fake_original(self, force=False):
        calls.append(("ftp", force))
        raise RuntimeError("ftp unavailable")

    def fake_fallback(self, force=False):
        calls.append(("http", force))
        return False

    monkeypatch.setattr(gael_ld, "_LDTK_ORIGINAL_DOWNLOAD_UNCACHED_FILES", fake_original)
    monkeypatch.setattr(gael_ld, "_download_ldtk_uncached_files_from_http", fake_fallback)

    assert Client.download_uncached_files(client, force=True) is False
    assert calls == [("ftp", True), ("http", True)]


def test_ldtk_file_list_wrapper_uses_http_fallback_when_ftp_listing_fails(monkeypatch):
    from ldtk.client import Client

    class ClientStub:
        edir = "SpecInt50FITS/PHOENIX-ACES-AGSS-COND-SPECINT-2011"

    calls = []
    requested_urls = []

    def fake_original(self):
        calls.append("ftp")
        raise TimeoutError("ftp listing timed out")

    def fake_get(url, timeout):
        requested_urls.append(url)
        if url.endswith("PHOENIX-ACES-AGSS-COND-SPECINT-2011/"):
            return DummyResponse(text="""
                <a href="../">../</a>
                <a href="README.txt">README.txt</a>
                <a href="Z%2B0.5/">Z+0.5/</a>
                <a href="Z-0.0/">Z-0.0/</a>
            """)
        if url.endswith("Z%2B0.5/"):
            return DummyResponse(text="""
                <a href="../">../</a>
                <a href="lte02300%2B0.00%2B0.5.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits">
                    lte02300+0.00+0.5.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
                </a>
            """)
        if url.endswith("Z-0.0/"):
            return DummyResponse(text="""
                <a href="../">../</a>
                <a href="lte02300-0.00-0.0.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits">
                    lte02300-0.00-0.0.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits
                </a>
            """)
        raise AssertionError(f"unexpected URL: {url}")

    monkeypatch.setenv(gael_ld._LDTK_HTTP_FALLBACK_ENV, "https://mirror.example/PHOENIX")
    monkeypatch.setattr(gael_ld, "_LDTK_ORIGINAL_GET_SERVER_FILE_LIST", fake_original)
    monkeypatch.setattr(gael_ld.requests, "get", fake_get)

    files = Client.get_server_file_list(ClientStub())

    assert calls == ["ftp"]
    assert requested_urls == [
        "https://mirror.example/PHOENIX/SpecInt50FITS/PHOENIX-ACES-AGSS-COND-SPECINT-2011/",
        "https://mirror.example/PHOENIX/SpecInt50FITS/PHOENIX-ACES-AGSS-COND-SPECINT-2011/Z%2B0.5/",
        "https://mirror.example/PHOENIX/SpecInt50FITS/PHOENIX-ACES-AGSS-COND-SPECINT-2011/Z-0.0/",
    ]
    assert files == {
        "Z+0.5": ["lte02300+0.00+0.5.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits"],
        "Z-0.0": ["lte02300-0.00-0.0.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits"],
    }
