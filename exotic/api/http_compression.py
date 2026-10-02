import gzip
import json

try:
    import zstandard
except ImportError:  # pragma: no cover - gzip fallback covers environments without zstandard
    zstandard = None


_GZIP_LEVEL = 6
_ZSTD_LEVEL = 3


def build_compressed_json_request(payload, content_encoding=None):
    raw_body = json.dumps(payload, separators=(",", ":"), ensure_ascii=False).encode("utf-8")

    normalized_encoding = None if content_encoding is None else str(content_encoding).strip().lower()
    if normalized_encoding not in (None, "", "gzip", "zstd"):
        raise ValueError(f"Unsupported content encoding: {content_encoding}")

    if normalized_encoding == "zstd":
        if zstandard is None:
            raise RuntimeError("zstandard compression requested, but the zstandard package is unavailable.")
        compressed_body = zstandard.ZstdCompressor(level=_ZSTD_LEVEL).compress(raw_body)
        encoding = "zstd"
    elif normalized_encoding == "gzip":
        compressed_body = gzip.compress(raw_body, compresslevel=_GZIP_LEVEL)
        encoding = "gzip"
    elif zstandard is not None:
        compressed_body = zstandard.ZstdCompressor(level=_ZSTD_LEVEL).compress(raw_body)
        encoding = "zstd"
    else:
        compressed_body = gzip.compress(raw_body, compresslevel=_GZIP_LEVEL)
        encoding = "gzip"

    headers = {
        "Content-Encoding": encoding,
        "Content-Type": "application/json",
    }
    return compressed_body, headers, encoding, len(raw_body), len(compressed_body)
