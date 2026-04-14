"""Lightweight BGZF block writer utilities (stdlib-only)."""

from __future__ import annotations

import struct
import zlib

BGZF_BLOCK_SIZE = 65280  # max uncompressed payload per BGZF block


def make_bgzf_block(data: bytes) -> bytes:
    """Compress *data* into a single BGZF block."""
    compressor = zlib.compressobj(zlib.Z_DEFAULT_COMPRESSION, zlib.DEFLATED, -15)
    compressed = compressor.compress(data) + compressor.flush()
    bsize = 18 + len(compressed) + 8 - 1
    header = (
        b"\x1f\x8b\x08\x04"
        b"\x00\x00\x00\x00"
        b"\x00\xff"
        + struct.pack("<H", 6)
        + b"BC"
        + struct.pack("<H", 2)
        + struct.pack("<H", bsize)
    )
    crc = zlib.crc32(data) & 0xFFFFFFFF
    trailer = struct.pack("<I", crc) + struct.pack("<I", len(data) & 0xFFFFFFFF)
    return header + compressed + trailer


class BgzfWriter:
    """Buffered BGZF text writer for producing block-gzipped files."""

    def __init__(self, path: str, encoding: str = "utf-8") -> None:
        self._encoding = encoding
        self._fout = open(path, "wb")
        self._buf = bytearray()

    def __enter__(self) -> "BgzfWriter":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def write(self, text: str) -> None:
        self._buf.extend(text.encode(self._encoding))
        while len(self._buf) >= BGZF_BLOCK_SIZE:
            chunk = bytes(self._buf[:BGZF_BLOCK_SIZE])
            self._buf = self._buf[BGZF_BLOCK_SIZE:]
            self._fout.write(make_bgzf_block(chunk))

    def close(self) -> None:
        if self._fout.closed:
            return
        if self._buf:
            self._fout.write(make_bgzf_block(bytes(self._buf)))
            self._buf.clear()
        self._fout.write(make_bgzf_block(b""))  # EOF block
        self._fout.close()
