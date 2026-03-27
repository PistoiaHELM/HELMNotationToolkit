"""Molfile encoding/decoding using gzip+base64 (replacing Pfizer sdlib)."""

from __future__ import annotations

import base64
import gzip


def encode(molfile: str) -> str | None:
    """Encode a molfile string using gzip compression and base64 encoding."""
    if molfile is None:
        return None
    compressed = gzip.compress(molfile.encode("utf-8"))
    return base64.b64encode(compressed).decode("ascii")


def decode(encoded_string: str) -> str | None:
    """Decode a gzip+base64 encoded molfile string."""
    if encoded_string is None:
        return None
    compressed = base64.b64decode(encoded_string.encode("ascii"))
    return gzip.decompress(compressed).decode("utf-8")
