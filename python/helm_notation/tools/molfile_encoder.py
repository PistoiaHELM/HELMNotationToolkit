"""Port of org.helm.notation.tools.MolfileEncoder.

The Java version uses sdlib (gzip-encoded SD blocks). We replicate the
same gzip+base64 wire format so existing MonomerDBGZEncoded.xml files
can be loaded without conversion.
"""

from __future__ import annotations

import base64
import gzip


def encode(text: str) -> str:
    """GZip-compress *text* and return it as a Base64 string."""
    if text is None:
        return None  # type: ignore[return-value]
    compressed = gzip.compress(text.encode("utf-8"))
    return base64.b64encode(compressed).decode("ascii")


def decode(encoded: str) -> str:
    """Decode a Base64+GZip string back to the original text."""
    if encoded is None:
        return None  # type: ignore[return-value]
    raw = base64.b64decode(encoded)
    return gzip.decompress(raw).decode("utf-8")
