"""Tests for molfile encoder."""

from __future__ import annotations

from helm_notation.tools import molfile_encoder


class TestEncode:
    """Tests for molfile encoding and decoding round-trip."""

    def test_encode_decode_round_trip(self):
        molfile = (
            "\n  Marvin  08200815002D          \n\n"
            " 10  9  0  0  0  0            999 V2000\n"
            "   -2.7541    2.1476    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "   -3.4686    0.9102    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "   -3.4686    1.7352    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "   -4.1830    2.1477    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "   -2.0397    1.7351    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n"
            "   -1.3250    2.1476    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "   -1.3250    2.9726    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "   -2.0397    0.9101    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "   -0.6105    1.7350    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n"
            "   -2.7542    0.4976    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n"
            "  3  1  1  0  0  0  0\n  5  1  1  1  0  0  0\n  3  2  1  0  0  0  0\n"
            "  3  4  1  0  0  0  0\n  5  8  1  0  0  0  0\n  5  6  1  0  0  0  0\n"
            "  6  7  2  0  0  0  0\n  6  9  1  0  0  0  0\n  8 10  1  0  0  0  0\n"
            "M  RGP  2   9   2  10   1\nM  END\n"
        )

        encoded = molfile_encoder.encode(molfile)
        assert encoded is not None
        assert len(encoded) > 0

        decoded = molfile_encoder.decode(encoded)
        assert decoded == molfile
