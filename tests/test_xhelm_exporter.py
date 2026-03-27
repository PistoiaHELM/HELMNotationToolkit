"""Tests for xHELM notation exporter."""

from __future__ import annotations

from helm_notation.tools import xhelm_parser
from helm_notation.tools import xhelm_exporter
from tests.conftest import load_test_resource


class TestExport:
    """Tests for xHELM export round-trip."""

    def test_export_round_trip(self, monomer_factory):
        """Test that importing and re-exporting produces consistent output.

        Load a .xhelm file, parse it, export it back to xHELM, and verify
        the output length matches the original.
        """
        xhelm_content = load_test_resource("PeptideLinkerNucleotide.xhelm")

        # Parse the xHELM content
        root_element = xhelm_parser.parse_xhelm_string(xhelm_content)
        helm_string = xhelm_parser.get_complex_notation_string(root_element)
        monomer_store = xhelm_parser.get_monomer_store(root_element)

        # Export back to xHELM
        output = xhelm_exporter.write_xhelm(helm_string, monomer_store)

        assert output is not None
        assert len(output) > 0
        # The round-trip output should be approximately the same size
        assert len(output) == len(xhelm_content)
