"""Tests for xHELM notation parser."""

from __future__ import annotations

import pytest

from helm_notation.tools import xhelm_parser
from helm_notation.tools import complex_notation_parser
from tests.conftest import load_test_resource, RESOURCES_DIR


def _get_xhelm_root_element(resource_name: str):
    """Load an xHELM resource and return the parsed root element."""
    content = load_test_resource(resource_name)
    return xhelm_parser.parse_xhelm_string(content)


class TestParseXHelmNotation:
    """Tests for parsing xHELM notation files."""

    def test_peptide_linker_nucleotide(self, monomer_factory, nucleotide_factory):
        root = _get_xhelm_root_element("PeptideLinkerNucleotide.xhelm")
        helm_string = xhelm_parser.get_complex_notation_string(root)
        store = xhelm_parser.get_monomer_store(root)

        assert helm_string == (
            "RNA1{[am6]P.R(C)P.R(U)P.R(U)P.R(G)P.R(A)P.R(G)P.R(G)}"
            "|PEPTIDE1{[aaa].C.G.K.E.D.K.R}|CHEM1{SMCC}"
            "$PEPTIDE1,CHEM1,2:R3-1:R2|RNA1,CHEM1,1:R1-1:R1$$$"
        )
        assert complex_notation_parser.validate_complex_notation(
            helm_string, store
        )

        canonical = complex_notation_parser.get_canonical_notation(
            helm_string, True, store
        )
        assert canonical == (
            "CHEM1{SMCC}|PEPTIDE1{[aaa].C.G.K.E.D.K.R}"
            "|RNA1{[am6]P.R(C)P.R(U)P.R(U)P.R(G)P.R(A)P.R(G)P.R(G)}"
            "$CHEM1,PEPTIDE1,1:R2-2:R3|CHEM1,RNA1,1:R1-1:R1$$$"
        )

    def test_simple_peptide_with_copy(self, monomer_factory, nucleotide_factory):
        root = _get_xhelm_root_element("simple.xhelm")
        helm_string = xhelm_parser.get_complex_notation_string(root)
        store = xhelm_parser.get_monomer_store(root)

        assert helm_string == "PEPTIDE1{G.K.A.[A_copy]}$$$$"
        assert complex_notation_parser.validate_complex_notation(
            helm_string, store
        )

    @pytest.mark.skipif(
        not (RESOURCES_DIR / "InlineSmiles.xhelm").exists(),
        reason="InlineSmiles.xhelm resource file not available",
    )
    def test_inline_smiles(self, monomer_factory, nucleotide_factory):
        root = _get_xhelm_root_element("InlineSmiles.xhelm")
        helm_string = xhelm_parser.get_complex_notation_string(root)
        store = xhelm_parser.get_monomer_store(root)

        assert helm_string == (
            "PEPTIDE1{A.C.A.C.G.K.E.E}|PEPTIDE2{A.C.A.C.G.K.E.E}"
            "|CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}|CHEM2{PEG2}"
            "$PEPTIDE1+PEPTIDE2,CHEM1,generic:K-1:R1"
            "|PEPTIDE1+PEPTIDE2,CHEM2,generic:Q1+Q2-1:R1$$$"
        )
        assert complex_notation_parser.validate_complex_notation(
            helm_string, store
        )

    def test_rna_with_inline(self, monomer_factory, nucleotide_factory):
        root = _get_xhelm_root_element("RNAWithInline.xhelm")
        helm_string = xhelm_parser.get_complex_notation_string(root)
        store = xhelm_parser.get_monomer_store(root)

        assert helm_string == (
            "RNA1{[C[C@@]1([*])O[C@H](CO[*])[C@@H](O[*])[C@H]1O "
            "|$;;_R3;;;;;_R1;;;_R2;;$|](A)P.RP."
            "[C[C@@]1([*])O[C@H](CO[*])[C@@H](O[*])[C@H]1O "
            "|$;;_R3;;;;;_R1;;;_R2;;$|](T)P."
            "R([Cc1nc2c(nc(N)[nH]c2=O)n1[*] |$;;;;;;;;;;;;_R1$|])P."
            "R([Cc1cc(N)nc(=O)n1[*] |$;;;;;;;;;_R1$|])}$$$$"
        )
        assert complex_notation_parser.validate_complex_notation(
            helm_string, store
        )


class TestQRPeptide:
    """Tests for QR peptide xHELM parsing."""

    def test_qr_peptide(self, monomer_factory, nucleotide_factory):
        root = _get_xhelm_root_element("qr_peptide.xhelm")
        helm_string = xhelm_parser.get_complex_notation_string(root)
        store = xhelm_parser.get_monomer_store(root)

        assert helm_string == "PEPTIDE1{[QR]}$$$$"
        assert complex_notation_parser.validate_complex_notation(
            helm_string, store
        )


class TestXHelmWithInlineSmiles:
    """Tests for xHELM with inline SMILES."""

    def test_simple_with_inline_smiles(self, monomer_factory, nucleotide_factory):
        root = _get_xhelm_root_element("simpleWithInlineSmiles.xhelm")
        helm_string = xhelm_parser.get_complex_notation_string(root)
        store = xhelm_parser.get_monomer_store(root)

        assert complex_notation_parser.validate_complex_notation(
            helm_string, store
        )
        smiles = complex_notation_parser.get_complex_polymer_smiles(
            helm_string, store
        )
        # Note: RDKit may produce different canonical SMILES than ChemAxon
        assert smiles is not None
        assert len(smiles) > 0


class TestInlineSmilesChem:
    """Tests for xHELM with inline SMILES in CHEM nodes."""

    def test_inline_smiles_chem(self, monomer_factory, nucleotide_factory):
        root = _get_xhelm_root_element("inlineSmilesChem.xhelm")
        helm_string = xhelm_parser.get_complex_notation_string(root)
        store = xhelm_parser.get_monomer_store(root)

        assert complex_notation_parser.validate_complex_notation(
            helm_string, store
        )


class TestPeptideLinkerWithChemInline:
    """Tests for peptide linker with inline CHEM xHELM."""

    def test_peptide_linker_with_chem_inline(self, monomer_factory, nucleotide_factory):
        root = _get_xhelm_root_element("PeptideLinkerWithChemInline.xhelm")
        helm_string = xhelm_parser.get_complex_notation_string(root)
        store = xhelm_parser.get_monomer_store(root)

        assert complex_notation_parser.validate_complex_notation(
            helm_string, store
        )


class TestInlineSmilesSinglePeptide:
    """Tests for xHELM with inline SMILES single peptide."""

    def test_inline_smiles_single_peptide(self, monomer_factory, nucleotide_factory):
        root = _get_xhelm_root_element("InlineSmilesSinglePeptide.xhelm")
        helm_string = xhelm_parser.get_complex_notation_string(root)
        store = xhelm_parser.get_monomer_store(root)

        assert complex_notation_parser.validate_complex_notation(
            helm_string, store
        )


class TestXHelmValidation:
    """Tests for xHELM validation of invalid documents."""

    def test_bad_xhelm_fails_validation(self, monomer_factory, nucleotide_factory):
        root = _get_xhelm_root_element("bad.xhelm")
        helm_string = xhelm_parser.get_complex_notation_string(root)
        store = xhelm_parser.get_monomer_store(root)

        with pytest.raises(Exception):
            complex_notation_parser.validate_complex_notation(
                helm_string, store
            )
