"""Tests for complex notation parser."""

from __future__ import annotations

import pytest

from helm_notation.tools import complex_notation_parser
from helm_notation.tools import simple_notation_parser
from helm_notation.tools import structure_parser
from helm_notation.model.monomer import Monomer


# --- Helper functions matching Java test helpers ---


def _test_complex_notation_validity(notation: str) -> bool:
    """Test that a complex notation is valid (returns True on success)."""
    try:
        complex_notation_parser.validate_notation_format(notation)
        complex_notation_parser.get_all_node_string(notation)
        complex_notation_parser.get_all_edge_string(notation)
        complex_notation_parser.get_all_base_pair_string(notation)
        complex_notation_parser.get_all_node_label_string(notation)
        complex_notation_parser.get_other_string(notation)
        complex_notation_parser.validate_complex_notation(notation)
        return True
    except Exception:
        return False


def _test_get_canonical_notation(notation: str) -> bool:
    """Test that canonical notation can be generated."""
    try:
        complex_notation_parser.get_canonical_notation(notation)
        return True
    except Exception:
        return False


def _test_get_molecule_info(notation: str) -> bool:
    """Test that molecule info can be calculated."""
    try:
        mi = complex_notation_parser.get_molecule_info(notation)
        assert mi.get_molecular_weight() > 0
        return True
    except Exception:
        return False


def _test_get_molecule_info_via_smiles(notation: str) -> bool:
    """Test that molecule info can be calculated via SMILES."""
    try:
        smiles = complex_notation_parser.get_complex_polymer_smiles(notation)
        mi = structure_parser.get_molecule_info(smiles)
        assert mi.get_molecular_weight() > 0
        return True
    except Exception:
        return False


def _create_peptide_notation(monomer_count: int) -> str:
    """Create a peptide notation with the given number of alanine monomers."""
    monomers = ".".join(["A"] * monomer_count)
    return f"PEPTIDE1{{{monomers}}}$$$$"


class TestGenericConnection:
    """Tests for generic connections in complex notation."""

    def test_peptide_chem_generic(self, monomer_factory, nucleotide_factory):
        notation = (
            "PEPTIDE1{A.C.A.C.G.K}|CHEM1{[*]OCCOCCOCCO[*] "
            "|$_R1;;;;;;;;;;;_R3$|}$PEPTIDE1,CHEM1,generic:K-1:R1$$$"
        )
        assert _test_complex_notation_validity(notation)
        assert _test_get_canonical_notation(notation)
        assert _test_get_molecule_info(notation)

    def test_peptide_peptide_generic_standard(self, monomer_factory, nucleotide_factory):
        notation = (
            "PEPTIDE1{A.C.A.C.G.K.E.E}|PEPTIDE2{A.K.D.C.A}"
            "$PEPTIDE1,PEPTIDE2,generic:K-3:R3$$$"
        )
        assert _test_complex_notation_validity(notation)
        assert _test_get_canonical_notation(notation)
        assert _test_get_molecule_info(notation)

    def test_peptide_peptide_two_generics(self, monomer_factory, nucleotide_factory):
        notation = (
            "PEPTIDE1{A.C.A.C.G.K.E.E}|PEPTIDE2{A.K.D.C.A}"
            "$PEPTIDE1,PEPTIDE2,generic:K-generic:E$$$"
        )
        assert _test_complex_notation_validity(notation)
        assert _test_get_canonical_notation(notation)
        assert _test_get_molecule_info(notation)

    def test_peptide_group_chem_generic(self, monomer_factory, nucleotide_factory):
        notation = (
            "PEPTIDE1{A.C.A.C.G.K.E.E}|PEPTIDE2{A.C.A.C.G.K.E.E}"
            "|CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}"
            "$PEPTIDE1+PEPTIDE2,CHEM1,generic:K2+K5-1:R1$$$"
        )
        assert _test_complex_notation_validity(notation)
        assert _test_get_canonical_notation(notation)
        assert _test_get_molecule_info(notation)

    def test_peptide_group_two_chem_generics(self, monomer_factory, nucleotide_factory):
        notation = (
            "PEPTIDE1{A.C.A.C.G.K.E.E}|PEPTIDE2{A.C.A.C.G.K.E.E}"
            "|CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}|CHEM2{PEG2}"
            "$PEPTIDE1+PEPTIDE2,CHEM1,generic:K-1:R1"
            "|PEPTIDE1+PEPTIDE2,CHEM2,generic:Q1+Q2-1:R1$$$"
        )
        assert _test_complex_notation_validity(notation)
        assert _test_get_canonical_notation(notation)
        assert _test_get_molecule_info(notation)

    def test_two_peptide_groups_generic(self, monomer_factory, nucleotide_factory):
        notation = (
            "PEPTIDE1{A.C.A.K.C.G.K.E.E}|PEPTIDE2{G.A.C.A.C.G.K.E.E}"
            "|PEPTIDE3{A.C.A.C.G.K.E.E}|PEPTIDE4{A.C.A.E.D.C.G.K.E.E}"
            "$PEPTIDE1+PEPTIDE4,PEPTIDE3+PEPTIDE2,generic:K-generic:D$$$"
        )
        assert _test_complex_notation_validity(notation)
        assert _test_get_canonical_notation(notation)
        assert _test_get_molecule_info(notation)

    def test_two_unknown_connections(self, monomer_factory, nucleotide_factory):
        notation = (
            "PEPTIDE1{A.A.A.A.A.A.A.A}|PEPTIDE2{G.G.G.G.G.G.G.G.G.G.G.G.G.G.G.G}"
            "|CHEM1{C1=CC=CC=C1 |c:0,2,4|}"
            "$PEPTIDE1+PEPTIDE2,CHEM1,generic:?-1:?$$$"
        )
        assert _test_complex_notation_validity(notation)
        assert _test_get_canonical_notation(notation)
        assert _test_get_molecule_info(notation)

    def test_one_unknown_connection(self, monomer_factory, nucleotide_factory):
        notation = (
            "PEPTIDE1{A.A.A.A.A.A.A.A}|PEPTIDE2{G.G.G.G.G.G.G.G.G.G.G.G.G.G.G.G}"
            "|CHEM1{[*]C1=CC=CC=C1 |$_R8;;;;;;$,c:3,5,t:1|}"
            "$PEPTIDE1+PEPTIDE2,CHEM1,generic:?-1:R8$$$"
        )
        assert _test_complex_notation_validity(notation)
        assert _test_get_canonical_notation(notation)
        assert _test_get_molecule_info(notation)


class TestCanonicalization:
    """Tests for canonical notation generation with various orderings."""

    def test_standard_order(self, monomer_factory, nucleotide_factory):
        notation = (
            "PEPTIDE1{H.H.E.E.E}|CHEM1{SS3}|CHEM2{EG}"
            "$PEPTIDE1,CHEM2,5:R2-1:R2|CHEM2,CHEM1,1:R1-1:R2"
            "|PEPTIDE1,CHEM1,1:R1-1:R1$$$"
        )
        assert _test_get_canonical_notation(notation)

    def test_changed_node_order(self, monomer_factory, nucleotide_factory):
        notation = (
            "CHEM1{SS3}|PEPTIDE1{H.H.E.E.E}|CHEM2{EG}"
            "$PEPTIDE1,CHEM2,5:R2-1:R2|CHEM2,CHEM1,1:R1-1:R2"
            "|PEPTIDE1,CHEM1,1:R1-1:R1$$$"
        )
        assert _test_get_canonical_notation(notation)

    def test_changed_edge_order(self, monomer_factory, nucleotide_factory):
        notation = (
            "CHEM1{SS3}|PEPTIDE1{H.H.E.E.E}|CHEM2{EG}"
            "$CHEM2,CHEM1,1:R1-1:R2|PEPTIDE1,CHEM2,5:R2-1:R2"
            "|PEPTIDE1,CHEM1,1:R1-1:R1$$$"
        )
        assert _test_get_canonical_notation(notation)

    def test_changed_edge_direction(self, monomer_factory, nucleotide_factory):
        notation = (
            "CHEM1{SS3}|PEPTIDE1{H.H.E.E.E}|CHEM2{EG}"
            "$CHEM2,CHEM1,1:R1-1:R2|CHEM2,PEPTIDE1,1:R2-5:R2"
            "|PEPTIDE1,CHEM1,1:R1-1:R1$$$"
        )
        assert _test_get_canonical_notation(notation)

    def test_changed_node_id(self, monomer_factory, nucleotide_factory):
        notation = (
            "CHEM4{SS3}|PEPTIDE2{H.H.E.E.E}|CHEM2{EG}"
            "$CHEM2,CHEM4,1:R1-1:R2|CHEM2,PEPTIDE2,1:R2-5:R2"
            "|PEPTIDE2,CHEM4,1:R1-1:R1$$$"
        )
        assert _test_get_canonical_notation(notation)

    def test_backbone_cyclic_peptide(self, monomer_factory, nucleotide_factory):
        notation = "PEPTIDE1{K.A.A.G.K}$PEPTIDE1,PEPTIDE1,1:R1-5:R2$$$"
        assert _test_get_canonical_notation(notation)

    def test_annotated_peptides(self, monomer_factory, nucleotide_factory):
        notation = (
            "PEPTIDE1{K.A.A.G.K}|PEPTIDE2{K.A.A.G.K}|RNA1{R(A)P.R(G)}|CHEM1{Alexa}"
            "$PEPTIDE1,PEPTIDE1,1:R1-5:R2$$PEPTIDE1{hc}|PEPTIDE2{lc}$"
        )
        assert _test_get_canonical_notation(notation)

    def test_unregistered_chem_fails(self, monomer_factory, nucleotide_factory):
        """CovX-2 is not a registered monomer, so canonicalization should fail."""
        notation = (
            "RNA1{R(A)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P.R(G)P."
            "R(G)P.R(G)P.R(C)P.R(C)P.R(C)P.R(C)}|CHEM1{CovX-2}"
            "$RNA1,CHEM1,1:R1-1:R1$$$"
        )
        assert not _test_get_canonical_notation(notation)


class TestMoleculeInfo:
    """Tests for molecule info calculation from complex notation."""

    def test_sirna_molecule_info(self, monomer_factory, nucleotide_factory):
        notation = (
            "RNA1{R(A)P.R(U)P.R(C)P.R(C)P.R(A)P.R(A)P.R(A)P.R(G)P.R(A)P."
            "R(U)P.R(A)P.R(C)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P."
            "R(G)P.R(C)P.R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)}"
            "|RNA2{R(U)P.R(U)P.R(C)P.R(U)P.R(G)P.R(C)P.R(A)P.R(A)P.R(A)P."
            "R(G)P.R(C)P.R(U)P.R(A)P.R(G)P.R(U)P.R(A)P.R(U)P.R(C)P.R(U)P."
            "R(U)P.R(U)P.R(G)P.R(G)P.[dR](A)P.[dR](T)}"
            "$$RNA1,RNA2,2:pair-74:pair|RNA1,RNA2,5:pair-71:pair"
            "|RNA1,RNA2,8:pair-68:pair|RNA1,RNA2,11:pair-65:pair"
            "|RNA1,RNA2,14:pair-62:pair|RNA1,RNA2,17:pair-59:pair"
            "|RNA1,RNA2,20:pair-56:pair|RNA1,RNA2,23:pair-53:pair"
            "|RNA1,RNA2,26:pair-50:pair|RNA1,RNA2,29:pair-47:pair"
            "|RNA1,RNA2,32:pair-44:pair|RNA1,RNA2,35:pair-41:pair"
            "|RNA1,RNA2,38:pair-38:pair|RNA1,RNA2,41:pair-35:pair"
            "|RNA1,RNA2,44:pair-32:pair|RNA1,RNA2,47:pair-29:pair"
            "|RNA1,RNA2,50:pair-26:pair|RNA1,RNA2,53:pair-23:pair"
            "|RNA1,RNA2,56:pair-20:pair|RNA1,RNA2,59:pair-17:pair"
            "|RNA1,RNA2,62:pair-14:pair|RNA1,RNA2,65:pair-11:pair"
            "|RNA1,RNA2,68:pair-8:pair|RNA1,RNA2,71:pair-5:pair"
            "|RNA1,RNA2,74:pair-2:pair$$"
        )
        assert _test_get_molecule_info(notation)
        assert _test_get_molecule_info_via_smiles(notation)

    def test_conjugate_molecule_info(self, monomer_factory, nucleotide_factory):
        notation = (
            "PEPTIDE1{A.G.G.G.C.C.K.K.K.K}|CHEM1{MCC}"
            "$PEPTIDE1,CHEM1,10:R3-1:R1$$$"
        )
        assert _test_get_molecule_info(notation)
        assert _test_get_molecule_info_via_smiles(notation)


class TestProteinMoleculeInfoCalculation:
    """Tests for molecule info calculation with large peptides."""

    @pytest.mark.parametrize("count", [100, 150, 199])
    def test_small_peptide_both_methods(self, monomer_factory, nucleotide_factory, count):
        notation = _create_peptide_notation(count)
        assert _test_get_molecule_info(notation)
        assert _test_get_molecule_info_via_smiles(notation)

    @pytest.mark.parametrize("count", [400, 800, 1600, 3200])
    def test_large_peptide_monomer_method(self, monomer_factory, nucleotide_factory, count):
        """Large peptides should work via monomer-based calculation.

        SMILES-based calculation may fail for very large molecules.
        """
        notation = _create_peptide_notation(count)
        assert _test_get_molecule_info(notation)


class TestNotationManipulation:
    """Tests for notation manipulation operations."""

    def test_rna_polymer_node_list(self, monomer_factory, nucleotide_factory):
        notation = (
            "RNA1{R(G)P.[fR](U)P.R(A)P.R(A)P.R(G)P.R(A)P.[fR](C)P.[fR](U)P."
            "[fR](A)P.R(C)P.R(U)P.R(C)P.R(A)P.[fR](U)P.R(G)P.R(A)P.[fR](U)P."
            "[fR](C)P.[fR](C)P.R(T)[sP].R(T)}"
            "|RNA2{R(G)P.R(G)P.R(A)P.[fR](U)P.[fR](C)P.R(A)P.[fR](U)P.[fR](G)P."
            "[fR](A)P.[fR](G)P.R(U)P.R(A)P.R(G)P.[fR](U)P.[fR](C)P.[fR](U)P."
            "[fR](U)P.R(A)P.[fR](C)P.R(T)[sP].R(T)}"
            "$$RNA2,RNA1,29:pair-29:pair|RNA2,RNA1,17:pair-41:pair"
            "|RNA2,RNA1,38:pair-20:pair|RNA2,RNA1,47:pair-11:pair"
            "|RNA2,RNA1,26:pair-32:pair|RNA2,RNA1,11:pair-47:pair"
            "|RNA2,RNA1,56:pair-2:pair|RNA2,RNA1,2:pair-56:pair"
            "|RNA2,RNA1,8:pair-50:pair|RNA2,RNA1,41:pair-17:pair"
            "|RNA2,RNA1,44:pair-14:pair|RNA2,RNA1,35:pair-23:pair"
            "|RNA2,RNA1,5:pair-53:pair|RNA2,RNA1,14:pair-44:pair"
            "|RNA2,RNA1,32:pair-26:pair|RNA2,RNA1,50:pair-8:pair"
            "|RNA2,RNA1,20:pair-38:pair|RNA2,RNA1,53:pair-5:pair"
            "|RNA2,RNA1,23:pair-35:pair$RNA1{as}|RNA2{ss}$"
        )
        rna_list = complex_notation_parser.get_rna_polymer_node_list(notation)
        assert len(rna_list) == 2

        node1 = rna_list[0]
        assert node1.get_anotation() == "as"
        assert node1.get_sequence() == "GUAAGACUACUCAUGAUCCTT"
        assert node1.get_modified_sequence() == "GfUAAGAfCfUfACUCAfUGAfUfCfCTT"

        node2 = rna_list[1]
        assert node2.get_anotation() == "ss"
        assert node2.get_sequence() == "GGAUCAUGAGUAGUCUUACTT"
        assert node2.get_modified_sequence() == "GGAfUfCAfUfGfAfGUAGfUfCfUfUAfCTT"

    def test_formatted_sirna_sequences(self, monomer_factory, nucleotide_factory):
        notation = (
            "RNA1{R(U)P.R(U)P.R(A)P.R(A)P.R(G)P.R(C)P.R(U)P.[dR](T)P.[dR](T)}"
            "|RNA2{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(A)P.R(A)P.[dR](T)P.[dR](T)}"
            "$$RNA2,RNA1,20:pair-2:pair|RNA2,RNA1,5:pair-17:pair"
            "|RNA2,RNA1,2:pair-20:pair|RNA2,RNA1,11:pair-11:pair"
            "|RNA2,RNA1,17:pair-5:pair|RNA2,RNA1,14:pair-8:pair"
            "|RNA2,RNA1,8:pair-14:pair$RNA1{ss}$"
        )
        seqs = complex_notation_parser.get_formated_sirna_sequences(
            notation, "*", "|"
        )
        assert seqs[0] == "**UUAAGCUTT"
        assert seqs[1] == "**|||||||**"
        assert seqs[2] == "TTAAUUCGA**"

    def test_formatted_sirna_sequences_default(self, monomer_factory, nucleotide_factory):
        notation = (
            "RNA1{R(U)P.R(U)P.R(A)P.R(A)P.R(G)P.R(C)P.R(U)P.[dR](T)P.[dR](T)}"
            "|RNA2{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(A)P.R(A)P.[dR](T)P.[dR](T)}"
            "$$RNA2,RNA1,20:pair-2:pair|RNA2,RNA1,5:pair-17:pair"
            "|RNA2,RNA1,2:pair-20:pair|RNA2,RNA1,11:pair-11:pair"
            "|RNA2,RNA1,17:pair-5:pair|RNA2,RNA1,14:pair-8:pair"
            "|RNA2,RNA1,8:pair-14:pair$RNA1{ss}$"
        )
        seqs = complex_notation_parser.get_formated_sirna_sequences(notation)
        assert seqs[0] == "  UUAAGCUTT"
        assert seqs[1] == "  |||||||  "
        assert seqs[2] == "TTAAUUCGA  "

    def test_combined_complex_notation(self, monomer_factory, nucleotide_factory):
        notation1 = (
            "RNA1{R(A)P.R(G)P.R(C)P.R(U)}|RNA2{R(A)P.R(G)P.R(C)P.R(U)}"
            "|RNA3{R(G)P.R(G)}|CHEM1{PEG2}"
            "$RNA3,CHEM1,4:R2-1:R2|RNA1,CHEM1,10:R2-1:R1$$$"
        )
        notation2 = (
            "CHEM1{PEG2}|RNA1{R(A)P.R(G)P.R(C)P.R(U)}|RNA2{R(A)P.R(G)P.R(C)P.R(U)}"
            "|RNA3{R(G)P.R(G)}"
            "$CHEM1,RNA1,1:R1-10:R2|CHEM1,RNA3,1:R2-4:R2$$$"
        )
        combined = complex_notation_parser.get_combined_complex_notation(
            notation1, notation2
        )
        assert combined is not None
        # Combined should have all nodes from both notations
        assert "RNA1{" in combined
        assert "CHEM1{" in combined
        assert "CHEM2{" in combined

    def test_hybridize(self, monomer_factory, nucleotide_factory):
        notation = (
            "RNA1{P.R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(U)}"
            "|RNA2{P.R(A)P.R(A)P.R(A)P.R(G)P.R(C)P.R(U)}$$$$"
        )
        result = complex_notation_parser.hybridize(notation)
        assert result is not None
        assert "pair" in result

    def test_decompose(self, monomer_factory, nucleotide_factory):
        notation = (
            "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(U)P}"
            "|RNA2{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(U)}"
            "$RNA1,RNA1,21:R2-1:R1$$$"
        )
        notations = complex_notation_parser.decompose(notation)
        assert len(notations) >= 2
        assert (
            notations[0]
            == "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(U)P}$RNA1,RNA1,21:R2-1:R1$$$"
        )
        assert (
            notations[1]
            == "RNA2{R(A)P.R(G)P.R(C)P.R(U)P.R(U)P.R(U)P.R(U)}$$$$"
        )

    def test_polymer_node_list(self, monomer_factory, nucleotide_factory):
        notation = (
            "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)}"
            "|RNA2{R(C)P.R(C)P.R(U)P.R(U)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)}"
            "$$RNA1,RNA2,2:pair-26:pair|RNA1,RNA2,14:pair-14:pair"
            "|RNA1,RNA2,5:pair-23:pair|RNA1,RNA2,26:pair-2:pair"
            "|RNA1,RNA2,23:pair-5:pair|RNA1,RNA2,20:pair-8:pair"
            "|RNA1,RNA2,11:pair-17:pair|RNA1,RNA2,17:pair-11:pair"
            "|RNA1,RNA2,8:pair-20:pair$RNA2{as}|RNA1{ss}$"
        )
        all_node = complex_notation_parser.get_all_node_string(notation)
        pn_list = complex_notation_parser.get_polymer_node_list(all_node)
        assert len(pn_list) == 2

        rna_list = complex_notation_parser.get_rna_polymer_node_list(notation)
        assert len(rna_list) == len(pn_list)

        rna = rna_list[0]
        annotation = rna.get_anotation()
        assert annotation == annotation.strip()
        assert rna.get_label() == (
            "R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)"
        )

    def test_replace_monomer_in_complex(self, monomer_factory, nucleotide_factory):
        notation = (
            "RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P.R(A)P.R(G)P.R(G)}"
            "|RNA2{R(C)P.R(C)P.R(U)P.R(U)P.R(U)P.R(A)P.R(G)P.R(C)P.R(U)}"
            "$$RNA1,RNA2,2:pair-26:pair|RNA1,RNA2,14:pair-14:pair"
            "|RNA1,RNA2,5:pair-23:pair|RNA1,RNA2,26:pair-2:pair"
            "|RNA1,RNA2,23:pair-5:pair|RNA1,RNA2,20:pair-8:pair"
            "|RNA1,RNA2,11:pair-17:pair|RNA1,RNA2,17:pair-11:pair"
            "|RNA1,RNA2,8:pair-20:pair$RNA2{as}|RNA1{ss}$"
        )
        result = complex_notation_parser.replace_monomer(
            notation, Monomer.NUCLIEC_ACID_POLYMER_TYPE, "P", "sP"
        )
        assert result is not None
        assert "[sP]" in result

    def test_has_nucleotide_modification_true(self, monomer_factory, nucleotide_factory):
        notation = (
            "RNA1{[5A6][sP].R(G)P.R(U)P.R(C)P.R(A)[sP].R(U)P.R(C)P.R(A)P."
            "R(C)P.R(A)P.R(C)P.R(U)P.R(G)P.R(A)P.R(A)P.R(U)P.R(A)P.R(C)P."
            "R(C)P.R(A)P.R(A)[sP].R(U)}"
            "|RNA2{R(A)P.R(U)P.R(U)P.R(G)P.R(G)P.R(U)P.R(A)P.R(U)P.R(U)P."
            "R(C)P.R(A)P.R(G)P.R(U)P.R(G)P.R(U)P.R(G)P.R(A)P.R(U)P.R(G)P."
            "R(A)P.[mR](C)[sP].[mR](A)[sP].R(C)}"
            "$$$RNA1{ss}|RNA2{as}$"
        )
        poly_node_list = complex_notation_parser.get_polymer_node_list(notation)
        assert complex_notation_parser.has_nucleotide_modification(poly_node_list)

    def test_has_nucleotide_modification_false(self, monomer_factory, nucleotide_factory):
        notation = (
            "RNA1{R(G)P.R(U)P.R(C)P.R(A)P.R(U)P.R(C)P.R(A)P.R(C)P.R(A)P."
            "R(C)P.R(U)P.R(G)P.R(A)P.R(A)P.R(U)P.R(A)P.R(C)P.R(C)P.R(A)P."
            "R(A)P.R(U)P}"
            "|RNA2{R(A)P.R(U)P.R(U)P.R(G)P.R(G)P.R(U)P.R(A)P.R(U)P.R(U)P."
            "R(C)P.R(A)P.R(G)P.R(U)P.R(G)P.R(U)P.R(G)P.R(A)P.R(U)P.R(G)P."
            "R(A)P.R(C)P}"
            "$$$RNA1{ss}|RNA2{as}$"
        )
        poly_node_list = complex_notation_parser.get_polymer_node_list(notation)
        assert not complex_notation_parser.has_nucleotide_modification(poly_node_list)

    def test_standardize(self, monomer_factory, nucleotide_factory):
        """Test standardize for various notation types."""
        test_cases = [
            "R(A)P.R(G)",
            "A.A.G.L",
            "PEG2",
            "RNA1{P.R(A)P.R(G)}$$$$",
        ]
        for notation in test_cases:
            result = complex_notation_parser.standardize(notation)
            assert result is not None


class TestInlineNotation:
    """Tests for inline SMILES in complex notation."""

    def test_inline_smiles_same_as_registered(self, monomer_factory, nucleotide_factory):
        """Replacing A with its inline SMILES should produce same result."""
        notation_registered = (
            "PEPTIDE1{A.G.G.G.C.C.K.K.K.K}|CHEM1{MCC}"
            "$PEPTIDE1,CHEM1,10:R3-1:R1$$$"
        )
        notation_inline = (
            "PEPTIDE1{[C[C@H](N[*])C([*])=O |$;;;_R1;;_R2;$|].G.G.G.C.C.K.K.K.K}"
            "|CHEM1{MCC}$PEPTIDE1,CHEM1,10:R3-1:R1$$$"
        )
        smiles_reg = complex_notation_parser.get_complex_polymer_smiles(
            notation_registered
        )
        smiles_inline = complex_notation_parser.get_complex_polymer_smiles(
            notation_inline
        )
        assert smiles_reg == smiles_inline

    def test_inline_smiles_validation(self, monomer_factory, nucleotide_factory):
        notation = (
            "PEPTIDE1{[C[C@H](N[*])C([*])=O |$;;;_R1;;_R2;$|].G.G.G.C.C.K.K.K.K}"
            "|CHEM1{MCC}$PEPTIDE1,CHEM1,10:R3-1:R1$$$"
        )
        assert complex_notation_parser.validate_complex_notation(notation)
        assert complex_notation_parser.validate_notation_format(notation)

    def test_modified_inline_smiles(self, monomer_factory, nucleotide_factory):
        """Slightly modified A should still be valid."""
        notation = (
            "PEPTIDE1{[C[C@H](N[*])C(=O)C[*] |$;;;_R1;;;;_R2$|].G.G.G.C.C.K.K.K.K}"
            "|CHEM1{MCC}$PEPTIDE1,CHEM1,10:R3-1:R1$$$"
        )
        assert complex_notation_parser.validate_notation_format(notation)
        assert complex_notation_parser.validate_complex_notation(notation)

    def test_modified_inline_smiles_produces_smiles(self, monomer_factory, nucleotide_factory):
        notation = (
            "PEPTIDE1{[C[C@H](N[*])C(=O)C[*] |$;;;_R1;;;;_R2$|].G.G.G.C.C.K.K.K.K}"
            "|CHEM1{MCC}$PEPTIDE1,CHEM1,10:R3-1:R1$$$"
        )
        smiles = complex_notation_parser.get_complex_polymer_smiles(notation)
        assert smiles is not None
        assert len(smiles) > 0


class TestReplaceSmiles:
    """Tests for replacing inline SMILES in complex notation."""

    def test_replace_smiles_in_complex(self, monomer_factory, monomer_store):
        notation = (
            "PEPTIDE1{[C[C@H](N[*])C(=O)C[*] |$;;;_R1;;;;_R2$|].G.G.G.C.C.K.K.K.K}"
            "|CHEM1{MCC}$PEPTIDE1,CHEM1,10:R3-1:R1$$$"
        )
        result = complex_notation_parser.get_notation_by_replacing_smiles(
            notation, monomer_store
        )
        assert result == (
            "PEPTIDE1{[PM#1].G.G.G.C.C.K.K.K.K}|CHEM1{MCC}"
            "$PEPTIDE1,CHEM1,10:R3-1:R1$$$"
        )
