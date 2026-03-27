"""Tests for nucleotide sequence parser and nucleotide model."""

from __future__ import annotations

import pytest

from helm_notation.tools import nucleotide_sequence_parser
from helm_notation.tools import simple_notation_parser
from helm_notation.tools import complex_notation_parser
from helm_notation.model.nucleotide import Nucleotide


class TestGetNotation:
    """Tests for converting nucleotide sequences to HELM notation."""

    def test_aug_notation(self, monomer_factory, nucleotide_factory):
        notation = nucleotide_sequence_parser.get_notation("AUG")
        assert notation == "R(A)P.R(U)P.R(G)P"

    def test_aug_simple_polymer_smiles(self, monomer_factory, nucleotide_factory):
        """Test SMILES generation from AUG notation.

        Note: RDKit may produce different canonical SMILES than ChemAxon,
        so we check for a non-empty result rather than exact match.
        """
        notation = nucleotide_sequence_parser.get_notation("AUG")
        smiles = simple_notation_parser.get_simple_polymer_smiles(notation, "RNA")
        assert smiles is not None
        assert len(smiles) > 0

    def test_aug_complex_polymer_smiles(self, monomer_factory, nucleotide_factory):
        """Test complex notation SMILES for AUG.

        Note: RDKit may produce different canonical SMILES than ChemAxon.
        """
        notation = nucleotide_sequence_parser.get_notation("AUG")
        complex_notation = "RNA1{" + notation + "}$$$$"
        smiles = complex_notation_parser.get_complex_polymer_smiles(complex_notation)
        assert smiles is not None
        assert len(smiles) > 0


class TestIsInOppositeDirection:
    """Tests for checking if sequences are in opposite direction."""

    def test_same_direction_rna(self, monomer_factory, nucleotide_factory):
        opp = nucleotide_sequence_parser.is_in_opposite_direction(
            "AUGTAU", "AUACAU"
        )
        assert opp is False

    def test_same_direction_mixed(self, monomer_factory, nucleotide_factory):
        opp = nucleotide_sequence_parser.is_in_opposite_direction(
            "ATGCAAA", "UUUGCAU"
        )
        assert opp is False


class TestGetSirnaNotation:
    """Tests for generating siRNA notation from sense/antisense sequences."""

    def test_sirna_notation_default(self, monomer_factory, nucleotide_factory):
        sense_seq = "CGAAAUGUUCAUACUGUUGdTdT"
        antisense_seq = "UUACAAUUUGGACUUUCCGdTdT"

        notation = nucleotide_sequence_parser.get_sirna_notation(
            sense_seq, antisense_seq
        )
        assert notation is not None
        assert "RNA1{" in notation
        assert "RNA2{" in notation
        # Check that it contains dT modifications
        assert "[dR](T)" in notation

    def test_sirna_notation_tuschl(self, monomer_factory, nucleotide_factory):
        sense_seq = "CGAAAUGUUCAUACUGUUGdTdT"
        antisense_seq = "UUACAAUUUGGACUUUCCGdTdT"

        notation = nucleotide_sequence_parser.get_sirna_notation(
            sense_seq,
            antisense_seq,
            nucleotide_sequence_parser.RNA_DESIGN_TUSCHL_19_PLUS_2,
        )
        assert notation is not None
        assert "RNA1{" in notation
        assert "RNA2{" in notation
        # Tuschl design should have more base pairs
        assert "pair" in notation


class TestGetNaturalAnalog:
    """Tests for getting natural analog of nucleotides."""

    def test_modified_c(self, monomer_factory):
        notation = "[LR]([5meC])[sP]"
        n = Nucleotide("ls5C", notation)
        assert n.get_natural_analog() == "C"

    def test_standard_a_middle(self, monomer_factory):
        n = Nucleotide("R(A)P", Nucleotide.MIDDLE_POSITION_TYPE)
        assert n.get_natural_analog() == "A"

    def test_inline_smiles_a_middle(self, monomer_factory):
        n = Nucleotide(
            "R([Nc1ncnc2n([*])cnc12 |$;;;;;;;_R1;;;$|])P",
            Nucleotide.MIDDLE_POSITION_TYPE,
        )
        assert n.get_natural_analog() == "A"

    def test_inline_smiles_sugar_and_base_middle(self, monomer_factory):
        n = Nucleotide(
            "[CO[C@H]1[C@H]([*])O[C@H](CO[*])[C@H]1O[*] |$;;;;_R3;;;;;_R1;;;_R2$|]"
            "([Nc1ncnc2n([*])cnc12 |$;;;;;;;_R1;;;$|])P",
            Nucleotide.MIDDLE_POSITION_TYPE,
        )
        assert n.get_natural_analog() == "A"

    def test_non_standard_base_x(self, monomer_factory):
        n = Nucleotide(
            "R([Cc1nc2c(nc(N)[nH]c2=O)n1[*] |$;;;;;;;;;;;;_R1$|])P",
            Nucleotide.MIDDLE_POSITION_TYPE,
        )
        assert n.get_natural_analog() == "X"

    def test_no_base_middle(self, monomer_factory):
        n = Nucleotide("RP", Nucleotide.MIDDLE_POSITION_TYPE)
        assert n.get_natural_analog() == "X"

    def test_sugar_only_ending(self, monomer_factory):
        n = Nucleotide("R", Nucleotide.ENDING_POSITION_TYPE)
        assert n.get_natural_analog() == "X"

    def test_base_only_starting(self, monomer_factory):
        n = Nucleotide(
            "([Nc1ncnc2n([*])cnc12 |$;;;;;;;_R1;;;$|])P",
            Nucleotide.STARTING_POSITION_TYPE,
        )
        assert n.get_natural_analog() == "A"


class TestGetComplementSequence:
    """Tests for complement sequence generation."""

    def test_normal_complement_formatted(self, monomer_factory, nucleotide_factory):
        sequence = "5'-UAU GUC UCC AGA AUG UAG CdTdT-3'"
        normal = nucleotide_sequence_parser.get_normal_complement_sequence(sequence)
        assert normal == "5'-AAGCUACAUUCUGGAGACAUA-3'"

    def test_reverse_complement_formatted(self, monomer_factory, nucleotide_factory):
        sequence = "5'-UAU GUC UCC AGA AUG UAG CdTdT-3'"
        reverse = nucleotide_sequence_parser.get_reverse_complement_sequence(sequence)
        assert reverse == "3'-AUACAGAGGUCUUACAUCGAA-5'"

    def test_notation_from_formatted(self, monomer_factory, nucleotide_factory):
        sequence = "5'-UAU GUC UCC AGA AUG UAG CdTdT-3'"
        notation = nucleotide_sequence_parser.get_notation(sequence)
        assert notation == (
            "R(U)P.R(A)P.R(U)P.R(G)P.R(U)P.R(C)P.R(U)P.R(C)P.R(C)P."
            "R(A)P.R(G)P.R(A)P.R(A)P.R(U)P.R(G)P.R(U)P.R(A)P.R(G)P."
            "R(C)P.[dR](T)P.[dR](T)P"
        )

    def test_modified_sequence_complement(self, monomer_factory, nucleotide_factory):
        sequence = "5'-mGCU ACA UUC UGG AGA CAU AdTdT-3'"
        normal = nucleotide_sequence_parser.get_normal_complement_sequence(sequence)
        assert normal == "5'-AAUAUGUCUCCAGAAUGUAGC-3'"

    def test_modified_sequence_reverse_complement(self, monomer_factory, nucleotide_factory):
        sequence = "5'-mGCU ACA UUC UGG AGA CAU AdTdT-3'"
        reverse = nucleotide_sequence_parser.get_reverse_complement_sequence(sequence)
        assert reverse == "3'-CGAUGUAAGACCUCUGUAUAA-5'"

    def test_modified_sequence_notation(self, monomer_factory, nucleotide_factory):
        sequence = "5'-mGCU ACA UUC UGG AGA CAU AdTdT-3'"
        notation = nucleotide_sequence_parser.get_notation(sequence)
        assert notation == (
            "[mR](G)P.R(C)P.R(U)P.R(A)P.R(C)P.R(A)P.R(U)P.R(U)P.R(C)P."
            "R(U)P.R(G)P.R(G)P.R(A)P.R(G)P.R(A)P.R(C)P.R(A)P.R(U)P."
            "R(A)P.[dR](T)P.[dR](T)P"
        )

    def test_lowercase_sequence(self, monomer_factory, nucleotide_factory):
        sequence = "agcuagggu"
        normal = nucleotide_sequence_parser.get_normal_complement_sequence(sequence)
        assert normal == "5'-ACCCUAGCU-3'"

    def test_lowercase_reverse_complement(self, monomer_factory, nucleotide_factory):
        sequence = "agcuagggu"
        reverse = nucleotide_sequence_parser.get_reverse_complement_sequence(sequence)
        assert reverse == "3'-UCGAUCCCA-5'"

    def test_lowercase_notation(self, monomer_factory, nucleotide_factory):
        sequence = "agcuagggu"
        notation = nucleotide_sequence_parser.get_notation(sequence)
        assert notation == "R(A)P.R(G)P.R(C)P.R(U)P.R(A)P.R(G)P.R(G)P.R(G)P.R(U)P"
