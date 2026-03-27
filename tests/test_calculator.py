"""Tests for extinction coefficient calculator."""

from __future__ import annotations

import pytest

from helm_notation.tools import extinction_coefficient_calculator


class TestCalculateFromAminoAcidSequence:
    """Tests for calculating extinction coefficient from amino acid sequences."""

    def test_no_absorbing_residues(self, monomer_factory):
        calc = extinction_coefficient_calculator.get_instance()
        result = calc.calculate_from_amino_acid_sequence(
            "AGGDDDDDDDDDDDDDDDDDDFFFFFFFFFFFFF"
        )
        assert result == pytest.approx(0.0, abs=1e-15)

    def test_cysteine_only(self, monomer_factory):
        calc = extinction_coefficient_calculator.get_instance()
        result = calc.calculate_from_amino_acid_sequence("AGGCFFFFFFFFFF")
        assert result == pytest.approx(0.0, abs=62.5)

    def test_tyrosine(self, monomer_factory):
        calc = extinction_coefficient_calculator.get_instance()
        result = calc.calculate_from_amino_acid_sequence(
            "AGGYEEEEEEEEEEEEEEEEEEE"
        )
        assert result == pytest.approx(1490.0, abs=1e-15)

    def test_tryptophan(self, monomer_factory):
        calc = extinction_coefficient_calculator.get_instance()
        result = calc.calculate_from_amino_acid_sequence(
            "AGGWEEEEEEEEEEEEEEEEEEE"
        )
        assert result == pytest.approx(5500.0, abs=1e-15)


class TestCalculateFromPeptidePolymerNotation:
    """Tests for calculating extinction coefficient from peptide polymer notation."""

    def test_peptide_notation(self, monomer_factory):
        calc = extinction_coefficient_calculator.get_instance()
        result = calc.calculate_from_peptide_polymer_notation(
            "A.G.G.W.E.E.E.E.E.W"
        )
        assert result == pytest.approx(11000.0, abs=1e-15)


class TestCalculateFromComplexNotation:
    """Tests for calculating extinction coefficient from complex notation."""

    def test_single_peptide(self, monomer_factory):
        calc = extinction_coefficient_calculator.get_instance()
        notation = "PEPTIDE1{A.G.G.W.E.E.E.E.E.W}$$$$"
        result = calc.calculate_from_complex_notation(
            notation,
            extinction_coefficient_calculator.PEPTIDE_UNIT_TYPE,
        )
        assert result == pytest.approx(11000.0, abs=1e-15)

    def test_two_peptides_default_unit(self, monomer_factory):
        calc = extinction_coefficient_calculator.get_instance()
        notation = "PEPTIDE1{A.G.G.W.E.E.E.E.E.W}|PEPTIDE2{A.G.G.W.E.Y.E.E.E.E.W}$$$$"
        result = calc.calculate_from_complex_notation(notation)
        assert result == pytest.approx(23.49, abs=0.01)

    def test_two_peptides_peptide_unit(self, monomer_factory):
        calc = extinction_coefficient_calculator.get_instance()
        notation = "PEPTIDE1{A.G.G.W.E.E.E.E.E.W}|PEPTIDE2{A.G.G.W.E.Y.E.E.E.E.W}$$$$"
        result = calc.calculate_from_complex_notation(
            notation,
            extinction_coefficient_calculator.PEPTIDE_UNIT_TYPE,
        )
        assert result == pytest.approx(23490.0, abs=1e-15)

    def test_rna_default_unit(self, monomer_factory):
        calc = extinction_coefficient_calculator.get_instance()
        notation = "RNA1{P.R(A)P.R([5meC])P.R(G)P.[mR](A)}$$$$"
        result = calc.calculate_from_complex_notation(notation)
        assert result == pytest.approx(46.200005, abs=0.01)

    def test_rna_peptide_unit(self, monomer_factory):
        calc = extinction_coefficient_calculator.get_instance()
        notation = "RNA1{P.R(A)P.R([5meC])P.R(G)P.[mR](A)}$$$$"
        result = calc.calculate_from_complex_notation(
            notation,
            extinction_coefficient_calculator.PEPTIDE_UNIT_TYPE,
        )
        assert result == pytest.approx(46200.004, abs=0.01)

    def test_mixed_rna_chem_peptide(self, monomer_factory):
        calc = extinction_coefficient_calculator.get_instance()
        notation = (
            "RNA1{P.R(A)P.R([5meC])P.R(G)P.[mR](A)}"
            "|CHEM1{PEG2}"
            "|PEPTIDE1{A.G.G.W.E.E.E.E.E.W}"
            "|PEPTIDE2{A.G.G.W.E.Y.E.E.E.E.W}$$$$"
        )
        result = calc.calculate_from_complex_notation(notation)
        assert result == pytest.approx(69.69, abs=0.01)


class TestCalculateFromNucleotideSequence:
    """Tests for calculating extinction coefficient from nucleotide sequences."""

    def test_nucleotide_sequence(self, monomer_factory):
        calc = extinction_coefficient_calculator.get_instance()
        result = calc.calculate_from_nucleotide_sequence("ACGTACGT")
        assert result == pytest.approx(81.119995, abs=0.01)


class TestCalculateFromModifiedNucleotideSequence:
    """Tests for calculating extinction coefficient from modified nucleotide sequences."""

    def test_modified_nucleotide_sequence(self, monomer_factory, nucleotide_factory):
        calc = extinction_coefficient_calculator.get_instance()
        result = calc.calculate_from_modified_nucleotide_sequence("ACGmTACmGT")
        assert result == pytest.approx(81.119995, abs=0.01)


class TestCalculateFromRnaPolymerNotation:
    """Tests for calculating extinction coefficient from RNA polymer notation."""

    def test_rna_polymer_notation(self, monomer_factory):
        calc = extinction_coefficient_calculator.get_instance()
        result = calc.calculate_from_rna_polymer_notation(
            "P.R(A)P.R(C)P.R(G)P.[mR](A)"
        )
        assert result == pytest.approx(46.200005, abs=0.01)
