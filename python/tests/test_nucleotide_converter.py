"""Tests for NucleotideConverter."""

import pytest

from helm_notation.tools.nucleotide_converter import NucleotideConverter


@pytest.fixture
def converter():
    """Get NucleotideConverter instance."""
    NucleotideConverter.reset()
    return NucleotideConverter.get_instance()


class TestGetNucleotideSequenceFromSimpleRNANotation:
    """Tests for get_nucleotide_sequence_from_simple_rna_notation."""

    def test_basic_rna_sequence(self, converter):
        """Test converting simple RNA notation to sequence."""
        result = converter.get_nucleotide_sequence_from_simple_rna_notation("R(A)P.R(U)P.R(G)P")
        expected = "AUG"
        print(f"Result: {result}, Expected: {expected}")
        assert result == expected

    def test_rna_with_cytosine(self, converter):
        """Test RNA with cytosine."""
        result = converter.get_nucleotide_sequence_from_simple_rna_notation("R(A)P.R(C)P.R(G)P.R(U)P")
        expected = "ACGU"
        print(f"Result: {result}, Expected: {expected}")
        assert result == expected

    def test_rna_with_all_four(self, converter):
        """Test RNA with all four nucleotides."""
        result = converter.get_nucleotide_sequence_from_simple_rna_notation("R(A)P.R(C)P.R(G)P.R(U)P")
        expected = "ACGU"
        print(f"Result: {result}, Expected: {expected}")
        assert result == expected

    def test_empty_notation(self, converter):
        """Test empty notation."""
        result = converter.get_nucleotide_sequence_from_simple_rna_notation("")
        expected = ""
        print(f"Result: {result}, Expected: {expected}")
        assert result == expected


class TestGetSimpleRNANotationFromSingleNucleotideSequence:
    """Tests for get_simple_rna_notation_from_single_nucleotide_sequence."""

    def test_basic_sequence_to_notation(self, converter):
        """Test converting nucleotide sequence to simple notation."""
        result = converter.get_simple_rna_notation_from_single_nucleotide_sequence("AUG")
        expected = "R(A)P...R(U)P...R(G)"
        print(f"Result: {result}, Expected: {expected}")
        assert result == expected

    def test_longer_sequence(self, converter):
        """Test longer sequence conversion."""
        result = converter.get_simple_rna_notation_from_single_nucleotide_sequence("AUGCUU")
        expected = "R(A)P...R(U)P...R(G)P...R(C)P...R(U)P...R(U)"
        print(f"Result: {result}, Expected: {expected}")
        assert result == expected

    def test_sequence_with_thymine(self, converter):
        """Test DNA sequence with thymine."""
        result = converter.get_simple_rna_notation_from_single_nucleotide_sequence("ATGC")
        expected = "R(A)P...R(T)P...R(G)P...R(C)"
        print(f"Result: {result}, Expected: {expected}")
        assert result == expected


class TestGetComplexNotationFromSingleNucleotideSequence:
    """Tests for get_complex_notation_from_single_nucleotide_sequence."""

    def test_single_sequence_to_complex(self, converter):
        """Test single sequence to complex notation."""
        result = converter.get_complex_notation_from_single_nucleotide_sequence("AUG")
        expected = "RNA1{R(A)P...R(U)P...R(G)}$$$"
        print(f"Result: {result}, Expected: {expected}")
        assert result == expected

    def test_longer_sequence_to_complex(self, converter):
        """Test longer sequence to complex notation."""
        result = converter.get_complex_notation_from_single_nucleotide_sequence("AUGCUU")
        expected = "RNA1{R(A)P...R(U)P...R(G)P...R(C)P...R(U)P...R(U)}$$$"
        print(f"Result: {result}, Expected: {expected}")
        assert result == expected


class TestGetComplexNotation:
    """Tests for get_complex_notation."""

    def test_single_sequence(self, converter):
        """Test single sequence conversion."""
        result = converter.get_complex_notation("AUGCUU")
        expected = "RNA1{R(A)P...R(U)P...R(G)P...R(C)P...R(U)P...R(U)}$$$"
        print(f"Result: {result}, Expected: {expected}")
        assert result == expected

    def test_empty_input(self, converter):
        """Test empty input."""
        result = converter.get_complex_notation("")
        expected = ""
        print(f"Result: {result}, Expected: {expected}")
        assert result == expected

    def test_whitespace_only(self, converter):
        """Test whitespace only input."""
        result = converter.get_complex_notation("   ")
        expected = ""
        print(f"Result: {result}, Expected: {expected}")
        assert result == expected


class TestGetNucleotideSequencesFromComplexNotation:
    """Tests for get_nucleotide_sequences_from_complex_notation."""

    def test_single_rna_node(self, converter):
        """Test single RNA node notation."""
        notation = "RNA1{R(A)P.R(U)P.R(G)P}$$$"
        result = converter.get_nucleotide_sequences_from_complex_notation(notation)
        expected = "AUG"
        print(f"Result: {result}, Expected: {expected}")
        assert result == expected


class TestRoundtrip:
    """Tests for sequence conversion roundtrip."""

    def test_sequence_to_notation_to_sequence(self, converter):
        """Test roundtrip: sequence -> notation -> sequence."""
        original = "AUGC"
        # Convert to simple notation
        simple = converter.get_simple_rna_notation_from_single_nucleotide_sequence(original)
        print(f"Simple notation: {simple}")
        # Convert back to sequence
        result = converter.get_nucleotide_sequence_from_simple_rna_notation(simple)
        expected = "AUGC"
        print(f"Result: {result}, Expected: {expected}")
        assert result == expected

    def test_full_roundtrip(self, converter):
        """Test full roundtrip conversion."""
        original = "AUG"
        # Sequence -> Complex
        complex_notation = converter.get_complex_notation(original)
        expected = "RNA1{R(A)P...R(U)P...R(G)}$$$"
        print(f"Result: {complex_notation}, Expected: {expected}")
        assert complex_notation == expected