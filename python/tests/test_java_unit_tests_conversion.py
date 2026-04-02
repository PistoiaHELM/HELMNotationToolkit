from __future__ import annotations

from pathlib import Path

import pytest

from helm_notation import ComplexNotationParser, SimpleNotationParser
from helm_notation.cli import _sequence_to_notation
from helm_notation.monomer_factory import MonomerFactory
from helm_notation.tools.molfile_encoder import decode, encode


JAVA_TEST_DIR = Path("test/org/helm/notation/tools")
JAVA_TEST_CLASSES = sorted(p.stem for p in JAVA_TEST_DIR.glob("*Test.java"))

CONVERTED_CLASSES = {
    "CalculatorTest",
    "ComplexNotationParserTest",
    "EncoderTest",
    "MonomerDBTest",
    "NucleotideTest",
    "SimpleNotationGroupIteratorTest",
    "SimpleNotationParserTest",
    "xHelmNotationExporterTest",
    "xHelmNotationParserTest",
}


def test_all_java_test_classes_are_represented_in_pytest_conversion():
    assert set(JAVA_TEST_CLASSES) == CONVERTED_CLASSES


# EncoderTest

def test_encoder_testencode_roundtrip():
    source = "HEL\nM"
    assert decode(encode(source)) == source


# MonomerDBTest

def test_monomerdb_initialize_monomer_db_smoke():
    factory = MonomerFactory.get_instance()
    db = factory.get_monomer_db()
    assert "PEPTIDE" in db
    assert "RNA" in db
    assert "CHEM" in db


# SimpleNotationParserTest (implemented subset)

def test_simple_notation_parser_get_smiles_subset():
    out = SimpleNotationParser.get_simple_polymer_smiles("A.G.C", "PEPTIDE")
    assert isinstance(out, str)
    assert "|$" in out


# ComplexNotationParserTest (implemented subset)

def test_complex_notation_parser_canonicalization_subset():
    notation = "PEPTIDE1{K.A.A.G.K}$PEPTIDE1,PEPTIDE1,1:R1-5:R2$$$"
    out = ComplexNotationParser.get_canonical_notation(notation)
    assert out in {
        "PEPTIDE1{A.A.G.K.K}$PEPTIDE1,PEPTIDE1,1:R1-5:R2$$$",
        "PEPTIDE1{A.A.G.K.K}$PEPTIDE1,PEPTIDE1,5:R2-1:R1$$$",
    }


# NucleotideTest (implemented subset)

def test_nucleotide_get_notation_subset():
    assert _sequence_to_notation("AUG") == "R(A)P.R(U)P.R(G)P"


def test_nucleotide_is_in_opposite_direction():
    from helm_notation.nucleotide_sequence_parser import is_in_opposite_direction

    assert is_in_opposite_direction("AUGTAU", "AUACAU") is False
    assert is_in_opposite_direction("ATGCAAA", "UUUGCAU") is False


def test_nucleotide_get_sirna_notation():
    from helm_notation.nucleotide_sequence_parser import (
        RNA_DESIGN_TUSCHL_19_PLUS_2,
        get_sirna_notation,
    )

    sense = "CGAAAUGUUCAUACUGUUGdTdT"
    antisense = "UUACAAUUUGGACUUUCCGdTdT"

    assert (
        get_sirna_notation(sense, antisense)
        == "RNA1{R(C)P.R(G)P.R(A)P.R(A)P.R(A)P.R(U)P.R(G)P.R(U)P.R(U)P.R(C)P.R(A)P.R(U)P.R(A)P.R(C)P.R(U)P.R(G)P.R(U)P.R(U)P.R(G)P.[dR](T)P.[dR](T)}|RNA2{R(U)P.R(U)P.R(A)P.R(C)P.R(A)P.R(A)P.R(U)P.R(U)P.R(U)P.R(G)P.R(G)P.R(A)P.R(C)P.R(U)P.R(U)P.R(U)P.R(C)P.R(C)P.R(G)P.[dR](T)P.[dR](T)}$$RNA1,RNA2,5:pair-50:pair|RNA1,RNA2,8:pair-47:pair|RNA1,RNA2,11:pair-44:pair|RNA1,RNA2,14:pair-41:pair$RNA1{ss}|RNA2{as}$"
    )

    assert (
        get_sirna_notation(sense, antisense, RNA_DESIGN_TUSCHL_19_PLUS_2)
        == "RNA1{R(C)P.R(G)P.R(A)P.R(A)P.R(A)P.R(U)P.R(G)P.R(U)P.R(U)P.R(C)P.R(A)P.R(U)P.R(A)P.R(C)P.R(U)P.R(G)P.R(U)P.R(U)P.R(G)P.[dR](T)P.[dR](T)}|RNA2{R(U)P.R(U)P.R(A)P.R(C)P.R(A)P.R(A)P.R(U)P.R(U)P.R(U)P.R(G)P.R(G)P.R(A)P.R(C)P.R(U)P.R(U)P.R(U)P.R(C)P.R(C)P.R(G)P.[dR](T)P.[dR](T)}$$RNA1,RNA2,2:pair-56:pair|RNA1,RNA2,5:pair-53:pair|RNA1,RNA2,11:pair-47:pair|RNA1,RNA2,14:pair-44:pair|RNA1,RNA2,20:pair-38:pair|RNA1,RNA2,23:pair-35:pair|RNA1,RNA2,29:pair-29:pair|RNA1,RNA2,32:pair-26:pair|RNA1,RNA2,38:pair-20:pair|RNA1,RNA2,44:pair-14:pair|RNA1,RNA2,47:pair-11:pair|RNA1,RNA2,50:pair-8:pair$RNA1{ss}|RNA2{as}$"
    )


def test_nucleotide_get_complement_sequence():
    from helm_notation.nucleotide_sequence_parser import (
        get_normal_complement_sequence,
        get_reverse_complement_sequence,
    )

    assert get_normal_complement_sequence("agcuagggu") == "5'-ACCCUAGCU-3'"
    assert get_reverse_complement_sequence("agcuagggu") == "3'-UCGAUCCCA-5'"


# CalculatorTest

def test_calculator_test_class():
    from helm_notation.tools.extinction_coefficient_calculator import (
        ExtinctionCoefficientCalculator,
    )

    calculator = ExtinctionCoefficientCalculator.get_instance()
    assert calculator.calculate_from_nucleotide_sequence("ACGTACGT") == pytest.approx(81.119995, abs=2e-5)
    assert calculator.calculate_from_modified_nucleotide_sequence("ACGmTACmGT") == pytest.approx(81.119995, abs=2e-5)
    assert calculator.calculate_from_amino_acid_sequence("AGGYEEEEEEEEEEEEEEEEEEE") == pytest.approx(1490.0, abs=1e-15)


# SimpleNotationGroupIteratorTest

def test_simple_notation_group_iterator():
    from helm_notation.simple_notation_group_iterator import SimpleNotationGroupIterator

    def collect(notation: str) -> list[str]:
        iterator = SimpleNotationGroupIterator(notation)
        groups: list[str] = []
        while iterator.has_next_group():
            groups.append(iterator.next_group())
        return groups

    assert collect("A.K.A") == ["A", "K", "A"]
    assert collect("P.R(A)[sP].RP.R(G)P.[LR]([5meC])") == [
        "P",
        "R(A)[sP]",
        "RP",
        "R(G)P",
        "[LR]([5meC])",
    ]
    assert collect("SMPEG2") == ["SMPEG2"]
    assert collect("G.G.K.A.A.[[SeH]C[C@H](N[*])C([*])=O |$;;;;_R1;;_R2;$|].[meC]") == [
        "G",
        "G",
        "K",
        "A",
        "A",
        "[[SeH]C[C@H](N[*])C([*])=O |$;;;;_R1;;_R2;$|]",
        "[meC]",
    ]
    assert collect("C[C@H](N[*])C([*])=O |$;;;;_R1;;_R2;$|") == [
        "C[C@H](N[*])C([*])=O |$;;;;_R1;;_R2;$|"
    ]
    assert collect(
        "P.[[*]OC[C@@H]1CN([*])C[C@H]([*])O1 |$_R1;;;;;;_R2;;;_R3;$|](A)[sP].RP.R(G)P.[LR]([5meC])"
    ) == [
        "P",
        "[[*]OC[C@@H]1CN([*])C[C@H]([*])O1 |$_R1;;;;;;_R2;;;_R3;$|](A)[sP]",
        "RP",
        "R(G)P",
        "[LR]([5meC])",
    ]


# xHelmNotationParserTest / xHelmNotationExporterTest

def test_xhelm_parser_port():
    from helm_notation.tools.xhelm_notation_parser import (
        get_complex_notation_string,
        get_monomer_store,
    )

    fixture = Path("test/org/helm/notation/tools/resources/simple.xhelm")
    notation = get_complex_notation_string(fixture)
    assert notation == "PEPTIDE1{G.K.A.[A_copy]}$$$$"

    store = get_monomer_store(fixture)
    peptide_db = store.get_monomer_db().get("PEPTIDE", {})
    assert "A" in peptide_db
    assert "A_copy" in peptide_db
    assert "G" in peptide_db
    assert "K" in peptide_db


def test_xhelm_exporter_port():
    from helm_notation.tools.xhelm_notation_exporter import write_xhelm
    from helm_notation.tools.xhelm_notation_parser import (
        get_complex_notation_string,
        get_monomer_store,
    )

    fixture = Path("test/org/helm/notation/tools/resources/simple.xhelm")
    notation = get_complex_notation_string(fixture)
    store = get_monomer_store(fixture)

    output = write_xhelm(notation, store)
    assert "<Xhelm>" in output
    assert f"<HelmNotation>{notation}</HelmNotation>" in output
    assert output.count("<Monomer>") >= 4
