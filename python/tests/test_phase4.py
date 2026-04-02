from __future__ import annotations

import json
from pathlib import Path

import pytest
from rdkit import Chem

from helm_notation import PeptideStructureParser, SimpleNotationParser, ComplexNotationParser
from helm_notation.exceptions import StructureException
from helm_notation.simple_notation_parser import _find_r_atom, cyclize_on_r_atoms


def _baseline_case(case_id: str) -> dict:
    data = json.loads(Path("migration/fixtures/java_baseline.json").read_text("utf-8"))
    for case in data["cases"]:
        if case["id"] == case_id:
            return case
    raise KeyError(case_id)


def test_smiles_to_notation_linear_peptide_baseline():
    case = _baseline_case("simple_peptide_smiles_1")
    parser = PeptideStructureParser.get_instance()

    notation = parser.smiles2notation(case["output"])

    assert notation == "A.G.G.G.C.C.K.K.K.K"


def test_roundtrip_simple_peptide_smiles_to_notation():
    parser = PeptideStructureParser.get_instance()

    input_notation = "A.G.C.K"
    smi = SimpleNotationParser.get_simple_polymer_smiles(input_notation, "PEPTIDE")
    out_notation = parser.smiles2notation(smi)

    assert out_notation == input_notation


def test_mixture_rejected():
    parser = PeptideStructureParser.get_instance()
    with pytest.raises(StructureException, match="mixture"):
        parser.smiles2notation("CCO.CCN")


def test_cyclic_peptide_backbone_to_complex_notation():
    parser = PeptideStructureParser.get_instance()
    structure = SimpleNotationParser.get_simple_polymer_structure("K.A.A.G.K", "PEPTIDE")
    mol = structure.molecule

    r1 = _find_r_atom(mol, monomer_position=1, r_label="R1")
    r2 = _find_r_atom(mol, monomer_position=5, r_label="R2")
    cyc = cyclize_on_r_atoms(mol, r1, r2)

    out = parser.molecule2notation(cyc)
    expected_1 = "PEPTIDE1{A.A.G.K.K}$PEPTIDE1,PEPTIDE1,1:R1-5:R2$$$"
    expected_2 = "PEPTIDE1{A.A.G.K.K}$PEPTIDE1,PEPTIDE1,5:R2-1:R1$$$"
    assert out in (expected_1, expected_2)


def test_disulfide_crosslink_to_complex_notation():
    parser = PeptideStructureParser.get_instance()
    structure = SimpleNotationParser.get_simple_polymer_structure("C.G.C", "PEPTIDE")
    mol = structure.molecule

    r3_1 = _find_r_atom(mol, monomer_position=1, r_label="R3")
    r3_3 = _find_r_atom(mol, monomer_position=3, r_label="R3")
    linked = cyclize_on_r_atoms(mol, r3_1, r3_3)

    out = parser.molecule2notation(linked)
    expected = "PEPTIDE1{C.G.C}$PEPTIDE1,PEPTIDE1,1:R3-3:R3$$$"
    assert out == expected
