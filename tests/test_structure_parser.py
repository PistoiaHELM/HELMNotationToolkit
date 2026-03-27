"""Tests for structure_parser (RDKit chemistry layer)."""

from __future__ import annotations

import pytest

from helm_notation.tools import structure_parser
from helm_notation.exceptions import StructureException
from helm_notation.model.molecule_info import MoleculeInfo


class TestGetMolecule:
    """Test SMILES/molfile parsing."""

    def test_simple_smiles(self):
        mol = structure_parser.get_molecule("CCO")
        assert mol is not None
        assert mol.GetNumAtoms() == 3  # C, C, O (no explicit H)

    def test_invalid_smiles_returns_none(self):
        mol = structure_parser.get_molecule("not_a_smiles!!!")
        assert mol is None

    def test_none_input_returns_none(self):
        mol = structure_parser.get_molecule(None)
        assert mol is None

    def test_cxsmiles_with_rgroups(self):
        """Parse ChemAxon extended SMILES with R-group labels."""
        cxsmi = "[*]OC1CC(O[*])C(CO[*])O1 |$_R1;;;;;_R2;;;;_R3;;$|"
        mol = structure_parser.get_molecule(cxsmi)
        assert mol is not None
        # Should have 3 dummy atoms with isotopes 1, 2, 3
        dummies = [
            a for a in mol.GetAtoms() if a.GetAtomicNum() == 0
        ]
        assert len(dummies) == 3
        isotopes = sorted(a.GetIsotope() for a in dummies)
        assert isotopes == [1, 2, 3]

    def test_cxsmiles_single_rgroup(self):
        cxsmi = "[*]O |$_R1;$|"
        mol = structure_parser.get_molecule(cxsmi)
        assert mol is not None
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert len(dummies) == 1
        assert dummies[0].GetIsotope() == 1


class TestValidateSmiles:
    def test_valid(self):
        assert structure_parser.validate_smiles("CCO") is True

    def test_invalid(self):
        assert structure_parser.validate_smiles("not_smiles!!!") is False


class TestMoleculeInfo:
    def test_water(self):
        mi = structure_parser.get_molecule_info("O")
        assert mi.molecular_formula == "H2O"
        assert abs(mi.molecular_weight - 18.02) < 0.1
        assert abs(mi.exact_mass - 18.01) < 0.1

    def test_ethanol(self):
        mi = structure_parser.get_molecule_info("CCO")
        assert "C2" in mi.molecular_formula
        assert abs(mi.molecular_weight - 46.07) < 0.1

    def test_invalid_raises(self):
        with pytest.raises(StructureException):
            structure_parser.get_molecule_info("invalid!!!")


class TestRgroupOperations:
    def _make_mol_with_rgroup(self, smiles_with_rgroup):
        """Helper to make a molecule with R-group from CXSMILES."""
        return structure_parser.get_molecule(smiles_with_rgroup)

    def test_find_rgroup_atom(self):
        mol = self._make_mol_with_rgroup("[*]O |$_R1;$|")
        idx = structure_parser.find_rgroup_atom(mol, 1)
        atom = mol.GetAtomWithIdx(idx)
        assert atom.GetAtomicNum() == 0
        assert atom.GetIsotope() == 1

    def test_find_nonexistent_rgroup_raises(self):
        mol = self._make_mol_with_rgroup("[*]O |$_R1;$|")
        with pytest.raises(StructureException, match="R99"):
            structure_parser.find_rgroup_atom(mol, 99)

    def test_remove_rgroup(self):
        from rdkit import Chem

        mol = self._make_mol_with_rgroup("[*]OCC |$_R1;;;$|")
        rw = Chem.RWMol(mol)
        neighbor_idx = structure_parser.remove_rgroup(rw, 1)
        # After removing the dummy, we should have O, C, C
        assert rw.GetNumAtoms() == 3
        # The neighbor should be an oxygen
        assert rw.GetAtomWithIdx(neighbor_idx).GetAtomicNum() == 8


class TestMerge:
    def test_merge_two_molecules(self):
        from rdkit import Chem

        # mol1: [1*]OC  (R1-O-C)
        # mol2: [2*]NC  (R2-N-C) -- use R2 but we merge at R2
        mol1 = structure_parser.get_molecule("[*]OC |$_R1;;$|")
        mol2 = structure_parser.get_molecule("[*]NC |$_R2;;$|")
        assert mol1 is not None and mol2 is not None

        rw1 = Chem.RWMol(mol1)
        rw2 = Chem.RWMol(mol2)

        r1_idx = structure_parser.find_rgroup_atom(rw1, 1)
        r2_idx = structure_parser.find_rgroup_atom(rw2, 2)

        merged = structure_parser.merge_ignore_stereo(rw1, r1_idx, rw2, r2_idx)
        # Should have O, C, N, C = 4 atoms
        assert merged.GetNumAtoms() == 4


class TestProcessMoleculeInfo:
    def test_add_and_subtract(self):
        chunk = MoleculeInfo()
        chunk.molecular_weight = 100.0
        chunk.exact_mass = 99.5
        chunk.molecular_formula = "C4H8O2"

        cap = MoleculeInfo()
        cap.molecular_weight = 18.0
        cap.exact_mass = 18.0
        cap.molecular_formula = "H2O"

        result = structure_parser.process_molecule_info([chunk], [cap])
        assert abs(result.molecular_weight - 82.0) < 0.01
        assert abs(result.exact_mass - 81.5) < 0.01
        # Formula should reflect C4H6O1
        assert "C4" in result.molecular_formula


class TestUniqueSmiles:
    def test_canonical(self):
        s1 = structure_parser.get_unique_smiles("OCC")
        s2 = structure_parser.get_unique_smiles("CCO")
        assert s1 == s2

    def test_invalid_raises(self):
        with pytest.raises(StructureException):
            structure_parser.get_unique_smiles("invalid!!!")


class TestGetSimpleExtendedSmiles:
    def test_strips_extra_extensions(self):
        # Format: SMILES |extra,$atoms$,extra|
        ext_smi = "[*]OC |r,$_R1;;$,c:0|"
        result = structure_parser.get_simple_extended_smiles(ext_smi)
        assert "_R1" in result
        assert "c:0" not in result
        assert "r," not in result


class TestDecomposeMixture:
    def test_single_component(self):
        # Need proper extended SMILES format with extra components
        smi = "[*]O |,$_R1;$,|"
        parts = structure_parser.decompose_mixture(smi)
        assert len(parts) >= 1


class TestGetRgroupsFromExtendedSmiles:
    def test_extracts_rgroups(self):
        smi = "[*]OC1CC(O[*])C(CO[*])O1 |$_R1;;;;;_R2;;;;_R3;;$|"
        rgroups = structure_parser.get_rgroups_from_extended_smiles(smi)
        assert "R1" in rgroups
        assert "R2" in rgroups
        assert "R3" in rgroups
