from helm_notation.chem import ChemEngine


def test_validate_smiles():
    assert ChemEngine.validate_smiles("CCO")
    assert not ChemEngine.validate_smiles("not-a-smiles")


def test_unique_smiles_stable():
    s1 = ChemEngine.get_unique_smiles("C(C)O")
    s2 = ChemEngine.get_unique_smiles("CCO")
    assert s1 == s2


def test_molecule_info_ethanol():
    info = ChemEngine.get_molecule_info("CCO")
    assert info.molecular_formula == "C2H6O"
    assert info.molecular_weight > 46.0
    assert info.exact_mass > 46.0


def test_decompose_extended_smiles_single_component():
    ext = "CCO|$;;;$|"
    parts = ChemEngine.decompose_mixture(ext)
    assert len(parts) == 1
    assert parts[0] == "CCO|$;;$|"
