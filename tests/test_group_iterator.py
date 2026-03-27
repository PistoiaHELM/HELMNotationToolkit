"""Tests for SimpleNotationGroupIterator."""

from __future__ import annotations

from helm_notation.tools.simple_notation_group_iterator import SimpleNotationGroupIterator


class TestSimplePeptide:
    """Tests for iterating over simple peptide notation."""

    def test_simple_peptide(self):
        iterator = SimpleNotationGroupIterator("A.K.A")
        groups = list(iterator)
        assert len(groups) == 3, "there should be three groups!"
        assert groups[0] == "A"
        assert groups[1] == "K"
        assert groups[2] == "A"


class TestSimpleRNA:
    """Tests for iterating over simple RNA notation."""

    def test_simple_rna(self):
        iterator = SimpleNotationGroupIterator("P.R(A)[sP].RP.R(G)P.[LR]([5meC])")
        groups = list(iterator)
        assert len(groups) == 5, "there should be five groups!"
        assert groups[0] == "P"
        assert groups[1] == "R(A)[sP]"
        assert groups[2] == "RP"
        assert groups[3] == "R(G)P"
        assert groups[4] == "[LR]([5meC])"


class TestSimpleChem:
    """Tests for iterating over simple chem notation."""

    def test_simple_chem(self):
        iterator = SimpleNotationGroupIterator("SMPEG2")
        groups = list(iterator)
        assert len(groups) == 1, "there should be one group!"
        assert groups[0] == "SMPEG2"


class TestInlineSmilesPeptide:
    """Tests for iterating over peptide notation with inline SMILES."""

    def test_inline_smiles_peptide(self):
        iterator = SimpleNotationGroupIterator(
            "G.G.K.A.A.[[SeH]C[C@H](N[*])C([*])=O |$;;;;_R1;;_R2;$|].[meC]"
        )
        groups = list(iterator)
        assert len(groups) == 7, "there should be 7 groups!"
        assert groups[0] == "G"
        assert groups[1] == "G"
        assert groups[2] == "K"
        assert groups[3] == "A"
        assert groups[4] == "A"
        assert groups[5] == "[[SeH]C[C@H](N[*])C([*])=O |$;;;;_R1;;_R2;$|]"
        assert groups[6] == "[meC]"


class TestInlineSmilesChem:
    """Tests for iterating over chem notation with inline SMILES."""

    def test_inline_smiles_chem(self):
        iterator = SimpleNotationGroupIterator(
            "C[C@H](N[*])C([*])=O |$;;;;_R1;;_R2;$|"
        )
        groups = list(iterator)
        assert len(groups) == 1, "there should be one group!"
        assert groups[0] == "C[C@H](N[*])C([*])=O |$;;;;_R1;;_R2;$|"


class TestInlineSmilesRNA:
    """Tests for iterating over RNA notation with inline SMILES."""

    def test_inline_smiles_rna(self):
        iterator = SimpleNotationGroupIterator(
            "P.[[*]OC[C@@H]1CN([*])C[C@H]([*])O1 |$_R1;;;;;;_R2;;;_R3;$|](A)[sP].RP.R(G)P.[LR]([5meC])"
        )
        groups = list(iterator)
        assert len(groups) == 5, "there should be five groups!"
        assert groups[0] == "P"
        assert groups[1] == (
            "[[*]OC[C@@H]1CN([*])C[C@H]([*])O1 |$_R1;;;;;;_R2;;;_R3;$|](A)[sP]"
        )
        assert groups[2] == "RP"
        assert groups[3] == "R(G)P"
        assert groups[4] == "[LR]([5meC])"
