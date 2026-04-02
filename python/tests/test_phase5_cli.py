from __future__ import annotations

from helm_notation.cli import main


def test_cli_seq2notation(capsys):
    code = main(["seq2notation", "AUGCUU"])
    captured = capsys.readouterr().out.strip()
    assert code == 0
    assert captured == "R(A)P.R(U)P.R(G)P.R(C)P.R(U)P.R(U)P"


def test_cli_simple_smiles(capsys):
    code = main(["simple-smiles", "PEPTIDE", "A.G"])
    captured = capsys.readouterr().out.strip()
    assert code == 0
    assert "|$" in captured
    assert "_R1" in captured
    assert "_R2" in captured


def test_cli_notation2property(capsys):
    notation = "PEPTIDE1{A.G.G.G.C.C.K.K.K.K}|CHEM1{MCC}$PEPTIDE1,CHEM1,10:R3-1:R1$$$"
    code = main(["notation2property", notation])
    captured = capsys.readouterr().out.strip()
    assert code == 0
    assert captured.startswith("MW=")
    assert ";EXACT=" in captured
    assert ";FORMULA=" in captured


def test_cli_peptide_from_smiles(capsys):
    smiles = "C[C@H](N[*])C(=O)NCC(=O)NCC(=O)NCC(=O)N[C@@H](CS[*])C(=O)N[C@@H](CS[*])C(=O)N[C@@H](CCCCN[*])C(=O)N[C@@H](CCCCN[*])C(=O)N[C@@H](CCCCN[*])C(=O)N[C@@H](CCCCN[*])C([*])=O |$;;;_R1;;;;;;;;;;;;;;;;;;;_R3;;;;;;;_R3;;;;;;;;;;_R3;;;;;;;;;;_R3;;;;;;;;;;_R3;;;;;;;;;;_R3;;_R2;$|"
    code = main(["peptide-from-smiles", smiles])
    captured = capsys.readouterr().out.strip()
    assert code == 0
    assert captured == "A.G.G.G.C.C.K.K.K.K"
