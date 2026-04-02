"""Phase 2 tests: data models, MonomerFactory, NucleotideFactory."""

import pytest

from helm_notation.model import (
    Attachment,
    ComplexPolymer,
    Monomer,
    MoleculeInfo,
    Nucleotide,
    PolymerEdge,
    PolymerNode,
    RNAPolymerNode,
    get_polymer_type,
    PEPTIDE_POLYMER_TYPE,
    NUCLIEC_ACID_POLYMER_TYPE,
    CHEMICAL_POLYMER_TYPE,
)


# ── Model smoke tests ────────────────────────────────────────────

class TestModels:
    def test_attachment_str(self):
        a = Attachment(label="R1", cap_group_name="H")
        assert str(a) == "R1: H"

    def test_monomer_attachment_lookup(self):
        m = Monomer(
            alternate_id="Ala",
            polymer_type=PEPTIDE_POLYMER_TYPE,
            monomer_type="Backbone",
            attachment_list=[
                Attachment(alternate_id="R1-H", label="R1", cap_group_name="H"),
                Attachment(alternate_id="R2-OH", label="R2", cap_group_name="OH"),
            ],
        )
        assert m.get_attachment("R1") is not None
        assert m.get_attachment("R1").cap_group_name == "H"
        assert m.get_attachment("R3") is None

    def test_nucleotide_is_modified(self):
        assert Nucleotide(notation="R(A)P").is_modified() is False
        assert Nucleotide(notation="[mR](A)P").is_modified() is True
        assert Nucleotide(notation="R(X)P").is_modified() is True
        assert Nucleotide(notation="R(A)").is_modified() is True  # no phosphate

    def test_polymer_node_type(self):
        assert get_polymer_type("RNA1") == "RNA"
        assert get_polymer_type("PEPTIDE2") == "PEPTIDE"
        assert get_polymer_type("CHEM3") == "CHEM"

    def test_rna_polymer_node(self):
        n = RNAPolymerNode(id="RNA1", label="R(A)P")
        assert n.type == NUCLIEC_ACID_POLYMER_TYPE

    def test_polymer_edge_properties(self):
        e = PolymerEdge(
            source_node="RNA1",
            target_node="CHEM1",
            connection="5:R2-1:R1",
        )
        assert e.source_connection == "5:R2"
        assert e.target_connection == "1:R1"

    def test_complex_polymer_default(self):
        cp = ComplexPolymer()
        assert cp.polymer_node_list == []
        assert cp.polymer_edge_list == []


# ── MonomerFactory ────────────────────────────────────────────────

class TestMonomerFactory:
    def test_loads_monomers(self):
        from helm_notation.monomer_factory import MonomerFactory

        MonomerFactory.finalize_monomer_cache()
        factory = MonomerFactory.get_instance()

        db = factory.get_monomer_db()
        assert "PEPTIDE" in db
        assert "RNA" in db
        assert "CHEM" in db

        # Check a well-known monomer
        assert "A" in db["RNA"]
        assert "Ala" in db["PEPTIDE"] or "A" in db["PEPTIDE"]

    def test_monomer_store_view(self):
        from helm_notation.monomer_factory import MonomerFactory

        factory = MonomerFactory.get_instance()
        store = factory.get_monomer_store()
        assert store.has_monomer("PEPTIDE", "G") or store.has_monomer("RNA", "A")

    def test_singleton(self):
        from helm_notation.monomer_factory import MonomerFactory

        f1 = MonomerFactory.get_instance()
        f2 = MonomerFactory.get_instance()
        assert f1 is f2

    def test_finalize_resets(self):
        from helm_notation.monomer_factory import MonomerFactory

        MonomerFactory.get_instance()
        MonomerFactory.finalize_monomer_cache()
        # After finalize, next get_instance creates a new one
        f = MonomerFactory.get_instance()
        assert f is not None


# ── NucleotideFactory ─────────────────────────────────────────────

class TestNucleotideFactory:
    def test_loads_templates(self):
        from helm_notation.nucleotide_factory import NucleotideFactory

        NucleotideFactory.reset()
        nf = NucleotideFactory.get_instance()
        templates = nf.get_nucleotide_templates()
        assert "HELM Notation" in templates
        helm = templates["HELM Notation"]
        assert helm["A"] == "R(A)P"
        assert helm["U"] == "R(U)P"
        assert helm["G"] == "R(G)P"
        assert helm["C"] == "R(C)P"
        assert helm["T"] == "R(T)P"
        assert helm["dA"] == "[dR](A)P"

    def test_reverse_map(self):
        from helm_notation.nucleotide_factory import NucleotideFactory

        nf = NucleotideFactory.get_instance()
        assert nf.get_symbol_for_notation("R(A)P") == "A"

    def test_notation_for_symbol(self):
        from helm_notation.nucleotide_factory import NucleotideFactory

        nf = NucleotideFactory.get_instance()
        assert nf.get_notation_for_symbol("A") == "R(A)P"
        assert nf.get_notation_for_symbol("NONEXISTENT") is None


# ── MonomerParser XML round-trip ──────────────────────────────────

class TestMonomerParser:
    def test_round_trip(self):
        import xml.etree.ElementTree as ET
        from helm_notation.tools.monomer_parser import monomer_to_element, parse_monomer

        m = Monomer(
            alternate_id="TestM",
            can_smiles="CCO",
            monomer_type="Backbone",
            polymer_type="PEPTIDE",
            natural_analog="X",
            name="test monomer",
            attachment_list=[
                Attachment(alternate_id="R1-H", label="R1", cap_group_name="H",
                           cap_group_smiles="[*][H] |$_R1;$|"),
            ],
        )
        el = monomer_to_element(m)
        m2 = parse_monomer(el)
        assert m2.alternate_id == "TestM"
        assert m2.can_smiles == "CCO"
        assert len(m2.attachment_list) == 1
        assert m2.attachment_list[0].label == "R1"


# ── MolfileEncoder ────────────────────────────────────────────────

class TestMolfileEncoder:
    def test_round_trip(self):
        from helm_notation.tools.molfile_encoder import encode, decode

        original = "Hello HELM World\nLine 2"
        encoded = encode(original)
        assert isinstance(encoded, str)
        decoded = decode(encoded)
        assert decoded == original

    def test_none(self):
        from helm_notation.tools.molfile_encoder import encode, decode

        assert encode(None) is None
        assert decode(None) is None
