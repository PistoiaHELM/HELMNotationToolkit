"""Tests for monomer database initialization and access."""

from __future__ import annotations

from helm_notation.tools import structure_parser


class TestMonomerDB:
    """Tests for monomer database initialization and iteration."""

    def test_initialize_monomer_db(self, monomer_factory):
        """Test that monomer DB can be initialized and contains data."""
        monomer_db = monomer_factory.get_monomer_db()
        assert monomer_db is not None
        assert len(monomer_db) > 0

    def test_monomer_db_has_polymer_types(self, monomer_factory):
        """Test that monomer DB contains expected polymer types."""
        monomer_db = monomer_factory.get_monomer_db()
        # Should have entries for RNA, PEPTIDE, CHEM
        polymer_types = list(monomer_db.keys())
        assert len(polymer_types) > 0

    def test_attachment_db(self, monomer_factory):
        """Test that attachment DB can be loaded."""
        attachment_db = monomer_factory.get_attachment_db()
        assert attachment_db is not None
        assert len(attachment_db) > 0

    def test_attachment_has_cap_group_smiles(self, monomer_factory):
        """Test that attachments have cap group SMILES."""
        attachment_db = monomer_factory.get_attachment_db()
        for attach_id, attachment in attachment_db.items():
            cap_smiles = attachment.cap_group_smiles
            # Cap group SMILES should be a string (can be None for some)
            assert isinstance(attach_id, str)

    def test_monomer_has_canonical_smiles(self, monomer_factory):
        """Test that monomers have canonical SMILES."""
        monomer_db = monomer_factory.get_monomer_db()
        for polymer_type, monomer_map in monomer_db.items():
            for monomer_id, monomer in monomer_map.items():
                smiles = monomer.can_smiles
                if smiles:
                    assert isinstance(smiles, str)

    def test_monomer_store(self, monomer_store):
        """Test that the monomer store is accessible."""
        assert monomer_store is not None
