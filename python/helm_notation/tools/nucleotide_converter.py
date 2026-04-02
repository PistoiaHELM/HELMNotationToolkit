"""Port of org.helm.notation.tools.NucleotideConverter.

Provides conversion between nucleotide sequences and HELM notations.
Note: Some functionality is provided by nucleotide_sequence_parser.py.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional

from helm_notation import constants
from helm_notation.exceptions import NotationException
from helm_notation.monomer_factory import MonomerFactory
from helm_notation.model import ID_P, NUCLIEC_ACID_POLYMER_TYPE, PolymerNode
from helm_notation.nucleotide_factory import NucleotideFactory
from helm_notation.nucleotide_sequence_parser import get_notation as seq_to_notation
from helm_notation.nucleotide_sequence_parser import get_normal_list, remove_last_p

logger = logging.getLogger(__name__)


class NucleotideConverter:
    """Converter for nucleotide sequences and HELM notations."""

    _instance: Optional["NucleotideConverter"] = None

    def __init__(self) -> None:
        pass

    @staticmethod
    def get_instance() -> "NucleotideConverter":
        """Get singleton instance, initializing factories if needed."""
        if NucleotideConverter._instance is None:
            MonomerFactory.get_instance()
            NucleotideFactory.get_instance()
            NucleotideConverter._instance = NucleotideConverter()
        return NucleotideConverter._instance

    @staticmethod
    def reset() -> None:
        """Reset singleton instance."""
        NucleotideConverter._instance = None

    def get_nucleotide_sequences_from_complex_notation(
        self, complex_notation: str
    ) -> str:
        """Get nucleotide sequences from a complex HELM notation.
        
        Args:
            complex_notation: HELM notation containing RNA polymers
            
        Returns:
            Nucleotide sequences separated by spaces
            
        Raises:
            NotationException: If notation contains non-nucleic acid polymers
        """
        polymer_nodes = self._get_polymer_node_list(complex_notation)

        parts: List[str] = []
        for node in polymer_nodes:
            polymer_type = node.type
            if polymer_type != NUCLIEC_ACID_POLYMER_TYPE:
                raise NotationException(
                    "Input complex notation contains non-nucleic acid polymer"
                )
            simple_notation = node.label
            notation = self.get_nucleotide_sequence_from_simple_rna_notation(
                simple_notation
            )
            parts.append(notation)

        return " ".join(parts)

    def _get_polymer_node_list(self, extended_notation: str) -> List[PolymerNode]:
        """Parse complex notation to get polymer nodes."""
        from helm_notation.complex_notation_parser import _parse_extended_notation

        nodes, _, _, _, _ = _parse_extended_notation(extended_notation)
        result: List[PolymerNode] = []
        for node in nodes:
            if node.polymer_type == NUCLIEC_ACID_POLYMER_TYPE:
                pn = PolymerNode(id=node.node_id, label=node.label)
                result.append(pn)
        return result

    def get_nucleotide_sequence_from_simple_rna_notation(
        self, simple_rna_notation: str
    ) -> str:
        """Get nucleotide sequence from simple RNA notation.
        
        Args:
            simple_rna_notation: Simple RNA notation like R(A)P.R(U)P.R(G)P
            
        Returns:
            Nucleotide sequence (e.g., "AUG")
            
        Raises:
            NotationException: If unknown nucleotide is found
        """
        nuc_list = self._parse_rna_notation(simple_rna_notation)

        reverse_nuc_map = NucleotideFactory.get_instance().get_reverse_nucleotide_map()

        parts: List[str] = []
        count = 0
        n = len(nuc_list)

        for nuc in nuc_list:
            nucleotide = nuc.notation
            nucleoside = nuc.natural_analog
            linker = nuc.linker_notation if hasattr(nuc, 'linker_notation') else ""

            # First nucleotide can skip nucleoside
            if count == 0 and not nucleoside:
                parts.append(nuc.symbol)
                count += 1
                continue

            # Last nucleotide can skip linker
            if count == n - 1 and not linker:
                nucleotide = nucleotide + ID_P

            if nucleotide in reverse_nuc_map:
                parts.append(reverse_nuc_map[nucleotide])
            else:
                raise NotationException(
                    f"Unknown nucleotide found: {nucleotide} - missing nucleotide template"
                )

            count += 1

        return "".join(parts)

    def _parse_rna_notation(self, notation: str) -> List["_NucToken"]:
        """Parse simple RNA notation into nucleotide tokens.
        
        Args:
            notation: Simple RNA notation like R(A)P.R(U)P.R(G)P
            
        Returns:
            List of nucleotide tokens
        """
        from dataclasses import dataclass
        
        @dataclass
        class _NucToken:
            symbol: str
            notation: str
            natural_analog: str
            linker_notation: str
        
        tokens = []
        parts = notation.split(".")
        
        for i, part in enumerate(parts):
            part = part.strip()
            if not part:
                continue
            
            # Match pattern like R(A)P, R(U)P, R(G)P, [mA]P, etc.
            import re
            match = re.match(r"(\[?[^\]]+\]?)\((\w+)\)(P?)", part)
            if match:
                monomer_id = match.group(1)
                natural_analog = match.group(2)
                has_p = match.group(3)
                
                # Look up the notation from factory
                nf = NucleotideFactory.get_instance()
                templates = nf.get_nucleotide_templates().get(constants.NOTATION_SOURCE, {})
                
                # Find the notation for this symbol
                notation_val = templates.get(monomer_id, "")
                if not notation_val:
                    # Try to build it
                    notation_val = f"R({natural_analog})P" if has_p else f"R({natural_analog})"
                
                tokens.append(_NucToken(
                    symbol=monomer_id,
                    notation=notation_val,
                    natural_analog=natural_analog,
                    linker_notation="P" if has_p else ""
                ))
            else:
                raise NotationException(f"Invalid nucleotide notation: {part}")
        
        return tokens

    def get_complex_notation(self, nucleotide_sequences: str) -> str:
        """Convert nucleotide sequences to HELM complex notation.
        
        Args:
            nucleotide_sequences: Space-separated nucleotide sequences (e.g., "AUGCUU AUGC")
            
        Returns:
            HELM complex notation
            
        Raises:
            NotationException: If conversion fails
        """
        if not nucleotide_sequences or not nucleotide_sequences.strip():
            return ""

        sequences = nucleotide_sequences.strip().split()
        
        if len(sequences) == 1:
            result = self.get_complex_notation_from_single_nucleotide_sequence(
                sequences[0]
            )
        else:
            result = self.get_complex_notation_from_single_nucleotide_sequence(
                sequences[0]
            )
            for seq in sequences[1:]:
                notation = self.get_complex_notation_from_single_nucleotide_sequence(seq)
                result = self._combine_complex_notations(result, notation)

        result = self._hybridize(result)
        return result

    def get_complex_notation_from_single_nucleotide_sequence(
        self, single_oligo_sequence: str
    ) -> str:
        """Convert single oligonucleotide sequence to complex notation.
        
        Args:
            single_oligo_sequence: Single sequence like "AUGCUU"
            
        Returns:
            Complex notation
        """
        simple_notation = self.get_simple_rna_notation_from_single_nucleotide_sequence(
            single_oligo_sequence
        )
        
        # Build complex notation from simple notation
        node_id = "RNA1"
        return f"{node_id}{{{simple_notation}}}$$$"

    def get_simple_rna_notation_from_single_nucleotide_sequence(
        self, single_oligo_sequence: str
    ) -> str:
        """Convert single oligonucleotide sequence to simple RNA notation.
        
        Args:
            single_oligo_sequence: Single sequence like "AUGCUU"
            
        Returns:
            Simple RNA notation (e.g., "R(A)P.R(U)P.R(G)P.R(C)P.R(U)P.R(U)P")
            
        Raises:
            NotationException: If unknown nucleotide is found
        """
        factory = MonomerFactory.get_instance()
        rna_monomers = factory.get_monomer_db().get(NUCLIEC_ACID_POLYMER_TYPE, {})

        # Find linker monomers (those with natural_analog = "P")
        # These are identified by their ID (like "P" for phosphate)
        linker_ids: List[str] = []
        for key, m in rna_monomers.items():
            if m.natural_analog == ID_P:
                linker_ids.append(key)

        # Find starting linker
        sb: List[str] = []
        seq = single_oligo_sequence.upper()
        
        for key in sorted(linker_ids, key=len, reverse=True):
            if seq.startswith(key):
                # Look up the notation from NucleotideFactory
                nf = NucleotideFactory.get_instance()
                templates = nf.get_nucleotide_templates().get(constants.NOTATION_SOURCE, {})
                notation_val = templates.get(key, "")
                if not notation_val:
                    notation_val = f"R(P)P"  # fallback
                sb.append(notation_val)
                seq = seq[len(key):]
                break

        # Process remaining nucleotides
        nucleotides = NucleotideFactory.get_instance().get_nucleotide_templates().get(
            constants.NOTATION_SOURCE, {}
        )

        while seq:
            if sb:
                sb.append(".")

            found = False
            for nuc in sorted(nucleotides.keys(), key=len, reverse=True):
                if seq.startswith(nuc):
                    sb.append(nucleotides[nuc])
                    seq = seq[len(nuc):]
                    found = True
                    break

            if not found:
                raise NotationException(f"Unknown nucleotide found starting at {seq}")

        result = ".".join(sb)
        
        # Remove trailing P
        return remove_last_p(result)

    def _combine_complex_notations(self, notation1: str, notation2: str) -> str:
        """Combine two complex notations."""
        # Parse the notations
        from helm_notation.complex_notation_parser import _parse_extended_notation
        
        nodes1, edges1, bp1, ann1, other1 = _parse_extended_notation(notation1)
        nodes2, edges2, bp2, ann2, other2 = _parse_extended_notation(notation2)
        
        # Renumber nodes in second notation
        offset = len(nodes1)
        
        new_nodes = list(nodes1)
        for node in nodes2:
            # Extract the type and number, then renumber
            import re
            match = re.match(r"([A-Za-z]+)(\d+)", node.node_id)
            if match:
                new_type = match.group(1)
                new_num = int(match.group(2)) + offset
                new_id = f"{new_type}{new_num}"
            else:
                new_id = node.node_id
            new_nodes.append(type(node)(node_id=new_id, polymer_type=node.polymer_type, label=node.label))
        
        # Simple combination - just concatenate nodes
        nodes_part = "|".join([f"{n.node_id}{{{n.label}}}" for n in new_nodes])
        
        # If there were edges in either, combine them
        edges_part = ""
        if edges1:
            edges_part = "|" + "|".join([f"{e.source},{e.target},{e.connection}" for e in edges1])
        if edges2:
            # Need to renumber edges from second notation
            if not edges_part:
                edges_part = "|"
            for edge in edges2:
                import re
                for src in edge.source.split("+"):
                    match = re.match(r"([A-Za-z]+)(\d+)", src)
                    if match:
                        new_num = int(match.group(2)) + offset
                        src = f"{match.group(1)}{new_num}"
                for tgt in edge.target.split("+"):
                    match = re.match(r"([A-Za-z]+)(\d+)", tgt)
                    if match:
                        new_num = int(match.group(2)) + offset
                        tgt = f"{match.group(1)}{new_num}"
                edges_part += f"|{src},{tgt},{edge.connection}"
        
        return f"{nodes_part}${edges_part}$$$"

    def _hybridize(self, notation: str) -> str:
        """Add base pair annotations for complementary sequences."""
        # This is a simplified version - full implementation would calculate
        # complementary pairs between strands
        return notation