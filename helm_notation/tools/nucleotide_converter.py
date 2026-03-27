"""Converter between nucleotide sequences and HELM notation.

Provides bidirectional conversion between nucleotide sequences
and HELM complex/simple notation.
"""

from __future__ import annotations

import logging

from helm_notation.constants import NOTATION_SOURCE
from helm_notation.exceptions import (
    MonomerException,
    NotationException,
    StructureException,
)
from helm_notation.model.monomer import (
    ID_P,
    NUCLEIC_ACID_POLYMER_TYPE,
)
from helm_notation.model.nucleotide import Nucleotide

logger = logging.getLogger(__name__)

# Module-level singleton
_instance: NucleotideConverter | None = None


class NucleotideConverter:
    """Converter between nucleotide sequences and HELM notation."""

    def __init__(self) -> None:
        pass

    def get_nucleotide_sequences_from_complex_notation(
        self, complex_notation: str
    ) -> str:
        """Return nucleotide sequences for nucleic acid polymer mixture.

        Args:
            complex_notation: HELM notation for nucleic acid polymer (may be mixture).

        Returns:
            Nucleotide sequence for nucleic acid, separated with space between strands.
        """
        from helm_notation.tools.complex_notation_parser import get_polymer_node_list

        polymer_nodes = get_polymer_node_list(complex_notation)

        parts: list[str] = []
        for node in polymer_nodes:
            polymer_type = node.type
            if polymer_type != NUCLEIC_ACID_POLYMER_TYPE:
                raise NotationException(
                    "Input complex notation contains non-nucleic acid polymer"
                )
            simple_notation = node.label
            notation = self.get_nucleotide_sequence_from_simple_rna_notation(
                simple_notation
            )
            parts.append(notation)
        return " ".join(parts)

    def get_nucleotide_sequence_from_simple_rna_notation(
        self, simple_rna_notation: str
    ) -> str:
        """Return nucleotide sequence for a single simple RNA notation.

        Args:
            simple_rna_notation: Simple HELM notation for RNA.

        Returns:
            Nucleotide sequence string.
        """
        from helm_notation.tools.simple_notation_parser import get_nucleotide_list

        nuc_list = get_nucleotide_list(simple_rna_notation, last_p=False)

        from helm_notation.nucleotide_factory import get_instance as get_nuc_factory

        count = 0
        sb: list[str] = []
        reverse_nuc_map = get_nuc_factory().get_reverse_nucleotide_template_map()

        for nuc in nuc_list:
            nucleotide = nuc.notation
            nucleoside = nuc.nucleoside_notation
            linker = nuc.linker_notation

            # It is ok for the first nucleotide not to have a nucleoside
            if count == 0 and len(nucleoside) == 0:
                phosphate_monomer = nuc.get_phosphate_monomer()
                if phosphate_monomer is not None:
                    sb.append(phosphate_monomer.alternate_id)
                count += 1
                continue

            # It is ok for the last nucleotide not to have a linker
            if count == len(nuc_list) - 1 and len(linker) == 0:
                nucleotide = nucleotide + ID_P

            if nucleotide in reverse_nuc_map:
                sb.append(reverse_nuc_map[nucleotide])
            else:
                raise NotationException(
                    f"Unknown nucleotide found for {nucleotide} : missing nucleotide template"
                )

            count += 1

        return "".join(sb)

    def get_complex_notation(self, nucleotide_sequences: str) -> str | None:
        """Return HELM complex notation for nucleotide sequences.

        Sequences may have more than one strand separated by space.

        Args:
            nucleotide_sequences: Nucleotide sequences (space-separated for multiple strands).

        Returns:
            Complex notation string, or None if input is empty.
        """
        if not nucleotide_sequences or len(nucleotide_sequences.strip()) == 0:
            return None

        from helm_notation.tools.complex_notation_parser import (
            get_combined_complex_notation,
            hybridize,
        )

        singles = nucleotide_sequences.split()
        result = self.get_complex_notation_from_single_nucleotide_sequence(
            singles[0]
        )

        for i in range(1, len(singles)):
            notation = self.get_complex_notation_from_single_nucleotide_sequence(
                singles[i]
            )
            result = get_combined_complex_notation(result, notation)

        result = hybridize(result)
        return result

    def get_complex_notation_from_single_nucleotide_sequence(
        self, single_oligo_sequence: str
    ) -> str:
        """Return complex notation for a single oligonucleotide sequence.

        Args:
            single_oligo_sequence: A single oligonucleotide sequence.

        Returns:
            Complex notation string.
        """
        simple_notation = (
            self.get_simple_rna_notation_from_single_nucleotide_sequence(
                single_oligo_sequence
            )
        )
        from helm_notation.tools.simple_notation_parser import (
            get_complex_notation_for_rna,
        )

        return get_complex_notation_for_rna(simple_notation)

    def get_simple_rna_notation_from_single_nucleotide_sequence(
        self, single_oligo_sequence: str
    ) -> str:
        """Return simple RNA notation for a single oligonucleotide sequence.

        Args:
            single_oligo_sequence: A single oligonucleotide sequence.

        Returns:
            Simple RNA notation string.
        """
        from helm_notation.monomer_factory import get_instance as get_monomer_factory
        from helm_notation.nucleotide_factory import get_instance as get_nuc_factory
        from helm_notation.tools.nucleotide_sequence_parser import remove_last_p

        rna_monomers = get_monomer_factory().monomer_db.get(
            NUCLEIC_ACID_POLYMER_TYPE, {}
        )
        linker_monomers: dict[str, object] = {}
        for key, m in rna_monomers.items():
            if m.natural_analog == ID_P:
                linker_monomers[key] = m

        sb: list[str] = []
        seq = single_oligo_sequence

        # Check if sequence starts with a linker monomer
        for key in linker_monomers:
            if seq.startswith(key):
                if len(key) > 1:
                    sb.append(f"[{key}]")
                else:
                    sb.append(key)
                seq = seq[len(key) :]
                break

        nucleotides = get_nuc_factory().get_nucleotide_templates().get(
            NOTATION_SOURCE, {}
        )
        nuc_set = list(nucleotides.keys())

        while len(seq) > 0:
            if len(sb) > 0:
                sb.append(".")

            found = False
            for nuc in nuc_set:
                if seq.startswith(nuc):
                    sb.append(nucleotides[nuc])
                    seq = seq[len(nuc) :]
                    found = True
                    break

            if not found:
                raise NotationException(
                    f"Unknown nucleotide found starting at {seq}"
                )

        return remove_last_p("".join(sb)) or ""


def get_instance() -> NucleotideConverter:
    """Get or create the singleton NucleotideConverter instance."""
    global _instance
    if _instance is None:
        # Ensure factories are initialized
        from helm_notation.monomer_factory import get_instance as get_monomer_factory
        from helm_notation.nucleotide_factory import get_instance as get_nuc_factory

        get_monomer_factory()
        get_nuc_factory()
        _instance = NucleotideConverter()
    return _instance
