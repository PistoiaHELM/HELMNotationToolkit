"""Extinction coefficient calculator for RNA and peptide polymers.

Calculates extinction coefficients from nucleotide sequences,
amino acid sequences, and HELM polymer notations.
"""

from __future__ import annotations

import importlib.resources
import logging
from typing import TYPE_CHECKING

from helm_notation.exceptions import (
    CalculationException,
    MonomerException,
    NotationException,
    StructureException,
)
from helm_notation.model.monomer import (
    NUCLEIC_ACID_POLYMER_TYPE,
    PEPTIDE_POLYMER_TYPE,
)

if TYPE_CHECKING:
    from helm_notation.monomer_store import MonomerStore

logger = logging.getLogger(__name__)

RNA_UNIT_TYPE = 1
PEPTIDE_UNIT_TYPE = 2

_mono_nucleotide_map: dict[str, float] = {}
_di_nucleotide_map: dict[str, float] = {}
_amino_acid_map: dict[str, float] = {}

_RNA_PROPERTY_FILE = "RNAExtinctionCoefficient.properties"
_PEPTIDE_PROPERTY_FILE = "PEPTIDEExtinctionCoefficient.properties"

# Module-level singleton
_instance: ExtinctionCoefficientCalculator | None = None


class ExtinctionCoefficientCalculator:
    """Calculator for extinction coefficients of RNA and peptide polymers."""

    def __init__(self) -> None:
        pass

    def get_unit(self, unit_type: int) -> str | None:
        """Return the unit string for the given unit type."""
        if unit_type == RNA_UNIT_TYPE:
            return self.get_rna_unit()
        elif unit_type == PEPTIDE_UNIT_TYPE:
            return self.get_peptide_unit()
        return None

    def get_rna_unit(self) -> str:
        return "mM-1cm-1"

    def get_peptide_unit(self) -> str:
        return "M-1cm-1"

    def get_default_unit(self) -> str:
        return self.get_rna_unit()

    def get_default_unit_type(self) -> int:
        return RNA_UNIT_TYPE

    def calculate_from_complex_notation(
        self,
        complex_notation: str,
        unit_type: int | None = None,
        monomer_store: MonomerStore | None = None,
    ) -> float:
        """Calculate extinction coefficient for complex polymer notation.

        Args:
            complex_notation: Complex polymer notation string.
            unit_type: RNA_UNIT_TYPE or PEPTIDE_UNIT_TYPE (default RNA).
            monomer_store: Optional MonomerStore (uses global if None).

        Returns:
            Extinction coefficient value.
        """
        if unit_type is None:
            unit_type = self.get_default_unit_type()

        if monomer_store is None:
            from helm_notation.monomer_factory import get_instance

            try:
                factory = get_instance()
            except Exception as ex:
                raise NotationException(
                    "Unable to initialize monomer factory"
                ) from ex
            monomer_store = factory.monomer_store

        from helm_notation.tools.complex_notation_parser import get_polymer_node_list

        result = 0.0
        polymer_nodes = get_polymer_node_list(complex_notation, monomer_store)

        for polymer_node in polymer_nodes:
            polymer_type = polymer_node.type
            notation = polymer_node.label
            ext = 0.0
            if polymer_type == NUCLEIC_ACID_POLYMER_TYPE:
                ext = self.calculate_from_rna_polymer_notation(
                    notation, monomer_store
                )
                if unit_type == PEPTIDE_UNIT_TYPE:
                    ext = ext * 1000
            elif polymer_type == PEPTIDE_POLYMER_TYPE:
                ext = self.calculate_from_peptide_polymer_notation(
                    notation, monomer_store
                )
                if unit_type == RNA_UNIT_TYPE:
                    ext = ext / 1000
            result += ext

        return result

    def calculate_from_rna_polymer_notation(
        self,
        simple_notation: str,
        monomer_store: MonomerStore | None = None,
    ) -> float:
        """Calculate extinction coefficient for RNA polymer notation.

        Args:
            simple_notation: RNA polymer notation string.
            monomer_store: Optional MonomerStore.

        Returns:
            Extinction coefficient value.
        """
        if monomer_store is None:
            from helm_notation.monomer_factory import get_instance

            monomer_store = get_instance().monomer_store

        from helm_notation.tools.simple_notation_parser import (
            get_trimmed_nucleotide_sequence,
        )

        natural_sequence = get_trimmed_nucleotide_sequence(
            simple_notation, monomer_store
        )
        return self.calculate_from_nucleotide_sequence(natural_sequence)

    def calculate_from_modified_nucleotide_sequence(
        self, sequence: str
    ) -> float:
        """Calculate extinction coefficient for modified nucleotide sequence.

        Args:
            sequence: Modified nucleotide sequence.

        Returns:
            Extinction coefficient value.
        """
        from helm_notation.tools.nucleotide_sequence_parser import (
            get_natural_analog_sequence,
        )

        natural_sequence = get_natural_analog_sequence(sequence)
        return self.calculate_from_nucleotide_sequence(natural_sequence)

    def calculate_from_nucleotide_sequence(self, sequence: str) -> float:
        """Calculate extinction coefficient for natural nucleotide sequence.

        Args:
            sequence: Natural nucleotide sequence.

        Returns:
            Extinction coefficient value.
        """
        if not sequence:
            raise CalculationException("Input sequence cannot be null")

        if len(sequence) == 1:
            if sequence in _mono_nucleotide_map:
                return _mono_nucleotide_map[sequence]
            else:
                raise CalculationException(
                    f"Unknown nucleotide [{sequence}] found"
                )

        mono_map: dict[str, int] = {}
        for i in range(1, len(sequence) - 1):
            mono = sequence[i]
            if mono in _mono_nucleotide_map:
                mono_map[mono] = mono_map.get(mono, 0) + 1
            else:
                raise CalculationException(
                    f"Unknown nucleotide [{mono}] found"
                )

        di_map: dict[str, int] = {}
        for i in range(len(sequence) - 1):
            di = sequence[i : i + 2]
            if di in _di_nucleotide_map:
                di_map[di] = di_map.get(di, 0) + 1
            else:
                raise CalculationException(
                    f"Unknown dinucleotide [{di}] found"
                )

        return _calculate_nucleotide(mono_map, di_map)

    def calculate_from_amino_acid_sequence(self, sequence: str) -> float:
        """Calculate extinction coefficient of single letter amino acid sequence.

        Args:
            sequence: Single letter amino acid sequence.

        Returns:
            Extinction coefficient value.
        """
        if not sequence:
            return 0.0

        ids = list(sequence)
        return _calculate_amino_acid(ids)

    def calculate_from_peptide_polymer_notation(
        self,
        simple_notation: str,
        monomer_store: MonomerStore | None = None,
    ) -> float:
        """Calculate extinction coefficient for peptide polymer notation.

        Args:
            simple_notation: Peptide polymer notation string.
            monomer_store: Optional MonomerStore.

        Returns:
            Extinction coefficient value.
        """
        if monomer_store is None:
            from helm_notation.monomer_factory import get_instance

            try:
                factory = get_instance()
            except Exception as ex:
                raise NotationException(
                    "Unable to initialize monomer factory"
                ) from ex
            monomer_store = factory.monomer_store

        from helm_notation.tools.simple_notation_parser import get_monomer_list

        monomers = get_monomer_list(
            simple_notation, PEPTIDE_POLYMER_TYPE, monomer_store
        )
        ids = [monomer.natural_analog for monomer in monomers]
        return _calculate_amino_acid(ids)


def _calculate_nucleotide(
    tmp_mono_nucleotide_map: dict[str, int],
    tmp_di_nucleotide_map: dict[str, int],
) -> float:
    """Calculate extinction coefficient from mono and di nucleotide counts."""
    mono_result = 0.0
    for key, count in tmp_mono_nucleotide_map.items():
        value = _mono_nucleotide_map[key]
        mono_result += value * count

    di_result = 0.0
    for key, count in tmp_di_nucleotide_map.items():
        value = _di_nucleotide_map[key]
        di_result += value * count

    return 2 * di_result - mono_result


def _calculate_amino_acid(aa_id_list: list[str] | None) -> float:
    """Calculate extinction coefficient from amino acid ID list."""
    if not aa_id_list:
        return 0.0

    count_map: dict[str, int] = {}
    for aa_id in aa_id_list:
        if aa_id in _amino_acid_map:
            count_map[aa_id] = count_map.get(aa_id, 0) + 1

    result = 0.0
    for key, count in count_map.items():
        factor = _amino_acid_map[key]
        result += factor * count

    return result


def _init_maps() -> None:
    """Initialize extinction coefficient maps from property files."""
    global _mono_nucleotide_map, _di_nucleotide_map, _amino_acid_map

    # Load RNA properties
    rna_data = (
        importlib.resources.files("helm_notation.resources")
        .joinpath(_RNA_PROPERTY_FILE)
        .read_text(encoding="utf-8")
    )
    for line in rna_data.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" in line:
            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()
            f = float(value)
            if len(key) == 1:
                _mono_nucleotide_map[key] = f
            elif len(key) == 2:
                _di_nucleotide_map[key] = f

    # Load peptide properties
    peptide_data = (
        importlib.resources.files("helm_notation.resources")
        .joinpath(_PEPTIDE_PROPERTY_FILE)
        .read_text(encoding="utf-8")
    )
    for line in peptide_data.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" in line:
            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()
            f = float(value)
            _amino_acid_map[key] = f


def get_instance() -> ExtinctionCoefficientCalculator:
    """Get or create the singleton ExtinctionCoefficientCalculator instance."""
    global _instance
    if _instance is None:
        _instance = ExtinctionCoefficientCalculator()
        try:
            _init_maps()
        except IOError:
            raise CalculationException(
                "Unable to initialize extinction coefficient property files"
            )
    return _instance
