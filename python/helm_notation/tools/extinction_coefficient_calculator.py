from __future__ import annotations

from dataclasses import dataclass
import struct
from threading import Lock

from helm_notation.exceptions import CalculationException
from helm_notation.nucleotide_sequence_parser import get_natural_analog_sequence
from helm_notation.resources import read_text


RNA_UNIT_TYPE = 1
PEPTIDE_UNIT_TYPE = 2


@dataclass
class _Maps:
    mono_nucleotide: dict[str, float]
    di_nucleotide: dict[str, float]
    amino_acid: dict[str, float]


class ExtinctionCoefficientCalculator:
    _instance: "ExtinctionCoefficientCalculator | None" = None
    _lock = Lock()

    def __init__(self) -> None:
        self._maps = _Maps(mono_nucleotide={}, di_nucleotide={}, amino_acid={})
        self._init_maps()

    @classmethod
    def get_instance(cls) -> "ExtinctionCoefficientCalculator":
        with cls._lock:
            if cls._instance is None:
                cls._instance = cls()
            return cls._instance

    def get_unit(self, unit_type: int) -> str | None:
        if unit_type == RNA_UNIT_TYPE:
            return self.get_rna_unit()
        if unit_type == PEPTIDE_UNIT_TYPE:
            return self.get_peptide_unit()
        return None

    @staticmethod
    def get_rna_unit() -> str:
        return "mM-1cm-1"

    @staticmethod
    def get_peptide_unit() -> str:
        return "M-1cm-1"

    def get_default_unit(self) -> str:
        return self.get_rna_unit()

    @staticmethod
    def get_default_unit_type() -> int:
        return RNA_UNIT_TYPE

    def calculate_from_modified_nucleotide_sequence(self, sequence: str) -> float:
        natural_sequence = get_natural_analog_sequence(sequence)
        return self.calculate_from_nucleotide_sequence(natural_sequence)

    def calculate_from_nucleotide_sequence(self, sequence: str) -> float:
        if not sequence:
            raise CalculationException("Input sequence cannot be null")

        sequence = sequence.upper()
        mono_map = self._maps.mono_nucleotide
        di_map = self._maps.di_nucleotide

        if len(sequence) == 1:
            if sequence in mono_map:
                return _f32(mono_map[sequence])
            raise CalculationException(f"Unknown nucleotide [{sequence}] found")

        internal_counts: dict[str, int] = {}
        for index in range(1, len(sequence) - 1):
            mono = sequence[index : index + 1]
            if mono not in mono_map:
                raise CalculationException(f"Unknown nucleotide [{mono}] found")
            internal_counts[mono] = internal_counts.get(mono, 0) + 1

        dinucleotide_counts: dict[str, int] = {}
        for index in range(0, len(sequence) - 1):
            di = sequence[index : index + 2]
            if di not in di_map:
                raise CalculationException(f"Unknown dinucleotide [{di}] found")
            dinucleotide_counts[di] = dinucleotide_counts.get(di, 0) + 1

        mono_result = _f32(0.0)
        for key, count in internal_counts.items():
            mono_result = _f32(mono_result + _f32(mono_map[key]) * _f32(count))

        di_result = _f32(0.0)
        for key, count in dinucleotide_counts.items():
            di_result = _f32(di_result + _f32(di_map[key]) * _f32(count))

        return _f32(_f32(2.0) * di_result - mono_result)

    def calculate_from_amino_acid_sequence(self, sequence: str) -> float:
        if not sequence:
            return 0.0

        counts: dict[str, int] = {}
        for residue in sequence.upper():
            if residue in self._maps.amino_acid:
                counts[residue] = counts.get(residue, 0) + 1

        return float(
            sum(self._maps.amino_acid[key] * count for key, count in counts.items())
        )

    def _init_maps(self) -> None:
        rna_text = read_text("RNAExtinctionCoefficient.properties")
        peptide_text = read_text("PEPTIDEExtinctionCoefficient.properties")

        for key, value in _parse_properties(rna_text).items():
            as_float = float(value)
            if len(key) == 1:
                self._maps.mono_nucleotide[key] = as_float
            elif len(key) == 2:
                self._maps.di_nucleotide[key] = as_float

        for key, value in _parse_properties(peptide_text).items():
            self._maps.amino_acid[key] = float(value)


def _parse_properties(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        out[key.strip()] = value.strip()
    return out


def _f32(value: float) -> float:
    return struct.unpack("!f", struct.pack("!f", float(value)))[0]
