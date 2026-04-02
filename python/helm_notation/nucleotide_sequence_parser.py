from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

from helm_notation import constants
from helm_notation.exceptions import NotationException
from helm_notation.monomer_factory import MonomerFactory
from helm_notation.nucleotide_factory import NucleotideFactory

MINIMUM_MATCH_FRAGMENT_LENGTH = 2

RNA_DESIGN_NONE = "NONE"
RNA_DESIGN_TUSCHL_19_PLUS_2 = "TUSCHL_19_PLUS_2"
RNA_DESIGN_DICER_27_R = "DICER_27_R"
RNA_DESIGN_DICER_27_L = "DICER_27_L"

SUPPORTED_DESIGN_LIST = [
    RNA_DESIGN_NONE,
    RNA_DESIGN_TUSCHL_19_PLUS_2,
    RNA_DESIGN_DICER_27_L,
    RNA_DESIGN_DICER_27_R,
]

COMPLEMENT_MAP: Dict[str, str] = {
    "A": "U",
    "G": "C",
    "C": "G",
    "U": "A",
    "T": "A",
    "X": "X",
}


@dataclass
class _NucleotideToken:
    symbol: str
    notation: str
    natural_analog: str


def is_normal_direction(sequence: str) -> bool:
    return not sequence.startswith("3")


def cleanup(sequence: str) -> str:
    condensed = "".join(sequence.split())
    tokens = condensed.split("-")

    if tokens and tokens[0].startswith(("5", "3")) and len(tokens) > 1:
        result = tokens[1]
    else:
        result = tokens[0] if tokens else ""

    if result == result.lower():
        result = result.upper()
    return result


def get_notation(sequence: str, notation_source: str = constants.NOTATION_SOURCE) -> str:
    normal = get_normal_list(sequence, notation_source)
    return ".".join(n.notation for n in normal)


def get_normal_list(
    sequence: str,
    notation_source: str = constants.NOTATION_SOURCE,
) -> List[_NucleotideToken]:
    in_normal = is_normal_direction(sequence)
    clean = cleanup(sequence)
    nuc_list = _get_nucleotide_list(clean, notation_source)
    return nuc_list if in_normal else list(reversed(nuc_list))


def get_reverse_complement_sequence(
    sequence: str,
    notation_source: str = constants.NOTATION_SOURCE,
) -> str:
    return _get_complement_sequence(sequence, notation_source, out_normal_direction=False)


def get_normal_complement_sequence(
    sequence: str,
    notation_source: str = constants.NOTATION_SOURCE,
) -> str:
    return _get_complement_sequence(sequence, notation_source, out_normal_direction=True)


def _get_complement_sequence(
    sequence: str,
    notation_source: str,
    out_normal_direction: bool,
) -> str:
    clean = cleanup(sequence)
    nuc_list = _get_nucleotide_list(clean, notation_source)
    comp_list = _get_complement_list(nuc_list)

    in_normal_direction = is_normal_direction(sequence)
    if in_normal_direction:
        if out_normal_direction:
            body = "".join(n.natural_analog for n in reversed(comp_list))
            return f"5'-{body}-3'"
        body = "".join(n.natural_analog for n in comp_list)
        return f"3'-{body}-5'"

    if out_normal_direction:
        body = "".join(n.natural_analog for n in comp_list)
        return f"5'-{body}-3'"
    body = "".join(n.natural_analog for n in reversed(comp_list))
    return f"3'-{body}-5'"


def get_reverse_sequence(sequence: str) -> str:
    nuc_list = _get_nucleotide_list(sequence, constants.NOTATION_SOURCE)
    return "".join(n.symbol for n in reversed(nuc_list))


def get_natural_analog_sequence(sequence: str) -> str:
    nuc_list = _get_nucleotide_list(sequence, constants.NOTATION_SOURCE)
    return "".join(n.natural_analog for n in nuc_list)


def is_in_opposite_direction(sense_seq: str, antisense_seq: str) -> bool:
    upper_ss = sense_seq.upper()
    upper_as = antisense_seq.upper()
    normal_comp_as = get_normal_complement_sequence(upper_as)
    reverse_comp_as = get_reverse_complement_sequence(upper_as)

    normal_max = get_max_match_fragment(upper_ss, normal_comp_as)
    reverse_max = get_max_match_fragment(upper_ss, reverse_comp_as)
    return len(normal_max) < len(reverse_max)


def get_max_match_fragment(
    seq1: str,
    seq2: str,
    min_length: int = MINIMUM_MATCH_FRAGMENT_LENGTH,
) -> str:
    if seq1 is None or seq2 is None:
        raise NotationException("Both sequences must not be null")

    if seq1 != seq1.upper() or seq2 != seq2.upper():
        raise NotationException(
            "Both sequences must be natural nucleotide sequence in upper case"
        )

    long_seq, short_seq = (seq1, seq2) if len(seq1) > len(seq2) else (seq2, seq1)
    long_seq = long_seq.replace("T", "U")
    short_seq = short_seq.replace("T", "U")

    minimum = max(MINIMUM_MATCH_FRAGMENT_LENGTH, min_length)
    for length in range(len(short_seq), minimum, -1):
        for idx in range(0, len(short_seq) - length + 1):
            fragment = short_seq[idx : idx + length]
            if fragment in long_seq:
                return fragment
    return ""


def remove_last_p(simple_rna_notation: str) -> str:
    if simple_rna_notation and simple_rna_notation.endswith("P"):
        return simple_rna_notation[:-1]
    return simple_rna_notation


def get_sirna_notation(
    sense_seq: str,
    antisense_seq: str,
    rna_design_type: str = RNA_DESIGN_NONE,
) -> str:
    _validate_sirna_design(sense_seq, antisense_seq, rna_design_type)

    tokens = ["", ""]
    count = 0

    if sense_seq:
        count += 1
        ss_notation = remove_last_p(get_notation(sense_seq))
        tokens[0] = f"RNA{count}{{{ss_notation}}}"
        tokens[1] = f"RNA{count}{{ss}}"

    if antisense_seq:
        count += 1
        as_notation = remove_last_p(get_notation(antisense_seq))
        tokens[0] = f"{tokens[0]}|RNA{count}{{{as_notation}}}" if tokens[0] else f"RNA{count}{{{as_notation}}}"
        tokens[1] = f"{tokens[1]}|RNA{count}{{as}}" if tokens[1] else f"RNA{count}{{as}}"

    base_pair = _hybridization(sense_seq, antisense_seq, rna_design_type)
    if not tokens[0]:
        return ""

    parts = [tokens[0], "", base_pair if base_pair else "", tokens[1], ""]
    return "$".join(parts)


def _hybridization(sense_seq: str, antisense_seq: str, rna_design_type: str) -> str:
    if not sense_seq or not antisense_seq:
        return ""

    analog_ss = get_natural_analog_sequence(sense_seq).replace("T", "U")
    analog_as = get_natural_analog_sequence(antisense_seq).replace("T", "U")

    if RNA_DESIGN_NONE.upper() == rna_design_type.upper():
        normal_comp_as = cleanup(get_normal_complement_sequence(analog_as))
        max_match = get_max_match_fragment(analog_ss, normal_comp_as)
        if not max_match:
            return ""

        ss_start = analog_ss.find(max_match)
        normal_comp_start = normal_comp_as.find(max_match)
        as_start = len(analog_as) - len(max_match) - normal_comp_start

        pairs: List[str] = []
        for i in range(len(max_match)):
            ss_pos = (i + ss_start) * 3 + 2
            as_pos = (as_start + len(max_match) - 1 - i) * 3 + 2
            pairs.append(f"RNA1,RNA2,{ss_pos}:pair-{as_pos}:pair")
        return "|".join(pairs)

    if RNA_DESIGN_TUSCHL_19_PLUS_2.upper() == rna_design_type.upper():
        return _hybridization_with_length_from_start(analog_ss, analog_as, 19)
    if RNA_DESIGN_DICER_27_R.upper() == rna_design_type.upper():
        return _hybridization_with_length_from_start(analog_ss, analog_as, 25)
    if RNA_DESIGN_DICER_27_L.upper() == rna_design_type.upper():
        return _hybridization_with_length_from_start(analog_ss, analog_as, 25)

    raise NotationException("Unknow RNA Design")


def _hybridization_with_length_from_start(
    sense_analog_seq: str,
    antisense_analog_seq: str,
    length_from_start: int,
) -> str:
    limit = min(length_from_start, len(sense_analog_seq), len(antisense_analog_seq))
    pairs: List[str] = []
    for i in range(limit):
        ss_pos = i * 3 + 2
        as_pos = (length_from_start - 1 - i) * 3 + 2
        ss_char = sense_analog_seq[i]
        as_char = antisense_analog_seq[length_from_start - 1 - i]

        if COMPLEMENT_MAP.get(ss_char, "X").upper() == as_char.upper():
            pairs.append(f"RNA1,RNA2,{ss_pos}:pair-{as_pos}:pair")
    return "|".join(pairs)


def _validate_sirna_design(sense_seq: str, antisense_seq: str, rna_design_type: str) -> bool:
    if RNA_DESIGN_NONE.upper() == rna_design_type.upper():
        return True

    if rna_design_type not in SUPPORTED_DESIGN_LIST:
        raise NotationException(f"Unsupported RNA Design Type '{rna_design_type}'")

    sense_nuc_list = _get_nucleotide_list(sense_seq, constants.NOTATION_SOURCE)
    antisense_nuc_list = _get_nucleotide_list(antisense_seq, constants.NOTATION_SOURCE)

    if rna_design_type == RNA_DESIGN_TUSCHL_19_PLUS_2:
        if len(sense_nuc_list) != 21:
            raise NotationException(
                "Sense strand for Tuschl 19+2 design must have 21 nucleotides"
            )
        if len(antisense_nuc_list) != 21:
            raise NotationException(
                "Antisense strand for Tuschl 19+2 design must have 21 nucleotides"
            )
    elif rna_design_type == RNA_DESIGN_DICER_27_R:
        if len(sense_nuc_list) != 25:
            raise NotationException(
                "Sense strand for Dicer 27R design must have 25 nucleotides"
            )
        if len(antisense_nuc_list) != 27:
            raise NotationException(
                "Antisense strand for Dicer 27R design must have 27 nucleotides"
            )
    elif rna_design_type == RNA_DESIGN_DICER_27_L:
        if len(sense_nuc_list) != 27:
            raise NotationException(
                "Sense strand for Dicer 27L design must have 27 nucleotides"
            )
        if len(antisense_nuc_list) != 25:
            raise NotationException(
                "Antisense strand for Dicer 27L design must have 25 nucleotides"
            )
    return True


def _get_complement_list(in_list: List[_NucleotideToken]) -> List[_NucleotideToken]:
    out: List[_NucleotideToken] = []
    for nucleotide in in_list:
        symbol = COMPLEMENT_MAP.get(nucleotide.natural_analog, "X")
        notation = f"R({symbol})P"
        out.append(_NucleotideToken(symbol=symbol, notation=notation, natural_analog=symbol))
    return out


def _get_nucleotide_list(
    non_direction_sequence: str,
    notation_source: str = constants.NOTATION_SOURCE,
) -> List[_NucleotideToken]:
    if non_direction_sequence is None:
        raise NotationException("Nucleotide Sequence must be specified")
    if notation_source is None:
        raise NotationException("Notation Source must be specified")

    templates = NucleotideFactory.get_instance().get_nucleotide_templates()
    if notation_source not in templates:
        raise NotationException(f"Unknown Notation Source [{notation_source}]")

    nucleotides = templates[notation_source]
    symbols = sorted(nucleotides.keys(), key=lambda s: (-len(s), s))

    out: List[_NucleotideToken] = []
    pos = 0
    while pos < len(non_direction_sequence):
        found = False
        for symbol in symbols:
            if non_direction_sequence.startswith(symbol, pos):
                notation = nucleotides[symbol]
                out.append(
                    _NucleotideToken(
                        symbol=symbol,
                        notation=notation,
                        natural_analog=_infer_natural_analog(symbol, notation),
                    )
                )
                pos += len(symbol)
                found = True
                break
        if not found:
            tail = non_direction_sequence[pos:]
            raise NotationException(
                f"Sequence contains unknown nucleotide starting at {tail}"
            )
    return out


def _infer_natural_analog(symbol: str, notation: str) -> str:
    if notation and "(" in notation and ")" in notation:
        base = notation[notation.find("(") + 1 : notation.find(")")].replace("[", "").replace("]", "")
        if base:
            monomer = (
                MonomerFactory.get_instance()
                .get_monomer_db()
                .get("RNA", {})
                .get(base)
            )
            if monomer and monomer.natural_analog:
                return monomer.natural_analog
            if len(base) == 1:
                return base.upper()

    if symbol:
        last = symbol[-1].upper()
        if last in COMPLEMENT_MAP:
            return last
    return "X"
