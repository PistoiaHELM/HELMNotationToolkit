"""Methods for handling nucleotide sequences.

Converts between nucleotide sequences and HELM notation,
provides complement sequence generation, and siRNA notation support.
"""

from __future__ import annotations

import re
from typing import Optional

from lxml import etree

from helm_notation.constants import NOTATION_SOURCE
from helm_notation.exceptions import NotationException
from helm_notation.model.nucleotide import Nucleotide

DEFAULT_NOTATION_SOURCE = "HELM Notation"  # deprecated
MINIMUM_MATCH_FRAGMENT_LENGTH = 2

complement_map: dict[str, str] = {
    "A": "U",
    "G": "C",
    "C": "G",
    "U": "A",
    "T": "A",
    "X": "X",
}

RNA_DESIGN_NONE = "NONE"
RNA_DESIGN_TUSCHL_19_PLUS_2 = "TUSCHL_19_PLUS_2"
RNA_DESIGN_DICER_27_R = "DICER_27_R"
RNA_DESIGN_DICER_27_L = "DICER_27_L"

SUPPORTED_DESIGN_LIST: list[str] = [
    RNA_DESIGN_NONE,
    RNA_DESIGN_TUSCHL_19_PLUS_2,
    RNA_DESIGN_DICER_27_L,
    RNA_DESIGN_DICER_27_R,
]

_NUCLEOTIDE_SYMBOL_ELEMENT = "SYMBOL"
_NUCLEOTIDE_MONOMER_NOTATION_ELEMENT = "MONOMER_NOTATION"
_NUCLEOTIDE_ELEMENT = "NUCLEOTIDE"
_TEMPLATE_ELEMENT = "TEMPLATE"
_TEMPLATE_NOTATION_SOURCE_ATTRIBUTE = "notationSource"


def get_reverse_complement_sequence(
    sequence: str, notation_source: str | None = None
) -> str:
    """Return complement sequence in reverse (3'-5') direction."""
    if notation_source is None:
        return _get_complement_sequence(sequence, NOTATION_SOURCE, False)
    return _get_complement_sequence(sequence, notation_source, False)


def get_normal_complement_sequence(
    sequence: str, notation_source: str | None = None
) -> str:
    """Return complement sequence in normal (5'-3') direction."""
    if notation_source is None:
        return _get_complement_sequence(sequence, NOTATION_SOURCE, True)
    return _get_complement_sequence(sequence, notation_source, True)


def _get_complement_sequence(
    sequence: str, notation_source: str, out_normal_direction: bool
) -> str:
    """Generate the complement sequence.

    Args:
        sequence: Input sequence, could be 5'-3' or 3'-5'.
        notation_source: The notation source identifier.
        out_normal_direction: True for 5'-3' output, False for 3'-5'.

    Returns:
        Complement sequence string.
    """
    clean_sequence = cleanup(sequence)
    nuc_list = _get_nucleotide_list(clean_sequence, notation_source)
    comp_list = _get_complement_list(nuc_list)

    sb: list[str] = []
    in_normal_direction = is_normal_direction(sequence)

    if in_normal_direction:
        if out_normal_direction:
            sb.append("5'-")
            for i in range(len(comp_list) - 1, -1, -1):
                sb.append(comp_list[i].get_natural_analog())
            sb.append("-3'")
        else:
            sb.append("3'-")
            for nuc in comp_list:
                sb.append(nuc.get_natural_analog())
            sb.append("-5'")
    else:
        if out_normal_direction:
            sb.append("5'-")
            for nuc in comp_list:
                sb.append(nuc.get_natural_analog())
            sb.append("-3'")
        else:
            sb.append("3'-")
            for i in range(len(comp_list) - 1, -1, -1):
                sb.append(comp_list[i].get_natural_analog())
            sb.append("-5'")

    return "".join(sb)


def get_notation(
    sequence: str, notation_source: str | None = None
) -> str:
    """Convert nucleotide sequence into HELM simple notation.

    Assumes sequence input in one of the following formats:
      - 5'-AGCUdT-3'
      - 3'-dTUCGA-5'
      - AFCUdT (default 5' to 3')

    Args:
        sequence: RNA nucleotide sequence.
        notation_source: Optional notation source.

    Returns:
        HELM simple notation string.
    """
    if notation_source is None:
        notation_source = NOTATION_SOURCE
    normal_nucleotide_list = get_normal_list(sequence, notation_source)

    parts: list[str] = []
    for nucleotide in normal_nucleotide_list:
        parts.append(nucleotide.notation)
    return ".".join(parts)


def is_normal_direction(sequence: str) -> bool:
    """Assume the sequence direction is from 5' to 3', unless started with 3."""
    return not sequence.startswith("3")


def _get_complement_list(in_list: list[Nucleotide]) -> list[Nucleotide]:
    """For non-natural ones, use 'X R(X)P'."""
    out_list: list[Nucleotide] = []
    for nucleotide in in_list:
        natural_analog = nucleotide.get_natural_analog()
        symbol = complement_map.get(natural_analog, "X")
        notation = f"R({symbol})P"
        comp_nu = Nucleotide(symbol=symbol, notation=notation)
        out_list.append(comp_nu)
    return out_list


def get_normal_list(
    sequence: str, notation_source: str | None = None
) -> list[Nucleotide]:
    """Get the normal nucleotide list regardless of input sequence direction."""
    if notation_source is None:
        notation_source = NOTATION_SOURCE
    in_normal_direction = is_normal_direction(sequence)
    clean_sequence = cleanup(sequence)
    nuc_list = _get_nucleotide_list(clean_sequence, notation_source)

    if in_normal_direction:
        return nuc_list
    else:
        return _get_reverse_list(nuc_list)


def cleanup(sequence: str) -> str:
    """Remove non-nucleotide related characters."""
    condensed_sequence = re.sub(r"\s", "", sequence)
    tokens = condensed_sequence.split("-")
    if tokens[0].startswith("5"):
        result = tokens[1]
    elif tokens[0].startswith("3"):
        result = tokens[1]
    else:
        result = tokens[0]

    if result == result.lower():
        result = result.upper()
    return result


def _get_reverse_list(in_list: list[Nucleotide]) -> list[Nucleotide]:
    """Reverse the order of list elements."""
    return list(reversed(in_list))


def _get_nucleotide_list(
    non_direction_sequence: str,
    notation_source: str | None = None,
) -> list[Nucleotide]:
    """Convert non-direction nucleotide sequence into a list of nucleotides.

    The sequence should contain letters only (after cleanup).
    """
    if notation_source is None:
        notation_source = NOTATION_SOURCE

    if non_direction_sequence is None:
        raise NotationException("Nucleotide Sequence must be specified")
    if notation_source is None:
        raise NotationException("Notation Source must be specified")

    from helm_notation.nucleotide_factory import get_instance

    templates = get_instance().get_nucleotide_templates()
    if notation_source not in templates:
        raise NotationException(
            f"Unknown Notation Source [{notation_source}]"
        )
    nucleotides = templates[notation_source]
    key_set = list(nucleotides.keys())

    result: list[Nucleotide] = []
    pos = 0
    while pos < len(non_direction_sequence):
        found = False
        for symbol in key_set:
            if non_direction_sequence[pos:].startswith(symbol):
                found = True
                notation = nucleotides[symbol]
                nuc = Nucleotide(symbol=symbol, notation=notation)
                result.append(nuc)
                pos += len(symbol)
                break
        if not found:
            raise NotationException(
                "Sequence contains unknown nucleotide starting at "
                + non_direction_sequence[pos:]
            )

    return result


def get_nucleotide_templates(
    templates_element: etree._Element,
) -> dict[str, dict[str, str]]:
    """Parse nucleotide templates from an XML element."""
    result: dict[str, dict[str, str]] = {}

    for template_element in templates_element:
        notation_source = template_element.get(_TEMPLATE_NOTATION_SOURCE_ATTRIBUTE, "")
        tmp_map: dict[str, str] = {}
        result[notation_source] = tmp_map
        for nucleotide_element in template_element:
            nucleotide = get_nucleotide_from_element(nucleotide_element)
            tmp_map[nucleotide.symbol] = nucleotide.notation

    return result


def get_nucleotide_from_element(nucleotide_element: etree._Element) -> Nucleotide:
    """Extract a Nucleotide from an XML element."""
    ns = nucleotide_element.nsmap.get(None, "")
    nsprefix = f"{{{ns}}}" if ns else ""

    symbol_el = nucleotide_element.find(f"{nsprefix}{_NUCLEOTIDE_SYMBOL_ELEMENT}")
    notation_el = nucleotide_element.find(
        f"{nsprefix}{_NUCLEOTIDE_MONOMER_NOTATION_ELEMENT}"
    )

    symbol_text = symbol_el.text if symbol_el is not None else ""
    notation_text = notation_el.text if notation_el is not None else ""
    return Nucleotide(symbol=symbol_text or "", notation=notation_text or "")


def get_nucleotide_from_xml(nucleotide_xml: str) -> Nucleotide | None:
    """Parse a Nucleotide from an XML string."""
    if not nucleotide_xml:
        return None
    root = etree.fromstring(nucleotide_xml.encode())
    return get_nucleotide_from_element(root)


def get_nucleotide_templates_xml(
    templates: dict[str, dict[str, str]],
) -> str:
    """Serialize nucleotide templates to XML string."""
    lines: list[str] = []
    lines.append(
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<NUCLEOTIDE_TEMPLATES xsi:schemaLocation="lmr NucleotideTemplateSchema.xsd" '
        'xmlns="lmr" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
    )

    for template_name, nuc_map in templates.items():
        template_element = etree.Element(_TEMPLATE_ELEMENT)
        template_element.set(_TEMPLATE_NOTATION_SOURCE_ATTRIBUTE, template_name)

        for symbol, notation in nuc_map.items():
            nucleotide_element = etree.SubElement(
                template_element, _NUCLEOTIDE_ELEMENT
            )
            symbol_element = etree.SubElement(
                nucleotide_element, _NUCLEOTIDE_SYMBOL_ELEMENT
            )
            symbol_element.text = symbol
            notation_element = etree.SubElement(
                nucleotide_element, _NUCLEOTIDE_MONOMER_NOTATION_ELEMENT
            )
            notation_element.text = notation

        template_string = etree.tostring(
            template_element, pretty_print=True, encoding="unicode"
        )
        lines.append(template_string)

    lines.append("\n</NUCLEOTIDE_TEMPLATES>")
    return "".join(lines)


def get_sirna_notation(
    sense_seq: str,
    anti_sense_seq: str,
    rna_design_type: str = RNA_DESIGN_NONE,
) -> str:
    """Convert nucleotide sequences into HELM notation for siRNA.

    Args:
        sense_seq: 5'-3' nucleotide sequence for sense strand.
        anti_sense_seq: 5'-3' nucleotide sequence for antisense strand.
        rna_design_type: RNA design type (default NONE).

    Returns:
        HELM notation string for siRNA.
    """
    _validate_sirna_design(sense_seq, anti_sense_seq, rna_design_type)

    sb: list[str] = []
    tokens = ["", ""]
    count = 0

    if sense_seq and len(sense_seq) > 0:
        count += 1
        ss_notation = get_notation(sense_seq)
        ss_notation = remove_last_p(ss_notation)
        tokens[0] = f"RNA{count}{{{ss_notation}}}"
        tokens[1] = f"RNA{count}{{ss}}"

    if anti_sense_seq and len(anti_sense_seq) > 0:
        count += 1
        as_notation = get_notation(anti_sense_seq)
        as_notation = remove_last_p(as_notation)
        if len(tokens[0]) > 0:
            tokens[0] += "|"
        tokens[0] += f"RNA{count}{{{as_notation}}}"
        if len(tokens[1]) > 0:
            tokens[1] += "|"
        tokens[1] += f"RNA{count}{{as}}"

    base_pair = _hybridization(sense_seq, anti_sense_seq, rna_design_type)

    if len(tokens[0]) > 0:
        sb.append(tokens[0])
        sb.append("$$")
        if len(base_pair) > 0:
            sb.append(base_pair)
        sb.append("$")
        sb.append(tokens[1])
        sb.append("$")

    return "".join(sb)


def _hybridization(
    sense_seq: str | None,
    anti_sense_seq: str | None,
    rna_design_type: str,
) -> str:
    """Compute hybridization base pair string."""
    base_pair = ""
    if (
        sense_seq
        and len(sense_seq) > 0
        and anti_sense_seq
        and len(anti_sense_seq) > 0
    ):
        analog_seq_ss = get_natural_analog_sequence(sense_seq).replace("T", "U")
        analog_seq_as = get_natural_analog_sequence(anti_sense_seq).replace("T", "U")

        if RNA_DESIGN_NONE.upper() == rna_design_type.upper():
            normal_comp_as = cleanup(get_normal_complement_sequence(analog_seq_as))
            max_match = get_max_match_fragment(analog_seq_ss, normal_comp_as)
            if len(max_match) > 0:
                ss_start = analog_seq_ss.index(max_match)
                normal_comp_start = normal_comp_as.index(max_match)
                as_start = (
                    len(analog_seq_as) - len(max_match) - normal_comp_start
                )

                parts: list[str] = []
                for i in range(len(max_match)):
                    ss_pos = (i + ss_start) * 3 + 2
                    as_pos = (as_start + len(max_match) - 1 - i) * 3 + 2
                    parts.append(
                        f"RNA1,RNA2,{ss_pos}:pair-{as_pos}:pair"
                    )
                base_pair = "|".join(parts)

        elif RNA_DESIGN_TUSCHL_19_PLUS_2.upper() == rna_design_type.upper():
            base_pair = _hybridization_with_length_from_start(
                analog_seq_ss, analog_seq_as, 19
            )
        elif RNA_DESIGN_DICER_27_R.upper() == rna_design_type.upper():
            base_pair = _hybridization_with_length_from_start(
                analog_seq_ss, analog_seq_as, 25
            )
        elif RNA_DESIGN_DICER_27_L.upper() == rna_design_type.upper():
            base_pair = _hybridization_with_length_from_start(
                analog_seq_ss, analog_seq_as, 25
            )
        else:
            raise NotationException("Unknow RNA Design")

    return base_pair


def _hybridization_with_length_from_start(
    sense_analog_seq: str,
    antisense_analog_seq: str,
    length_from_start: int,
) -> str:
    """Compute base pair string with fixed length from start."""
    parts: list[str] = []
    for i in range(length_from_start):
        ss_pos = i * 3 + 2
        as_pos = (length_from_start - 1 - i) * 3 + 2
        ss_char = sense_analog_seq[i]
        as_char = antisense_analog_seq[length_from_start - 1 - i]

        comp = complement_map.get(ss_char, "")
        if comp.upper() == as_char.upper():
            parts.append(f"RNA1,RNA2,{ss_pos}:pair-{as_pos}:pair")
    return "|".join(parts)


def _validate_sirna_design(
    sense_seq: str | None,
    antisense_seq: str | None,
    rna_design_type: str,
) -> bool:
    """Validate siRNA design parameters."""
    if RNA_DESIGN_NONE.upper() == rna_design_type.upper():
        return True

    if rna_design_type not in SUPPORTED_DESIGN_LIST:
        raise NotationException(
            f"Unsupported RNA Design Type '{rna_design_type}'"
        )

    sense_nuc_list = _get_nucleotide_list(sense_seq or "")
    antisense_nuc_list = _get_nucleotide_list(antisense_seq or "")

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


def remove_last_p(simple_rna_notation: str | None) -> str | None:
    """Clean up the last phosphate from simple RNA notation."""
    if simple_rna_notation is not None and simple_rna_notation.endswith("P"):
        return simple_rna_notation[: len(simple_rna_notation) - 1]
    return simple_rna_notation


def get_reverse_sequence(sequence: str) -> str:
    """Return the sequence in reverse direction (e.g. AAGCU => UCGAA)."""
    nuc_list = _get_nucleotide_list(sequence)
    rev_list = _get_reverse_list(nuc_list)
    return "".join(nuc.symbol for nuc in rev_list)


def get_natural_analog_sequence(sequence: str) -> str:
    """Convert modified sequence to natural analog sequence."""
    nuc_list = _get_nucleotide_list(sequence)
    return "".join(nuc.get_natural_analog() for nuc in nuc_list)


def is_in_opposite_direction(sense_seq: str, anti_sense_seq: str) -> bool:
    """Check if antisense is in opposite direction relative to sense.

    Args:
        sense_seq: siRNA sense strand sequence (generally 5'-3').
        anti_sense_seq: siRNA antisense strand sequence.

    Returns:
        True if opposite, False if same direction.
    """
    upper_ss = sense_seq.upper()
    upper_as = anti_sense_seq.upper()
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
    """Return the largest matched fragment between two sequences.

    Replaces T with U before matching.

    Args:
        seq1: Single letter, all upper case nucleotide sequence.
        seq2: Single letter, all upper case nucleotide sequence.
        min_length: Minimum fragment length.

    Returns:
        Largest match fragment string.
    """
    if seq1 is None or seq2 is None:
        raise NotationException("Both sequences must not be null ")

    if seq1 != seq1.upper() or seq2 != seq2.upper():
        raise NotationException(
            "Both sequences must be natural nucleotide sequence in upper case "
        )

    if len(seq1) > len(seq2):
        long_seq = seq1
        short_seq = seq2
    else:
        long_seq = seq2
        short_seq = seq1

    # Replace T with U
    long_seq = long_seq.replace("T", "U")
    short_seq = short_seq.replace("T", "U")

    minimum = MINIMUM_MATCH_FRAGMENT_LENGTH
    if min_length > minimum:
        minimum = min_length

    for length in range(len(short_seq), minimum, -1):
        for i in range(len(short_seq) - length + 1):
            tmp = short_seq[i : i + length]
            if tmp in long_seq:
                return tmp

    return ""
