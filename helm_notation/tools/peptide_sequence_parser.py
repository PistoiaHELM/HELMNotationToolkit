"""Methods for handling peptide sequences.

Converts between peptide amino acid sequences and HELM simple notation.
"""

from __future__ import annotations

import re
from typing import Optional

from helm_notation.exceptions import MonomerException, NotationException
from helm_notation.model.monomer import PEPTIDE_POLYMER_TYPE

SEQUENCE_DELIMITERS = [".", ":", ",", ";"]


def get_notation(
    sequence: str, delimiter: str | None = None
) -> str:
    """Convert peptide sequence into simple notation for peptide polymer.

    Args:
        sequence: Peptide sequence (with optional delimiter).
        delimiter: Optional delimiter in the sequence.

    Returns:
        Simple notation for peptide polymer.
    """
    sb: list[str] = []
    aa_list = get_amino_acid_list(sequence, delimiter)

    for aa in aa_list:
        token = f"[{aa}]" if len(aa) > 1 else aa
        sb.append(token)
    return ".".join(sb)


def get_amino_acid_list(
    peptide_sequence: str, delimiter: str | None = None
) -> list[str]:
    """Convert peptide sequence into a list of amino acid identifiers.

    Args:
        peptide_sequence: Input peptide sequence.
        delimiter: Optional delimiter in the input sequence.

    Returns:
        List of amino acid identifier strings.
    """
    if peptide_sequence is None:
        raise NotationException("Peptide Sequence must be specified")

    if delimiter is None or len(delimiter) == 0:
        return _get_amino_acid_list_no_delimiter(peptide_sequence)
    else:
        if not _is_valid_delimiter(delimiter):
            raise NotationException(
                f"Invalid sequence delimiter [{delimiter}], "
                f"only the following are supported: [{_get_valid_delimiters()}]"
            )

    blocks = peptide_sequence.split(delimiter)
    result: list[str] = []
    for block in blocks:
        if block:
            tmp_list = _get_amino_acid_list_no_delimiter(block)
            result.extend(tmp_list)
    return result


def _get_amino_acid_list_no_delimiter(peptide_sequence: str) -> list[str]:
    """Convert peptide sequence (without delimiter) into amino acid list."""
    if peptide_sequence is None:
        raise NotationException("Peptide Sequence must be specified")

    clean_seq = cleanup(peptide_sequence)

    from helm_notation.monomer_factory import get_instance

    peptide_map = get_instance().monomer_db.get(PEPTIDE_POLYMER_TYPE, {})
    key_set = list(peptide_map.keys())

    result: list[str] = []
    pos = 0
    while pos < len(clean_seq):
        found = False
        for symbol in key_set:
            if clean_seq[pos:].startswith(symbol):
                found = True
                result.append(symbol)
                pos += len(symbol)
                break
        if not found:
            raise NotationException(
                "Sequence contains unknown amino acid starting at "
                + clean_seq[pos:]
            )

    return result


def cleanup(sequence: str) -> str:
    """Remove white space, and convert all lower case to upper case."""
    result = re.sub(r"\s", "", sequence)
    if result == result.lower():
        result = result.upper()
    return result


def _is_valid_delimiter(delimiter: str) -> bool:
    """Check if delimiter is in the supported list."""
    return delimiter in SEQUENCE_DELIMITERS


def _get_valid_delimiters() -> str:
    """Return string of valid delimiters."""
    return " ".join(SEQUENCE_DELIMITERS)
