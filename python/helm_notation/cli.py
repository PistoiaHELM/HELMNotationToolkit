from __future__ import annotations

import argparse
from typing import Iterable, List

from helm_notation import ComplexNotationParser, PeptideStructureParser, SimpleNotationParser
from helm_notation.nucleotide_factory import NucleotideFactory


def _sequence_to_notation(sequence: str) -> str:
    seq = sequence.strip().upper().replace("T", "U")
    if not seq:
        raise ValueError("Sequence is empty")

    templates = NucleotideFactory.get_instance().get_nucleotide_templates().get("HELM Notation", {})
    parts: List[str] = []
    for base in seq:
        notation = templates.get(base)
        if notation is None:
            raise ValueError(f"Unsupported base in sequence: {base}")
        parts.append(notation)
    return ".".join(parts)


def _notation_to_property(notation: str) -> str:
    info = ComplexNotationParser.get_molecule_info(notation)
    return f"MW={info.molecular_weight:.8f};EXACT={info.exact_mass:.9f};FORMULA={info.molecular_formula}"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="helm-notation",
        description="HELM Notation Toolkit (Python/RDKit)",
    )

    sub = parser.add_subparsers(dest="command", required=True)

    seq2notation = sub.add_parser("seq2notation", help="Convert RNA/DNA sequence to HELM simple notation")
    seq2notation.add_argument("sequence", help="Sequence string, e.g. AUGCUU or ATGCTT")

    notation2property = sub.add_parser(
        "notation2property",
        help="Compute molecule info (MW/EXACT/FORMULA) from HELM complex notation",
    )
    notation2property.add_argument("notation", help="HELM complex notation string")

    simple_smiles = sub.add_parser("simple-smiles", help="Generate CXSMILES from simple notation")
    simple_smiles.add_argument("polymer_type", choices=["PEPTIDE", "RNA", "CHEM"], help="Polymer type")
    simple_smiles.add_argument("notation", help="Simple notation string")

    complex_smiles = sub.add_parser("complex-smiles", help="Generate canonical smiles from complex notation")
    complex_smiles.add_argument("notation", help="Complex HELM notation string")

    peptide_notation = sub.add_parser(
        "peptide-from-smiles",
        help="Infer peptide HELM notation from peptide smiles",
    )
    peptide_notation.add_argument("smiles", help="Peptide smiles string")

    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)

    if args.command == "seq2notation":
        print(_sequence_to_notation(args.sequence))
        return 0

    if args.command == "notation2property":
        print(_notation_to_property(args.notation))
        return 0

    if args.command == "simple-smiles":
        print(SimpleNotationParser.get_simple_polymer_smiles(args.notation, args.polymer_type))
        return 0

    if args.command == "complex-smiles":
        print(ComplexNotationParser.get_complex_polymer_smiles(args.notation))
        return 0

    if args.command == "peptide-from-smiles":
        print(PeptideStructureParser.get_instance().smiles2notation(args.smiles))
        return 0

    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
