"""Command-line interface for HELM Notation Toolkit."""

from __future__ import annotations

import argparse
import json
import sys


def main(args: list[str] | None = None) -> None:
    """Main entry point for the helm-toolkit CLI."""
    parser = argparse.ArgumentParser(
        prog="helm-toolkit",
        description="HELM Notation Toolkit - Parse and manipulate HELM notation using RDKit",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {_get_version()}"
    )
    parser.add_argument(
        "--json", action="store_true", dest="output_json",
        help="Output results as JSON",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Validate command
    validate_parser = subparsers.add_parser(
        "validate", help="Validate a HELM notation string"
    )
    validate_parser.add_argument("notation", help="HELM notation string to validate")

    # SMILES command
    smiles_parser = subparsers.add_parser(
        "smiles", help="Convert HELM notation to canonical SMILES"
    )
    smiles_parser.add_argument("notation", help="HELM notation string")

    # Info command
    info_parser = subparsers.add_parser(
        "info", help="Get molecular info (MW, formula, exact mass) from HELM notation"
    )
    info_parser.add_argument("notation", help="HELM notation string")

    # Canonicalize command
    canon_parser = subparsers.add_parser(
        "canonicalize", help="Canonicalize a HELM notation string"
    )
    canon_parser.add_argument("notation", help="HELM notation string")

    # Monomers command
    monomers_parser = subparsers.add_parser(
        "monomers", help="List available monomers in the database"
    )
    monomers_parser.add_argument(
        "--type", choices=["PEPTIDE", "RNA", "CHEM"],
        help="Filter by polymer type",
    )

    # Extinction coefficient command
    ext_parser = subparsers.add_parser(
        "extinction", help="Calculate extinction coefficient"
    )
    ext_parser.add_argument("notation", help="HELM notation string")

    # Sequence command
    seq_parser = subparsers.add_parser(
        "sequence", help="Extract sequence from HELM notation"
    )
    seq_parser.add_argument("notation", help="HELM notation string")
    seq_parser.add_argument(
        "--type", choices=["peptide", "rna"], default="peptide",
        help="Sequence type to extract (default: peptide)",
    )

    parsed = parser.parse_args(args)

    if parsed.command is None:
        parser.print_help()
        sys.exit(1)

    try:
        if parsed.command == "validate":
            _cmd_validate(parsed.notation, parsed.output_json)
        elif parsed.command == "smiles":
            _cmd_smiles(parsed.notation, parsed.output_json)
        elif parsed.command == "info":
            _cmd_info(parsed.notation, parsed.output_json)
        elif parsed.command == "canonicalize":
            _cmd_canonicalize(parsed.notation, parsed.output_json)
        elif parsed.command == "monomers":
            _cmd_monomers(parsed.type, parsed.output_json)
        elif parsed.command == "extinction":
            _cmd_extinction(parsed.notation, parsed.output_json)
        elif parsed.command == "sequence":
            _cmd_sequence(parsed.notation, parsed.type, parsed.output_json)
    except Exception as e:
        if parsed.output_json:
            print(json.dumps({"error": str(e)}))
        else:
            print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def _get_version() -> str:
    from helm_notation import __version__
    return __version__


def _output(data: dict, as_json: bool) -> None:
    """Output data as JSON or human-readable text."""
    if as_json:
        print(json.dumps(data, indent=2))
    else:
        for key, value in data.items():
            print(f"{key}: {value}")


def _cmd_validate(notation: str, as_json: bool) -> None:
    from helm_notation.tools.complex_notation_parser import validate_complex_notation

    try:
        validate_complex_notation(notation)
        if as_json:
            print(json.dumps({"valid": True, "notation": notation}))
        else:
            print("Valid HELM notation")
    except Exception as e:
        if as_json:
            print(json.dumps({"valid": False, "error": str(e)}))
        else:
            print(f"Invalid: {e}")
        sys.exit(1)


def _cmd_smiles(notation: str, as_json: bool) -> None:
    from helm_notation.tools.complex_notation_parser import get_complex_polymer_smiles

    smiles = get_complex_polymer_smiles(notation)
    if as_json:
        print(json.dumps({"smiles": smiles, "notation": notation}))
    else:
        print(smiles)


def _cmd_info(notation: str, as_json: bool) -> None:
    from helm_notation.tools.complex_notation_parser import get_molecule_info

    mi = get_molecule_info(notation)
    data = {
        "Molecular Weight": f"{mi.molecular_weight:.2f}",
        "Molecular Formula": mi.molecular_formula,
        "Exact Mass": f"{mi.exact_mass:.2f}",
    }
    if as_json:
        print(json.dumps({
            "molecular_weight": mi.molecular_weight,
            "molecular_formula": mi.molecular_formula,
            "exact_mass": mi.exact_mass,
            "notation": notation,
        }))
    else:
        _output(data, False)


def _cmd_canonicalize(notation: str, as_json: bool) -> None:
    from helm_notation.tools.complex_notation_parser import get_canonical_notation

    canonical = get_canonical_notation(notation)
    if as_json:
        print(json.dumps({"canonical": canonical, "original": notation}))
    else:
        print(canonical)


def _cmd_monomers(polymer_type: str | None, as_json: bool) -> None:
    from helm_notation.monomer_factory import get_instance

    factory = get_instance()
    db = factory.get_monomer_db()

    if polymer_type:
        types_to_show = [polymer_type]
    else:
        types_to_show = sorted(db.keys())

    if as_json:
        result = {}
        for pt in types_to_show:
            if pt in db:
                result[pt] = sorted(db[pt].keys())
        print(json.dumps(result, indent=2))
    else:
        for pt in types_to_show:
            if pt in db:
                monomers = sorted(db[pt].keys())
                print(f"\n{pt} ({len(monomers)} monomers):")
                print(f"  {', '.join(monomers)}")


def _cmd_extinction(notation: str, as_json: bool) -> None:
    from helm_notation.tools.extinction_coefficient_calculator import get_instance

    calc = get_instance()
    result = calc.calculate_from_complex_notation(notation)
    unit = calc.get_default_unit()

    if as_json:
        print(json.dumps({
            "extinction_coefficient": result,
            "unit": unit,
            "notation": notation,
        }))
    else:
        print(f"Extinction Coefficient: {result:.2f} {unit}")


def _cmd_sequence(notation: str, seq_type: str, as_json: bool) -> None:
    from helm_notation.tools.complex_notation_parser import get_all_node_string
    from helm_notation.tools import simple_notation_parser as snp

    # Extract the first node of the matching type
    node_str = get_all_node_string(notation)
    target_prefix = "PEPTIDE" if seq_type == "peptide" else "RNA"

    for node in node_str.split("|"):
        node = node.strip()
        if not node:
            continue
        label, simple = node.split("{", 1)
        simple = simple.rstrip("}")
        if label.startswith(target_prefix):
            if seq_type == "peptide":
                sequence = snp.get_peptide_sequence(simple)
            else:
                sequence = snp.get_nucleotide_sequence(simple)

            if as_json:
                print(json.dumps({
                    "sequence": sequence,
                    "type": seq_type,
                    "polymer": label,
                    "notation": notation,
                }))
            else:
                print(sequence)
            return

    msg = f"No {target_prefix} polymer found in notation"
    if as_json:
        print(json.dumps({"error": msg}))
    else:
        print(msg, file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
