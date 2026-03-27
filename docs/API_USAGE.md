# HELM Notation Toolkit — Python API Usage Guide

A Python library for parsing, validating, and manipulating HELM (Hierarchical Editing Language for Macromolecules) notation. Uses **RDKit** for chemistry and **lxml** for XML processing.

## Installation

```bash
pip install rdkit lxml
```

## Quick Start

```python
from helm_notation.monomer_factory import get_instance as get_monomer_factory
from helm_notation.nucleotide_factory import get_instance as get_nucleotide_factory

# Initialize factories (required before any operations)
mf = get_monomer_factory()
nf = get_nucleotide_factory()
```

---

## Table of Contents

1. [Initialization](#1-initialization)
2. [Monomer Database](#2-monomer-database)
3. [Structure Parsing (RDKit)](#3-structure-parsing-rdkit)
4. [Simple Notation Parsing](#4-simple-notation-parsing)
5. [Complex Notation Parsing](#5-complex-notation-parsing)
6. [Nucleotide Sequences](#6-nucleotide-sequences)
7. [Peptide Sequences](#7-peptide-sequences)
8. [Extinction Coefficient Calculator](#8-extinction-coefficient-calculator)
9. [XHELM Import/Export](#9-xhelm-importexport)
10. [Molfile Encoding](#10-molfile-encoding)

---

## 1. Initialization

### MonomerFactory

Loads the monomer database from `MonomerDBGZEncoded.xml`. Must be initialized before using any parser.

```python
from helm_notation.monomer_factory import get_instance, finalize

# Initialize (loads ~193 monomers: 119 PEPTIDE, 63 RNA, 11 CHEM)
factory = get_instance()

# Access the monomer store
store = factory.monomer_store

# Get the raw monomer database
monomer_db = factory.get_monomer_db()
for polymer_type, monomers in monomer_db.items():
    print(f"{polymer_type}: {len(monomers)} monomers")

# Get attachment database
attachment_db = factory.get_attachment_db()

# Clean up when done
finalize()
```

### NucleotideFactory

Loads nucleotide templates from `NucleotideTemplates.xml`. Required for modified nucleotide operations.

```python
from helm_notation.nucleotide_factory import get_instance as get_nf

nf = get_nf()
templates = nf.get_nucleotide_templates()
# templates: {"A": {"symbol": "A", "notation": "R(A)P"}, ...}

# Reverse lookup: notation → symbol
reverse_map = nf.get_reverse_nucleotide_template_map()
```

### MonomerStore

Access and query monomers programmatically.

```python
from helm_notation.model.monomer import PEPTIDE_POLYMER_TYPE, NUCLEIC_ACID_POLYMER_TYPE

store = factory.monomer_store

# Look up a monomer
ala = store.get_monomer(PEPTIDE_POLYMER_TYPE, "A")
print(f"Alanine SMILES: {ala.can_smiles}")
print(f"Natural analog: {ala.natural_analog}")

# Check existence
assert store.has_monomer(PEPTIDE_POLYMER_TYPE, "G")

# Get all monomers for a polymer type
rna_monomers = store.get_monomers(NUCLEIC_ACID_POLYMER_TYPE)

# Look up by SMILES
monomer = store.get_monomer_by_smiles("C(C(=O)O)N |$;;;_R2;;_R1$|")
```

---

## 2. Monomer Database

### Monomer Model

```python
from helm_notation.model.monomer import Monomer

monomer = store.get_monomer("PEPTIDE", "W")
print(monomer.alternate_id)     # "W"
print(monomer.name)             # "Tryptophan"
print(monomer.can_smiles)       # canonical SMILES
print(monomer.natural_analog)   # "W"
print(monomer.polymer_type)     # "PEPTIDE"
print(monomer.monomer_type)     # "Backbone"
print(monomer.attachments)      # list of Attachment objects
```

### Attachment Model

```python
from helm_notation.model.attachment import Attachment

for att in monomer.attachments:
    print(f"  {att.alternate_id}: {att.label} - cap: {att.cap_group_smiles}")
```

---

## 3. Structure Parsing (RDKit)

All chemistry operations go through `structure_parser`.

### SMILES Validation

```python
from helm_notation.tools import structure_parser

# Validate a SMILES string
assert structure_parser.validate_smiles("CCO") is True
assert structure_parser.validate_smiles("invalid") is False
```

### Molecule Info (MW, Formula, Exact Mass)

```python
info = structure_parser.get_molecule_info("CC(=O)OC1=CC=CC=C1C(=O)O")
print(f"MW: {info.molecular_weight}")
print(f"Formula: {info.molecular_formula}")
print(f"Exact mass: {info.exact_mass}")
```

### Parse Molecules

```python
from rdkit import Chem

# From SMILES
mol = structure_parser.get_molecule("CCO")

# Canonical SMILES
unique = structure_parser.get_unique_smiles("OCC")
print(unique)  # "CCO"

# From RDKit Mol
smiles = structure_parser.get_unique_smiles_from_mol(mol)
```

### Extended SMILES (CXSMILES) with R-groups

```python
# Monomer SMILES with R-group labels
ext_smiles = "[*]OC1CC(O[*])C(CO[*])O1 |$_R1;;;;;_R2;;;;_R3;;$|"

# Get simplified extended SMILES
simple = structure_parser.get_simple_extended_smiles(ext_smiles)

# Get unique (canonical) extended SMILES
unique = structure_parser.get_unique_extended_smiles(ext_smiles)

# List R-groups
rgroups = structure_parser.get_rgroups_from_extended_smiles(ext_smiles)
print(rgroups)  # ["R1", "R2", "R3"]
```

### Molecule Merging at R-groups

```python
# Merge two molecules at their R-group attachment points
mol1 = structure_parser.get_molecule("[1*]OCC")
mol2 = structure_parser.get_molecule("[2*]NCC")

r1_idx = structure_parser.find_rgroup_atom(mol1, 1)
r2_idx = structure_parser.find_rgroup_atom(mol2, 2)

merged = structure_parser.merge(mol1, r1_idx, mol2, r2_idx)
print(Chem.MolToSmiles(merged))
```

### Decompose Mixtures

```python
components = structure_parser.decompose_mixture(mixture_smiles)
# Returns list of individual component SMILES
```

---

## 4. Simple Notation Parsing

Handles single-polymer HELM notation (e.g., `A.G.G.W.E` for peptides, `R(A)P.R(C)P.R(G)` for RNA).

### Validate Notation

```python
from helm_notation.tools import simple_notation_parser as snp

# Validate peptide notation
assert snp.validate_simple_notation_for_peptide("A.G.G.W.E")

# Validate RNA notation
assert snp.validate_simple_notation_for_rna("R(A)P.R(C)P.R(G)P.[mR](A)")

# Validate CHEM notation
assert snp.validate_simple_notation_for_chem("PEG2")

# Generic validation
assert snp.validate_simple_notation("A.G.W", "PEPTIDE")
```

### Get SMILES from Notation

```python
smiles = snp.get_simple_polymer_smiles("A.G.G.W.E.E", "PEPTIDE")
print(f"Peptide SMILES: {smiles}")

rna_smiles = snp.get_simple_polymer_smiles(
    "R(A)P.R(C)P.R(G)P.[mR](A)", "RNA"
)
```

### Get Polymer Structure (RDKit Mol)

```python
from helm_notation.model.rgroup_structure import RgroupStructure

structure = snp.get_simple_polymer_structure("A.G.W", "PEPTIDE")
# structure.molecule is an RDKit Mol
# structure.rgroup_map maps R-group IDs to atom indices
```

### Monomer Operations

```python
# Get list of Monomer objects
monomers = snp.get_monomer_list("A.G.W.E", "PEPTIDE")
for m in monomers:
    print(f"{m.alternate_id}: {m.name}")

# Get monomer IDs only
ids = snp.get_monomer_id_list("A.G.W.E", "PEPTIDE")
print(ids)  # ["A", "G", "W", "E"]

# Count monomers
count = snp.get_monomer_count_for_peptide("A.G.W.E")
print(count)  # 4

# Look up a single monomer
ala = snp.get_monomer("A", "PEPTIDE")
```

### Nucleotide Operations

```python
# Get nucleotide list from RNA notation
nucleotides = snp.get_nucleotide_list("R(A)P.R(C)P.R(G)P.[mR](A)")
for nuc in nucleotides:
    print(f"{nuc.symbol}: {nuc.notation}")

# Get nucleotide sequence (natural analogs)
seq = snp.get_nucleotide_sequence(polymer_notation="R(A)P.R(C)P.R(G)P.[mR](A)")
print(seq)  # "ACGA"

# Get modified nucleotide sequence
mod_seq = snp.get_modified_nucleotide_sequence(
    polymer_notation="R(A)P.R(C)P.R(G)P.[mR](A)"
)
```

### Peptide Sequence

```python
seq = snp.get_peptide_sequence("A.G.G.W.E.E")
print(seq)  # "AGGWEE"

mod_seq = snp.get_modified_peptide_sequence("A.G.[meA].W")
```

### Canonical Notation

```python
canonical = snp.get_simple_canonical_notation("G.A.W", "PEPTIDE")
```

### Convert to Complex Notation

```python
complex_notation = snp.get_complex_notation_for_peptide("A.G.W.E")
print(complex_notation)  # "PEPTIDE1{A.G.W.E}$$$$"

complex_rna = snp.get_complex_notation_for_rna("R(A)P.R(C)P.R(G)")
```

### Replace Monomers

```python
new_notation = snp.replace_monomer("A.G.W.E", "PEPTIDE", "A", "meA")
```

### Molecule Info for Polymer

```python
info = snp.get_molecule_info_for_polymer("A.G.W.E", "PEPTIDE")
print(f"MW: {info.molecular_weight}")
print(f"Formula: {info.molecular_formula}")
```

---

## 5. Complex Notation Parsing

Handles full HELM notation with multiple polymers, connections, base pairs, and annotations.

### HELM Notation Format

```
PEPTIDE1{A.G.W}|RNA1{R(A)P.R(C)P}$PEPTIDE1,RNA1,3:R2-1:R1$$$
 └─ polymers ─┘                     └── connections ──┘
```

### Parse Complex Notation

```python
from helm_notation.tools import complex_notation_parser as cnp

notation = "PEPTIDE1{A.G.G.W.E.E.E.E.E.W}$$$$"

# Validate
assert cnp.validate_complex_notation(notation)

# Parse into ComplexPolymer object
cp = cnp.parse(notation)
print(f"Polymer nodes: {len(cp.polymer_node_list)}")
print(f"Edges: {len(cp.polymer_edge_list)}")
```

### Get SMILES

```python
smiles = cnp.get_complex_polymer_smiles(
    "PEPTIDE1{A.G.G.W.E.E.E.E.E.W}$$$$"
)

canonical_smiles = cnp.get_complex_polymer_canonical_smiles(
    "PEPTIDE1{A.G.G.W.E.E.E.E.E.W}$$$$"
)
```

### Molecule Info

```python
info = cnp.get_molecule_info("PEPTIDE1{A.G.G.W.E.E.E.E.E.W}$$$$")
print(f"MW: {info.molecular_weight}")
print(f"Formula: {info.molecular_formula}")
```

### Extract Notation Components

```python
notation = "PEPTIDE1{A.G.W}|RNA1{R(A)P.R(C)P}$PEPTIDE1,RNA1,3:R2-1:R1$$V1$"

nodes = cnp.get_all_node_string(notation)       # "PEPTIDE1{A.G.W}|RNA1{R(A)P.R(C)P}"
edges = cnp.get_all_edge_string(notation)        # "PEPTIDE1,RNA1,3:R2-1:R1"
pairs = cnp.get_all_base_pair_string(notation)   # ""
labels = cnp.get_all_node_label_string(notation)  # "V1"
```

### Polymer Node Lists

```python
node_list = cnp.get_polymer_node_list(nodes)
for node in node_list:
    print(f"{node.label}: {node.type} - {node.notation}")

# RNA-specific nodes
rna_nodes = cnp.get_rna_polymer_node_list(notation)
```

### Canonical Notation

```python
canonical = cnp.get_canonical_notation(
    "PEPTIDE1{A.G.W}|PEPTIDE2{E.E.E}$$$$"
)
```

### Combine Notations

```python
combined = cnp.get_combined_complex_notation(
    "PEPTIDE1{A.G.W}$$$$",
    "RNA1{R(A)P.R(C)P}$$$$"
)
```

### siRNA Operations

```python
# Hybridize RNA strands
hybridized = cnp.hybridize("RNA1{R(A)P.R(C)P}|RNA2{R(G)P.R(U)P}$$$$")

# Get formatted siRNA sequences
sequences = cnp.get_formated_sirna_sequences(hybridized)
```

### Decompose and Standardize

```python
# Decompose into independent components
components = cnp.decompose("PEPTIDE1{A.G}|PEPTIDE2{W.E}$$$$")

# Standardize notation
std = cnp.standardize("PEPTIDE1{A.G.W}$$$$")
```

### Monomer Count

```python
count = cnp.get_total_monomer_count("PEPTIDE1{A.G.W}|RNA1{R(A)P.R(C)P}$$$$")
```

---

## 6. Nucleotide Sequences

Convert between nucleotide sequences and HELM notation.

```python
from helm_notation.tools import nucleotide_sequence_parser as nsp

# Convert sequence to HELM RNA notation
notation = nsp.get_notation("ACGT")
print(notation)  # "R(A)P.R(C)P.R(G)P.R(T)"

# Get complement sequences
rev_comp = nsp.get_reverse_complement_sequence("ACGT")
norm_comp = nsp.get_normal_complement_sequence("ACGT")

# Check direction
is_normal = nsp.is_normal_direction("5'-ACGT-3'")

# Clean up sequence string
cleaned = nsp.cleanup("  5'-ACGT-3'  ")

# Get natural analog sequence
analog = nsp.get_natural_analog_sequence("ACGmT")

# Generate siRNA notation from sense/antisense
sirna = nsp.get_sirna_notation("ACGUACGU", "ACGUACGU")

# Find maximum matching fragment between sequences
match = nsp.get_max_match_fragment("ACGTACGT", "CGTACGTA")
```

---

## 7. Peptide Sequences

Convert between peptide sequences and HELM notation.

```python
from helm_notation.tools import peptide_sequence_parser as psp

# Convert amino acid sequence to HELM peptide notation
notation = psp.get_notation("AGGWEE")
print(notation)  # "A.G.G.W.E.E"

# Parse amino acids from sequence
aa_list = psp.get_amino_acid_list("AGGWEE")
print(aa_list)  # ["A", "G", "G", "W", "E", "E"]

# With delimiter for modified residues
aa_list = psp.get_amino_acid_list("A-G-[meA]-W", delimiter="-")

# Clean up sequence
cleaned = psp.cleanup("  AGGWEE  ")
```

---

## 8. Extinction Coefficient Calculator

Calculate extinction coefficients at 280nm for peptides and nucleotides.

```python
from helm_notation.tools import extinction_coefficient_calculator as ecc

calc = ecc.get_instance()

# From amino acid sequence
ec = calc.calculate_from_amino_acid_sequence("AGGWEEEEEEEEEEEEEEEEE")
print(f"ε₂₈₀ = {ec}")  # 5500.0

# From peptide polymer notation
ec = calc.calculate_from_peptide_polymer_notation("A.G.G.W.E.E.E.E.E.W")
print(f"ε₂₈₀ = {ec}")  # 11000.0

# From nucleotide sequence
ec = calc.calculate_from_nucleotide_sequence("ACGTACGT")

# From RNA polymer notation
ec = calc.calculate_from_rna_polymer_notation("P.R(A)P.R(C)P.R(G)P.[mR](A)")

# From modified nucleotide sequence (requires NucleotideFactory)
ec = calc.calculate_from_modified_nucleotide_sequence("ACGmTACmGT")

# From full HELM complex notation
ec = calc.calculate_from_complex_notation(
    "PEPTIDE1{A.G.G.W.E.E.E.E.E.W}$$$$",
    unit_type=ecc.PEPTIDE_UNIT_TYPE
)

# Unit types
print(ecc.PEPTIDE_UNIT_TYPE)   # 1 — returns raw molar extinction coefficient
print(ecc.DEFAULT_UNIT_TYPE)   # 0 — returns mg/mL extinction coefficient

# Get unit strings
print(calc.get_peptide_unit())  # "M⁻¹cm⁻¹" or similar
print(calc.get_default_unit())
```

---

## 9. XHELM Import/Export

Parse and export XHELM XML format, which embeds monomer definitions alongside HELM notation.

### Parse XHELM

```python
from helm_notation.tools import xhelm_parser
from lxml import etree

# Parse from XML string
xhelm_content = open("molecule.xhelm").read()
root = xhelm_parser.parse_xhelm_string(xhelm_content)

# Extract HELM notation string
helm_notation = xhelm_parser.get_complex_notation_string(root)
print(helm_notation)

# Extract embedded monomer store (includes ad-hoc monomers)
store = xhelm_parser.get_monomer_store(root)
```

### Export XHELM

```python
from helm_notation.tools import xhelm_exporter

store = get_monomer_factory().monomer_store
xhelm_xml = xhelm_exporter.write_xhelm(
    "PEPTIDE1{A.G.W}$$$$",
    store
)
# Returns XML string with HELM notation + monomer definitions
```

### Round-Trip

```python
# Parse XHELM → get notation → export back to XHELM
root = xhelm_parser.parse_xhelm_string(xhelm_xml)
notation = xhelm_parser.get_complex_notation_string(root)
store = xhelm_parser.get_monomer_store(root)
new_xhelm = xhelm_exporter.write_xhelm(notation, store)
```

---

## 10. Molfile Encoding

Encode/decode molfiles using gzip + base64 (used in MonomerDB XML).

```python
from helm_notation.tools import molfile_encoder

# Encode a molfile string
molfile = """
  Mrv1718 ...
  3  2  0  0  0  ...
"""
encoded = molfile_encoder.encode(molfile)

# Decode back
decoded = molfile_encoder.decode(encoded)
assert decoded.strip() == molfile.strip()
```

---

## Complete Example: Analyze a HELM Molecule

```python
from helm_notation.monomer_factory import get_instance as get_mf
from helm_notation.nucleotide_factory import get_instance as get_nf
from helm_notation.tools import complex_notation_parser as cnp
from helm_notation.tools import extinction_coefficient_calculator as ecc

# Initialize
mf = get_mf()
nf = get_nf()

# Define a HELM notation with two peptide chains
notation = (
    "PEPTIDE1{A.G.G.W.E.E.E.E.E.W}"
    "|PEPTIDE2{A.G.G.W.E.Y.E.E.E.E.W}"
    "$$$$"
)

# Validate
assert cnp.validate_complex_notation(notation)

# Parse
cp = cnp.parse(notation)
print(f"Number of chains: {len(cp.polymer_node_list)}")

# Molecule info
info = cnp.get_molecule_info(notation)
print(f"Molecular weight: {info.molecular_weight:.2f}")
print(f"Molecular formula: {info.molecular_formula}")

# Extinction coefficient
calc = ecc.get_instance()
ec = calc.calculate_from_complex_notation(notation, ecc.PEPTIDE_UNIT_TYPE)
print(f"ε₂₈₀ = {ec:.1f} M⁻¹cm⁻¹")

# Canonical form
canonical = cnp.get_canonical_notation(notation)
print(f"Canonical: {canonical}")

# Total monomer count
count = cnp.get_total_monomer_count(notation)
print(f"Total monomers: {count}")
```

---

## Complete Example: RNA siRNA Design

```python
from helm_notation.monomer_factory import get_instance as get_mf
from helm_notation.nucleotide_factory import get_instance as get_nf
from helm_notation.tools import nucleotide_sequence_parser as nsp
from helm_notation.tools import complex_notation_parser as cnp

# Initialize
mf = get_mf()
nf = get_nf()

# Define sense and antisense sequences
sense = "ACGUACGU"
antisense = nsp.get_reverse_complement_sequence(sense)
print(f"Sense:      {sense}")
print(f"Antisense:  {antisense}")

# Generate siRNA HELM notation
sirna = nsp.get_sirna_notation(sense, antisense)
print(f"siRNA HELM: {sirna}")

# Validate and analyze
assert cnp.validate_complex_notation(sirna)
info = cnp.get_molecule_info(sirna)
print(f"MW: {info.molecular_weight:.2f}")
```

---

## Error Handling

```python
from helm_notation.exceptions import (
    HelmError,
    MonomerException,
    NotationException,
    StructureException,
    NucleotideLoadingException,
)

try:
    cnp.validate_complex_notation("INVALID{notation}$$$$")
except NotationException as e:
    print(f"Invalid notation: {e}")

try:
    snp.get_monomer("UNKNOWN", "PEPTIDE")
except MonomerException as e:
    print(f"Monomer not found: {e}")
```
