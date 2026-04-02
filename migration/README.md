# Java -> Python Migration Workspace

This folder tracks staged migration from Java/MarvinBeans to Python/RDKit.

## Phase 0 (complete)
- Add Java baseline exporter for deterministic parity fixtures.
- Add script to regenerate fixture JSON from current Java implementation.
- Baseline output file: `migration/fixtures/java_baseline.json`.

### Regenerate baseline

```bash
./migration/scripts/generate_java_baseline.sh
```

Or write to custom path:

```bash
./migration/scripts/generate_java_baseline.sh /tmp/java_baseline.json
```

## Phase 1 (complete) - StructureParser Parity
- Added missing methods to `chem/engine.py`:
  - `get_rgroups_from_extended_smiles()` - Extract R-group labels from extended SMILES
  - `process_molecule_info()` - Combine MoleculeInfo from chunks and caps

## Phase 2 (complete) - EdgeParser Port
- Ported `EdgeParser.java` → `tools/edge_parser.py`
- Implemented functions:
  - `parse()` - Parse edge notation into PolymerEdge
  - `canonicalize()` - Canonicalize edge based on sorted node IDs
  - `sort()` - Sort canonical polymer edges

## Phase 3 (complete) - Baseline Fixture Expansion
- Expanded `migration/fixtures/java_baseline.json` with 13 new test cases:
  - Edge parsing (standard, pair, generic connections)
  - SMILES validation (valid/invalid)
  - Extended SMILES operations
  - Mixture decomposition
  - Peptide conversion
  - Complex notation validity

## Phase 4 (complete) - Test Fixes
- Updated RNA parser test to check robust traits rather than exact match
- RDKit and ChemAxon produce different canonical RNA forms (expected)

## Phase 5 (complete)
- Added Python CLI entrypoint (`helm-notation`) via `pyproject.toml` console script.
- Added command parity for Java command-line use-cases:
	- `seq2notation`
	- `notation2property`
- Added additional utility commands:
	- `simple-smiles`
	- `complex-smiles`
	- `peptide-from-smiles`
- Added validation script: `migration/scripts/run_phase5_checks.sh`.

## Phase 6 (complete) - Utility Classes
- Ported `DeepCopy.java` → `tools/deep_copy.py`
- Ported `PermutationAndExpansion.java` → `tools/permutation_and_expansion.py`

### CLI examples

```bash
# sequence -> HELM simple notation
helm-notation seq2notation AUGCUU

# simple notation -> CXSMILES
helm-notation simple-smiles PEPTIDE A.G.C

# complex notation -> canonical smiles
helm-notation complex-smiles "PEPTIDE1{A.G.C}$$$"

# complex notation -> molecular properties
helm-notation notation2property "PEPTIDE1{A.G.C}$$$"

# peptide structure -> notation
helm-notation peptide-from-smiles "C[C@H](N[*])C(=O)NCC([*])=O |$;;;_R1;;;;_R2;$|"
```

### Run checks

```bash
./migration/scripts/run_phase5_checks.sh
```

---

## Conversion Status

| Category | Count |
|----------|-------|
| Java Classes Total | 42 |
| Converted | 31 (~74%) |
| Unused (Not Needed) | 5 |

### Run Tests

```bash
pip install -e .
python -m pytest python/tests/ -v
```