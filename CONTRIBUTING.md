# Developer Contributing Guide

This guide covers how to set up the HELM Notation Toolkit Python project for local development.

## Prerequisites

- **Python 3.10+** 
- **uv** (package manager) - Install via: `brew install uv` or `curl -LsSf https://astral.sh/uv/install.sh | sh`
- **RDKit** - Auto-installed as dependency

## Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/PistoiaHELM/HELMNotiaionToolkit.git
cd HELMNotiaionToolkit

# 2. Create and activate virtual environment
uv venv
source .venv/bin/activate

# 3. Install package in editable mode
make install-dev

# 4. Run tests to verify setup
make test
```

## Development Setup

### Using uv (Recommended)

```bash
# Create virtual environment
uv venv

# Activate (on macOS/Linux)
source .venv/bin/activate

# Install with dev dependencies
make install-dev
```

### Using Make Commands

```bash
make test          # Run all tests
make test-s        # Run tests with stdout output
make test-coverage # Run tests with coverage report
make build         # Build package with uv
make install       # Install package in editable mode
make install-dev   # Install with dev dependencies
make clean-build   # Remove build artifacts
make help          # Show available commands
```

### Verify Installation

```bash
# Check CLI is available
helm-notation --help

# Run all tests
python -m pytest python/tests/ -v

# Run specific test file
python -m pytest python/tests/test_nucleotide_converter.py -v -s
```

## Project Structure

```
HELMNotiaionToolkit/
├── python/
│   ├── helm_notation/          # Main package
│   │   ├── chem/               # Chemistry engine (RDKit)
│   │   ├── peptide/            # Peptide parsing
│   │   ├── tools/              # Ported Java tools
│   │   └── resources/          # XML templates
│   └── tests/                  # Test suite
├── migration/
│   ├── fixtures/               # Baseline test fixtures
│   └── README.md               # Migration documentation
└── dist/                       # Built distributions
```

## Running Tests

### Run All Tests
```bash
python -m pytest python/tests/ -v
```

### Run with Output (for debugging)
```bash
python -m pytest python/tests/ -v -s
```

### Run Specific Test Class
```bash
python -m pytest python/tests/test_nucleotide_converter.py::TestGetNucleotideSequenceFromSimpleRNANotation -v -s
```

### Run Tests Matching Pattern
```bash
python -m pytest python/tests/ -v -k "nucleotide"
```

## Building the Package

### Build with uv
```bash
# Clean previous builds
rm -rf dist/ build/ python/*.egg-info/

# Build sdist and wheel
uv build

# View built files
ls -la dist/
```

### Install Built Package
```bash
uv pip install dist/helm_notation_toolkit_py-0.1.0-py3-none-any.whl
```

## Code Style

This project follows Python conventions:

- **Type hints** for all function signatures
- **Docstrings** for public APIs
- **pytest** for testing

### Running Linting (if configured)

```bash
# Check for issues (flake8, mypy, etc.)
ruff check python/
```

## Common Tasks

### Adding a New Ported Java Class

1. Create module in `python/helm_notation/tools/`
2. Add tests in `python/tests/test_<module>.py`
3. Run tests to verify parity
4. Update baseline fixtures if needed

### Updating Baseline Fixtures

Edit `migration/fixtures/java_baseline.json` and add test cases:

```json
{
  "id": "your_test_case",
  "operation": "operation_name",
  "input": "test input",
  "output": "expected output",
  "error": null
}
```

### Running CLI Commands

```bash
# Sequence to notation
helm-notation seq2notation AUGCUU

# Simple notation to SMILES
helm-notation simple-smiles PEPTIDE A.G.C

# Complex notation to SMILES
helm-notation complex-smiles "PEPTIDE1{A.G.C}$$$"

# Complex notation to molecule info
helm-notation notation2property "PEPTIDE1{A.G.C}$$$"
```

## Troubleshooting

### RDKit Installation Issues

```bash
# If RDKit fails to install, try
uv pip install rdkit --no-cache
```

### Test Failures

- Check Python version: `python --version` (should be 3.10+)
- Verify rdkit installed: `python -c "from rdkit import Chem; print('OK')"`
- Reinstall: `uv pip install -e ".[dev]"`

### Import Errors

```bash
# Ensure package is installed
uv pip install -e "."

# Check PYTHONPATH
echo $PYTHONPATH
```

## Additional Resources

- [Migration README](migration/README.md) - Detailed phase documentation
- [RDKit Documentation](https://rdkit.org/docs/)
- [HELM Notation Specification](https://www.helm-notation.org/)