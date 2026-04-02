#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

PYTHON_BIN="python"
if [[ -x "$ROOT_DIR/.venv/bin/python" ]]; then
  PYTHON_BIN="$ROOT_DIR/.venv/bin/python"
elif command -v python3 >/dev/null 2>&1; then
  PYTHON_BIN="python3"
fi

"$PYTHON_BIN" -m pip install -e ".[dev]"

if ! "$PYTHON_BIN" -c "import pytest" >/dev/null 2>&1; then
  echo "pytest is not installed for: $PYTHON_BIN" >&2
  echo "Install with: $PYTHON_BIN -m pip install -e '.[dev]'" >&2
  exit 1
fi

echo "Running phase 5 CLI tests..."
"$PYTHON_BIN" -m pytest python/tests/test_phase5_cli.py -q

echo "Running full python test suite..."
"$PYTHON_BIN" -m pytest python/tests -q

echo "Phase 5 checks passed."
