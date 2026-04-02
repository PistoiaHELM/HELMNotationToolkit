#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
OUTPUT_FILE="${1:-$ROOT_DIR/migration/fixtures/java_baseline.json}"

cd "$ROOT_DIR"

ant compile-demo >/dev/null

java -cp "build/classes:build/demo:lib/*" \
  org.helm.notation.demo.tools.BaselineGenerator \
  "$OUTPUT_FILE"

echo "Baseline generated at: $OUTPUT_FILE"
