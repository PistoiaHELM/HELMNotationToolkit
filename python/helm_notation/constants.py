"""Port of org.helm.notation.NotationConstant."""

import os
from pathlib import Path

NOTATION_SOURCE: str = "HELM Notation"
NOTATION_DIRECTORY: str = str(Path.home() / ".helm")

# Notation→SMILES conversion not performed when monomer count exceeds this
MONOMER_COUNT_THRESHOLD: int = 200

# Break long sequences into chunks of this size
MONOMER_GROUP_COUNT_INTERVAL: int = 20
