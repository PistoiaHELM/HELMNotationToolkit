"""Constants for HELM notation."""

import os
from pathlib import Path

NOTATION_SOURCE = "HELM Notation"
NOTATION_DIRECTORY = Path.home() / ".helm"

# Notation to SMILES conversion will not be performed when total monomer
# count is above this threshold
MONOMER_COUNT_THRESHOLD = 200

# Break long sequences to small chunks using this value
MONOMER_GROUP_COUNT_INTERVAL = 20
