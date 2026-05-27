#!/usr/bin/env bash
set -euo pipefail
python -m venv .venv-alphagenome
source .venv-alphagenome/bin/activate
python -m pip install --upgrade pip
pip install alphagenome alphagenome-research huggingface_hub pandas numpy matplotlib
