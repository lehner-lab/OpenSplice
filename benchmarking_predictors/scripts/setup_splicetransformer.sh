#!/usr/bin/env bash
set -euo pipefail
python -m venv .venv-splicetransformer
source .venv-splicetransformer/bin/activate
python -m pip install --upgrade pip
pip install -r requirements_splicetransformer.txt
