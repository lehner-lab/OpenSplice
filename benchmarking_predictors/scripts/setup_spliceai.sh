#!/usr/bin/env bash
set -euo pipefail
python -m venv .venv-spliceai
source .venv-spliceai/bin/activate
python -m pip install --upgrade pip
pip install -r requirements_spliceai.txt
