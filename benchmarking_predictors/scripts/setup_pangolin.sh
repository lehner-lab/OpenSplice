#!/usr/bin/env bash
set -euo pipefail
python -m venv .venv-pangolin
source .venv-pangolin/bin/activate
python -m pip install --upgrade pip
pip install -r requirements_pangolin.txt
