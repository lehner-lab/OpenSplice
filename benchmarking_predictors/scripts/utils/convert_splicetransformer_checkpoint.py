#@title Convert SpliceTransformer Lightning checkpoint to plain PyTorch state_dict
"""
Convert a SpliceTransformer Lightning checkpoint into a plain PyTorch state_dict.

This is a one-time setup utility for users who have downloaded the original
SpliceTransformer checkpoint. It is kept separate from the inference script so
that inference remains a simple, reproducible analysis step.

The original manuscript-analysis notebook used:

    /content/SpliceTransformer/model/weights/SpTransformer_pytorch (3).ckpt

and wrote:

    /content/SpliceTransformer/model/weights/model_splice_transformer.pt

Adjust the paths below to match where the SpliceTransformer repository and
weights are stored locally.
"""
from __future__ import annotations

from pathlib import Path

import torch


# ============================================================
# INPUT / OUTPUT PATHS
# ============================================================

# /content/SpliceTransformer/model/weights/SpTransformer_pytorch (3).ckpt
CKPT_PATH = Path("external/SpliceTransformer/model/weights/SpTransformer_pytorch.ckpt")

# Output path used by the inference environment.
PT_PATH = Path("external/SpliceTransformer/model/weights/model_splice_transformer.pt")
PT_PATH.parent.mkdir(parents=True, exist_ok=True)


# ============================================================
# CONVERT CHECKPOINT
# ============================================================

print(f"Loading checkpoint: {CKPT_PATH}")

# weights_only=False is required for Lightning checkpoints.
ckpt = torch.load(CKPT_PATH, map_location="cpu", weights_only=False)

# Lightning checkpoints usually store model weights under "state_dict".
state_dict = ckpt["state_dict"] if "state_dict" in ckpt else ckpt

# Remove "model." prefix if present.
new_state_dict = {}
for key, value in state_dict.items():
    new_key = key.replace("model.", "") if key.startswith("model.") else key
    new_state_dict[new_key] = value

torch.save(new_state_dict, PT_PATH)

print(f"Saved converted state_dict to: {PT_PATH}")
