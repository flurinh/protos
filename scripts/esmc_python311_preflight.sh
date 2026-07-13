#!/usr/bin/env bash
set -euo pipefail

# Disposable, no-weight production registration check.  Config metadata is
# fetched, but no model/tokenizer weights are downloaded or allocated.
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
VENV=$(mktemp -d "${TMPDIR:-/tmp}/protos-esmc-preflight.XXXXXX")
trap 'rm -rf "$VENV"' EXIT
PYTHON311=${PYTHON311:-python3.11}
command -v "$PYTHON311" >/dev/null || PYTHON311=/data/fast/envs/lambda/bin/python3.11
"$PYTHON311" -m venv --system-site-packages "$VENV"
"$VENV/bin/python" -m pip install --upgrade pip
"$VENV/bin/python" -m pip install \
  --ignore-installed \
  'transformers @ git+https://github.com/Biohub/transformers.git@ef32577f55da19a4989cd7b22e004dc43a4998cb' \
  'huggingface_hub>=0.20'
"$VENV/bin/python" -m pip install --no-deps 'accelerate>=1.0.0'
PYTHONPATH="$ROOT/src" ROOT="$ROOT" "$VENV/bin/python" - <<'PY'
import json, os
import sys
import torch
from importlib.metadata import version
from transformers import AutoConfig, AutoModelForMaskedLM
import importlib.util, sys
spec = importlib.util.spec_from_file_location("esmc_adapter", os.path.join(os.environ["ROOT"], "src/protos/processing/embedding/esmc_adapter.py"))
module = importlib.util.module_from_spec(spec); sys.modules[spec.name] = module; spec.loader.exec_module(module)
ESMCModelProvenance, ESMCLoadPolicy, HuggingFaceESMCBackend = module.ESMCModelProvenance, module.ESMCLoadPolicy, module.HuggingFaceESMCBackend

model = ESMCModelProvenance()
config = AutoConfig.from_pretrained(model.model_id, revision=model.model_revision, trust_remote_code=False)
assert config.model_type == "esmc" and config.d_model == 2560 and config.n_layers == 80
model_class = AutoModelForMaskedLM._model_mapping[type(config)]
assert model_class.__name__
kwargs = HuggingFaceESMCBackend.model_load_kwargs(model, ESMCLoadPolicy(), "cpu", torch_dtype="bfloat16")
assert kwargs["revision"] == model.model_revision and kwargs["trust_remote_code"] is False
print(json.dumps({"python": os.sys.version.split()[0], "transformers": version("transformers"), "accelerate": version("accelerate"), "model_type": config.model_type, "registered_model": model_class.__name__, "weights": "not loaded", "loader_kwargs": kwargs}, sort_keys=True))
PY
